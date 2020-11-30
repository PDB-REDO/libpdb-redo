/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "config.hpp"

#include <atomic>
#include <mutex>

#include "cif++/CifUtils.hpp"

#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/Symmetry-2.hpp"

//#define DEBUG_VOOR_BART

namespace mmcif
{

// --------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& os, const Atom& a)
{
	os << a.labelAsymID() << ':' << a.labelSeqID() << '/' << a.labelAtomID() << a.labelAltID();
	
	return os;
}

// --------------------------------------------------------------------

std::vector<clipper::RTop_orth> DistanceMap::AlternativeSites(const clipper::Spacegroup& spacegroup,
	const clipper::Cell& cell)
{
	std::vector<clipper::RTop_orth> result;
	
	// to make the operation at index 0 equal to identity
	result.push_back(clipper::RTop_orth::identity());

	for (int i = 0; i < spacegroup.num_symops(); ++i)
	{
		const auto& symop = spacegroup.symop(i);
		
		for (int u: { -1, 0, 1})
			for (int v: { -1, 0, 1})
				for (int w: { -1, 0, 1})
				{
					if (i == 0 and u == 0 and v == 0 and w == 0)
						continue;

					auto rtop = clipper::RTop_frac(
							symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w)
						).rtop_orth(cell);

					result.push_back(std::move(rtop));
				}
	}
	
	return result;
}

// --------------------------------------------------------------------

DistanceMap::DistanceMap(const Structure& p, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell,
		float maxDistance)
	: structure(p), cell(cell), spacegroup(spacegroup), dim(0), mMaxDistance(maxDistance), mMaxDistanceSQ(maxDistance * maxDistance)
{
	auto& atoms = p.atoms();
	dim = atoms.size();
	
	std::vector<clipper::Coord_orth> locations(dim);
	
	// bounding box
	Point pMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
		  pMax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
	
	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
		rIndex[ix] = atom.id();
		
		locations[ix] = toClipper(atom.location());
		
		auto p = atom.location();

		if (pMin.mX > p.mX)
			pMin.mX = p.mX;
		if (pMin.mY > p.mY)
			pMin.mY = p.mY;
		if (pMin.mZ > p.mZ)
			pMin.mZ = p.mZ;

		if (pMax.mX < p.mX)
			pMax.mX = p.mX;
		if (pMax.mY < p.mY)
			pMax.mY = p.mY;
		if (pMax.mZ < p.mZ)
			pMax.mZ = p.mZ;
	};
	
	// correct locations so that the median of x, y and z are inside the cell
	std::vector<float> c(locations.size());
	auto median = [&]()
	{
		return dim % 1 == 0
			? c[dim / 2]
			: (c[dim / 2 - 1] + c[dim / 2]) / 2;
	};
	
	transform(locations.begin(), locations.end(), c.begin(), [](auto& l) { return l[0]; });
	sort(c.begin(), c.end());
	float mx = median();
	
	transform(locations.begin(), locations.end(), c.begin(), [](auto& l) { return l[1]; });
	sort(c.begin(), c.end());
	float my = median();
	
	transform(locations.begin(), locations.end(), c.begin(), [](auto& l) { return l[2]; });
	sort(c.begin(), c.end());
	float mz = median();

	if (cif::VERBOSE > 1)
		std::cerr << "median position of atoms: " << Point(mx, my, mz) << std::endl;
	
	auto calculateD = [&](float m, float c)
	{
		float d = 0;
		while (m + d < -(c / 2))
			d += c;
		while (m + d > (c / 2))
			d -= c;
		return d;
	};

	mD.mX = calculateD(mx, cell.a());
	mD.mY = calculateD(my, cell.b());
	mD.mZ = calculateD(mz, cell.c());
	
	clipper::Coord_orth D = toClipper(mD);
	
	if (mD.mX != 0 or mD.mY != 0 or mD.mZ != 0)
	{
		if (cif::VERBOSE)
			std::cerr << "moving coorinates by " << mD.mX << ", " << mD.mY << " and " << mD.mZ << std::endl;
		
		for_each(locations.begin(), locations.end(), [&](auto& p) { p += toClipper(mD); });
	}
	
	pMin -= mMaxDistance;	// extend bounding box
	pMax += mMaxDistance;

	mRtOrth = AlternativeSites(spacegroup, cell);
	
	DistMap dist;
	
	std::vector<const Residue*> residues;

	for (auto& poly: p.polymers())
	{
		for (auto& m: poly)
		{
			residues.emplace_back(&m);

			// Add distances for atoms in this residue		
			AddDistancesForAtoms(m, m, dist, 0);
		}
	}

	for (auto& r: p.getNonPolymers())
	{
		residues.emplace_back(&r);

		// Add distances for atoms in this residue		
		AddDistancesForAtoms(r, r, dist, 0);
	}

	for (auto& r: p.getBranchResidues())
	{
		residues.emplace_back(&r);

		// Add distances for atoms in this residue		
		AddDistancesForAtoms(r, r, dist, 0);
	}

	cif::Progress progress(residues.size() * residues.size(), "Creating distance std::map");
	
	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		auto& ri = *residues[i];
		
		Point centerI;
		float radiusI;
		std::tie(centerI, radiusI) = ri.centerAndRadius();
		
		for (size_t j = i + 1; j < residues.size(); ++j)
		{
			progress.consumed(1);

			auto& rj = *residues[j];
			
			// first case, no symmetry operations
			
			Point centerJ;
			float radiusJ;
			std::tie(centerJ, radiusJ) = rj.centerAndRadius();
			
			auto d = Distance(centerI, centerJ) - radiusI - radiusJ;
			if (d < mMaxDistance)
			{
				AddDistancesForAtoms(ri, rj, dist, 0);
				continue;
			}
			
			// now try all symmetry operations to see if we can move rj close to ri
			
			clipper::Coord_orth cI = toClipper(centerI);
			clipper::Coord_orth cJ = toClipper(centerJ);
			
			auto minR2 = d;
			
			int32_t kbest = 0;
			for (int32_t k = 1; k < static_cast<int32_t>(mRtOrth.size()); ++k)
			{
				auto& rt = mRtOrth[k];
				
				auto pJ = (cJ + D).transform(rt) - D;
				double r2 = sqrt((cI - pJ).lengthsq()) - radiusI - radiusJ;

				if (minR2 > r2)
				{
					minR2 = r2;
					kbest = k;
				}
			}
			
			if (minR2 < mMaxDistance)
				AddDistancesForAtoms(ri, rj, dist, kbest);
		}
	}
	
	// Store as a sparse CSR compressed matrix
	
	size_t nnz = dist.size();
	mA.reserve(nnz);
	mIA.reserve(dim + 1);
	mJA.reserve(nnz);
	
	size_t lastR = 0;
	mIA.push_back(0);
	
	for (auto& di: dist)
	{
		size_t c, r;
		std::tie(r, c) = di.first;

		if (r != lastR)	// new row
		{
			for (size_t ri = lastR; ri < r; ++ri)
				mIA.push_back(mA.size());
			lastR = r;
		}

		mA.push_back(di.second);
		mJA.push_back(c);
	}

	for (size_t ri = lastR; ri < dim; ++ri)
		mIA.push_back(mA.size());
}

// --------------------------------------------------------------------

void DistanceMap::AddDistancesForAtoms(const Residue& a, const Residue& b, DistMap& dm, int32_t rtix)
{
	for (auto& aa: a.atoms())
	{
		clipper::Coord_orth pa = toClipper(aa.location());
		size_t ixa = index[aa.id()];
		
		for (auto& bb: b.atoms())
		{
			if (aa.id() == bb.id())
				continue;
			
			clipper::Coord_orth pb = toClipper(bb.location());
			
			if (rtix)
				pb = pb.transform(mRtOrth[rtix]);
			
			auto d = (pa - pb).lengthsq();
			if (d > mMaxDistanceSQ)
				continue;

			d = sqrt(d);

			size_t ixb = index[bb.id()];

			dm[std::make_tuple(ixa, ixb)] = std::make_tuple(d, rtix);
			dm[std::make_tuple(ixb, ixa)] = std::make_tuple(d, -rtix);
		}
	}
}

float DistanceMap::operator()(const Atom& a, const Atom& b) const
{
	size_t ixa, ixb;
	
	try
	{
		ixa = index.at(a.id());
	}
	catch (const std::out_of_range& ex)
	{
		throw std::runtime_error("atom " + a.id() + " not found in distance std::map");
	}
		
	try
	{
		ixb = index.at(b.id());
	}
	catch (const std::out_of_range& ex)
	{
		throw std::runtime_error("atom " + b.id() + " not found in distance std::map");
	}
	
//	if (ixb < ixa)
//		std::swap(ixa, ixb);

	size_t L = mIA[ixa];
	size_t R = mIA[ixa + 1] - 1;
	
	while (L <= R)
	{
		size_t i = (L + R) / 2;

		if (mJA[i] == ixb)
			return std::get<0>(mA[i]);

		if (mJA[i] < ixb)
			L = i + 1;
		else
			R = i - 1;
	}
	
	return 100.f;
}

std::vector<Atom> DistanceMap::near(const Atom& a, float maxDistance) const
{
	assert(maxDistance <= mMaxDistance);
	if (maxDistance > mMaxDistance)
		throw std::runtime_error("Invalid max distance in DistanceMap::near");
	
	size_t ixa;
	try
	{
		ixa = index.at(a.id());
	}
	catch (const std::out_of_range& ex)
	{
		throw std::runtime_error("atom " + a.id() + " not found in distance map");
	}

	std::vector<Atom> result;
	auto alta = a.labelAltID();
	
	for (size_t i = mIA[ixa]; i < mIA[ixa + 1]; ++i)
	{
		float d;
		int32_t rti;
		std::tie(d, rti) = mA[i];

		if (d > maxDistance)
			continue;

		size_t ixb = mJA[i];
		Atom b = structure.getAtomByID(rIndex.at(ixb));

		auto altb = b.labelAltID();
		if (altb != alta and not altb.empty() and not alta.empty())
			continue;
		
		if (rti > 0)
			result.emplace_back(symmetryCopy(b, mD, spacegroup, cell, mRtOrth.at(rti)));
		else if (rti < 0)
			result.emplace_back(symmetryCopy(b, mD, spacegroup, cell, mRtOrth.at(-rti).inverse()));
		else
			result.emplace_back(b);
	}
	
	return result;
}

}
