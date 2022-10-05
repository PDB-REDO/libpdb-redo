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

#include <atomic>
#include <mutex>

#include <cif++/utilities.hpp>

#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/Symmetry-2.hpp"

namespace pdb_redo
{

using cif::point;

// --------------------------------------------------------------------

std::vector<clipper::RTop_orth> DistanceMap::AlternativeSites(const clipper::Spacegroup &spacegroup,
	const clipper::Cell &cell)
{
	std::vector<clipper::RTop_orth> result;

	// to make the operation at index 0 equal to identity
	result.push_back(clipper::RTop_orth::identity());

	for (int i = 0; i < spacegroup.num_symops(); ++i)
	{
		const auto &symop = spacegroup.symop(i);

		for (int u : { -1, 0, 1 })
			for (int v : { -1, 0, 1 })
				for (int w : { -1, 0, 1 })
				{
					if (i == 0 and u == 0 and v == 0 and w == 0)
						continue;

					auto rtop = clipper::RTop_frac(
						symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w))
					                .rtop_orth(cell);

					result.push_back(std::move(rtop));
				}
	}

	return result;
}

// --------------------------------------------------------------------

std::tuple<point, float> calculateCenterAndRadius(const std::vector<std::tuple<size_t,point>> &atoms)
{
	std::vector<point> pts;
	for (const auto &[ix, pt] : atoms)
		pts.emplace_back(pt);

	auto center = centroid(pts);
	float radius = 0;

	for (auto &pt : pts)
	{
		float d = static_cast<float>(distance(pt, center));
		if (radius < d)
			radius = d;
	}

	return std::make_tuple(center, radius);
}

// --------------------------------------------------------------------

DistanceMap::DistanceMap(const cif::mm::structure &p, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	float maxDistance)
	: mStructure(p)
	, cell(cell)
	, spacegroup(spacegroup)
	, dim(0)
	, mMaxDistance(maxDistance)
	, mMaxDistanceSQ(maxDistance * maxDistance)
{
	using namespace cif::literals;

	// First collect the atoms from the datablock
	std::vector<cif::row_handle> atoms;

	auto &db = p.get_datablock();

	for (auto rh : db["atom_site"].find("pdbx_model_nr"_key == p.get_model_nr() or "pdbx_model_nr"_key == cif::null))
		atoms.push_back(rh);

	dim = uint32_t(atoms.size());

	std::vector<point> locations(dim);

	// bounding box
	point pMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
		pMax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());

	for (auto &atom : atoms)
	{
		const auto &[id, x, y, z] = atom.get<std::string, float, float, float>("id", "Cartn_x", "Cartn_y", "Cartn_z");

		size_t ix = index.size();
		index[id] = ix;
		rIndex[ix] = id;

		cif::point pt{ x, y, z };
		locations[ix] = pt;

		if (pMin.m_x > pt.m_x)
			pMin.m_x = pt.m_x;
		if (pMin.m_y > pt.m_y)
			pMin.m_y = pt.m_y;
		if (pMin.m_z > pt.m_z)
			pMin.m_z = pt.m_z;

		if (pMax.m_x < pt.m_x)
			pMax.m_x = pt.m_x;
		if (pMax.m_y < pt.m_y)
			pMax.m_y = pt.m_y;
		if (pMax.m_z < pt.m_z)
			pMax.m_z = pt.m_z;
	};

	// correct locations so that the median of x, y and z are inside the cell
	std::vector<float> c(locations.size());
	auto median = [&]()
	{
		return dim % 1 == 0
		           ? c[dim / 2]
		           : (c[dim / 2 - 1] + c[dim / 2]) / 2;
	};

	transform(locations.begin(), locations.end(), c.begin(), [](point &l)
		{ return static_cast<float>(l.m_x); });
	sort(c.begin(), c.end());
	float mx = median();

	transform(locations.begin(), locations.end(), c.begin(), [](point &l)
		{ return static_cast<float>(l.m_y); });
	sort(c.begin(), c.end());
	float my = median();

	transform(locations.begin(), locations.end(), c.begin(), [](point &l)
		{ return static_cast<float>(l.m_z); });
	sort(c.begin(), c.end());
	float mz = median();

	if (cif::VERBOSE > 1)
		std::cerr << "median position of atoms: " << point(mx, my, mz) << std::endl;

	auto calculateD = [&](float m, float c)
	{
		float d = 0;
		while (m + d < -(c / 2))
			d += c;
		while (m + d > (c / 2))
			d -= c;
		return d;
	};

	mD.m_x = calculateD(mx, static_cast<float>(cell.a()));
	mD.m_y = calculateD(my, static_cast<float>(cell.b()));
	mD.m_z = calculateD(mz, static_cast<float>(cell.c()));

	clipper::Coord_orth D = toClipper(mD);

	if (mD.m_x != 0 or mD.m_y != 0 or mD.m_z != 0)
	{
		if (cif::VERBOSE > 1)
			std::cerr << "moving coorinates by " << mD.m_x << ", " << mD.m_y << " and " << mD.m_z << std::endl;

		for_each(locations.begin(), locations.end(), [&](point &p)
			{ p += mD; });
	}

	pMin -= mMaxDistance; // extend bounding box
	pMax += mMaxDistance;

	mRtOrth = AlternativeSites(spacegroup, cell);

	DistMap dist;

	std::vector<std::tuple<cif::point, float, std::vector<std::tuple<size_t,point>>>> residues;

	// loop over poly_seq_scheme
	for (const auto &[asymID, seqID] : db["pdbx_poly_seq_scheme"].rows<std::string, int>("asym_id", "seq_id"))
	{
		std::vector<std::tuple<size_t,point>> rAtoms;
		for (size_t i = 0; i < dim; ++i)
		{
			if (atoms[i]["label_asym_id"] == asymID and atoms[i]["label_seq_id"] == seqID)
				rAtoms.emplace_back(i, locations[i]);
		}
		
		AddDistancesForAtoms(rAtoms, rAtoms, dist, 0);

		auto &&[center, radius] = calculateCenterAndRadius(rAtoms);
		residues.emplace_back(center, radius, std::move(rAtoms));
	}

	// loop over pdbx_nonpoly_scheme
	for (const auto& [ asymID ] : db["pdbx_nonpoly_scheme"].rows<std::string>("asym_id"))
	{
		std::vector<std::tuple<size_t,point>> rAtoms;
		for (size_t i = 0; i < dim; ++i)
		{
			if (atoms[i]["label_asym_id"] == asymID)
				rAtoms.emplace_back(i, locations[i]);
		}
		
		AddDistancesForAtoms(rAtoms, rAtoms, dist, 0);

		auto &&[center, radius] = calculateCenterAndRadius(rAtoms);
		residues.emplace_back(center, radius, std::move(rAtoms));
	}

	// loop over pdbx_branch_scheme
	for (const auto &[asym_id, pdb_seq_num] : db["pdbx_branch_scheme"].rows<std::string, std::string>("asym_id", "pdb_seq_num"))
	{
		std::vector<std::tuple<size_t,point>> rAtoms;
		for (size_t i = 0; i < dim; ++i)
		{
			if (atoms[i]["label_asym_id"] == asym_id and atoms[i]["label_seq_id"] == pdb_seq_num)
				rAtoms.emplace_back(i, locations[i]);
		}
		
		AddDistancesForAtoms(rAtoms, rAtoms, dist, 0);

		auto &&[center, radius] = calculateCenterAndRadius(rAtoms);
		residues.emplace_back(center, radius, std::move(rAtoms));
	}

	cif::Progress progress(residues.size() * (residues.size() - 1), "Creating distance map");

	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		const auto &[centerI, radiusI, atomsI] = residues[i];

		for (size_t j = i + 1; j < residues.size(); ++j)
		{
			progress.consumed(1);

			const auto &[centerJ, radiusJ, atomsJ] = residues[j];

			// first case, no symmetry operations

			auto d = distance(centerI, centerJ) - radiusI - radiusJ;
			if (d < mMaxDistance)
			{
				AddDistancesForAtoms(atomsI, atomsJ, dist, 0);
				continue;
			}

			// now try all symmetry operations to see if we can move rj close to ri

			clipper::Coord_orth cI = toClipper(centerI);
			clipper::Coord_orth cJ = toClipper(centerJ);

			auto minR2 = d;

			int32_t kbest = 0;
			for (int32_t k = 1; k < static_cast<int32_t>(mRtOrth.size()); ++k)
			{
				auto &rt = mRtOrth[k];

				auto pJ = (cJ + D).transform(rt) - D;
				double r2 = std::sqrt((cI - pJ).lengthsq()) - radiusI - radiusJ;

				if (minR2 > r2)
				{
					minR2 = static_cast<float>(r2);
					kbest = k;
				}
			}

			if (minR2 < mMaxDistance)
				AddDistancesForAtoms(atomsI, atomsJ, dist, kbest);
		}
	}

	// Store as a sparse CSR compressed matrix

	size_t nnz = dist.size();
	mA.reserve(nnz);
	mIA.reserve(dim + 1);
	mJA.reserve(nnz);

	size_t lastR = 0;
	mIA.push_back(0);

	for (const auto &[key, value] : dist)
	{
		size_t col, row;
		std::tie(row, col) = key;

		if (row != lastR) // new row
		{
			for (size_t ri = lastR; ri < row; ++ri)
				mIA.push_back(mA.size());
			lastR = row;
		}

		mA.push_back(value);
		mJA.push_back(col);
	}

	for (size_t ri = lastR; ri < dim; ++ri)
		mIA.push_back(mA.size());
}

// --------------------------------------------------------------------

void DistanceMap::AddDistancesForAtoms(const std::vector<std::tuple<size_t,point>> &a, const std::vector<std::tuple<size_t,point>> &b, DistMap &dm, int32_t rtix)
{
	for (const auto &[ixa, loc_a] : a)
	{
		clipper::Coord_orth pa = toClipper(loc_a);

		for (const auto &[ixb, loc_b] : b)
		{
			if (ixa == ixb)
				continue;

			clipper::Coord_orth pb = toClipper(loc_b);

			if (rtix)
				pb = pb.transform(mRtOrth[rtix]);

			auto d = static_cast<float>((pa - pb).lengthsq());
			if (d > mMaxDistanceSQ)
				continue;

			d = std::sqrt(d);

			dm[std::make_tuple(ixa, ixb)] = std::make_tuple(d, rtix);
			dm[std::make_tuple(ixb, ixa)] = std::make_tuple(d, -rtix);
		}
	}
}

float DistanceMap::operator()(const std::string &a, const std::string &b) const
{
	size_t ixa, ixb;

	try
	{
		ixa = index.at(a);
	}
	catch (const std::out_of_range &ex)
	{
		throw std::runtime_error("atom " + a + " not found in distance map");
	}

	try
	{
		ixb = index.at(b);
	}
	catch (const std::out_of_range &ex)
	{
		throw std::runtime_error("atom " + b + " not found in distance map");
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

// std::vector<std::string> DistanceMap::near(cif::row_handle atom, float maxDistance) const
// {
// 	using namespace cif::literals;

// 	assert(maxDistance <= mMaxDistance);
// 	if (maxDistance > mMaxDistance)
// 		throw std::runtime_error("Invalid max distance in DistanceMap::near");

// 	const auto &[a_id, alta] = atom.get<std::string, std::string>("id", "label_alt_id");

// 	size_t ixa;
// 	try
// 	{
// 		ixa = index.at(a_id);
// 	}
// 	catch (const std::out_of_range &ex)
// 	{
// 		throw std::runtime_error("atom " + a_id + " not found in distance map");
// 	}

// 	std::vector<std::string> result;
// 	auto &atom_site = db["atom_site"];

// 	for (size_t i = mIA[ixa]; i < mIA[ixa + 1]; ++i)
// 	{
// 		float d;
// 		int32_t rti;
// 		std::tie(d, rti) = mA[i];

// 		if (d > maxDistance)
// 			continue;

// 		size_t ixb = mJA[i];

// 		std::string b_id = rIndex.at(ixb);
// 		auto altb = atom_site.find1<std::string>("id"_key == b_id, "label_alt_id");

// 		if (altb != alta and not altb.empty() and not alta.empty())
// 			continue;

// 		if (rti > 0)
// 			result.emplace_back(symmetryCopy(b_id, mD, spacegroup, cell, mRtOrth.at(rti)));
// 		else if (rti < 0)
// 			result.emplace_back(symmetryCopy(b_id, mD, spacegroup, cell, mRtOrth.at(-rti).inverse()));
// 		else
// 			result.emplace_back(b_id);
// 	}

// 	return result;
// }

} // namespace pdb_redo
