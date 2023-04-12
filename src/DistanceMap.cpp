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

	for (auto a : p.atoms())
		atoms.push_back(a.get_row());

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

	pMin -= mMaxDistance; // extend bounding box
	pMax += mMaxDistance;

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
		
		AddDistancesForAtoms(rAtoms, rAtoms, dist);

		auto &&[center, radius] = calculateCenterAndRadius(rAtoms);
		residues.emplace_back(center, radius, std::move(rAtoms));
	}

	// treat waters special
	auto water_entity_id = db["entity"].find1<std::optional<std::string>>("type"_key == "water", "id");
	if (water_entity_id.has_value())
	{
		for (size_t i = 0; i < dim; ++i)
		{
			if (atoms[i]["label_entity_id"] == *water_entity_id)
			{
				auto pt = locations[i];
				residues.emplace_back(pt, 0.f, std::vector<std::tuple<size_t,point>>{ { i, pt } });
			}
		}
	}

	// loop over pdbx_nonpoly_scheme
	for (const auto &[asymID, entityID] : db["pdbx_nonpoly_scheme"].rows<std::string, std::string>("asym_id", "entity_id"))
	{
		if (water_entity_id.has_value() and entityID == *water_entity_id)
			continue;

		std::vector<std::tuple<size_t,point>> rAtoms;
		for (size_t i = 0; i < dim; ++i)
		{
			if (atoms[i]["label_asym_id"] == asymID)
				rAtoms.emplace_back(i, locations[i]);
		}
		
		AddDistancesForAtoms(rAtoms, rAtoms, dist);

		auto &&[center, radius] = calculateCenterAndRadius(rAtoms);
		residues.emplace_back(center, radius, std::move(rAtoms));
	}

	// loop over pdbx_branch_scheme
	for (const auto &[asym_id, pdb_seq_num] : db["pdbx_branch_scheme"].rows<std::string, std::string>("asym_id", "num"))
	{
		std::vector<std::tuple<size_t,point>> rAtoms;
		for (size_t i = 0; i < dim; ++i)
		{
			if (atoms[i]["label_asym_id"] == asym_id and atoms[i]["auth_seq_id"] == pdb_seq_num)
				rAtoms.emplace_back(i, locations[i]);
		}
		
		AddDistancesForAtoms(rAtoms, rAtoms, dist);

		auto &&[center, radius] = calculateCenterAndRadius(rAtoms);
		residues.emplace_back(center, radius, std::move(rAtoms));
	}

	cif::Progress progress((residues.size() - 1) * (residues.size() - 2), "Creating distance map");

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
				AddDistancesForAtoms(atomsI, atomsJ, dist);
				continue;
			}

			const auto &[ds, symop] = closestSymmetryCopy(spacegroup, cell, centerI, centerJ);
			if (ds - radiusI - radiusJ < mMaxDistance)
				AddDistancesForAtoms(atomsI, atomsJ, dist, symop);
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

cif::point DistanceMap::offsetToOrigin(const cif::point &p) const
{
	cif::point d{};

	while (p.m_x + d.m_x < (cell.a() / 2))
		d.m_x += cell.a();
	while (p.m_x + d.m_x > (cell.a() / 2))
		d.m_x -= cell.a();

	while (p.m_y + d.m_y < (cell.b() / 2))
		d.m_y += cell.b();
	while (p.m_y + d.m_y > (cell.b() / 2))
		d.m_y -= cell.b();

	while (p.m_z + d.m_z < (cell.c() / 2))
		d.m_z += cell.c();
	while (p.m_z + d.m_z > (cell.c() / 2))
		d.m_z -= cell.c();

	return d;
};

// --------------------------------------------------------------------

void DistanceMap::AddDistancesForAtoms(const std::vector<std::tuple<size_t,point>> &a, const std::vector<std::tuple<size_t,point>> &b, DistMap &dm)
{
	for (const auto &[ixa, loc_a] : a)
	{
		for (const auto &[ixb, loc_b] : b)
		{
			if (ixa == ixb)
				continue;

			float d = cif::distance_squared(loc_a, loc_b);

			if (d > mMaxDistanceSQ)
				continue;

			d = std::sqrt(d);

			dm[std::make_tuple(ixa, ixb)] = std::make_tuple(d, sym_op{}, false);
			dm[std::make_tuple(ixb, ixa)] = std::make_tuple(d, sym_op{}, false);
		}
	}
}

void DistanceMap::AddDistancesForAtoms(const std::vector<std::tuple<size_t,point>> &a, const std::vector<std::tuple<size_t,point>> &b,
	DistMap &dm, sym_op symop)
{
	for (const auto &[ixa, loc_a] : a)
	{
		for (const auto &[ixb, loc_b] : b)
		{
			if (ixa == ixb)
				continue;

			float d = cif::distance_squared(loc_a, symmetryCopy(loc_b, spacegroup, cell, symop));

			if (d > mMaxDistanceSQ)
				continue;

			d = std::sqrt(d);

			dm[std::make_tuple(ixa, ixb)] = std::make_tuple(d, symop, false);
			dm[std::make_tuple(ixb, ixa)] = std::make_tuple(d, symop, true);
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
		throw std::out_of_range("atom " + a + " not found in distance map");
	}

	try
	{
		ixb = index.at(b);
	}
	catch (const std::out_of_range &ex)
	{
		throw std::out_of_range("atom " + b + " not found in distance map");
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

std::vector<cif::mm::atom> DistanceMap::near(const cif::mm::atom &atom, float maxDistance) const
{
	using namespace cif::literals;

	assert(maxDistance <= mMaxDistance);
	if (maxDistance > mMaxDistance)
		throw std::runtime_error("Invalid max distance in DistanceMap::near");

	std::string a_id = atom.id();
	std::string alta = atom.get_label_alt_id();

	size_t ixa;
	try
	{
		ixa = index.at(a_id);
	}
	catch (const std::out_of_range &ex)
	{
		throw std::runtime_error("atom " + a_id + " not found in distance map");
	}

	std::vector<cif::mm::atom> result;

	auto &atom_site = atom.get_row().get_category();

	for (size_t i = mIA[ixa]; i < mIA[ixa + 1]; ++i)
	{
		const auto &[d, symop, inverse] = mA[i];

		if (d > maxDistance)
			continue;

		size_t ixb = mJA[i];

		std::string b_id = rIndex.at(ixb);
		auto altb = atom_site.find1<std::string>("id"_key == b_id, "label_alt_id");

		if (altb != alta and not altb.empty() and not alta.empty())
			continue;

		auto atom_b = mStructure.get_atom_by_id(b_id);

		if (symop)
		{
			if (inverse)
			{
				auto rt = symop.toClipperOrth(spacegroup, cell);
				result.emplace_back(symmetryCopy(atom_b, spacegroup, cell, rt.inverse()));
			}
			else
				result.emplace_back(symmetryCopy(atom_b, spacegroup, cell, symop));
		}
		else
			result.emplace_back(atom_b);
	}

	return result;
}

} // namespace pdb_redo
