/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "pdb-redo/BlobFinder.hpp"

#include <algorithm>
#include <cif++/point.hpp>
#include <cif++/symmetry.hpp>
#include <limits>
#include <pdb-redo/Restraints.hpp>
#include <ranges>
#include <stdexcept>

namespace pdb_redo
{

BlobFinder::BlobFinder(clipper::Xmap<float> &xmm, cif::mm::structure &structure, float growingPercentile)
	: mXmap(xmm)
	, mProteinAtoms(structure.atoms())
	, mCrystal(structure.get_datablock())
{
	// To make sure we iterate through density around the protein
	// we intend to make a cuboid around the protein by taking the min and max coordinate in x,y and z
	// and extend with 6 angstrom in each direction

	cif::point min, max;
	std::vector<cif::point> pts;
	pts.reserve(mProteinAtoms.size());

	for (bool first = true; auto &a : mProteinAtoms)
	{
		auto loc = a.get_location();

		pts.emplace_back(loc);

		if (std::exchange(first, false))
			min = max = loc;
		else
		{
			if (min.x > loc.x)
				min.x = loc.x;
			if (min.y > loc.y)
				min.y = loc.y;
			if (min.z > loc.z)
				min.z = loc.z;

			if (max.x < loc.x)
				max.x = loc.x;
			if (max.y < loc.y)
				max.y = loc.y;
			if (max.z < loc.z)
				max.z = loc.z;
		}
	}

	std::tie(mProteinCenter, mProteinRadius) = cif::smallest_sphere_around_points(pts);

	// Store all residue spheres as well
	for (auto &poly : structure.polymers())
	{
		for (auto &res : poly)
			mResidueSpheres.emplace_back(res.center_and_radius());
	}

	// use radius
	float max_r_sq = mProteinRadius * mProteinRadius;

	using namespace clipper;

	float extend = 6;
	cif::point extending{ extend, extend, extend };
	cif::point pMin = min - extending, pMax = max + extending;
	Coord_orth oMin{ pMin.x, pMin.y, pMin.z }, oMax{ pMax.x, pMax.y, pMax.z };
	Coord_frac fMin = oMin.coord_frac(mXmap.cell()), fMax = oMax.coord_frac(mXmap.cell());
	Coord_map mMin = fMin.coord_map(mXmap.grid_sampling()), mMax = fMax.coord_map(mXmap.grid_sampling());
	Coord_grid gMin = mMin.floor(), gMax = mMax.ceil();

	// Set starting points and initialize vector of potential interesting gridpoints
	auto i0 = clipper::Xmap_base::Map_reference_coord(xmm, gMin);

	// Create vector with density heights for all values >0
	for (auto iu = i0; iu.coord().u() <= gMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gMax[2]; iw.next_w())
			{
				double dens_height = xmm[iw];
				auto iw_op = iw.coord_orth();
				if (dens_height > 0 and cif::distance_squared(cif::point{ iw_op.x(), iw_op.y(), iw_op.z() }, mProteinCenter) < max_r_sq + extend * extend)
					mPotentialGridPoints.emplace_back(iw);
			}

	// Check if vector not empty
	if (mPotentialGridPoints.empty())
		throw std::runtime_error("No gridpoints with density height above 0");

	// Sort vector on density height (from high to low numbers)
	std::ranges::sort(mPotentialGridPoints, [this](GridPoint a, GridPoint b)
		{ return mXmap[a] < mXmap[b]; });

	auto ix = static_cast<size_t>(std::ceil(growingPercentile * mPotentialGridPoints.size()));
	mGrowingThreshold = mXmap[mPotentialGridPoints.at(ix)];
	if (mGrowingThreshold == 0)
		mGrowingThreshold = 1e-6;

	mPotentialGridPoints.erase(mPotentialGridPoints.begin(), mPotentialGridPoints.begin() + ix);
}

std::vector<cif::point> BlobFinder::next(float minimalVolume)
{
	auto cellVolume = mXmap.cell().volume();
	auto gridSize = mXmap.grid_sampling().size();
	auto gridPointVolume = cellVolume / gridSize;

	while (mPotentialGridPoints.size() > 0)
	{
		auto newblob = pop();

		std::ranges::sort(newblob, [](auto &a, auto &b)
			{ return a.index() < b.index(); });

		newblob.erase(std::ranges::unique(newblob, [](auto &a, auto &b)
			{ return a.index() == b.index(); }).begin(), newblob.end());

		std::vector<cif::point> result;

		for (auto &gp : newblob)
		{
			auto op = gp.coord_orth();
			result.emplace_back(op.x(), op.y(), op.z());
		}

		// std::ranges::sort(result);
		// result.erase(std::ranges::unique(result).begin(), result.end());

		if (result.size() < 30 or gridPointVolume * result.size() < minimalVolume)
			continue;

		auto [blobCenter, blobRadius] = cif::smallest_sphere_around_points(result);

		float bestD = std::numeric_limits<float>::max();
		cif::sym_op bestSO{};

		for (auto &[c, r] : mResidueSpheres)
		{
			if (auto d = distance(c, blobCenter); bestD > d)
				bestD = d;
		}

		if (bestD > 3.0f)
		{
			for (auto &[c, r] : mResidueSpheres)
			{
				auto [d, p, so] = mCrystal.closest_symmetry_copy(c, blobCenter);

				if (bestD > d)
				{
					bestD = d;
					bestSO = so;
				}
			}

			if (bestSO)
			{
				for (auto &bp : result)
					bp = mCrystal.symmetry_copy(bp, bestSO);

				std::tie(blobCenter, blobRadius) = cif::smallest_sphere_around_points(result);
			}
		}

		// Check if found blob is in proximity of protein atoms
		if (not mProteinAtoms.empty())
		{
			auto max_d = (10 + blobRadius) * (10 + blobRadius);

			if (std::ranges::find_if(mProteinAtoms, [=](const cif::mm::atom &a)
					{ return cif::distance_squared(a.get_location(), blobCenter) < max_d; }) == mProteinAtoms.end())
				continue;
		}

		return result;
	}

	return {};
}

std::vector<BlobFinder::GridPoint> BlobFinder::pop()
{
	std::stack<GridPoint> stack{ { mPotentialGridPoints.back() } };
	mPotentialGridPoints.pop_back();

	std::vector<GridPoint> blob;

	while (not stack.empty())
	{
		auto gridpoint = stack.top();
		stack.pop();

		blob.emplace_back(gridpoint);

		// Find neighbouring gridpoints to get a 3*3*3 cube (excluding the center gridpoint)
		for (auto u = -1; u <= 1; u++)
			for (auto v = -1; v <= 1; v++)
				for (auto w = -1; w <= 1; w++)
				{
					if (u == 0 and v == 0 and w == 0)
						continue;

					auto gp = gridpoint;
					if (u < 0)
						gp = gp.prev_u();
					else if (u > 0)
						gp = gp.next_u();

					if (v < 0)
						gp = gp.prev_v();
					else if (v > 0)
						gp = gp.next_v();

					if (w < 0)
						gp = gp.prev_w();
					else if (w > 0)
						gp = gp.next_w();

					if (std::ranges::find_if(blob, [ix = gp.index()](const GridPoint &p)
							{ return p.index() == ix; }) != blob.end())
						continue;

					if (mXmap[gp] < mGrowingThreshold)
						continue;

					stack.push(gp);
				};
	}

	return blob;
}

} // namespace pdb_redo