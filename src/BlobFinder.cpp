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
#include <pdb-redo/Restraints.hpp>
#include <stdexcept>

namespace pdb_redo
{

// BlobFinder::BlobFinder(clipper::Xmap<float> &xmm, float growingPercentile)
// 	: mXmap(xmm)
//     , mCrystal(structure.get_datablock())
// {
// 	// Create vector with density heights for all values >0
// 	for (auto i = clipper::Xmap_base::Map_reference_coord(xmm); not i.last(); i.next())
// 	{
// 		double dens_height = xmm[i];
// 		if (dens_height > 0)
// 			mPotentialGridPoints.emplace_back(i);
// 	}

// 	// Check if vector not empty
// 	if (mPotentialGridPoints.empty())
// 		throw std::runtime_error("No gridpoints with density height above 0");

// 	// Sort vector on density height (from high to low numbers)
// 	std::ranges::sort(mPotentialGridPoints, [this](GridPoint a, GridPoint b)
// 		{ return mXmap[a] < mXmap[b]; });

// 	auto ix = static_cast<size_t>(std::ceil(growingPercentile * mPotentialGridPoints.size()));
// 	mGrowingThreshold = mXmap[mPotentialGridPoints.at(ix)];
// 	if (mGrowingThreshold == 0)
// 		mGrowingThreshold = 1e-6;

// 	mPotentialGridPoints.erase(mPotentialGridPoints.begin(), mPotentialGridPoints.begin() + ix);
// }

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
			if (min.m_x > loc.m_x)
				min.m_x = loc.m_x;
			if (min.m_y > loc.m_y)
				min.m_y = loc.m_y;
			if (min.m_z > loc.m_z)
				min.m_z = loc.m_z;

			if (max.m_x < loc.m_x)
				max.m_x = loc.m_x;
			if (max.m_y < loc.m_y)
				max.m_y = loc.m_y;
			if (max.m_z < loc.m_z)
				max.m_z = loc.m_z;
		}
	}

    std::tie(mProteinCenter, mProteinRadius) = cif::smallest_sphere_around_points(pts);

	// use radius
	float max_r_sq = mProteinRadius * mProteinRadius;

	using namespace clipper;

	float extend = 6;
	cif::point extending{ extend, extend, extend };
	cif::point pMin = min - extending, pMax = max + extending;
	Coord_orth oMin = pMin, oMax = pMax;
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
				if (dens_height > 0 and cif::distance_squared(cif::point{ iw.coord_orth() }, mProteinCenter) < max_r_sq + extend * extend)
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

		// Erase all gridpoints from the potential list, including symmetry copies
		std::erase_if(mPotentialGridPoints,
				[&newblob](const GridPoint &gp)
				{
					return std::ranges::find_if(newblob, [gp](const GridPoint &p)
							   { return gp.index() == p.index(); }) != newblob.end();
				});

		std::vector<cif::point> result;

		for (auto &gp : newblob)
			result.emplace_back(gp.coord_orth());

		std::ranges::sort(result);
		result.erase(std::ranges::unique(result).begin(), result.end());

		if (result.size() < 30 or gridPointVolume * result.size() < minimalVolume)
			continue;

        auto [blobCenter, blobRadius] = cif::smallest_sphere_around_points(result);

        if (cif::distance(mProteinCenter, blobCenter) > blobRadius + mProteinRadius)
        {
            auto [d, p, so] = mCrystal.closest_symmetry_copy(mProteinCenter, blobCenter);

            if (so)
            {
                for (auto &bp : result)
                    bp = mCrystal.symmetry_copy(bp, so);
            }
        }

		// // Check if found blob is in proximity of protein atoms
		// if (not mProteinAtoms.empty())
		// {
		// 	auto [center, radius] = cif::smallest_sphere_around_points(result);
		// 	auto max_d = (10 + radius) * (10 + radius);

		// 	if (std::ranges::find_if(mProteinAtoms, [=](const cif::mm::atom &a)
		// 			{ return cif::distance_squared(a.get_location(), center) < max_d; }) == mProteinAtoms.end())
		// 		continue;
		// }

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

		// Define coordinates of the starting gridpoint
		auto u0 = gridpoint.coord().u();
		auto v0 = gridpoint.coord().v();
		auto w0 = gridpoint.coord().w();

		// Find neighbouring gridpoints to get a 3*3*3 cube (excluding the center gridpoint)
		for (auto u = u0 - 1; u <= u0 + 1; u++)
			for (auto v = v0 - 1; v <= v0 + 1; v++)
				for (auto w = w0 - 1; w <= w0 + 1; w++)
				{
					if (u == u0 and v == v0 and w == w0)
						continue;

					auto n_gp = clipper::Xmap_base::Map_reference_coord(mXmap, { u, v, w });

					if (std::ranges::find_if(blob, [ix = n_gp.index()](const GridPoint &p)
							{ return p.index() == ix; }) != blob.end())
						continue;

					if (mXmap[n_gp] < mGrowingThreshold)
						continue;

					stack.push(n_gp);
				};
	}

	return blob;
}

} // namespace pdb_redo