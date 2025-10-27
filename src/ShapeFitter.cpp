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

#include "pdb-redo/ShapeFitter.hpp"

#include "cif++/point.hpp"
#include "clipper/core/coords.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace pdb_redo
{

double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap)
{
	auto cellVolume = xmap.cell().volume();
	auto gridSize = xmap.grid_sampling().size();
	auto gridPointVolume = cellVolume / gridSize;
	auto gridPointRadius = std::pow((3 * gridPointVolume) / (4 * cif::kPI), 1 / 3.0);

	// Minimal blob finding algo
	std::vector<std::pair<clipper::Coord_grid, bool>> blobPoints;

	for (clipper::Xmap<float>::Map_reference_index i(xmap); not i.last(); i.next())
	{
		if (xmap[i] <= 0)
			continue;

		blobPoints.emplace_back(i.coord(), true);
		break;
	}

	if (blobPoints.empty())
		throw std::runtime_error("No density found in map");

	std::vector<cif::point> blob;

	for (;;)
	{
		auto i = std::find_if(blobPoints.begin(), blobPoints.end(), [](auto &bp)
			{ return bp.second; });

		if (i == blobPoints.end())
			break;

		i->second = false;

		clipper::Xmap<float>::Map_reference_coord iw(xmap, i->first);
		int symNr = iw.sym();

		blob.emplace_back(iw.coord_orth());

		for (int du : { -1, 0, 1 })
			for (int dv : { -1, 0, 1 })
				for (int dw : { -1, 0, 1 })
				{
					if (du == 0 and dv == 0 and dw == 0)
						continue;

					clipper::Coord_grid g(i->first.u() + du, i->first.v() + dv, i->first.w() + dw);
					clipper::Xmap<float>::Map_reference_coord iw(xmap, g);
					if (iw.sym() != symNr)
						continue;

					if (xmap[iw] == 0)
						continue;

					auto j = std::find_if(blobPoints.begin(), blobPoints.end(), [g](auto &bp)
						{ return bp.first == g; });

					if (j == blobPoints.end())
						blobPoints.emplace_back(g, true);
				}
	}

	// // Minimal blob finding algo

	// std::vector<std::pair<clipper::Coord_grid, bool>> blobPoints;

	// for (clipper::Xmap<float>::Map_reference_index i(xmap); not i.last(); i.next())
	// {
	// 	if (xmap[i] <= 0)
	// 		continue;

	// 	blobPoints.emplace_back(i.coord(), true);
	// 	break;
	// }

	// if (blobPoints.empty())
	// 	throw std::runtime_error("No density found in map");

	// std::vector<cif::point> blob;
	// auto range = xmap.grid_asu();
	// auto gMin = range.min();
	// auto gMax = range.max();

	// for (;;)
	// {
	// 	auto i = std::find_if(blobPoints.begin(), blobPoints.end(), [](auto &bp)
	// 		{ return bp.second; });

	// 	if (i == blobPoints.end())
	// 		break;

	// 	i->second = false;

	// 	clipper::Xmap<float>::Map_reference_index iw(xmap, i->first);

	// 	blob.emplace_back(iw.coord_orth());

	// 	for (int du : { -1, 0, 1 })
	// 		for (int dv : { -1, 0, 1 })
	// 			for (int dw : { -1, 0, 1 })
	// 			{
	// 				if (du == 0 and dv == 0 and dw == 0)
	// 					continue;

	// 				auto ix = iw.index_offset(du, dv, dw);

	// 				if (xmap.get_data(ix) == 0)
	// 					continue;

	// 				clipper::Coord_grid g(i->first.u() + du, i->first.v() + dv, i->first.w() + dw);

	// 				if (g.u() < gMin.u() or g.u() > gMax.u() or
	// 					g.v() < gMin.v() or g.v() > gMax.v() or
	// 					g.w() < gMin.w() or g.w() > gMax.w())
	// 				{
	// 					continue;
	// 				}

	// 				// clipper::Coord_grid g(i->first.u() + du, i->first.v() + dv, i->first.w() + dw);

	// 				auto j = std::find_if(blobPoints.begin(), blobPoints.end(), [g](auto &bp)
	// 					{ return bp.first == g; });

	// 				if (j == blobPoints.end())
	// 					blobPoints.emplace_back(g, true);
	// 			}
	// }


	// std::vector<cif::point> blob;
	// auto range = xmap.grid_asu();
	// auto gMin = range.min();
	// auto gMax = range.max();

	// for (int u = gMin.u(); u < gMax.u(); ++u)
	// {
	// 	for (int v = gMin.v(); v < gMax.v(); ++v)
	// 	{
	// 		for (int w = gMin.w(); w < gMax.w(); ++w)
	// 		{
	// 			clipper::Coord_grid g(u, v, w);

	// 			clipper::Xmap<float>::Map_reference_coord ix(xmap, g);

	// 			if (xmap[ix] == 0)
	// 				continue;

	// 			blob.emplace_back(ix.coord_orth());
	// 		}
	// 	}
	// }

	cif::point center = cif::center_points(blob);

	using list_of_spheres = std::vector<std::tuple<cif::point, float>>;

	list_of_spheres mapSpheres;
	mapSpheres.reserve(blob.size());
	for (auto p : blob)
	{
		mapSpheres.emplace_back(p, gridPointRadius);
		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f} }},\n", p.get_x(), p.get_y(), p.get_z());
	}

	auto dots = cif::spherical_dots<12>::instance();

	std::vector<std::tuple<cif::point, float, float>> surface;
	surface.reserve(dots.size());
	for (auto p : dots)
		surface.emplace_back(p, std::numeric_limits<float>::max(), std::numeric_limits<float>::min());

	for (const auto &[sp_c, radius] : mapSpheres)
	{
		for (auto &[line, l1, l2] : surface)
		{
			auto a = cif::dot_product(line, line);
			auto b = 2 * cif::dot_product(line, -sp_c);
			auto c = cif::dot_product(-sp_c, -sp_c) - radius * radius;

			auto D = b * b - 4 * a * c;
			if (D < 0)
				continue;

			float d1, d2;

			if (D == 0)
				d1 = d2 = -b / (2 * a);
			else
			{
				d1 = (-b + std::sqrt(D)) / (2 * a);
				d2 = (-b - std::sqrt(D)) / (2 * a);

				if (d1 > d2)
					std::swap(d1, d2);
			}

			if (l1 > d1)
				l1 = d1;
			if (l2 < d2)
				l2 = d2;
		}
	}

	for (auto &[line, l1, l2] : surface)
	{
		auto p1 = center + l1 * line;
		auto p2 = center + l2 * line;

		std::cout << "p1: " << p1 << ", p2: " << p2 << "\n";
	}

	return 0;
}

} // namespace pdb_redo