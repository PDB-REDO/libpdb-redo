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

#include "cif++/atom_type.hpp"
#include "cif++/point.hpp"
#include "clipper/core/coords.h"

#include <algorithm>
#include <clipper/core/clipper_types.h>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pdb_redo
{

// Locate the single blob in the xmap
std::vector<clipper::Coord_grid> findSingleBlob(clipper::Xmap<float> &xmap)
{
	auto &sg = xmap.spacegroup();

	struct Vec3Less
	{
		bool operator()(const clipper::Vec3<int> &a, const clipper::Vec3<int> &b) const
		{
			int d = a[0] - b[0];
			if (d == 0)
				d = a[1] - b[1];
			if (d == 0)
				d = a[2] - b[2];
			return d < 0;
		}
	};

	// Minimal blob finding algo
	std::set<clipper::Coord_grid, Vec3Less> result;
	std::stack<clipper::Coord_grid> stack;

	for (clipper::Xmap<float>::Map_reference_coord i(xmap); not i.last(); i.next())
	{
		if (xmap[i] <= 0)
			continue;

		
		if (int symNr = i.sym(); symNr == 0)
			stack.push(i.coord());
		else
		{
			clipper::Coord_map cm(i.coord());
			cm = clipper::Coord_map(sg.symop(symNr).inverse() * cm);
			stack.push(cm.coord_grid());
		}

		break;
	}


	while (not stack.empty())
	{
		auto p = stack.top();
		stack.pop();

		clipper::Xmap<float>::Map_reference_coord iw(xmap, p);

		result.insert(p);

		for (int du : { -1, 0, 1 })
			for (int dv : { -1, 0, 1 })
				for (int dw : { -1, 0, 1 })
				{
					if (du == 0 and dv == 0 and dw == 0)
						continue;

					clipper::Coord_grid g(p.u() + du, p.v() + dv, p.w() + dw);
					clipper::Xmap<float>::Map_reference_coord iw(xmap, g);

					if (xmap[iw] == 0)
						continue;

					if (iw.sym() != 0)
					{
						// Move point into the correct symmetry
						clipper::Coord_map cm(g);
						cm = clipper::Coord_map(sg.symop(iw.sym()).inverse() * cm);
						g = cm.coord_grid();
					}

					if (not result.contains(g))
						stack.push(g);
				}
	}

	return { result.begin(), result.end() };
}

double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap)
{
	auto cellVolume = xmap.cell().volume();
	auto gridSize = xmap.grid_sampling().size();
	auto gridPointVolume = cellVolume / gridSize;
	auto gridPointRadius = std::pow((3 * gridPointVolume) / (4 * cif::kPI), 1 / 3.0);

	// Minimal blob finding algo
	std::vector<std::pair<clipper::Coord_grid, bool>> blobPoints;

	std::vector<cif::point> blob;
	for (auto p : findSingleBlob(xmap))
		blob.emplace_back(xmap.coord_orth(p.coord_map()));
	
	cif::point center = cif::center_points(blob);

	using list_of_spheres = std::vector<std::tuple<cif::point, float>>;

	list_of_spheres mapSpheres;
	mapSpheres.reserve(blob.size());
	for (auto p : blob)
	{
		mapSpheres.emplace_back(p, gridPointRadius * 1.5);
		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f}, r: {:.3f} }},\n", p.get_x(), p.get_y(), p.get_z(), gridPointRadius);
	}

	std::cout << "\n\n";

	auto dots = cif::spherical_dots<15>::instance();

	std::vector<std::tuple<cif::point, float, float>> blobSurface;
	blobSurface.reserve(dots.size());
	for (auto p : dots)
		blobSurface.emplace_back(p, std::numeric_limits<float>::max(), std::numeric_limits<float>::min());

	for (const auto &[sp_c, radius] : mapSpheres)
	{
		for (auto &[line, l1, l2] : blobSurface)
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

	for (auto &[line, l1, l2] : blobSurface)
	{
		auto p1 = /* center + */ l1 * line;
		auto p2 = /* center + */ l2 * line;

		// std::cout << "p1: " << p1 << ", p2: " << p2 << "\n";
		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f} }},\n", p1.get_x(), p1.get_y(), p1.get_z())
				  << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f} }},\n", p2.get_x(), p2.get_y(), p2.get_z());
	}

	// Same for the ligand

	auto &ligand = structure.get_residue(asym_id);

	std::vector<cif::point> atomLocations;
	for (auto a : ligand.atoms())
		atomLocations.emplace_back(a.get_location());
	center = cif::center_points(atomLocations);

	list_of_spheres atomSpheres;
	atomSpheres.reserve(dots.size());
	for (size_t ix = 0; auto a : ligand.atoms())
	{
		cif::point loc = atomLocations[ix++];
		auto radius = cif::atom_type_traits(a.get_type()).radius();
		atomSpheres.emplace_back(loc, radius);
		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f}, r: {:.4f} }},\n", loc.get_x(), loc.get_y(), loc.get_z(), radius);
	}

	std::cout << "\n\n";

	std::vector<std::tuple<cif::point, float, float>> ligandSurface;
	ligandSurface.reserve(dots.size());
	for (auto p : dots)
		ligandSurface.emplace_back(p, std::numeric_limits<float>::max(), std::numeric_limits<float>::min());

	for (const auto &[sp_c, radius] : atomSpheres)
	{
		for (auto &[line, l1, l2] : ligandSurface)
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

	for (auto &[line, l1, l2] : ligandSurface)
	{
		auto p1 = /* center + */ l1 * line;
		auto p2 = /* center + */ l2 * line;

		// std::cout << "p1: " << p1 << ", p2: " << p2 << "\n";
		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f} }},\n", p1.get_x(), p1.get_y(), p1.get_z())
				  << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f} }},\n", p2.get_x(), p2.get_y(), p2.get_z());
	}

	if (ligandSurface.size() != blobSurface.size())
		throw std::runtime_error("Internal error fitting shape");

	const int N = 2 * ligandSurface.size();

	struct Score
	{
		cif::quaternion q;
		double v;

		bool operator<(const Score &rhs) const
		{
			return v < rhs.v;
		}
	};

	std::vector<Score> best;

	for (size_t i = 0; i < dots.size(); ++i)
	{
		std::vector<cif::point> blobDots(N), ligandDots(N);

		for (auto li = ligandDots.begin(); auto &[line, l1, l2] : ligandSurface)
		{
			*li++ = line * l1;
			*li++ = line * l2;
		}

		for (auto li = blobDots.begin() + i; auto &[line, l1, l2] : blobSurface)
		{
			if (li == blobDots.end())
				li = blobDots.begin();

			*li++ = line * l1;
			*li++ = line * l2;
		}

		auto q = cif::align_points(ligandDots, blobDots);

		for (auto &p : ligandDots)
			p.rotate(q);

		auto s = cif::RMSd(ligandDots, blobDots);

		best.emplace_back(q, s);
		std::push_heap(best.begin(), best.end());
	}

	std::sort_heap(best.begin(), best.end());

	for (auto [q, v] : best)
	{
		std::cout << "q: " << q << ", v: " << v << "\n";
	}

	return 0;
}

} // namespace pdb_redo