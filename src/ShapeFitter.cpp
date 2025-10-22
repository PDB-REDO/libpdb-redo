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

#include <cmath>
#include <limits>

namespace pdb_redo
{

double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap)
{
	std::vector<cif::point> blob;
	cif::point center;

	for (auto i = xmap.first(); not i.last(); i.next())
	{
		if (xmap[i] > 0)
		{
			cif::point p{ i.coord_orth() };
			center += p;
			blob.emplace_back(p);
		}
	}

	center /= blob.size();

	auto cellVolume = xmap.cell().volume();
	auto gridSize = xmap.grid_sampling().size();
	auto gridPointVolume = cellVolume / gridSize;
	auto gridPointRadius = std::pow((3 * gridPointVolume) / (4 * cif::kPI), 1 / 3.0);

	using list_of_spheres = std::vector<std::tuple<cif::point, float>>;

	list_of_spheres mapSpheres;
	mapSpheres.reserve(blob.size());
	for (auto p : blob)
		mapSpheres.emplace_back(p - center, gridPointRadius);

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