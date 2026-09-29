// SPDX-FileCopyrightText: Maarten L. Hekkelman, 2026
// SPDX-License-Identifier: BSD-2-Clause

#include "pdb-redo/DensityCalculations.hpp"

#include "pdb-redo/ShapeFitter.hpp"

// --------------------------------------------------------------------

namespace pdb_redo
{

HullDensity calculateDensityAroundLigand(const cif::mm::residue &lig, const clipper::Xmap<float> &xmap, float r)
{
	const auto dots = pdb_redo::create_spherical_dots(15);

	std::vector<cif::point> pts;
	pts.reserve(lig.atoms().size() * dots.size());

	std::vector<std::tuple<cif::point, float>> locationsAndRadii;
	for (auto aa : lig.atoms())
		locationsAndRadii.emplace_back(aa.get_location(), cif::atom_type_traits(aa.get_type()).radius() + r);

	for (const auto &[loc, r] : locationsAndRadii)
	{
		for (auto dot : dots)
			pts.emplace_back(loc + r * dot);
	}

	for (const auto [loc, r] : locationsAndRadii)
	{
		std::erase_if(pts, [loc, r](const cif::point &p)
			{ return distance(p, loc) < r; });
	};

	float minD = 100, maxD = -100, sumD = 0;

	for (auto pt : pts)
	{
		clipper::Coord_orth cp{ pt.m_x, pt.m_y, pt.m_z };
		clipper::Coord_frac pf = cp.coord_frac(xmap.cell());
		auto dp = xmap.interp<clipper::Interp_cubic>(pf);

		if (minD > dp)
			minD = dp;
		if (maxD < dp)
			maxD = dp;
		sumD += dp;
	}

	float avgD = sumD / pts.size();
	float sumSqD = 0;
	for (auto pt : pts)
	{
		clipper::Coord_orth cp{ pt.m_x, pt.m_y, pt.m_z };
		clipper::Coord_frac pf = cp.coord_frac(xmap.cell());
		auto dp = xmap.interp<clipper::Interp_cubic>(pf);

		sumSqD += (dp - avgD) * (dp - avgD);
	}

	auto variance = sumSqD / pts.size();
	auto sd = std::sqrt(variance);

	return { minD, maxD, avgD, sd };
}

} // namespace pdb_redo