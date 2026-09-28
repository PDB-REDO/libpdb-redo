// SPDX-FileCopyrightText: Maarten L. Hekkelman, 2026
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cif++/model.hpp>
#include <clipper/clipper.h>
#include <tuple>

namespace pdb_redo
{

/// Result type for calculateDensityAroundLigand
struct HullDensity
{
	float minDensity, maxDensity, avgDensity, sd;
};

/// Calculate the average density for the hull around \a lig with a distance of \a r to the atoms
HullDensity calculateDensityAroundLigand(const cif::mm::residue &lig, const clipper::Xmap<float> &xmap, float r);

} // namespace pdb_redo
