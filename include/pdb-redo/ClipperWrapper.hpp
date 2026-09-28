// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cif++/cif++.hpp>
#include <clipper/core/coords.h>

namespace pdb_redo
{

clipper::Atom toClipper(const cif::mm::atom &atom);
clipper::Atom toClipper(cif::const_row_handle atom, cif::const_row_handle aniso_row);

// --------------------------------------------------------------------

clipper::Spacegroup getSpacegroup(const cif::datablock &db);
clipper::Cell getCell(const cif::datablock &db);

// --------------------------------------------------------------------

cif::symop_data GetSymOpDataForRTop_frac(const clipper::RTop_frac &rt);
int getSpacegroupNumber(const clipper::Spacegroup &sg);

} // namespace pdb_redo