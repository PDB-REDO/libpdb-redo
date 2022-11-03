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

#pragma once

#include <cif++.hpp>
#include <clipper/core/coords.h>

namespace pdb_redo
{

clipper::Atom toClipper(const cif::mm::atom &atom);
clipper::Atom toClipper(cif::row_handle atom, cif::row_handle aniso_row);

inline clipper::Coord_orth toClipper(const cif::point &pt)
{
	return { pt.m_x, pt.m_y, pt.m_z };
}

inline cif::point toPoint(const clipper::Coord_orth &pt)
{
	return { static_cast<float>(pt.x()), static_cast<float>(pt.y()), static_cast<float>(pt.z()) };
}

// --------------------------------------------------------------------

clipper::Spacegroup getSpacegroup(const cif::datablock &db);
clipper::Cell getCell(const cif::datablock &db);

// --------------------------------------------------------------------

cif::symop_data GetSymOpDataForRTop_frac(const clipper::RTop_frac &rt);
int getSpacegroupNumber(const clipper::Spacegroup &sg);

} // namespace pdb_redo