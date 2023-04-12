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

#include "pdb-redo/ClipperWrapper.hpp"

namespace pdb_redo
{


namespace literals
{

	inline sym_op operator""_so(const char *text, size_t length)
	{
		return sym_op({ text, length });
	}

}


// --------------------------------------------------------------------
// Functions to use when working with symmetry stuff

std::vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell);
// int GetSpacegroupNumber(std::string spacegroup);	// alternative for clipper's parsing code
// std::string SpacegroupToHall(std::string spacegroup);

// cif::mm::atom symmetryCopy(const cif::mm::atom &atom, const cif::point &d,
// 	const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);

cif::mm::atom symmetryCopy(const cif::mm::atom &atom, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);
cif::mm::atom symmetryCopy(const cif::mm::atom &atom, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, sym_op symop);
cif::point symmetryCopy(const cif::point &loc, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);
cif::point symmetryCopy(const cif::point &loc, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, sym_op symop);

std::string describeRToperation(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);

/// Return the closest RTop and distance. The rtop should be applied to \a b to get the actual point nearest to \a a.
std::tuple<float,sym_op> closestSymmetryCopy(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	cif::point a, cif::point b);

std::tuple<float,cif::mm::atom> closestSymmetryCopy(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	const cif::mm::atom &a, const cif::mm::atom &b);

std::tuple<float,cif::mm::atom> closestSymmetryCopy(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	const cif::point &pt, const cif::mm::atom &b);


} // namespace pdb_redo