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

// AtomShape, analogue to the similarly named code in clipper

#pragma once

#include <cif++.hpp>

namespace pdb_redo
{

// --------------------------------------------------------------------
// Class used in calculating radii

class AtomShape
{
  public:
	AtomShape(const cif::datablock &db, std::string_view atom_id, float resHigh, float resLow,
		bool electronScattering)
		: AtomShape(db["atom_site"].find1(cif::key("id") == atom_id), db["atom_site_anisotrop"].find1(cif::key("id") == atom_id),
			  resHigh, resLow, electronScattering)
	{
	}

	AtomShape(const cif::datablock &db, std::string_view atom_id, float resHigh, float resLow,
		bool electronScattering, float bFactor)
		: AtomShape(db["atom_site"].find1(cif::key("id") == atom_id), db["atom_site_anisotrop"].find1(cif::key("id") == atom_id),
			  resHigh, resLow, electronScattering, bFactor)
	{
	}

	AtomShape(cif::row_handle atom, cif::row_handle atom_aniso, float resHigh, float resLow,
		bool electronScattering, std::optional<float> bFactor = {});

	AtomShape(const cif::mm::atom &atom, float resHigh, float resLow, bool electronScattering, std::optional<float> bFactor = {})
		: AtomShape(atom.get_row(), atom.get_row_aniso(), resHigh, resLow, electronScattering, bFactor)
	{
	}

	~AtomShape();

	AtomShape(const AtomShape &) = delete;
	AtomShape &operator=(const AtomShape &) = delete;

	float radius() const;
	float calculatedDensity(float r) const;
	float calculatedDensity(cif::point p) const;

  private:
	struct AtomShapeImpl *mImpl;
};

} // namespace pdb_redo
