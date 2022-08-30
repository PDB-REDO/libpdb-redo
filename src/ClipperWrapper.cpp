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

#include "pdb-redo/ClipperWrapper.hpp"
#include <pdbx++/Symmetry.hpp>
// #include "cif++/Symmetry.hpp"

namespace pdb_redo
{

// --------------------------------------------------------------------

clipper::Atom toClipper(cif::row_handle atom, cif::row_handle aniso_row)
{
	const double kPI = pdbx::kPI;

	clipper::Atom result;

	pdbx::Point location = atom.get<float, float, float>("Cartn_x", "Cartn_y", "Cartn_z");
	result.set_coord_orth({ location.mX, location.mY, location.mZ });

	if (atom["occupancy"].empty())
		result.set_occupancy(1.0);
	else
		result.set_occupancy(atom["occupancy"].as<float>());

	std::string element = atom["type_symbol"].as<std::string>();
	if (not atom["pdbx_formal_charge"].empty())
	{
		int charge = atom["pdbx_formal_charge"].as<int>();
		if (std::abs(charge) > 1)
			element += std::to_string(charge);
		if (charge < 0)
			element += '-';
		else
			element += '+';
	}
	result.set_element(element);

	if (not atom["U_iso_or_equiv"].empty())
		result.set_u_iso(atom["U_iso_or_equiv"].as<float>());
	else if (not atom["B_iso_or_equiv"].empty())
		result.set_u_iso(atom["B_iso_or_equiv"].as<float>() / (8 * kPI * kPI));
	else
		throw std::runtime_error("Missing B_iso or U_iso");

	if (aniso_row.empty())
		result.set_u_aniso_orth(clipper::U_aniso_orth(nan("0"), 0, 0, 0, 0, 0));
	else
	{
		const auto &[u11, u12, u13, u22, u23, u33] =
			aniso_row.get<float, float, float, float, float, float>("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
		result.set_u_aniso_orth(clipper::U_aniso_orth(u11, u22, u33, u12, u13, u23));
	}

	return result;
}

// --------------------------------------------------------------------

clipper::Spacegroup getSpacegroup(const cif::datablock &db)
{
	std::string spacegroup = db["symmetry"].find1<std::string>(cif::key("entry_id") == db.name(), "space_group_name_H-M");

	if (spacegroup == "P 1-")
		spacegroup = "P -1";
	else if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw std::runtime_error("No spacegroup, cannot continue");

	try
	{
		return clipper::Spacegroup{ clipper::Spgr_descr(pdbx::GetSpacegroupNumber(spacegroup)) };
	}
	catch (const clipper::Message_fatal &m)
	{
		// std::cout << m.text() << std::endl;
	}

	try
	{
		return clipper::Spacegroup{ clipper::Spgr_descr(spacegroup) };
	}
	catch (const clipper::Message_fatal &e)
	{
		std::cerr << e.text() << std::endl;
	}

	throw std::runtime_error("Unsupported spacegroup: " + spacegroup);
}

clipper::Cell getCell(const cif::datablock &db)
{
	const auto &[a, b, c, alpha, beta, gamma] =
		db["cell"].find1<float, float, float, float, float, float>(cif::key("entry_id") == db.name(),
			"length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma");

	return clipper::Cell{ clipper::Cell_descr(a, b, c, alpha, beta, gamma) };
}

} // namespace pdb_redo