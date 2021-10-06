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

#include "cif++/Symmetry.hpp"
#include "pdb-redo/ClipperWrapper.hpp"

// --------------------------------------------------------------------

clipper::Atom toClipper(const mmcif::Atom& atom)
{
	using mmcif::kPI;

	auto row = atom.getRow();

	clipper::Atom result;

	auto location = atom.location();
	result.set_coord_orth({location.mX, location.mY, location.mZ});
	
	if (row["occupancy"].empty())
		result.set_occupancy(1.0);
	else
		result.set_occupancy(row["occupancy"].as<float>());
	
	std::string element = row["type_symbol"].as<std::string>();
	if (not row["pdbx_formal_charge"].empty())
	{
		int charge = row["pdbx_formal_charge"].as<int>();
		if (std::abs(charge) > 1)
			element += std::to_string(charge);
		if (charge < 0)
			element += '-';
		else
			element += '+';
	}
	result.set_element(element);
	
	if (not row["U_iso_or_equiv"].empty())
		result.set_u_iso(row["U_iso_or_equiv"].as<float>());
	else if (not row["B_iso_or_equiv"].empty())
		result.set_u_iso(row["B_iso_or_equiv"].as<float>() / (8 * kPI * kPI));
	else
		throw std::runtime_error("Missing B_iso or U_iso");	
	
	auto r = atom.getRowAniso();
	if (r.empty())
		result.set_u_aniso_orth(clipper::U_aniso_orth(nan("0"), 0, 0, 0, 0, 0));
	else
	{
		float u11, u12, u13, u22, u23, u33;
		cif::tie(u11, u12, u13, u22, u23, u33) =
			r.get("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
		
		result.set_u_aniso_orth(clipper::U_aniso_orth(u11, u22, u33, u12, u13, u23));
	}
	
	return result;
}

// --------------------------------------------------------------------

clipper::Spacegroup getSpacegroup(const mmcif::Structure& structure)
{
	auto &db = structure.getFile().data();

	auto refine = db["refine"][cif::Key("entry_id") == db.getName()];
	if (refine.empty())
		throw std::runtime_error("No refinement data found");
	
	double reso, hires = 99, lowres = 0;
	cif::tie(reso, hires, lowres) = refine.get("ls_d_res_high", "ls_d_res_high", "ls_d_res_low");
	
	std::string spacegroup = db["symmetry"]
		[cif::Key("entry_id") == db.getName()]
		["space_group_name_H-M"].as<std::string>();
	
	if (spacegroup == "P 1-")
		spacegroup = "P -1";
	else if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw std::runtime_error("No spacegroup, cannot continue");
	
	try
	{
		return clipper::Spacegroup{clipper::Spgr_descr(mmcif::GetSpacegroupNumber(spacegroup))};
	}
	catch (const clipper::Message_fatal& m)
	{
		// std::cout << m.text() << std::endl;
	}

	try
	{
		return clipper::Spacegroup{clipper::Spgr_descr(spacegroup)};
	}
	catch (const clipper::Message_fatal& e)
	{
		std::cerr << e.text() << std::endl;
	}
	
	throw std::runtime_error("Unsupported spacegroup: " + spacegroup);

	// // reconstruct a clipper spacegroup.

	// int sg = mmcif::GetSpacegroupNumber(spacegroup);

	// // locate the symmetry operations
	// int L = 0, R = mmcif::kSymopNrTableSize - 1;
	// while (R >= L)
	// {
	// 	int i = (L + R) / 2;
	// 	if (mmcif::kSymopNrTable[i].spacegroup() < sg)
	// 		L = i + 1;
	// 	else
	// 		R = i - 1;
	// }

	// if (mmcif::kSymopNrTable[L].spacegroup() != sg)
	// 	throw std::runtime_error("Could not locate spacegroup " + spacegroup);
	
	// clipper::Spgr_descr::Symop_codes codes;

	// for (int i = L; i < mmcif::kSymopNrTableSize and mmcif::kSymopNrTable[i].spacegroup() == sg; ++i)
	// {
	// 	auto data = mmcif::kSymopNrTable[i].symop().data();

	// 	clipper::Mat33<> m(
	// 			data[0], data[1], data[2],
	// 			data[3], data[4], data[5],
	// 			data[6], data[7], data[8]);

	// 	clipper::Vec3<> v(
	// 		1.0f * data[ 9] / data[10],
	// 		1.0f * data[11] / data[12],
	// 		1.0f * data[13] / data[14]);

	// 	clipper::RTop<> rtop(m, v);

	// 	codes.emplace_back(clipper::Symop{rtop});
	// }

	// return clipper::Spacegroup{clipper::Spgr_descr{codes}};
}

clipper::Cell getCell(const mmcif::Structure& structure)
{
	auto &db = structure.getFile().data();

	double a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = db["cell"][cif::Key("entry_id") == db.getName()]
		.get("length_a", "length_b", "length_c",
			 "angle_alpha", "angle_beta", "angle_gamma");

	return clipper::Cell{clipper::Cell_descr(a, b, c, alpha, beta, gamma)};
}


	// AtomImpl(const AtomImpl& impl, const Point& d, const clipper::RTop_orth& rt)
	// 	: mFile(impl.mFile), mID(impl.mID), mType(impl.mType), mAtomID(impl.mAtomID)
	// 	, mCompID(impl.mCompID), mAsymID(impl.mAsymID), mSeqID(impl.mSeqID)
	// 	, mAltID(impl.mAltID), mLocation(impl.mLocation), mRefcount(1)
	// 	, mRow(impl.mRow), mCompound(impl.mCompound), mRadius(impl.mRadius)
	// 	, mCachedProperties(impl.mCachedProperties)
	// 	, mSymmetryCopy(true), mRTop(rt), mD(d)
	// {
	// 	mLocation += d;
	// 	mLocation = ((clipper::Coord_orth)mLocation).transform(rt);
	// 	mLocation -= d;
	// }
