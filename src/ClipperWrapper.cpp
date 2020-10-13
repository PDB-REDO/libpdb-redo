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

#include "pdb-redo.hpp"

#include "ClipperWrapper.hpp"

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
		if (abs(charge) > 1)
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

mmcif::Atom symmetryCopy(const mmcif::Atom& atom, const mmcif::Point& d, const clipper::RTop_orth& rt)
{
	auto loc = atom.location();
	loc += d;
	loc = ((clipper::Coord_orth)loc).transform(rt);

	return mmcif::Atom(atom, loc);
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
