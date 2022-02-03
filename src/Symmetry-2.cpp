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

#include <atomic>
#include <mutex>

#include "cif++/CifUtils.hpp"
#include "cif++/Symmetry.hpp"

#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/Symmetry-2.hpp"

namespace c = cif;
namespace m = mmcif;

namespace pdb_redo
{

// --------------------------------------------------------------------

clipper::Coord_orth CalculateOffsetForCell(const m::Structure &p, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell)
{
	auto &atoms = p.atoms();
	size_t dim = atoms.size();

	std::vector<clipper::Coord_orth> locations;
	locations.reserve(dim);

	// bounding box
	m::Point pMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
		pMax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());

	for (auto &atom : atoms)
	{
		auto pt = atom.location();
		locations.push_back(toClipper(pt));

		if (pMin.mX > pt.mX)
			pMin.mX = pt.mX;
		if (pMin.mY > pt.mY)
			pMin.mY = pt.mY;
		if (pMin.mZ > pt.mZ)
			pMin.mZ = pt.mZ;

		if (pMax.mX < pt.mX)
			pMax.mX = pt.mX;
		if (pMax.mY < pt.mY)
			pMax.mY = pt.mY;
		if (pMax.mZ < pt.mZ)
			pMax.mZ = pt.mZ;
	};

	// correct locations so that the median of x, y and z are inside the cell
	std::vector<float> c(dim);
	auto median = [&]()
	{
		return dim % 1 == 0
		           ? c[dim / 2]
		           : (c[dim / 2 - 1] + c[dim / 2]) / 2;
	};

	std::transform(locations.begin(), locations.end(), c.begin(), [](auto &l)
		{ return static_cast<float>(l[0]); });
	std::sort(c.begin(), c.end());
	float mx = median();

	std::transform(locations.begin(), locations.end(), c.begin(), [](auto &l)
		{ return static_cast<float>(l[1]); });
	std::sort(c.begin(), c.end());
	float my = median();

	std::transform(locations.begin(), locations.end(), c.begin(), [](auto &l)
		{ return static_cast<float>(l[2]); });
	std::sort(c.begin(), c.end());
	float mz = median();

	if (cif::VERBOSE > 1)
		std::cerr << "median position of atoms: " << m::Point(mx, my, mz) << std::endl;

	auto calculateD = [&](float m, float c)
	{
		float d = 0;
		assert(c != 0);
		if (c != 0)
		{
			while (m + d < -(c / 2))
				d += c;
			while (m + d > (c / 2))
				d -= c;
		}
		return d;
	};

	if (cell.a() == 0 or cell.b() == 0 or cell.c() == 0)
		throw std::runtime_error("Invalid cell, contains a dimension that is zero");

	m::Point D;

	D.mX = calculateD(mx, static_cast<float>(cell.a()));
	D.mY = calculateD(my, static_cast<float>(cell.b()));
	D.mZ = calculateD(mz, static_cast<float>(cell.c()));

	if (D.mX != 0 or D.mY != 0 or D.mZ != 0)
	{
		if (cif::VERBOSE)
			std::cerr << "moving coorinates by " << D.mX << ", " << D.mY << " and " << D.mZ << std::endl;
	}

	return toClipper(D);
}

// --------------------------------------------------------------------

std::vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup &spacegroup,
	const clipper::Cell &cell)
{
	std::vector<clipper::RTop_orth> result;

	// to make the operation at index 0 equal to identity
	result.push_back(clipper::RTop_orth::identity());

	for (int i = 0; i < spacegroup.num_symops(); ++i)
	{
		const auto &symop = spacegroup.symop(i);

		for (int u : {-1, 0, 1})
			for (int v : {-1, 0, 1})
				for (int w : {-1, 0, 1})
				{
					if (i == 0 and u == 0 and v == 0 and w == 0)
						continue;

					auto rtop = clipper::RTop_frac(
						symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w))
					                .rtop_orth(cell);

					result.push_back(std::move(rtop));
				}
	}

	return result;
}

// --------------------------------------------------------------------
// Unfortunately, clipper has a different numbering scheme than PDB
// for rotation numbers. So we created a table to map those.
// Perhaps a bit over the top, but hey....

int32_t GetRotationalIndexNumber(int spacegroup, const clipper::RTop_frac &rt)
{
	auto &rot = rt.rot();
	auto &trn = rt.trn();

	auto rte = [&rot](int i, int j)
	{ return static_cast<int8_t>(lrint(rot(i, j))); };

	std::array<int, 15> krt{
		rte(0, 0), rte(0, 1), rte(0, 2),
		rte(1, 0), rte(1, 1), rte(1, 2),
		rte(2, 0), rte(2, 1), rte(2, 2)};

	for (int i = 0; i < 3; ++i)
	{
		int n = lrint(trn[i] * 24);
		int d = 24;

		if (n == 0 or std::abs(n) == 24)
			continue; // is 0, 0 in our table

		for (int j = 5; j > 1; --j)
			if (n % j == 0 and d % j == 0)
			{
				n /= j;
				d /= j;
			}

		n = (n + d) % d;

		switch (i)
		{
			case 0:
				krt[9] = n;
				krt[10] = d;
				break;
			case 1:
				krt[11] = n;
				krt[12] = d;
				break;
			case 2:
				krt[13] = n;
				krt[14] = d;
				break;
		}
	}

	m::SymopData k(krt);

	const size_t N = m::kSymopNrTableSize;
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;
		if (m::kSymopNrTable[i].spacegroup() < spacegroup)
			L = i + 1;
		else
			R = i - 1;
	}

	for (size_t i = L; i < N and m::kSymopNrTable[i].spacegroup() == spacegroup; ++i)
	{
		if (m::kSymopNrTable[i].symop() == k)
			return m::kSymopNrTable[i].rotational_number();
	}

	throw std::runtime_error("Symmetry operation was not found in table, cannot find rotational number");
}

// -----------------------------------------------------------------------

std::string SpacegroupToHall(std::string spacegroup)
{
	int nr = m::GetSpacegroupNumber(spacegroup);

	// yeah, sucks, I know, might be looping three times this way
	std::string result;
	for (size_t i = 0; i < m::kNrOfSpaceGroups; ++i)
	{
		auto &sp = m::kSpaceGroups[i];
		if (sp.nr == nr)
		{
			result = sp.Hall;
			break;
		}
	}

	if (result.empty())
		throw std::runtime_error("Spacegroup name " + spacegroup + " was not found in table");

	return result;
}

clipper::Spgr_descr GetCCP4SpacegroupDescr(int nr)
{
	for (size_t i = 0; i < m::kNrOfSpaceGroups; ++i)
	{
		auto &sg = m::kSpaceGroups[i];
		if (sg.nr == nr)
			return clipper::Spgr_descr(sg.Hall, clipper::Spgr_descr::Hall);
	}

	throw std::runtime_error("Invalid spacegroup number: " + std::to_string(nr));
}

// --------------------------------------------------------------------

std::string describeRToperation(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt)
{
	auto spacegroup_nr = mmcif::GetSpacegroupNumber(spacegroup.symbol_xhm(), mmcif::SpacegroupName::xHM);

	if (not(rt.is_null() or rt.equals(clipper::RTop_orth::identity(), 0.0001f)))
	{
		for (int i = 0; i < spacegroup.num_symops(); ++i)
		{
			const auto &symop = spacegroup.symop(i);

			for (int u : {-1, 0, 1})
				for (int v : {-1, 0, 1})
					for (int w : {-1, 0, 1})
					{
						// if (i == 0 and u == 0 and v == 0 and w == 0)
						// 	continue;

						auto rtop = clipper::RTop_frac(
							symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w))
						                .rtop_orth(cell);

						if (rtop.rot().equals(rt.rot(), 0.00001) and rtop.trn().equals(rt.trn(), 0.000001))
						{
							// gotcha

							auto rtop_f = rtop.rtop_frac(cell);

							int rnr = GetRotationalIndexNumber(spacegroup_nr, rtop_f);

							uint32_t t[3] = {
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[0]))),
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[1]))),
								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[2])))};

							if (t[0] > 9 or t[1] > 9 or t[2] > 9)
								throw std::runtime_error("Symmetry operation has an out-of-range translation.");

							return std::to_string(rnr) + "_" + std::to_string(t[0]) + std::to_string(t[1]) + std::to_string(t[2]);
						}
					}
		}
	}

	return "1_555";
}

// --------------------------------------------------------------------

mmcif::Atom symmetryCopy(const mmcif::Atom &atom, const mmcif::Point &d,
	const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt)
{
	auto loc = atom.location();
	loc += d;
	loc = toPoint(toClipper(loc).transform(rt));
	loc -= d;

	std::string rt_operation = describeRToperation(spacegroup, cell, rt);

	return mmcif::Atom(atom, loc, rt_operation);
}

// -----------------------------------------------------------------------

SymmetryAtomIteratorFactory::SymmetryAtomIteratorFactory(const m::Structure &p, int spacegroupNr, const clipper::Cell &cell)
	: mSpacegroup(GetCCP4SpacegroupDescr(spacegroupNr))
	, mD(toPoint(CalculateOffsetForCell(p, mSpacegroup, cell)))
	, mRtOrth(AlternativeSites(mSpacegroup, cell))
	, mCell(cell)
{
}

// string SymmetryAtomIteratorFactory::symop_mmcif(const m::Atom& a) const
// {
// 	string result;

// 	if (not a.isSymmetryCopy())
// 		result = "1_555";
// 	else
// 	{
// 		auto rtop_o = a.symop();

// 		for (int i = 0; i < mSpacegroup.num_symops(); ++i)
// 		{
// 			const auto& symop = mSpacegroup.symop(i);

// 			for (int u: { -1, 0, 1})
// 				for (int v: { -1, 0, 1})
// 					for (int w: { -1, 0, 1})
// 					{
// 						// if (i == 0 and u == 0 and v == 0 and w == 0)
// 						// 	continue;

// 						auto rtop = clipper::RTop_frac(
// 								symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w)
// 							).rtop_orth(mCell);

// 						if (rtop.rot().equals(rtop_o.rot(), 0.00001) and rtop.trn().equals(rtop_o.trn(), 0.000001))
// 						{
// 							// gotcha

// 							auto rtop_f = rtop.rtop_frac(mCell);

// 							int rnr = GetRotationalIndexNumber(mSpacegroupNr, rtop_f);

// 							uint32_t t[3] = {
// 								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[0]))),
// 								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[1]))),
// 								static_cast<uint32_t>(5 + static_cast<int>(rint(rtop_f.trn()[2])))
// 							};

// 							if (t[0] > 9 or t[1] > 9 or t[2] > 9)
// 								throw runtime_error("Symmetry operation has an out-of-range translation.");

// 							result += to_string(rnr) + "_"
// 								   + to_string(t[0])
// 								   + to_string(t[1])
// 								   + to_string(t[2]);

// 							return result;
// 						}
// 					}
// 		}
// 	}

// 	return result;
// }

} // namespace pdb_redo