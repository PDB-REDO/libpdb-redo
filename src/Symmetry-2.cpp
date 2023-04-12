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

#include "cif++.hpp"

#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/Symmetry-2.hpp"

namespace c = cif;

namespace pdb_redo
{

int32_t GetRotationalIndexNumber(int spacegroup, const clipper::RTop_frac &rt);

// --------------------------------------------------------------------

sym_op::sym_op(std::string_view s)
{
	auto b = s.data();
	auto e = b + s.length();

	int rnri;

	auto r = std::from_chars(b, e, rnri);
	
	rnr = rnri;
	t[0] = r.ptr[1] - '0';
	t[1] = r.ptr[2] - '0';
	t[2] = r.ptr[3] - '0';

	if (r.ec != std::errc() or rnri > 192 or r.ptr[0] != '_' or t[0] > 9 or t[1] > 9 or t[2] > 9)
		throw std::invalid_argument("Could not convert string into sym_op");
}

sym_op::sym_op(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt)
	: sym_op(describeRToperation(spacegroup, cell, rt))
{
}

std::string sym_op::string() const
{
	std::ostringstream os;
	os << (int)rnr << '_' << (int)t[0] << (int)t[1] << (int)t[2];
	return os.str();
}

clipper::RTop_frac sym_op::toClipperFrac(const clipper::Spacegroup &spacegroup) const
{
	auto spacegroup_nr = getSpacegroupNumber(spacegroup);

	const size_t N = cif::kSymopNrTableSize;
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;
		if (cif::kSymopNrTable[i].spacegroup() < spacegroup_nr)
			L = i + 1;
		else
			R = i - 1;
	}

	for (size_t i = L; i < N and cif::kSymopNrTable[i].spacegroup() == spacegroup_nr; ++i)
	{
		if (cif::kSymopNrTable[i].rotational_number() != rnr)
			continue;
		
		auto d = cif::kSymopNrTable[i].symop().data();

		clipper::Mat33<> rot(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);

		clipper::Vec3<> trn(
			(d[9] == 0 ? 0 : 1.0 * d[9] / d[10]),
			(d[11] == 0 ? 0 : 1.0 * d[11] / d[12]),
			(d[13] == 0 ? 0 : 1.0 * d[13] / d[14])
		);

		return clipper::RTop_frac(rot, trn + clipper::Vec3<clipper::ftype>{
			t[0] - 5.0f,
			t[1] - 5.0f,
			t[2] - 5.0f
		});
	}

	throw std::invalid_argument("symmetry operator not found!");
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

std::tuple<float,sym_op> closestSymmetryCopy(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	cif::point a, cif::point b)
{
	// move a as close as possible to the origin, and move b in the same way

	auto calculateD = [&](float c, float a)
	{
		float d = 0;
		assert(a != 0);
		if (a != 0)
		{
			while (c + d < -(a / 2))
				d += a;
			while (c + d > (a / 2))
				d -= a;
		}
		return d;
	};

	if (cell.a() == 0 or cell.b() == 0 or cell.c() == 0)
		throw std::runtime_error("Invalid cell, contains a dimension that is zero");

	cif::point d;

	d.m_x = calculateD(a.m_x, static_cast<float>(cell.a()));
	d.m_y = calculateD(a.m_y, static_cast<float>(cell.b()));
	d.m_z = calculateD(a.m_z, static_cast<float>(cell.c()));

	a += d;
	b += d;

	d.m_x = calculateD(b.m_x, static_cast<float>(cell.a()));
	d.m_y = calculateD(b.m_y, static_cast<float>(cell.b()));
	d.m_z = calculateD(b.m_z, static_cast<float>(cell.c()));

	auto ca = toClipper(a);
	auto cb = toClipper(b);

	auto cfa = ca.coord_frac(cell);
	auto cfb = cb.coord_frac(cell);

	float result_d = std::numeric_limits<float>::max();
	sym_op result_s;

	for (int i = 0; i < spacegroup.num_symops(); ++i)
	{
		auto rt = spacegroup.symop(i);
		sym_op s(i + 1);	//GetRotationalIndexNumber(getSpacegroupNumber(spacegroup), rt)

		auto scfb = cfb.transform(rt);

		for (int j = 0; j < 3; ++j)
		{
			while (scfb[j] - 0.5f > cfa[j])
			{
				s.t[j] -= 1;
				scfb[j] -= 1;
			}

			while (scfb[j] + 0.5f < cfa[j])
			{
				s.t[j] += 1;
				scfb[j] += 1;
			}
		}

		auto dsq = cell.metric_real().lengthsq(cfa - scfb);
		if (result_d > dsq)
		{
			result_d = dsq;
			result_s = s;
		}
	}

	result_s.rnr = GetRotationalIndexNumber(getSpacegroupNumber(spacegroup), spacegroup.symop(result_s.rnr - 1));

	return { std::sqrt(result_d), result_s };
}

std::tuple<float,cif::mm::atom> closestSymmetryCopy(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	const cif::mm::atom &a, const cif::mm::atom &b)
{
	auto &&[d, symop] = closestSymmetryCopy(spacegroup, cell, a.get_location(), b.get_location());
	return { d, symmetryCopy(b, spacegroup, cell, symop) };
}

// --------------------------------------------------------------------
// Unfortunately, clipper has a different numbering scheme than PDB
// for rotation numbers. So we created a table to map those.
// Perhaps a bit over the top, but hey....

int32_t GetRotationalIndexNumber(int spacegroup, const clipper::RTop_frac &rt)
{
	auto k = GetSymOpDataForRTop_frac(rt);

	const size_t N = cif::kSymopNrTableSize;
	int32_t L = 0, R = static_cast<int32_t>(N - 1);
	while (L <= R)
	{
		int32_t i = (L + R) / 2;
		if (cif::kSymopNrTable[i].spacegroup() < spacegroup)
			L = i + 1;
		else
			R = i - 1;
	}

	for (size_t i = L; i < N and cif::kSymopNrTable[i].spacegroup() == spacegroup; ++i)
	{
		if (cif::kSymopNrTable[i].symop() == k)
			return cif::kSymopNrTable[i].rotational_number();
	}

	throw std::runtime_error("Symmetry operation was not found in table, cannot find rotational number");
}

// -----------------------------------------------------------------------

std::string SpacegroupToHall(std::string spacegroup)
{
	int nr = cif::get_space_group_number(spacegroup);

	// yeah, sucks, I know, might be looping three times this way
	std::string result;
	for (size_t i = 0; i < cif::kNrOfSpaceGroups; ++i)
	{
		auto &sp = cif::kSpaceGroups[i];
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
	for (size_t i = 0; i < cif::kNrOfSpaceGroups; ++i)
	{
		auto &sg = cif::kSpaceGroups[i];
		if (sg.nr == nr)
			return clipper::Spgr_descr(sg.Hall, clipper::Spgr_descr::Hall);
	}

	throw std::runtime_error("Invalid spacegroup number: " + std::to_string(nr));
}

// --------------------------------------------------------------------

std::string describeRToperation(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt)
{
	if (not(rt.is_null() or rt.equals(clipper::RTop_orth::identity(), 0.0001f)))
	{
		auto spacegroup_nr = getSpacegroupNumber(spacegroup);

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

							std::ostringstream os;
							os << rnr << '_' << (u + 5) << (v + 5) << (w + 5);
							return os.str();
						}
					}
		}
	}

	return "1_555";
}

// --------------------------------------------------------------------

cif::point offsetToOrigin(const clipper::Cell &cell, const cif::point &p)
{
	cif::point d{};

	while (p.m_x + d.m_x < (cell.a() / 2))
		d.m_x += cell.a();
	while (p.m_x + d.m_x > (cell.a() / 2))
		d.m_x -= cell.a();

	while (p.m_y + d.m_y < (cell.b() / 2))
		d.m_y += cell.b();
	while (p.m_y + d.m_y > (cell.b() / 2))
		d.m_y -= cell.b();

	while (p.m_z + d.m_z < (cell.c() / 2))
		d.m_z += cell.c();
	while (p.m_z + d.m_z > (cell.c() / 2))
		d.m_z -= cell.c();

	return d;
};

std::tuple<int,int,int> offsetToOriginInt(const clipper::Cell &cell, const cif::point &p)
{
	auto o = offsetToOrigin(cell, p);
	return {
		std::rintf(o.m_x / cell.a()),
		std::rintf(o.m_y / cell.b()),
		std::rintf(o.m_z / cell.c())
	};
}

cif::mm::atom symmetryCopy(const cif::mm::atom &atom, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt)
{
	auto loc = atom.get_location();
	auto d = offsetToOrigin(cell, loc);

	loc += d;
	loc = toPoint(toClipper(loc).transform(rt));
	loc -= d;

	std::string rt_operation = describeRToperation(spacegroup, cell, rt);

	return cif::mm::atom(atom, loc, rt_operation);
}

cif::mm::atom symmetryCopy(const cif::mm::atom &atom, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, sym_op symop)
{
	auto loc = symmetryCopy(atom.get_location(), spacegroup, cell, symop.toClipperOrth(spacegroup, cell));
	return cif::mm::atom(atom, loc, symop.string());
}

cif::point symmetryCopy(const cif::point &loc, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt)
{
	auto o = offsetToOrigin(cell, loc);

	if (o.m_x or o.m_y or o.m_z)
	{
		clipper::RTop_orth rt_o(clipper::Mat33<>::identity(), toClipper(o));

		auto c1 = toClipper(loc);
		auto c2 = c1.transform(rt_o);
		auto c3 = c2.transform(rt);
		auto c4 = c3.transform(rt_o.inverse());
		return c4;

		// return toClipper(loc).transform(clipper::RTop_orth(rt_o * rt * rt_o.inverse()));
	}
	else
		return toClipper(loc).transform(rt);
}

cif::point symmetryCopy(const cif::point &loc, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, sym_op symop)
{
	return symmetryCopy(loc, spacegroup, cell, symop.toClipperOrth(spacegroup, cell));
}

} // namespace pdb_redo