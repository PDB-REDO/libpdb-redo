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

struct sym_op
{
	uint8_t rnr;
	uint8_t t[3];

	sym_op(uint8_t rnri = 1, uint8_t tx = 5, uint8_t ty = 5, uint8_t tz = 5)
	{
		rnr = rnri;
		t[0] = tx;
		t[1] = ty;
		t[2] = tz;
	}

	sym_op(std::string_view s);

	sym_op(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);

	sym_op(const sym_op &) = default;
	sym_op(sym_op &&) = default;
	sym_op &operator=(const sym_op &) = default;
	sym_op &operator=(sym_op &&) = default;

	constexpr bool is_identity() const
	{
		return rnr == 1 and t[0] == 5 and t[1] == 5 and t[2] == 5;
	}

	explicit constexpr operator bool() const
	{
		return not is_identity();
	}

	std::string string() const;

	clipper::RTop_frac toClipperFrac(const clipper::Spacegroup &spacegroup) const;

	clipper::RTop_orth toClipperOrth(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell) const
	{
		return toClipperFrac(spacegroup).rtop_orth(cell);
	}
};

static_assert(sizeof(sym_op) == 4, "Sym_op should be four bytes");

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
std::tuple<float,clipper::RTop_orth> closestSymmetryCopy(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
	cif::point a, cif::point b);

// --------------------------------------------------------------------
// To iterate over all symmetry copies of an atom

class SymmetryAtomIteratorFactory
{
  public:
	SymmetryAtomIteratorFactory(const cif::mm::structure &structure, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell);
	SymmetryAtomIteratorFactory(const cif::mm::structure &structure, int spacegroupNr, const clipper::Cell &cell);

	SymmetryAtomIteratorFactory(const cif::mm::structure &structure)
		: SymmetryAtomIteratorFactory(structure, getSpacegroup(structure.get_datablock()), getCell(structure.get_datablock()))
	{
	}

	SymmetryAtomIteratorFactory(const SymmetryAtomIteratorFactory &) = delete;
	SymmetryAtomIteratorFactory &operator=(const SymmetryAtomIteratorFactory &) = delete;

	template <typename ConditionFunc>
	class SymmetryAtomIterator
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = cif::mm::atom;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type*;
		using reference = value_type &;

		SymmetryAtomIterator(const SymmetryAtomIteratorFactory &factory, const cif::mm::atom &atom, ConditionFunc &cond)
			: m_f(&factory)
			, m_i(0)
			, m_a(atom)
			, m_c(atom)
			, m_l(m_a.get_location())
			, m_cond(cond)
		{
			// Find the offset first, that needs to be applied
			auto calc_offset = [&](float c, float e)
			{
				float d = 0;
				assert(e != 0);
				if (e != 0)
				{
					while (c + d < -(e / 2))
						d += e;
					while (c + d > (e / 2))
						d -= e;
				}
				return d;
			};

			m_o.m_x = calc_offset(m_l.m_x, m_f->mCell.a());
			m_o.m_y = calc_offset(m_l.m_y, m_f->mCell.b());
			m_o.m_z = calc_offset(m_l.m_z, m_f->mCell.c());

			while (not test() and m_i < m_f->mRtOrth.size())
				++m_i;
		}

		SymmetryAtomIterator(const SymmetryAtomIteratorFactory &factory, const cif::mm::atom &atom, ConditionFunc &cond, int)
			: SymmetryAtomIterator(factory, atom, cond)
		{
			m_i = m_f->mRtOrth.size();
		}

		SymmetryAtomIterator(const SymmetryAtomIterator &iter)
			: m_f(iter.m_f)
			, m_i(iter.m_i)
			, m_a(iter.m_a)
			, m_c(iter.m_c)
			, m_cond(iter.m_cond)
		{
		}

		SymmetryAtomIterator &operator=(const SymmetryAtomIterator &iter)
		{
			if (this != &iter)
			{
				m_f = iter.m_f;
				m_i = iter.m_i;
				m_a = iter.m_a;
				m_c = iter.m_c;
			}
			return *this;
		}

		reference operator*() { return m_c; }
		pointer operator->() { return &m_c; }

		SymmetryAtomIterator operator++()
		{
			while (m_i < m_f->mRtOrth.size())
			{
				++m_i;
				if (test())
					break;
			}

			return *this;
		}

		SymmetryAtomIterator operator++(int)
		{
			SymmetryAtomIterator result(*this);
			this->operator++();
			return result;
		}

		bool operator==(const SymmetryAtomIterator &iter) const
		{
			return m_f == iter.m_f and m_i == iter.m_i;
		}

		bool operator!=(const SymmetryAtomIterator &iter) const
		{
			return m_f != iter.m_f or m_i != iter.m_i;
		}

	  private:
		bool test()
		{
			bool result = false;

			if (m_i < m_f->mRtOrth.size())
			{
				auto &rt = m_f->mRtOrth[m_i];
				auto loc = m_l;

				loc += m_o;
				loc = toPoint(toClipper(loc).transform(rt));
				loc -= m_o;

				if (m_cond(loc))
				{
					std::string rt_operation = describeRToperation(m_f->mSpacegroup, m_f->mCell, rt);
					m_c = cif::mm::atom(m_a, loc, rt_operation);
					result = true;
				}
			}

			return result;
		}

		const SymmetryAtomIteratorFactory *m_f;
		size_t m_i;
		cif::mm::atom m_a, m_c;
		cif::point m_l, m_o;
		ConditionFunc &m_cond;
	};

	template <typename ConditionFunc>
	class SymmetryAtomIteratorRange
	{
	  public:
		SymmetryAtomIteratorRange(const SymmetryAtomIteratorFactory &f, const cif::mm::atom &a, ConditionFunc &&cond)
			: m_f(f)
			, m_a(a)
			, m_cond(std::move(cond))
		{
		}

		SymmetryAtomIterator<ConditionFunc> begin()
		{
			return SymmetryAtomIterator(m_f, m_a, m_cond);
		}

		SymmetryAtomIterator<ConditionFunc> end()
		{
			return SymmetryAtomIterator(m_f, m_a, m_cond, 1);
		}

	  private:
		const SymmetryAtomIteratorFactory &m_f;
		cif::mm::atom m_a;
		ConditionFunc m_cond;
	};

	template <typename ConditionFunc>
	SymmetryAtomIteratorRange<ConditionFunc> operator()(const cif::mm::atom &a, ConditionFunc &&cond) const
	{
		return SymmetryAtomIteratorRange(*this, a, std::forward<ConditionFunc>(cond));
	}

	auto operator()(const cif::mm::atom &a, const cif::point &loc, float maxDistance) const
	{
		return this->operator()(a, [loc, dsq = maxDistance * maxDistance](const cif::point &p) { return distance_squared(p, loc) <= dsq; });
	}

	auto operator()(const cif::mm::atom &a) const
	{
		return this->operator()(a, [](const cif::point &p) { return true; });
	}

  private:
	clipper::Spacegroup mSpacegroup;
	std::vector<clipper::RTop_orth> mRtOrth;
	clipper::Cell mCell;
};

} // namespace pdb_redo