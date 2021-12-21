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

#include "cif++/Structure.hpp"

#include "pdb-redo/ClipperWrapper.hpp"

namespace pdb_redo
{

// --------------------------------------------------------------------
// Functions to use when working with symmetry stuff

clipper::Coord_orth CalculateOffsetForCell(const mmcif::Structure &p, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell);
std::vector<clipper::RTop_orth> AlternativeSites(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell);
// int GetSpacegroupNumber(std::string spacegroup);	// alternative for clipper's parsing code
// std::string SpacegroupToHall(std::string spacegroup);

mmcif::Atom symmetryCopy(const mmcif::Atom &atom, const mmcif::Point &d,
	const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);

std::string describeRToperation(const clipper::Spacegroup &spacegroup, const clipper::Cell &cell, const clipper::RTop_orth &rt);

// --------------------------------------------------------------------
// To iterate over all symmetry copies of an atom

class SymmetryAtomIteratorFactory
{
  public:
	SymmetryAtomIteratorFactory(const mmcif::Structure &p, int spacegroupNr, const clipper::Cell &cell);

	SymmetryAtomIteratorFactory(const SymmetryAtomIteratorFactory &) = delete;
	SymmetryAtomIteratorFactory &operator=(const SymmetryAtomIteratorFactory &) = delete;

	template <typename ConditionFunc>
	class SymmetryAtomIterator : public std::iterator<std::forward_iterator_tag, const mmcif::Atom>
	{
	  public:
		typedef std::iterator<std::forward_iterator_tag, const mmcif::Atom> baseType;
		typedef typename baseType::pointer pointer;
		typedef typename baseType::reference reference;

		SymmetryAtomIterator(const SymmetryAtomIteratorFactory &factory, const mmcif::Atom &atom, ConditionFunc &cond)
			: m_f(&factory)
			, m_i(0)
			, m_a(atom)
			, m_c(atom)
			, m_cond(cond)
		{
			while (not test() and m_i < m_f->mRtOrth.size())
				++m_i;
		}

		SymmetryAtomIterator(const SymmetryAtomIteratorFactory &factory, const mmcif::Atom &atom, ConditionFunc &cond, int)
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
				auto loc = m_a.location();

				loc += m_f->mD;
				loc = toPoint(toClipper(loc).transform(rt));
				loc -= m_f->mD;

				if (m_cond(loc))
				{
					std::string rt_operation = describeRToperation(m_f->mSpacegroup, m_f->mCell, rt);
					m_c = mmcif::Atom(m_a, loc, rt_operation);
					result = true;
				}
			}

			return result;
		}

		const SymmetryAtomIteratorFactory *m_f;
		size_t m_i;
		mmcif::Atom m_a, m_c;
		ConditionFunc &m_cond;
	};

	template <typename ConditionFunc>
	class SymmetryAtomIteratorRange
	{
	  public:
		SymmetryAtomIteratorRange(const SymmetryAtomIteratorFactory &f, const mmcif::Atom &a, ConditionFunc &&cond)
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
		mmcif::Atom m_a;
		ConditionFunc m_cond;
	};

	template <typename ConditionFunc>
	SymmetryAtomIteratorRange<ConditionFunc> operator()(const mmcif::Atom &a, ConditionFunc &&cond) const
	{
		return SymmetryAtomIteratorRange(*this, a, std::forward<ConditionFunc>(cond));
	}

	// std::string symop_mmcif(const mmcif::Atom& a) const;

  private:
	clipper::Spacegroup mSpacegroup;
	mmcif::Point mD; // needed to move atoms to center
	std::vector<clipper::RTop_orth> mRtOrth;
	clipper::Cell mCell;
};

} // namespace pdb_redo