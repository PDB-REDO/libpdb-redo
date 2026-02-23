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

#include <cif++/cif++.hpp>

#include <unordered_map>

#ifdef near
#undef near
#endif

namespace pdb_redo
{

class DistanceMap
{
  public:
	DistanceMap(const cif::mm::structure &p, cif::crystal crystal, float maxDistance);

	DistanceMap(const cif::mm::structure &p, float maxDistance)
		: DistanceMap(p, cif::crystal(p.get_datablock()), maxDistance)
	{
	}

	DistanceMap(const DistanceMap &) = delete;
	DistanceMap &operator=(const DistanceMap &) = delete;

	float operator()(const std::string &a, const std::string &b) const;

	std::vector<cif::mm::atom> near(const cif::mm::atom &atom, float maxDistance = 3.5f) const;

  private:
	const cif::mm::structure &m_structure;
	cif::crystal m_crystal;
	float m_grid_spacing;

	struct key_type
	{
		int x, y, z;

		constexpr bool operator<=>(const key_type &) const noexcept = default;
	};

	struct key_type_hash
	{
		std::size_t operator()(const key_type &s) const noexcept
		{
			auto h0 = std::hash<int>{}(s.x);
			auto h1 = std::hash<int>{}(s.y);
			auto h2 = std::hash<int>{}(s.z);

			return h0 ^ (h1 << 1) ^ (h2 << 2);
		}
	};

	struct entry
	{
		std::string id;
		cif::sym_op symop;
	};

	std::unordered_multimap<key_type, entry, key_type_hash> m_index;
};

} // namespace pdb_redo
