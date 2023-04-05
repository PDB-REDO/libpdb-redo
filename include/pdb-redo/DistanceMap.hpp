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

#include <unordered_map>

#include <clipper/clipper.h>

#include <cif++.hpp>

#include <pdb-redo/ClipperWrapper.hpp>
#include <pdb-redo/Symmetry-2.hpp>

#ifdef near
#undef near
#endif

namespace pdb_redo
{

class DistanceMap
{
  public:

	DistanceMap(const cif::mm::structure &p, const clipper::Spacegroup &spacegroup, const clipper::Cell &cell,
		float maxDistance);

	DistanceMap(const cif::mm::structure &p, float maxDistance)
		: DistanceMap(p, getSpacegroup(p.get_datablock()), getCell(p.get_datablock()), maxDistance)
	{
	}

	DistanceMap(const DistanceMap &) = delete;
	DistanceMap &operator=(const DistanceMap &) = delete;

	float operator()(const std::string &a, const std::string &b) const;

	std::vector<cif::mm::atom> near(const cif::mm::atom &atom, float maxDistance = 3.5f) const;

  private:
	using DistKeyType = std::tuple<size_t, size_t>;
	using DistValueType = std::tuple<float, sym_op>;
	using DistMap = std::map<DistKeyType, DistValueType>;

	void AddDistancesForAtoms(const std::vector<std::tuple<size_t,cif::point>> &a,
		const std::vector<std::tuple<size_t,cif::point>> &b, DistMap &dm, sym_op rtop = {});

	cif::point offsetToOrigin(const cif::point &p) const;

	const cif::mm::structure &mStructure;
	clipper::Cell cell;
	clipper::Spacegroup spacegroup;
	size_t dim;
	std::unordered_map<std::string, size_t> index;
	std::map<size_t, std::string> rIndex;

	float mMaxDistance, mMaxDistanceSQ;

	std::vector<std::tuple<float, sym_op>> mA;
	std::vector<size_t> mIA, mJA;
};

} // namespace pdb_redo
