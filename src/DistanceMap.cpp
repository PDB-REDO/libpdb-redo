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

#include "pdb-redo/DistanceMap.hpp"

#include <cif++/utilities.hpp>

namespace pdb_redo
{

using cif::point;

// --------------------------------------------------------------------

std::tuple<point, float> calculateCenterAndRadius(const std::vector<std::tuple<std::size_t, point>> &atoms)
{
	std::vector<point> pts;
	for (const auto &[ix, pt] : atoms)
		pts.emplace_back(pt);

	auto center = centroid(pts);
	float radius = 0;

	for (auto &pt : pts)
	{
		auto d = static_cast<float>(distance(pt, center));
		if (radius < d)
			radius = d;
	}

	return std::make_tuple(center, radius);
}

// --------------------------------------------------------------------

DistanceMap::DistanceMap(std::vector<cif::mm::atom> atoms, cif::crystal crystal, float maxDistance)
	: mAtoms(std::move(atoms))
	, mCrystal(std::move(crystal))
{
	std::ranges::sort(mAtoms, [](auto &a, auto &b) { return a.id().compare(b.id()) < 0; });

	std::vector<std::tuple<cif::point, std::string>> pts;

	pts.reserve(mAtoms.size());

	KeyType k1{}, k2{};

	for (auto a : mAtoms)
	{
		pts.emplace_back(a.get_location(), a.id());

		auto p = a.get_location();
		KeyType k{
			static_cast<int16_t>(std::rint(p.m_x)),
			static_cast<int16_t>(std::rint(p.m_y)),
			static_cast<int16_t>(std::rint(p.m_z))
		};

		if (mIndex.empty())
		{
			k1.x = k2.x = k.x;
			k1.y = k2.y = k.y;
			k1.z = k2.z = k.z;
		}
		else
		{
			if (k1.x > k.x)
				k1.x = k.x;
			else if (k2.x < k.x)
				k2.x = k.x;

			if (k1.y > k.y)
				k1.y = k.y;
			else if (k2.y < k.y)
				k2.y = k.y;

			if (k1.z > k.z)
				k1.z = k.z;
			else if (k2.z < k.z)
				k2.z = k.z;
		}

		mIndex.emplace(k, Entry{ a.id(), cif::sym_op{} });
	}

	int d = static_cast<int>(std::rint(maxDistance));
	k1.x -= d;
	k2.x += d;
	k1.y -= d;
	k2.y += d;
	k1.z -= d;
	k2.z += d;

	auto &sg = mCrystal.get_spacegroup();
	auto &cell = mCrystal.get_cell();

	cif::progress_bar progress(sg.size() * 9 * 9 * 9 - 1, "Creating distancemap");

	for (uint8_t i = 1; std::cmp_less(i, sg.size() + 1); ++i)
	{
		for (uint8_t tx = 1; tx <= 9; ++tx)
		{
			for (uint8_t ty = 1; ty <= 9; ++ty)
			{
				for (uint8_t tz = 1; tz <= 9; ++tz)
				{
					cif::sym_op symop(i, tx, ty, tz);

					if (not symop) // skip the identity symop
						continue;

					for (auto &[pt, id] : pts)
					{
						auto ap = sg(pt, cell, symop);

						KeyType k{
							static_cast<int16_t>(std::rint(ap.m_x)),
							static_cast<int16_t>(std::rint(ap.m_y)),
							static_cast<int16_t>(std::rint(ap.m_z))
						};

						if (k.x >= k1.x and k.x <= k2.x and
							k.y >= k1.y and k.y <= k2.y and
							k.y >= k1.z and k.z <= k2.z)
						{
							mIndex.emplace(k, Entry{ id, symop });
						}
					}

					progress.consumed(1);
				}
			}
		}
	}
}

cif::mm::atom DistanceMap::getAtomByID(const std::string &id) const
{
	cif::mm::atom result;

	ssize_t L = 0z, R = mAtoms.size();
	while (L <= R)
	{
		auto i = (L + R) / 2;

		auto d = mAtoms[i].id().compare(id);
		if (d == 0)
		{
			result = mAtoms[i];
			break;
		}
		else if (d < 0)
			L = i + 1;
		else
			R = i - 1;
	}

	return result;
}

std::vector<cif::mm::atom> DistanceMap::near(const cif::mm::atom &atom, float maxDistance) const
{
	std::vector<cif::mm::atom> result;

	auto maxDistanceSq = maxDistance * maxDistance;

	auto p = atom.get_location();

	KeyType k{
		static_cast<int16_t>(std::rint(p.m_x)),
		static_cast<int16_t>(std::rint(p.m_y)),
		static_cast<int16_t>(std::rint(p.m_z))
	};

	KeyType k1 = k, k2 = k;

	int d = static_cast<int>(std::ceil(maxDistance));
	k1.x -= d;
	k1.y -= d;
	k1.z -= d;

	k2.x += d;
	k2.y += d;
	k2.z += d;

	for (k.x = k1.x; k.x <= k2.x; ++k.x)
	{
		for (k.y = k1.y; k.y <= k2.y; ++k.y)
		{
			for (k.z = k1.z; k.z <= k2.z; ++k.z)
			{
				auto r = mIndex.equal_range(k);
				for (auto &[k, e] : std::ranges::subrange(r.first, r.second))
				{
					auto [id, symop] = e;
					if (id == atom.id())
						continue;

					auto a = getAtomByID(id);
					auto loc = mCrystal.symmetry_copy(a.get_location(), symop);

					if (auto d = cif::distance_squared(p, loc); d <= maxDistanceSq)
					{
						if (symop)
							result.emplace_back(a, loc, symop.string());
						else
							result.emplace_back(a);
					}
				}
			}
		}
	}

	std::ranges::sort(result, [](auto &a, auto &b)
		{ return a.id().compare(b.id()) < 0; });
	auto r = std::ranges::unique(result);

	if (r.begin() != r.end())
		result.erase(r.begin(), r.end());

	return result;
}

} // namespace pdb_redo
