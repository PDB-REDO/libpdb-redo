/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2026 NKI/AVL, Netherlands Cancer Institute
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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cif++/condition.hpp>
#include <cif++/model.hpp>
#include <cif++/point.hpp>
#include <cif++/symmetry.hpp>
#include <cif++/utilities.hpp>
#include <cif++/validate.hpp>
#include <ranges>
#include <unordered_map>
#define CATCH_CONFIG_RUNNER

#include "pdb-redo/DistanceMap.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cif++/cif++.hpp>
#include <filesystem>
#include <utility>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

std::filesystem::path gTestDir;

int main(int argc, char *argv[])
{
	gTestDir = std::filesystem::current_path();

	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
#if CATCH22
	using namespace Catch::clara;
#else
	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;
#endif

	auto cli = session.cli()                                // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                  // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                   // the option names it will respond to
	           ("The directory containing the data files"); // description string for the help output

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	if (fs::exists(gTestDir / "minimal-components.cif"))
		cif::compound_factory::instance().push_dictionary(gTestDir / "minimal-components.cif");

	return session.run();
}

// --------------------------------------------------------------------

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

DistanceMap::DistanceMap(const cif::mm::structure &structure, cif::crystal crystal, float maxDistance)
	: m_structure(structure)
	, m_crystal(std::move(crystal))
	, m_grid_spacing(1)
{
	std::vector<std::tuple<cif::point, std::string>> pts;

	pts.reserve(structure.atoms().size());

	key_type k1{}, k2{};

	for (auto a : structure.atoms())
	{
		pts.emplace_back(a.get_location(), a.id());

		auto p = a.get_location();
		key_type k{
			static_cast<int>(std::rint(p.x / m_grid_spacing)),
			static_cast<int>(std::rint(p.y / m_grid_spacing)),
			static_cast<int>(std::rint(p.z / m_grid_spacing))
		};

		if (m_index.empty())
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

		m_index.emplace(k, entry{ a.id(), cif::sym_op{} });
	}

	int d = static_cast<int>(std::rint(maxDistance / m_grid_spacing));
	k1.x -= d;
	k2.x += d;
	k1.y -= d;
	k2.y += d;
	k1.z -= d;
	k2.z += d;

	auto &sg = m_crystal.get_spacegroup();
	auto &cell = m_crystal.get_cell();

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

						key_type k{
							static_cast<int>(std::rint(ap.x / m_grid_spacing)),
							static_cast<int>(std::rint(ap.y / m_grid_spacing)),
							static_cast<int>(std::rint(ap.z / m_grid_spacing))
						};

						if (k.x >= k1.x and k.x <= k2.x and
							k.y >= k1.y and k.y <= k2.y and
							k.y >= k1.z and k.z <= k2.z)
						{
							m_index.emplace(k, entry{ id, symop });
						}
					}
				}
			}
		}
	}
}

float DistanceMap::operator()(const std::string &a, const std::string &b) const
{
	return cif::distance(m_structure.get_atom_by_id(a).get_location(), m_structure.get_atom_by_id(b).get_location());
}

std::vector<cif::mm::atom> DistanceMap::near(const cif::mm::atom &atom, float maxDistance) const
{
	std::vector<cif::mm::atom> result;

	auto maxDistanceSq = maxDistance * maxDistance;

	auto p = atom.get_location();

	key_type k{
		static_cast<int>(std::rint(p.x / m_grid_spacing)),
		static_cast<int>(std::rint(p.y / m_grid_spacing)),
		static_cast<int>(std::rint(p.z / m_grid_spacing))
	};

	key_type k1 = k, k2 = k;

	int d = static_cast<int>(maxDistance / m_grid_spacing) + 1;
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
				auto r = m_index.equal_range(k);
				for (auto &[k, e] : std::ranges::subrange(r.first, r.second))
				{
					auto [id, symop] = e;
					if (id == atom.id())
						continue;

					auto a = m_structure.get_atom_by_id(id);
					auto loc = m_crystal.symmetry_copy(a.get_location(), symop);

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

// --------------------------------------------------------------------

TEST_CASE("test_0")
{
	std::cout << "reading file..." << std::flush;
	// cif::file f(gTestDir / "8p4v_final.cif");
	// cif::file f(gTestDir / "4hea_final.cif");
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	std::cout << " loading dictionary..." << std::flush;
	f.front().load_dictionary("mmcif_pdbx.dic");
	std::cout << " building structure..." << std::flush;
	cif::mm::structure s(f);
	std::cout << " done\n";

	auto n1 = std::chrono::system_clock::now();
	pdb_redo::DistanceMap dm1(s, 3.5f);
	std::cout << "dm1 took: " << std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now() - n1) << "\n";

	auto n0 = std::chrono::system_clock::now();
	DistanceMap dm2(s, 3.5f);
	std::cout << "dm2 took: " << std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now() - n0) << "\n";

	for (auto a : s.atoms())
	{
		auto n1 = dm1.near(a);
		auto n2 = dm2.near(a);

		std::erase_if(n1, [a](const cif::mm::atom &b)
			{ return cif::distance(a.get_location(), b.get_location()) > 3.5f; });
		std::ranges::sort(n1, [](auto &a, auto &b)
			{ return a.id().compare(b.id()) < 0; });

		std::vector<cif::mm::atom> oi1, oi2;

		auto b1 = n1.begin(), b2 = n2.begin();
		while (b1 != n1.end() and b2 != n2.end())
		{
			if (*b1 == *b2)
				++b1, ++b2;
			else if (b1->id().compare(b2->id()) < 0)
				oi1.emplace_back(*b1++);
			else
				oi2.emplace_back(*b1++);
		}

		while (b1 != n1.end())
			oi1.emplace_back(*b1++);
		while (b2 != n2.end())
			oi2.emplace_back(*b2++);

		if (oi1.empty() and oi2.empty())
			continue;

		CHECK(oi1.size() == 0);

		if (oi1.size() == 0)
			continue;

		std::cout << "For a = " << a.id() << ": " << a << " @ " << a.symmetry() << " d: " << cif::distance(a.get_location(), a.get_location()) << ' ' << a.get_location() << "\n";

		std::cout << "only in n1:\n";

		for (auto ai : oi1)
			std::cout << ai.id() << ": " << ai << " @ " << ai.symmetry() << " d: " << cif::distance(a.get_location(), ai.get_location()) << ' ' << ai.get_location() << "\n";

		// std::cout << "only in n2:\n";
		// for (auto ai : oi2)
		// 	std::cout << ai.id() << ": " << ai << " @ " << ai.symmetry() << " d: " << cif::distance(a.get_location(), ai.get_location()) << ' ' << ai.get_location() << "\n";
	}

	// std::cout << "test: " << N << " less: " << M << " more: "
	// std::println(std::cout, "test: {}, less: {}, more: {}\n", N, M, O);
	// std::cout << std::format("test: {}, less: {}, more: {}\n", N, M, O);
}
