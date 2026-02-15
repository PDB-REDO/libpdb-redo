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
#include <cif++/model.hpp>
#include <cif++/point.hpp>
#include <cif++/symmetry.hpp>
#include <cif++/utilities.hpp>
#include <cif++/validate.hpp>
#include <ranges>
#include <unordered_map>
#define CATCH_CONFIG_RUNNER

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Minimizer.hpp"
#include "pdb-redo/Statistics.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cif++.hpp>
#include <filesystem>
#include <stdexcept>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

std::filesystem::path gTestDir = std::filesystem::current_path();

int main(int argc, char *argv[])
{
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
	DistanceMap(const cif::mm::structure &p, const cif::crystal &crystal, float maxDistance);

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

	using key_type = std::tuple<int, int, int>;

	struct key_type_hash
	{
		std::size_t operator()(const key_type &s) const noexcept
		{
			auto h0 = std::hash<int>{}(std::get<0>(s));
			auto h1 = std::hash<int>{}(std::get<1>(s));
			auto h2 = std::hash<int>{}(std::get<2>(s));

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

DistanceMap::DistanceMap(const cif::mm::structure &structure, const cif::crystal &crystal, float maxDistance)
	: m_structure(structure)
	, m_crystal(crystal)
	, m_grid_spacing(1)
{
	std::vector<cif::point> pts;
	pts.reserve(structure.atoms().size());

	for (auto a : structure.atoms())
		pts.emplace_back(a.get_location());

	auto [center, radius] = cif::smallest_sphere_around_points(pts);

	// radius += maxDistance;

	// auto maxDistanceSq = maxDistance * maxDistance;
	auto maxAtomDistance = radius + maxDistance + m_grid_spacing;
	maxAtomDistance *= maxAtomDistance;

	// std::vector<std::tuple<std::string, cif::point, float>> asymSpheres;
	std::map<std::string, std::vector<cif::point>> asyms;
	std::map<std::string, std::vector<std::string>> atomIDs;

	cif::progress_bar progress_bar(structure.atoms().size(), "Creating distance map");

	for (auto a : structure.atoms())
	{
		auto p = a.get_location();
		key_type k{
			static_cast<int>(std::rint(p.m_x / m_grid_spacing)),
			static_cast<int>(std::rint(p.m_y / m_grid_spacing)),
			static_cast<int>(std::rint(p.m_z / m_grid_spacing))
		};
		m_index.emplace(k, entry{ a.id(), cif::sym_op{} });

		asyms[a.get_label_asym_id()].emplace_back(p);
		atomIDs[a.get_label_asym_id()].emplace_back(a.id());

		progress_bar.consumed(1);
	}

	std::cout << "Number of asyms: " << asyms.size() << "\n";

	int N = 0;
	for (auto &[asym_id, points] : asyms)
	{
		// auto [asym_center, asym_radius] = cif::smallest_sphere_around_points(points);

		// auto dsq = asym_radius + radius + maxDistance + m_grid_spacing;
		// dsq *= dsq;

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

						// auto p = sg(asym_center, cell, symop);

// if (i == 4 and tx == 5 and ty == 5 and tz == 4 and asym_id == "C")
// 						std::cout << "letop\n";

						// if (cif::distance_squared(p, center) > dsq)
						// 	continue;

						++N;
						auto &ids = atomIDs[asym_id];

						for (size_t ix = 0; ix < points.size(); ++ix)
						{
							auto ap = sg(points[ix], cell, symop);

// if (i == 4 and tx == 5 and ty == 5 and tz == 4 and asym_id == "C" and ids[ix] == "1132")
// 						std::cout << "letop\n";


							if (cif::distance_squared(ap, center) <= maxAtomDistance)
							{
								key_type k{
									static_cast<int>(std::rint(ap.m_x / m_grid_spacing)),
									static_cast<int>(std::rint(ap.m_y / m_grid_spacing)),
									static_cast<int>(std::rint(ap.m_z / m_grid_spacing))
								};
								m_index.emplace(k, entry{ ids[ix], symop });
							}
						}
					}
				}
			}
		}
	}

	std::cout << "N: " << N << "\n";
}

float DistanceMap::operator()(const std::string &a, const std::string &b) const
{
	return cif::distance(m_structure.get_atom_by_id(a).get_location(), m_structure.get_atom_by_id(b).get_location());
}

std::vector<cif::mm::atom> DistanceMap::near(const cif::mm::atom &atom, float maxDistance) const
{
	std::vector<cif::mm::atom> result;

	auto p = atom.get_location();

	for (float x = p.m_x - maxDistance; x <= p.m_x + maxDistance; x += m_grid_spacing)
	{
		for (float y = p.m_y - maxDistance; y <= p.m_y + maxDistance; y += m_grid_spacing)
		{
			for (float z = p.m_z - maxDistance; z <= p.m_z + maxDistance; z += m_grid_spacing)
			{
				key_type k{
					static_cast<int>(std::rint(x / m_grid_spacing)),
					static_cast<int>(std::rint(y / m_grid_spacing)),
					static_cast<int>(std::rint(z / m_grid_spacing))
				};

				auto r = m_index.equal_range(k);
				for (auto &[k, e] : std::ranges::subrange(r.first, r.second))
				{
					auto [id, symop] = e;
					if (id == atom.id())
						continue;

					auto a = m_structure.get_atom_by_id(id);
					auto loc = m_crystal.symmetry_copy(a.get_location(), symop);

					if (auto d = cif::distance(p, loc); d <= maxDistance)
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
	cif::file f(gTestDir / "4hea_final.cif");
	// cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	std::cout << " loading dictionary..." << std::flush;
	f.front().load_dictionary("mmcif_pdbx.dic");
	std::cout << " building structure..." << std::flush;
	cif::mm::structure s(f);
	std::cout << " done\n";

	auto n1 = std::chrono::system_clock::now();
	pdb_redo::DistanceMap dm1(s, 3.5f);

	std::cout << "dm1 took: " << std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now() - n1) << "\n";

	DistanceMap dm2(s, 3.5f);

	std::vector<int> N1(s.atoms().size()), N2(s.atoms().size());

	{
		cif::progress_bar p1(s.atoms().size(), "near in 1");
		for (size_t ix = 0; auto a : s.atoms())
		{
			auto n1 = dm1.near(a);
	
			std::erase_if(n1, [a](const cif::mm::atom &b)
				{ return cif::distance(a.get_location(), b.get_location()) >= 3.5f; });
	
			N1[ix++] = n1.size();
			p1.consumed(1);
		}
	}

	{
		cif::progress_bar p2(s.atoms().size(), "near in 2");
		for (size_t ix = 0; auto a : s.atoms())
		{
			auto n2 = dm2.near(a);
	
			N2[ix++] = n2.size();
			p2.consumed(1);
		}
	}


	int N = N1.size();
	int M = 0, O = 0;
	for (size_t ix = 0; ix < N; ++ix)
	{
		if (N1[ix] > N2[ix])
			++M;
		if (N1[ix] < N2[ix])
			++O;
	}

	// 	auto n2 = dm2.near(a);

	// 	std::erase_if(n1, [a](const cif::mm::atom &b)
	// 		{ return cif::distance(a.get_location(), b.get_location()) >= 3.5f; });
	// 	// std::erase_if(n2, [](const cif::mm::atom &a)
	// 	// 	{ return a.is_symmetry_copy(); });

	// 	++N;
		
	// 	if (n1.size() > n2.size())
	// 		++M;

	// 	if (n2.size() < n2.size())
	// 		++O;

	// 	// CHECK(n1.size() == n2.size());

	// 	// if (n1.size() == n2.size())
	// 	// 	continue;

	// 	// std::cout << "a: " << a.id() << ": " << a << ' ' << a.get_location() << "\n";

	// 	// std::cout << "n1:\n";
	// 	// for (auto ai : n1)
	// 	// 	std::cout << ai.id() << ": " << ai << " @ " << ai.symmetry() << " d: " << cif::distance(a.get_location(), ai.get_location()) << ' ' << ai.get_location() << "\n";
	// 	// std::cout << "n2:\n";
	// 	// for (auto ai : n2)
	// 	// 	std::cout << ai.id() << ": " << ai << " @ " << ai.symmetry() << " d: " << cif::distance(a.get_location(), ai.get_location()) << ' ' << ai.get_location() << "\n";
	// }

	// std::cout << "test: " << N << " less: " << M << " more: "
	// std::println(std::cout, "test: {}, less: {}, more: {}\n", N, M, O);
	std::cout << std::format("test: {}, less: {}, more: {}\n", N, M, O);
}
