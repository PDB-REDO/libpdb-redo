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
#include <cif++/compound.hpp>
#include <cif++/condition.hpp>
#include <cif++/model.hpp>
#include <cif++/point.hpp>
#include <cif++/symmetry.hpp>
#include <cif++/utilities.hpp>
#include <cif++/validate.hpp>
#include <cstddef>
#include <ranges>
#include <unordered_map>
#define CATCH_CONFIG_RUNNER

#include "pdb-redo/BondMap.hpp"

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

class BondMapException : public std::runtime_error
{
  public:
	BondMapException(const std::string &msg)
		: runtime_error(msg)
	{
	}
};

class BondMap
{
  public:
	BondMap(const cif::datablock &db, std::optional<std::tuple<cif::point, float>> around = {}, std::size_t model_nr = 1);

	BondMap(const BondMap &) = delete;
	BondMap &operator=(const BondMap &) = delete;

	BondMap(BondMap &&);
	BondMap &operator=(BondMap &&);

	bool operator()(const std::string &atom_1, const std::string &atom_2) const
	{
		auto aix1 = index.find(atom_1);
		auto aix2 = index.find(atom_2);
		return aix1 != index.end() and aix2 != index.end() and isBonded(aix1->second, aix2->second);
	}

	bool operator()(const cif::mm::atom &atom_1, const cif::mm::atom &atom_2) const
	{
		return operator()(atom_1.id(), atom_2.id());
	}

	bool is1_4(const std::string &atom_1, const std::string &atom_2) const
	{
		uint32_t ixa = index.at(atom_1);
		uint32_t ixb = index.at(atom_2);

		return bond_1_4.count(key(ixa, ixb));
	}

	bool is1_4(const cif::mm::atom &atom_1, const cif::mm::atom &atom_2) const
	{
		return is1_4(atom_1.id(), atom_2.id());
	}

	// links coming from the struct_conn records:
	std::vector<std::string> linked(const std::string &atom) const;

	// This list of atomID's is comming from either CCD or the CCP4 dictionaries loaded
	static std::vector<std::string> atomIDsForCompound(const std::string &compoundID);

	//   private:
	constexpr std::tuple<uint32_t, uint32_t> key(uint32_t a, uint32_t b) const
	{
		if (a > b)
			std::swap(a, b);
		return { a, b };
	}

	constexpr bool isBonded(uint32_t ai, uint32_t bi) const
	{
		if (ai > bi)
			std::swap(ai, bi);

		return bond.count({ ai, bi }) != 0;
	}

	uint32_t dim;
	std::unordered_map<std::string, uint32_t> index;
	std::set<std::tuple<uint32_t, uint32_t>> bond, bond_1_4;

	std::map<std::string, std::set<std::string>> link;
};

// --------------------------------------------------------------------

struct CompoundBondInfo
{
	std::string mID;
	std::set<std::tuple<uint32_t, uint32_t>> mBonded;

	[[nodiscard]] bool bonded(uint32_t a1, uint32_t a2) const
	{
		return mBonded.count({ a1, a2 }) > 0;
	}
};

// --------------------------------------------------------------------

class CompoundBondMap
{
  public:
	static CompoundBondMap &instance()
	{
		static std::unique_ptr<CompoundBondMap> s_instance(new CompoundBondMap);
		return *s_instance;
	}

	bool bonded(const std::string &compoundID, const std::string &atomID1, const std::string &atomID2);

  private:
	CompoundBondMap() = default;

	uint32_t getAtomID(const std::string &atomID)
	{
		std::string id(atomID);

		uint32_t result;

		auto i = mAtomIDIndex.find(id);
		if (i == mAtomIDIndex.end())
		{
			result = static_cast<uint32_t>(mAtomIDIndex.size());
			mAtomIDIndex[id] = result;
		}
		else
			result = i->second;

		return result;
	}

	std::map<std::string, uint32_t> mAtomIDIndex;
	std::vector<CompoundBondInfo> mCompounds;
	std::mutex mMutex;
};

bool CompoundBondMap::bonded(const std::string &compoundID, const std::string &atomID1, const std::string &atomID2)
{
	std::scoped_lock lock(mMutex);

	using namespace std::literals;

	std::string id(compoundID);
	uint32_t a1 = getAtomID(atomID1);
	uint32_t a2 = getAtomID(atomID2);
	if (a1 > a2)
		std::swap(a1, a2);

	for (auto &bi : mCompounds)
	{
		if (bi.mID != id)
			continue;

		return bi.bonded(a1, a2);
	}

	bool result = false;

	// not found in our cache, calculate
	CompoundBondInfo bondInfo{ id };

	auto compound = cif::compound_factory::instance().create(compoundID);
	if (not compound)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Missing compound bond info for " << compoundID << '\n';
	}
	else
	{
		for (auto &atom : compound->bonds())
		{
			uint32_t ca1 = getAtomID(atom.atom_id[0]);
			uint32_t ca2 = getAtomID(atom.atom_id[1]);
			if (ca1 > ca2)
				std::swap(ca1, ca2);

			bondInfo.mBonded.insert({ ca1, ca2 });
			result = result or (a1 == ca1 and a2 == ca2);
		}
	}

	mCompounds.push_back(bondInfo);

	return result;
}

// --------------------------------------------------------------------

BondMap::BondMap(BondMap &&bm)
	: dim(bm.dim)
	, index(std::move(bm.index))
	, bond(std::move(bm.bond))
	, bond_1_4(std::move(bm.bond_1_4))
	, link(std::move(bm.link))
{
}

BondMap &BondMap::operator=(BondMap &&bm)
{
	dim = bm.dim;
	index.swap(bm.index);
	bond.swap(bm.bond);
	bond_1_4.swap(bm.bond_1_4);
	link.swap(bm.link);

	return *this;
}

BondMap::BondMap(const cif::datablock &db, std::optional<std::tuple<cif::point, float>> around, std::size_t model_nr)
{
	using namespace cif::literals;

	auto &compoundBondInfo = CompoundBondMap::instance();

	// First collect the atoms from the datablock
	std::vector<cif::const_row_handle> atoms;

	cif::crystal crystal(db);

	for (auto rh : db["atom_site"].find("pdbx_PDB_model_num"_key == model_nr or "pdbx_PDB_model_num"_key == cif::null))
	{
		if (around)
		{
			const auto &[p, r] = *around;
			const auto &[x, y, z] = rh.get<float, float, float>("Cartn_x", "Cartn_y", "Cartn_z");

			const auto &[d, pt, op] = crystal.closest_symmetry_copy(p, { x, y, z });

			if (d <= r)
				atoms.push_back(rh);
		}
		else
			atoms.push_back(rh);
	}

	dim = static_cast<uint32_t>(atoms.size());

	for (auto &atom : atoms)
		index[atom["id"].get<std::string>()] = static_cast<uint32_t>(index.size());

	auto bindAtoms = [this](const std::string &a, const std::string &b)
	{
		uint32_t ixa = index[a];
		uint32_t ixb = index[b];

		bond.insert(key(ixa, ixb));
	};

	auto linkAtoms = [this, &bindAtoms](const std::string &a, const std::string &b)
	{
		bindAtoms(a, b);

		link[a].insert(b);
		link[b].insert(a);
	};

	// collect all compounds first
	std::set<std::string> compounds;
	for (const auto &comp_id : db["chem_comp"].rows<std::string>("id"))
		compounds.insert(comp_id);

	// make sure we also have all residues in the polyseq
	for (const auto &mon_id : db["entity_poly_seq"].rows<std::string>("mon_id"))
	{
		if (compounds.count(mon_id))
			continue;

		if (cif::VERBOSE > 1)
			std::cerr << "Warning: mon_id " << mon_id << " is missing in the chem_comp category\n";
		compounds.insert(mon_id);
	}

	cif::progress_bar progress_bar(compounds.size(), "Creating bond map");

	// some helper indices to speed things up a bit
	using atom_map_key_type = std::tuple<std::string, int, std::string, std::string>;
	std::map<atom_map_key_type, std::string> atomMapByAsymSeqAndAtom;
	for (auto a : atoms)
	{
		atom_map_key_type key = a.get("label_asym_id", "label_seq_id", "label_atom_id", "auth_seq_id");
		atomMapByAsymSeqAndAtom[key] = a.get<std::string>("id");
	}

	// first link all residues in a polyseq

	std::string lastAsymID, lastAuthSeqID;
	int lastSeqID = 0;
	for (const auto &[asymID, seqID, authSeqID] : db["pdbx_poly_seq_scheme"].rows<std::string, int, std::string>("asym_id", "seq_id", "pdb_seq_num"))
	{
		if (asymID != lastAsymID) // first in a new sequece
		{
			lastAsymID = asymID;
			lastSeqID = seqID;
			lastAuthSeqID = authSeqID;
			continue;
		}

		auto kc = make_tuple(asymID, lastSeqID, "C", lastAuthSeqID);
		auto kn = make_tuple(asymID, seqID, "N", authSeqID);

		if (atomMapByAsymSeqAndAtom.count(kc) and atomMapByAsymSeqAndAtom.count(kn))
		{
			auto c = atomMapByAsymSeqAndAtom.at(kc);
			auto n = atomMapByAsymSeqAndAtom.at(kn);

			bindAtoms(c, n);
		}
		// if (not(c.empty() or n.empty()))

		lastSeqID = seqID;
		lastAuthSeqID = authSeqID;
	}

	for (auto l : db["struct_conn"])
	{
		atom_map_key_type ka = l.get("ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id", "ptnr1_auth_seq_id");
		atom_map_key_type kb = l.get("ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id", "ptnr2_auth_seq_id");

		if (atomMapByAsymSeqAndAtom.count(ka) and atomMapByAsymSeqAndAtom.count(kb))
		{
			auto a = atomMapByAsymSeqAndAtom.at(ka);
			auto b = atomMapByAsymSeqAndAtom.at(kb);

			linkAtoms(a, b);
		}
	}

	// then link all atoms in the compounds

	// std::vector<cif::const_row_handle> rAtoms;
	// std::string lastEntityID, lastCompID;

	// lastAsymID.clear();
	// lastAuthSeqID.clear();

	// for (auto atom : atoms)
	// {
	// 	const auto [asym_id, seq_id, auth_seq_id, entity_id, comp_id] = atom.get<std::string, int, std::string, std::string, std::string>(
	// 		"label_asym_id", "label_seq_id", "auth_seq_id", "entity_id", "label_comp_id");

	// 	if (asym_id == lastAsymID and entity_id == lastEntityID and lastCompID == comp_id and lastSeqID == seq_id and lastAuthSeqID == auth_seq_id)
	// 	{
	// 		rAtoms.emplace_back(atom);
	// 		continue;
	// 	}

	// 	if (not rAtoms.empty())
	// 	{
	// 		for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
	// 		{
	// 			auto id_i = rAtoms[i].get<std::string>("id");
	// 			auto atom_id_i = rAtoms[i].get<std::string>("label_atom_id");

	// 			for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
	// 			{
	// 				auto atom_id_j = rAtoms[j].get<std::string>("label_atom_id");

	// 				if (compoundBondInfo.bonded(comp_id, atom_id_i, atom_id_j))
	// 					bindAtoms(id_i, rAtoms[j].get<std::string>("id"));
	// 			}
	// 		}
	// 	}

	// 	rAtoms = { atom };
	// 	lastAsymID = asym_id;
	// 	lastEntityID = entity_id;
	// 	lastCompID = comp_id;
	// 	lastSeqID = seq_id;
	// 	lastAuthSeqID = auth_seq_id;
	// }

	// for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
	// {
	// 	auto id_i = rAtoms[i].get<std::string>("id");
	// 	auto atom_id_i = rAtoms[i].get<std::string>("label_atom_id");

	// 	for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
	// 	{
	// 		auto atom_id_j = rAtoms[j].get<std::string>("label_atom_id");

	// 		if (compoundBondInfo.bonded(lastCompID, atom_id_i, atom_id_j))
	// 			bindAtoms(id_i, rAtoms[j].get<std::string>("id"));
	// 	}
	// }

	for (auto c : compounds)
	{
		progress_bar.consumed(1);

		if (c == "HOH" or c == "H2O" or c == "WAT")
		{
			if (cif::VERBOSE > 1)
				std::cerr << "skipping water in bond map calculation\n";
			continue;
		}

		auto bonded = [c, &compoundBondInfo](cif::const_row_handle a, cif::const_row_handle b)
		{
			auto label_a = a.get<std::string>("label_atom_id");
			auto label_b = b.get<std::string>("label_atom_id");

			return compoundBondInfo.bonded(c, label_a, label_b);
		};

		// loop over poly_seq_scheme
		for (const auto &[asymID, seqID] : db["pdbx_poly_seq_scheme"].find<std::string, int>(cif::key("mon_id") == c, "asym_id", "seq_id"))
		{
			std::vector<cif::const_row_handle> rAtoms;
			std::ranges::copy_if(atoms, back_inserter(rAtoms),
				[asymID = asymID, seqID = seqID](cif::const_row_handle a)
				{ return a["label_asym_id"] == asymID and a["label_seq_id"] == seqID; });

			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
						bindAtoms(rAtoms[i].get<std::string>("id"), rAtoms[j].get<std::string>("id"));
				}
			}
		}

		// loop over pdbx_nonpoly_scheme
		for (auto r : db["pdbx_nonpoly_scheme"].find(cif::key("mon_id") == c))
		{
			std::string asymID;
			cif::tie(asymID) = r.get("asym_id");

			std::vector<cif::const_row_handle> rAtoms;
			std::ranges::copy_if(atoms, back_inserter(rAtoms),
				[&](cif::const_row_handle a)
				{ return a["label_asym_id"] == asymID; });

			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
					{
						uint32_t ixa = index[rAtoms[i].get<std::string>("id")];
						uint32_t ixb = index[rAtoms[j].get<std::string>("id")];

						bond.insert(key(ixa, ixb));
					}
				}
			}
		}

		// loop over pdbx_branch_scheme
		for (const auto &[asym_id, pdb_seq_num] : db["pdbx_branch_scheme"].find<std::string, std::string>(cif::key("mon_id") == c, "asym_id", "pdb_seq_num"))
		{
			std::vector<cif::const_row_handle> rAtoms;
			std::ranges::copy_if(atoms, back_inserter(rAtoms),
				[id = asym_id, nr = pdb_seq_num](cif::const_row_handle a)
				{ return a["label_asym_id"] == id and a["auth_seq_id"] == nr; });

			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
					{
						uint32_t ixa = index[rAtoms[i].get<std::string>("id")];
						uint32_t ixb = index[rAtoms[j].get<std::string>("id")];

						bond.insert(key(ixa, ixb));
					}
				}
			}
		}
	}

	// start by creating an index for single bonds

	std::multimap<uint32_t, uint32_t> b1_2;
	for (auto [a, b] : bond)
	{
		b1_2.insert({ a, b });
		b1_2.insert({ b, a });
	}

	std::multimap<uint32_t, uint32_t> b1_3;
	for (uint32_t i = 0; i < dim; ++i)
	{
		auto a = b1_2.equal_range(i);

		std::vector<uint32_t> s;
		for (auto j = a.first; j != a.second; ++j)
			s.push_back(j->second);

		for (std::size_t si1 = 0; si1 + 1 < s.size(); ++si1)
		{
			for (std::size_t si2 = si1 + 1; si2 < s.size(); ++si2)
			{
				uint32_t x = s[si1];
				uint32_t y = s[si2];

				if (isBonded(x, y))
					continue;

				b1_3.insert({ x, y });
				b1_3.insert({ y, x });
			}
		}
	}

	for (uint32_t i = 0; i < dim; ++i)
	{
		auto a1 = b1_2.equal_range(i);
		auto a2 = b1_3.equal_range(i);

		for (auto ai1 = a1.first; ai1 != a1.second; ++ai1)
		{
			for (auto ai2 = a2.first; ai2 != a2.second; ++ai2)
			{
				uint32_t b1 = ai1->second;
				uint32_t b2 = ai2->second;

				if (isBonded(b1, b2))
					continue;

				bond_1_4.insert(key(b1, b2));
			}
		}
	}
}

std::vector<std::string> BondMap::linked(const std::string &atom_id) const
{
	auto i = link.find(atom_id);

	std::vector<std::string> result;

	if (i != link.end())
		result = std::vector<std::string>(i->second.begin(), i->second.end());

	return result;
}

std::vector<std::string> BondMap::atomIDsForCompound(const std::string &compoundID)
{
	std::vector<std::string> result;

	auto *compound = cif::compound_factory::instance().create(compoundID);

	if (compound == nullptr)
		throw BondMapException("Missing bond information for compound " + compoundID);

	for (auto &compAtom : compound->atoms())
		result.push_back(compAtom.id);

	return result;
}

// --------------------------------------------------------------------

TEST_CASE("test_0")
{
	std::cout << "reading file..." << std::flush;
	// cif::file f(gTestDir / "8p4v_final.cif");
	// cif::file f(gTestDir / "4hea_final.cif");
	// cif::file f(gTestDir / "1dex.cif.gz");
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	std::cout << " loading dictionary..." << std::flush;
	f.front().load_dictionary("mmcif_pdbx.dic");
	std::cout << " building structure..." << std::flush;
	cif::mm::structure s(f);
	std::cout << " done\n";

	auto n1 = std::chrono::system_clock::now();
	pdb_redo::BondMap bm1(f.front());
	std::cout << "bm1 took: " << std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now() - n1) << "\n";

	auto n0 = std::chrono::system_clock::now();
	BondMap bm2(f.front());
	std::cout << "bm2 took: " << std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now() - n0) << "\n";

	// CHECK(bm1.dim == bm2.dim);
	// CHECK(bm1.index == bm2.index);
	// // CHECK(bm1.bond_1_4 == bm2.bond_1_4);
	// CHECK(bm1.link == bm2.link);
}

// TEST_CASE("test_1")
// {
// 	using namespace cif::literals;

// 	std::cout << "reading file..." << std::flush;
// 	// cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
// 	cif::file f(gTestDir / "2b8h.cif.gz");
// 	std::cout << " loading dictionary..." << std::flush;
// 	f.front().load_dictionary("mmcif_pdbx.dic");
// 	std::cout << " building structure..." << std::flush;
// 	cif::mm::structure s(f);
// 	std::cout << " done\n";

// 	auto &struct_conn = f.front()["struct_conn"];

// 	pdb_redo::BondMap bm(f.front());

// 	cif::progress_bar pb(s.atoms().size() * (s.atoms().size() - 1), "testing");

// 	for (size_t ixa = 0; ixa + 1 < s.atoms().size(); ++ixa)
// 	{
// 		auto a = s.atoms()[ixa];

// 		auto compound_a = cif::compound_factory::instance().create(a.get_label_comp_id());
// 		REQUIRE(compound_a != nullptr);

// 		for (size_t ixb = ixa + 1; ixb < s.atoms().size(); ++ixb)
// 		{
// 			auto b = s.atoms()[ixb];

// 			pb.consumed(1);

// 			auto compound_b = cif::compound_factory::instance().create(a.get_label_comp_id());
// 			REQUIRE(compound_b != nullptr);

// 			bool bonded = false;

// 			// std::cout << "a: " << a << "b: " << b << "\n";

// 			// Same asym?
// 			if (a.get_label_asym_id() == b.get_label_asym_id())
// 			{
// 				// same residue?
// 				if (a.get_label_seq_id() == b.get_label_seq_id() and
// 					a.get_label_comp_id() == b.get_label_comp_id() and
// 					a.get_auth_seq_id() == b.get_auth_seq_id())
// 				{
// 					CHECK(compound_a == compound_b);
// 					bonded = compound_a->atoms_bonded(a.get_label_atom_id(), b.get_label_atom_id());
// 				}
// 				else if (a.get_label_seq_id() == b.get_label_seq_id() - 1)
// 					bonded = a.get_label_atom_id() == "C" and b.get_label_atom_id() == "N";
// 				else if (b.get_label_seq_id() == a.get_label_seq_id() - 1)
// 					bonded = b.get_label_atom_id() == "C" and a.get_label_atom_id() == "N";
// 			}

// 			if (not bonded)
// 			{
// 				if (auto r = struct_conn.find_first(

// 						"ptnr1_label_asym_id"_key == a.get_label_asym_id() and
// 						"ptnr1_label_comp_id"_key == a.get_label_comp_id() and
// 						"ptnr1_label_seq_id"_key == a.get_label_seq_id() and
// 						"ptnr1_label_atom_id"_key == a.get_label_atom_id() and

// 						"ptnr2_label_asym_id"_key == b.get_label_asym_id() and
// 						"ptnr2_label_comp_id"_key == b.get_label_comp_id() and
// 						"ptnr2_label_seq_id"_key == b.get_label_seq_id() and
// 						"ptnr2_label_atom_id"_key == b.get_label_atom_id()))
// 				{
// 					bonded = true;
// 				}
// 				else if (auto r = struct_conn.find_first(

// 						"ptnr1_label_asym_id"_key == b.get_label_asym_id() and
// 						"ptnr1_label_comp_id"_key == b.get_label_comp_id() and
// 						"ptnr1_label_seq_id"_key == b.get_label_seq_id() and
// 						"ptnr1_label_atom_id"_key == b.get_label_atom_id() and

// 						"ptnr2_label_asym_id"_key == a.get_label_asym_id() and
// 						"ptnr2_label_comp_id"_key == a.get_label_comp_id() and
// 						"ptnr2_label_seq_id"_key == a.get_label_seq_id() and
// 						"ptnr2_label_atom_id"_key == a.get_label_atom_id()))
// 				{
// 					bonded = true;
// 				}
// 			}

// 			CHECK(bm(a.id(), b.id()) == bonded);

// 			if (bm(a.id(), b.id()) != bonded)
// 			{
// 				std::cout << "a: " << a << " and b: " << b << " should " << (bonded ? "" : "not ") << "be bonded\n";
// 			}
// 		}
// 	}
// }

TEST_CASE("test-iter")
{
	using namespace cif::literals;

	std::cout << "reading file..." << std::flush;
	cif::file f(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	// cif::file f(gTestDir / "2b8h.cif.gz");
	std::cout << " loading dictionary..." << std::flush;
	f.front().load_dictionary("mmcif_pdbx.dic");
	// std::cout << " building structure..." << std::flush;
	// cif::mm::structure s(f);
	std::cout << " done\n";

	auto &struct_conn = f.front()["struct_conn"];

	pdb_redo::BondMap bm(f.front());

	for (const auto &[one, two] : bm)
	{
		CHECK(bm(one, two));
		if (not bm(one, two))
		{
			println(std::cout, "not bonded? {} and {}", one, two);
			break;
		}
	}
	
}