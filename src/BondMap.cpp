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

#include <cassert>

#include <algorithm>
#include <fstream>
#include <mutex>

#include <cif++.hpp>

#include <pdb-redo/BondMap.hpp>

namespace pdb_redo
{

// --------------------------------------------------------------------

struct CompoundBondInfo
{
	std::string mID;
	std::set<std::tuple<uint32_t, uint32_t>> mBonded;

	bool bonded(uint32_t a1, uint32_t a2) const
	{
		return mBonded.count({a1, a2}) > 0;
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
	CompoundBondMap() {}

	uint32_t getAtomID(const std::string &atomID)
	{
		std::string id(atomID);

		uint32_t result;

		auto i = mAtomIDIndex.find(id);
		if (i == mAtomIDIndex.end())
		{
			result = uint32_t(mAtomIDIndex.size());
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
	std::lock_guard lock(mMutex);

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
	CompoundBondInfo bondInfo{id};

	auto compound = cif::compound_factory::instance().create(compoundID);
	if (not compound)
	{
		if (cif::VERBOSE >= 0)
			std::cerr << "Missing compound bond info for " << compoundID << std::endl;
	}
	else
	{
		for (auto &atom : compound->bonds())
		{
			uint32_t ca1 = getAtomID(atom.atom_id[0]);
			uint32_t ca2 = getAtomID(atom.atom_id[1]);
			if (ca1 > ca2)
				std::swap(ca1, ca2);

			bondInfo.mBonded.insert({ca1, ca2});
			result = result or (a1 == ca1 and a2 == ca2);
		}
	}

	mCompounds.push_back(bondInfo);

	return result;
}

// --------------------------------------------------------------------

BondMap::BondMap(const cif::datablock &db, int model_nr)
{
	using namespace cif::literals;

	auto &compoundBondInfo = CompoundBondMap::instance();

	// First collect the atoms from the datablock
	std::vector<cif::row_handle> atoms;

	for (auto rh : db["atom_site"].find("pdbx_PDB_model_num"_key == model_nr or "pdbx_PDB_model_num"_key == cif::null))
		atoms.push_back(rh);

	dim = uint32_t(atoms.size());

	for (auto &atom : atoms)
		index[atom["id"].as<std::string>()] = uint32_t(index.size());

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
			std::cerr << "Warning: mon_id " << mon_id << " is missing in the chem_comp category" << std::endl;
		compounds.insert(mon_id);
	}

	cif::Progress progress(compounds.size(), "Creating bond map");

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

	for (auto c : compounds)
	{
		if (c == "HOH" or c == "H2O" or c == "WAT")
		{
			if (cif::VERBOSE > 0)
				std::cerr << "skipping water in bond map calculation" << std::endl;
			continue;
		}

		auto bonded = [c, &compoundBondInfo](cif::row_handle a, cif::row_handle b)
		{
			auto label_a = a.get<std::string>("label_atom_id");
			auto label_b = b.get<std::string>("label_atom_id");

			return compoundBondInfo.bonded(c, label_a, label_b);
		};

		// loop over poly_seq_scheme
		for (const auto &[asymID, seqID] : db["pdbx_poly_seq_scheme"].find<std::string, int>(cif::key("mon_id") == c, "asym_id", "seq_id"))
		{
			std::vector<cif::row_handle> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[asymID=asymID,seqID=seqID](cif::row_handle a)
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

			std::vector<cif::row_handle> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](cif::row_handle a)
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
			std::vector<cif::row_handle> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[id = asym_id, nr = pdb_seq_num](cif::row_handle a)
				{ return a["label_asym_id"] == id and a["label_seq_id"] == nr; });

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
	for (auto &bk : bond)
	{
		uint32_t a, b;
		std::tie(a, b) = dekey(bk);

		b1_2.insert({a, b});
		b1_2.insert({b, a});
	}

	std::multimap<uint32_t, uint32_t> b1_3;
	for (uint32_t i = 0; i < dim; ++i)
	{
		auto a = b1_2.equal_range(i);

		std::vector<uint32_t> s;
		for (auto j = a.first; j != a.second; ++j)
			s.push_back(j->second);

		for (size_t si1 = 0; si1 + 1 < s.size(); ++si1)
		{
			for (size_t si2 = si1 + 1; si2 < s.size(); ++si2)
			{
				uint32_t x = s[si1];
				uint32_t y = s[si2];

				if (isBonded(x, y))
					continue;

				b1_3.insert({x, y});
				b1_3.insert({y, x});
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

} // namespace pdbx
