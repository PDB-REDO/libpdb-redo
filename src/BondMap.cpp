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

#include "config.hpp"

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/CifUtils.hpp"

#include "pdb-redo/BondMap.hpp"

namespace mmcif
{

// --------------------------------------------------------------------

BondMap::BondMap(const Structure& p)
{
	auto atoms = p.atoms();
	dim = atoms.size();

//	bond = std::vector<bool>(dim * (dim - 1), false);

	for (auto& atom: atoms)
	{
		size_t ix = index.size();
		index[atom.id()] = ix;
	};
	
	auto bindAtoms = [this](const std::string& a, const std::string& b)
	{
		uint32_t ixa = index[a];
		uint32_t ixb = index[b];
		
		bond.insert(key(ixa, ixb));
	};

	auto linkAtoms = [this,&bindAtoms](const std::string& a, const std::string& b)
	{
		bindAtoms(a, b);

		link[a].insert(b);
		link[b].insert(a);
	};

	cif::Datablock& db = p.getFile().data();

	// collect all compounds first
	std::set<std::string> compounds;
	for (auto c: db["chem_comp"])
		compounds.insert(c["id"].as<std::string>());
	
	// make sure we also have all residues in the polyseq
	for (auto m: db["entity_poly_seq"])
	{
		std::string c = m["mon_id"].as<std::string>();
		if (compounds.count(c))
			continue;
		
		if (cif::VERBOSE > 1)
			std::cerr << "Warning: mon_id " << c << " is missing in the chem_comp category" << std::endl;
		compounds.insert(c);
	}

	cif::Progress progress(compounds.size(), "Creating bond map");

	// some helper indices to speed things up a bit
	std::map<std::tuple<std::string,int,std::string>,std::string> atomMapByAsymSeqAndAtom;
	for (auto& a: p.atoms())
	{
		auto key = make_tuple(a.labelAsymID(), a.labelSeqID(), a.labelAtomID());
		atomMapByAsymSeqAndAtom[key] = a.id();
	}
	
	// first link all residues in a polyseq
	
	std::string lastAsymID;
	int lastSeqID = 0;
	for (auto r: db["pdbx_poly_seq_scheme"])
	{
		std::string asymID;
		int seqID;

		cif::tie(asymID, seqID) = r.get("asym_id", "seq_id");

		if (asymID != lastAsymID)		// first in a new sequece
		{
			lastAsymID = asymID;
			lastSeqID = seqID;
			continue;
		}
		
		auto c = atomMapByAsymSeqAndAtom[make_tuple(asymID, lastSeqID, "C")];
		auto n = atomMapByAsymSeqAndAtom[make_tuple(asymID, seqID, "N")];

		if (not (c.empty() or n.empty()))
			bindAtoms(c, n);
		
		lastSeqID = seqID;
	}

	for (auto l: db["struct_conn"])
	{
		std::string asym1, asym2, atomId1, atomId2;
		int seqId1 = 0, seqId2 = 0;
		cif::tie(asym1, asym2, atomId1, atomId2, seqId1, seqId2) =
			l.get("ptnr1_label_asym_id", "ptnr2_label_asym_id",
				  "ptnr1_label_atom_id", "ptnr2_label_atom_id",
				  "ptnr1_label_seq_id", "ptnr2_label_seq_id");

		std::string a = atomMapByAsymSeqAndAtom[make_tuple(asym1, seqId1, atomId1)];
		std::string b = atomMapByAsymSeqAndAtom[make_tuple(asym2, seqId2, atomId2)];
			
		if (not (a.empty() or b.empty()))
			linkAtoms(a, b);
	}

	// then link all atoms in the compounds
	
	for (auto c: compounds)
	{
		auto* compound = mmcif::Compound::create(c);
		if (not compound)
		{
			if (cif::VERBOSE)
				std::cerr << "Missing compound information for " << c << std::endl;
			continue;
		}
		
		if (compound->isWater())
		{
			if (cif::VERBOSE)
				std::cerr << "skipping water in bond map calculation" << std::endl;
			continue;
		}
		
		// loop over poly_seq_scheme
		for (auto r: db["pdbx_poly_seq_scheme"].find(cif::Key("mon_id") == c))
		{
			std::string asymID;
			int seqID;
			cif::tie(asymID, seqID) = r.get("asym_id", "seq_id");
			
			std::vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymID() == asymID and a.labelSeqID() == seqID; });
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (compound->atomsBonded(rAtoms[i].labelAtomID(), rAtoms[j].labelAtomID()))
						bindAtoms(rAtoms[i].id(), rAtoms[j].id());
				}
			}
		}

		// loop over pdbx_nonpoly_scheme
		for (auto r: db["pdbx_nonpoly_scheme"].find(cif::Key("mon_id") == c))
		{
			std::string asymID;
			cif::tie(asymID) = r.get("asym_id");
			
			std::vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymID() == asymID; });
//			for (auto a: db["atom_site"].find(cif::Key("label_asym_id") == asymID))
//				rAtoms.push_back(p.getAtomByID(a["id"].as<std::string>()));
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (compound->atomsBonded(rAtoms[i].labelAtomID(), rAtoms[j].labelAtomID()))
					{
						uint32_t ixa = index[rAtoms[i].id()];
						uint32_t ixb = index[rAtoms[j].id()];
						
						bond.insert(key(ixa, ixb));
					}
				}
			}
		}

		// loop over pdbx_branch_scheme
		for (auto r: db["pdbx_branch_scheme"].find(cif::Key("mon_id") == c))
		{
			std::string asymID;
			cif::tie(asymID) = r.get("asym_id");
			
			std::vector<Atom> rAtoms;
			copy_if(atoms.begin(), atoms.end(), back_inserter(rAtoms),
				[&](auto& a) { return a.labelAsymID() == asymID; });
//			for (auto a: db["atom_site"].find(cif::Key("label_asym_id") == asymID))
//				rAtoms.push_back(p.getAtomByID(a["id"].as<std::string>()));
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (compound->atomsBonded(rAtoms[i].labelAtomID(), rAtoms[j].labelAtomID()))
					{
						uint32_t ixa = index[rAtoms[i].id()];
						uint32_t ixb = index[rAtoms[j].id()];
						
						bond.insert(key(ixa, ixb));
					}
				}
			}
		}
	}
	
	
	// start by creating an index for single bonds

//std::cout << "Maken van b1_2 voor " << bond.size() << " bindingen" << std::endl;
	
	std::multimap<uint32_t,uint32_t> b1_2;
	for (auto& bk: bond)
	{
		uint32_t a, b;
		std::tie(a, b) = dekey(bk);
		
		b1_2.insert({ a, b });
		b1_2.insert({ b, a });
	}
	
//std::cout << "Afmeting b1_2: " << b1_2.size() << std::endl;

	std::multimap<uint32_t,uint32_t> b1_3;
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
				
				b1_3.insert({ x, y });
				b1_3.insert({ y, x });
			}
		}
	}

//std::cout << "Afmeting b1_3: " << b1_3.size() << std::endl;

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
//std::cout << "Afmeting b1_4: " << bond_1_4.size() << std::endl;
}

std::vector<std::string> BondMap::linked(const Atom& a) const
{
	auto i = link.find(a.id());
	
	std::vector<std::string> result;
	
	if (i != link.end())
		result = std::vector<std::string>(i->second.begin(), i->second.end());

	return result;
}


}
