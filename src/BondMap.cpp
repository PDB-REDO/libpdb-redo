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

#include <fstream>
#include <algorithm>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/CifUtils.hpp"

#include "pdb-redo/BondMap.hpp"
#if USE_RSRC
#include "mrsrc.hpp"
#endif

#ifndef DATADIR
#define DATADIR "/rsrc/"
#endif

namespace fs = std::filesystem;
namespace io = boost::iostreams;

namespace mmcif
{

// --------------------------------------------------------------------

uint32_t to_id(const std::string& id)
{
	assert(id.length() <= 4);

	uint32_t result = 0;

	for (uint8_t ch: id)
		result = result << 8 | ch;

	return result;
}

std::string from_id(uint32_t id)
{
	std::string result;

	while (id)
	{
		result.insert(result.begin(), static_cast<char>(id & 0x0ff));
		id >>= 8;
	}

	return result;
}

// --------------------------------------------------------------------

struct CompoundBondInfoFileHeader
{
	char		signature[4] = { 'C', 'B', 'I', 'H'  };
	uint32_t	headerSize = sizeof(CompoundBondInfoFileHeader);
	uint32_t	indexEntries;
	uint32_t	atomEntries;
	uint64_t	dataSize;
};

struct CompoundBondInfo
{
	uint32_t	id;
	uint32_t	nrOfAtoms;
	int64_t		dataOffset;
};

void createBondInfoFile(const fs::path& components, const fs::path& infofile)
{
	std::ofstream outfile(infofile, std::ios::binary);
	if (not outfile.is_open())
		throw BondMapException("Could not create bond info file " + infofile.string());

	cif::File infile(components);

	std::set<uint32_t> atomIDs;
	std::vector<uint32_t> compoundIDs;

	for (auto& db: infile)
	{
		auto chem_comp_bond = db.get("chem_comp_bond");
		if (not chem_comp_bond)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "Missing chem_comp_bond category in data block " << db.getName() << std::endl;
			continue;
		}

		for (const auto& [atom_id_1, atom_id_2]: chem_comp_bond->rows<std::string,std::string>({"atom_id_1", "atom_id_2"}))
		{
			atomIDs.insert(to_id(atom_id_1));
			atomIDs.insert(to_id(atom_id_2));
		}

		compoundIDs.push_back(to_id(db.getName()));
	}

	if (cif::VERBOSE)
		std::cout << "Number of unique atom names is " << atomIDs.size() << std::endl
				  << "Number of unique residue names is " << compoundIDs.size() << std::endl;

	CompoundBondInfoFileHeader header = {};
	header.indexEntries = compoundIDs.size();
	header.atomEntries = atomIDs.size();

	outfile.write(reinterpret_cast<char*>(&header), sizeof(header));
	
	for (auto atomID: atomIDs)
		outfile.write(reinterpret_cast<char*>(&atomID), sizeof(uint32_t));

	auto dataOffset = outfile.tellp();

	std::vector<CompoundBondInfo> entries;
	entries.reserve(compoundIDs.size());

	std::map<uint32_t, uint16_t> atomIDMap;
	for (auto& atomID: atomIDs)
		atomIDMap[atomID] = atomIDMap.size();

	for (auto& db: infile)
	{
		auto chem_comp_bond = db.get("chem_comp_bond");
		if (not chem_comp_bond)
			continue;

		std::set<uint16_t> bondedAtoms;

		for (const auto& [atom_id_1, atom_id_2]: chem_comp_bond->rows<std::string,std::string>({"atom_id_1", "atom_id_2"}))
		{
			bondedAtoms.insert(atomIDMap[to_id(atom_id_1)]);
			bondedAtoms.insert(atomIDMap[to_id(atom_id_2)]);
		}

		std::map<uint16_t, int32_t> bondedAtomMap;
		for (auto id: bondedAtoms)
			bondedAtomMap[id] = static_cast<int32_t>(bondedAtomMap.size());
		
		CompoundBondInfo info = {
			to_id(db.getName()),
			static_cast<uint32_t>(bondedAtomMap.size()),
			outfile.tellp() - dataOffset
		};

		entries.push_back(info);

		// An now first write the array of atom ID's in this compound
		std::vector<uint16_t> atoms(bondedAtoms.begin(), bondedAtoms.end());
		outfile.write(reinterpret_cast<char*>(atoms.data()), sizeof(uint16_t) * atoms.size());

		// And then the symmetric matrix with bonds
		size_t N = atoms.size();
		size_t M = (N * (N - 1)) / 2;

		size_t K = M / 8;
		if (M % 8)
			K += 1;
		
		std::vector<uint8_t> m(K);

		for (const auto& [atom_id_1, atom_id_2]: chem_comp_bond->rows<std::string,std::string>({"atom_id_1", "atom_id_2"}))
		{
			auto a = bondedAtomMap[atomIDMap[to_id(atom_id_1)]];
			auto b = bondedAtomMap[atomIDMap[to_id(atom_id_2)]];

			assert(a != b);
			assert((int)b < (int)N);

			if (a > b)
				std::swap(a, b);
			
			size_t ix = ((b - 1) * b) / 2 + a;
			assert(ix < M);

			auto Bix = ix / 8;
			auto bix = ix % 8;

			m[Bix] |= 1 << bix;
		}

		outfile.write(reinterpret_cast<char*>(m.data()), m.size());
	}

	header.dataSize = outfile.tellp() - dataOffset;

	std::sort(entries.begin(), entries.end(), [](CompoundBondInfo& a, CompoundBondInfo& b)
	{
		return a.id < b.id;
	});

	outfile.write(reinterpret_cast<char*>(entries.data()), sizeof(CompoundBondInfo) * entries.size());

	outfile.seekp(0);
	outfile.write(reinterpret_cast<char*>(&header), sizeof(header));
}

// --------------------------------------------------------------------

class CompoundBondMap
{
  public:

	~CompoundBondMap()
	{
		delete[] m_data;
	}

	static CompoundBondMap& instance()
	{
		if (not s_instance)
			s_instance.reset(new CompoundBondMap());
		return *s_instance;
	}

	bool bonded(const std::string& compoundID, const std::string& atomID1, const std::string& atomID2) const
	{
		return bonded(compound_index_nr(compoundID), atomID1, atomID2);
	}

	bool bonded(int32_t compoundNr, const std::string& atomID1, const std::string& atomID2) const;

	int32_t compound_index_nr(const std::string& compoundID) const
	{
		int32_t L = 0, R = static_cast<int32_t>(m_header.indexEntries) - 1;
		uint32_t compID = to_id(compoundID);

		while (L <= R)
		{
			auto i = (L + R) / 2;

			if (m_compounds[i].id <= compID)
				L = i + 1;
			else
				R = i - 1;
		}

		return R >= 0 and m_compounds[R].id == compID ? R : -1;
	}

	std::vector<std::string> atomIDsForCompound(const std::string& compoundID) const;

  private:

	static std::unique_ptr<CompoundBondMap> s_instance;

	CompoundBondMap();

	void init(std::istream&& is);

	int32_t atom_nr(const std::string& atomID) const
	{
		auto id = to_id(atomID);

		int L = 0, R = m_atom_ids.size() - 1;
		while (L <= R)
		{
			int i = (L + R) / 2;
			
			if (m_atom_ids[i] <= id)
				L = i + 1;
			else
				R = i - 1;
		}

		return R >= 0 and m_atom_ids[R] == id ? R : -1;
	}

	CompoundBondInfoFileHeader m_header;
	std::vector<CompoundBondInfo> m_compounds;
	std::vector<uint32_t> m_atom_ids;
	const uint8_t* m_data = nullptr;
};

std::unique_ptr<CompoundBondMap> CompoundBondMap::s_instance;

CompoundBondMap::CompoundBondMap()
{
#if USE_RSRC
	mrsrc::rsrc rsrc("bond-info.bin");
	if (rsrc)
	{
		init(mrsrc::istream(rsrc));
		return;
	}
#endif

	fs::path bondInfoFile;

	if (getenv("BOND_INFO_FILE") != nullptr)
		bondInfoFile = fs::path(getenv("BOND_INFO_FILE"));
	else
		bondInfoFile = fs::path(DATADIR) / "bond-info.bin";

	if (fs::exists(bondInfoFile))
		init(std::ifstream(bondInfoFile, std::ios::binary));
	else
		throw BondMapException("Missing bond info file");
}

void CompoundBondMap::init(std::istream&& is)
{
	io::filtering_stream<io::input> in;
	in.push(io::gzip_decompressor());
	in.push(is);

	in.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));

	if (m_header.signature[0] != 'C' or
		m_header.signature[1] != 'B' or
		m_header.signature[2] != 'I' or
		m_header.signature[3] != 'H' or
		m_header.headerSize != sizeof(CompoundBondInfoFileHeader))
		throw BondMapException("Invalid bond info file");

	m_atom_ids.resize(m_header.atomEntries);
	in.read(reinterpret_cast<char*>(m_atom_ids.data()), m_header.atomEntries * sizeof(uint32_t));

	m_data = new uint8_t[m_header.dataSize];
	in.read(reinterpret_cast<char*>(const_cast<uint8_t*>(m_data)), m_header.dataSize);

	m_compounds.resize(m_header.indexEntries);
	in.read(reinterpret_cast<char*>(m_compounds.data()), m_header.indexEntries * sizeof(CompoundBondInfo));
}

bool CompoundBondMap::bonded(int32_t compoundNr, const std::string& atomID1, const std::string& atomID2) const
{
	bool result = false;

	if (compoundNr >= 0 and compoundNr < static_cast<int32_t>(m_header.indexEntries))
	{
		auto a1 = atom_nr(atomID1);		if (a1 < 0) throw BondMapException("Unknown atom ID " + atomID1);
		auto a2 = atom_nr(atomID2);		if (a2 < 0) throw BondMapException("Unknown atom ID " + atomID2);

		auto& comp = m_compounds[compoundNr];
		auto comp_data = m_data + comp.dataOffset;

		// memory might be unaligned (yes, I'm old...)
		std::vector<uint16_t> comp_atoms(comp.nrOfAtoms);
		std::memcpy(comp_atoms.data(), comp_data, comp.nrOfAtoms * sizeof(uint16_t));

		auto matrix = comp_data + comp.nrOfAtoms * sizeof(uint16_t);

		auto i = std::find(comp_atoms.begin(), comp_atoms.end(), a1);
		if (i == comp_atoms.end())
			throw BondMapException("Compound " + from_id(compoundNr) + " does not have an atom named " + atomID1);
		auto ix_1 = i - comp_atoms.begin();

		i = std::find(comp_atoms.begin(), comp_atoms.end(), a2);
		if (i == comp_atoms.end())
			throw BondMapException("Compound " + from_id(compoundNr) + " does not have an atom named " + atomID2);
		auto ix_2 = i - comp_atoms.begin();

		if (ix_1 > ix_2)
			std::swap(ix_1, ix_2);
		

		size_t ix = ((ix_2 - 1) * ix_2) / 2 + ix_1;

		auto Bix = ix / 8;
		auto bix = ix % 8;

		result = matrix[Bix] & (1 << bix);
	}

	return result;
}

std::vector<std::string> CompoundBondMap::atomIDsForCompound(const std::string& compoundID) const
{
	std::vector<std::string> result;

	auto comp_nr = compound_index_nr(compoundID);
	{
		auto& comp = m_compounds[comp_nr];
		auto comp_data = m_data + comp.dataOffset;

		// memory might be unaligned (yes, I'm old...)
		std::vector<uint16_t> comp_atoms(comp.nrOfAtoms);
		std::memcpy(comp_atoms.data(), comp_data, comp.nrOfAtoms * sizeof(uint16_t));

		for (auto& compAtom: comp_atoms)
			result.push_back(from_id(m_atom_ids[compAtom]));
	}

	return result;
}

// --------------------------------------------------------------------

BondMap::BondMap(const Structure& p)
{
	auto& compoundBondInfo = CompoundBondMap::instance();

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
		if (c == "HOH" or c == "H2O" or c == "WAT")
		{
			if (cif::VERBOSE)
				std::cerr << "skipping water in bond map calculation" << std::endl;
			continue;
		}

		auto compound_nr = compoundBondInfo.compound_index_nr(c);
		auto* compound = compound_nr < 0 ? mmcif::Compound::create(c) : nullptr;

		if (compound_nr < 0 and compound == nullptr)
			std::cerr << "Missing bond information for compound " << c << std::endl;

		auto bonded = [compound_nr, compound, &compoundBondInfo](const Atom& a, const Atom& b)
		{
			auto label_a = a.labelAtomID();
			auto label_b = b.labelAtomID();

			return compound_nr >= 0
				? compoundBondInfo.bonded(compound_nr, label_a, label_b)
				: (compound ? compound->atomsBonded(label_a, label_b) : false);
		};

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
					if (bonded(rAtoms[i], rAtoms[j]))
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
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
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
			
			for (uint32_t i = 0; i + 1 < rAtoms.size(); ++i)
			{
				for (uint32_t j = i + 1; j < rAtoms.size(); ++j)
				{
					if (bonded(rAtoms[i], rAtoms[j]))
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
	
	std::multimap<uint32_t,uint32_t> b1_2;
	for (auto& bk: bond)
	{
		uint32_t a, b;
		std::tie(a, b) = dekey(bk);
		
		b1_2.insert({ a, b });
		b1_2.insert({ b, a });
	}
	
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

std::vector<std::string> BondMap::linked(const Atom& a) const
{
	auto i = link.find(a.id());
	
	std::vector<std::string> result;
	
	if (i != link.end())
		result = std::vector<std::string>(i->second.begin(), i->second.end());

	return result;
}

std::vector<std::string> BondMap::atomIDsForCompound(const std::string& compoundID)
{
	auto& compoundBondInfo = CompoundBondMap::instance();

	auto result = compoundBondInfo.atomIDsForCompound(compoundID);

	if (result.empty())
	{
		auto* compound = mmcif::Compound::create(compoundID);

		if (compound == nullptr)
			throw BondMapException("Missing bond information for compound " + compoundID);

		for (auto& compAtom: compound->atoms())
			result.push_back(compAtom.id);
	}

	return result;
}

}
