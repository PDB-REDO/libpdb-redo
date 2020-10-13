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

/* 
   Created by: Maarten L. Hekkelman
   Date: maandag 07 januari, 2019
*/

#include "cif++/Config.hpp"

#include <fstream>
#include <filesystem>

#include <zeep/xml/serialize.hpp>
#include <zeep/xml/document.hpp>

#include "skiplist.h"

using namespace std;
namespace fs = std::filesystem;
namespace zx = zeep::xml;

// --------------------------------------------------------------------

struct SkipResidueLabel
{
	std::string	asymID;
	int			seqID;
	std::string	monID;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		ar & zx::make_attribute_nvp("asym-id", asymID)
		   & zx::make_attribute_nvp("seq-id", seqID)
		   & zx::make_attribute_nvp("mon-id", monID);
	}
};

// --------------------------------------------------------------------

struct SkipResiduePDB
{
	std::string	strandID;
	int			seqNum;
	std::string	monID;
	std::string	insCode;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		ar & zx::make_attribute_nvp("strand-id", strandID)
		   & zx::make_attribute_nvp("seq-num", seqNum)
		   & zx::make_attribute_nvp("mon-id", monID)
		   & zx::make_attribute_nvp("ins-code", insCode);
	}
};

//struct SkipResidueAuth
//{
//	std::string	asymID;
//	int			seqID;
//	std::string	compID;
//
//	template<class Archive>
//	void serialize(Archive& ar, const unsigned int version)
//	{
//		ar & zx::make_attribute_nvp("asym-id", asymID)
//		   & zx::make_attribute_nvp("seq-id", seqID)
//		   & zx::make_attribute_nvp("comp-id", compID);
//	}
//};

// --------------------------------------------------------------------

template<typename T>
struct SkipResidueList
{
	vector<T>	res;
	
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		ar & zx::make_element_nvp("residue", res);
	}
};

// --------------------------------------------------------------------

struct SkipListContainer
{
	SkipListNumberingScheme scheme;
	SkipResidueList<SkipResidueLabel> label;
	SkipResidueList<SkipResiduePDB> pdb;
//	vector<SkipResidueAuth> auth;

	template<class Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		ar & zx::make_attribute_nvp("scheme", scheme);
		
		switch (scheme)
		{
			case sl_Label:
				ar & zx::make_element_nvp("residues", label);
				break;
			
			case sl_PDB:
				ar & zx::make_element_nvp("residues", pdb);
				break;
				
			default:
				throw runtime_error("Unimplemented skiplist scheme");
		}
	}
};

// --------------------------------------------------------------------

class SkipListInitializer
{
  public:
	SkipListInitializer();
};

SkipListInitializer::SkipListInitializer()
{
	zeep::value_serializer<SkipListNumberingScheme>::init("SkipListNumberingScheme", {
		{ sl_Label, "label" },
		{ sl_PDB, "pdb" }
	});
}

SkipListInitializer sInitSkipListTypes;

// --------------------------------------------------------------------

SkipList readSkipList(const fs::path& file, const mmcif::Structure& structure)
{
	std::ifstream is(file);
	if (not is.is_open())
		throw runtime_error("Could not open skip list file " + file.string());
	
	zx::document doc(is);
	
	SkipListContainer list;
	doc.deserialize("skip-list", list);
	
	SkipList result;
	
	switch (list.scheme)
	{
		case sl_Label:
			for (auto& r: list.label.res)
				result.push_back({r.asymID, r.seqID, r.monID});
			break;
		
		case sl_PDB:
			for (auto& r: list.pdb.res)
			{
				SkipResidue lr;
				tie(lr.asymID, lr.seqID, lr.monID) =
					structure.MapPDBToLabel(r.strandID, r.seqNum, r.monID, r.insCode);
				
				if (lr.asymID.empty())
					throw runtime_error("Could not map residue specified in skip list (" +
						r.monID + ' ' + r.strandID + to_string(r.seqNum) + r.insCode + ')');
				
				result.emplace_back(move(lr));
			}
			break;
		
		default:
			throw runtime_error("Unknown skip list scheme");
	}
	
	return result;
}

void writeSkipList(const fs::path& file, const mmcif::Structure& structure,
	const SkipList& list, SkipListNumberingScheme scheme)
{
	
}
