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

#include "pdb-redo.hpp"

#include <set>
#include <regex>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "cif++/Cif++.hpp"

#include "HBondTraits.h"

namespace fs = std::filesystem;
namespace ba = boost::algorithm;

NotAHBondSet* NotAHBondSet::Create()
{
	const char* clibdMon = getenv("CLIBD_MON");
	if (clibdMon == nullptr)
		throw std::runtime_error("Cannot locate peptide list, please source the CCP4 environment");
	
	fs::path dir(clibdMon);
	
	const std::regex kNotAHBondEnergyTypeRX("NR(?:5|55|56|6|66)");
	NotAHBondAtomSet atomSet;
	
	for (auto d = fs::directory_iterator(dir); d != fs::directory_iterator(); ++d)
	{
		if (d->path().filename().string().length() != 1)
			continue;
		
		for (auto e = fs::directory_iterator(d->path()); e != fs::directory_iterator(); ++e)
		{
			if (e->path().extension().string() != ".cif")
				continue;
			
			std::ifstream file(e->path());
			if (not file.is_open())
				throw std::runtime_error("Could not open file: " + e->path().string());
			
			cif::File data(file);
			
			for (auto& compName: data["comp_list"]["chem_comp"])
			{
				std::string compId = compName["id"].as<std::string>();
				
				try
				{
					for (auto& compAtom: data["comp_" + compId]["chem_comp_atom"])
					{
						std::string atomType = compAtom["type_symbol"].as<std::string>();
						std::string typeEnergy = compAtom["type_energy"].as<std::string>();
					
						if (atomType == "N" and (typeEnergy == "N" or regex_match(typeEnergy, kNotAHBondEnergyTypeRX)))
							atomSet.insert({ compId, compAtom["atom_id"].as<std::string>() });
					}
				}
				catch (const std::exception& ex)
				{
					std::cerr << "Error reading " << e->path() << ": " << ex.what() << std::endl;
				}
			}
		}
	}
	
	return new NotAHBondSet(move(atomSet));
}

void NotAHBondSet::Save(zeep::xml::container& inNode)
{
	using zeep::xml::element;
	
	auto setNode = new element("not-a-hbond");
	
	for (auto& atom: mData)
	{
		auto atomNode = new element("atom");
		atomNode->set_attribute("monomer", atom.monomer);
		atomNode->set_attribute("atom_id", atom.atomId);
		setNode->append(atomNode);
	}
	
	inNode.append(setNode);
}

NotAHBondSet* NotAHBondSet::Load(zeep::xml::container& inNode)
{
	NotAHBondAtomSet data;
	
	auto setNode = inNode.find_first("not-a-hbond");
	if (setNode != nullptr)
	{
		for (auto atom: *setNode)
		{
			if (atom->name() != "atom")
				continue;
			data.insert({ atom->get_attribute("monomer"), atom->get_attribute("atom_id") });
		}
	}
	
	return new NotAHBondSet(move(data));
}
