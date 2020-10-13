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
   Date: maandag 19 februari, 2018
*/

#include "pdb-redo.hpp"

#include <iostream>
#include <fstream>
#include <regex>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <ctgmath>

#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>

#include <zeep/xml/document.hpp>
#include <zeep/xml/serialize.hpp>
// #include <zeep/xml/writer.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Structure.hpp"
#include "cif++/Compound.hpp"
#include "cif++/CifUtils.hpp"

#include "AtomShape.hpp"
#include "BondMap.hpp"
#include "MapMaker.hpp"
#include "ResolutionCalculator.hpp"

#include "HBondTraits.h"

#include "svm++.h"

using namespace clipper;
using data32::Flag;
using data32::F_phi;
using data32::F_sigF;
using data32::Phi_fom;
using mmcif::kPI;

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace zx = zeep::xml;
namespace c = mmcif;

// --------------------------------------------------------------------

class ChironHelper
{
  public:
	ChironHelper(const mmcif::Compound& compound, const std::string& centre)
		: mCompound(compound), mCentre(centre)
		, mAtoms(mCompound.atoms()), mBonds(mCompound.bonds())
	{
	}
	
	bool canSwapChains(std::vector<std::tuple<std::string,std::string>>& swapAtoms);
	bool canSwapIsomers(const mmcif::Residue& res, std::string& swapCompound, std::vector<std::tuple<std::string,std::string>>& swapAtoms);
	bool canPullThroughPlane(const mmcif::Residue& r, mmcif::Point& newCoordinates);
	
  private:

	bool canSwapChains(const std::string& atom1, const std::string& atom2,
		std::vector<std::tuple<std::string,std::string>>& swapAtoms);

	bool canSwapSubChain(std::string atom1, std::string atom2,
		std::set<std::string> visited1, std::set<std::string> visited2,
		std::vector<std::tuple<std::string,std::string>>& swapAtoms);

	const mmcif::CompoundAtom& atomForID(const std::string& id) const
	{
		auto i = find_if(mAtoms.begin(), mAtoms.end(),
			[&](const mmcif::CompoundAtom& a) { return a.id == id; });
		if (i == mAtoms.end())
			throw std::runtime_error("Could not find atom " + id + " in compound " + mCompound.id());
		return *i;
	}

	size_t countChiralErrors(const mmcif::Residue& res, const mmcif::Compound& c) const
	{
		std::vector<std::tuple<std::string,std::string>> remapped;
		return countChiralErrors(res, c, remapped);
	}
	
	size_t countChiralErrors(const mmcif::Residue& res, const mmcif::Compound& c, const std::vector<std::tuple<std::string,std::string>>& remapped) const;

	const mmcif::Compound&			mCompound;
	std::string							mCentre;
	std::vector<mmcif::CompoundAtom>	mAtoms;
	std::vector<mmcif::CompoundBond>	mBonds;
};

bool ChironHelper::canPullThroughPlane(const mmcif::Residue& r, mmcif::Point& newCoordinates)
{
	bool result = false;

	// see if one of the bonds is to a hydrogen
	std::vector<std::string> l;
	for (auto& b: mBonds)
	{
		std::string linked;
		
		if (b.atomID[0] == mCentre)
			linked = b.atomID[1];
		else if (b.atomID[1] == mCentre)
			linked = b.atomID[0];
		
		if (not linked.empty() and atomForID(linked).typeSymbol != mmcif::H)
			l.push_back(linked);
	}
	
	// Only when we have three atoms other than H bonded we can pull
	// the centre through the plane
	if (l.size() == 3)
	{
		result = true;
		
		// Construct plane based on three points
		// a.x + b.y + c.z + d = 0

		auto centre = r.atomByID(mCentre).location();
		
		// Normal std::vector on the plane is:
		auto p1 = r.atomByID(l[0]).location();
		auto p2 = r.atomByID(l[1]).location();
		auto p3 = r.atomByID(l[2]).location();
		
		auto n = CrossProduct(p2 - p1, p3 - p1);
		
		float a = n.getX();
		float b = n.getY();
		float c = n.getZ();
		float d = p1.getX() * a + p1.getY() * b + p1.getZ() * c;
		
		d -= a * centre.getX() + b * centre.getY() + c * centre.getZ();
		
		// point closest to centre for this plane:
		mmcif::Point o {
			a * d / (a * a + b * b + c * c),
			b * d / (a * a + b * b + c * c),
			c * d / (a * a + b * b + c * c)
		};
		o += centre;
		
		auto l = o - centre;
		
		newCoordinates = o + l / 100.f;
	}
	
	return result;
}

bool ChironHelper::canSwapChains(std::vector<std::tuple<std::string,std::string>>& swapAtoms)
{
	std::vector<std::string> l;
	
	for (auto& b: mBonds)
	{
		if (b.atomID[0] == mCentre)
			l.push_back(b.atomID[1]);
		else if (b.atomID[1] == mCentre)
			l.push_back(b.atomID[0]);
	}
	
	bool result = false;

	for (size_t i = 0; i + 1 < l.size(); ++i)
	{
		for (size_t j = i + 1; j < l.size(); ++j)
		{
			if (canSwapChains(l[i], l[j], swapAtoms))
			{
				result = true;
				break;
			}
		}
	}

	return result;
}

bool ChironHelper::canSwapChains(const std::string& atom1, const std::string& atom2,
		std::vector<std::tuple<std::string,std::string>>& swapAtoms)
{
	std::set<std::string> v1, v2;
	
	v1.insert(mCentre);
	v2.insert(mCentre);
	
	return canSwapSubChain(atom1, atom2, v1, v2, swapAtoms);
}

bool ChironHelper::canSwapIsomers(const mmcif::Residue& res, std::string& swapCompound, std::vector<std::tuple<std::string,std::string>>& swapAtoms)
{
	bool result = false;
	
	if (cif::VERBOSE or not mCompound.isSugar())
	{
		auto isomers = mCompound.isomers();
	
		if (not isomers.empty())
		{
			auto chiralErrInitial = countChiralErrors(res, mCompound);
			auto chiralErr = chiralErrInitial;
			
			if (cif::VERBOSE > 1)
				std::cerr << "Trying to swap isomers, initial error count is " << chiralErrInitial << std::endl;
			
			for (auto i: isomers)
			{
				auto c = mmcif::Compound::create(i);
				
				std::vector<std::tuple<std::string,std::string>> m = c->mapToIsomer(mCompound);
				
				if (cif::VERBOSE > 2)
				{
					for (auto& a: m)
						std::cerr << "  " << std::get<0>(a) << " => " << std::get<1>(a) << std::endl;
				}
				
				auto err = countChiralErrors(res, *c, m);
	
				if (cif::VERBOSE > 1)
					std::cerr << "isomer " << i << " has " << err << " errors" << std::endl;
	
				if (chiralErr > err)
				{
					chiralErr = err;
					result = true;
					swapCompound = i;
					swap(swapAtoms, m);
					
					if (cif::VERBOSE > 1)
						std::cerr << "Err count decreased to " << chiralErr << std::endl; 
	
					if (chiralErr == 0)	 // we're done.
						break;
				}
			}
	
			if (result and mCompound.isSugar())
			{
				if (cif::VERBOSE)
					std::cerr << "Since residue " << res.compoundID() << " is a sugar, it will not be swapped with " << swapCompound << std::endl;
				result = false;
			}
		}
	}
	
	return result;
}

size_t ChironHelper::countChiralErrors(const mmcif::Residue& res, const mmcif::Compound& c, const std::vector<std::tuple<std::string,std::string>>& mapping) const
{
	size_t result = 0;
	
	for (auto cc: c.chiralCentres())
	{
		if (cc.volumeSign == mmcif::both)
			continue;
		
		auto rename = [&](const std::string& name) -> std::string
		{
			std::string result;
			
			auto i = find_if(mapping.begin(), mapping.end(), [&](auto& m) { return std::get<0>(m) == name; });
			if (i == mapping.end())
			{
//				if (cif::VERBOSE > 1 and not mapping.empty())
//					std::cerr << "no mapping found for atom " << name << " in " << c.id() << std::endl;
				 result = name;
			}
			else
				result = std::get<1>(*i);

			return result;
		};
		
		try
		{
			auto centre = res.atomByID(rename(cc.atomIDCentre));
			auto atom1 = res.atomByID(rename(cc.atomID[0]));
			auto atom2 = res.atomByID(rename(cc.atomID[1]));
			auto atom3 = res.atomByID(rename(cc.atomID[2]));

			auto chiralVolume = DotProduct(atom1.location() - centre.location(),
				CrossProduct(atom2.location() - centre.location(), atom3.location() - centre.location()));

			if ((chiralVolume < 0 and cc.volumeSign == mmcif::positiv) or
				(chiralVolume > 0 and cc.volumeSign == mmcif::negativ))
			{
				if (cif::VERBOSE > 1)
					std::cerr << "chiral error in " << c.id() << " around " << cc.atomIDCentre << " with volume: " << chiralVolume << std::endl;
				
				++result;
			}
		}
		catch (const std::exception& ex)
		{
			if (cif::VERBOSE)
				std::cerr << "Missing atom in counting chiral errors: " << ex.what() << std::endl;
		}
	}
	
	return result;
}

bool ChironHelper::canSwapSubChain(std::string atom1, std::string atom2,
		std::set<std::string> visited1, std::set<std::string> visited2,
		std::vector<std::tuple<std::string,std::string>>& swapAtoms)
{
	auto& a1 = atomForID(atom1);
	auto& a2 = atomForID(atom2);
	
	bool result = false;
	for (;;)
	{
		if (a1.typeSymbol != a2.typeSymbol)
			break;
		
		if (a1.typeSymbol == mmcif::H)
		{
			result = true;
			break;
		}

		std::vector<std::pair<std::string,mmcif::BondType>> l1, l2;
		
		visited1.insert(atom1);
		visited2.insert(atom2);
		
		for (auto& b: mBonds)
		{
			if (atom1 == b.atomID[0] and not visited1.count(b.atomID[1]))
				l1.push_back(std::make_pair(b.atomID[1], b.type));
			else if (atom1 == b.atomID[1] and not visited1.count(b.atomID[0]))
				l1.push_back(std::make_pair(b.atomID[0], b.type));

			if (atom2 == b.atomID[0] and not visited2.count(b.atomID[1]))
				l2.push_back(std::make_pair(b.atomID[1], b.type));
			else if (atom2 == b.atomID[1] and not visited2.count(b.atomID[0]))
				l2.push_back(std::make_pair(b.atomID[0], b.type));
		}
		
		if (l1.size() != l2.size())
			break;
		
		std::vector<std::tuple<std::string,std::string>> subSwap;		
		
		auto test = [&](int a, int b) -> bool
		{
			bool result = l1[a].second == l2[b].second and
				(l1[a].first == l2[a].first or canSwapSubChain(l1[a].first, l2[b].first, visited1, visited2, subSwap));
			if (not result)
				subSwap.clear();
			return result;
		};
		
		auto add = [&](std::initializer_list<std::pair<int,int>> l)
		{
			if (atom1 != atom2)
				swapAtoms.emplace_back(atom1, atom2);
			swapAtoms.insert(swapAtoms.end(), subSwap.begin(), subSwap.end());
			result = true;
		};
		
		switch (l1.size())
		{
			case 0:
				add({});
				break;
			
			case 1:
				if (test(0, 0))
					add({ std::make_pair(0, 0) });
				break;
			
			case 2:
				if (test(0, 0) and test(1, 1))
					add({ std::make_pair(0, 0), std::make_pair(1, 1) });
				else if (test(0, 1) and test(1, 0))
					add({ std::make_pair(0, 1), std::make_pair(1, 0) });
				break;
			
			case 3:
				if (test(0, 0))
				{
					if (test(1, 1) and test(2, 2)) 			add({ std::make_pair(0, 0), std::make_pair(1, 1), std::make_pair(2, 2) });
					else if (test(1, 2) and test(2, 1))		add({ std::make_pair(0, 0), std::make_pair(1, 2), std::make_pair(2, 1) });
				}
				else if (test(0, 1))
				{
					if (test(1, 2) and test(2, 0))			add({ std::make_pair(0, 1), std::make_pair(1, 2), std::make_pair(2, 0) });
					else if (test(1, 0) and test(2, 2))		add({ std::make_pair(0, 1), std::make_pair(1, 0), std::make_pair(2, 2) });
				}
				else if (test(0, 2))
				{
					if (test(1, 0) and test(2, 1))			add({ std::make_pair(0, 2), std::make_pair(1, 0), std::make_pair(2, 1) });
					else if (test(1, 1) and test(2, 0))		add({ std::make_pair(0, 2), std::make_pair(1, 1), std::make_pair(2, 0) });
				}
	
				break;
			
			default:
				throw std::runtime_error("unimplemented number of bonds");
		}

		break;
	}

	if (cif::VERBOSE > 2)
		std::cerr << "canSwap(" << atom1 << ", " << atom2 << ") => " << std::boolalpha << result << std::endl;

	return result;
}

// --------------------------------------------------------------------

int Process(mmcif::Structure& structure, const mmcif::Residue& res, const mmcif::Compound& compound)
{
	int result = 0;
	
	if (cif::VERBOSE > 1)
		std::cerr << "Process " << res.compoundID() << " " << res.asymID() << res.seqID() << std::endl;
	
	for (auto cc: compound.chiralCentres())
	{
		if (cc.volumeSign == mmcif::both)
			continue;
		
		try
		{
			auto centre = res.atomByID(cc.atomIDCentre);
			auto atom1 = res.atomByID(cc.atomID[0]);
			auto atom2 = res.atomByID(cc.atomID[1]);
			auto atom3 = res.atomByID(cc.atomID[2]);

			auto chiralVolume = DotProduct(atom1.location() - centre.location(),
				CrossProduct(atom2.location() - centre.location(), atom3.location() - centre.location()));

			if (cif::VERBOSE)
			{
				std::cerr << "chiral volume for " << res.compoundID() << " " << res.asymID() << res.seqID()
					 << " with centre " << cc.atomIDCentre
					 << " is " << chiralVolume << " and should be "
					 << (cc.volumeSign == mmcif::positiv ? "positive" : "negative") << std::endl;
			}

			if ((chiralVolume < 0 and cc.volumeSign == mmcif::positiv) or
				(chiralVolume > 0 and cc.volumeSign == mmcif::negativ))
			{
				std::cerr << "Error in chiral volume for " << res.compoundID() << " " << res.asymID() << res.seqID()
					 << " with centre " << cc.atomIDCentre;
				
				if (chiralVolume < 0)
					std::cerr << " volume should be positive but is negative: " << chiralVolume << std::endl;
				else
					std::cerr << " volume should be negative but is positive: " << chiralVolume << std::endl;

				ChironHelper test(compound, cc.atomIDCentre);
				
				std::vector<std::tuple<std::string,std::string>> swapAtoms;
				std::string swapCompound;
				mmcif::Point newCoordinates;
				
				if (test.canSwapChains(swapAtoms))
				{
					++result;
					
					std::cerr << "Flipping labels: ";
					
					for (auto p: swapAtoms)
					{
						std::string na1, na2;
						std::tie(na1, na2) = p;
						
						std::cerr << "{ " << na1 << " and " << na2 << " }, ";

						auto a1 = res.atomByID(na1);
						auto a2 = res.atomByID(na2);

						structure.swapAtoms(a1, a2);
					}

					std::cerr << std::endl;
				}
				else if (test.canSwapIsomers(res, swapCompound, swapAtoms))
				{
					++result;
					std::cerr << "Replacing with isomer " << swapCompound;

					if (swapAtoms.empty())
						std::cerr << std::endl;
					else
					{
						std::cerr << ", swapping atom labels: ";

						for (auto p: swapAtoms)
						{
							std::string na1, na2;
							std::tie(na1, na2) = p;
							
							std::cerr << "{ " << na1 << " -> " << na2 << " }, ";
						}

						std::cerr << std::endl;
					}
					
					structure.changeResidue(res, swapCompound, swapAtoms);
				}
				else if (test.canPullThroughPlane(res, newCoordinates))
				{
					++result;

					structure.moveAtom(centre, newCoordinates);
					
					chiralVolume = DotProduct(atom1.location() - newCoordinates,
						CrossProduct(atom2.location() - newCoordinates, atom3.location() - newCoordinates));
		
					std::cerr << "Pulling it through the plane from " << centre.location() << " to " << newCoordinates << std::endl
						 << "  new chiral volume is " << chiralVolume << std::endl;
				}
				else
					std::cerr << "Cannot fix this problem" << std::endl;
			}
		}
		catch (const std::runtime_error& ex)
		{
			if (cif::VERBOSE)
				std::cerr << ex.what() << std::endl;
			
			continue;
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

int Process(mmcif::Structure& structure)
{
	int result = 0;
	
	for (auto& res: structure.nonPolymers())
	{
		// skip waters...
		if (res.isWater())
			continue;
		
		auto& compound = res.compound();

		if (compound.chiralCentres().empty())
			continue;
		
		result += Process(structure, res, compound);
	}

	for (auto& poly: structure.polymers())
	{
		for (auto& m: poly)
		{
			auto& compound = m.compound();

			if (compound.chiralCentres().empty())
				continue;
			
			result += Process(structure, m, compound);
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("xyzin",				po::value<std::string>(),	"coordinates file")
		("output,o",			po::value<std::string>(),	"Write output to this file instead of stdout")
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("xyzin") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}
	
	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	mmcif::File f(vm["xyzin"].as<std::string>());
	mmcif::Structure structure(f);

	int result = Process(structure);

	if (result)
	{
		if (vm.count("output"))
			structure.getFile().save(vm["output"].as<std::string>());
		else
			structure.getFile().file().save(std::cout);
	}
	
	return result;
}
