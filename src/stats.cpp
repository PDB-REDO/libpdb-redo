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
   Date: woensdag 27 december, 2017
*/

#include "pdb-redo.hpp"

#include <fstream>
#include <iomanip>
#include <filesystem>

#include <boost/program_options.hpp>

#include "zeep/json/element.hpp"

#include "BondMap.hpp"
#include "Statistics.hpp"


namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " [options] <mtzfile> <coordinatesfile> [<output>]");
	visible_options.add_options()
		("hklin",				po::value<std::string>(),	"mtz file")
		("recalc",										"Recalculate Fc from FP/SIGFP in mtz file")
		("aniso-scaling",		po::value<std::string>(),	"Anisotropic scaling (none/observed/calculated)")
		("no-bulk",										"No bulk correction")
		("xyzin",				po::value<std::string>(),	"coordinates file")
		("fomap",				po::value<std::string>(),	"Fo map file -- 2mFo - DFc")
		("dfmap",				po::value<std::string>(),	"difference map file -- 2(mFo - DFc)")
		("reshi",				po::value<float>(),		"High resolution")
		("reslo",				po::value<float>(),		"Low resolution")
		("sampling-rate",		po::value<float>(),		"Sampling rate")
		("electron-scattering",							"Use electron scattering factors")
		("no-edia",										"Skip EDIA score calculation")
		("output,o",			po::value<std::string>(),	"Write output to this file instead of stdout")
		("output-format",		po::value<std::string>(),	"Output format, can be either 'edstats' or 'json'")
		("use-auth-ids",								"Write auth_ identities instead of label_")
		("dict",				po::value<std::string>(),	"Dictionary file containing restraints for residues in this specific target")
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
	p.add("hklin", 1);
	p.add("xyzin", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

//	po::variables_map vm;
//	po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
//	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cerr << visible_options << std::endl;
		return 0;
	}
	
	if (vm.count("xyzin") == 0 or
		(vm.count("hklin") == 0 and (vm.count("fomap") == 0 or vm.count("dfmap") == 0)))
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}
	
	const std::set<std::string> kAnisoOptions{ "none", "calculated", "observed" };
	if (vm.count("aniso-scaling") and kAnisoOptions.count(vm["aniso-scaling"].as<std::string>()) == 0)
	{
		std::cerr << "Invalid option for aniso-scaling, allowed values are none, observed and calculated" << std::endl;
		exit(1);
	}
	
	if (vm.count("fomap") and (vm.count("reshi") == 0 or vm.count("reslo") == 0))
	{
		std::cerr << "The reshi and reslo parameters are required when using std::map files" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// Load dict, if any
	
	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<std::string>());

	mmcif::File f(vm["xyzin"].as<std::string>());
	mmcif::Structure structure(f);

	bool electronScattering = vm.count("electron-scattering") > 0;
	if (not electronScattering)
	{
		auto& exptl = f.data()["exptl"];
		electronScattering = not exptl.empty() and exptl.front()["method"] == "ELECTRON CRYSTALLOGRAPHY";
	}
	
	c::MapMaker<float> mm;
	
	if (vm.count("hklin"))
	{
		float samplingRate = 0.75;
		if (vm.count("sampling-rate"))
			samplingRate = vm["sampling-rate"].as<float>();
	
		if (vm.count("recalc"))
		{
			auto aniso = c::MapMaker<float>::as_None;
			if (vm.count("aniso-scaling"))
			{
				if (vm["aniso-scaling"].as<std::string>() == "observed")
					aniso = c::MapMaker<float>::as_Observed;
				else if (vm["aniso-scaling"].as<std::string>() == "calculated")
					aniso = c::MapMaker<float>::as_Calculated;
			}
			
			mm.calculate(
				vm["hklin"].as<std::string>(), structure, vm.count("no-bulk"), aniso, samplingRate, electronScattering);
		}
		else
			mm.loadMTZ(vm["hklin"].as<std::string>(), samplingRate);
	}
	else
	{
		float reshi = vm["reshi"].as<float>();
		float reslo = vm["reslo"].as<float>();
		
		mm.loadMaps(vm["fomap"].as<std::string>(), vm["dfmap"].as<std::string>(), reshi, reslo);
	}
	
	std::vector<mmcif::ResidueStatistics> r;
	
	if (vm.count("no-edia"))
	{
		mmcif::StatsCollector collector(mm, structure, electronScattering);
		r = collector.collect();
	}
	else
	{
		mmcif::BondMap bm(structure);

		mmcif::EDIAStatsCollector collector(mm, structure, electronScattering, bm);
		r = collector.collect();
	}

	std::ofstream of;
	if (vm.count("output"))
	{
		of.open(vm["output"].as<std::string>());
		if (not of.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}
		std::cout.rdbuf(of.rdbuf());
	}

	if (vm.count("output-format") and vm["output-format"].as<std::string>() == "json")
	{
		using object = zeep::json::element;
		
		// std::vector<object> rs;
		
		// for (auto i: r)
		// {
		// 	object res;
		// 	res["asymID"] = i.asymID;
		// 	res["seqID"] = i.seqID;
		// 	res["compID"] = i.compID;
			
		// 	std::tuple<std::string,int,std::string,std::string> pdbID = structure.MapLabelToPDB(i.asymID, i.seqID, i.compID, i.authSeqID);
			
		// 	object pdb;
		// 	pdb["strandID"] = get<0>(pdbID);
		// 	pdb["seqNum"] = get<1>(pdbID);
		// 	pdb["compID"] = get<2>(pdbID);
		// 	pdb["insCode"] = get<3>(pdbID);
		// 	res["pdb"] = pdb;
			
		// 	res["RSR"] = i.RSR;
		// 	res["SRSR"] = i.SRSR;
		// 	res["RSCCS"] = i.RSCCS;
		// 	res["NGRID"] = i.ngrid;
		// 	res["EDIAm"] = i.EDIAm;
		// 	res["OPIA"] = i.OPIA;
			
		// 	rs.push_back(move(res));
		// }
		
		// object stats(rs);

		object stats;
		
		for (auto i: r)
		{
			std::tuple<std::string,int,std::string,std::string> pdbID = structure.MapLabelToPDB(i.asymID, i.seqID, i.compID, i.authSeqID);

			stats.emplace_back(object{
				{ "asymID", i.asymID },
				{ "seqID", i.seqID },
				{ "compID", i.compID },
				{
					"pdb", {
						{ "strandID", std::get<0>(pdbID) },
						{ "seqNum", std::get<1>(pdbID) },
						{ "compID", std::get<2>(pdbID) },
						{ "insCode", std::get<3>(pdbID) }
					}
				},
				{ "RSR", i.RSR },
				{ "SRSR", i.SRSR },
				{ "RSCCS", i.RSCCS },
				{ "NGRID", i.ngrid },
				{ "EDIAm", i.EDIAm },
				{ "OPIA", i.OPIA }
			});
		}
		
		std::cout << stats << std::endl;
	}
	else
	{
		std::cout << "RESIDUE" << '\t'
			 << "RSR" << '\t'
			 << "SRSR" << '\t'
			 << "RSCCS" << '\t'
			 << "NGRID" << '\t'
			 << "EDIAm" << '\t'
			 << "OPIA" << std::endl;
	
		bool writeAuth = vm.count("use-auth-ids");
	
		for (auto i: r)
		{
			std::string id;
			
			if (writeAuth)
			{
				std::tuple<std::string,int,std::string,std::string> pdbID = structure.MapLabelToPDB(i.asymID, i.seqID, i.compID, i.authSeqID);
				id = std::get<2>(pdbID) + '_' + std::get<0>(pdbID) + '_' + std::to_string(std::get<1>(pdbID)) + std::get<3>(pdbID);
			}
			else if (i.compID == "HOH")
				id = i.compID + '_' + i.asymID + '_' + i.authSeqID;
			else
				id = i.compID + '_' + i.asymID + '_' + std::to_string(i.seqID);
			
			std::cout << std::fixed << std::setprecision(3)
				 << id << '\t'
				 << i.RSR << '\t'
				 << i.SRSR << '\t'
				 << i.RSCCS << '\t'
				 << i.ngrid << '\t'
				 << i.EDIAm << '\t'
				 << std::setprecision(1) << i.OPIA << std::endl;
		}
	}
	
	return 0;
}
