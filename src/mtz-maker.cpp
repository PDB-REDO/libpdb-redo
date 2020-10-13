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

// my humble attempt to create mtz files

#include "pdb-redo.hpp"

#include <filesystem>

#include <boost/program_options.hpp>

#include "MapMaker.hpp"

using namespace std;
namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("make-mtz " + VERSION_STRING + " options");
	visible_options.add_options()
		("help,h",										"Display help message")
		("hklin",			po::value<string>(),		"Input file (either mtz or cif reflections file)")
		("xyzin",			po::value<string>(),		"Input coordinates file")
		("hklout",			po::value<string>(),		"Output MTZ file")
		("status-ccp4",									"In case the status flag contains '0' and '1', interpret status flag as CCP4 ('0' is free), the default is XPLOR ('1' is free)")
		("no-bulk",										"No bulk ")
		("electron-scattering",							"Use electron scattering factors")
		("aniso-obs",									"Anisotropic scaling of observed")
		("aniso-cal",									"Anisotropic scaling of calculated")
//		("num-reflns",		po::value<int>(),			"Number of reflections to use in structure factor weighting")
//		("num-params",		po::value<int>(),			"Number of spline parameters to use in structure factor weighting")
//		("fo-labels",		po::value<string>(),		"Comma seperated list of the (two) labels containing the observed reflections, default is 'FP,SIGFP'")
//		("free-label",		po::value<string>(),		"Label for the free flag")
		("version",										"Print version")
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("hklin") == 0 or vm.count("xyzin") == 0 or vm.count("hklout") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// --------------------------------------------------------------------

//	MTZMaker mm(vm["hklin"].as<string>(), vm["xyzin"].as<string>(), vm);

	fs::path hklin = vm["hklin"].as<string>();
	fs::path xyzin = vm["xyzin"].as<string>();

	mmcif::File f(xyzin);
	mmcif::Structure structure(f);
	
	c::MapMaker<float> mm;

	bool electronScattering = vm.count("electron-scattering") > 0;
	if (not electronScattering)
	{
		auto& exptl = f.data()["exptl"];
		electronScattering = not exptl.empty() and exptl.front()["method"] == "ELECTRON CRYSTALLOGRAPHY";
	}

	auto aniso = c::MapMaker<float>::as_None;
	if (vm.count("aniso-scaling"))
	{
		if (vm["aniso-scaling"].as<string>() == "observed")
			aniso = c::MapMaker<float>::as_Observed;
		else if (vm["aniso-scaling"].as<string>() == "calculated")
			aniso = c::MapMaker<float>::as_Calculated;
	}

	float samplingRate = 4.5;

//	int nRefln = 1000;
//	if (vm.count("num-reflns"))
//		nRefln = vm["num-reflns"].as<int>();
//	
//	int nParam = 20;
//	if (vm.count("num-params"))
//		nParam = vm["num-param"].as<int>();
		
	mm.calculate(hklin, structure, vm.count("no-bulk"), aniso, samplingRate, electronScattering);

	string name = f.data().getName();

	mm.writeMTZ(vm["hklout"].as<string>(), name, name);	

	return 0;
}
