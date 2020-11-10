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

// #include <iomanip>
#include <fstream>
#include <filesystem>
// #include <unordered_set>

#include <boost/program_options.hpp>
// #include <boost/iostreams/concepts.hpp>    // output_filter
// #include <boost/iostreams/operations.hpp>  // put
// #include <boost/iostreams/filtering_stream.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/Structure.hpp"
// #include "cif++/Symmetry.hpp"

// #include "Symmetry-2.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	int result = 0;

	po::options_description visible_options("carbonanza " + VERSION_STRING + " options file]" );
	visible_options.add_options()
		("output,o",	po::value<std::string>(),	"The output file, default is stdout")
		("help,h",									"Display help message")
		("version",									"Print version")
		("verbose,v",								"Verbose output")
		("dict",		po::value<std::string>(),	"Dictionary file containing restraints for residues in this specific target")
		;

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",		po::value<std::string>(),	"Input files")
		("test",									"Run test-suite")
		("debug,d",		po::value<int>(),			"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "platonyzer.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "platonyzer.conf";

	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}

	po::notify(vm);

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		std::cerr << visible_options << std::endl;
		exit(1);
	}

	// Load dict, if any

	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<std::string>());

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	if (cif::VERBOSE)
		std::cerr << "Loading data...";

	fs::path input = vm["input"].as<std::string>();
	c::File pdb(input);
	c::Structure structure(pdb);

	if (cif::VERBOSE)
		std::cerr << " done" << std::endl;

	auto& db = pdb.data();

	// -----------------------------------------------------------------------

	std::string entryId = db["entry"].front()["id"].as<std::string>();
	if (entryId.empty())
		throw std::runtime_error("Missing _entry.id in coordinates file");

	auto& nonpoly_scheme = db["pdbx_nonpoly_scheme"];
	auto& struct_conn = db["struct_conn"];
	auto& atom_site = db["atom_site"];

	for (auto& r: nonpoly_scheme.find(cif::Key("mon_id") == "NAG" or cif::Key("mon_id") == "NDG"))
	{
		std::string mon_id = r["mon_id"].as<>();
		std::string asym_id = r["asym_id"].as<>();

		auto l = struct_conn.find(
				(cif::Key("ptnr1_label_asym_id") == asym_id and cif::Key("ptnr2_label_comp_id") == "ASN") or
				(cif::Key("ptnr2_label_asym_id") == asym_id and cif::Key("ptnr1_label_comp_id") == "ASN"));

		if (l)	// already linked to an ASN
		{
			// see if the O1 is gone

			auto o = atom_site.find(cif::Key("label_asym_id") == asym_id and cif::Key("label_atom_id") == "O1");
			if (o)
			{
				if (cif::VERBOSE)
					std::cerr << "Deleting O1 atom from " << mon_id << " in asym " << asym_id << std::endl;
				
				atom_site.erase(cif::Key("label_asym_id") == asym_id and cif::Key("label_atom_id") == "O1");
			}

			// move 

			continue;
		}

		std::cout << r << std::endl;

	}

	db.add_software("carbonanza", "other", get_version_nr(), get_version_date());


	return 0;
}