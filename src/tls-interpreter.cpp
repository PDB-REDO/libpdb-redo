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

#include <termios.h>
#include <sys/ioctl.h>

#include <iostream>
#include <iomanip>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "cif++/CifUtils.hpp"
#include "cif++/Structure.hpp"
#include "cif++/TlsParser.hpp"

using namespace std;
namespace po = boost::program_options;
namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace c = mmcif;

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options("available tls-interpreter options");
	visible_options.add_options()
		("input,i",		po::value<string>(),	"Input PDB or mmCIF file. Can be compressed.")
		("output,o",	po::value<string>(),	"Output file, will write mmCIF unless extension is .pdb. Add .gz or .bz2 extension to create compressed file.")

		("mmcif-namespace",						"By default the TLS specification is interpreted as being in PDB namespace, with this option you can interpret the TLS specification as begin in mmCIF namespace")

		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output");
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cout << "Usage: tls-interpreter [options] input-file [output-file]" << endl
			 << endl
			 << visible_options << endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	bool cifNamespace = vm.count("mmcif-namespace") != 0;

	fs::path input = vm["input"].as<string>();
	c::File pdb(input);
//	c::structure structure(pdb, 1);
	
	auto& db = pdb.data();
	
	string program;
	
	for (auto p: db["software"].find(cif::Key("classification") == "refinement"))
	{
		program = p["name"].as<string>();
		break;
	}

	if (program.empty())
		cout << "No program specified when looking for _software.name with classification == refinement" << endl;
	else
		cout << "Parsing selection details for program " << program << endl;
	
	uint32_t group_id = 1;

	for (auto tls: db["pdbx_refine_tls"])
	{
		string refine_id, refine_tls_id;
		
		cif::tie(refine_id, refine_tls_id) = tls.get("pdbx_refine_id", "id");

		std::vector<tuple<std::string,int,int>> ranges;
		
		for (auto t: db["pdbx_refine_tls_group"].find(cif::Key("refine_tls_id") == refine_tls_id))
		{
			string selection = t["selection_details"].as<string>();
			
			if (selection.empty())
				continue;

			if (t["beg_auth_asym_id"].empty() == false)
			{
				cout << "beg_auth_asym_id already filled in, ignoring" << endl;
				continue;
			}

			cout << " --> " << selection << endl;
		
			auto s = cif::ParseSelectionDetails(program, selection);
	
			if (s == nullptr)
			{
				cerr << "Unable to parse" << endl
					 << endl;
				
				continue;
			}

			for (auto r: s->GetRanges(db, not cifNamespace))
			{
				string chainID;
				int from, to;
				tie(chainID, from, to) = r;
				
				cout << "RANGE  '" << chainID;
				
				if (from == cif::kNoSeqNum)
					cout << "    ";
				else
					cout << setw(4) << from;
				
				cout << ".' '" << chainID;
				
				if (to == cif::kNoSeqNum)
					cout << "    ";
				else
					cout << setw(4) << to;
				cout << ".' ALL" << endl;
				
				ranges.push_back(r);
			}
			
			cout << endl;
		}
		
		if (ranges.empty())
		{
			cout << "Empty range for TLS " << refine_tls_id << endl;
			continue;
		}

		db["pdbx_refine_tls_group"].erase(cif::Key("refine_tls_id") == refine_tls_id);
		
		for (auto r: ranges)
		{
			string asym_id;
			int from, to;
			
			tie(asym_id, from, to) = r;

			string from_s, to_s;
			if (from != cif::kNoSeqNum)
				from_s = to_string(from);
			else	
				from_s = ".";
			
			if (to != cif::kNoSeqNum)
				to_s = to_string(to);
			else
				to_s = ".";

			db["pdbx_refine_tls_group"].emplace({
				{ "pdbx_refine_id", refine_id },
				{ "id", "id_" + to_string(group_id++) },
				{ "refine_tls_id", refine_tls_id },
				{ "beg_label_asym_id", asym_id },
				{ "beg_label_seq_id", from_s },
				{ "end_label_asym_id", asym_id },
				{ "end_label_seq_id", to_s }
			});
		}
	}
	
	if (vm.count("output"))
		pdb.save(vm["output"].as<string>());
	
	return 0;
}
