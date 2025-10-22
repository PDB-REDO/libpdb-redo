/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/compound.hpp"
#include "cif++/datablock.hpp"
#include "cif++/model.hpp"
#include "cif++/validate.hpp"

#include <catch2/catch_test_macros.hpp>
#include <clipper/core/xmap.h>
#define CATCH_CONFIG_RUNNER

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Minimizer.hpp"
#include "pdb-redo/ShapeFitter.hpp"
#include "pdb-redo/Statistics.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cif++.hpp>
#include <filesystem>
#include <stdexcept>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

std::filesystem::path gTestDir = std::filesystem::current_path();

int main(int argc, char *argv[])
{
	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
#if CATCH22
	using namespace Catch::clara;
#else
	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;
#endif

	auto cli = session.cli()                                // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                  // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                   // the option names it will respond to
	           ("The directory containing the data files"); // description string for the help output

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	if (fs::exists(gTestDir / "minimal-components.cif"))
		cif::compound_factory::instance().push_dictionary(gTestDir / "minimal-components.cif");

	return session.run();
}

// --------------------------------------------------------------------

TEST_CASE("sf-1")
{
	pdb_redo::Map<float> map;
	map.read(gTestDir / "1cbs-REA-blob.map");
	clipper::Xmap<float> &xmap = map.get();

	std::vector<cif::point> blob;
	for (auto i = xmap.first(); not i.last(); i.next())
	{
		if (xmap[i] > 0)
			blob.emplace_back(i.coord_orth());
	}

	CHECK(blob.size() == 781);

	// Create a ligand
	cif::datablock db; // empty
    db.set_validator(&cif::validator_factory::instance().get("mmcif_pdbx.dic"));
	cif::mm::structure s(db);

	std::vector<cif::row_initializer> atoms;
    auto compound = cif::compound_factory::instance().create("REA");

	for (auto a : compound->atoms())
	{
		// We skip H-atoms, as fitting without H-atoms works better and we avoid conflicts in protonation states between CCD and MONLIB
		if (cif::atom_type_traits(a.type_symbol).symbol() == "H")
			continue;

		auto ax = a.get_location().get_x();
		auto ay = a.get_location().get_y();
		auto az = a.get_location().get_z();

		atoms.emplace_back(cif::row_initializer{
			{ "type_symbol", cif::atom_type_traits(a.type_symbol).symbol() },
			{ "label_atom_id", a.id },
			{ "auth_atom_id", a.id },
			{ "Cartn_x", ax },
			{ "Cartn_y", ay },
			{ "Cartn_z", az },
			{ "B_iso_or_equiv", 30.00 } });
	}
	
    auto ligand_entity_id = s.create_non_poly_entity("REA");
    auto ligand_asym_id = s.create_non_poly(ligand_entity_id, atoms);

    std::cout << "ligand: " << ligand_asym_id << " created\n";

    pdb_redo::fitShape(s, ligand_asym_id, xmap);
}