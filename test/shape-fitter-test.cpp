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
#include "pdb-redo/BlobFinder.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cif++/pdb.hpp>
#include <clipper/core/xmap.h>
#define CATCH_CONFIG_RUNNER

#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/ShapeFitter.hpp"

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cif++.hpp>
#include <filesystem>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

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

	auto cf = R"(
data_1CBS
# 
_entry.id   1CBS
# 
_cell.entry_id           1CBS
_cell.length_a           45.650
_cell.length_b           47.560
_cell.length_c           77.610
_cell.angle_alpha        90.00
_cell.angle_beta         90.00
_cell.angle_gamma        90.00
_cell.Z_PDB              4
_cell.pdbx_unique_axis   ?
# 
_symmetry.entry_id                         1CBS
_symmetry.space_group_name_H-M             'P 21 21 21'
_symmetry.pdbx_full_space_group_name_H-M   ?
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                19
# 
	)"_cf;

	// Create a ligand
	cif::datablock &db = cf.front(); // almost empty
	db.set_validator(cif::validator_factory::instance().get("mmcif_pdbx.dic"));
	cif::mm::structure s(db);

	pdb_redo::BlobFinder bf(xmap, 0);

	auto blob = bf.next();

	CHECK(blob.size() == 831);

	auto ligand_asym_id = s.create_non_poly("REA", true);

	std::cout << "ligand: " << ligand_asym_id << " created\n";

	auto score = pdb_redo::fitShape(s, ligand_asym_id, xmap, blob);

	CHECK(score < 0);

	// std::ofstream file(std::filesystem::temp_directory_path() / "test.cif");
	// cf.save(file);
}

// --------------------------------------------------------------------

TEST_CASE("sf-2")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());
	cif::mm::structure s(file);
	s.remove_residue(s.get_residue("B"));

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

    auto &mm_fb = mm.fb();
    auto maskedmap = mm_fb.masked(s, s.atoms());

	pdb_redo::BlobFinder blobFinder(maskedmap, s);

	auto ligand_asym_id = s.create_non_poly("REA", true);

	for (;;)
	{
		auto blob = blobFinder.next();

		auto score = pdb_redo::fitShape(s, ligand_asym_id, mm_fb, blob);
	
		CHECK(score < 0);
	
		// std::ofstream of(std::filesystem::temp_directory_path() / "test-2.cif");
		// file.save(of);

		break;
	}

}

// // --------------------------------------------------------------------

// TEST_CASE("sf-3")
// {
// 	const fs::path example(gTestDir / "3aba_besttls.cif.gz");
// 	cif::file file = cif::pdb::read(example.string());

// 	cif::pdb::reconstruct_pdbx(file);

// 	cif::mm::structure s(file);
// 	// auto &db = s.get_datablock();

// 	for (std::string asymID : { "H", "I", "J", "K", "L"})
// 	{
// 		s.remove_residue(s.get_residue(asymID));
	
// 		pdb_redo::MapMaker<float> mm;
// 		float samplingRate = 0.75;
// 		mm.loadMTZ(gTestDir / "3aba_loopwhole.mtz", samplingRate);
	
// 		auto &mm_fb = mm.fb();
// 		auto maskedmap = mm_fb.masked(s, s.atoms());
	
// 		pdb_redo::BlobFinder blobFinder(maskedmap, s);
	
// 		auto ligand_asym_id = s.create_non_poly("GOL", true);
	
// 		for (int i = 0;; ++i)
// 		{
// 			auto blob = blobFinder.next();
// 			if (blob.empty())
// 				break;
	
// 			auto score = pdb_redo::fitShape(s, ligand_asym_id, mm_fb, blob);
		
// 			if (score < 0)
// 			{
// 				std::ofstream of(std::filesystem::temp_directory_path() / std::format("{}-{}-{}.cif", "3aba", asymID, i));
// 				file.save(of);
// 			}
// 		}
// 	}
// }
