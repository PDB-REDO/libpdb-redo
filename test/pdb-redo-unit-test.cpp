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

#define CATCH_CONFIG_RUNNER

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/SkipList.hpp"

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <filesystem>

#include <cif++/cif++.hpp>
#include <cif++/text.hpp>

namespace fs = std::filesystem;

using namespace pdb_redo;

// --------------------------------------------------------------------

std::filesystem::path gTestDir;

int main(int argc, char *argv[])
{
	gTestDir = std::filesystem::current_path();

	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;

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
// skip list test

TEST_CASE("skip_1")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure structure(file);

	auto& chain = structure.polymers().front();

	SkipList skiplist;

	int n = 10;
	for (auto& res: chain)
	{
		skiplist.emplace_back(res);
		if (--n == 0)
			break;
	}

	for (const auto fmt: { SkipListFormat::OLD, SkipListFormat::CIF })
	{
		std::stringstream ss;
		writeSkipList(ss, skiplist, fmt);

		CHECK(skiplist.size() == 10);

		// std::cout << ss.str() << '\n';

		SkipList list2 = readSkipList(ss);
		
		CHECK(list2.size() == skiplist.size());

		if (list2.size() != skiplist.size())
			continue;

		for (std::size_t i = 0; i < skiplist.size(); ++i)
		{
			CHECK(skiplist[i].pdb_asym_id == list2[i].pdb_asym_id);
			CHECK(skiplist[i].pdb_seq_id == list2[i].pdb_seq_id);
			CHECK(skiplist[i].pdbx_PDB_ins_code == list2[i].pdbx_PDB_ins_code);
		}
	}
}

// --------------------------------------------------------------------

// atom radii calculated with AtomShape and a minimizer from GSL
struct TestRadius {
	std::string		type;
	float			radius;
} kTestRadii[] = {
	{ "N", 1.073270559310913086f },
	{ "C", 1.077472805976867676f },
	{ "C", 1.060930848121643066f },
	{ "O", 1.037933468818664551f },
	{ "C", 1.080411434173583984f },
	{ "C", 1.080696582794189453f },
	{ "C", 1.090956211090087891f },
	{ "N", 1.02884829044342041f },
	{ "C", 1.017064213752746582f },
	{ "C", 0.985809326171875f },
	{ "O", 0.9498787522315979004f },
	{ "C", 1.032562017440795898f },
	{ "C", 1.043723225593566895f },
	{ "O", 1.000524282455444336f },
	{ "N", 1.057830214500427246f },
	{ "N", 0.9573949575424194336f },
	{ "C", 0.9534692764282226562f },
	{ "C", 0.9520707130432128906f },
	{ "O", 0.9262598156929016113f },
	{ "C", 0.938775181770324707f },
	{ "C", 0.9474387764930725098f },
	{ "C", 0.9400410056114196777f },
	{ "C", 0.9585416316986083984f },
	{ "C", 0.9261589646339416504f },
	{ "C", 0.9467949271202087402f },
	{ "C", 0.935360252857208252f },
	{ "N", 0.930846095085144043f },
	{ "C", 0.9431300163269042969f },
	{ "C", 0.9361689090728759766f },
	{ "O", 0.9119053483009338379f },
	{ "C", 0.9605298042297363281f },
	{ "O", 0.9615512490272521973f },
	{ "N", 0.9203097224235534668f },
	{ "C", 0.9396781325340270996f },
	{ "C", 0.9424930214881896973f },
	{ "O", 0.9095469117164611816f },
	{ "N", 0.9243852496147155762f },
	{ "C", 0.9413107037544250488f },
	{ "C", 0.9356296658515930176f },
	{ "O", 0.9198570847511291504f },
	{ "C", 0.9559983015060424805f },
	{ "C", 0.9961333274841308594f },
	{ "O", 0.9828038811683654785f },
	{ "N", 1.007446646690368652f },
	{ "N", 0.9110313653945922852f },
	{ "C", 0.9231967926025390625f },
	{ "C", 0.9301134943962097168f },
	{ "O", 0.9059692621231079102f },
	{ "C", 0.917514503002166748f },
	{ "C", 0.9176003932952880859f }
};

TEST_CASE("atom_shape_1")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");

	cif::file file(example.string());
	cif::mm::structure structure(file);

	const float kResHi = 1.80009, kResLo = 7.99918; 

	const std::size_t N = sizeof(kTestRadii) / sizeof(TestRadius);
	std::size_t i = 0;

	for (auto& atom: structure.atoms())
	{
		if (i >= N)
			break;

		pdb_redo::AtomShape shape(atom, kResHi, kResLo, false);

		CHECK(cif::atom_type_traits(atom.get_type()).symbol() == kTestRadii[i].type);

		float radius = shape.radius();
		float test = kTestRadii[i].radius;

		CHECK_THAT(radius, Catch::Matchers::WithinAbs(test, 0.001f));

		++i;
	}
}

// --------------------------------------------------------------------

cif::file operator""_cf(const char* text, std::size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, std::size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::file(is);
}

// --------------------------------------------------------------------

TEST_CASE("map_maker_1")
{
	pdb_redo::MapMaker<float> mm;
	
	float samplingRate = 0.75;

	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	CHECK_THAT(mm.resHigh(), Catch::Matchers::WithinAbs(1.8, 0.01));
	CHECK_THAT(mm.resLow(), Catch::Matchers::WithinAbs(8.0, 0.01));
}

TEST_CASE("map_maker_2")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");

	cif::file file(example.string());
	cif::mm::structure structure(file);

	MapMaker<float> mm;
	
	float samplingRate = 0.75;

	auto aniso = MapMaker<float>::as_None;
	// 	if (vm.count("aniso-scaling"))
	// 	{
	// 		if (vm["aniso-scaling"].get<std::string>() == "observed")
	// 			aniso = MapMaker<float>::as_Observed;
	// 		else if (vm["aniso-scaling"].get<std::string>() == "calculated")
	// 			aniso = MapMaker<float>::as_Calculated;
	// 	}
		
	mm.calculate(gTestDir / ".." / "examples" / "1cbs_map.mtz", structure, false, aniso, samplingRate, false);

	CHECK_THAT(mm.resHigh(), Catch::Matchers::WithinAbs(1.8, 0.01));
	CHECK_THAT(mm.resLow(), Catch::Matchers::WithinAbs(8.0, 0.01));
}

// --------------------------------------------------------------------

// First residues from 1cbs

struct TestResidue
{
	std::string		compID;
	std::string		asymID;
	int				seqID;
	double			RSR, SRSR, RSCCS;
	int				NGRID;
	double			EDIAm, OPIA;
};

TEST_CASE("stats_1")
{
	cif::VERBOSE = 2;

	// read test data first (output from previous stats version)

	std::vector<TestResidue> test;
	std::ifstream testFile(gTestDir / "1cbs-test.eds");

	REQUIRE(testFile.is_open());

	std::string line;
	std::getline(testFile, line);
	REQUIRE(line.substr(0, 7) == "RESIDUE");

	while (std::getline(testFile, line))
	{
		auto items = cif::split<std::string>(line, "\t");

		REQUIRE(items.size() == 7);

		auto id = cif::split<std::string>(items[0], "_");
		REQUIRE(id.size() == 3);

		test.push_back({
			id[0],
			id[1],
			std::stoi(id[2]),
			std::stod(items[1]),
			std::stod(items[2]),
			std::stod(items[3]),
			std::stoi(items[4]),
			std::stod(items[5]),
			std::stod(items[6])
		});
	}

	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());
	file.front().load_dictionary();
	// cif::mm::structure structure(file);

	MapMaker<float> mm;
	float samplingRate = 1.5;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	pdb_redo::EDIAStatsCollector collector(mm, file.front(), 1, false);
	auto r = collector.collect();

	auto ti = test.begin();
	for (auto& ri: r)
	{
		REQUIRE(ti != test.end());

		auto t = *ti++;

		if (std::isnan(t.RSCCS))
			continue;

		CHECK(ri.asymID == t.asymID);
		CHECK(ri.compID == t.compID);

		CHECK_THAT(std::abs(ri.RSR), Catch::Matchers::WithinAbs(t.RSR, 0.01));
		CHECK_THAT(std::abs(ri.SRSR), Catch::Matchers::WithinAbs(t.SRSR, 0.01));

		if (not (std::isnan(ri.RSCCS) and std::isnan(t.RSCCS)))
			CHECK_THAT(std::abs(ri.RSCCS), Catch::Matchers::WithinAbs(t.RSCCS, 0.1));
		else
			CHECK(std::isnan(ri.RSCCS) == std::isnan(t.RSCCS));

		if (not (std::isnan(ri.EDIAm) or std::isnan(t.EDIAm)))
		{
			CHECK_THAT(std::abs(ri.EDIAm), Catch::Matchers::WithinAbs(t.EDIAm, 0.1));

			// if (std::abs(ri.EDIAm - t.EDIAm) > 0.1)
			// 	std::cerr << ri << '\n';

			CHECK_THAT(std::abs(ri.OPIA), Catch::Matchers::WithinAbs(t.OPIA, 0.1));
		}
		else
		{
			CHECK(std::isnan(ri.EDIAm) == std::isnan(t.EDIAm));
		}

		CHECK(ri.ngrid == t.NGRID);
	}
}

// --------------------------------------------------------------------
// test stats on a file with an unkown residue, should give nan's for EDIA

TEST_CASE("stats_2")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file = cif::pdb::read(example.string());
	
	CHECK(file.is_valid());

	auto &db = file.front();
	
	// Rename a compound to an unknown ID

	for (auto r: db["chem_comp"].find(cif::key("id") == "ALA"))
	{
		r["id"] = "U_K";
		break;
	}

	for (auto r: db["pdbx_poly_seq_scheme"].find(cif::key("mon_id") == "ALA"))
		r["mon_id"] = "U_K";

	file.save("/tmp/1cbs-test.cif");

	// and load this into a structure (note, structure caches data from the file, so order is important)

	// cif::mm::structure structure(file);
	file.front().load_dictionary();

	MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	pdb_redo::EDIAStatsCollector collector(mm, file.front(), 1, false);
	auto r = collector.collect();

	for (auto& ri: r)
	{
		CHECK(ri.compID != "ALA");

		if (ri.compID != "U_K")
			continue;

		CHECK(not std::isnan(ri.RSR));
		CHECK(not std::isnan(ri.SRSR));
		CHECK(not std::isnan(ri.RSCCS));

		CHECK(std::isnan(ri.EDIAm));
		CHECK(std::isnan(ri.OPIA));
	}
}

// --------------------------------------------------------------------

TEST_CASE("bond_map_1")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file = cif::pdb::read(example.string());
	
	CHECK(file.is_valid());

	auto &db = file.front();
	auto &atom_site = db["atom_site"];

	std::vector<std::string> atom_ids;
	std::vector<cif::point> atom_locs;
	for (const auto &[id, x, y, z] : atom_site.rows<std::string,float,float,float>("id", "cartn_x", "cartn_y", "cartn_z"))
	{
		atom_ids.emplace_back(id);
		atom_locs.emplace_back(x, y, z);
	}

	BondMap bm(db);

	using key_type = std::tuple<std::string,std::string>;
	std::map<key_type,bool> bonded;
	
	for (std::size_t i = 0; i + 1 < atom_ids.size(); ++i)
	{
		auto a = atom_ids[i];

		for (std::size_t j = i + 1; j < atom_ids.size(); ++j)
		{
			auto b = atom_ids[j];
			
			bonded.emplace(std::make_tuple(a, b), bm(a, b));
		}
	}

	cif::point c = atom_locs.front();
	BondMap bm2(db, std::make_tuple(c, 6.0f));

	for (std::size_t i = 0; i + 1 < atom_ids.size(); ++i)
	{
		auto a = atom_ids[i];
		auto pa = atom_locs[i];

		for (std::size_t j = i + 1; j < atom_ids.size(); ++j)
		{
			auto b = atom_ids[j];
			auto pb = atom_locs[j];

			if (distance(pa, c) < 6 and distance(pb, c) < 6)
				CHECK(bm2(a, b) == bonded[std::make_tuple(a, b)]);
			else
				// CHECK_THROW(bm2(a, b), std::out_of_range);
				CHECK(bm2(a, b) == false);
		}
	}
}