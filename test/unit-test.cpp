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

#define BOOST_TEST_MODULE Libpdb_redo_Test
#include <boost/test/included/unit_test.hpp>

#include <boost/algorithm/string.hpp>

#include <stdexcept>
#include <filesystem>

#include <cif++/Structure.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/SkipList.hpp"

namespace fs = std::filesystem;
namespace tt = boost::test_tools;
namespace utf = boost::unit_test;
namespace ba = boost::algorithm;


// --------------------------------------------------------------------
// skip list test

BOOST_AUTO_TEST_CASE(skip_1)
{
	namespace c = mmcif;

	const fs::path example("../examples/1cbs.cif.gz");
	mmcif::File file(example);

	c::Structure structure(file);

	auto& chain = structure.polymers().front();

	SkipList skiplist;

	int n = 10;
	for (auto& res: chain)
	{
		skiplist.emplace_back(res);
		if (--n == 0)
			break;
	}

	for (const auto fmt: { SkipListFormat::OLD, SkipListFormat::JSON, SkipListFormat::CIF })
	{
		std::stringstream ss;
		writeSkipList(ss, skiplist, fmt);

		BOOST_CHECK(skiplist.size() == 10);

		// std::cout << ss.str() << std::endl;

		SkipList list2 = readSkipList(ss);
		
		BOOST_CHECK(list2.size() == skiplist.size());

		if (list2.size() != skiplist.size())
			continue;

		for (std::size_t i = 0; i < skiplist.size(); ++i)
		{
			BOOST_CHECK(skiplist[i].auth_asym_id == list2[i].auth_asym_id);
			BOOST_CHECK(skiplist[i].auth_seq_id == list2[i].auth_seq_id);
			BOOST_CHECK(skiplist[i].pdbx_PDB_ins_code == list2[i].pdbx_PDB_ins_code);
		}
	}
}

// --------------------------------------------------------------------

// atom radii calculated with AtomShape and NEWUOA
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

BOOST_AUTO_TEST_CASE(atom_shape_1, *utf::tolerance(0.0001f))
{
	const fs::path example("../examples/1cbs.cif.gz");

	mmcif::File file(example);
	mmcif::Structure structure(file);

	const float kResHi = 1.80009, kResLo = 7.99918; 

	const size_t N = sizeof(kTestRadii) / sizeof(TestRadius);
	size_t i = 0;

	for (auto& atom: structure.atoms())
	{
		if (i >= N)
			break;

		mmcif::AtomShape shape(atom, kResHi, kResLo, false);

		BOOST_CHECK(mmcif::AtomTypeTraits(atom.type()).symbol() == kTestRadii[i].type);

		float radius = shape.radius();
		float test = kTestRadii[i].radius;

		BOOST_TEST(radius == test);

		++i;
	}
}

// --------------------------------------------------------------------

cif::File operator""_cf(const char* text, size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::File(is);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(map_maker_1)
{
	namespace c = mmcif;

	c::MapMaker<float> mm;
	
	float samplingRate = 0.75;

	mm.loadMTZ("../examples/1cbs_map.mtz", samplingRate);

	BOOST_TEST(mm.resHigh() == 1.8, tt::tolerance(0.01));
	BOOST_TEST(mm.resLow() == 8.0, tt::tolerance(0.01));
}

BOOST_AUTO_TEST_CASE(map_maker_2)
{
	namespace c = mmcif;

	const fs::path example("../examples/1cbs.cif.gz");

	mmcif::File file(example);
	mmcif::Structure structure(file);

	c::MapMaker<float> mm;
	
	float samplingRate = 0.75;

	auto aniso = c::MapMaker<float>::as_None;
	// 	if (vm.count("aniso-scaling"))
	// 	{
	// 		if (vm["aniso-scaling"].as<std::string>() == "observed")
	// 			aniso = c::MapMaker<float>::as_Observed;
	// 		else if (vm["aniso-scaling"].as<std::string>() == "calculated")
	// 			aniso = c::MapMaker<float>::as_Calculated;
	// 	}
		
	mm.calculate("../examples/1cbs_map.mtz", structure, false, aniso, samplingRate, false);

	BOOST_TEST(mm.resHigh() == 1.8, tt::tolerance(0.01));
	BOOST_TEST(mm.resLow() == 8.0, tt::tolerance(0.01));
}

// --------------------------------------------------------------------

// First residues from 1cbs

struct TestResidue
{
	std::string		compID;
	std::string		asymID;
	int				seqID;
	double			RSR, SRSR, RSCCS;
	size_t			NGRID;
	double			EDIAm, OPIA;
};

BOOST_AUTO_TEST_CASE(stats_1)
{
	namespace c = mmcif;

	// read test data first (output from previous stats version)

	std::vector<TestResidue> test;
	std::ifstream testFile("1cbs-test.eds");

	BOOST_ASSERT(testFile.is_open());

	std::string line;
	std::getline(testFile, line);
	BOOST_ASSERT(line.substr(0, 7) == "RESIDUE");

	while (std::getline(testFile, line))
	{
		std::vector<std::string> fields;
		ba::split(fields, line, ba::is_any_of("\t"));

		BOOST_ASSERT(fields.size() == 7);

		std::vector<std::string> id;
		ba::split(id, fields[0], ba::is_any_of("_"));
		BOOST_ASSERT(id.size() == 3);

		test.push_back({
			id[0],
			id[1],
			std::stoi(id[2]),
			std::stod(fields[1]),
			std::stod(fields[2]),
			std::stod(fields[3]),
			static_cast<size_t>(std::stoi(fields[4])),
			std::stod(fields[5]),
			std::stod(fields[6])
		});
	}

	const fs::path example("../examples/1cbs.cif.gz");
	mmcif::File file(example);
	mmcif::Structure structure(file);

	c::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ("../examples/1cbs_map.mtz", samplingRate);

	mmcif::BondMap bm(structure);

	mmcif::EDIAStatsCollector collector(mm, structure, false, bm);
	auto r = collector.collect();

	auto ti = test.begin();
	for (auto& ri: r)
	{
		BOOST_ASSERT(ti != test.end());

		auto t = *ti++;

		if (std::isnan(t.RSCCS))
			continue;

		BOOST_TEST(ri.asymID == t.asymID);
		BOOST_TEST(ri.compID == t.compID);

		BOOST_TEST(std::abs(ri.RSR - t.RSR) <= 0.01, tt::tolerance(0.01));
		BOOST_TEST(std::abs(ri.SRSR - t.SRSR) <= 0.01, tt::tolerance(0.01));

		if (not (std::isnan(ri.RSCCS) and std::isnan(t.RSCCS)))
			BOOST_TEST(std::abs(ri.RSCCS - t.RSCCS) <= 0.1, tt::tolerance(0.1));
		else
			BOOST_CHECK(std::isnan(ri.RSCCS) == std::isnan(t.RSCCS));

		if (not (std::isnan(ri.EDIAm) or std::isnan(t.EDIAm)))
		{
			BOOST_TEST(std::abs(ri.EDIAm - t.EDIAm) <= 0.1, tt::tolerance(0.1));
			BOOST_TEST(std::abs(ri.OPIA - t.OPIA) <= 0.1, tt::tolerance(0.1));
		}
		else
		{
			BOOST_CHECK(std::isnan(ri.EDIAm) == std::isnan(t.EDIAm));
		}

		BOOST_TEST(ri.ngrid == t.NGRID);
	}
}

// --------------------------------------------------------------------
// test stats on a file with an unkown residue, should give nan's for EDIA

BOOST_AUTO_TEST_CASE(stats_2)
{
	namespace c = mmcif;

	const fs::path example("../examples/1cbs.cif.gz");
	mmcif::File file(example);

	BOOST_CHECK(file.file().isValid());
	
	// Rename a compound to an unknown ID

	for (auto r: file.data()["chem_comp"].find(cif::Key("id") == "ALA"))
	{
		r["id"] = "U_K";
		break;
	}

	// and load this into a structure (note, structure caches data from the file, so order is important)

	mmcif::Structure structure(file);

	c::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ("../examples/1cbs_map.mtz", samplingRate);

	mmcif::BondMap bm(structure);

	mmcif::EDIAStatsCollector collector(mm, structure, false, bm);
	auto r = collector.collect();

	for (auto& ri: r)
	{
		BOOST_CHECK(ri.compID != "ALA");

		if (ri.compID != "U_K")
			continue;

		BOOST_CHECK(not std::isnan(ri.RSR));
		BOOST_CHECK(not std::isnan(ri.SRSR));
		BOOST_CHECK(not std::isnan(ri.RSCCS));

		BOOST_CHECK(std::isnan(ri.EDIAm));
		BOOST_CHECK(std::isnan(ri.OPIA));
	}
}
