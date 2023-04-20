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

#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/SkipList.hpp"

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <stdexcept>
#include <filesystem>

#include <cif++.hpp>
#include <cif++/text.hpp>

namespace fs = std::filesystem;
namespace tt = boost::test_tools;
namespace utf = boost::unit_test;

using namespace pdb_redo;

// --------------------------------------------------------------------

fs::path gTestDir = fs::current_path();

bool init_unit_test()
{
	// not a test, just initialize test dir

	if (boost::unit_test::framework::master_test_suite().argc == 2)
	{
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

		if (fs::exists(gTestDir / "minimal-components.cif"))
			cif::compound_factory::instance().push_dictionary(gTestDir / "minimal-components.cif");
	}

	return true;
}


// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(q_4, *utf::tolerance(0.1f))
{
	std::vector<cif::point> pA{
		{ -0.188163757, 1.20824051, -0.0742692947 },
		{ -0.258220673, -0.184207916, -0.430510521 },
		{ 0.446380615, -1.02402496, 0.504780769 }
	};
	
	std::vector<cif::point> pB{
		{ -0.266670227, 0.961334228, 0.659666061 },
		{ 0.223335266, -0.319665909, 0.366666794 },
		{ 0.0433349609, -0.641666412, -1.02633286 }
	};

	auto cA = center_points(pA);
	auto cB = center_points(pB);

	auto q = cif::align_points(pB, pA);

	BOOST_TEST(q.get_a() == 0.137875929f);
	BOOST_TEST(q.get_b() == 0.15067713f);
	BOOST_TEST(q.get_c() == -0.942455589f);
	BOOST_TEST(q.get_d() == -0.264696091f);
}

// --------------------------------------------------------------------


// --------------------------------------------------------------------
// skip list test

BOOST_AUTO_TEST_CASE(skip_1)
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
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");

	cif::file file(example.string());
	cif::mm::structure structure(file);

	const float kResHi = 1.80009, kResLo = 7.99918; 

	const size_t N = sizeof(kTestRadii) / sizeof(TestRadius);
	size_t i = 0;

	for (auto& atom: structure.atoms())
	{
		if (i >= N)
			break;

		pdb_redo::AtomShape shape(atom, kResHi, kResLo, false);

		BOOST_CHECK_EQUAL(cif::atom_type_traits(atom.get_type()).symbol(), kTestRadii[i].type);

		float radius = shape.radius();
		float test = kTestRadii[i].radius;

		BOOST_TEST(radius == test);

		++i;
	}
}

// --------------------------------------------------------------------

cif::file operator""_cf(const char* text, size_t length)
{
    struct membuf : public std::streambuf
    {
        membuf(char* text, size_t length)
        {
            this->setg(text, text, text + length);
        }
    } buffer(const_cast<char*>(text), length);

    std::istream is(&buffer);
    return cif::file(is);
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(map_maker_1)
{
	pdb_redo::MapMaker<float> mm;
	
	float samplingRate = 0.75;

	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	BOOST_TEST(mm.resHigh() == 1.8, tt::tolerance(0.01));
	BOOST_TEST(mm.resLow() == 8.0, tt::tolerance(0.01));
}

BOOST_AUTO_TEST_CASE(map_maker_2)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");

	cif::file file(example.string());
	cif::mm::structure structure(file);

	MapMaker<float> mm;
	
	float samplingRate = 0.75;

	auto aniso = MapMaker<float>::as_None;
	// 	if (vm.count("aniso-scaling"))
	// 	{
	// 		if (vm["aniso-scaling"].as<std::string>() == "observed")
	// 			aniso = MapMaker<float>::as_Observed;
	// 		else if (vm["aniso-scaling"].as<std::string>() == "calculated")
	// 			aniso = MapMaker<float>::as_Calculated;
	// 	}
		
	mm.calculate(gTestDir / ".." / "examples" / "1cbs_map.mtz", structure, false, aniso, samplingRate, false);

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
	cif::VERBOSE = 2;

	// read test data first (output from previous stats version)

	std::vector<TestResidue> test;
	std::ifstream testFile(gTestDir / "1cbs-test.eds");

	BOOST_ASSERT(testFile.is_open());

	std::string line;
	std::getline(testFile, line);
	BOOST_ASSERT(line.substr(0, 7) == "RESIDUE");

	while (std::getline(testFile, line))
	{
		auto fields = cif::split<std::string>(line, "\t");

		BOOST_ASSERT(fields.size() == 7);

		auto id = cif::split<std::string>(fields[0], "_");
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

	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());
	cif::mm::structure structure(file);

	MapMaker<float> mm;
	float samplingRate = 1.5;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	BondMap bm(structure);

	pdb_redo::EDIAStatsCollector collector(mm, structure, false, bm);
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

			if (std::abs(ri.EDIAm - t.EDIAm) > 0.1)
				std::cerr << ri << std::endl;

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
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file = cif::pdb::read(example.string());
	
	BOOST_CHECK(file.is_valid());

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

	cif::mm::structure structure(file);

	MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	BondMap bm(structure);

	pdb_redo::EDIAStatsCollector collector(mm, structure, false, bm);
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

BOOST_AUTO_TEST_CASE(eigen_1, *utf::tolerance(0.1f))
{
	clipper::Matrix<float> m(4, 4);

	m(0, 0) = 4;
	m(1, 0) = m(0, 1) = -30;
	m(2, 0) = m(0, 2) = 60;
	m(3, 0) = m(0, 3) = -35;
	m(1, 1) = 300;
	m(2, 1) = m(1, 2) = -675;
	m(3, 1) = m(1, 3) = 420;
	m(2, 2) = 1620;
	m(3, 2) = m(2, 3) = -1050;
	m(3, 3) = 700;

	auto ev = m.eigen(true);

	BOOST_TEST(ev[0] == 0.1666428611718905f);
	BOOST_TEST(ev[1] == 1.4780548447781369f);
	BOOST_TEST(ev[2] == 37.1014913651276582f);
	BOOST_TEST(ev[3] == 2585.25381092892231f);

	BOOST_TEST(m(0, 0) == 0.792608291163763585f);
	BOOST_TEST(m(1, 0) == 0.451923120901599794f);
	BOOST_TEST(m(2, 0) == 0.322416398581824992f);
	BOOST_TEST(m(3, 0) == 0.252161169688241933f);

	BOOST_TEST(m(0, 1) == -0.582075699497237650f);
	BOOST_TEST(m(1, 1) == 0.370502185067093058f);
	BOOST_TEST(m(2, 1) == 0.509578634501799626f);
	BOOST_TEST(m(3, 1) == 0.514048272222164294f);

	// BOOST_TEST(m(0, 2) == -0.179186290535454826f);
	// BOOST_TEST(m(1, 2) == 0.741917790628453435f);
	// BOOST_TEST(m(2, 2) == -0.100228136947192199f);
	// BOOST_TEST(m(3, 2) == -0.638282528193614892f);

	BOOST_TEST(m(0, 3) == 0.0291933231647860588f);
	BOOST_TEST(m(1, 3) == -0.328712055763188997f);
	BOOST_TEST(m(2, 3) == 0.791411145833126331f);
	BOOST_TEST(m(3, 3) == -0.514552749997152907f);


}


BOOST_AUTO_TEST_CASE(eigen_2, *utf::tolerance(0.1f))
{
	cif::symmetric_matrix4x4<float> m;

	m(0, 0) = 4;
	m(0, 1) = -30;
	m(0, 2) = 60;
	m(0, 3) = -35;
	m(1, 1) = 300;
	m(1, 2) = -675;
	m(1, 3) = 420;
	m(2, 2) = 1620;
	m(2, 3) = -1050;
	m(3, 3) = 700;

	cif::matrix4x4<float> m2;
	m2 = m;

	const auto &[ev, em] = cif::eigen(m2, true);

	BOOST_TEST(ev[0] == 0.1666428611718905f);
	BOOST_TEST(ev[1] == 1.4780548447781369f);
	BOOST_TEST(ev[2] == 37.1014913651276582f);
	BOOST_TEST(ev[3] == 2585.25381092892231f);

	BOOST_TEST(em(0, 0) == 0.792608291163763585f);
	BOOST_TEST(em(1, 0) == 0.451923120901599794f);
	BOOST_TEST(em(2, 0) == 0.322416398581824992f);
	BOOST_TEST(em(3, 0) == 0.252161169688241933f);

	BOOST_TEST(em(0, 1) == -0.582075699497237650f);
	BOOST_TEST(em(1, 1) == 0.370502185067093058f);
	BOOST_TEST(em(2, 1) == 0.509578634501799626f);
	BOOST_TEST(em(3, 1) == 0.514048272222164294f);

	// BOOST_TEST(em(0, 2) == -0.179186290535454826f);
	// BOOST_TEST(em(1, 2) == 0.741917790628453435f);
	// BOOST_TEST(em(2, 2) == -0.100228136947192199f);
	// BOOST_TEST(em(3, 2) == -0.638282528193614892f);

	BOOST_TEST(em(0, 3) == 0.0291933231647860588f);
	BOOST_TEST(em(1, 3) == -0.328712055763188997f);
	BOOST_TEST(em(2, 3) == 0.791411145833126331f);
	BOOST_TEST(em(3, 3) == -0.514552749997152907f);
}


