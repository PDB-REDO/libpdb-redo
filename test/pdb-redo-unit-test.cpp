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

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <stdexcept>
#include <filesystem>

#include <cif++.hpp>
#include <cif++/text.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/SkipList.hpp"
#include "pdb-redo/Symmetry-2.hpp"

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

BOOST_AUTO_TEST_CASE(symm_1)
{
	if (not fs::exists(gTestDir / "1mdn.mtz"))
		return;

	const std::string spacegroup = "I 1 21 1";

	int nr = cif::get_space_group_number(spacegroup);
	BOOST_TEST(nr == 5005);

	// clipper::Spacegroup sgr(clipper::Spgr_descr{ spacegroup, clipper::Spgr_descr::HM });

	// BOOST_TEST(sgr.symbol_hm() == spacegroup);

	using clipper::CCP4MTZfile;

	CCP4MTZfile mtzin;
	mtzin.open_read((gTestDir / "1mdn.mtz").string());

	HKL_info mHKLInfo;
	mtzin.import_hkl_info(mHKLInfo);

	nr = pdb_redo::getSpacegroupNumber(mHKLInfo.spacegroup());

	BOOST_TEST(nr == 5005);
}

BOOST_AUTO_TEST_CASE(symm_2)
{
	cif::file f(gTestDir / "3bwh.cif.gz");
	cif::mm::structure s(f);

	auto sg = getSpacegroup(f.front());
	auto c = getCell(f.front());

	for (int i = 0; i < sg.num_symops(); ++i)
	{
		auto csymop = sg.symop(i);

		auto rt = clipper::RTop_frac(csymop.rot(), csymop.trn()).rtop_orth(c);

		auto d = describeRToperation(sg, c, rt);

		std::cout << d << std::endl;

		sym_op so(d);

		BOOST_TEST(so.string() == d);

		auto rt2 = so.toClipperOrth(sg, c);

		BOOST_TEST(rt.equals(rt2, 0.00001));
	}
}

BOOST_AUTO_TEST_CASE(symm_4, *utf::tolerance(0.1f))
{
	// based on 2b8h
	auto sg = clipper::Spacegroup(clipper::Spgr_descr(154)); // p 32 2 1
	auto c = clipper::Cell(clipper::Cell_descr(107.516, 107.516, 338.487, 90.00, 90.00, 120.00));
	
	// cif::point a{-23.016, 47.514, -8.469};	// C1 NAG F 1
	// cif::point b{-16.707, 59.364, 24.320};	// CG1 VAL A 41-122

	// cif::point sb(59.76,15.21,-24.32);		// 4_555 copy of b

	cif::point a{   -8.688,  79.351, 10.439 }; // O6 NAG A 500
	cif::point b{  -35.356,  33.693, -3.236 }; // CG2 THR D 400
	cif::point sb(  -6.916,   79.34,   3.236); // 4_565 copy of b

	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_455"))) == static_cast<float>(c.a()));
	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_545"))) == static_cast<float>(c.b()));
	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_554"))) == static_cast<float>(c.c()));


	auto sa2 = symmetryCopy(a, sg, c, sym_op("4_565"));

	auto sb2 = symmetryCopy(b, sg, c, sym_op("4_565"));
	BOOST_TEST(sb.m_x == sb2.m_x);
	BOOST_TEST(sb.m_y == sb2.m_y);
	BOOST_TEST(sb.m_z == sb2.m_z);

	// BOOST_TEST(distance(a, b) == 35.43f);
	// BOOST_TEST(distance(a, symmetryCopy(b, sg, c, sym_op("1_555"))) == 35.43f);
	BOOST_TEST(distance(a, sb2) == 7.42f);
}

BOOST_AUTO_TEST_CASE(symm_3)
{
	auto sg = clipper::Spacegroup(clipper::Spgr_descr(195)); // P 2 3, cubic
	auto c = clipper::Cell(clipper::Cell_descr(1, 2, 3));

	cif::point a{ 1, 0, 0 }, b{ 0, 1, 0 };
	auto kD = distance(a, b);

	BOOST_TEST(distance(a, symmetryCopy(b, sg, c, sym_op("1_555"))) == kD);

	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_455"))) == c.a());
	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_545"))) == c.b());
	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_554"))) == c.c());

	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_455"))) == distance(a, { 0, 0, 0 }));
	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_545"))) == distance(a, { 1, -2, 0 }));
	BOOST_TEST(distance(a, symmetryCopy(a, sg, c, sym_op("1_554"))) == distance(a, { 1, 0, -3 }));

}

// --------------------------------------------------------------------
// more symmetry tests

BOOST_AUTO_TEST_CASE(symm_2bi3_1, *utf::tolerance(0.1f))
{
	cif::file f(gTestDir / "2bi3.cif.gz");
	cif::mm::structure s(f);

	auto &db = f.front();

	auto sg = getSpacegroup(db);
	auto c = getCell(db);

	auto struct_conn = db["struct_conn"];
	for (const auto &[
			asym1, seqid1, authseqid1, atomid1, symm1,
			asym2, seqid2, authseqid2, atomid2, symm2,
			dist] : struct_conn.find<
				std::string,int,std::string,std::string,std::string,
				std::string,int,std::string,std::string,std::string,
				float>(
			cif::key("ptnr1_symmetry") != "1_555" or cif::key("ptnr2_symmetry") != "1_555",
			"ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id", "ptnr1_label_atom_id", "ptnr1_symmetry", 
			"ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_auth_seq_id", "ptnr2_label_atom_id", "ptnr2_symmetry", 
			"pdbx_dist_value"
		))
	{
		auto &r1 = s.get_residue(asym1, seqid1, authseqid1);
		auto &r2 = s.get_residue(asym2, seqid2, authseqid2);

		auto a1 = r1.get_atom_by_atom_id(atomid1);
		auto a2 = r2.get_atom_by_atom_id(atomid2);

		auto sa1 = symmetryCopy(a1, sg, c, sym_op(symm1));
		auto sa2 = symmetryCopy(a2, sg, c, sym_op(symm2));

		BOOST_TEST(distance(sa1, sa2) == dist);
	}
}

BOOST_AUTO_TEST_CASE(symm_3bwh_1)
{
	cif::file f(gTestDir / "3bwh.cif.gz");
	cif::mm::structure s(f);

	auto a = s.get_residue("B", 0, "6").get_atom_by_atom_id("O2");
	auto b = s.get_residue("A", 1, "").get_atom_by_atom_id("O");

	// {
	// 	pdb_redo::SymmetryAtomIteratorFactory saif(s);

	// 	// auto a = s.atoms().front();
	// 	for (auto sa : saif(b))
	// 		std::cout << sa << " " << sa.get_location() << " " << sa.symmetry() << std::endl;
	// }

	const auto &[d, rtop] = closestSymmetryCopy(getSpacegroup(s.get_datablock()), getCell(s.get_datablock()), a.get_location(), b.get_location());

	BOOST_CHECK(d < 3);
}

BOOST_AUTO_TEST_CASE(symm_3bwh_2, *utf::tolerance(0.1f))
{
	cif::file f(gTestDir / "3bwh.cif.gz");
	cif::mm::structure s(f);

	auto a = s.get_residue("B", 0, "6").get_atom_by_atom_id("O2");
	auto b = s.get_residue("A", 1, "").get_atom_by_atom_id("O");

	auto sb = symmetryCopy(b, getSpacegroup(s.get_datablock()), getCell(s.get_datablock()), sym_op("4_945"));

	BOOST_TEST(distance(a, sb) == 2.54f);

	pdb_redo::DistanceMap dm(s, 3);

	BOOST_TEST(dm(a.id(), b.id()) == 2.54f);
	BOOST_TEST(dm(b.id(), a.id()) == 2.54f);

	std::cout << "near a " << a << " " << a.get_location() << " " << a.symmetry() << std::endl;
	for (auto n : dm.near(a, 3))
		std::cout << n << " " << n.get_location() << " " << n.symmetry() << std::endl;

	std::cout << "near b " << b << " " << b.get_location() << " " << b.symmetry() << std::endl;
	for (auto n : dm.near(b, 3))
		std::cout << n << " " << n.get_location() << " " << n.symmetry() << std::endl;
}

BOOST_AUTO_TEST_CASE(symm_2b8h_1, *utf::tolerance(0.1f))
{
	cif::file f(gTestDir / "2b8h.cif.gz");
	cif::mm::structure s(f);

	pdb_redo::DistanceMap dm(s, 3);

	auto a = s.get_residue("F", 0, "1").get_atom_by_atom_id("C1");
	auto b = s.get_residue("A", 120, "").get_atom_by_atom_id("ND2");

	BOOST_TEST(dm(a.id(), b.id()) == 1.44f);
	BOOST_TEST(dm(b.id(), a.id()) == 1.44f);

	std::cout << "near a " << a << " " << a.get_location() << " " << a.symmetry() << std::endl;
	for (auto n : dm.near(a, 3))
		std::cout << n << " " << n.get_location() << " " << n.symmetry() << std::endl;

	std::cout << "near b " << b << " " << b.get_location() << " " << b.symmetry() << std::endl;
	for (auto n : dm.near(b, 3))
		std::cout << n << " " << n.get_location() << " " << n.symmetry() << std::endl;
}




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
