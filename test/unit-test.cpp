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

#include <stdexcept>
#include <filesystem>

#include <cif++/Structure.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"

namespace fs = std::filesystem;
namespace tt = boost::test_tools;
namespace utf = boost::unit_test;

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

BOOST_AUTO_TEST_CASE(bondmap_1)
{
	// sections taken from CCD compounds.cif
	auto components = R"(
data_ASN
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
ASN N   CA   SING N N 1
ASN N   H    SING N N 2
ASN N   H2   SING N N 3
ASN CA  C    SING N N 4
ASN CA  CB   SING N N 5
ASN CA  HA   SING N N 6
ASN C   O    DOUB N N 7
ASN C   OXT  SING N N 8
ASN CB  CG   SING N N 9
ASN CB  HB2  SING N N 10
ASN CB  HB3  SING N N 11
ASN CG  OD1  DOUB N N 12
ASN CG  ND2  SING N N 13
ASN ND2 HD21 SING N N 14
ASN ND2 HD22 SING N N 15
ASN OXT HXT  SING N N 16
data_PHE
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
PHE N   CA  SING N N 1
PHE N   H   SING N N 2
PHE N   H2  SING N N 3
PHE CA  C   SING N N 4
PHE CA  CB  SING N N 5
PHE CA  HA  SING N N 6
PHE C   O   DOUB N N 7
PHE C   OXT SING N N 8
PHE CB  CG  SING N N 9
PHE CB  HB2 SING N N 10
PHE CB  HB3 SING N N 11
PHE CG  CD1 DOUB Y N 12
PHE CG  CD2 SING Y N 13
PHE CD1 CE1 SING Y N 14
PHE CD1 HD1 SING N N 15
PHE CD2 CE2 DOUB Y N 16
PHE CD2 HD2 SING N N 17
PHE CE1 CZ  DOUB Y N 18
PHE CE1 HE1 SING N N 19
PHE CE2 CZ  SING Y N 20
PHE CE2 HE2 SING N N 21
PHE CZ  HZ  SING N N 22
PHE OXT HXT SING N N 23
data_PRO
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
PRO N   CA  SING N N 1
PRO N   CD  SING N N 2
PRO N   H   SING N N 3
PRO CA  C   SING N N 4
PRO CA  CB  SING N N 5
PRO CA  HA  SING N N 6
PRO C   O   DOUB N N 7
PRO C   OXT SING N N 8
PRO CB  CG  SING N N 9
PRO CB  HB2 SING N N 10
PRO CB  HB3 SING N N 11
PRO CG  CD  SING N N 12
PRO CG  HG2 SING N N 13
PRO CG  HG3 SING N N 14
PRO CD  HD2 SING N N 15
PRO CD  HD3 SING N N 16
PRO OXT HXT SING N N 17
)"_cf;

	const fs::path example("/usr/share/doc/libcifpp-dev/examples/1cbs.cif.gz");
	mmcif::File file(example);
	mmcif::Structure structure(file);

	mmcif::BondMap bm(structure);

	// Test the bonds of the first three residues, that's PRO A 1, ASN A 2, PHE A 3	

	for (const auto& [compound, seqnr]: std::initializer_list<std::tuple<std::string,int>>{ { "PRO", 1 }, { "ASN", 2 }, { "PHE", 3 } })
	{
		auto& res = structure.getResidue("A", compound, seqnr);
		auto atoms = res.atoms();

		auto dc = components.get(compound);
		BOOST_ASSERT(dc != nullptr);

		auto cc = dc->get("chem_comp_bond");
		BOOST_ASSERT(cc != nullptr);

		std::set<std::tuple<std::string,std::string>> bonded;

		for (const auto& [atom_id_1, atom_id_2]: cc->rows<std::string,std::string>({ "atom_id_1", "atom_id_2" }))
		{
			if (atom_id_1 > atom_id_2)
				bonded.insert({ atom_id_2, atom_id_1 });
			else
				bonded.insert({ atom_id_1, atom_id_2 });
		}

		for (size_t i = 0; i + 1 < atoms.size(); ++i)
		{
			auto label_i = atoms[i].labelAtomID();

			for (size_t j = i + 1; j < atoms.size(); ++j)
			{
				auto label_j = atoms[j].labelAtomID();

				bool bonded_1 = bm(atoms[i], atoms[j]);
				bool bonded_1_i = bm(atoms[j], atoms[i]);

				bool bonded_t = label_i > label_j
					? bonded.count({ label_j, label_i })
					: bonded.count({ label_i, label_j });

				BOOST_CHECK(bonded_1 == bonded_t);
				BOOST_CHECK(bonded_1_i == bonded_t);
			}
		}
	}
}

// --------------------------------------------------------------------

// atom radii calculated with AtomShape and NEWUOA
struct TestRadius {
	std::string		type;
	float			radius;
} kTestRadii[] = {
	{ "N", 1.07327 },
	{ "C", 1.07747 },
	{ "C", 1.06093 },
	{ "O", 1.03793 },
	{ "C", 1.08041 },
	{ "C", 1.0807 },
	{ "C", 1.09096 },
	{ "N", 1.02885 },
	{ "C", 1.01706 },
	{ "C", 0.98581 },
	{ "O", 0.949879 },
	{ "C", 1.03256 },
	{ "C", 1.04372 },
	{ "O", 1.00052 },
	{ "N", 1.05783 },
	{ "N", 0.957396 },
	{ "C", 0.95347 },
	{ "C", 0.952071 },
	{ "O", 0.926261 },
	{ "C", 0.938776 },
	{ "C", 0.947439 },
	{ "C", 0.940042 },
	{ "C", 0.958542 },
	{ "C", 0.92616 },
	{ "C", 0.946796 },
	{ "C", 0.935361 },
	{ "N", 0.930847 },
	{ "C", 0.943131 },
	{ "C", 0.93617 },
	{ "O", 0.911906 },
	{ "C", 0.96053 },
	{ "O", 0.961552 },
	{ "N", 0.92031 },
	{ "C", 0.939679 },
	{ "C", 0.942494 },
	{ "O", 0.909548 },
	{ "N", 0.924386 },
	{ "C", 0.941311 },
	{ "C", 0.93563 },
	{ "O", 0.919858 },
	{ "C", 0.955999 },
	{ "C", 0.996134 },
	{ "O", 0.982804 },
	{ "N", 1.00745 },
	{ "N", 0.911032 },
	{ "C", 0.923198 },
	{ "C", 0.930114 },
	{ "O", 0.90597 },
	{ "C", 0.917515 },
	{ "C", 0.917601 }
};

BOOST_AUTO_TEST_CASE(atom_shape_1, *utf::tolerance(0.001f))
{
	const fs::path example("/usr/share/doc/libcifpp-dev/examples/1cbs.cif.gz");

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

	const fs::path example("/usr/share/doc/libcifpp-dev/examples/1cbs.cif.gz");

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

struct TestResidue {
	std::string		compID;
	std::string		asymID;
	int				seqID;
	double			RSR, SRSR, RSCCS;
	size_t			NGRID;
	double			EDIAm, OPIA;
} kTestResidues[] = {
	{ "PRO", "A", 1, 	0.116,	0.056,	0.849,	54,		0.008,		14.3 },
	{ "ASN", "A", 2, 	0.079,	0.036,	0.873,	46,		0.864,		87.5 },
	{ "PHE", "A", 3, 	0.027,	0.022,	0.980,	50,		0.863,		72.7 },
	{ "SER", "A", 4, 	0.087,	0.038,	0.877,	36,		0.668,		50.0 },
	{ "GLY", "A", 5, 	0.019,	0.035,	0.983,	18,		0.902,		100.0 },
	{ "ASN", "A", 6, 	0.065,	0.030,	0.874,	40,		0.641,		37.5 },
	{ "TRP", "A", 7, 	0.030,	0.018,	0.980,	56,		0.719,		57.1 },
	{ "LYS", "A", 8, 	0.080,	0.031,	0.911,	54,		0.057,		66.7 },
	{ "ILE", "A", 9, 	0.050,	0.026,	0.964,	36,		0.932,		87.5 },
	{ "ILE", "A", 10, 	0.054,	0.025,	0.938,	42,		0.841,		75.0 },
	{ "ARG", "A", 11, 	0.039,	0.023,	0.957,	50,		0.844,		54.5 },
	{ "SER", "A", 12, 	0.051,	0.030,	0.856,	28,		0.849,		83.3 },
	{ "GLU", "A", 13, 	0.047,	0.028,	0.959,	40,		0.476,		55.6 },
	{ "ASN", "A", 14, 	0.035,	0.025,	0.970,	36,		0.821,		75.0 },
	{ "PHE", "A", 15, 	0.036,	0.025,	0.973,	54,		0.925,		100.0 },
	{ "GLU", "A", 16, 	0.066,	0.028,	0.940,	52,		0.589,		55.6 },
	{ "GLU", "A", 17, 	0.103,	0.032,	0.897,	56,		0.137,		44.4 },
	{ "LEU", "A", 18, 	0.028,	0.025,	0.962,	30,		0.749,		37.5 },
	{ "LEU", "A", 19, 	0.039,	0.027,	0.960,	40,		0.670,		62.5 },
	{ "LYS", "A", 20, 	0.108,	0.040,	0.911,	60,		0.221,		44.4 },
	{ "VAL", "A", 21, 	0.058,	0.035,	0.951,	34,		0.752,		42.9 },
	{ "LEU", "A", 22, 	0.042,	0.033,	0.958,	40,		0.658,		50.0 },
	{ "GLY", "A", 23, 	0.078,	0.046,	0.720,	28,		1.013,		100.0 },
	{ "VAL", "A", 24, 	0.044,	0.028,	0.954,	28,		0.654,		28.6 },
	{ "ASN", "A", 25, 	0.080,	0.039,	0.922,	42,		0.827,		75.0 },
	{ "VAL", "A", 26, 	0.084,	0.039,	0.871,	44,		0.602,		28.6 },
	{ "MET", "A", 27, 	0.064,	0.033,	0.911,	52,		0.488,		62.5 },
	{ "LEU", "A", 28, 	0.057,	0.028,	0.929,	42,		0.502,		50.0 },
	{ "ARG", "A", 29, 	0.100,	0.027,	0.818,	68,		0.673,		36.4 },
	{ "LYS", "A", 30, 	0.117,	0.051,	0.868,	62,		0.299,		44.4 },
	{ "ILE", "A", 31, 	0.057,	0.029,	0.898,	44,		0.690,		50.0 },
	{ "ALA", "A", 32, 	0.056,	0.040,	0.935,	20,		0.871,		80.0 },
	{ "VAL", "A", 33, 	0.042,	0.035,	0.965,	36,		0.724,		28.6 },
	{ "ALA", "A", 34, 	0.061,	0.041,	0.900,	28,		0.845,		60.0 },
	{ "ALA", "A", 35, 	0.039,	0.030,	0.968,	26,		0.688,		20.0 },
	{ "ALA", "A", 36, 	0.032,	0.037,	0.975,	28,		0.952,		100.0 },
	{ "SER", "A", 37, 	0.065,	0.035,	0.926,	34,		0.816,		33.3 },
	{ "LYS", "A", 38, 	0.130,	0.037,	0.868,	54,		0.006,		55.6 },
	{ "PRO", "A", 39, 	0.064,	0.030,	0.953,	36,		0.909,		85.7 },
	{ "ALA", "A", 40, 	0.058,	0.038,	0.930,	28,		1.009,		100.0 },
	{ "VAL", "A", 41, 	0.056,	0.026,	0.929,	40,		0.776,		71.4 },
	{ "GLU", "A", 42, 	0.080,	0.028,	0.932,	54,		0.410,		55.6 },
	{ "ILE", "A", 43, 	0.036,	0.029,	0.968,	40,		0.715,		25.0 },
	{ "LYS", "A", 44, 	0.091,	0.033,	0.904,	48,		0.022,		66.7 },
	{ "GLN", "A", 45, 	0.052,	0.028,	0.932,	40,		0.589,		22.2 },
	{ "GLU", "A", 46, 	0.152,	0.041,	0.802,	66,		0.051,		55.6 },
	{ "GLY", "A", 47, 	0.049,	0.038,	0.956,	18,		0.712,		75.0 },
	{ "ASP", "A", 48, 	0.056,	0.031,	0.972,	50,		0.669,		50.0 },
	{ "THR", "A", 49, 	0.056,	0.032,	0.963,	40,		0.610,		28.6 },
	{ "PHE", "A", 50, 	0.047,	0.025,	0.946,	54,		0.735,		54.5 },
	{ "TYR", "A", 51, 	0.037,	0.024,	0.983,	52,		0.838,		75.0 },
	{ "ILE", "A", 52, 	0.046,	0.023,	0.934,	40,		0.821,		75.0 },
	{ "LYS", "A", 53, 	0.070,	0.027,	0.945,	46,		0.486,		55.6 },
	{ "THR", "A", 54, 	0.031,	0.026,	0.972,	26,		0.908,		71.4 },
	{ "SER", "A", 55, 	0.051,	0.030,	0.969,	26,		0.690,		66.7 },
	{ "THR", "A", 56, 	0.032,	0.024,	0.957,	30,		0.942,		100.0 },
	{ "THR", "A", 57, 	0.051,	0.025,	0.965,	34,		0.851,		57.1 },
	{ "VAL", "A", 58, 	0.042,	0.027,	0.948,	30,		0.620,		42.9 },
	{ "ARG", "A", 59, 	0.066,	0.030,	0.946,	54,		0.828,		63.6 },
	{ "THR", "A", 60, 	0.047,	0.028,	0.952,	32,		0.688,		71.4 },
	{ "THR", "A", 61, 	0.057,	0.030,	0.968,	32,		0.811,		85.7 },
	{ "GLU", "A", 62, 	0.091,	0.029,	0.882,	56,		0.841,		55.6 },
	{ "ILE", "A", 63, 	0.051,	0.029,	0.955,	36,		0.865,		75.0 },
	{ "ASN", "A", 64, 	0.057,	0.029,	0.931,	38,		0.830,		62.5 },
	{ "PHE", "A", 65, 	0.040,	0.022,	0.977,	52,		0.757,		36.4 },
	{ "LYS", "A", 66, 	0.130,	0.035,	0.815,	60,		0.525,		44.4 },
	{ "VAL", "A", 67, 	0.048,	0.038,	0.921,	34,		0.791,		71.4 },
	{ "GLY", "A", 68, 	0.068,	0.041,	0.893,	24,		0.877,		100.0 },
};

const size_t kNResidues = sizeof(kTestResidues) / sizeof(TestResidue);

BOOST_AUTO_TEST_CASE(stats_1)
{
	namespace c = mmcif;

	const fs::path example("/usr/share/doc/libcifpp-dev/examples/1cbs.cif.gz");
	mmcif::File file(example);
	mmcif::Structure structure(file);

	c::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ("../examples/1cbs_map.mtz", samplingRate);

	mmcif::BondMap bm(structure);

	mmcif::EDIAStatsCollector collector(mm, structure, false, bm);
	auto r = collector.collect();

	size_t i = 0;
	for (auto& ri: r)
	{
		if (i >= kNResidues)
			break;

		auto& t = kTestResidues[i++];

		BOOST_TEST(ri.asymID == t.asymID);
		BOOST_TEST(ri.compID == t.compID);
		BOOST_TEST(ri.seqID == t.seqID);

		BOOST_TEST(std::abs(ri.RSR - t.RSR) <= 0.01, tt::tolerance(0.01));
		BOOST_TEST(std::abs(ri.SRSR - t.SRSR) <= 0.01, tt::tolerance(0.01));
		BOOST_TEST(std::abs(ri.RSCCS - t.RSCCS) <= 0.01, tt::tolerance(0.01));

		BOOST_TEST(std::abs(ri.EDIAm - t.EDIAm) <= 0.1, tt::tolerance(0.1));
		BOOST_TEST(std::abs(ri.OPIA - t.OPIA) <= 0.1, tt::tolerance(0.1));

		BOOST_TEST(ri.ngrid == t.NGRID);
	}
}
