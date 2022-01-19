#define BOOST_TEST_MODULE Libpdb_redo_Test
#include <boost/test/included/unit_test.hpp>

#include <boost/algorithm/string.hpp>

#include <stdexcept>
#include <filesystem>

#include <cif++/Structure.hpp>
#include <cif++/CifUtils.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/Minimizer.hpp"

namespace fs = std::filesystem;
namespace tt = boost::test_tools;
namespace utf = boost::unit_test;
namespace ba = boost::algorithm;
namespace c = mmcif;

// --------------------------------------------------------------------

fs::path gTestDir = fs::current_path();

BOOST_AUTO_TEST_CASE(init)
{
	// not a test, just initialize test dir

	if (boost::unit_test::framework::master_test_suite().argc == 2)
	{
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

		if (fs::exists(gTestDir / "minimal-components.cif"))
			mmcif::CompoundFactory::instance().pushDictionary(gTestDir / "minimal-components.cif");
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(refine_0)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	mmcif::File file(example.string());

	c::Structure structure(file);

	auto &chain = structure.polymers().front();

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	c::BondMap bonds(structure);

	auto minimizer = pdb_redo::Minimizer::create(chain, 3, 3, bonds, mm.fb());

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << std::endl;

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	mmcif::File refFile(example.string());
	c::Structure reference(refFile);

	auto &refChain = reference.polymers().front();

	auto &atoms3 = chain.at(2).atoms();
	auto &refAtoms3 = refChain.at(2).atoms();

	BOOST_ASSERT(atoms3.size() == refAtoms3.size());

	for (size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		std::cout << a1 << ": " << a1.location() << " => " << a2.location() << "  distance: " << Distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '=') << std::endl;
}

BOOST_AUTO_TEST_CASE(refine_1)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	mmcif::File file(example.string());

	c::Structure structure(file);

	c::BondMap bonds(structure);

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	auto &chain = structure.polymers().front();
	auto &atoms3 = chain.at(2).atoms();

	auto minimizer = pdb_redo::Minimizer::create(chain, 3, 3, bonds, mm.fb());

	minimizer->printStats();

	// nudge the atoms of the third residue

	auto &res3 = chain.at(2);
	for (auto a : res3.atoms())
		structure.moveAtom(a, Nudge(a.location(), 0.5f));

	mmcif::File refFile(example.string());
	c::Structure reference(refFile);

	auto &refChain = reference.polymers().front();

	auto &refAtoms3 = refChain.at(2).atoms();

	BOOST_ASSERT(atoms3.size() == refAtoms3.size());

	for (size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		std::cout << a1 << ": " << a1.location() << " => " << a2.location() << "  distance: " << Distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << std::endl;

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	for (size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		std::cout << a1 << ": " << a1.location() << " => " << a2.location() << "  distance: " << Distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '=') << std::endl;
}

BOOST_AUTO_TEST_CASE(refine_2)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	mmcif::File file(example.string());

	c::Structure structure(file);

	const float kNearBy = 3;

	pdb_redo::DistanceMap dm(structure, kNearBy);

	// Move the REA residue somewhat

	auto &rea = structure.getResidue("B");

	std::vector<mmcif::Atom> atoms;
	for (auto a : rea.atoms())
	{
		for (auto b : dm.near(a, kNearBy))
		{
			if (find(atoms.begin(), atoms.end(), b) != atoms.end())
				continue;

			atoms.push_back(b);
		}
	}

	// translate by { 0.1, 0.1, 0.1 } and then
	// rotate around 1, 0, 0 for 5 degrees

	const float angle = 5 * (c::kPI / 180); // 5 degrees
	c::Quaternion q(
		std::cos(angle / 2), std::sin(angle / 2), 0, 0
	);

	for (auto a : rea.atoms())
		a.translateAndRotate({ 0.1, 0.1, 0.1 }, q);

	mmcif::File refFile(example.string());
	c::Structure reference(refFile);

	auto &atomsRea = rea.atoms();
	auto &refAtomsRea = reference.getResidue("B").atoms();

	BOOST_ASSERT(atomsRea.size() == refAtomsRea.size());

	for (size_t i = 0; i < atomsRea.size(); ++i)
	{
		auto a1 = atomsRea.at(i);
		auto a2 = refAtomsRea.at(i);

		std::cout << a1 << ": " << a1.location() << " => " << a2.location() << "  distance: " << Distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	c::BondMap bonds(structure);

	auto minimizer = pdb_redo::Minimizer::create(structure, atoms, bonds, mm.fb());

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << std::endl;

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	for (size_t i = 0; i < atomsRea.size(); ++i)
	{
		auto a1 = atomsRea.at(i);
		auto a2 = refAtomsRea.at(i);

		std::cout << a1 << ": " << a1.location() << " => " << a2.location() << "  distance: " << Distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;
}