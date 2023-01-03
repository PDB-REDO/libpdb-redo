#define BOOST_TEST_MODULE Libpdb_redo_Test
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>
#include <filesystem>

#include <cif++.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/Minimizer.hpp"

namespace fs = std::filesystem;
namespace tt = boost::test_tools;
namespace utf = boost::unit_test;

// --------------------------------------------------------------------

fs::path gTestDir = fs::current_path();

BOOST_AUTO_TEST_CASE(init)
{
	cif::VERBOSE = 1;

	// not a test, just initialize test dir

	if (boost::unit_test::framework::master_test_suite().argc == 2)
	{
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

		if (fs::exists(gTestDir / "minimal-components.cif"))
			cif::compound_factory::instance().push_dictionary(gTestDir / "minimal-components.cif");
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(refine_0)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure structure(file);

	auto &chain = structure.polymers().front();

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	pdb_redo::BondMap bonds(structure);

	auto minimizer = pdb_redo::Minimizer::create(chain, 3, 3, bonds, mm.fb());

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << std::endl;

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	cif::file refFile(example.string());
	cif::mm::structure reference(refFile);

	auto &refChain = reference.polymers().front();

	auto &atoms3 = chain.at(2).atoms();
	auto &refAtoms3 = refChain.at(2).atoms();

	BOOST_ASSERT(atoms3.size() == refAtoms3.size());

	double d_sum = 0;

	for (size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		auto d = distance(a1, a2);
		d_sum += d * d;

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << d << std::endl;
	}

	auto rmsd = std::sqrt(d_sum / atoms3.size());
	std::cout << "RMSd: " << rmsd << std::endl;
	BOOST_CHECK(rmsd < 0.2);

	std::cout << std::string(cif::get_terminal_width(), '=') << std::endl;
}

BOOST_AUTO_TEST_CASE(refine_1)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure structure(file);

	pdb_redo::BondMap bonds(structure);

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
		structure.move_atom(a, cif::nudge(a.get_location(), 0.5f));

	cif::file refFile(example.string());
	cif::mm::structure reference(refFile);

	auto &refChain = reference.polymers().front();

	auto &refAtoms3 = refChain.at(2).atoms();

	BOOST_ASSERT(atoms3.size() == refAtoms3.size());

	for (size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << std::endl;

	minimizer->printStats();

	double d_sum = 0;

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	for (size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		auto d = distance(a1, a2);
		d_sum += d * d;

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << d << std::endl;
	}

	auto rmsd = std::sqrt(d_sum / atoms3.size());
	std::cout << "RMSd: " << rmsd << std::endl;
	BOOST_CHECK(rmsd < 0.2);

	std::cout << std::string(cif::get_terminal_width(), '=') << std::endl;
}

BOOST_AUTO_TEST_CASE(refine_2)
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure structure(file);

	const float kNearBy = 3;

	pdb_redo::DistanceMap dm(structure, kNearBy);

	// Move the REA residue somewhat

	auto &rea = structure.get_residue("B");

	std::vector<cif::mm::atom> atoms;
	for (auto a : rea.atoms())
	{
		for (auto b : dm.near(a, kNearBy))
		{
			if (find(atoms.begin(), atoms.end(), b) != atoms.end())
				continue;

			atoms.push_back(b);
		}
	}

	// // translate by { 0.1, 0.1, 0.1 } and then
	// // rotate around 1, 0, 0 for 5 degrees

	// const float angle = 5 * (cif::kPI / 180);
	// cif::quaternion q(
	// 	std::cos(angle / 2), std::sin(angle / 2), 0, 0
	// );

	// for (auto a : rea.atoms())
	// 	a.translate_and_rotate({ 0.1, 0.1, 0.1 }, q);

	for (auto a : rea.atoms())
	{
		auto l = a.get_location();
		l = nudge(l, 0.5);
		a.set_location(l);
	}

	cif::file refFile(example.string());
	cif::mm::structure reference(refFile);

	auto &atomsRea = rea.atoms();
	auto &refAtomsRea = reference.get_residue("B").atoms();

	BOOST_ASSERT(atomsRea.size() == refAtomsRea.size());

	for (size_t i = 0; i < atomsRea.size(); ++i)
	{
		auto a1 = atomsRea.at(i);
		auto a2 = refAtomsRea.at(i);

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << distance(a1, a2) << std::endl;
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	pdb_redo::BondMap bonds(structure);

	auto minimizer = pdb_redo::Minimizer::create(structure, atoms, bonds, mm.fb());

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << std::endl;

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;

	double d_sum = 0;

	for (size_t i = 0; i < atomsRea.size(); ++i)
	{
		auto a1 = atomsRea.at(i);
		auto a2 = refAtomsRea.at(i);

		auto d = distance(a1, a2);
		d_sum += d * d;

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << d << std::endl;
	}

	file.save("/tmp/rsr-test-3.cif");

	auto rmsd = std::sqrt(d_sum / atomsRea.size());
	std::cout << "RMSd: " << rmsd << std::endl;
	BOOST_CHECK(rmsd < 0.2);

	std::cout << std::string(cif::get_terminal_width(), '-') << std::endl;
}
