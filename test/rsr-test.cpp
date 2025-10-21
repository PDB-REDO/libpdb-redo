/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2024 NKI/AVL, Netherlands Cancer Institute
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

#include <stdexcept>
#include <filesystem>

#include <cif++.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/Minimizer.hpp"

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

TEST_CASE("refine_0")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure structure(file);
	cif::crystal crystal(structure.get_datablock());

	auto &chain = structure.polymers().front();

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	auto minimizer = pdb_redo::Minimizer::create(crystal, chain, 3, 3, mm.fb());

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << '\n';

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << '\n';

	cif::file refFile(example.string());
	cif::mm::structure reference(refFile);

	auto &refChain = reference.polymers().front();

	auto &atoms3 = chain.at(2).atoms();
	auto &refAtoms3 = refChain.at(2).atoms();

	CHECK(atoms3.size() == refAtoms3.size());

	double d_sum = 0;

	for (std::size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		auto d = distance(a1, a2);
		d_sum += d * d;

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << d << '\n';
	}

	auto rmsd = std::sqrt(d_sum / atoms3.size());
	std::cout << "RMSd: " << rmsd << '\n';
	CHECK_THAT(rmsd, Catch::Matchers::WithinAbs(0.125, 0.125));

	std::cout << std::string(cif::get_terminal_width(), '=') << '\n';
}

TEST_CASE("refine_1")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());

	cif::mm::structure structure(file);
	cif::crystal crystal(structure.get_datablock());

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	auto &chain = structure.polymers().front();
	auto &atoms3 = chain.at(2).atoms();

	auto minimizer = pdb_redo::Minimizer::create(crystal, chain, 3, 3, mm.fb());

	minimizer->printStats();

	// nudge the atoms of the third residue

	auto &res3 = chain.at(2);
	for (auto a : res3.atoms())
		structure.move_atom(a, cif::nudge(a.get_location(), 0.5f));

	cif::file refFile(example.string());
	cif::mm::structure reference(refFile);

	auto &refChain = reference.polymers().front();

	auto &refAtoms3 = refChain.at(2).atoms();

	CHECK(atoms3.size() == refAtoms3.size());

	for (std::size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << distance(a1, a2) << '\n';
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << '\n';

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << '\n';

	minimizer->printStats();

	double d_sum = 0;

	std::cout << std::string(cif::get_terminal_width(), '-') << '\n';

	for (std::size_t i = 0; i < atoms3.size(); ++i)
	{
		auto a1 = atoms3.at(i);
		auto a2 = refAtoms3.at(i);

		auto d = distance(a1, a2);
		d_sum += d * d;

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << d << '\n';
	}

	auto rmsd = std::sqrt(d_sum / atoms3.size());
	std::cout << "RMSd: " << rmsd << '\n';
	CHECK(rmsd < 0.3);

	std::cout << std::string(cif::get_terminal_width(), '=') << '\n';
}

TEST_CASE("refine_2")
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

	CHECK(atomsRea.size() == refAtomsRea.size());

	for (std::size_t i = 0; i < atomsRea.size(); ++i)
	{
		auto a1 = atomsRea.at(i);
		auto a2 = refAtomsRea.at(i);

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << distance(a1, a2) << '\n';
	}

	std::cout << std::string(cif::get_terminal_width(), '-') << '\n';

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	cif::crystal crystal(structure.get_datablock());
	auto minimizer = pdb_redo::Minimizer::create(crystal, structure, atoms, mm.fb());

	minimizer->printStats();

	auto score = minimizer->refine(true);

	std::cout << "minimizer score: " << score << '\n';

	minimizer->printStats();

	std::cout << std::string(cif::get_terminal_width(), '-') << '\n';

	double d_sum = 0;

	for (std::size_t i = 0; i < atomsRea.size(); ++i)
	{
		auto a1 = atomsRea.at(i);
		auto a2 = refAtomsRea.at(i);

		auto d = distance(a1, a2);
		d_sum += d * d;

		std::cout << a1 << ": " << a1.get_location() << " => " << a2.get_location() << "  distance: " << d << '\n';
	}

	file.save("/tmp/rsr-test-3.cif");

	auto rmsd = std::sqrt(d_sum / atomsRea.size());
	std::cout << "RMSd: " << rmsd << '\n';
	CHECK_THAT(rmsd, Catch::Matchers::WithinAbs(0.35, 0.35 / 2));

	std::cout << std::string(cif::get_terminal_width(), '-') << '\n';
}
