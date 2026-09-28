// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2024
// SPDX-License-Identifier: BSD-2-Clause

#define CATCH_CONFIG_RUNNER

#include <pdb-redo/Compound.hpp>

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cif++/cif++.hpp>
#include <filesystem>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

std::filesystem::path gTestDir;

int main(int argc, char *argv[])
{
	gTestDir = std::filesystem::current_path();

	cif::VERBOSE = 1;

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

TEST_CASE("een")
{
	auto &cf = pdb_redo::CompoundFactory::instance();

	for (fs::directory_iterator i(gTestDir / "restr2ccd"); i != fs::directory_iterator(); ++i)
	{
		auto fn = i->path().filename().string();

		if (not fn.ends_with(".cif"))
			continue;

		auto comp_id = i->path().filename().stem().string();

		cf.pushDictionary(i->path());

		auto compound = cf.create(comp_id);

		REQUIRE(compound != nullptr);
		CHECK(cif::iequals(compound->id(), comp_id));

		auto ccompound = cif::compound_factory::instance().create(comp_id);

		REQUIRE(ccompound != nullptr);
		CHECK(ccompound->id() == comp_id);

		cf.popDictionary();
	}
}