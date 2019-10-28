#include "pdb-redo.h"

// #include <sys/wait.h>

// #include <fstream>
// #include <chrono>

#include <boost/program_options.hpp>
// #include <boost/algorithm/string.hpp>
// #include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
// #include <boost/iostreams/filter/bzip2.hpp>
// #include <boost/iostreams/filter/gzip.hpp>
// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/device/file_descriptor.hpp>
// #include <boost/iostreams/copy.hpp>

// #include <zeep/xml/document.hpp>

#include "cif++/Cif++.h"
#include "cif++/Compound.h"
#include "cif++/Structure.h"

using namespace std;
namespace po = boost::program_options;
// namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
// namespace io = boost::iostreams;
namespace c = mmcif;
// namespace zx = zeep::xml;

int pr_main(int argc, char* argv[])
{
	int result = 0;
	
	po::options_description visible_options("platonyzer " + VERSION + " options file]" );
	visible_options.add_options()
		("output,o",	po::value<string>(),	"The output file, default is stdout")
		("help,h",								"Display help message")
		("version",								"Print version")
		("verbose,v",							"Verbose output")
		// ("pdb-redo-data", po::value<string>(),	"The PDB-REDO dat file" /*, default is the built in one"*/)
		("dict",		po::value<string>(),	"Dictionary file containing restraints for residues in this specific target")
		;

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("input,i",		po::value<string>(),	"Input files")
		("test",								"Run test-suite")
		("debug,d",		po::value<int>(),		"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 2);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "platonyzer.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "platonyzer.conf";
	
	if (fs::exists(configFile))
	{
		fs::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}
	
	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION << endl;
		exit(0);
	}

	if (vm.count("help") or vm.count("input") == 0)
	{
		cerr << visible_options << endl;
		exit(1);
	}

	// Load dict, if any
	
	if (vm.count("dict"))
		c::CompoundFactory::instance().pushDictionary(vm["dict"].as<string>());

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();
	
	fs::path input = vm["input"].as<string>();
	c::File pdb(input);
	
	auto& db = pdb.data();

	ofstream output;
	streambuf* buf = nullptr;
	
	if (vm.count("output"))
	{
		output.open(vm["output"].as<string>());
		if (not output.is_open())
			throw runtime_error("Could not open output file " + vm["output"].as<string>());

		buf = cout.rdbuf(output.rdbuf());
	}
	
		 
	db.add_software("platonyzer", "other", get_version_nr(), get_version_date());
	pdb.save(vm["output"].as<string>());
	
	if (buf != nullptr)
	{
		output.close();
		cout.rdbuf(buf);
	}

	return result;
}
