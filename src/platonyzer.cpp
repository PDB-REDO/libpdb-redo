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
#include "cif++/DistanceMap.h"

using namespace std;
namespace po = boost::program_options;
// namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
// namespace io = boost::iostreams;
namespace c = mmcif;
// namespace zx = zeep::xml;

// -----------------------------------------------------------------------

const set<string> kBackBone = {
	"N", "CA", "C", "O", "OXT"
};

const float
	kMaxZnHisDistanceInCluster = 3.8f,
	kMaxZnCysDistanceInCluster = 4.8f,
	kBoundaryClosenessAtomToZn = 2.9f,
	
	kDixonQTest95Perc5Points = 0.710f;

// -----------------------------------------------------------------------

struct ZincSite
{
	c::Atom zn;
	vector<tuple<c::Atom,float>> lig;
};

// -----------------------------------------------------------------------

vector<ZincSite> findZincSites(c::Structure& structure, cif::Datablock& db, const c::DistanceMap& distance)
{
	// vector<c::Atom> zn, his, cys, cysCA, SO;

	vector<ZincSite> result;

	for (auto atom: structure.atoms())
	{
	 	if (atom.labelCompId() != "ZN")
		 	continue;

		ZincSite zs = { atom };

		for (auto a: distance.near(atom, 5.0f))
		{
			float d = distance(atom, a);

			if (a.labelCompId() == "HIS" and (a.labelAtomId() == "ND1" or a.labelAtomId() == "NE2") and
				d <= kMaxZnHisDistanceInCluster)
			{
				zs.lig.emplace_back(a, d);
				continue;
			}

			if (a.labelCompId() == "CYS" and a.labelAtomId() == "SG" and
				d <= kMaxZnCysDistanceInCluster)
			{
				zs.lig.emplace_back(a, d);
				continue;
			}
		}

		// sort ligands on distance
		sort(zs.lig.begin(), zs.lig.end(), [](auto& a, auto& b) { return get<1>(a) < get<1>(b); });

		// take the nearest atom of a residue, this takes care of NE2/ND1 from a single HIS
		// and also the alt cases
		for (auto a = zs.lig.begin(); a != zs.lig.end() and next(a) != zs.lig.end(); ++a)
		{
			auto& aa = get<0>(*a);

			if (aa.labelCompId() != "HIS")
				continue;

			auto ad = get<1>(*a);

			for (auto b = next(a); b != zs.lig.end(); ++b)
			{
				auto& ba = get<0>(*b);

				if (ba.labelCompId() != aa.labelCompId() or aa.labelSeqId() != ba.labelSeqId() or aa.labelAsymId() != ba.labelAsymId())
					continue;

				auto bd = get<1>(*b);
				assert(bd > ad);

				zs.lig.erase(b);
				break;
			}
		}

		// -----------------------------------------------------------------------
		
		if (cif::VERBOSE)
		{
			cerr << "preliminary cluster: " << endl
				 << " zn: " << zs.zn.id() << endl;

			for (auto& l: zs.lig)
			{
				auto& a = get<0>(l);
				cerr << " " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << endl;
			}
		}

		// -----------------------------------------------------------------------
		// if there are more than five atoms, give up. If there are exactly five
		// see if we can use a dixon q-test to assign the last to be an outlier.

		if (zs.lig.size() > 5)
		{
			if (cif::VERBOSE)
				cerr << "Rejecting cluster since there are more than 5 near atoms" << endl;

			continue;
		}

		if (zs.lig.size() == 5)
		{
			auto gap = get<1>(zs.lig[4]) - get<1>(zs.lig[3]);
			auto range = get<1>(zs.lig[4]) - get<1>(zs.lig[0]);
			if ((gap / range) < kDixonQTest95Perc5Points)
			{
				if (cif::VERBOSE)
					cerr << "Rejecting cluster since there are 5 atoms near by and none is considered to be an outlier" << endl;
				continue;
			}

			if (cif::VERBOSE)
			{
				auto& a = get<0>(zs.lig[4]);
				cerr << "Atom " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << " was considered to be an outlier" << endl;
			}

			zs.lig.erase(prev(zs.lig.end()));
		}

	}

	return result;
}


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
	c::Structure structure(pdb);
	
	auto& db = pdb.data();

	// -----------------------------------------------------------------------
	// first build a distance map, makes life easier

	string entryId = db["entry"].front()["id"].as<string>();
	if (entryId.empty())
		throw runtime_error("Missing _entry.id in coordinates file");
	
	double a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = db["cell"][cif::Key("entry_id") == entryId]
		.get("length_a", "length_b", "length_c",
			 "angle_alpha", "angle_beta", "angle_gamma");
	
	clipper::Cell cell(clipper::Cell_descr(a, b, c, alpha, beta, gamma));

	string spacegroup = db["symmetry"]
		[cif::Key("entry_id") == entryId]
		["space_group_name_H-M"].as<string>();
	
	if (spacegroup == "P 1-")
		spacegroup = "P -1";
	else if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw runtime_error("No spacegroup, cannot continue");
	
	c::DistanceMap dm(structure, clipper::Spacegroup(clipper::Spgr_descr(spacegroup)), cell, 5.0f);

	// -----------------------------------------------------------------------
	
	findZincSites(structure, db, dm);

	// -----------------------------------------------------------------------
	
	db.add_software("platonyzer", "other", get_version_nr(), get_version_date());

	if (vm.count("output"))
		pdb.save(vm["output"].as<string>());
	else
		pdb.file().save(cout);

	return result;
}
