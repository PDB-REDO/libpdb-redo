#include "pdb-redo.h"

// #include <sys/wait.h>

// #include <fstream>
// #include <chrono>
#include <iomanip>

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

// enum class RestraintType { Dist, Angle }; //TODO implement: TORSION (and others?) when necessary

// string to_string(RestraintType type)
// {
// 	switch (type)
// 	{
// 		case RestraintType::Angle:	return "angle";
// 		case RestraintType::Dist:	return "dist";
// 	}
// }

class RestraintGenerator
{
  public:

	RestraintGenerator(const string& file)
		: m_file(file)
	{
		if (not m_file.is_open())
			throw runtime_error("Could not open restraint file " + file);
	}

	~RestraintGenerator() = default;

	size_t writeTetrahedral(const c::Atom& ion, const vector<c::Atom>& ligands);

  private:

	void writeAngleRestraint(float target, float sd, const c::Atom& a, const c::Atom& b, const c::Atom& c);

	struct AtomPart
	{
		const c::Atom& m_a;

		AtomPart(const c::Atom& a) : m_a(a) {}

		friend ostream& operator<<(ostream& os, const AtomPart& aw)
		{
			os << "chain " << aw.m_a.authAsymId()
			   << " resi " << aw.m_a.authSeqId()
			   << " ins " << (aw.m_a.pdbxAuthInsCode().empty() ? "." : aw.m_a.pdbxAuthInsCode())
			   << " atom " << aw.m_a.authAtomId();
			if (not aw.m_a.authAltId().empty())
				os << " alt " + aw.m_a.authAltId();
			os << " symm " << (aw.m_a.isSymmetryCopy() ? "Y" : "N");

			return os;
		}
	};

	ofstream m_file;
};

size_t RestraintGenerator::writeTetrahedral(const c::Atom& ion, const vector<c::Atom>& ligands)
{
	const float
		kTarget = 109.5,
		kSD = 3;

	size_t n = 0;
	for (auto a = ligands.begin(); next(a) != ligands.end(); ++a)
	{
		for (auto b = a + 1; b != ligands.end(); ++b)
		{
			writeAngleRestraint(kTarget, kSD, *a, ion, *b);
			++n;
		}
	}

	return n;
}

void RestraintGenerator::writeAngleRestraint(float target, float sd, const c::Atom& a, const c::Atom& b, const c::Atom& c)
{
	m_file << "exte angle first " << AtomPart(a) << " next " << AtomPart(b) << " next " << AtomPart(c)
		   << fixed << setprecision(3) << " value " << target << " sigma " << sd << endl;
}

// ------------------------------------------------------------------------

struct ZincSite
{
	c::Atom zn;
	vector<tuple<c::Atom,float>> lig;
};

// -----------------------------------------------------------------------

vector<ZincSite> findZincSites(c::Structure& structure, cif::Datablock& db, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
{
	vector<ZincSite> result;

    // factory for symmetry atom iterators
	c::SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom: structure.atoms())
	{
	 	if (atom.labelCompId() != "ZN")
		 	continue;

		ZincSite zs = { atom };

		for (auto a: structure.atoms())
		{
			if (a.labelCompId() == "HIS" and (a.labelAtomId() == "ND1" or a.labelAtomId() == "NE2"))
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxZnHisDistanceInCluster)
						zs.lig.emplace_back(sa, d);
				}
				continue;
			}

			if (a.labelCompId() == "CYS" and a.labelAtomId() == "SG")
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxZnCysDistanceInCluster)
						zs.lig.emplace_back(sa, d);
				}
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

			auto ad = get<1>(*a);

			for (auto b = next(a); b != zs.lig.end(); ++b)
			{
				auto& ba = get<0>(*b);

				if (ba.labelCompId() != aa.labelCompId() or aa.labelSeqId() != ba.labelSeqId() or aa.labelAsymId() != ba.labelAsymId() or aa.symmetry() != ba.symmetry())
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
				 << " zn: " << zs.zn.labelAsymId() << '/' << zs.zn.labelAtomId() << endl;

			for (auto& l: zs.lig)
			{
				auto& a = get<0>(l);
				cerr << " " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << ' ' << a.symmetry() << " @ " << get<1>(l) << endl;
			}
		}

		// -----------------------------------------------------------------------
		// if there are more than five atoms, give up. If there are exactly five
		// see if we can use a dixon q-test to assign the last to be an outlier.

		if (zs.lig.size() >= 5)
		{
			auto gap = get<1>(zs.lig[4]) - get<1>(zs.lig[3]);
			auto range = get<1>(zs.lig[4]) - get<1>(zs.lig[0]);
			if ((gap / range) < kDixonQTest95Perc5Points)
			{
				if (cif::VERBOSE)
					cerr << "Rejecting cluster since there are 5 or more atoms near by and nr 5 is not considered to be an outlier" << endl;
				continue;
			}

			if (cif::VERBOSE)
			{
                for (size_t i = 4; i < zs.lig.size(); ++i)
                {
                    auto& a = get<0>(zs.lig[i]);
                    cerr << "Atom " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << " was considered to be an outlier" << endl;
                }
			}

			zs.lig.erase(zs.lig.begin() + 4, zs.lig.end());
		}

		if (zs.lig.size() != 4)
		{
			if (cif::VERBOSE)
				cerr << "Rejecting cluster since there not 4 atoms near by" << endl;
			continue;
		}

		result.emplace_back(move(zs));
	}

	return result;
}

int pr_main(int argc, char* argv[])
{
	int result = 0;

	po::options_description visible_options("platonyzer " + VERSION + " options file]" );
	visible_options.add_options()
		("output,o",	po::value<string>(),	"The output file, default is stdout")
		("restraints-file,r",
						po::value<string>(),	"Restraint file name, default is {id}_platonyze.restraints")
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

	// c::DistanceMap dm(structure, clipper::Spacegroup(clipper::Spgr_descr(spacegroup)), cell, 5.0f);

	// -----------------------------------------------------------------------

	auto zincSites = findZincSites(structure, db, clipper::Spacegroup(clipper::Spgr_descr(spacegroup)), cell);

	string restraintFileName = entryId + "_platonyzer.restraints";
	if (vm.count("restraints-file"))
		restraintFileName = vm["restraints-file"].as<string>();

	auto& structConn = db["struct_conn"];

	size_t removedLinks = 0;
	size_t platonyzerLinkId = 1;

	RestraintGenerator rg(restraintFileName);
	for (auto zs: zincSites)
	{
		// write restraints
		vector<c::Atom> ligands;
		transform(zs.lig.begin(), zs.lig.end(), back_inserter(ligands), [](auto& l) { return get<0>(l); });
		rg.writeTetrahedral(zs.zn, ligands);

		// replace LINK/struct_conn records
		size_t n = structConn.size();
		structConn.erase(
			(cif::Key("ptnr1_label_asym_id") == zs.zn.labelAsymId() and cif::Key("ptnr1_label_atom_id") == zs.zn.labelAtomId()) or
			(cif::Key("ptnr2_label_asym_id") == zs.zn.labelAsymId() and cif::Key("ptnr2_label_atom_id") == zs.zn.labelAtomId()) or
			(cif::Key("ptnr3_label_asym_id") == zs.zn.labelAsymId() and cif::Key("ptnr3_label_atom_id") == zs.zn.labelAtomId()));
		removedLinks += (n - structConn.size());

		for (auto&& [atom, distance] : zs.lig)
		{
			string id = "metalc_p-" + to_string(platonyzerLinkId++);

			structConn.emplace({
				{ "id", id },
				{ "conn_type_id", "metalc" },
				{ "ptnr1_label_asym_id", atom.labelAsymId() },
				{ "ptnr1_label_comp_id", atom.labelCompId() },
				{ "ptnr1_label_seq_id", atom.labelSeqId() },
				{ "ptnr1_label_atom_id", atom.labelAtomId() },
				{ "pdbx_ptnr1_label_alt_id", atom.labelAltId() },
				{ "pdbx_ptnr1_PDB_ins_code", atom.pdbxAuthInsCode() },
				// { "ptnr1_symmetry", atom.symmetry() },
				{ "ptnr1_auth_asym_id", atom.authAsymId() },
				{ "ptnr1_auth_comp_id", atom.authCompId() },
				{ "ptnr1_auth_seq_id", atom.authSeqId() },

				{ "ptnr2_label_asym_id", zs.zn.labelAsymId() },
				{ "ptnr2_label_comp_id", zs.zn.labelCompId() },
				{ "ptnr2_label_seq_id", "?" },
				{ "ptnr2_label_atom_id", zs.zn.labelAtomId() },
				{ "pdbx_ptnr2_label_alt_id", zs.zn.labelAltId() },
				{ "pdbx_ptnr2_PDB_ins_code", zs.zn.pdbxAuthInsCode() },
				{ "ptnr2_auth_asym_id", zs.zn.authAsymId() },
				{ "ptnr2_auth_comp_id", zs.zn.authCompId() },
				{ "ptnr2_auth_seq_id", zs.zn.authSeqId() },
				// { "ptnr2_symmetry", zs.zn.symmetry() },

				{ "pdbx_dist_value", distance }
			});
		}
	}

	if (cif::VERBOSE)
		cerr << "Removed " << removedLinks << " link records" << endl;

	// -----------------------------------------------------------------------

	db.add_software("platonyzer", "other", get_version_nr(), get_version_date());

	if (vm.count("output"))
		pdb.save(vm["output"].as<string>());
	else
		pdb.file().save(cout);

	return result;
}
