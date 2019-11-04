#include "pdb-redo.h"

#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>

#include "cif++/Cif++.h"
#include "cif++/Compound.h"
#include "cif++/Structure.h"
#include "cif++/Symmetry.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace c = mmcif;

// -----------------------------------------------------------------------

const set<string> kBackBone = {
	"N", "CA", "C", "O", "OXT"
};

const float
	kMaxZnHisDistanceInCluster = 3.8f,
	kMaxZnCysDistanceInCluster = 4.8f,
	kBoundaryClosenessAtomToZn = 2.9f,

	kDixonQTest95Perc5Points = 0.710f,

	kMaxMetalLigandDistance = 3.5f;

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
	size_t writeOctahedral(const c::Atom& ion, const vector<c::Atom>& ligands,
		const vector<tuple<size_t,size_t>>& opposing);

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
			if (not aw.m_a.pdbxAuthAltId().empty())
				os << " alt " + aw.m_a.pdbxAuthAltId();
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

size_t RestraintGenerator::writeOctahedral(const c::Atom& ion, const vector<c::Atom>& ligands,
	const vector<tuple<size_t,size_t>>& opposing)
{
	const float kSD = 3;

	size_t n = 0, ai = 0;
	for (auto a = ligands.begin(); next(a) != ligands.end(); ++a, ++ai)
	{
		size_t bi = ai + 1;
		for (auto b = a + 1; b != ligands.end(); ++b, ++bi)
		{
			if (find(opposing.begin(), opposing.end(), make_tuple(ai, bi)) == opposing.end())
				writeAngleRestraint(90, kSD, *a, ion, *b);
			else
				writeAngleRestraint(180, kSD, *a, ion, *b);
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

struct IonSite
{
	c::Atom ion;
	vector<tuple<c::Atom,float,string>> lig;
	vector<tuple<size_t,size_t>> opposing;

	bool isOctaHedral();
};

bool IonSite::isOctaHedral()
{
	const float kMaxAllowedAngleDeviation = 30.0f;

	bool result = lig.size() == 6;

	// check angles
	for (size_t a = 0; result and a + 1 < 6; ++a)
	{
		auto& la = get<0>(lig[a]);
		size_t opposing_la = 0;

		for (size_t b = a + 1; result and b < 6; ++b)
		{
			auto& lb = get<0>(lig[b]);
			float angle = c::Angle(la.location(), ion.location(), lb.location());
			
			if (abs(angle - 180) < kMaxAllowedAngleDeviation)	// opposing?
			{
				if (opposing_la++ > 0)
					result = false;
				else
					opposing.emplace_back(a, b);
			}
			else if (abs(angle - 90) > kMaxAllowedAngleDeviation) // should be 90 degrees then
				result = false;
		}
	}

	return result and opposing.size() == 3;
}

// -----------------------------------------------------------------------

vector<IonSite> findZincSites(c::Structure& structure, cif::Datablock& db, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
{
	vector<IonSite> result;

    // factory for symmetry atom iterators
	c::SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom: structure.atoms())
	{
	 	if (atom.labelCompId() != "ZN")
		 	continue;

		IonSite zs = { atom };

		for (auto a: structure.atoms())
		{
			if (a.labelCompId() == "HIS" and (a.labelAtomId() == "ND1" or a.labelAtomId() == "NE2"))
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxZnHisDistanceInCluster)
						zs.lig.emplace_back(sa, d, saif.symop_mmcif(sa));
				}
				continue;
			}

			if (a.labelCompId() == "CYS" and a.labelAtomId() == "SG")
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxZnCysDistanceInCluster)
						zs.lig.emplace_back(sa, d, saif.symop_mmcif(sa));
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
				 << " zn: " << zs.ion.labelAsymId() << '/' << zs.ion.labelAtomId() << " (" << zs.ion.pdbID() << ')' << endl;

			for (auto& l: zs.lig)
			{
				auto& a = get<0>(l);
				cerr << " " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << " (" << a.pdbID() << ')'  << ' ' << saif.symop_mmcif(a) << " @ " << get<1>(l) << endl;
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
                    cerr << "Atom " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << " (" << a.pdbID() << ')'  << " was considered to be an outlier" << endl;
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

// -----------------------------------------------------------------------

constexpr float get_t_90(size_t N)
{
	const float t_dist_90[4] = {
		1.638,	// 3
		1.533,	// 4
		1.476,	// 5
		1.440,	// 6
	};

	assert(N > 5 and N < 10);
	return t_dist_90[N - 3];
}

vector<IonSite> findOctahedralSites(c::Structure& structure, cif::Datablock& db, const clipper::Spacegroup& spacegroup, const clipper::Cell& cell)
{
	vector<IonSite> result;

    // factory for symmetry atom iterators
	c::SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom: structure.atoms())
	{
	 	if (atom.labelCompId() != "NA" and atom.labelCompId() != "MG")
		 	continue;

		IonSite is = { atom };

		for (auto a: structure.atoms())
		{
			if (a.type() == c::AtomType::S or
				a.type() == c::AtomType::O or
				a.type() == c::AtomType::N)
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxMetalLigandDistance)
						is.lig.emplace_back(sa, d, saif.symop_mmcif(sa));
				}
			}
		}

		// sort ligands on distance
		sort(is.lig.begin(), is.lig.end(), [](auto& a, auto& b) { return get<1>(a) < get<1>(b); });

		// take the nearest atom of a residue, this takes care of NE2/ND1 from a single HIS
		// and also the alt cases
		for (auto a = is.lig.begin(); a != is.lig.end() and next(a) != is.lig.end(); ++a)
		{
			auto& aa = get<0>(*a);

			auto ad = get<1>(*a);

			for (auto b = next(a); b != is.lig.end(); ++b)
			{
				auto& ba = get<0>(*b);

				// mmCIF...
				if (aa.isWater())
				{
					if (aa.authSeqId() != ba.authSeqId() or aa.authAsymId() != ba.authAsymId() or aa.symmetry() != ba.symmetry())
						continue;
				}
				else if (ba.labelCompId() != aa.labelCompId() or aa.labelSeqId() != ba.labelSeqId() or aa.labelAsymId() != ba.labelAsymId() or aa.symmetry() != ba.symmetry())
					continue;

				auto bd = get<1>(*b);
				assert(bd > ad);

				is.lig.erase(b);
				break;
			}
		}

		// -----------------------------------------------------------------------

		if (cif::VERBOSE)
		{
			cerr << "preliminary cluster: " << endl
				 << " metal: " << is.ion.labelAsymId() << '/' << is.ion.labelAtomId() << " (" << is.ion.pdbID() << ')' << endl;

			for (auto& l: is.lig)
			{
				auto& a = get<0>(l);
				cerr << " " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << " (" << a.pdbID() << ')' << ' ' << saif.symop_mmcif(a) << " @ " << get<1>(l) << endl;
			}
		}

		// -----------------------------------------------------------------------
		
		for (auto a = is.lig.begin(); a != is.lig.end() and next(a) != is.lig.end(); ++a)
		{
			auto& aa = get<0>(*a);

			if (aa.type() != c::AtomType::N)
				continue;
				
			if (aa.labelCompId() == "GLN")
			{
				assert(aa.labelAtomId() == "NE2");
				auto o = structure.getAtomByLabel("OE1", aa.labelAsymId(), "GLN", aa.labelSeqId(), aa.labelAltId());

				if (cif::VERBOSE)
					cerr << "Flipping side chain for GLN " << " " << aa.labelAsymId() << aa.labelSeqId() << " (" << aa.pdbID() << ')' << endl;

				structure.swapAtoms(aa, o);
				continue;
			}
			
			if (aa.labelCompId() == "ASN")
			{
				assert(aa.labelAtomId() == "ND2");
				auto o = structure.getAtomByLabel("OD1", aa.labelAsymId(), "GLN", aa.labelSeqId(), aa.labelAltId());

				if (cif::VERBOSE)
					cerr << "Flipping side chain for ASN " << " " << aa.labelAsymId() << aa.labelSeqId() << " (" << aa.pdbID() << ')' << endl;

				structure.swapAtoms(aa, o);
				continue;
			}
		}

		// -----------------------------------------------------------------------
		// if there are less than six atoms or more than nine, give up.

		if (is.lig.size() < 6 or is.lig.size() > 9)
		{
			if (cif::VERBOSE)
				cerr << "Rejecting cluster since the number of atoms is unreasonable" << endl;
			continue;
		}
		
		// However, if we have more than six, try to remove outliers with a Grubbs test
		// See: https://en.wikipedia.org/wiki/Grubbs%27s_test_for_outliers

		while (is.lig.size() > 6)
		{
			double sum = accumulate(is.lig.begin(), is.lig.end(), 0.0, [](double s, auto& l) { return s + get<1>(l); });
			double avg = sum / is.lig.size();
			double stddev = sqrt(accumulate(is.lig.begin(), is.lig.end(), 0.0, [avg](double s, auto& l) { return s + (get<1>(l) - avg) * (get<1>(l) - avg); }) / (is.lig.size() - 1));

			// only test if max distance is outlier
			double G = (get<1>(is.lig.back()) - avg) / stddev;

			// extracted from student t-distribution table, with one-sided confidence level 90%
			// and degrees of freedom 
			
			if (G < get_t_90(is.lig.size()))
				break;
			
			if (cif::VERBOSE)
			{
				auto& a = get<0>(is.lig.back());
				cerr << "Removing outlier " << a.labelAsymId() << a.labelSeqId() << '/' << a.labelAtomId() << " (" << a.pdbID() << ')' << ' ' << saif.symop_mmcif(a) << " @ " << get<1>(is.lig.back()) << endl;
			}

			is.lig.erase(is.lig.begin() + is.lig.size() - 1);
		}

		if (not is.isOctaHedral())
		{
			if (cif::VERBOSE)
				cerr << "Rejecting cluster since it is not an octahedral" << endl;
			continue;
		}

		result.emplace_back(move(is));
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

	// -----------------------------------------------------------------------

	string restraintFileName = entryId + "_platonyzer.restraints";
	if (vm.count("restraints-file"))
		restraintFileName = vm["restraints-file"].as<string>();
	RestraintGenerator rg(restraintFileName);

	auto& structConn = db["struct_conn"];

	size_t removedLinks = 0;
	size_t platonyzerLinkId = 1;

	for (auto& ionSites: {
		findZincSites(structure, db, clipper::Spacegroup(clipper::Spgr_descr(spacegroup)), cell),
		findOctahedralSites(structure, db, clipper::Spacegroup(clipper::Spgr_descr(spacegroup)), cell)
		})
	{
		for (auto ionSite: ionSites)
		{
			// write restraints
			vector<c::Atom> ligands;
			transform(ionSite.lig.begin(), ionSite.lig.end(), back_inserter(ligands), [](auto& l) { return get<0>(l); });

			switch (ionSite.lig.size())
			{
				case 4: 
					rg.writeTetrahedral(ionSite.ion, ligands);
					break;
				case 6:
					rg.writeOctahedral(ionSite.ion, ligands, ionSite.opposing);
					break;
			};

			// replace LINK/struct_conn records
			size_t n = structConn.size();
			structConn.erase(
				(cif::Key("ptnr1_label_asym_id") == ionSite.ion.labelAsymId() and cif::Key("ptnr1_label_atom_id") == ionSite.ion.labelAtomId()) or
				(cif::Key("ptnr2_label_asym_id") == ionSite.ion.labelAsymId() and cif::Key("ptnr2_label_atom_id") == ionSite.ion.labelAtomId()));
			removedLinks += (n - structConn.size());

			for (auto&& [atom, distance, symop] : ionSite.lig)
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
					{ "ptnr1_symmetry", "1_555" },
					{ "ptnr1_auth_asym_id", atom.authAsymId() },
					{ "ptnr1_auth_comp_id", atom.authCompId() },
					{ "ptnr1_auth_seq_id", atom.authSeqId() },

					{ "ptnr2_label_asym_id", ionSite.ion.labelAsymId() },
					{ "ptnr2_label_comp_id", ionSite.ion.labelCompId() },
					{ "ptnr2_label_seq_id", "?" },
					{ "ptnr2_label_atom_id", ionSite.ion.labelAtomId() },
					{ "pdbx_ptnr2_label_alt_id", ionSite.ion.labelAltId() },
					{ "pdbx_ptnr2_PDB_ins_code", ionSite.ion.pdbxAuthInsCode() },
					{ "ptnr2_auth_asym_id", ionSite.ion.authAsymId() },
					{ "ptnr2_auth_comp_id", ionSite.ion.authCompId() },
					{ "ptnr2_auth_seq_id", ionSite.ion.authSeqId() },
					{ "ptnr2_symmetry", symop },

					{ "pdbx_dist_value", distance }
				});
			}
		}
	}

	// -----------------------------------------------------------------------

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
