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

#include "pdb-redo.hpp"

#include <iomanip>
#include <fstream>
#include <filesystem>

#include <boost/program_options.hpp>

#include "cif++/Cif++.hpp"
#include "cif++/Compound.hpp"
#include "cif++/Structure.hpp"
#include "cif++/Symmetry.hpp"

#include "Symmetry-2.hpp"

using namespace std;
namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace c = mmcif;

// -----------------------------------------------------------------------

const set<string> kBackBone = {
	"N", "CA", "C", "O", "OXT"
};

const float
	kMaxZnHisDistanceInCluster = 3.8f,
	kMaxZnCysDistanceInCluster = 4.8f,
	// kBoundaryClosenessAtomToZn = 2.9f,

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

	RestraintGenerator(const string& file, bool deleteVDWRestraints)
		: m_file(file), m_deleteVDWRestraints(deleteVDWRestraints)
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
			os << "chain " << aw.m_a.authAsymID()
			   << " resi " << aw.m_a.authSeqID()
			   << " ins " << (aw.m_a.pdbxAuthInsCode().empty() ? "." : aw.m_a.pdbxAuthInsCode())
			   << " atom " << aw.m_a.authAtomID();
			if (not aw.m_a.pdbxAuthAltID().empty())
				os << " alt " + aw.m_a.pdbxAuthAltID();
			os << " symm " << (aw.m_a.isSymmetryCopy() ? "Y" : "N");

			return os;
		}
	};

	ofstream m_file;
	bool m_deleteVDWRestraints;
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

	if (m_deleteVDWRestraints)
		m_file << "vdwr exclude " << AtomPart(ion) << endl;

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

vector<IonSite> findZincSites(c::Structure& structure, cif::Datablock& db, int spacegroup, const clipper::Cell& cell)
{
	vector<IonSite> result;

    // factory for symmetry atom iterators
	SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom: structure.atoms())
	{
	 	if (atom.labelCompID() != "ZN")
		 	continue;

		IonSite zs = { atom };

		for (auto a: structure.atoms())
		{
			if (a.labelCompID() == "HIS" and (a.labelAtomID() == "ND1" or a.labelAtomID() == "NE2"))
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxZnHisDistanceInCluster)
						zs.lig.emplace_back(sa, d, sa.symmetry());
				}
				continue;
			}

			if (a.labelCompID() == "CYS" and a.labelAtomID() == "SG")
			{
				for (auto sa: saif(a))
				{
					float d = Distance(atom, sa);
					if (d <= kMaxZnCysDistanceInCluster)
						zs.lig.emplace_back(sa, d, sa.symmetry());
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

				if (ba.labelCompID() != aa.labelCompID() or aa.labelSeqID() != ba.labelSeqID() or aa.labelAsymID() != ba.labelAsymID() or aa.symmetry() != ba.symmetry())
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
				 << " zn: " << zs.ion.labelAsymID() << '/' << zs.ion.labelAtomID() << " (" << zs.ion.pdbID() << ')' << endl;

			for (auto& l: zs.lig)
			{
				auto& a = get<0>(l);
				cerr << " " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')'  << ' ' << a.symmetry() << " @ " << get<1>(l) << endl;
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
                    cerr << "Atom " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')'  << " was considered to be an outlier" << endl;
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
	const float t_dist_90[] = {
		1.638,	// 3
		1.533,	// 4
		1.476,	// 5
		1.440,	// 6
		1.415,	// 7
		1.397,	// 8
		1.383,	// 9
		1.372,	// 10
	};

	assert(N >= 3 and N < sizeof(t_dist_90) / sizeof(float) + 3);
	return t_dist_90[N - 3];
}

vector<IonSite> findOctahedralSites(c::Structure& structure, cif::Datablock& db, int spacegroup, const clipper::Cell& cell)
{
	vector<IonSite> result;

    // factory for symmetry atom iterators
	SymmetryAtomIteratorFactory saif(structure, spacegroup, cell);

	for (auto atom: structure.atoms())
	{
	 	if (atom.labelCompID() != "NA" and atom.labelCompID() != "MG")
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
						is.lig.emplace_back(sa, d, sa.symmetry());
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
					if (aa.authSeqID() != ba.authSeqID() or aa.authAsymID() != ba.authAsymID() or aa.symmetry() != ba.symmetry())
						continue;
				}
				else if (ba.labelAtomID() != aa.labelAtomID() or ba.labelCompID() != aa.labelCompID() or aa.labelSeqID() != ba.labelSeqID() or aa.labelAsymID() != ba.labelAsymID() or aa.symmetry() != ba.symmetry())
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
				 << " metal: " << is.ion.labelAsymID() << '/' << is.ion.labelAtomID() << " (" << is.ion.pdbID() << ')' << endl;

			for (auto& l: is.lig)
			{
				auto& a = get<0>(l);
				cerr << " " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')' << ' ' << a.symmetry() << " @ " << get<1>(l) << endl;
			}
		}

		// -----------------------------------------------------------------------
		
		for (auto a = is.lig.begin(); a != is.lig.end() and next(a) != is.lig.end(); ++a)
		{
			auto& aa = get<0>(*a);

			if (aa.type() != c::AtomType::N)
				continue;
				
			try
			{
				if (aa.labelCompID() == "GLN")
				{
					assert(aa.labelAtomID() == "NE2");
					auto o = structure.getAtomByLabel("OE1", aa.labelAsymID(), "GLN", aa.labelSeqID(), aa.labelAltID());

					if (cif::VERBOSE)
						cerr << "Flipping side chain for GLN " << " " << aa.labelAsymID() << aa.labelSeqID() << " (" << aa.pdbID() << ')' << endl;

					structure.swapAtoms(aa, o);
					continue;
				}
				
				if (aa.labelCompID() == "ASN")
				{
					assert(aa.labelAtomID() == "ND2");
					auto o = structure.getAtomByLabel("OD1", aa.labelAsymID(), "GLN", aa.labelSeqID(), aa.labelAltID());

					if (cif::VERBOSE)
						cerr << "Flipping side chain for ASN " << " " << aa.labelAsymID() << aa.labelSeqID() << " (" << aa.pdbID() << ')' << endl;

					structure.swapAtoms(aa, o);
					continue;
				}
			}
			catch (const std::out_of_range& ex)
			{
				if (cif::VERBOSE)
					cerr << "Could not flip " << aa.labelCompID() << ": " << ex.what() << endl;

				// is.lig.clear();	// give up
				// break;
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
				cerr << "Removing outlier " << a.labelAsymID() << a.labelSeqID() << '/' << a.labelAtomID() << " (" << a.pdbID() << ')' << ' ' << a.symmetry() << " @ " << get<1>(is.lig.back()) << endl;
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

	po::options_description visible_options("platonyzer " + VERSION_STRING + " options file]" );
	visible_options.add_options()
		("output,o",	po::value<string>(),	"The output file, default is stdout")
		("restraints-file,r",
						po::value<string>(),	"Restraint file name, default is {id}_platonyze.restraints")
		("delete-vdw-rest",						"Delete vanderWaals restraints for octahedral ions in the external for Refmac")
		("create-na-mg-links",					"Create links for Na/Mg ion sites that were found")
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
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}

	po::notify(vm);

	if (vm.count("version"))
	{
		cout << argv[0] << " version " << VERSION_STRING << endl;
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

	if (cif::VERBOSE)
		cerr << "Loading data...";

	fs::path input = vm["input"].as<string>();
	c::File pdb(input);
	c::Structure structure(pdb);

	if (cif::VERBOSE)
		cerr << " done" << endl;

	auto& db = pdb.data();

	// -----------------------------------------------------------------------

	string entryId = db["entry"].front()["id"].as<string>();
	if (entryId.empty())
		throw runtime_error("Missing _entry.id in coordinates file");

	double a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = db["cell"][cif::Key("entry_id") == entryId]
		.get("length_a", "length_b", "length_c",
			 "angle_alpha", "angle_beta", "angle_gamma");

	clipper::Cell cell(clipper::Cell_descr(a, b, c, alpha, beta, gamma));

	string spacegroupName = db["symmetry"]
		[cif::Key("entry_id") == entryId]
		["space_group_name_H-M"].as<string>();

	int spacegroupNr = mmcif::GetSpacegroupNumber(spacegroupName);

	// -----------------------------------------------------------------------

	string restraintFileName = entryId + "_platonyzer.restraints";
	if (vm.count("restraints-file"))
		restraintFileName = vm["restraints-file"].as<string>();
	RestraintGenerator rg(restraintFileName, vm.count("delete-vdw-rest"));

	auto& structConn = db["struct_conn"];

	size_t removedLinks = 0, createdLinks = 0;
	size_t platonyzerLinkId = 1;

	bool createNaMgLinks = vm.count("create-na-mg-links");

	for (auto& ionSites: {
		findZincSites(structure, db, spacegroupNr, cell),
		findOctahedralSites(structure, db, spacegroupNr, cell)
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
				(cif::Key("ptnr1_label_asym_id") == ionSite.ion.labelAsymID() and cif::Key("ptnr1_label_atom_id") == ionSite.ion.labelAtomID()) or
				(cif::Key("ptnr2_label_asym_id") == ionSite.ion.labelAsymID() and cif::Key("ptnr2_label_atom_id") == ionSite.ion.labelAtomID()));
			removedLinks += (n - structConn.size());

			if (not createNaMgLinks and (ionSite.ion.type() == c::AtomType::Na or ionSite.ion.type() == c::AtomType::Mg))
				continue;

			createdLinks += ionSite.lig.size();

			for (auto&& [atom, distance, symop] : ionSite.lig)
			{
				string id = "metalc_p-" + to_string(platonyzerLinkId++);

				structConn.emplace({
					{ "id", id },
					{ "conn_type_id", "metalc" },
					{ "ptnr1_label_asym_id", atom.labelAsymID() },
					{ "ptnr1_label_comp_id", atom.labelCompID() },
					{ "ptnr1_label_seq_id", atom.labelSeqID() },
					{ "ptnr1_label_atom_id", atom.labelAtomID() },
					{ "pdbx_ptnr1_label_alt_id", atom.labelAltID() },
					{ "pdbx_ptnr1_PDB_ins_code", atom.pdbxAuthInsCode() },
					{ "ptnr1_symmetry", "1_555" },
					{ "ptnr1_auth_asym_id", atom.authAsymID() },
					{ "ptnr1_auth_comp_id", atom.authCompID() },
					{ "ptnr1_auth_seq_id", atom.authSeqID() },

					{ "ptnr2_label_asym_id", ionSite.ion.labelAsymID() },
					{ "ptnr2_label_comp_id", ionSite.ion.labelCompID() },
					{ "ptnr2_label_seq_id", "?" },
					{ "ptnr2_label_atom_id", ionSite.ion.labelAtomID() },
					{ "pdbx_ptnr2_label_alt_id", ionSite.ion.labelAltID() },
					{ "pdbx_ptnr2_PDB_ins_code", ionSite.ion.pdbxAuthInsCode() },
					{ "ptnr2_auth_asym_id", ionSite.ion.authAsymID() },
					{ "ptnr2_auth_comp_id", ionSite.ion.authCompID() },
					{ "ptnr2_auth_seq_id", ionSite.ion.authSeqID() },
					{ "ptnr2_symmetry", symop },

					{ "pdbx_dist_value", distance }
				});
			}
		}
	}

	// -----------------------------------------------------------------------

	if (cif::VERBOSE)
		cerr << "Removed " << removedLinks << " link records" << endl
			 << "Created " << createdLinks << " link records" << endl;

	// -----------------------------------------------------------------------

	db.add_software("platonyzer", "other", get_version_nr(), get_version_date());

	if (vm.count("output"))
		pdb.save(vm["output"].as<string>());
	else
		pdb.file().save(cout);

	return result;
}
