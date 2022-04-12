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

/*
   Created by: Maarten L. Hekkelman
   Date: dinsdag 22 mei, 2018
*/

#include <filesystem>
#include <future>
#include <iomanip>
#include <random>
#include <regex>

#include <boost/format.hpp>

#include "pdb-redo/Minimizer.hpp"
#include "pdb-redo/Symmetry-2.hpp"

using namespace mmcif;
namespace fs = std::filesystem;

namespace pdb_redo
{

// --------------------------------------------------------------------

const uint32_t kRefSentinel = std::numeric_limits<uint32_t>::max();

const double
	kNonBondedContactDistanceSq = 11.0 * 11.0,
	kMaxPeptideBondLengthSq = 3.5 * 3.5;

// --------------------------------------------------------------------

struct lessAtom
{
	bool operator()(const Atom &a, const Atom &b) const { return a.id().compare(b.id()) < 0; }
};

typedef std::set<Atom, lessAtom> AtomSet;

// --------------------------------------------------------------------

DPoint AtomLocationProvider::operator[](AtomRef atomID) const
{
	if (atomID >= mAtoms.size())
		throw std::range_error("Unknown atom " + std::to_string(atomID));
	return mAtoms[atomID].location();
}

std::string AtomLocationProvider::atom(AtomRef atomID) const
{
	if (atomID >= mAtoms.size())
		throw std::range_error("Unknown atom " + std::to_string(atomID));
	auto &a = mAtoms[atomID];
	return std::to_string(a.labelSeqID()) + ' ' + a.labelAtomID();
}

// --------------------------------------------------------------------

Minimizer::Minimizer(const mmcif::Structure &structure, const mmcif::BondMap &bonds, float plane5AtomsESD)
	: mStructure(structure)
	, mBonds(bonds)
	, mPlane5ESD(plane5AtomsESD)
{
}

void Minimizer::addResidue(const mmcif::Residue &res)
{
	auto compound = Compound::create(res.compoundID()); // r.compound();
	if (not compound)
		throw std::runtime_error("Missing compound information for " + res.compoundID());

	for (auto a : res.atoms())
	{
		(void)ref(a);
		mAtoms.push_back(a);
	}

	for (auto &b : compound->bonds())
	{
		try
		{
			if (compound->getAtomByID(b.atomID[0]).typeSymbol == H or
				compound->getAtomByID(b.atomID[1]).typeSymbol == H)
			{
				continue;
			}

			Atom a1 = res.atomByID(b.atomID[0]);
			Atom a2 = res.atomByID(b.atomID[1]);

			if (not(a1 and a2))
				continue;

			mBondRestraints.emplace_back(ref(a1), ref(a2), b.distance, b.esd);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "While processing bond restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	for (auto &a : compound->angles())
	{
		try
		{
			if (compound->getAtomByID(a.atomID[0]).typeSymbol == H or
				compound->getAtomByID(a.atomID[1]).typeSymbol == H or
				compound->getAtomByID(a.atomID[2]).typeSymbol == H)
			{
				continue;
			}

			Atom a1 = res.atomByID(a.atomID[0]);
			Atom a2 = res.atomByID(a.atomID[1]);
			Atom a3 = res.atomByID(a.atomID[2]);

			if (not(a1 and a2 and a3))
				continue;

			mAngleRestraints.emplace_back(
				ref(a1), ref(a2), ref(a3), a.angle, a.esd);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "While processing angle restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	for (auto &a : compound->torsions())
	{
		if (a.esd == 0)
			continue;

		try
		{
			if (compound->getAtomByID(a.atomID[0]).typeSymbol == H or
				compound->getAtomByID(a.atomID[1]).typeSymbol == H or
				compound->getAtomByID(a.atomID[2]).typeSymbol == H or
				compound->getAtomByID(a.atomID[3]).typeSymbol == H)
			{
				continue;
			}

			Atom a1 = res.atomByID(a.atomID[0]);
			Atom a2 = res.atomByID(a.atomID[1]);
			Atom a3 = res.atomByID(a.atomID[2]);
			Atom a4 = res.atomByID(a.atomID[3]);

			mTorsionRestraints.emplace_back(ref(a1), ref(a2), ref(a3), ref(a4), a.angle, a.esd, a.period);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "While processing torsion restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	for (auto &cv : compound->chiralCentres())
	{
		try
		{
			if (compound->getAtomByID(cv.atomID[0]).typeSymbol == H or
				compound->getAtomByID(cv.atomID[1]).typeSymbol == H or
				compound->getAtomByID(cv.atomID[2]).typeSymbol == H)
			{
				continue;
			}

			Atom cc = res.atomByID(cv.atomIDCentre);
			Atom a1 = res.atomByID(cv.atomID[0]);
			Atom a2 = res.atomByID(cv.atomID[1]);
			Atom a3 = res.atomByID(cv.atomID[2]);

			if (not(cc and a1 and a2 and a3))
				continue;

			auto volume = compound->chiralVolume(cv.id);

			mChiralVolumeRestraints.emplace_back(ref(cc), ref(a1),
				ref(a2), ref(a3), volume * 6);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "While processing chiral volume restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	for (auto &p : compound->planes())
	{
		try
		{
			std::vector<AtomRef> atoms;

			for (auto a : p.atomID)
			{
				if (compound->getAtomByID(a).typeSymbol == H)
					continue;

				auto a1 = res.atomByID(a);
				if (not a1)
					continue;

				atoms.push_back(ref(a1));
			}

			if (atoms.size() > 3)
				mPlanarityRestraints.emplace_back(move(atoms), p.esd);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "While processing planarity restraints: " << ex.what() << std::endl;
			continue;
		}
	}
}

void Minimizer::addPolySection(const mmcif::Polymer &poly, int first, int last)
{
	const Monomer *prev = nullptr; // used to link residues

	for (auto &r : poly)
	{
		if (r.seqID() < first)
		{
			prev = &r;
			continue;
		}

		if (r.seqID() <= last)
			addResidue(r);

		if (prev != nullptr)
			try
			{
				Atom c = prev->atomByID("C"), n = r.atomByID("N");

				if (c and n and DistanceSquared(c, n) < kMaxPeptideBondLengthSq)
				{
					bool trans = not Monomer::isCis(*prev, r);

					if (trans)
						addLinkRestraints(*prev, r, r.compoundID() == "PRO" ? "PTRANS" : "TRANS");
					else
						addLinkRestraints(*prev, r, r.compoundID() == "PRO" ? "PCIS" : "CIS");

					if (trans)
					{
						Atom ca1 = prev->atomByID("CA");
						Atom ca2 = r.atomByID("CA");

						mTransPeptideRestraints.emplace_back(TransPeptideRestraint{
							ref(ca1),
							ref(c),
							ref(n),
							ref(ca2)});
					}

					// add planar restraints
					std::vector<AtomRef> atoms = {
						ref(prev->atomByID("CA")),
						ref(c),
						ref(prev->atomByID("O")),
						ref(n),
						ref(r.atomByID("CA"))};

					mPlanarityRestraints.emplace_back(PlanarityRestraint{move(atoms), mPlane5ESD});
				}
			}
			catch (const std::exception &ex)
			{
				if (cif::VERBOSE > 1)
					std::cerr << "While processing plane-5-atoms restraints: " << ex.what() << std::endl;
				//			continue;
			}

		if (r.seqID() > last)
			break;

		prev = &r;
	}

	// Add link BondRestraints
	//	for (auto& a1: mAtoms)
	//	{
	//		// TODO: implement based on struct_conn information and radii from ener_lib?
	//	}
}

void Minimizer::addDensityMap(const XMap &xMap, float mapWeight)
{
	std::vector<std::pair<AtomRef, double>> densityAtoms;
	densityAtoms.reserve(mAtoms.size());

	transform(mAtoms.begin(), mAtoms.end(), back_inserter(densityAtoms),
		[this](const Atom &a)
		{
			double z = static_cast<int>(a.type());
			double weight = 1;
			double occupancy = a.occupancy();

			if (occupancy > 1)
				occupancy = 1;

			// TODO: cryo_em support

			return std::make_pair(ref(a), z * weight * occupancy);
		});

	mDensityRestraint.reset(new DensityRestraint(move(densityAtoms), xMap, mapWeight));
}

void Minimizer::Finish()
{
	if (mAtoms.empty())
		throw std::runtime_error("No atoms to refine");

	fs::path enerLibFilePath(getenv("CLIBD_MON"));
	enerLibFilePath /= "ener_lib.cif";

	cif::File enerLibFile(enerLibFilePath);
	auto &db = enerLibFile["energy"];
	auto &libAtom = db["lib_atom"];

	const std::regex donorRx("B|D|H"), acceptorRx("B|A|H");

	std::set<std::tuple<AtomRef, AtomRef>> nbc;

	SymmetryAtomIteratorFactory saif(mStructure, getSpacegroup(mStructure), getCell(mStructure));

	// now add the non-bonded restraints
	for (auto &a1 : mAtoms)
	{
		AtomRef ra1 = ref(a1);

		for (auto sAtom : mStructure.atoms())
		{
			for (auto a2 : saif(sAtom, [l = a1.location()](const Point &p) { return DistanceSquared(p, l) <= kNonBondedContactDistanceSq; }))
			{
				if (a1 == a2)
					continue;

				if (mBonds(a1, a2))
					continue;

				if (DistanceSquared(a1, a2) > kNonBondedContactDistanceSq)
					continue;

				AtomRef ra2 = ref(a2);

				if (nbc.count(std::make_tuple(ra1, ra2)))
					continue;

				if (find_if(mAngleRestraints.begin(), mAngleRestraints.end(),
						[&](auto &ar)
						{ return (ar.mA == ra1 and ar.mC == ra2) or (ar.mA == ra2 and ar.mC == ra1); }) != mAngleRestraints.end())
					continue;

				if ((a1.labelCompID() == "PRO" or a1.labelCompID() == "HYP") and
					a1.labelSeqID() == a2.labelSeqID() + 1 and
					a1.labelAtomID() == "CD")
				{
					continue;
				}

				if ((a2.labelCompID() == "PRO" or a2.labelCompID() == "HYP") and
					a2.labelSeqID() == a1.labelSeqID() + 1 and
					a2.labelAtomID() == "CD")
				{
					continue;
				}

				if ((a1.labelCompID() == "ASN" or a2.labelCompID() == "NAG") and
					a1.labelAtomID() == "OD1" and a2.labelAtomID() == "C1")
				{
					continue;
				}

				if ((a1.labelCompID() == "NAG" or a2.labelCompID() == "ASN") and
					a1.labelAtomID() == "C1" and a2.labelAtomID() == "OD1")
				{
					continue;
				}

				double minDist = 2.8;

				if (mBonds.is1_4(a1, a2))
				{
					if (cif::VERBOSE > 1)
						std::cerr << "1_4 for " << a1 << " and " << a2 << std::endl;
					minDist = 2.64;
				}
				else if ((a1.labelSeqID() + 1 == a2.labelSeqID() and a1.labelAtomID() == "O" and a2.labelAtomID() == "C") or
						(a2.labelSeqID() + 1 == a1.labelSeqID() and a2.labelAtomID() == "O" and a1.labelAtomID() == "C"))
				{
					minDist = 2.84;
				}
				else
				{
					try
					{
						auto c1 = Compound::create(a1.labelCompID());
						auto c2 = Compound::create(a2.labelCompID());

						std::string et1 = c1->getAtomByID(a1.labelAtomID()).typeEnergy;
						std::string et2 = c2->getAtomByID(a2.labelAtomID()).typeEnergy;

						if (not(et1.empty() or et2.empty()))
						{
							auto r1 = libAtom.find(cif::Key("type") == et1);
							auto r2 = libAtom.find(cif::Key("type") == et2);

							if (not(r1.empty() or r2.empty()))
							{
								if (AtomTypeTraits(a1.type()).isMetal())
									minDist = r1.front()["ion_radius"].as<float>();
								else
									minDist = r1.front()["vdw_radius"].as<float>();

								if (AtomTypeTraits(a2.type()).isMetal())
									minDist += r2.front()["ion_radius"].as<float>();
								else
									minDist += r2.front()["vdw_radius"].as<float>();

								// OK, now that we're here, see if the atoms are in the same residue...

								if (a1.labelAsymID() == a2.labelAsymID() and a1.labelSeqID() == a2.labelSeqID())
									minDist *= 0.84;

								std::string hbType1 = r1.front()["hb_type"].as<std::string>(),
											hbType2 = r2.front()["hb_type"].as<std::string>();

								if (std::regex_match(hbType1, donorRx) and regex_match(hbType2, acceptorRx))
								{
									minDist -= 0.5;
									if (hbType1 == "H")
										minDist -= 0.3;
								}

								if (regex_match(hbType2, donorRx) and regex_match(hbType1, acceptorRx))
								{
									minDist -= 0.5;
									if (hbType2 == "H")
										minDist -= 0.3;
								}
							}
						}

						// so-called strange exceptions in coot code

						if (find(mAtoms.begin(), mAtoms.end(), a2) == mAtoms.end())
						{
							switch (std::abs(a1.labelSeqID() - a2.labelSeqID()))
							{
								case 1:
									if ((a1.labelAtomID() == "O" and a2.labelAtomID() == "CA") or
										(a1.labelAtomID() == "CA" and a2.labelAtomID() == "O") or
										(a1.labelAtomID() == "N" and a2.labelAtomID() == "CB") or
										(a1.labelAtomID() == "CB" and a2.labelAtomID() == "N") or
										(a1.labelAtomID() == "C" and a2.labelAtomID() == "CB") or
										(a1.labelAtomID() == "CB" and a2.labelAtomID() == "C"))
									{
										minDist = 2.7;
									}
									break;

								case 2:
									if ((a1.labelAtomID() == "C" and a2.labelAtomID() == "N") or
										(a1.labelAtomID() == "N" and a2.labelAtomID() == "C"))
									{
										minDist = 2.7;
									}
									break;
							}
						}
					}
					catch (const std::exception &ex)
					{
						if (cif::VERBOSE)
							std::cerr << "err calculating nbc distance: " << ex.what() << std::endl;
						minDist = 2.8;
					}
				}

				mNonBondedContactRestraints.emplace_back(ra1, ra2, minDist, 0.02);
				nbc.insert(std::make_tuple(ra1, ra2));
				nbc.insert(std::make_tuple(ra2, ra1));
			}
		}
	}

	// create reverse index (for dfcollector)
	mRef2AtomIndex = std::vector<size_t>(mReferencedAtoms.size(), kRefSentinel);
	for (size_t i = 0; i < mAtoms.size(); ++i)
	{
		AtomRef ar = ref(mAtoms[i]);
		assert(ar < mRef2AtomIndex.size());
		mRef2AtomIndex[ar] = i;
	}

	// collect the restraints

	for (auto &r : mBondRestraints)
		mRestraints.push_back(&r);
	for (auto &r : mAngleRestraints)
		mRestraints.push_back(&r);
	for (auto &r : mTransPeptideRestraints)
		mRestraints.push_back(&r);
	for (auto &r : mTorsionRestraints)
		mRestraints.push_back(&r);
	for (auto &r : mPlanarityRestraints)
		mRestraints.push_back(&r);
	for (auto &r : mChiralVolumeRestraints)
		mRestraints.push_back(&r);
	for (auto &r : mNonBondedContactRestraints)
		mRestraints.push_back(&r);
	if (mDensityRestraint)
		mRestraints.push_back(mDensityRestraint.get());

	// report

	if (cif::VERBOSE > 1)
		std::cout << "created " << mBondRestraints.size() << " bond restraints" << std::endl
				  << "created " << mAngleRestraints.size() << " angle restraints" << std::endl
				  << "created " << mTorsionRestraints.size() << " torsion restraints" << std::endl
				  << "created " << mPlanarityRestraints.size() << " plane restraints" << std::endl
				  << "created " << mTransPeptideRestraints.size() << " trans peptide restraints" << std::endl
				  << "created " << mChiralVolumeRestraints.size() << " chiral vol restraints" << std::endl
				  << "created " << mNonBondedContactRestraints.size() << " non-bonded-contact restraints" << std::endl
				  << std::endl;

	AtomLocationProvider loc(mReferencedAtoms);

	if (cif::VERBOSE > 2)
		for (auto r : mRestraints)
			r->print(loc);
}

AtomRef Minimizer::ref(const mmcif::Atom &atom)
{
	std::string atomID = atom.id();
	AtomRef result;

	auto k = mRefIndex.find(atomID);
	if (k != mRefIndex.end())
		result = k->second;
	else
	{
		result = mReferencedAtoms.size();
		mReferencedAtoms.push_back(atom);
		mRefIndex[atomID] = result;
	}

	return result;
}

void Minimizer::addLinkRestraints(const Residue &a, const Residue &b, const Link &link)
{
	auto &c1 = a.compound();
	auto &c2 = b.compound();

	auto getCompoundAtom = [&](const LinkAtom &la)
	{
		return la.compID == 1 ? c1.getAtomByID(la.atomID) : c2.getAtomByID(la.atomID);
	};

	auto getAtom = [&](const LinkAtom &la)
	{
		const Residue &r = la.compID == 1 ? a : b;
		return r.atomByID(la.atomID);
	};

	for (auto &b : link.bonds())
	{
		try
		{
			if (getCompoundAtom(b.atom[0]).typeSymbol == H or
				getCompoundAtom(b.atom[1]).typeSymbol == H)
			{
				continue;
			}

			Atom a1 = getAtom(b.atom[0]);
			Atom a2 = getAtom(b.atom[1]);

			mBondRestraints.emplace_back(ref(a1), ref(a2), b.distance, b.esd);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE)
				std::cerr << "While processing bond restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	for (auto &a : link.angles())
	{
		try
		{
			if (getCompoundAtom(a.atom[0]).typeSymbol == H or
				getCompoundAtom(a.atom[1]).typeSymbol == H or
				getCompoundAtom(a.atom[2]).typeSymbol == H)
			{
				continue;
			}

			Atom a1 = getAtom(a.atom[0]);
			Atom a2 = getAtom(a.atom[1]);
			Atom a3 = getAtom(a.atom[2]);

			mAngleRestraints.emplace_back(ref(a1), ref(a2), ref(a3), a.angle, a.esd);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 1)
				std::cerr << "While processing angle restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	//	for (auto& a: link.torsions())
	//	{
	//		if (a.esd == 0)
	//			continue;
	//
	//		try
	//		{
	//			if (getCompoundAtom(a.atom[0]).typeSymbol == H or
	//				getCompoundAtom(a.atom[1]).typeSymbol == H or
	//				getCompoundAtom(a.atom[2]).typeSymbol == H or
	//				getCompoundAtom(a.atom[3]).typeSymbol == H)
	//			{
	//				continue;
	//			}
	//
	//			Atom a1 = getAtom(a.atom[0]);
	//			Atom a2 = getAtom(a.atom[1]);
	//			Atom a3 = getAtom(a.atom[2]);
	//			Atom a4 = getAtom(a.atom[3]);
	//
	//			mTorsionRestraints.emplace_back(ref(a1), ref(a2), ref(a3), ref(a4), a.angle, a.esd, a.period);
	//		}
	//		catch (const exception& ex)
	//		{
	//			if (cif::VERBOSE)
	//				std::cerr << "While processing torsion restraints: " << ex.what() << std::endl;
	//			continue;
	//		}
	//	}

	for (auto &cv : link.chiralCentres())
	{
		try
		{
			if (getCompoundAtom(cv.atom[0]).typeSymbol == H or
				getCompoundAtom(cv.atom[1]).typeSymbol == H or
				getCompoundAtom(cv.atom[2]).typeSymbol == H)
			{
				continue;
			}

			Atom cc = getAtom(cv.atomCentre);
			Atom a1 = getAtom(cv.atom[0]);
			Atom a2 = getAtom(cv.atom[1]);
			Atom a3 = getAtom(cv.atom[2]);

			auto volume = link.chiralVolume(cv.id, a.compoundID(), b.compoundID());

			mChiralVolumeRestraints.emplace_back(ref(cc), ref(a1), ref(a2), ref(a3), volume);
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE)
				std::cerr << "While processing chiral volume restraints: " << ex.what() << std::endl;
			continue;
		}
	}

	for (auto &p : link.planes())
	{
		try
		{
			std::vector<AtomRef> atoms;

			for (auto a : p.atoms)
			{
				if (getCompoundAtom(a).typeSymbol == H)
					continue;

				atoms.push_back(ref(getAtom(a)));
			}

			if (atoms.size() > 3)
				mPlanarityRestraints.emplace_back(PlanarityRestraint{move(atoms), p.esd});
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE)
				std::cerr << "While processing planarity restraints: " << ex.what() << std::endl;
			continue;
		}
	}
}

void Minimizer::printStats()
{
	AtomLocationProvider loc(mReferencedAtoms);

	// for (auto &r : mBondRestraints)
	// 	std::cout << mReferencedAtoms[r.mA] << " -> " << mReferencedAtoms[r.mB] << " = " << r.f(loc) << std::endl;

	double bondScore = rmsz(loc, mBondRestraints);
	double angleScore = rmsz(loc, mAngleRestraints);
	double torsionScore = rmsz(loc, mTorsionRestraints);
	double chiralityVolumeScore = rmsz(loc, mChiralVolumeRestraints);
	double planarityScore = rmsz(loc, mPlanarityRestraints);
	double transpeptideScore = rmsz(loc, mTransPeptideRestraints);
	double nbcScore = rmsz(loc, mNonBondedContactRestraints);
	double densityScore = mDensityRestraint ? mDensityRestraint->f(loc) : 0;

	std::cerr << "  Bonds:              " << bondScore << std::endl
			  << "  Angles:             " << angleScore << std::endl
			  << "  Torsion:            " << torsionScore << std::endl
			  << "  Chirality:          " << chiralityVolumeScore << std::endl
			  << "  Planarity:          " << planarityScore << std::endl
			  << "  Transpeptide:       " << transpeptideScore << std::endl
			  << "  Non-Bonded-Contact: " << nbcScore << std::endl
			  << "  Density:            " << densityScore << std::endl;
}

double Minimizer::score()
{
	AtomLocationProvider loc(mReferencedAtoms);
	return score(loc);
}

double Minimizer::score(const AtomLocationProvider &loc)
{
	double result = 0;
	for (auto r : mRestraints)
		result += r->f(loc);

	if (cif::VERBOSE > 3)
		std::cout << "score: " << result << std::endl;

	return result;
}

// --------------------------------------------------------------------

#include <gsl/gsl_blas.h> // for debugging norm of gradient
#include <gsl/gsl_multimin.h>

class GSLAtomLocation : public AtomLocationProvider
{
  public:
	GSLAtomLocation(std::vector<Atom> &atoms, const std::vector<Point> &fixedAtoms,
		const std::vector<size_t> &index, const gsl_vector *v)
		: AtomLocationProvider(atoms)
		, mFixedLocations(fixedAtoms)
		, mIndex(index)
		, mV(v)
	{
		assert(mIndex.size() == mFixedLocations.size());

		if (cif::VERBOSE > 2)
		{
			for (size_t i = 0; i < mIndex.size(); ++i)
			{
				size_t ri = mIndex[i];
				if (ri == kRefSentinel)
					continue;

				DPoint p = {
					gsl_vector_get(mV, ri * 3),
					gsl_vector_get(mV, ri * 3 + 1),
					gsl_vector_get(mV, ri * 3 + 2)};

				std::cout << mAtoms[i] << p << std::endl;
			}
		}
	}

	virtual DPoint operator[](AtomRef atom) const;

	void storeLocations();

  private:
	const std::vector<Point> &mFixedLocations;
	const std::vector<size_t> &mIndex;
	const gsl_vector *mV;
};

DPoint GSLAtomLocation::operator[](AtomRef atomID) const
{
	assert(atomID < mIndex.size());

	size_t ix = mIndex.at(atomID);

	if (ix == kRefSentinel)
		return mFixedLocations.at(atomID);

	return DPoint(
		gsl_vector_get(mV, ix * 3 + 0),
		gsl_vector_get(mV, ix * 3 + 1),
		gsl_vector_get(mV, ix * 3 + 2));
}

void GSLAtomLocation::storeLocations()
{
	for (size_t i = 0; i < mIndex.size(); ++i)
	{
		size_t ri = mIndex[i];
		if (ri == kRefSentinel)
			continue;

		DPoint p = {
			gsl_vector_get(mV, ri * 3),
			gsl_vector_get(mV, ri * 3 + 1),
			gsl_vector_get(mV, ri * 3 + 2)};

		mAtoms[i].location(p);
	}
}

// --------------------------------------------------------------------

class GSLDFCollector : public DFCollector
{
  public:
	GSLDFCollector(const std::vector<mmcif::Atom> &atoms, const std::vector<size_t> &index, gsl_vector *df)
		: mAtoms(atoms)
		, mIndex(index)
		, mDF(df)
	{
		for (size_t ix : mIndex)
		{
			if (ix == kRefSentinel)
				continue;

			gsl_vector_set(mDF, ix * 3 + 0, 0.0);
			gsl_vector_set(mDF, ix * 3 + 1, 0.0);
			gsl_vector_set(mDF, ix * 3 + 2, 0.0);
		}
	}

	~GSLDFCollector();

	virtual void add(AtomRef atom, double dx, double dy, double dz);

  private:
	// for debugging
	std::string label(AtomRef atom) const
	{
		std::string atomName = " " + mAtoms[atom].labelAtomID();
		atomName += std::string(5 - atomName.length(), ' ');
		return std::to_string(mAtoms[atom].labelSeqID()) + atomName;
	}

	const std::vector<mmcif::Atom> &mAtoms;
	const std::vector<size_t> &mIndex;
	gsl_vector *mDF;
};

GSLDFCollector::~GSLDFCollector()
{
	if (cif::VERBOSE > 1)
	{
		std::cerr << std::string(cif::get_terminal_width(), '-') << std::endl
				  << "Collected gradient: " << std::endl;

		for (size_t i = 0; i < mAtoms.size(); ++i)
		{
			size_t ix = mIndex[i];
			if (ix == kRefSentinel)
				continue;

			double dx = gsl_vector_get(mDF, ix * 3 + 0);
			double dy = gsl_vector_get(mDF, ix * 3 + 1);
			double dz = gsl_vector_get(mDF, ix * 3 + 2);

			std::cerr << "atom: " << label(i) << " d: " << std::setprecision(10) << dx << " " << dy << " " << dz << std::endl;
		}

		std::cerr << std::string(cif::get_terminal_width(), '-') << std::endl;
	}
}

void GSLDFCollector::add(AtomRef atom, double dx, double dy, double dz)
{
	assert(atom < mIndex.size());

	size_t ix = mIndex[atom];
	if (ix != kRefSentinel)
	{
		gsl_vector_set(mDF, ix * 3 + 0, gsl_vector_get(mDF, ix * 3 + 0) + dx);
		gsl_vector_set(mDF, ix * 3 + 1, gsl_vector_get(mDF, ix * 3 + 1) + dy);
		gsl_vector_set(mDF, ix * 3 + 2, gsl_vector_get(mDF, ix * 3 + 2) + dz);

		if (cif::VERBOSE > 4)
			std::cerr << "atom: " << label(atom) << " d: " << std::setprecision(10) << dx << ", " << dy << ", " << dz << std::endl;
	}
}

// --------------------------------------------------------------------

class GSLMinimizer : public Minimizer
{
  public:
	GSLMinimizer(const mmcif::Structure &structure, const mmcif::BondMap &bm, float plane5AtomsESD)
		: Minimizer(structure, bm, plane5AtomsESD)
	{
	}

	virtual void Finish()
	{
		Minimizer::Finish();

		for (auto &a : mReferencedAtoms)
			mFixedLocations.push_back(a.location());
	}

	~GSLMinimizer()
	{
		if (m_s != nullptr)
			gsl_multimin_fdfminimizer_free(m_s);
	}

	virtual double refine(bool storeAtoms);
	virtual std::vector<std::pair<std::string, mmcif::Point>> getAtoms() const;
	virtual void storeAtomLocations();

  private:
	static double F(const gsl_vector *v, void *params);
	static void Df(const gsl_vector *v, void *params, gsl_vector *df);
	static void Fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

	double F(const gsl_vector *v);
	void Df(const gsl_vector *v, gsl_vector *df);
	void Fdf(const gsl_vector *x, double *f, gsl_vector *df);

	std::vector<Point> mFixedLocations;
	gsl_multimin_fdfminimizer *m_s = nullptr;
};

double GSLMinimizer::refine(bool storeAtoms)
{
	const size_t iterations = 4000;

	gsl_multimin_function_fdf fdf = {};
	fdf.f = &GSLMinimizer::F;
	fdf.df = &GSLMinimizer::Df;
	fdf.fdf = &GSLMinimizer::Fdf;
	fdf.n = mAtoms.size() * 3;
	fdf.params = this;

	auto T = gsl_multimin_fdfminimizer_conjugate_pr;
	auto x = gsl_vector_alloc(3 * mAtoms.size());

	size_t ix = 0;
	for (auto &a : mAtoms)
	{
		auto l = a.location();
		gsl_vector_set(x, ix++, l.mX);
		gsl_vector_set(x, ix++, l.mY);
		gsl_vector_set(x, ix++, l.mZ);
	}

	m_s = gsl_multimin_fdfminimizer_alloc(T, 3 * mAtoms.size());

	float tolerance = 0.06f;
	double stepSize = 0.1 * gsl_blas_dnrm2(x);

	gsl_multimin_fdfminimizer_set(m_s, &fdf, x, stepSize, tolerance);

	double gradLim = std::sqrt(mRestraints.size()) * 0.15;
	if (gradLim < 0.3)
		gradLim = 0.3;

	for (size_t i = 0; i < iterations; ++i)
	{
		int status = gsl_multimin_fdfminimizer_iterate(m_s);

		if (cif::VERBOSE > 1)
		{
			size_t ix = 0;
			for (auto &a : mAtoms)
			{
				auto l = a.location();

				Point p{
					static_cast<float>(gsl_vector_get(m_s->x, ix + 0)),
					static_cast<float>(gsl_vector_get(m_s->x, ix + 1)),
					static_cast<float>(gsl_vector_get(m_s->x, ix + 2))};

				ix += 3;

				std::cerr << a << " l: " << l << " => p: " << p << " d = " << (p - l) << std::endl;
			}
		}

		if (status != 0)
		{
			if (status != GSL_ENOPROG)
				std::cerr << "Unexpected result from gsl_multimin_fdfminimizer_iterate: " << status << std::endl;
			else if (cif::VERBOSE)
				std::cerr << "Minimizer stopped at iteration " << i << " at " << m_s->f << std::endl;
			break;
		}

		status = gsl_multimin_test_gradient(m_s->gradient, gradLim);

		if (cif::VERBOSE > 1)
		{
			double norm = gsl_blas_dnrm2(m_s->gradient);
			std::cout << "iteration number " << i << " with f: " << m_s->f
					  << " status from gsl_multimin_test_gradient() " << status << " for norm "
					  << norm << std::endl;
		}

		if (status == GSL_SUCCESS)
		{
			if (cif::VERBOSE)
				std::cerr << "Minimum found at iteration " << i << " at " << m_s->f << std::endl;
			break;
		}

		if (status != GSL_CONTINUE)
			break;
	}

	gsl_vector_free(x);

	if (storeAtoms)
		storeAtomLocations();

	return m_s->f;
}

std::vector<std::pair<std::string, mmcif::Point>> GSLMinimizer::getAtoms() const
{
	std::vector<std::pair<std::string, mmcif::Point>> result;

	for (size_t i = 0; i < mRef2AtomIndex.size(); ++i)
	{
		size_t ri = mRef2AtomIndex[i];
		if (ri == kRefSentinel)
			continue;

		DPoint p = {
			gsl_vector_get(m_s->x, ri * 3),
			gsl_vector_get(m_s->x, ri * 3 + 1),
			gsl_vector_get(m_s->x, ri * 3 + 2)};

		result.emplace_back(mReferencedAtoms[i].id(), p);
	}

	return result;
}

void GSLMinimizer::storeAtomLocations()
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, m_s->x);
	loc.storeLocations();
}

double GSLMinimizer::F(const gsl_vector *v, void *params)
{
	GSLMinimizer *self = reinterpret_cast<GSLMinimizer *>(params);
	return self->F(v);
}

void GSLMinimizer::Df(const gsl_vector *v, void *params, gsl_vector *df)
{
	GSLMinimizer *self = reinterpret_cast<GSLMinimizer *>(params);
	self->Df(v, df);
}

void GSLMinimizer::Fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
	GSLMinimizer *self = reinterpret_cast<GSLMinimizer *>(params);
	self->Fdf(v, f, df);

	if (cif::VERBOSE > 1)
		std::cout << "FDF => " << std::setprecision(10) << *f << std::endl;
}

double GSLMinimizer::F(const gsl_vector *v)
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, v);

	//	return score(loc);
	auto F = score(loc);
	if (cif::VERBOSE > 1)
		std::cout << "F => " << std::setprecision(10) << F << std::endl;
	return F;
}

void GSLMinimizer::Df(const gsl_vector *v, gsl_vector *df)
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, v);

	GSLDFCollector c(mReferencedAtoms, mRef2AtomIndex, df);

	for (auto r : mRestraints)
		r->df(loc, c);
}

void GSLMinimizer::Fdf(const gsl_vector *x, double *f, gsl_vector *df)
{
	GSLAtomLocation loc(mReferencedAtoms, mFixedLocations, mRef2AtomIndex, x);
	*f = score(loc);

	GSLDFCollector c(mReferencedAtoms, mRef2AtomIndex, df);

	for (auto r : mRestraints)
		r->df(loc, c);
}

// --------------------------------------------------------------------

Minimizer *Minimizer::create(const Polymer &poly, int first, int last, const BondMap &bonds,
	const XMap &xMap, float mapWeight, float plane5AtomsESD)
{
	std::unique_ptr<Minimizer> result(new GSLMinimizer(*poly.structure(), bonds, plane5AtomsESD));
	result->addPolySection(poly, first, last);
	result->addDensityMap(xMap, mapWeight);
	result->Finish();
	return result.release();
}

Minimizer *Minimizer::create(mmcif::Structure &structure, const std::vector<mmcif::Atom> &atoms,
	const mmcif::BondMap &bm, float plane5AtomsESD, const XMap *xMap, float mapWeight)
{
	std::unique_ptr<Minimizer> result(new GSLMinimizer(structure, bm, plane5AtomsESD));

	std::vector<const mmcif::Residue *> residues;

	for (auto atom : atoms)
	{
		auto &res = structure.getResidue(atom);

		auto ri = std::find_if(residues.begin(), residues.end(), [rp = &res](const Residue *r)
			{ return r == rp; });
		if (ri != residues.end())
			continue;

		residues.emplace_back(&res);
	}

	// sort by asym, seq_id

	sort(residues.begin(), residues.end(), [](const mmcif::Residue *a, const mmcif::Residue *b)
		{
		int d = a->asymID().compare(b->asymID());
		if (d == 0)
			d = a->seqID() - b->seqID();
		return d < 0; });

	auto &polymers = structure.polymers();

	for (auto ri = residues.begin(); ri != residues.end(); ++ri)
	{
		auto res = *ri;
		auto monomer = dynamic_cast<const mmcif::Monomer *>(res);

		if (monomer == nullptr)
		{
			result->addResidue(*res);
			continue;
		}

		int startSeqID = monomer->seqID();
		int endSeqID = startSeqID;

		while (ri != residues.end())
		{
			if ((*ri)->seqID() != endSeqID + 1 or (*ri)->asymID() != monomer->asymID())
				break;
			++endSeqID;
			++ri;
		}

		auto pi = find_if(polymers.begin(), polymers.end(), [id = monomer->asymID()](mmcif::Polymer &poly)
			{ return poly.asymID() == id; });
		if (pi == polymers.end())
			throw std::runtime_error("Polymer not found for asym ID " + monomer->asymID());

		result->addPolySection(*pi, startSeqID, endSeqID);
	}

	// Add any residue that might be bonded to our list of residues via a struct_conn record

	auto &db = structure.datablock();
	auto &struct_conn = db["struct_conn"];
	std::vector<std::tuple<const Residue *, const Residue *, std::string>> linked;

	for (const auto &[ptnr1_label_asym_id, ptnr1_label_comp_id, ptnr1_label_seq_id, ptnr1_auth_seq_id,
			 ptnr2_label_asym_id, ptnr2_label_comp_id, ptnr2_label_seq_id, ptnr2_auth_seq_id,
			 link_id] : struct_conn.rows<std::string, std::string, int, std::string, std::string, std::string, int, std::string, std::string>("ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id", "ptnr1_auth_seq_id",
			 "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id", "ptnr2_auth_seq_id",
			 "ccp4_link_id"))
	{
		auto ai = find_if(residues.begin(), residues.end(),
			[asym_id = ptnr1_label_asym_id, seq_id = ptnr1_label_seq_id, auth_seq_id = ptnr1_auth_seq_id](const Residue *res)
			{ return res->asymID() == asym_id and res->seqID() == seq_id and res->authSeqID() == auth_seq_id; });
		auto bi = find_if(residues.begin(), residues.end(),
			[asym_id = ptnr2_label_asym_id, seq_id = ptnr2_label_seq_id, auth_seq_id = ptnr2_auth_seq_id](const Residue *res)
			{ return res->asymID() == asym_id and res->seqID() == seq_id and res->authSeqID() == auth_seq_id; });

		if (ai == residues.end() and bi == residues.end())
			continue;

		const Residue *ra = *ai;
		const Residue *rb = *bi;

		if (ai != residues.end() and bi != residues.end())
		{
			linked.emplace_back(ra, rb, link_id);
			continue;
		}

		if (ai != residues.end())
		{
			residues.emplace_back(&structure.getResidue(ptnr2_label_asym_id, ptnr2_label_comp_id, ptnr2_label_seq_id, ptnr2_auth_seq_id));
			linked.emplace_back(ra, residues.back(), link_id);
		}
		else
		{
			residues.emplace_back(&structure.getResidue(ptnr1_label_asym_id, ptnr1_label_comp_id, ptnr1_label_seq_id, ptnr1_auth_seq_id));
			linked.emplace_back(rb, residues.back(), link_id);
		}
	}

	// The struct conn records
	for (const auto &[a, b, link_id] : linked)
	{
		if (not link_id.empty())
		{
			result->addLinkRestraints(*b, *a, link_id);
			continue;
		}

		try
		{
			result->addLinkRestraints(*a, *b, a->compoundID() + "-" + b->compoundID());
		}
		catch (const std::exception &e)
		{
			result->addLinkRestraints(*b, *a, b->compoundID() + "-" + a->compoundID());
		}
	}

	if (xMap != nullptr and mapWeight != 0)
		result->addDensityMap(*xMap, mapWeight);

	result->Finish();

	return result.release();
}

} // namespace pdb_redo