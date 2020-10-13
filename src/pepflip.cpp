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
   Date: maandag 19 februari, 2018
*/

// test 3fvl

#include "pdb-redo.hpp"

#include <fcntl.h>

#include <fstream>
#include <iomanip>
#include <numeric>
#include <future>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "cif++/Secondary.hpp"
#include "cif++/CifUtils.hpp"

#include "Statistics.hpp"
#include "minimizer.h"
#include "ramachandran.h"
#include "skiplist.h"
#include "utils.hpp"


namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace ba = boost::algorithm;

using mmcif::Atom;
using mmcif::Point;
using mmcif::Structure;
using mmcif::Monomer;
using mmcif::BondMap;
using mmcif::Polymer;
using mmcif::StatsCollector;

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Coord_map;
using clipper::Coord_frac;

typedef mmcif::MapMaker<float> MapMaker;

// --------------------------------------------------------------------

size_t
	NTHREADS = std::thread::hardware_concurrency();

float
	gMinAngle		= 90,
	gMinRSCC		= 0.6,
	gMaxRSCCDrop	= 0.95;

const double
	kMinRSCC = 0.6;

// --------------------------------------------------------------------

struct PFRes
{
	PFRes(int seqID, std::string compoundID, std::string entityID, std::string asymID, std::string authID, bool tryFlip)
		: mSeqID(seqID), mCompoundID(compoundID), mEntityID(entityID), mAsymID(asymID), mAuthID(authID), mTryFlip(tryFlip) {}
	
	int			mSeqID;
	std::string		mCompoundID, mEntityID, mAsymID, mAuthID;
	float		mTryAngle = 180;
	bool		mTrust = false;
	bool		mTryFlip;
	bool		mFlipped = false;
};

std::ostream& operator<<(std::ostream& os, const PFRes& r)
{
	os << r.mCompoundID << ' ' << r.mAsymID << r.mSeqID << " (" << r.mAuthID << ')';
	return os;
}

// --------------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const Atom& a);

std::ostream& operator<<(std::ostream& os, const RamachandranScore& s)
{
	switch (s)
	{
		case rsNotAllowed:		os << "N"; break;
		case rsAllowed:			os << "A"; break;
		case rsFavoured:		os << "F"; break;
	}

	return os;
}

// --------------------------------------------------------------------

class AtomLocationSaver
{
  public:

	AtomLocationSaver(mmcif::Residue& r)
		: mCommitted(false)
	{
		for (auto& a: r.atoms())
			mAtoms.emplace_back(std::make_pair(a, a.location()));
	}

	AtomLocationSaver(mmcif::Structure& s, const std::string& asymID, int seqFirst, int seqLast)
		: mCommitted(false)
	{
		for (auto& a: s.atoms())
		{
			if (a.labelAsymID() != asymID or a.labelSeqID() < seqFirst or a.labelSeqID() > seqLast)
				continue;
			
			mAtoms.emplace_back(std::make_pair(a, a.location()));
		}
	}
	
	~AtomLocationSaver()
	{
		if (not mCommitted)
		{
			for (auto& a: mAtoms)
				a.first.location(a.second);
		}
	}
	
	void commit() { mCommitted = true; }
	
	void storeAtomPositions(std::vector<std::pair<std::string,Point>>& atoms) const
	{
		atoms.reserve(mAtoms.size());
		for (auto& a: mAtoms)
			atoms.push_back(make_pair(a.first.id(), a.first.location()));
	}
	
  private:
	std::vector<std::pair<mmcif::Atom,mmcif::Point>> mAtoms;
	bool mCommitted;
};

// --------------------------------------------------------------------

void TrustBasedOnSS(Structure& structure, std::vector<PFRes>& residues)
{
	mmcif::DSSP dssp(structure, 3, false);
	
	const std::set<mmcif::SecondaryStructureType> kTrustedDSSPTypes = {
			mmcif::ssAlphahelix, mmcif::ssHelix_3, mmcif::ssHelix_5, mmcif::ssStrand };
	const int kBufferSS = 1;
	
	enum { LOOP, SS_RUN } state = LOOP;
	mmcif::SecondaryStructureType cur = mmcif::ssLoop;
	
	size_t i = 0, b = 0;
	while (i < residues.size())
	{
		auto ss = dssp(residues[i].mAsymID, residues[i].mSeqID);
		++i;
		
		switch (state)
		{
			case LOOP:
				if (kTrustedDSSPTypes.count(ss))
				{
					state = SS_RUN;
					cur = ss;
					b = i - 1;
				}
				break;
			
			case SS_RUN:
				if (ss != cur)
				{
					// trust range, if long enough
					if (i - b > 2 * kBufferSS + 2)
					{
						for (auto j = b + kBufferSS; j < i - kBufferSS - 2; ++j)
						{
							if (cif::VERBOSE)
								std::cerr << residues[j] << "  - trusted based on secondary structure (" << char(cur) << ')' << std::endl;
							if (not dssp.isAlphaHelixEndBeforeStart(residues[j].mAsymID, residues[j].mSeqID))
								residues[j].mTrust = true;
						}

						if (cif::VERBOSE)
							std::cerr << residues[i - kBufferSS - 2] << "  - n-term of secondary structure (" << char(cur) << "). It will be considered for flips" << std::endl;
					}
					
					state = LOOP;
					--i;	// re-start
				}
				break;
		}
	}
}

class Rotation
{
  public:
	Rotation(Point ca1, Point ca2, float angle)
	{
		mT = ca1;
		
		angle = angle * mmcif::kPI / 180;
		
		auto s = sin(angle / 2);
		auto v = ca2 - ca1;
		v.normalize();
		
		mQ = mmcif::quaternion(cos(angle / 2), v.mX * s, v.mY * s, v.mZ * s);
	}

	Point operator()(const Point& p) const
	{
		Point result = p - mT;
		result.rotate(mQ);
		return result + mT;
	}

  private:
	Point mT;
	mmcif::quaternion mQ;
};

std::tuple<double,double,float,float> ScorePotentialFlips(MapMaker& mm, const Monomer& r1, const Monomer& r2, std::initializer_list<float> angles)
{
	double bestScoreB = 0, bestScoreD = 0;
	float bestAngleB = 0, bestAngleD = 0;
	
	auto& fb = static_cast<clipper::Xmap<float>&>(mm.fb());
	auto& fd = static_cast<clipper::Xmap<float>&>(mm.fd());
	auto& cell = mm.cell();
	
	for (auto angle: angles)
	{
		auto ca1 = r1.atomByID("CA");
		auto ca2 = r2.atomByID("CA");
		
		Rotation rt(ca1.location(), ca2.location(), angle);

		auto c = r1.atomByID("C");
		auto o = r1.atomByID("O");
		auto n = r2.atomByID("N");
//		auto h = r2.atomByID("H");

		double scoreB = 0, scoreD = 0;
		
		for (auto& a: { c, o, n })
		{
			auto l = rt(a.location());

			clipper::Coord_orth p = l;
			clipper::Coord_frac pf = p.coord_frac(cell);
			
			scoreB += static_cast<int>(a.type()) * a.occupancy() * fb.interp<clipper::Interp_cubic>(pf);
			scoreD += static_cast<int>(a.type()) * a.occupancy() * fd.interp<clipper::Interp_cubic>(pf);
		}

		if (bestScoreB < scoreB)
		{
			bestScoreB = scoreB;
			bestAngleB = angle;
		}

		if (bestScoreD < scoreD)
		{
			bestScoreD = scoreD;
			bestAngleD = angle;
		}
	}
	
	return std::make_tuple(bestScoreB, bestScoreD, bestAngleB, bestAngleD);
}

// --------------------------------------------------------------------

void CheckDensity(Structure& structure, std::vector<PFRes>& residues, MapMaker& mm)
{
	const double kMinRatioOrigFlip = 0.9;

	size_t potentialFlipCount = 0;
	
	auto& polymers = structure.polymers();
	
	for (size_t i = 0; i + 1 < residues.size(); ++i)
	{
		auto& cur = residues[i];
		auto& next = residues[i + 1];

		if (next.mAsymID != cur.mAsymID)
		{
			cur.mTryFlip = false;
			continue;
		}

		std::string asymID = cur.mAsymID;
		auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
		assert(poly != polymers.end());
		
		auto& r1 = poly->getBySeqID(cur.mSeqID);
		auto& r2 = poly->getBySeqID(next.mSeqID);

		if (cur.mTrust or not cur.mTryFlip)
		{
			std::cout << cur << "  - residue excluded from peptide flips" << std::endl;
			cur.mTryFlip = false;
			continue;
		}

		cur.mTryFlip = false;
		
		if (not Monomer::areBonded(r1, r2))
		{
			std::cout << cur << "  - Distance between this CA and the next CA doesn't correspond with a peptide bond" << std::endl;
			continue;
		}
		
		try
		{
			double scoreBBefore, scoreBAfter, scoreDBefore, scoreDAfter;
			float angleBAfter, angleDAfter;

			std::tie(scoreBBefore, scoreDBefore, std::ignore, std::ignore) = ScorePotentialFlips(mm, r1, r2, { 330, 0, 30 });
			std::tie(scoreBAfter, scoreDAfter, angleBAfter, angleDAfter) = ScorePotentialFlips(mm, r1, r2, { 150, 180, 210 });
			
			if (scoreBAfter > kMinRatioOrigFlip * scoreBBefore)
			{
				std::cout << cur << "  - might need to be flipped based on the density fit, before: " << scoreBBefore << " after: " << scoreBAfter << std::endl;
				cur.mTryFlip = true;
				cur.mTryAngle = angleBAfter;
				++potentialFlipCount;
			}
			else if (scoreDAfter > scoreDBefore)
			{
				std::cout << cur << "  - might need to be flipped based on the difference density fit, before: " << scoreDBefore << " after: " << scoreDAfter << std::endl;
				cur.mTryFlip = true;
				cur.mTryAngle = angleDAfter;
				++potentialFlipCount;
			}
			else
				std::cout << cur << "  - no reason to flip this peptide based on the density" << std::endl;
		}
		catch (const std::exception& ex)
		{
			std::cerr << "Error processing " << cur << std::endl
				 << ex.what() << std::endl;	
		}
	}
}

// --------------------------------------------------------------------

//const regex kNotAHBondEnergyTypeRX("N(?:R(?:5|55|56|6|66))?");	// <-- N, NR5, NR55, NR56, NR6 or NR66
//
//int CalculateHBondsForAtoms(Structure& s, initializer_list<Atom> atoms)
//{
//	int result = 0;
//	
//	for (auto& a: atoms)
//	{
//		Polymer p(s, a.labelAsymId());
//		auto m = p.getBySeqID(a.labelSeqId());
//		
//		// all atoms near 3.5
//		for (auto& a2: dm.near(a, 3.5f))
//		{
//			// but only O or N
//			if (a2.type() != mmcif::O and a2.type() != mmcif::N)
//				continue;
//			
//			// And if N, not with energytype N, NR5, etc
//			if (a2.type() == mmcif::N)
//			{
//				auto& c = a2.comp();
//				const auto& ca = c.getAtomById(a2.labelAtomId());
//				if (regex_match(ca.typeEnergy, kNotAHBondEnergyTypeRX))
//					continue;
//			}
//			
//			// Not on same residue of course
//			if (a.labelSeqId() == a2.labelSeqId() and a.labelAsymId() == a2.labelAsymId())
//				continue;
//			
//			// if on same chain and both backbone, at least 3 residues apart
//			if (a.labelAsymId() == a2.labelAsymId() and a2.isBackBone())
//			{
//				auto m2 = p.getBySeqID(a2.labelSeqId());
//				if (p.Distance(m, m2) < 3)
//					continue;
//			}
//
//			if (cif::VERBOSE)
//				std::cerr << "HBond between " << a << " and " << a2 << std::endl;
//			
//			result += 1;
//		}
//	}
//	
//	return result;
//}

// --------------------------------------------------------------------

struct PepFlipScore
{
	size_t			fpix;
	std::string			id;
	float			angle;
	mmcif::Point	c;
	bool			refined;
	struct {
		double score;					// from minimizer
		double rsccsRes;
		float densityDifferenceO;
		RamachandranScore z1, z2;
//		int hbonds;
	} orig, flip;
	std::vector<std::pair<std::string,Point>> atoms;

	bool improved() const
	{
		return
				refined
			and
				abs(angle) > gMinAngle
			and
				flip.score < orig.score
			and
				flip.z1 >= rsAllowed and flip.z2 >= rsAllowed
			and
				(flip.z1 + flip.z2) >= (orig.z1 + orig.z2)
			and
				flip.rsccsRes > gMinRSCC
			and
				flip.rsccsRes >= gMaxRSCCDrop * orig.rsccsRes
			and	
				flip.densityDifferenceO > 0
			and
				flip.densityDifferenceO >= orig.densityDifferenceO
//			and
//				flip.hbonds >= orig.hbonds
			;
	}
	
	bool operator<(const PepFlipScore& rhs) const { return fpix < rhs.fpix; }
};

template<typename T>
inline auto best(const T& v, bool best, int prec)
{
	std::stringstream s;
	s << std::fixed << std::setprecision(prec) << std::setw(8) << v;
	std::string vs = s.str();
	return best ?
		coloured(vs.c_str(), cif::scCYAN, cif::scBLUE) :
		coloured(vs.c_str(), cif::scNONE, cif::scNONE, false);
}

template<>
inline auto best(const RamachandranScore& v, bool best, int prec)
{
	std::stringstream s;
	s << v;
	std::string vs = s.str();
	return best ?
		coloured(vs.c_str(), cif::scCYAN, cif::scBLUE) :
		coloured(vs.c_str(), cif::scNONE, cif::scNONE, false);
}

template<>
inline auto best(const int& v, bool best, int prec)
{
	std::string s = std::to_string(v);
	return best ?
		coloured(s.c_str(), cif::scCYAN, cif::scBLUE) :
		coloured(s.c_str(), cif::scNONE, cif::scNONE, false);
}

struct headers1 {};
struct headers2 {};

std::ostream& operator<<(std::ostream& os, const headers1&)
{
	os << "Angle    "
	   << "Score             "
	   << "RSCCS residues    "
	   << "Diff.Density O    "
	   << "RamaZ  "
	   ;
	return os;
}

std::ostream& operator<<(std::ostream& os, const headers2)
{
	os << "         "
	   << "Original Flipped  "
	   << "Original Flipped  "
	   << "Original Flipped  "
	   << "O F O F "
	   ;
	return os;
}

std::ostream& operator<<(std::ostream& os, const PepFlipScore& rhs)
{
	os << best(rhs.angle, abs(rhs.angle) > 90.0f, 1) << ' '

	   << best(rhs.orig.score, false, 1) << ' '
	   << best(rhs.flip.score, rhs.flip.score < rhs.orig.score, 1) << ' '

	   << best(rhs.orig.rsccsRes, false, 3) << ' '
	   << best(rhs.flip.rsccsRes, rhs.flip.rsccsRes > rhs.orig.rsccsRes, 3) << ' '

	   << best(rhs.orig.densityDifferenceO, false, 3) << ' '
	   << best(rhs.flip.densityDifferenceO, rhs.flip.densityDifferenceO > 0 and rhs.flip.densityDifferenceO > rhs.orig.densityDifferenceO, 3) << ' '

	   << best(rhs.orig.z1, false, 0) << ' '
	   << best(rhs.flip.z1, rhs.flip.z1 > rhs.orig.z1, 0) << ' '

	   << best(rhs.orig.z2, false, 0) << ' '
	   << best(rhs.flip.z2, rhs.flip.z2 > rhs.orig.z2, 0) << ' '

//	   << best(rhs.orig.hbonds, false, 0) << ' '
//	   << best(rhs.flip.hbonds, rhs.flip.hbonds > rhs.orig.hbonds, 1)
	   
	   ;

	return os;
}

PepFlipScore FlipPeptide(Structure& structure, const Polymer& poly,
	const std::vector<PFRes>& residues, size_t i, StatsCollector& sc, MapMaker& mm,
	mmcif::BondMap& bm, const std::string& algorithm, float mapWeight, float plane5AtomsESD,
	bool testMode)
{
	auto& r1 = poly.getBySeqID(residues[i].mSeqID);
	auto& r2 = poly.getBySeqID(residues[i + 1].mSeqID);
	
	const clipper::Xmap<float>& fd = mm.fd();
	
	auto getDf = [&fd](const Point& p)
	{
		Coord_orth po = p;
		Coord_frac pf = po.coord_frac(fd.cell());
		return fd.interp<clipper::Interp_cubic>(pf);
	};
	
	AtomLocationSaver s(structure, r1.asymID(), r1.seqID() - 1, r2.seqID() + 1);

	PepFlipScore result = { i, boost::lexical_cast<std::string>(residues[i]) };

	auto prePro = [&](size_t i)
	{
		return (i + 1 < residues.size() and residues[i + 1].mCompoundID == "PRO" and residues[i].mCompoundID != "PRO");
	};

	auto ca1 = r1.atomByID("CA");
	auto ca2 = r2.atomByID("CA");
	
	auto c = r1.atomByID("C");
	auto o = r1.atomByID("O");
	auto n = r2.atomByID("N");
//	auto h = r2.atomByID("H");

	result.c = c.location();

	std::unique_ptr<Minimizer> origMinimizer(Minimizer::create(algorithm, poly, r1.seqID() - 1, r2.seqID() + 1, bm, mm.fb(), mapWeight, plane5AtomsESD));
	result.orig.score = origMinimizer->refine(false);

	result.orig.rsccsRes = sc.collect({ &r1, &r2 }).RSCCS;
	result.orig.densityDifferenceO = getDf(o.location());
	result.orig.z1 = calculateRamachandranScore(r1.compoundID(), prePro(i), r1.phi(), r1.psi());
	result.orig.z2 = calculateRamachandranScore(r2.compoundID(), prePro(i + 1), r2.phi(), r2.psi());
//	result.orig.hbonds = CalculateHBondsForAtoms(structure, { o, n }, dm);

//	auto cl = c.location();
	auto ol = o.location();
//	auto nl = n.location();

	auto angle = residues[i].mTryAngle;
	Rotation rt(ca1.location(), ca2.location(), angle);

	c.location(rt(c.location()));
	o.location(rt(o.location()));
	n.location(rt(n.location()));

	if (testMode)
		structure.getFile().save("unrefined-" + r1.labelID() + ".pdb");

	std::unique_ptr<Minimizer> flipMinimizer(Minimizer::create(algorithm, poly, r1.seqID() - 1, r2.seqID() + 1, bm, mm.fb(), mapWeight, plane5AtomsESD));

	double unrefinedScore = flipMinimizer->score();
	result.flip.score = flipMinimizer->refine(true);
	result.refined = result.flip.score < unrefinedScore;

	if (testMode)
		structure.getFile().save("refined-" + r1.labelID() + ".pdb");

//	auto rcl = c.location();
	auto rol = o.location();
//	auto rnl = n.location();

	result.angle = DihedralAngle(ol, ca1.location(), ca2.location(), rol);

	result.flip.rsccsRes = sc.collect({ &r1, &r2 }).RSCCS;
	result.flip.densityDifferenceO = getDf(o.location());
	result.flip.z1 = calculateRamachandranScore(r1.compoundID(), prePro(i), r1.phi(), r1.psi());
	result.flip.z2 = calculateRamachandranScore(r2.compoundID(), prePro(i + 1), r2.phi(), r2.psi());
//	result.flip.hbonds = CalculateHBondsForAtoms(structure, { o, n });

//std::cout << residues[i] << std::endl;
	if (result.improved())
		result.atoms = flipMinimizer->getAtoms();

	return result;
}

std::vector<size_t> CombineTrustBasedOnNCS(const std::list<Polymer>& polymers, std::vector<PFRes>& residues)
{
	typedef std::vector<PFRes>::iterator PFResIter;
	
	std::set<std::string> entities;
	for (auto& p: polymers)
		entities.insert(p.entityID());
	
	std::vector<size_t> result;
	
	for (const std::string& entityID: entities)
	{
		std::vector<std::pair<PFResIter,PFResIter>> pri;
		
		for (auto ri = residues.begin(); ri != residues.end(); ++ri)
		{
			if (ri->mEntityID == entityID)
			{
				auto bri = ri, eri = ri + 1;
				
				while (eri != residues.end() and eri->mEntityID == entityID and eri->mAsymID == bri->mAsymID)
					++eri;
				
				pri.emplace_back(bri, eri);
				ri = eri - 1;
			}
		}

		int seqID = 0;
		std::vector<PFResIter> tri;

		while (pri.size() > 1)
		{
			tri.clear();
			for (auto& ri: pri)
			{
				if (ri.first != ri.second and ri.first->mSeqID == seqID)
				{
					tri.push_back(ri.first);
					++ri.first;
				}
			}
			
			if (tri.size() < 2)
			{
				++seqID;
				continue;
			}
			
			size_t flipCount = accumulate(tri.begin(), tri.end(), 0, [](int sum, auto ri) { return sum + (ri->mFlipped ? 1 : 0); });
			if (flipCount < tri.size() and flipCount > 0)
			{
				for (auto ri: tri)
				{
					if (ri->mTryFlip)
						continue;
					
					result.push_back(ri - residues.begin());
				}
			}
			
			pri.erase(remove_if(pri.begin(), pri.end(), [](auto i) { return i.first == i.second; }), pri.end());
		}
	}
	
	return result;
}

void JoinFlips(Structure& structure, const Polymer& poly, MapMaker& mm, mmcif::BondMap& bm,
	const std::string& algorithm, float mapWeight, float plane5AtomsESD,
	const std::vector<PFRes>& residues, std::vector<PepFlipScore>& flipped, size_t f1, size_t f2)
{
	size_t ri1 = flipped[f1].fpix;
	size_t ri2 = flipped[f2].fpix;
	
	int s1 = residues[ri1].mSeqID - 1;
	int s2 = residues[ri2].mSeqID + 1;
	
	std::cout << "Joining flips from " << residues[ri1] << " to " << residues[ri2] << std::endl;
	
	// Eerste poging, gewoon alles flippen en dan verfijnen...
	
	for (size_t fi = f1; fi <= f2; ++fi)
	{
		size_t ri = flipped[fi].fpix;
		
		auto& r1 = poly.getBySeqID(residues[ri].mSeqID);
		auto& r2 = poly.getBySeqID(residues[ri + 1].mSeqID);

		auto ca1 = r1.atomByID("CA");
		auto ca2 = r2.atomByID("CA");
		
		auto c = r1.atomByID("C");
		auto o = r1.atomByID("O");
		auto n = r2.atomByID("N");
		
		auto angle = residues[ri].mTryAngle;
		Rotation rt(ca1.location(), ca2.location(), angle);
	
		c.location(rt(c.location()));
		o.location(rt(o.location()));
		n.location(rt(n.location()));
		
		flipped[fi].atoms.clear();
	}
	
	std::unique_ptr<Minimizer> flipMinimizer(Minimizer::create(algorithm, poly,
		s1, s2, bm, mm.fb(), mapWeight, plane5AtomsESD));

	double unrefinedScore = flipMinimizer->score();
	double refinedScore = flipMinimizer->refine(true);

	if (refinedScore > unrefinedScore)
	{
		std::cerr << "Oeps, deze verfijnde niet" << std::endl;
		
		for (size_t fi = f1; fi <= f2; ++fi)
			flipped[fi].refined = false;
	}
}

void FlipPeptides(Structure& structure, const std::string& asymID,
	int resFirst, int resLast, const SkipList& skip,
	MapMaker& mm, bool trustDSSP,
	const std::string& algorithm, float mapWeight, float plane5AtomsESD,
	std::ofstream& cootScript, bool testMode)
{
	StatsCollector sc(mm, structure, false /*electronScattering*/);

	auto& polymers = structure.polymers();
	
	std::vector<PFRes> residues;
	for (auto& p: polymers)
	{
		for (auto& m: p)
		{
			bool tryFlip = asymID.empty() or (m.asymID() == asymID and m.seqID() >= resFirst and m.seqID() <= resLast);
			
			if (tryFlip)
				tryFlip = find_if(skip.begin(), skip.end(), [&](auto s) { return s.asymID == m.asymID() and s.seqID == m.seqID(); }) == skip.end();
			
			residues.emplace_back(m.seqID(), m.compoundID(), p.entityID(), m.asymID(), m.authID(), tryFlip);
		}
	}
	
	if (trustDSSP)
	{
		if (cif::VERBOSE)
			std::cerr << std::endl
				 << "Calculating DSSP" << std::endl;
	
		TrustBasedOnSS(structure, residues);

		if (cif::VERBOSE)
			std::cerr << "DSSP done" << std::endl;
	}

	if (cif::VERBOSE)
		std::cerr << std::endl
			 << "Initial check" << std::endl;

	auto df = std::async(std::launch::async, std::bind(CheckDensity, std::ref(structure), std::ref(residues), std::ref(mm)));

	mmcif::BondMap bm(structure);
	
	df.wait();

	std::cout << std::endl
	     << "Potential number of flips: " << accumulate(residues.begin(), residues.end(), 0, [](int sum, auto& r) { return sum + (r.mTryFlip ? 1 : 0); }) << std::endl
	     << std::endl;

	std::vector<PepFlipScore> flipped;
	flipped.reserve(residues.size());
	
	std::unique_ptr<cif::Progress> p;
	if (not cif::VERBOSE)
		p.reset(new cif::Progress(residues.size(), "Flipping, 1st round"));
	
	if (NTHREADS == 1)
	{
		for (size_t i = 0; i + 1 < residues.size(); ++i)
		{
			if (p)
				p->progress(i);
	
			if (residues[i].mTrust or not residues[i].mTryFlip)
				continue;
	
			std::string asymID = residues[i].mAsymID;
			auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
			assert(poly != polymers.end());
	
			try
			{
				flipped.emplace_back(
					FlipPeptide(structure, *poly, residues, i, sc, mm, bm, algorithm, mapWeight, plane5AtomsESD, testMode));
				residues[i].mFlipped = flipped.back().improved();
			}
			catch (const std::exception& ex)
			{
				std::cerr << "Error processing " << residues[i] << ": " << ex.what() << std::endl;
				continue;
			}
		}
	}
	else
	{
		std::list<std::thread> t;
		std::atomic<size_t> next(0);
		std::mutex m;
		
		for (size_t ti = 0; ti < NTHREADS; ++ti)
			t.emplace_back([&]()
			{
				Structure s(structure);
				StatsCollector sc(mm, s, false /*electronScattering*/);
				auto& polymers = s.polymers();

				for (;;)
				{
					size_t i = next++;
					
					if (i + 1 >= residues.size())
						break;
					
					if (p)
						p->progress(i);
					
					if (residues[i].mTrust or not residues[i].mTryFlip)
						continue;
			
					std::string asymID = residues[i].mAsymID;
					auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
					assert(poly != polymers.end());
			
					try
					{
						auto flipResult = FlipPeptide(s, *poly, residues, i, sc, mm, bm, algorithm, mapWeight, plane5AtomsESD, testMode);
						residues[i].mFlipped = flipResult.improved();
						
						std::unique_lock lock(m);
						flipped.emplace_back(std::move(flipResult));
					}
					catch (const std::exception& ex)
					{
						std::cerr << "Error processing " << residues[i] << ": " << ex.what() << std::endl;
						continue;
					}
				}
			});
	
		for (auto& ti: t)
			ti.join();
		
		sort(flipped.begin(), flipped.end());
		
	}

	if (p)
		p->progress(residues.size());

	std::cout << std::endl;

	// second round, flip any residue whose NCS copy was flipped
	if (asymID.empty())
	{
		auto residuesRound2 = CombineTrustBasedOnNCS(polymers, residues);
		if (not residuesRound2.empty())
		{
			std::cout << std::endl
				 << "Potential number of flips in second round (NCS copies): " << residuesRound2.size() << std::endl
				 << std::endl;
			
			std::unique_ptr<cif::Progress> p;
			if (not cif::VERBOSE)
				p.reset(new cif::Progress(residuesRound2.size(), "Flipping, 2nd round"));
			
			for (size_t i: residuesRound2)
			{
				std::string asymID = residues[i].mAsymID;
				auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
				assert(poly != polymers.end());
		
				try
				{
					flipped.emplace_back(
						FlipPeptide(structure, *poly, residues, i, sc, mm, bm, algorithm, mapWeight, plane5AtomsESD, testMode));
					residues[i].mFlipped = flipped.back().improved();
				}
				catch (const std::exception& ex)
				{
					std::cerr << "Error processing " << residues[i] << ": " << ex.what() << std::endl;
					continue;
				}
	
				if (p)
					p->consumed(1);
			}
		}
	}
	
	// third step, some of the flips might overlap, process them separately
	for (size_t i = 0; i + 1 < flipped.size(); ++i)
	{
		if (not flipped[i].improved())
			continue;
		
		std::string asymID = residues[flipped[i].fpix].mAsymID;
		int firstSeq = residues[flipped[i].fpix].mSeqID - 1;
		int lastSeq = firstSeq + 3;
		
		size_t j = i;
		while (j + 1 < flipped.size() and flipped[j + 1].improved() and lastSeq + 1 >= residues[flipped[j + 1].fpix].mSeqID)
		{
			lastSeq = residues[flipped[j].fpix].mSeqID + 2;
			++j;
		}
		
		if (j > i)
		{
			auto poly = find_if(polymers.begin(), polymers.end(), [&](auto& p) { return p.asymID() == asymID; });
			JoinFlips(structure, *poly, mm, bm, algorithm, mapWeight, plane5AtomsESD, residues, flipped, i, j);
		}
	}

	size_t idLen = 12;
	for (auto& s: flipped)
		if (idLen < s.id.length())
			idLen = s.id.length();
	
	std::cout << std::endl
		 << std::string(idLen + 8, ' ') << headers1() << std::endl
		 << std::string(idLen + 8, ' ') << headers2() << std::endl
		 << std::string(cif::get_terminal_width(), '-') << std::endl;
	
	std::vector<std::string> flippedIDs;
	for (auto& score: flipped)
	{
		std::string id = score.id + std::string(idLen - score.id.length(), ' ');

		auto red = [](const char* s){ return cif::coloured(s, cif::scWHITE, cif::scRED); };
		
		if (not score.improved())
		{
			std::cout << id << " " << "noflip" << " " << score << std::endl;
			continue;
		}

		std::cout << id << " " << red("flip") << "   " << score << std::endl;

		if (cootScript.is_open())
			cootScript << "(list \"Flipped " << score.id << "\" "
					   << score.c.mX << ' ' << score.c.mY << ' ' << score.c.mZ << " )" << std::endl;
		
		for (auto a: score.atoms)
			structure.getAtomByID(a.first).location(a.second);
		
		flippedIDs.push_back(score.id);
	}

	std::cout << std::endl
		 << std::string(cif::get_terminal_width(), '-') << std::endl
		 << "Summary: Flipped " << flippedIDs.size() << " peptides" << std::endl
		 << ba::join(flippedIDs, ", ") << std::endl;
}

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	po::options_description visible_options(fs::path(argv[0]).filename().string() + " options");
	visible_options.add_options()
		("hklin",				po::value<std::string>(),	"reflections file")
		("xyzin",				po::value<std::string>(),	"coordinates file")
		("xyzout,o",			po::value<std::string>(),	"output coordinates to this file")
		
		("log",					po::value<std::string>(),	"Write verbose output to log file")
		
		("asym-id",				po::value<std::string>(),	"Asymetric unit to pepflip")
		("algorithm",			po::value<std::string>(),	"Refinement algorithm to use, can be sa or gsl")
		("res-first",			po::value<int>(),		"Sequence number for first residue to pepflip, default = 1")
		("res-last",			po::value<int>(),		"Sequence number for last residue to pepflip, default is last in sequence")
		("dict",				po::value<std::vector<std::string>>(),
														"Dictionary file containing restraints for residues in this specific target, can be specified multiple times.")

		("no-dssp",										"By default, residues in helix and strand are trusted, use this switch to turn this off")

		("skip",				po::value<std::string>(),	"File containing the skip list: the residues that should not be flipped")

		("std::map-weight",			po::value<float>(),		"Map weight in minimisation, default is 60")

		("minimal-angle",		po::value<float>(),		"Minimal angle for flip in degrees")
		("minimal-rscc",		po::value<float>(),		"Minimal RSCC score for the two residues")
		("max-rscc-drop",		po::value<float>(),		"The flipped rscc should be at least this factor times the orignal rscc")

		("plane-5-atoms-esd",	po::value<float>(),		"ESD for the atoms in the 5 atom peptide bond plane, default is 0.11")
		("coot-script",									"Write a Coot script for the changed peptides")

		("help,h",										"Display help message")
		("version",										"Print version")
		
		("sampling-rate",		po::value<float>(),		"Sampling rate")
		("max-threads",			po::value<uint16_t>(),	"Max number of threads")
		
		("verbose,v",									"Verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),		"Debug level (for even more verbose output)")
		("test",										"Do a test run, flipping all potential residues and writing intermediate files for further debugging")
		;

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("hklin", 1);
	p.add("xyzin", 1);
	p.add("xyzout", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "pepflip.conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / "pepflip.conf";
	
	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}
	
	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cerr << visible_options << std::endl;
		exit(0);
	}
	
	if (vm.count("xyzin") == 0 or vm.count("hklin") == 0)
	{
		std::cerr << "Input files not specified" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	if (vm.count("log"))
	{
		std::string logFile = vm["log"].as<std::string>();
		
	    // open the log file
	    int fd = open(logFile.c_str(), O_CREAT|O_APPEND|O_RDWR, 0644);
	    if (fd < 0)
	        throw std::runtime_error("Opening log file " + logFile + " failed: " + strerror(errno));
	
	    // redirect stdout and stderr to the log file
	    dup2(fd, STDOUT_FILENO);
	    dup2(fd, STDERR_FILENO);
	    close(fd);
	}

	if (vm.count("max-threads"))
		NTHREADS = vm["max-threads"].as<uint16_t>();
	if (NTHREADS > std::thread::hardware_concurrency())
		NTHREADS = std::thread::hardware_concurrency();

	if (cif::VERBOSE or NTHREADS < 1)
		NTHREADS = 1;

	if (vm.count("dict"))
	{
		for (auto dict: vm["dict"].as<std::vector<std::string>>())
			mmcif::CompoundFactory::instance().pushDictionary(dict);
	}

	mmcif::File f(vm["xyzin"].as<std::string>());
	Structure structure(f);
	
	fs::path mtzFile = vm["hklin"].as<std::string>();

	float samplingRate = 1.5;
	if (vm.count("sampling-rate"))
		samplingRate = vm["sampling-rate"].as<float>();

	MapMaker mm;
	mm.loadMTZ(mtzFile, samplingRate);

	std::string asymID;
	if (vm.count("asym-id"))
		asymID = vm["asym-id"].as<std::string>();
		
	int resFirst = 1;
	if (vm.count("res-first"))
		resFirst = vm["res-first"].as<int>();
	
	int resLast = std::numeric_limits<int>::max();
	if (vm.count("res-last"))
		resLast = vm["res-last"].as<int>();

	std::string algorithm = "gsl";
	if (vm.count("algorithm"))
		algorithm = vm["algorithm"].as<std::string>();

	if (vm.count("minimal-angle"))	gMinAngle = vm["minimal-angle"].as<float>();
	if (vm.count("minimal-rscc"))	gMinRSCC = vm["minimal-rscc"].as<float>();
	if (vm.count("max-rscc-drop"))	gMaxRSCCDrop = vm["max-rscc-drop"].as<float>();


	bool trustDSSP = vm.count("no-dssp") == 0;
	
	float mapWeight = 60;
	if (vm.count("map-weight"))
		mapWeight = vm["map-weight"].as<float>();
	
	float plane5AtomsESD = 0.11;
	if (vm.count("plane-5-atoms-esd"))
		plane5AtomsESD = vm["plane-5-atoms-esd"].as<float>();
	
	fs::path xyzout;
	if (vm.count("xyzout"))
		xyzout = vm["xyzout"].as<std::string>();
	else
	{
		xyzout = vm["xyzin"].as<std::string>();
		if (xyzout.extension() == ".gz")
			xyzout = xyzout.stem();
		xyzout = xyzout.parent_path() / (xyzout.filename().stem().string() + "-flipped.cif");
	}
	
	xyzout = fs::canonical(xyzout);
	
	std::ofstream cootScript;
	if (vm.count("coot-script"))
	{
		fs::path cs = xyzout.extension() == ".gz" ?
			xyzout.parent_path() / (xyzout.filename().stem().stem().string() + ".scm") :
			xyzout.parent_path() / (xyzout.filename().stem().string() + ".scm");
		cootScript.open(cs);
		
		if (not cootScript.is_open())	
			throw std::runtime_error("Failed to open coot script " + cs.string());
		
		cootScript << "(interesting-things-gui \"Pepflips\"" << std::endl
				   << "(list" << std::endl;
	}
	
	bool testMode = false;
	if (vm.count("test"))
	{
		testMode = true;
		
		fs::path testOutputDir = fs::current_path() / (xyzout.filename().stem().string() + "-test");
		
		if (not fs::exists(testOutputDir))
			fs::create_directory(testOutputDir);
		
		fs::current_path(testOutputDir);
		
		NTHREADS = 1;
	}
	
	SkipList skip;
	if (vm.count("skip"))
		skip = readSkipList(vm["skip"].as<std::string>(), structure);
	
	FlipPeptides(structure, asymID,
		resFirst, resLast, skip, mm, trustDSSP, algorithm, mapWeight, plane5AtomsESD, cootScript, testMode);
	
	if (cootScript.is_open())
	{
		cootScript << "))" << std::endl;
		cootScript.close();
	}

	f.save(xyzout);
	
	return 0;
}
