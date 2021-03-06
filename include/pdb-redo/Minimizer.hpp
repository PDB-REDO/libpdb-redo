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

#pragma once

#include <boost/optional.hpp>

#include "cif++/Structure.hpp"

#include "pdb-redo/MapMaker.hpp"
#include "cif++/BondMap.hpp"
#include "pdb-redo/Restraints.hpp"

// --------------------------------------------------------------------

class AtomLocationProvider
{
  public:
	AtomLocationProvider(const AtomLocationProvider&) = delete;
	AtomLocationProvider& operator=(const AtomLocationProvider&) = delete;

	AtomLocationProvider(std::vector<mmcif::Atom>& atoms)
		: mAtoms(atoms) {}
	virtual ~AtomLocationProvider() {}
	
	virtual mmcif::DPoint operator[](AtomRef atomID) const;
	virtual std::string atom(AtomRef atomID) const;

  protected:
	std::vector<mmcif::Atom>& mAtoms;
};

// --------------------------------------------------------------------

class DFCollector
{
  public:
	DFCollector(const DFCollector&) = delete;
	DFCollector& operator=(const DFCollector&) = delete;

	DFCollector() {}
	virtual ~DFCollector() {}
	
	virtual void add(AtomRef atom, double dx, double dy, double dz) = 0;
	void add(AtomRef atom, mmcif::DPoint&& d)
	{
		add(atom, d.mX, d.mY, d.mZ);
	}
};

// --------------------------------------------------------------------

class Minimizer
{
  public:
	typedef clipper::Xmap<float> XMap;

	Minimizer(const Minimizer&) = delete;
	Minimizer& operator=(const Minimizer&) = delete;

	virtual ~Minimizer() {}

	// factory method:
	static Minimizer* create(const std::string& algorithm,
		const mmcif::Polymer& poly, int first, int last, const mmcif::BondMap& bm,
		const XMap& xMap, float mapWeight = 60, float plane5AtomsESD = 0.11);

	void printStats();

	virtual double refine(bool storeAtoms) = 0;
	double score();
	virtual std::vector<std::pair<std::string,mmcif::Point>> getAtoms() const = 0;
	virtual void storeAtomLocations() = 0;

  protected:

	Minimizer(const mmcif::Polymer& poly, int first, int last,
		const mmcif::BondMap& bm, const XMap& xMap, float mapWeight,
		float plane5AtomsESD);

	double score(const AtomLocationProvider& loc);

	void addLinkRestraints(const mmcif::Monomer& a, const mmcif::Monomer& b, const std::string& linkName)
	{
		addLinkRestraints(a, b, mmcif::Link::create(linkName));
	}
	
	void addLinkRestraints(const mmcif::Monomer& a, const mmcif::Monomer& b, const mmcif::Link& link);

	template<typename R>
	double rmsz(const AtomLocationProvider& atoms, const std::vector<R>& a) const
	{
		double result = 0;
	
		if (not a.empty())
		{
			double sumZ = accumulate(a.begin(), a.end(),
				0.0, [&atoms](double sum, const R& r) { double z = r.f(atoms); return sum + z; });
			
			result = sqrt(sumZ / a.size());
		}
		
		return result;
	}

	AtomRef ref(const mmcif::Atom& atom);

	bool mElectronScattering = false;	// TODO: use!

	std::vector<mmcif::Atom> mAtoms, mReferencedAtoms;
	std::vector<size_t> mRef2AtomIndex;
	std::map<std::string,AtomRef> mRefIndex;
	
	std::vector<BondRestraint> mBondRestraints;
	std::vector<AngleRestraint> mAngleRestraints;
	std::vector<TorsionRestraint> mTorsionRestraints;
	std::vector<TransPeptideRestraint> mTransPeptideRestraints;
	std::vector<ChiralVolumeRestraint> mChiralVolumeRestraints;
	std::vector<PlanarityRestraint> mPlanarityRestraints;
	std::vector<NonBondedContactRestraint> mNonBondedContactRestraints;
	std::unique_ptr<DensityRestraint> mDensityRestraint;
	
	std::vector<Restraint*> mRestraints;
};

