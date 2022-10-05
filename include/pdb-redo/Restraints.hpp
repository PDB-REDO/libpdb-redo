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

#include <set>

#include <cif++/point.hpp>

#include <pdb-redo/MapMaker.hpp>

namespace pdb_redo
{

using DPoint = cif::point_type<double>;

// --------------------------------------------------------------------

class AtomLocationProvider;
class DFCollector;

// --------------------------------------------------------------------

typedef uint32_t AtomRef;
typedef typename Map<float>::Xmap Xmap;

// --------------------------------------------------------------------

struct Restraint
{
	virtual ~Restraint() {}

	virtual double f(const AtomLocationProvider &atoms) const = 0;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const = 0;

	virtual void print(const AtomLocationProvider &atoms) const = 0;
};

struct BondRestraint : public Restraint
{
	BondRestraint(AtomRef a, AtomRef b, double distance, double esd)
		: mA(a)
		, mB(b)
		, mDist(distance)
		, mDistESD(esd)
	{
	}

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	AtomRef mA, mB;
	double mDist, mDistESD;
};

struct AngleRestraint : public Restraint
{
	AngleRestraint(AtomRef a, AtomRef b, AtomRef c, double angle, double esd)
		: mA(a)
		, mB(b)
		, mC(c)
		, mAngle(angle)
		, mESD(esd)
	{
	}

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	AtomRef mA, mB, mC;
	double mAngle, mESD;
};

struct TorsionRestraint : public Restraint
{
	TorsionRestraint(AtomRef a, AtomRef b, AtomRef c, AtomRef d, double target, double esd, int periodicity)
		: mA(a)
		, mB(b)
		, mC(c)
		, mD(d)
		, mPeriodicity(periodicity)
		, mTarget(target)
		, mESD(esd)
	{
	}

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	AtomRef mA, mB, mC, mD;
	int mPeriodicity;
	double mTarget, mESD;

  private:
	std::tuple<DPoint, DPoint, DPoint, DPoint> CalculateTorsionGradients(float theta, DPoint p[4]) const;
};

struct TransPeptideRestraint : public TorsionRestraint
{
	TransPeptideRestraint(AtomRef a, AtomRef b, AtomRef c, AtomRef d, double esd = 2.0)
		: TorsionRestraint(a, b, c, d, 180.0, esd, 2)
	{
	}
};

struct ChiralVolumeRestraint : public Restraint
{
	ChiralVolumeRestraint(AtomRef c, AtomRef a1, AtomRef a2, AtomRef a3, double volume)
		: mCentre(c)
		, mA1(a1)
		, mA2(a2)
		, mA3(a3)
		, mVolume(volume)
	{
	}

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	AtomRef mCentre, mA1, mA2, mA3;
	double mVolume;
};

struct PlanarityRestraint : public Restraint
{
	PlanarityRestraint(std::vector<AtomRef> &&atoms, double esd)
		: mAtoms(std::move(atoms))
		, mESD(esd)
	{
		if (mAtoms.size() < 3)
			throw std::runtime_error("Insufficient number of atoms in planar restraint");
	}

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	void calculatePlaneFunction(const AtomLocationProvider &atoms, double abcd[4]) const;

	std::vector<AtomRef> mAtoms;
	double mESD;
};

struct NonBondedContactRestraint : public Restraint
{
	NonBondedContactRestraint(AtomRef a, AtomRef b, double minDist, double esd)
		: mA(a)
		, mB(b)
		, mMinDist(minDist)
		, mMinDistSq(minDist * minDist)
		, mDistESD(esd)
	{
	}

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	AtomRef mA, mB;
	double mMinDist, mMinDistSq, mDistESD;
};

struct DensityRestraint : public Restraint
{
	DensityRestraint(std::vector<std::pair<AtomRef, double>> &&atoms,
		const Xmap &xMap, double mapWeight = 60);

	virtual double f(const AtomLocationProvider &atoms) const;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const;
	virtual void print(const AtomLocationProvider &atoms) const;

	std::vector<std::pair<AtomRef, double>> mAtoms;
	const Xmap &mXMap;
	double mMapWeight;
	bool mElectronScattering = false;
};

} // namespace pdb_redo