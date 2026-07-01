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

#include <cif++/point.hpp>
#include <pdb-redo/MapMaker.hpp>

namespace pdb_redo
{

using DPoint = cif::point_type<double>;

// --------------------------------------------------------------------

class AtomLocationProvider;
class DFCollector;

// --------------------------------------------------------------------

using AtomRef = std::size_t;
using Xmap = typename Map<float>::Xmap;

// --------------------------------------------------------------------

struct Restraint
{
	virtual ~Restraint() = default;

	[[nodiscard]] virtual double f(const AtomLocationProvider &atoms) const = 0;
	virtual void df(const AtomLocationProvider &atoms, DFCollector &d) const = 0;
	[[nodiscard]] virtual std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const = 0;

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

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

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

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

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

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

	AtomRef mA, mB, mC, mD;
	int mPeriodicity;
	double mTarget, mESD;

  private:
	std::tuple<DPoint, DPoint, DPoint, DPoint> CalculateTorsionGradients(double theta, DPoint p[4]) const;
};

struct TransPeptideRestraint : public TorsionRestraint
{
	TransPeptideRestraint(AtomRef a, AtomRef b, AtomRef c, AtomRef d, double esd = 2.0)
		: TorsionRestraint(a, b, c, d, 180.0, esd, 2)
	{
	}
};

const double kChiralVolumeESD = 0.2; // according to coot that's a reasonable value...

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

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

	AtomRef mCentre, mA1, mA2, mA3;
	double mVolume, mESD = kChiralVolumeESD;
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

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

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

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

	AtomRef mA, mB;
	double mMinDist, mMinDistSq, mDistESD;
};

struct DensityRestraint : public Restraint
{
	DensityRestraint(std::vector<std::pair<AtomRef, double>> &&atoms,
		const Xmap &xMap, double mapWeight = 60);

	[[nodiscard]] double f(const AtomLocationProvider &atoms) const override;
	void df(const AtomLocationProvider &atoms, DFCollector &d) const override;
	[[nodiscard]] std::tuple<double, double> distortion(const AtomLocationProvider &atoms) const override;
	void print(const AtomLocationProvider &atoms) const override;

	std::vector<std::pair<AtomRef, double>> mAtoms;
	const Xmap &mXMap;
	double mMapWeight;
	bool mElectronScattering = false;
};

} // namespace pdb_redo