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

#pragma once

#include <pdb-redo/BondMap.hpp>
#include "pdb-redo/DistanceMap.hpp"
#include "pdb-redo/MapMaker.hpp"

namespace pdb_redo
{

// --------------------------------------------------------------------

struct AtomData;
class BoundingBox;

struct ResidueStatistics
{
	std::string asymID;
	int seqID;
	std::string compID;
	std::string authSeqID;

	double RSR, SRSR, RSCCS, EDIAm, OPIA;
	int ngrid;
};

std::ostream &operator<<(std::ostream &os, const ResidueStatistics &st);

// --------------------------------------------------------------------

template <typename F, typename FTYPE = float>
void iterateGrid(const clipper::Coord_orth &p, float r, const clipper::Xmap<FTYPE> &m, F &&func)
{
	using namespace clipper;

	Coord_frac fp = p.coord_frac(m.cell());

	Coord_frac o = Coord_orth(r, r, r).coord_frac(m.cell());
	o[0] = std::abs(o[0]);
	o[1] = std::abs(o[1]);
	o[2] = std::abs(o[2]);

	Coord_frac fMin = fp - o, fMax = fp + o;
	Coord_map mMin = fMin.coord_map(m.grid_sampling()), mMax = fMax.coord_map(m.grid_sampling());
	Coord_grid gMin = mMin.floor(), gMax = mMax.ceil();

	auto i0 = Xmap_base::Map_reference_coord(m, gMin);
	for (auto iu = i0; iu.coord().u() <= gMax[0]; iu.next_u())
		for (auto iv = iu; iv.coord().v() <= gMax[1]; iv.next_v())
			for (auto iw = iv; iw.coord().w() <= gMax[2]; iw.next_w())
				func(iw);
}

// --------------------------------------------------------------------

class StatsCollector
{
  public:
	StatsCollector(const StatsCollector &) = delete;
	StatsCollector &operator=(const StatsCollector &) = delete;

	StatsCollector(const MapMaker<float> &mm,
		cif::mm::structure &structure, bool electronScattering);

	virtual std::vector<ResidueStatistics> collect() const;

	virtual std::vector<ResidueStatistics> collect(const std::string &asymID,
		int resFirst, int resLast, bool authNameSpace = true) const;

	virtual ResidueStatistics collect(std::initializer_list<const cif::mm::residue *> residues) const;

	virtual ResidueStatistics collect(std::initializer_list<cif::mm::atom> atoms) const;

	virtual ResidueStatistics collect(const std::vector<cif::mm::atom> &atoms) const;

  protected:
	using residue_list = std::vector<std::tuple<std::string, int, std::string>>;

	// asym-seqid-compid
	std::vector<ResidueStatistics> collect(const residue_list &residues, BoundingBox &bbox, bool addWaters) const;

	void initialize();

	virtual void calculate(std::vector<AtomData> &atomData) const;

	struct cmpGPt
	{
		bool operator()(const clipper::Coord_grid &a, const clipper::Coord_grid &b) const
		{
			int d = a.u() - b.u();
			if (d == 0)
				d = a.v() - b.v();
			if (d == 0)
				d = a.w() - b.w();
			return d < 0;
		}
	};

	typedef std::map<clipper::Coord_grid, double, cmpGPt> GridPtDataMap;

	cif::mm::structure &mStructure;
	const MapMaker<float> &mMapMaker;

	clipper::Spacegroup mSpacegroup;
	clipper::Cell mCell;
	clipper::Grid_sampling mGrid;
	float mResHigh, mResLow;
	bool mElectronScattering;

	std::map<std::string, std::pair<double, double>> mRmsScaled;

	void collectSums(std::vector<AtomData> &atomData, GridPtDataMap &gridPointDensity) const;
	void sumDensity(std::vector<AtomData> &atomData,
		GridPtDataMap &gridPointDensity, std::map<std::string, std::vector<double>> &zScoresPerAsym) const;

	// Other variables we cache

	double mMeanDensityFb, mRMSDensityFb, mRMSDensityFd;
	double mSZ; // average electron density in cell
	double mVF; // degrees of freedom
	double mVC; // cell volume?
};

// --------------------------------------------------------------------

class EDIAStatsCollector : public StatsCollector
{
  public:
	EDIAStatsCollector(MapMaker<float> &mm,
		cif::mm::structure &structure, bool electronScattering,
		const BondMap &bondMap);

  protected:
	virtual void calculate(std::vector<AtomData> &atomData) const;

	DistanceMap mDistanceMap;
	const BondMap &mBondMap;
	std::map<cif::atom_type, float> mRadii;
};

} // namespace pdb_redo
