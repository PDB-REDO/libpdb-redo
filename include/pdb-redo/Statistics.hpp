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

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/MapMaker.hpp"

#include <cif++/datablock.hpp>
#include <pdb-redo/BondMap.hpp>

namespace pdb_redo
{

// --------------------------------------------------------------------

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

	StatsCollector(const MapMaker<float> &mm, cif::datablock &db, int modelNr, bool electronScattering);

	[[nodiscard]] virtual std::vector<ResidueStatistics> collect() const;

	[[nodiscard]] virtual std::vector<ResidueStatistics> collect(const std::string &asymID) const;

	// [[nodiscard]] virtual std::vector<ResidueStatistics> collect(const std::string &asymID,
	// 	int resFirst, int resLast, bool authNameSpace = false) const;

	// [[nodiscard]] virtual ResidueStatistics collect(std::initializer_list<const cif::mm::residue *> residues) const;

	// [[nodiscard]] virtual ResidueStatistics collect(std::initializer_list<cif::mm::atom> atoms) const;

	// [[nodiscard]] virtual ResidueStatistics collect(const std::vector<cif::mm::atom> &atoms) const;

  protected:
  // asym-seqid-authseqid-compid
	using residue_list = std::vector<std::tuple<std::string, int, std::string, std::string>>;

	std::vector<ResidueStatistics> collect(const residue_list &residues, BoundingBox &bbox, bool addWaters) const;

	void initialize();

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

	using GridPtDataMap = std::map<clipper::Coord_grid, double, cmpGPt>;

	struct AtomGridData
	{
		AtomGridData(const clipper::Coord_grid &gp, double density)
			: p(gp)
			, density(density)
		{
		}

		clipper::Coord_grid p;
		double density;
	};

	struct AtomDataSums
	{
		std::size_t ngrid = 0;
		double rfSums[2] = {}; // sums for R-Factor
		double edSums[2] = {}; // Sums for ED1 and ED3
		double ccSums[3] = {}; // Sums for CC calculation
		double rgSums[2] = {};
		double swSums[3] = {}; // Sums used for sample CC calculation

		AtomDataSums &operator+=(const AtomDataSums &rhs)
		{
			ngrid += rhs.ngrid;
			rfSums[0] += rhs.rfSums[0];
			rfSums[1] += rhs.rfSums[1];
			edSums[0] += rhs.edSums[0];
			edSums[1] += rhs.edSums[1];
			ccSums[0] += rhs.ccSums[0];
			ccSums[1] += rhs.ccSums[1];
			ccSums[2] += rhs.ccSums[2];
			rgSums[0] += rhs.rgSums[0];
			rgSums[1] += rhs.rgSums[1];
			swSums[0] += rhs.swSums[0];
			swSums[1] += rhs.swSums[1];
			swSums[2] += rhs.swSums[2];
			return *this;
		}

		[[nodiscard]] double cc() const
		{
			double s = (ccSums[1] - (edSums[0] * edSums[0]) / ngrid) * (ccSums[2] - (edSums[1] * edSums[1]) / ngrid);
			return (ccSums[0] - edSums[0] * edSums[1] / ngrid) / std::sqrt(s);
		}

		[[nodiscard]] double srg() const
		{
			double rgsq = rgSums[0] / rgSums[1];
			double rg = std::sqrt(rgsq);

			return std::sqrt(swSums[0] - rgsq * swSums[1] + 0.5 * rgsq * rgsq * swSums[2]) / (rg * rgSums[1]);
		}
	};

	struct AtomData
	{
		AtomData(cif::mm::atom atom, AtomShape shape)
			: atom(atom)
			, asymID(atom.get_label_asym_id())
			, seqID(atom.get_label_seq_id())
			, authSeqID(atom.get_auth_seq_id())
			, shape(std::move(shape))
			, radius(this->shape.radius())
			, occupancy(atom.get_occupancy())
		{
		}

		cif::mm::atom atom;
		std::string asymID;
		int seqID;
		std::string authSeqID; // required for waters
		AtomShape shape;
		float radius;
		float occupancy;
		std::vector<AtomGridData> points;
		double averageDensity = 0;
		double edia = 0;
		AtomDataSums sums;
	};

	// cif::mm::structure &mStructure;
	cif::datablock &mDb;
	int mModelNr;
	const MapMaker<float> &mMapMaker;

	clipper::Spacegroup mSpacegroup;
	clipper::Cell mCell;
	clipper::Grid_sampling mGrid;
	float mResHigh, mResLow;
	bool mElectronScattering;

	std::map<std::string, std::pair<double, double>> mRmsScaled;
	GridPtDataMap mGridPointDensity;
	std::map<std::string, std::vector<double>> mZScoresPerAsym;
	std::vector<AtomData> mAtomData;

	virtual void calculate(std::vector<AtomData> &atomData) const;
	void collectSums(std::vector<AtomData> &atomData, const GridPtDataMap &gridPointDensity) const;
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
	EDIAStatsCollector(const MapMaker<float> &mm, cif::datablock &db, int modelNr, bool electronScattering);

  protected:
	void calculate(std::vector<AtomData> &atomData) const override;

	BondMap createBondMap(std::vector<AtomData> &atomData) const;

	std::map<cif::atom_type, float> mRadii;
};

} // namespace pdb_redo
