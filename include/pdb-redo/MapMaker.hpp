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

#include <clipper/clipper.h>

#include <cif++/cif++.hpp>

#include <filesystem>

// My apologies, but this code is emitting way too many warnings...
#if defined(_MSC_VER)
# pragma warning(disable : 4244) // possible loss of data (in conversion to smaller type)
#endif

namespace pdb_redo
{

template <typename FTYPE = float>
class Map
{
  public:
	using ftype = FTYPE;
	using Xmap = typename clipper::Xmap<ftype>;

	Map();
	Map(const Map &rhs) = default;
	~Map();

	Map &operator=(const Map &rhs) = default;

	void calculateStats();

	[[nodiscard]] double rmsDensity() const { return mRMSDensity; }
	[[nodiscard]] double meanDensity() const { return mMeanDensity; }

	operator Xmap &() { return mMap; }
	operator const Xmap &() const { return mMap; }
	Xmap &get() { return mMap; }
	[[nodiscard]] const Xmap &get() const { return mMap; }

	// These routines work with CCP4 map files
	void read(const std::filesystem::path &f);
	void write(const std::filesystem::path &f);

	void write_masked(std::ostream &os, clipper::Grid_range range);
	void write_masked(const std::filesystem::path &f,
		clipper::Grid_range range);

	[[nodiscard]] clipper::Spacegroup spacegroup() const { return mMap.spacegroup(); }
	[[nodiscard]] clipper::Cell cell() const { return mMap.cell(); }

	/// \brief Create a masked map blotting out the density for all \a atom_ids in the structure contained in \a db
	[[nodiscard]] Map masked(const cif::mm::structure &structure, const std::vector<cif::mm::atom> &atom_ids) const;

	/// \brief Return the z-weighted density sum for the atoms \a atom_ids in the structure contained in \a db
	[[nodiscard]] float z_weighted_density(const cif::mm::structure &structure, const std::vector<cif::mm::atom> &atom_ids) const;

  private:
	Xmap mMap;
	double mMinDensity, mMaxDensity;
	double mRMSDensity, mMeanDensity;
};

// --------------------------------------------------------------------

bool IsMTZFile(const std::filesystem::path &p);

// --------------------------------------------------------------------

template <typename FTYPE = float>
class MapMaker
{
  public:
	using MapType = Map<FTYPE>;
	using Xmap = typename MapType::Xmap;

	enum AnisoScalingFlag
	{
		as_None,
		as_Observed,
		as_Calculated
	};

	MapMaker();
	~MapMaker();

	MapMaker(const MapMaker &) = delete;
	MapMaker &operator=(const MapMaker &) = delete;

	void loadMTZ(const std::filesystem::path &mtzFile,
		float samplingRate,
		std::initializer_list<std::string> fbLabels = { "FWT", "PHWT" },
		std::initializer_list<std::string> fdLabels = { "DELFWT", "PHDELWT" },
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> fcLabels = { "FC_ALL", "PHIC_ALL" },
		std::initializer_list<std::string> faLabels = { "FAN", "PHAN" });

	void loadMaps(
		const std::filesystem::path &fbMapFile,
		const std::filesystem::path &fdMapFile,
		float reshi, float reslo);

	// following works on both mtz files and structure factor files in CIF format
	void calculate(const std::filesystem::path &hklin,
		const cif::mm::structure &structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate, bool electronScattering = false,
		std::initializer_list<std::string> foLabels = { "FP", "SIGFP" },
		std::initializer_list<std::string> freeLabels = { "FREE" });

	void recalc(const cif::mm::structure &structure,
		bool noBulk, AnisoScalingFlag anisoScaling,
		float samplingRate, bool electronScattering = false);

	void printStats();

	void writeMTZ(const std::filesystem::path &file,
		const std::string &project, const std::string &crystal);

	MapType &fb() { return mFb; }
	MapType &fd() { return mFd; }
	MapType &fa() { return mFa; }

	[[nodiscard]] const MapType &fb() const { return mFb; }
	[[nodiscard]] const MapType &fd() const { return mFd; }
	[[nodiscard]] const MapType &fa() const { return mFa; }

	[[nodiscard]] float resLow() const { return mResLow; }
	[[nodiscard]] float resHigh() const { return mResHigh; }

	[[nodiscard]] const clipper::Spacegroup &spacegroup() const { return mHKLInfo.spacegroup(); }
	[[nodiscard]] const clipper::Cell &cell() const { return mHKLInfo.cell(); }
	[[nodiscard]] const clipper::Grid_sampling &gridSampling() const { return mGrid; }

  private:
	void loadFoFreeFromReflectionsFile(const std::filesystem::path &hklin);
	void loadFoFreeFromMTZFile(const std::filesystem::path &hklin,
		std::initializer_list<std::string> foLabels,
		std::initializer_list<std::string> freeLabels);

	void fixMTZ();

	MapType mFb, mFd, mFa;
	clipper::Grid_sampling mGrid;
	float mResLow, mResHigh;
	int mNumRefln = 1000, mNumParam = 20;

	// Cached raw data
	clipper::HKL_info mHKLInfo;
	clipper::HKL_data<clipper::data32::F_sigF> mFoData;
	clipper::HKL_data<clipper::data32::Flag> mFreeData;
	clipper::HKL_data<clipper::data32::F_phi> mFcData, mFbData, mFdData, mFaData;
	clipper::HKL_data<clipper::data32::Phi_fom> mPhiFomData;
};

} // namespace pdb_redo
