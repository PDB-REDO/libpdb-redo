// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2025
// SPDX-License-Identifier: BSD-2-Clause

#pragma once


#include <cif++/cif++.hpp>
#include <cif++/symmetry.hpp>
#include <clipper/clipper.h>
#include <clipper/core/coords.h>

namespace pdb_redo
{

class BlobFinder
{
  public:
	// /// \brief Find all blobs in the map
	// BlobFinder(clipper::Xmap<float> &xmm, float growingPercentile = 0.95f);

	/// \brief Find only blobs near the molecule(s) in @a structure
	BlobFinder(clipper::Xmap<float> &xmm, cif::mm::structure &structure,
		float growingPercentile = 0.95f);

	std::vector<cif::point> next(float minimalVolume = 17);

  private:
	using GridPoint = clipper::Xmap<float>::Map_reference_coord;

	std::vector<GridPoint> pop();

	[[nodiscard]] bool blobIsInProximityOfAtoms(const std::vector<GridPoint> &blob) const;

	const clipper::Xmap<float> &mXmap;
	std::vector<GridPoint> mPotentialGridPoints;
	std::vector<cif::mm::atom> mProteinAtoms;
	std::vector<std::tuple<cif::point,float>> mResidueSpheres;
    cif::point mProteinCenter;
    cif::crystal mCrystal;
    float mProteinRadius;
	float mGrowingThreshold;
};

} // namespace pdb_redo