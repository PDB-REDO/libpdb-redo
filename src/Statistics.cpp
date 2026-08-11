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

#include "pdb-redo/Statistics.hpp"

#include "pdb-redo/BondMap.hpp"
#include "pdb-redo/DistanceMap.hpp"

#include <algorithm>
#include <cif++/cif++.hpp>
#include <cif++/datablock.hpp>
#include <cif++/symmetry.hpp>
#include <cif++/utilities.hpp>
#include <cif++/validate.hpp>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

// --------------------------------------------------------------------

namespace pdb_redo
{

using cif::atom_type_traits;
using std::unary_function;

// --------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const ResidueStatistics &st)
{
	if (st.compID == "HOH")
		os << st.asymID << '_' << st.authSeqID << '_' << st.compID << '\t';
	else
		os << st.asymID << '_' << st.seqID << '_' << st.compID << '\t';

	os << st.RSR << '\t'
	   << st.SRSR << '\t'
	   << st.RSCCS << '\t'
	   << st.ngrid << '\t'
	   << st.EDIAm << '\t'
	   << st.OPIA;

	return os;
}

// --------------------------------------------------------------------

double anorm(double x)
{
	return 0.5 * erfc(-x * std::sqrt(0.5));
}

double phinvs(double p)
{
	//
	// ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3.
	//
	// Produces the normal deviate Z corresponding to a given lower tail
	// area of P; Z is accurate to about 1 part in 10**16.

	// Coefficients for P close to 0.5
	const double A[8] = {
		3.3871328727963666080, 1.3314166789178437745e+2, 1.9715909503065514427e+3, 1.3731693765509461125e+4,
		4.5921953931549871457e+4, 6.7265770927008700853e+4, 3.3430575583588128105e+4, 2.5090809287301226727e+3
	},
				 B[8] = { 0, 4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3, 2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4, 5.2264952788528545610e+3 };

	// Coefficients for P not close to 0, 0.5 or 1.
	const double C[8] = {
		1.42343711074968357734e0,
		4.63033784615654529590e0,
		5.76949722146069140550e0,
		3.64784832476320460504e0,
		1.27045825245236838258e0,
		2.41780725177450611770e-1,
		2.27238449892691845833e-2,
		7.74545014278341407640e-4,
	},
				 D[8] = { 0, 2.05319162663775882187e0, 1.67638483018380384940e0, 6.89767334985100004550e-1, 1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9 };

	// Coefficients for P near 0 or 1.
	const double E[8] = {
		6.65790464350110377720e0,
		5.46378491116411436990e0,
		1.78482653991729133580e0,
		2.96560571828504891230e-1,
		2.65321895265761230930e-2,
		1.24266094738807843860e-3,
		2.71155556874348757815e-5,
		2.01033439929228813265e-7,
	},
				 F[8] = { 0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2, 7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15 };

	if (p < 0 or p > 1)
		throw std::runtime_error("P should be >=0 and <=1");

	double q = p - 0.5;
	double result;

	if (std::abs(q) < 0.425)
	{
		double r = 0.180625e0 - q * q;
		result =
			q * (((((((A[7] * r + A[6]) * r + A[5]) * r + A[4]) * r + A[3]) * r + A[2]) * r + A[1]) * r + A[0]) / (((((((B[7] * r + B[6]) * r + B[5]) * r + B[4]) * r + B[3]) * r + B[2]) * r + B[1]) * r + 1);
	}
	else
	{
		double r;

		if (q < 0)
			r = p;
		else
			r = 1 - p;

		r = std::sqrt(-std::log(r));

		if (r <= 5)
		{
			r -= 1.6;
			result = (((((((C[7] * r + C[6]) * r + C[5]) * r + C[4]) * r + C[3]) * r + C[2]) * r + C[1]) * r + C[0]) / (((((((D[7] * r + D[6]) * r + D[5]) * r + D[4]) * r + D[3]) * r + D[2]) * r + D[1]) * r + 1);
		}
		else
		{
			r -= 0.5;
			result = (((((((E[7] * r + E[6]) * r + E[5]) * r + E[4]) * r + E[3]) * r + E[2]) * r + E[1]) * r + E[0]) / (((((((F[7] * r + F[6]) * r + F[5]) * r + F[4]) * r + F[3]) * r + F[2]) * r + F[1]) * r + 1);
		}

		if (q < 0)
			result = -result;
	}
	return result;
}

double errsol(double a)
{
	auto c = std::sqrt(2.0 / std::numbers::pi);
	auto b = std::abs(a);

	double result = 0;
	if (b > 3 / c)
	{
		auto x = std::abs(std::pow(b, 1 / 3.0) - 2 * std::pow(std::numbers::pi / b, 2));
		if (a < 0)
			x = -x;

		for (;;)
		{
			auto xx = x * x;
			auto y = c * std::exp(-0.5 * xx);
			auto d = (b * (2 * anorm(x) - 1 - x * y) / xx - x) / (b * y - 3);
			x -= d;

			if (std::abs(d) <= 1e-4)
				break;
		}

		result = x;
	}

	return result;
}

// --------------------------------------------------------------------

class PointWeightFunction
{
  public:
	PointWeightFunction(cif::point center, float atomRadius)
		: m_Center(center)
		, m_Radius(atomRadius)
	{
		m_P[0] = P{ -1.0f, 0, 1.0f, 1.0822f };
		m_P[1] = P{ 5.1177f, 1.29366f, -0.4f, 1.4043f };
		m_P[2] = P{ -0.9507f, 2, 0, 2 };
	}

	float operator()(cif::point p) const
	{
		float d = distance(m_Center, p);
		d /= m_Radius;

		float result = 0;

		for (auto &pi : m_P)
		{
			if (d > pi.x)
				continue;

			result = pi.m * (d - pi.c) * (d - pi.c) + pi.b;

			//			assert(result != 0);
			if (result == 0)
				result = std::numeric_limits<float>::epsilon();

			break;
		}

		return result;
	}

  private:
	struct P
	{
		float m, c, b, x;
	};

	cif::point m_Center;
	float m_Radius;
	P m_P[3];
};

// --------------------------------------------------------------------

std::tuple<float, float> CalculateMapStatistics(const clipper::Xmap<float> &f)
{
	double sum = 0, sum2 = 0;
	int count = 0;

	for (auto ix = f.first(); not ix.last(); ix.next())
	{
		auto v = f[ix];

		if (std::isnan(v))
			throw std::runtime_error("map contains NaN values");

		++count;
		sum += v;
		sum2 += v * v;
	}

	auto meanDensity = static_cast<float>(sum / count);
	auto rmsDensity = static_cast<float>(std::sqrt((sum2 / count) - (meanDensity * meanDensity)));

	return std::make_tuple(meanDensity, rmsDensity);
}

// --------------------------------------------------------------------

class BoundingBox
{
  public:
	BoundingBox(float margin = 5.f)
		: mMargin(margin)
		, mMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max())
		, mMax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min())
	{
	}

	void extend(cif::point pt)
	{
		if (mMin.m_x > pt.m_x - mMargin)
			mMin.m_x = pt.m_x - mMargin;
		if (mMin.m_y > pt.m_y - mMargin)
			mMin.m_y = pt.m_y - mMargin;
		if (mMin.m_z > pt.m_z - mMargin)
			mMin.m_z = pt.m_z - mMargin;
		if (mMax.m_x < pt.m_x + mMargin)
			mMax.m_x = pt.m_x + mMargin;
		if (mMax.m_y < pt.m_y + mMargin)
			mMax.m_y = pt.m_y + mMargin;
		if (mMax.m_z < pt.m_z + mMargin)
			mMax.m_z = pt.m_z + mMargin;
	}

	[[nodiscard]] bool contains(const cif::point &p) const
	{
		return p.m_x >= mMin.m_x and p.m_x <= mMax.m_x and p.m_y >= mMin.m_y and p.m_y <= mMax.m_y and p.m_z >= mMin.m_z and p.m_z <= mMax.m_z;
	}

  private:
	float mMargin;
	cif::point mMin, mMax;
};
// --------------------------------------------------------------------

StatsCollector::StatsCollector(const MapMaker<float> &mm, cif::datablock &db, int modelNr, bool electronScattering)
	: mDb(db)
	, mModelNr(modelNr)
	, mMapMaker(mm)
	, mElectronScattering(electronScattering)
{
	mSpacegroup = mm.spacegroup();
	mCell = mm.cell();
	mGrid = mm.gridSampling();
	mResHigh = static_cast<float>(mm.resHigh());
	mResLow = static_cast<float>(mm.resLow());

	initialize();
}

void StatsCollector::initialize()
{
	// easiest way to prime this map:
	for (auto asym_id : mDb["struct_asym"].rows<std::string>("id"))
		mRmsScaled[asym_id] = { 1, 1 };

	mMeanDensityFb = mMapMaker.fb().meanDensity();
	mRMSDensityFb = mMapMaker.fb().rmsDensity();
	mRMSDensityFd = mMapMaker.fd().rmsDensity();

	// calculate degrees of freedom
	auto omcd = mCell.matrix_orth();

	mVF = 1;
	mVC = 1;

	for (int i = 0; i < 3; ++i)
	{
		mVC *= omcd(i, i);
		mVF *= omcd(i, i) / mGrid[i];
	}

	mVF *= std::pow(2 / mResHigh, 3);

	mSZ = 0;
	//	double so = 0;
	//	const double C = std::sqrt(2.0 / std::numbers::pi);

	using namespace cif::literals;

	{
		cif::progress_bar progress(mDb["atom_site"].size(), "Initialising stats collector");

		for (auto rh : mDb["atom_site"].find("pdbx_PDB_model_num"_key == mModelNr or "pdbx_PDB_model_num"_key == cif::null))
		{
			const auto [id, type, occupancy] = rh.get<std::string, std::string, float>("id", "type_symbol", "occupancy");

			progress.consumed(1);

			auto t = cif::atom_type_traits(type).type();
			if (t <= cif::atom_type::He)
				continue;

			float w = occupancy * static_cast<int>(t);

			if (w <= 0)
				continue;

			mSZ += w;

			cif::mm::atom atom(const_cast<cif::datablock &>(mDb), rh);
			AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);
			mAtomData.emplace_back(std::move(atom), std::move(shape));
		}
	}

	//	auto bo = mSZ;

	mSZ = mSZ * mSpacegroup.num_symops() / mVC;
	//	mMeanBIso = std::pow(mResHigh * errsol(bo / so), 2);

	// Calculate overall rms data
	std::map<std::string, std::vector<double>> zScoresPerAsym;
	sumDensity(mAtomData, mGridPointDensity, zScoresPerAsym);

	// Now that we have the density data, we can calculate the correction/rescale factors
	for (auto zsc : zScoresPerAsym)
	{
		// collect array of z-scores
		std::vector<double> &zdca0 = zsc.second;

		auto &z = zdca0;
		auto vf = mVF;

		std::ranges::sort(z);

		double qa = 0, qb = 1;

		std::size_t nd = z.size();
		auto n = static_cast<std::size_t>(round(vf * nd));

		if (n > 100)
		{
			std::size_t i1 = static_cast<std::size_t>((n + 1) * anorm(-1.5)) + 1;
			auto i2 = static_cast<std::size_t>((n + 1) * anorm(1.5));

			std::size_t ns = i2 - i1 + 1;

			double vr = (nd - 1) / (n - 1.0);
			double sw = 0, swx = 0, swxs = 0, swy = 0, swxy = 0, swys = 0;

			for (auto i = i1; i <= i2; ++i)
			{
				double qx = phinvs(static_cast<double>(i) / (n + 1));
				double x = vr * i;
				auto j = static_cast<std::size_t>(x);
				x -= j;

				//		assert(j < z.size());
				if (j < 1 or j >= z.size())
					continue;

				auto qyd = (1.0 - x) * z[j - 1] + x * z[j] - qx;

				auto wx = std::exp(-0.5 * qx * qx);
				sw += wx;
				swx += wx * qx;
				swxs += wx * qx * qx;
				swy += wx * qyd;
				swxy += wx * qx * qyd;
				swys += wx * qyd * qyd;
			}

			double dd = 1.0 / (sw * swxs - swx * swx);
			qa = dd * (swxs * swy - swx * swxy);
			qb = dd * (sw * swxy - swx * swy);

			if (cif::VERBOSE > 1)
			{
				swys = dd * (swys - (qa * swy + qb * swxy)) / (ns - 2);
				std::cerr << '\n'
						  << "Intercept & gradient before LS: " << qa << " (" << std::sqrt(swys * swxs) << ") " << qb << " (" << std::sqrt(swys * sw) << ')' << '\n';
			}

			qb += 1.0;

			if (cif::VERBOSE > 1)
			{
				std::cerr << '\n'
						  << "Rescale SD(delta-rho) using Q-Q plot for asym " << zsc.first << ':' << '\n'
						  << std::string(54, '=') << '\n'
						  << "Input & updated SD(delta-rho): " << mRMSDensityFd << " ; " << qb * mRMSDensityFd << '\n'
						  << '\n';
			}
		}

		mRmsScaled[zsc.first] = std::make_pair(qa * mRMSDensityFd, qb * mRMSDensityFd);
	}
}

std::vector<ResidueStatistics> StatsCollector::collect() const
{
	ResidueList residues;
	BoundingBox bbox;

	for (auto atom : mAtomData | std::views::transform(&AtomData::atom))
	{
		if (atom.is_water())
			continue;

		PerResidueInfo pr{ atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id(), atom.get_label_comp_id() };

		if (residues.empty() or residues.back() != pr)
			residues.emplace_back(std::move(pr));

		if (auto altID = atom.get_label_alt_id(); not altID.empty())
			residues.back().altIDs.insert(altID);

		bbox.extend(atom.get_location());
	}

	return collect(residues, bbox, true);
}

std::vector<ResidueStatistics> StatsCollector::collect(const std::string &asymID) const
{
	using namespace std::literals;

	ResidueList residues;
	BoundingBox bbox;

	for (auto atom : mAtomData                                    //
						 | std::views::transform(&AtomData::atom) //
#if 0
						 | std::views::filter([asymID](auto &a)
							   { return a.get_label_asym_id() == asymID; })
#endif

	)
	{
		if (atom.get_label_asym_id() != asymID)
			continue;

		PerResidueInfo pr{ atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id(), atom.get_label_comp_id() };

		if (residues.empty() or residues.back() != pr)
			residues.emplace_back(std::move(pr));

		if (auto altID = atom.get_label_alt_id(); not altID.empty())
			residues.back().altIDs.insert(altID);

		bbox.extend(atom.get_location());
	}

	return collect(residues, bbox, false);
}

// std::vector<ResidueStatistics> StatsCollector::collect(const std::string &asymID, int resFirst, int resLast, bool authNameSpace) const
// {
// 	residue_list residues;
// 	std::vector<cif::mm::atom> atoms;

// 	// for (auto atom : mStructure.atoms())
// 	// {
// 	// 	if (atom.is_water())
// 	// 		continue;

// 	// 	if (authNameSpace)
// 	// 	{
// 	// 		int auth_seq_id = stoi(atom.get_auth_seq_id());

// 	// 		if (atom.get_auth_asym_id() != asymID or auth_seq_id < resFirst or auth_seq_id > resLast)
// 	// 			continue;
// 	// 	}
// 	// 	else
// 	// 	{
// 	// 		if (atom.get_label_asym_id() != asymID or atom.get_label_seq_id() < resFirst or atom.get_label_seq_id() > resLast)
// 	// 			continue;
// 	// 	}

// 	// 	auto k = std::make_tuple(atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());

// 	// 	if (residues.empty() or residues.back() != k)
// 	// 	{
// 	// 		residues.emplace_back(move(k));
// 	// 		atoms.emplace_back(std::move(atom));
// 	// 	}
// 	// }

// 	for (auto atom : mStructure.atoms())
// 	{
// 		if (atom.is_water())
// 			continue;

// 		if (authNameSpace)
// 		{
// 			int auth_seq_id = stoi(atom.get_auth_seq_id());

// 			if (atom.get_auth_asym_id() != asymID or auth_seq_id < resFirst or auth_seq_id > resLast)
// 				continue;
// 		}
// 		else
// 		{
// 			if (atom.get_label_asym_id() != asymID or atom.get_label_seq_id() < resFirst or atom.get_label_seq_id() > resLast)
// 				continue;
// 		}

// 		auto k = std::make_tuple(atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());

// 		if (residues.empty() or residues.back() != k)
// 			residues.emplace_back(std::move(k));
// 	}

// 	for (const auto &[asymID, seqID, authSeqID] : residues)
// 	{
// 		auto &res = mStructure.get_residue(asymID, seqID, authSeqID);

// 		for (auto atom : res.unique_atoms())
// 			atoms.push_back(atom);
// 	}

// 	BoundingBox bbox(mStructure, atoms, 5.0f);
// 	return collect(residues, bbox, false);
// }

std::vector<ResidueStatistics> StatsCollector::collect(const ResidueList &residues, BoundingBox &bbox, bool addWaters) const
{
	std::vector<AtomData> atomData;

	for (auto ad : mAtomData)
	{
		if (ad.atom.is_water())
		{
			if (not addWaters)
				continue;
		}

		if (not bbox.contains(ad.atom.get_location()))
			continue;

		atomData.emplace_back(ad);
	}

	calculate(atomData);

	std::set<std::string> missing;
	std::vector<ResidueStatistics> result;

	cif::progress_bar progress(residues.size(), "Collecting per residue");

	// And now collect the per residue information
	for (const auto &[asymID, seqID, authSeqID, compID, altIDs] : residues)
	{
		// TODO: Need to do something with hetero residues (alternate compound types)
		// auto &res = mStructure.get_residue(asymID, seqID, authSeqID);

		AtomDataSums sums;

		std::vector<const AtomData *> resAtomData;
		for (const auto &d : atomData)
		{
			if (d.asymID == asymID and d.seqID == seqID and d.authSeqID == authSeqID)
				resAtomData.push_back(&d);
		}

		std::vector<std::string> atomIDs;

		try
		{
			if (not missing.count(compID))
			{
				atomIDs = BondMap::atomIDsForCompound(compID);

				for (auto &compAtom : atomIDs)
				{
					if (compAtom == "OXT")
						continue;

					for (auto d : resAtomData)
					{
						if (d->atom.get_label_atom_id() != compAtom)
							continue;

						// We used to factor in the occupancy here, but that seems to be incorrect
						// since occupancy was already used.
						// sums += d->sums * d->occupancy;

						sums += d->sums;
						break;
					}

					std::erase_if(resAtomData, [id = compAtom](const AtomData *d)
						{ return d->atom.get_label_atom_id() == id; });
				}
			}
		}
		catch (const BondMapException &ex)
		{
			std::cerr << "Missing information for compound '" << compID << '\'' << '\n'
					  << ex.what() << '\n';

			missing.insert(compID);
		}

		// atoms that were present but not part of the Compound
		for (auto d : resAtomData)
			sums += d->sums;

		// EDIA
		std::set<std::string> alts = altIDs;
		if (alts.empty())
			alts.insert("");

		double EDIAm = 0, OPIA = 0, OCC = 0;

		for (auto alt : alts)
		{
			double ediaSum = 0, occSum = 0;
			std::size_t n = 0, m = 0;

			for (const auto &d : atomData)
			{
				if (d.asymID != asymID or d.seqID != seqID or d.authSeqID != authSeqID)
					continue;

				if (alt.empty())
				{
					resAtomData.push_back(&d);
					continue;
				}

				auto altd = d.atom.get_label_alt_id();
				if (altd.empty() or altd == alt)
					resAtomData.push_back(&d);
			}

			if (atomIDs.empty())
			{
				for (const auto &d : resAtomData)
				{
					occSum += d->occupancy;
					ediaSum += std::pow(d->edia + 0.1, -2);
					++n;
					if (d->edia >= 0.8)
						++m;
				}
			}
			else
			{
				for (auto &compAtom : atomIDs)
				{
					if (compAtom == "OXT")
						continue;

					for (auto d : resAtomData)
					{
						if (d->atom.get_label_atom_id() != compAtom)
							continue;

						occSum += d->occupancy;
						ediaSum += std::pow(d->edia + 0.1, -2);

						++n;
						if (d->edia >= 0.8)
							++m;
					}
				}
			}

			if (n == 0) // I'm paranoid
				continue;

			OCC += occSum;
			EDIAm += occSum * (1 / std::sqrt(ediaSum / n) - 0.1);
			OPIA += occSum * (100. * m / n);
		}

		if (atomIDs.empty())
		{
			EDIAm = std::nan("0");
			OPIA = std::nan("0");
		}
		else
		{
			EDIAm /= OCC;
			OPIA /= OCC;
		}

		result.emplace_back(ResidueStatistics{ asymID, seqID, compID,
			authSeqID,
			(sums.rfSums[0] / sums.rfSums[1]),            // rsr
			sums.srg(),                                   // srsr
			sums.cc(),                                    // rsccs
			EDIAm,                                        // ediam
			OPIA,                                         // opia
			static_cast<int>(round(mVF * sums.ngrid)) }); // ngrid

		progress.consumed(1);
	}

	if (addWaters)
	{
		for (const auto &d : atomData)
		{
			const cif::mm::atom &atom = d.atom;

			if (not atom.is_water())
				continue;

			result.emplace_back(ResidueStatistics{ d.asymID, d.seqID, "HOH", d.authSeqID,
				(d.sums.rfSums[0] / d.sums.rfSums[1]),          // rsr
				d.sums.srg(),                                   // srsr
				d.sums.cc(),                                    // rsccs
				d.edia,                                         // ediam
				(d.edia > 0.8 ? 100. : 0.),                     // opia
				static_cast<int>(round(mVF * d.sums.ngrid)) }); // ngrid
		}
	}

	return result;
}

// ResidueStatistics StatsCollector::collect(std::initializer_list<cif::mm::atom> atoms) const
// {
// 	std::vector<cif::mm::atom> v(atoms);
// 	return collect(v);
// }

ResidueStatistics StatsCollector::collect(const std::vector<cif::mm::atom> &atoms) const
{
	std::vector<AtomData> atomData;

	for (auto &ad : mAtomData)
	{
		if (std::ranges::contains(atoms, ad.atom))
			atomData.emplace_back(ad);
	}

	calculate(atomData);

	AtomDataSums sums;
	std::size_t n = 0, m = 0;
	double ediaSum = 0;

	for (auto &ad : atomData)
	{
		++n;

		sums += ad.sums;
		ediaSum += std::pow(ad.edia + 0.1, -2);

		if (ad.edia >= 0.8)
			++m;
	}

	ResidueStatistics result{
		"", 0, "", "",
		(sums.rfSums[0] / sums.rfSums[1]),        // rsr
		sums.srg(),                               // srsr
		sums.cc(),                                // rsccs
		1 / std::sqrt(ediaSum / n) - 0.1,         // ediam
		100. * m / n,                             // opia
		static_cast<int>(round(mVF * sums.ngrid)) // ngrid
	};

	return result;
}

[[nodiscard]] ResidueStatistics StatsCollector::collectSum(const std::string &asymID) const
{
	std::vector<cif::mm::atom> atoms;

	for (auto &ad : mAtomData)
	{
		if (ad.atom.get_label_asym_id() == asymID)
			atoms.emplace_back(ad.atom);
	}

	return collect(atoms);
}

// ResidueStatistics StatsCollector::collect(const std::vector<cif::mm::atom> &atoms) const
// {
// 	AtomDataSums sums;
// 	std::size_t n = 0, m = 0;
// 	double ediaSum = 0;

// 	for (auto &ad : mAtomData)
// 	{
// 		if (not std::ranges::contains(atoms, ad.atom))
// 			continue;

// 		sums += ad.sums;
// 		ediaSum += std::pow(ad.edia + 0.1, -2);
		
// 		++n;
// 		if (ad.edia >= 0.8)
// 			++m;
// 	}

// 	ResidueStatistics result{
// 		"", 0, "", "",
// 		(sums.rfSums[0] / sums.rfSums[1]),        // rsr
// 		sums.srg(),                               // srsr
// 		sums.cc(),                                // rsccs
// 		1 / std::sqrt(ediaSum / n) - 0.1,         // ediam
// 		100. * m / n,                             // opia
// 		static_cast<int>(round(mVF * sums.ngrid)) // ngrid
// 	};

// 	return result;
// }

void StatsCollector::sumDensity(std::vector<AtomData> &atomData,
	GridPtDataMap &gridPointDensity, std::map<std::string, std::vector<double>> &zScoresPerAsym) const
{
	using namespace clipper;

	const Xmap<float> &Fb = mMapMaker.fb();
	const Xmap<float> &Fd = mMapMaker.fd();

	// First step, iterate over atoms, then over grid points covered by this atom
	// collecting per gridpoint statistics

	cif::progress_bar progress(atomData.size(), "Summarising density");

	for (auto &data : atomData)
	{
		auto &atom = data.atom;

		progress.consumed(1);

		if (atom.get_occupancy() == 0)
			continue;

		std::string asymID = data.asymID;
		if (atom.is_water())
			asymID = "0";

		auto radius = data.radius;
		double sumDensity = 0;

		iterateGrid(atom.get_location(), radius, Fb, [&, radius_sq = radius * radius](Xmap_base::Map_reference_coord &iw)
			{
			cif::point p = iw.coord_orth();
			
			double d = distance_squared(p, atom.get_location());

			if (d <= radius_sq)
			{
				double density = data.shape.calculatedDensity(p);
				
				if (std::isnan(density))
					return;

				gridPointDensity[iw.coord()] += density;
				data.points.emplace_back(iw.coord(), density);
				
				sumDensity += density;
				
				zScoresPerAsym[data.asymID].push_back(Fd[iw] / (Fd.multiplicity(iw.coord()) * mRMSDensityFd));
			} });

		data.averageDensity = sumDensity / data.points.size();
	}
}

void StatsCollector::collectSums(std::vector<AtomData> &atomData, const GridPtDataMap &gridPointDensity) const
{
	using namespace clipper;

	const Xmap<float> &Fb = mMapMaker.fb();
	const Xmap<float> &Fd = mMapMaker.fd();

	cif::progress_bar progress_bar(atomData.size(), "Stats calculation");

	// Iterate over the atom data to collect the sums
	for (auto &d : atomData)
	{
		auto rmsi = mRmsScaled.find(d.asymID);

		if (rmsi == mRmsScaled.end())
			continue;

		auto rmsScaledF = rmsi->second;

		for (auto gp : d.points)
		{
			++d.sums.ngrid;

			auto gpd = gridPointDensity.at(gp.p);
			if (gpd == 0)
				continue;

			double e = gp.density / gpd;
			double t = e * mSZ / rmsScaledF.second;

			clipper::Xmap_base::Map_reference_coord ix(Fb, gp.p);

			double fb = Fb[ix];
			double fd = Fd[ix];

			double ed1 = e * (fb - rmsScaledF.first) / rmsScaledF.second + t;
			double ed2 = e * (fd - rmsScaledF.first) / rmsScaledF.second;
			double ed3 = ed1 - ed2;

			d.sums.rfSums[0] += std::abs(ed2);
			d.sums.rfSums[1] += std::abs(ed1 + ed3);

			double w = gp.density / d.averageDensity;
			if (w < 0)
				w = 0;
			if (w > 1)
				w = 1;

			d.sums.rgSums[0] += w * ed2 * ed2;
			d.sums.rgSums[1] += w * ed1 * ed1;

			d.sums.swSums[0] += (w * ed2) * (w * ed2);
			d.sums.swSums[1] += (w * ed1) * (w * ed2);
			d.sums.swSums[2] += (w * ed1) * (w * ed1);

			ed1 -= t;
			ed3 -= t;

			d.sums.ccSums[0] += ed1 * ed3;
			d.sums.ccSums[1] += ed1 * ed1;
			d.sums.ccSums[2] += ed3 * ed3;

			d.sums.edSums[0] += ed1;
			d.sums.edSums[1] += ed3;
		}

		progress_bar.consumed(1);
	}
}

void StatsCollector::calculate(std::vector<AtomData> &atomData) const
{
	collectSums(atomData, mGridPointDensity);
}

// --------------------------------------------------------------------

EDIAStatsCollector::EDIAStatsCollector(const MapMaker<float> &mm, cif::datablock &db, int modelNr, bool electronScattering)
	: StatsCollector(mm, db, modelNr, electronScattering)
{
	// create a atom radius map, for EDIA

	const float kResolutions[] = {
		0.5, 1.0, 1.5, 2.0, 2.5
	};

	// The following numbers were harvested with the application collect-b-factors
	const float kAverageBFactors[] = {
		6.31912, // 0.5
		14.4939, // 1.0
		20.8827, // 1.5
		27.7075, // 2.0
		55.6378  // 2.5
	};
	const int kAverageBFactorCount = sizeof(kAverageBFactors) / sizeof(float);

	int i = static_cast<int>(floor(mResHigh / 0.5)) - 1;
	if (i > kAverageBFactorCount - 1)
		i = kAverageBFactorCount - 1;
	if (i < 0)
		i = 0;

	float ediaBFactor;
	if (i < kAverageBFactorCount - 1)
		ediaBFactor = kAverageBFactors[i] +
		              ((kAverageBFactors[i + 1] - kAverageBFactors[i]) * (mResHigh - kResolutions[i]) / (kResolutions[i + 1] - kResolutions[i]));
	else
		ediaBFactor = kAverageBFactors[i];

	if (cif::VERBOSE > 1)
		std::cerr << "Calculating radii with B Factor " << ediaBFactor << '\n';

	for (auto atom : mAtomData | std::views::transform(&AtomData::atom))
	{
		if (mRadii.count(atom.get_type()))
			continue;

		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering, ediaBFactor);
		mRadii[atom.get_type()] = shape.radius();

		if (cif::VERBOSE > 1)
			std::cerr << "Radius for atom with type " << atom_type_traits(atom.get_type()).symbol() << " is " << mRadii[atom.get_type()] << '\n';
	}
}

void EDIAStatsCollector::calculate(std::vector<AtomData> &atomData) const
{
	StatsCollector::calculate(atomData);

	const clipper::Xmap<float> &Fb = mMapMaker.fb();
	//	Xmap<float>& fd = mMapMaker.fd();

	struct lessAtom
	{
		bool operator()(const cif::mm::atom &a, const cif::mm::atom &b) const { return a.id().compare(b.id()) < 0; }
	};

	using atomSet = std::set<cif::mm::atom, lessAtom>;

	// Calculate EDIA scores

#if __cpp_lib_ranges_to_container >= 202202L and __GNUC__ >= 16
	DistanceMap dm(mAtomData | std::views::transform(&AtomData::atom) | std::ranges::to<std::vector>(),
		cif::crystal(mDb), 3.5f);
#else
	std::vector<cif::mm::atom> dAtoms;
	for (auto &da : mAtomData)
		dAtoms.emplace_back(da.atom);
	DistanceMap dm(std::move(dAtoms), cif::crystal(mDb), 3.5f);
#endif
	BondMap bm{ mDb, std::nullopt, mModelNr };

	cif::progress_bar progress_bar(atomData.size(), "EDIA calculation");

	for (auto &data : atomData)
	{
		auto &atom = data.atom;
		float radius = mRadii.at(atom.get_type());

		PointWeightFunction w(atom.get_location(), radius);

		std::vector<cif::mm::atom> atomsNearBy = dm.near(atom, 3.5f);

		std::vector<PointWeightFunction> wn;
		for (auto a : atomsNearBy)
			wn.emplace_back(a.get_location(), mRadii.at(a.get_type()));

		float ediaSum[2] = {};

		iterateGrid(atom.get_location(), radius, Fb, [&](auto iw)
			{
			cif::point p = iw.coord_orth();
			
			// EDIA calculations
			auto fb = Fb[iw];

			// fix z calculation, thanks to Dmytro Guzenko for spotting the error
			auto z = static_cast<float>((fb - mMeanDensityFb) / mRMSDensityFb);

			if (z < 0)
				z = 0;
			
			if (z > 1.2)
				z = 1.2f;
			
			float wp = w(p);
			
			// And divide the ownership
			
			atomSet S, D, I;
			
			if (wp != 0)
			{
				if (wp < 0)
					D.insert(atom);
				else
				{
					S.insert(atom);
					I.insert(atom);
				}
			}
			
			for (std::size_t i = 0; i < atomsNearBy.size(); ++i)
			{
				float wpi = wn[i](p);
				if (wpi == 0)
					continue;
				
				if (wpi < 0)
					D.insert(atomsNearBy[i]);
				else if (wpi > 0)
				{
					S.insert(atomsNearBy[i]);
					
					if (not bm(atomsNearBy[i], atom))
						I.insert(atomsNearBy[i]);
				}
			}
			
			float o = 0;
			if (wp > 0)
			{
				if (I.size() == 1)
					o = 1;
				else
				{
					auto sumpb = std::ranges::fold_left(I, 0.f, [p](float s, auto &a) {return s + cif::distance(p, a.get_location()); });
					o = 1 - distance(atom.get_location(), p) / sumpb;
				}
			}
			else if (D.count(atom) and S.empty())
			{
				if (D.size() == 1)
					o = 1;
				else
				{
					auto sumpb = std::ranges::fold_left(D, 0.f, [p](float s, auto &a) {return s + cif::distance(p, a.get_location()); });
					o = 1 - distance(atom.get_location(), p) / sumpb;
				}
			}

			ediaSum[0] += z * wp * o;
			if (wp > 0)
				ediaSum[1] += wp; });

		data.edia = ediaSum[0] / ediaSum[1];
		if (data.edia < 0)
			data.edia = 0;

		progress_bar.consumed(1);
	}
}

BondMap EDIAStatsCollector::createBondMap(std::vector<AtomData> &atomData) const
{
	std::vector<cif::point> pts;
	for (auto a : atomData)
		pts.emplace_back(a.atom.get_location());

	auto [center, radius] = cif::smallest_sphere_around_points(pts);

	return { mDb, std::make_tuple(center, radius + 3.5f) };
}

} // namespace pdb_redo
