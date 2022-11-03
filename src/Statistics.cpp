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

#include <fstream>
#include <numeric>

#include <cif++.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/BondMap.hpp"
#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/Statistics.hpp"

// --------------------------------------------------------------------

namespace pdb_redo
{

using clipper::Coord_grid;
using clipper::Coord_orth;
using clipper::Xmap;

using cif::atom_type_traits;

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
		4.5921953931549871457e+4, 6.7265770927008700853e+4, 3.3430575583588128105e+4, 2.5090809287301226727e+3},
				 B[8] = {0, 4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3, 2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4, 5.2264952788528545610e+3};

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
				 D[8] = {0, 2.05319162663775882187e0, 1.67638483018380384940e0, 6.89767334985100004550e-1, 1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9};

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
				 F[8] = {0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2, 7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15};

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
	auto c = std::sqrt(2.0 / kPI);
	auto b = std::abs(a);

	double result = 0;
	if (b > 3 / c)
	{
		auto x = std::abs(std::pow(b, 1 / 3.0) - 2 * std::pow(kPI / b, 2));
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
		m_P[0] = P{-1.0f, 0, 1.0f, 1.0822f};
		m_P[1] = P{5.1177f, 1.29366f, -0.4f, 1.4043f};
		m_P[2] = P{-0.9507f, 2, 0, 2};
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

struct AtomGridData
{
	AtomGridData(const Coord_grid &gp, double density)
		: p(gp)
		, density(density)
	{
	}

	Coord_grid p;
	double density;
};

struct AtomDataSums
{
	size_t ngrid = 0;
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

	double cc() const
	{
		double s = (ccSums[1] - (edSums[0] * edSums[0]) / ngrid) * (ccSums[2] - (edSums[1] * edSums[1]) / ngrid);
		return (ccSums[0] - edSums[0] * edSums[1] / ngrid) / std::sqrt(s);
	}

	double srg() const
	{
		double rgsq = rgSums[0] / rgSums[1];
		double rg = std::sqrt(rgsq);

		return std::sqrt(swSums[0] - rgsq * swSums[1] + 0.5 * rgsq * rgsq * swSums[2]) / (rg * rgSums[1]);
	}
};

struct AtomData
{
	AtomData(cif::mm::atom atom, float radius)
		: atom(atom)
		, asymID(atom.get_label_asym_id())
		, seqID(atom.get_label_seq_id())
		, auth_seq_id(atom.get_auth_seq_id())
		, radius(radius)
		, occupancy(atom.get_occupancy())
	{
	}

	cif::mm::atom atom;
	std::string asymID;
	int seqID;
	std::string auth_seq_id; // required for waters
	float radius;
	float occupancy;
	std::vector<AtomGridData> points;
	double averageDensity = 0;
	double edia = 0;
	AtomDataSums sums;
};

// --------------------------------------------------------------------

std::tuple<float, float> CalculateMapStatistics(const Xmap<float> &f)
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

	float meanDensity = static_cast<float>(sum / count);
	float rmsDensity = static_cast<float>(std::sqrt((sum2 / count) - (meanDensity * meanDensity)));

	return std::make_tuple(meanDensity, rmsDensity);
}

// --------------------------------------------------------------------

class BoundingBox
{
  public:

	template <class List>
	BoundingBox(const cif::mm::structure &structure, List atoms, float margin)
	{
		mXMin = mYMin = mZMin = std::numeric_limits<float>::max();
		mXMax = mYMax = mZMax = std::numeric_limits<float>::min();

		for (auto &atom : atoms)
		{
			auto l = atom.get_location();
			if (mXMin > l.m_x)
				mXMin = l.m_x;
			if (mXMax < l.m_x)
				mXMax = l.m_x;
			if (mYMin > l.m_y)
				mYMin = l.m_y;
			if (mYMax < l.m_y)
				mYMax = l.m_y;
			if (mZMin > l.m_z)
				mZMin = l.m_z;
			if (mZMax < l.m_z)
				mZMax = l.m_z;
		}

		mXMin -= margin;
		mXMax += margin;
		mYMin -= margin;
		mYMax += margin;
		mZMin -= margin;
		mZMax += margin;
	}

	bool contains(const cif::point &p) const
	{
		return p.m_x >= mXMin and p.m_x <= mXMax and p.m_y >= mYMin and p.m_y <= mYMax and p.m_z >= mZMin and p.m_z <= mZMax;
	}

  private:
	float mXMin, mXMax, mYMin, mYMax, mZMin, mZMax;
};
// --------------------------------------------------------------------

StatsCollector::StatsCollector(const MapMaker<float> &mm, cif::mm::structure &structure, bool electronScattering)
	: mStructure(structure)
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
	for (auto &asym_id : mStructure.get_datablock()["struct_asym"].rows<std::string>("id"))
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
	//	const double C = std::sqrt(2.0 / kPI);

	for (auto &a : mStructure.atoms())
	{
		auto t = a.get_type();
		if (t <= cif::atom_type::He)
			continue;

		float w = a.get_occupancy() * t;

		if (w <= 0)
			continue;

		mSZ += w;

		//		float bIso = Util::u2b(a.uIso());
		//		if (bIso < 4)
		//			bIso = 4;
		//		float x = std::sqrt(bIso) / mResHigh;
		//		x = w * (2 * anorm(x) - 1 - C * x * std::exp(-0.5 * std::pow(x, 2))) / std::pow(x, 3);
		//
		//		so += x;
	}

	//	auto bo = mSZ;

	mSZ = mSZ * mSpacegroup.num_symops() / mVC;
	//	mMeanBIso = std::pow(mResHigh * errsol(bo / so), 2);

	// Calculate overall rms data
	std::vector<AtomData> atomData;

	for (auto atom : mStructure.atoms())
	{
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);

		float radius = shape.radius();

		if (cif::VERBOSE > 2)
			std::cerr << (atomData.size() + 1) << '\t'
					  << atom_type_traits(atom.get_type()).symbol() << '\t'
					  << radius << std::endl;

		atomData.emplace_back(atom, radius);
	}

	GridPtDataMap gridPointDensity;
	std::map<std::string, std::vector<double>> zScoresPerAsym;
	sumDensity(atomData, gridPointDensity, zScoresPerAsym);

	// Now that we have the density data, we can calculate the correction/rescale factors
	for (auto zsc : zScoresPerAsym)
	{
		// collect array of z-scores
		std::vector<double> &zdca0 = zsc.second;

		auto &z = zdca0;
		auto vf = mVF;

		sort(z.begin(), z.end());

		double qa = 0, qb = 1;

		size_t nd = z.size();
		size_t n = static_cast<size_t>(round(vf * nd));

		if (n > 100)
		{
			size_t i1 = static_cast<size_t>((n + 1) * anorm(-1.5)) + 1;
			size_t i2 = static_cast<size_t>((n + 1) * anorm(1.5));

			size_t ns = i2 - i1 + 1;

			double vr = (nd - 1) / (n - 1.0);
			double sw = 0, swx = 0, swxs = 0, swy = 0, swxy = 0, swys = 0;

			for (auto i = i1; i <= i2; ++i)
			{
				double qx = phinvs(static_cast<double>(i) / (n + 1));
				double x = vr * i;
				size_t j = static_cast<size_t>(x);
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
				std::cerr << std::endl
						  << "Intercept & gradient before LS: " << qa << " (" << std::sqrt(swys * swxs) << ") " << qb << " (" << std::sqrt(swys * sw) << ')' << std::endl;
			}

			qb += 1.0;

			if (cif::VERBOSE > 1)
			{
				std::cerr << std::endl
						  << "Rescale SD(delta-rho) using Q-Q plot for asym " << zsc.first << ':' << std::endl
						  << std::string(54, '=') << std::endl
						  << "Input & updated SD(delta-rho): " << mRMSDensityFd << " ; " << qb * mRMSDensityFd << std::endl
						  << std::endl;
			}
		}

		mRmsScaled[zsc.first] = std::make_pair(qa * mRMSDensityFd, qb * mRMSDensityFd);
	}
}

std::vector<ResidueStatistics> StatsCollector::collect() const
{
	std::vector<std::tuple<std::string, int, std::string>> residues;
	std::vector<cif::mm::atom> atoms;

	for (auto atom : mStructure.atoms())
	{
		if (atom.is_water())
			continue;

		auto k = std::make_tuple(atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());

		if (residues.empty() or residues.back() != k)
		{
			residues.emplace_back(move(k));
			atoms.emplace_back(std::move(atom));
		}
	}

	BoundingBox bbox(mStructure, atoms, 5.0f);
	return collect(residues, bbox, true);
}

std::vector<ResidueStatistics> StatsCollector::collect(const std::string &asymID, int resFirst, int resLast, bool authNameSpace) const
{
	residue_list residues;
	std::vector<cif::mm::atom> atoms;

	// for (auto atom : mStructure.atoms())
	// {
	// 	if (atom.is_water())
	// 		continue;

	// 	if (authNameSpace)
	// 	{
	// 		int auth_seq_id = stoi(atom.get_auth_seq_id());

	// 		if (atom.get_auth_asym_id() != asymID or auth_seq_id < resFirst or auth_seq_id > resLast)
	// 			continue;
	// 	}
	// 	else
	// 	{
	// 		if (atom.get_label_asym_id() != asymID or atom.get_label_seq_id() < resFirst or atom.get_label_seq_id() > resLast)
	// 			continue;
	// 	}

	// 	auto k = std::make_tuple(atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());

	// 	if (residues.empty() or residues.back() != k)
	// 	{
	// 		residues.emplace_back(move(k));
	// 		atoms.emplace_back(std::move(atom));
	// 	}
	// }

	for (auto atom : mStructure.atoms())
	{
		if (atom.is_water())
			continue;

		if (authNameSpace)
		{
			int auth_seq_id = stoi(atom.get_auth_seq_id());

			if (atom.get_auth_asym_id() != asymID or auth_seq_id < resFirst or auth_seq_id > resLast)
				continue;
		}
		else
		{
			if (atom.get_label_asym_id() != asymID or atom.get_label_seq_id() < resFirst or atom.get_label_seq_id() > resLast)
				continue;
		}

		auto k = std::make_tuple(atom.get_label_asym_id(), atom.get_label_seq_id(), atom.get_auth_seq_id());

		if (residues.empty() or residues.back() != k)
			residues.emplace_back(move(k));
	}

	for (const auto &[asymID, seqID, authSeqID] : residues)
	{
		auto &res = mStructure.get_residue(asymID, seqID, authSeqID);

		for (auto atom : res.unique_atoms())
			atoms.push_back(atom);
	}

	BoundingBox bbox(mStructure, atoms, 5.0f);
	return collect(residues, bbox, false);
}

std::vector<ResidueStatistics> StatsCollector::collect(const residue_list &residues, BoundingBox &bbox, bool addWaters) const
{
	std::vector<AtomData> atomData;

	//	BoundingBox bb(mStructure, residues, 5.0f);

	for (auto atom : mStructure.atoms())
	{
		if (atom.is_water())
		{
			if (not addWaters)
				continue;
		}
		else if (not bbox.contains(atom.get_location()))
			continue;

		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);

		float radius = shape.radius();

		if (cif::VERBOSE > 2)
			std::cerr << (atomData.size() + 1) << '\t'
					  << atom_type_traits(atom.get_type()).symbol() << '\t'
					  << radius << std::endl;

		atomData.emplace_back(atom, radius);
	}

	calculate(atomData);

	std::set<std::string> missing;
	std::vector<ResidueStatistics> result;

	// And now collect the per residue information
	for (const auto &[asymID, seqID, auth_seq_id] : residues)
	{
		// TODO: Need to do something with hetero residues (alternate compound types)
		auto &res = mStructure.get_residue(asymID, seqID, auth_seq_id);
		auto compID = res.get_compound_id();

		AtomDataSums sums;

		std::vector<const AtomData *> resAtomData;
		for (const auto &d : atomData)
		{
			if (d.asymID == asymID and d.seqID == seqID and d.auth_seq_id == auth_seq_id)
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

					resAtomData.erase(
						std::remove_if(resAtomData.begin(), resAtomData.end(), [id = compAtom](const AtomData *d)
							{ return d->atom.get_label_atom_id() == id; }),
						resAtomData.end());
				}
			}
		}
		catch (const BondMapException &ex)
		{
			std::cerr << "Missing information for compound '" << compID << '\'' << std::endl
					  << ex.what() << std::endl;

			missing.insert(compID);
		}

		// atoms that were present but not part of the Compound
		for (auto d : resAtomData)
			sums += d->sums;

		// EDIA

		std::set<std::string> alts;

		if (not missing.count(compID))
			alts = res.get_alternate_ids();

		if (alts.empty())
			alts.insert("");

		double EDIAm = 0, OPIA = 0, OCC = 0;

		for (auto alt : alts)
		{
			double ediaSum = 0, occSum = 0;
			size_t n = 0, m = 0;

			for (const auto &d : atomData)
			{
				if (d.asymID != asymID or d.seqID != seqID or d.auth_seq_id != auth_seq_id)
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

		result.emplace_back(ResidueStatistics{asymID, seqID, compID,
			auth_seq_id,
			(sums.rfSums[0] / sums.rfSums[1]),           // rsr
			sums.srg(),                                  // srsr
			sums.cc(),                                   // rsccs
			EDIAm,                                       // ediam
			OPIA,                                        // opia
			static_cast<int>(round(mVF * sums.ngrid))}); // ngrid
	}

	if (addWaters)
	{
		for (const auto &d : atomData)
		{
			const cif::mm::atom &atom = d.atom;

			if (not atom.is_water())
				continue;

			result.emplace_back(ResidueStatistics{d.asymID, d.seqID, "HOH",
				atom.get_auth_seq_id(),
				(d.sums.rfSums[0] / d.sums.rfSums[1]),         // rsr
				d.sums.srg(),                                  // srsr
				d.sums.cc(),                                   // rsccs
				d.edia,                                        // ediam
				(d.edia > 0.8 ? 100. : 0.),                    // opia
				static_cast<int>(round(mVF * d.sums.ngrid))}); // ngrid
		}
	}

	return result;
}

ResidueStatistics StatsCollector::collect(std::initializer_list<const cif::mm::residue *> residues) const
{
	std::vector<cif::mm::atom> atoms;
	for (auto &r : residues)
		for (auto a : r->atoms())
			atoms.push_back(a);

	return collect(atoms);
}

ResidueStatistics StatsCollector::collect(std::initializer_list<cif::mm::atom> atoms) const
{
	std::vector<cif::mm::atom> v(atoms);
	return collect(v);
}

ResidueStatistics StatsCollector::collect(const std::vector<cif::mm::atom> &atoms) const
{
	std::vector<AtomData> atomData;

	BoundingBox bb(mStructure, atoms, 4.f);

	for (auto atom : mStructure.atoms())
	{
		if (not bb.contains(atom.get_location()))
			continue;

		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);

		float radius = shape.radius();

		if (cif::VERBOSE > 2)
			std::cerr << (atomData.size() + 1) << '\t'
					  << atom_type_traits(atom.get_type()).symbol() << '\t'
					  << radius << std::endl;

		atomData.emplace_back(atom, radius);
	}

	calculate(atomData);

	AtomDataSums sums;
	size_t n = 0, m = 0;
	double ediaSum = 0;

	for (auto &atom : atoms)
	{
		++n;

		auto ci = find_if(atomData.begin(), atomData.end(),
			[=](auto &d)
			{ return d.asymID == atom.get_label_asym_id() and d.seqID == atom.get_label_seq_id() and d.atom.get_label_atom_id() == atom.get_label_atom_id(); });

		if (ci == atomData.end())
			continue;

		sums += ci->sums;
		ediaSum += std::pow(ci->edia + 0.1, -2);

		if (ci->edia >= 0.8)
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

void StatsCollector::sumDensity(std::vector<AtomData> &atomData,
	GridPtDataMap &gridPointDensity, std::map<std::string, std::vector<double>> &zScoresPerAsym) const
{
	using namespace clipper;

	const Xmap<float> &Fb = mMapMaker.fb();
	const Xmap<float> &Fd = mMapMaker.fd();

	// First step, iterate over atoms, then over grid points covered by this atom
	// collecting per gridpoint statistics

	for (auto &data : atomData)
	{
		auto &atom = data.atom;

		if (atom.get_occupancy() == 0)
			continue;
		
		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering);

		std::string asymID = data.asymID;
		if (atom.is_water())
			asymID = "0";
		
		auto radius = data.radius;
		double sumDensity = 0;

		iterateGrid(toClipper(atom.get_location()), radius, Fb, [&](Xmap_base::Map_reference_coord &iw)
			{
			auto p = toPoint(iw.coord_orth());
			
			double d = distance(p, atom.get_location());

			if (d <= radius)
			{
				double density = shape.calculatedDensity(p);
				
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

void StatsCollector::collectSums(std::vector<AtomData> &atomData, GridPtDataMap &gridPointDensity) const
{
	using namespace clipper;

	const Xmap<float> &Fb = mMapMaker.fb();
	const Xmap<float> &Fd = mMapMaker.fd();

	cif::Progress progress(atomData.size(), "Stats calculation");

	// Iterate over the atom data to collect the sums
	for (auto &d : atomData)
	{
		auto rmsScaledF = mRmsScaled.at(d.asymID);

		for (auto gp : d.points)
		{
			++d.sums.ngrid;

			auto gpd = gridPointDensity[gp.p];
			if (gpd == 0)
				continue;

			double e = gp.density / gpd;
			double t = e * mSZ / rmsScaledF.second;

			Xmap_base::Map_reference_coord ix(Fb, gp.p);

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

		progress.consumed(1);
	}
}

void StatsCollector::calculate(std::vector<AtomData> &atomData) const
{
	GridPtDataMap gridPointDensity;
	std::map<std::string, std::vector<double>> zScoresPerAsym;

	sumDensity(atomData, gridPointDensity, zScoresPerAsym);
	collectSums(atomData, gridPointDensity);
}

// --------------------------------------------------------------------

EDIAStatsCollector::EDIAStatsCollector(MapMaker<float> &mm,
	cif::mm::structure &structure, bool electronScattering, const BondMap &bondMap)
	: StatsCollector(mm, structure, electronScattering)
	, mDistanceMap(structure, mm.spacegroup(), mm.cell(), 3.5f)
	, mBondMap(bondMap)
{
	// create a atom radius map, for EDIA

	const double kResolutions[] =
		{
			0.5, 1.0, 1.5, 2.0, 2.5};

	// The following numbers were harvested with the application collect-b-factors
	const double kAverageBFactors[] =
		{
			6.31912, // 0.5
			14.4939, // 1.0
			20.8827, // 1.5
			27.7075, // 2.0
			55.6378  // 2.5
		};
	const int kAverageBFactorCount = sizeof(kAverageBFactors) / sizeof(double);

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

	if (cif::VERBOSE)
		std::cerr << "Calculating radii with B Factor " << ediaBFactor << std::endl;

	for (auto atom : mStructure.atoms())
	{
		if (mRadii.count(atom.get_type()))
			continue;

		AtomShape shape(atom, mResHigh, mResLow, mElectronScattering, ediaBFactor);
		mRadii[atom.get_type()] = shape.radius();

		if (cif::VERBOSE)
			std::cerr << "Radius for atom with type " << atom_type_traits(atom.get_type()).symbol() << " is " << mRadii[atom.get_type()] << std::endl;
	}
}

void EDIAStatsCollector::calculate(std::vector<AtomData> &atomData) const
{
	StatsCollector::calculate(atomData);

	const Xmap<float> &Fb = mMapMaker.fb();
	//	Xmap<float>& fd = mMapMaker.fd();

	struct lessAtom
	{
		bool operator()(const cif::mm::atom &a, const cif::mm::atom &b) const { return a.id().compare(b.id()) < 0; }
	};

	typedef std::set<cif::mm::atom, lessAtom> atomSet;

	// Calculate EDIA scores

	cif::Progress progress(atomData.size(), "EDIA calculation");

	for (auto &data : atomData)
	{
		auto &atom = data.atom;
		float radius = mRadii.at(atom.get_type());

		//		if (cif::VERBOSE > 2)
		//			std::cerr << (atomData.size() + 1) << '\t'
		//				 << atom_type_traits(atom.get_type()).symbol() << '\t'
		//				 << radius << std::endl;
		//
		PointWeightFunction w(atom.get_location(), radius);

		std::vector<cif::mm::atom> atomsNearBy = mDistanceMap.near(atom, 3.5f);

		std::vector<PointWeightFunction> wn;
		for (auto a : atomsNearBy)
			wn.emplace_back(a.get_location(), mRadii.at(a.get_type()));

		float ediaSum[2] = {};

		iterateGrid(toClipper(atom.get_location()), radius, Fb, [&](auto iw)
		{
			cif::point p = toPoint(iw.coord_orth());
			
			// EDIA calculations
			auto fb = Fb[iw];

			float z = 0;
			if (fb >= mMeanDensityFb + mRMSDensityFb)
				z = static_cast<float>((fb - mMeanDensityFb) / mRMSDensityFb);
			
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
			
			for (size_t i = 0; i < atomsNearBy.size(); ++i)
			{
				float wpi = wn[i](p);
				if (wpi == 0)
					continue;
				
				if (wpi < 0)
					D.insert(atomsNearBy[i]);
				else if (wpi > 0)
				{
					S.insert(atomsNearBy[i]);
					
					if (not mBondMap(atomsNearBy[i], atom))
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
					float sumpb = accumulate(I.begin(), I.end(), 0.f,
						[p](float s, const cif::mm::atom& b) -> float
						{
							return s + distance(p, b.get_location());
						});

					o = 1 - distance(atom.get_location(), p) / sumpb;
				}
			}
			else if (D.count(atom) and S.empty())
			{
				if (D.size() == 1)
					o = 1;
				else
				{
					float sumpb = accumulate(D.begin(), D.end(), 0.f,
						[p](float s, const cif::mm::atom& b) -> float
						{
							return s + distance(p, b.get_location());
						});

					o = 1 - distance(atom.get_location(), p) / sumpb;
				}
			}

			ediaSum[0] += z * wp * o;
			if (wp > 0)
				ediaSum[1] += wp;
		});

		data.edia = ediaSum[0] / ediaSum[1];
		if (data.edia < 0)
			data.edia = 0;

		progress.consumed(1);
	}
}

} // namespace pdb_redo
