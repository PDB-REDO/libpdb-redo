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

#include <mutex>
#include <thread>

#include <cif++.hpp>

#include "pdb-redo/AtomShape.hpp"
#include "pdb-redo/ClipperWrapper.hpp"

#if __has_include(<newuoa.h>)

#include <newuoa.h>
#define HAVE_NEWUOA 1

#elif __has_include(<dlib/global_optimization.h>)

#include <dlib/global_optimization.h>

#define HAVE_DLIB 1

#else
#error "Should have either newuoa or dlib"
#endif

namespace pdb_redo
{

using cif::atom_type;
using cif::atom_type_traits;
using cif::point;

using cif::kPI;

// --------------------------------------------------------------------

// sine integration function based on the description in
//	GalSim: The modular galaxy image simulation toolkit
//	B.T.P. Rowe et. al.
//	DOI: 10.1016/j.ascom.2015.02.002

double sineIntegration(double x)
{
	struct P
	{
		double n, d;
	};
	double result = 0;

	if (x <= 4)
	{
		const P kP[] = {
			{
				1,
				1,
			},
			{ -4.54393409816329991e-2, 1.01162145739225565e-2 },
			{ 1.15457225751016682e-3, 4.99175116169755106e-5 },
			{ -1.41018536821330254e-5, 1.55654986308745614e-7 },
			{ 9.43280809438713025e-8, 3.28067571055789734e-10 },
			{ -3.53201978997168357e-10, 4.5049097575386581e-13 },
			{ 7.08240282274875911e-13, 3.21107051193712168e-16 },
			{ -6.05338212010422477e-16, 0 }
		};

		double xs = x * x;
		double xi = 1;
		double sn = 0, sd = 0;

		for (auto p : kP)
		{
			sn += xi * p.n;
			sd += xi * p.d;
			xi *= xs;
		}

		result = x * (sn / sd);
	}
	else
	{
		const P
			kF[] = {
				{
					1,
					1,
				},
				{
					7.44437068161936700618e2,
					7.46437068161927678031e2,
				},
				{
					1.96396372895146869801e5,
					1.97865247031583951450e5,
				},
				{
					2.37750310125431834034e7,
					2.41535670165126845144e7,
				},
				{
					1.43073403821274636888e9,
					1.47478952192985464958e9,
				},
				{
					4.33736238870432522765e10,
					4.58595115847765779830e10,
				},
				{
					6.40533830574022022911e11,
					7.08501308149515401563e11,
				},
				{
					4.20968180571076940208e12,
					5.06084464593475076774e12,
				},
				{
					1.00795182980368574617e13,
					1.43468549171581016479e13,
				},
				{
					4.94816688199951963482e12,
					1.11535493509914254097e13,
				},
				{ -4.94701168645415959931e11, 0 }
			},
			kG[] = { {
						 1,
						 1,
					 },
				{
					8.1359520115168615e2,
					8.19595201151451564e2,
				},
				{
					2.35239181626478200e5,
					2.40036752835578777e5,
				},
				{
					3.12557570795778731e7,
					3.26026661647090822e7,
				},
				{
					2.06297595146763354e9,
					2.23355543278099360e9,
				},
				{
					6.83052205423625007e10,
					7.87465017341829930e10,
				},
				{
					1.09049528450362786e12,
					1.39866710696414565e12,
				},
				{
					7.57664583257834349e12,
					1.17164723371736605e13,
				},
				{
					1.81004487464664575e13,
					4.01839087307656620e13,
				},
				{
					6.43291613143049485e12,
					3.99653257887490811e13,
				},
				{ -1.36517137670871689e12, 0 } };

		double xs = std::pow(x, -2);
		double xi = 1;

		double sn = 0, sd = 0;
		for (auto p : kF)
		{
			sn += xi * p.n;
			sd += xi * p.d;
			xi *= xs;
		}

		double fx = (sn / sd) / x;

		sn = 0;
		sd = 0;
		xi = 1;
		for (auto p : kG)
		{
			sn += xi * p.n;
			sd += xi * p.d;
			xi *= xs;
		}

		double gx = (sn / sd) / (x * x);

		result = cif::kPI / 2 - fx * std::cos(x) - gx * std::sin(x);
	}

	return result;
}

// --------------------------------------------------------------------
// Internal class to cache some common data

class DensityIntegration
{
  public:
	DensityIntegration(float resolutionLow, float resolutionHigh);

	static DensityIntegration &instance(float resolutionLow, float resolutionHigh);

	double integrateRadius(float perc, float occupancy, double yi, const std::vector<double> &fst) const;
	double integrateDensity(double r, int ks, const std::vector<double> &fst) const;

	float a() const { return mA; }
	float b() const { return mB; }
	const std::vector<double> &st() const { return mST; }
	const std::vector<double> &sts() const { return mSTS; }
	const std::vector<double> &wa() const { return mWA; }

  private:
	float mA, mB;
	int mM;

	// Gauss-Legendre quadrature weights and abscissae
	std::vector<double> mWA, mST, mSTS;
	static std::list<DensityIntegration> sInstances;
};

std::list<DensityIntegration> DensityIntegration::sInstances;

DensityIntegration &DensityIntegration::instance(float resolutionLow, float resolutionHigh)
{
	static std::mutex m;
	std::lock_guard lock(m);

	float a = 0.5f / resolutionLow, b = 0.5f / resolutionHigh;

	auto i = find_if(sInstances.begin(), sInstances.end(), [=](const DensityIntegration &di)
		{ return di.mA == a and di.mB == b; });

	if (i == sInstances.end())
	{
		sInstances.emplace_back(resolutionLow, resolutionHigh);
		i = prev(sInstances.end());
	}

	return *i;
}

DensityIntegration::DensityIntegration(float resolutionLow, float resolutionHigh)
{
	mA = 0.5f / resolutionLow;
	mB = 0.5f / resolutionHigh;

	mM = static_cast<int>(12.0 * std::sqrt(1.2 * mB) + 1);
	if (mM < 3)
		mM = 3;

	int N = 2 * mM;
	//	int J = N;
	double xr = mB - mA;
	double xh = .5 * xr;
	double xm = 0.5 * (mB + mA);

	mWA = std::vector<double>(N, 0);
	mST = std::vector<double>(N, 0);
	mSTS = std::vector<double>(N, 0);

	for (int i = 1, j = N; i <= mM; ++i, --j)
	{
		double z, zo, dp;

		z = std::cos(cif::kPI * (i - 0.25) / (N + 0.5));

		do
		{
			double p1 = 1;
			double p2 = 0;

			for (int k = 1; k <= N; ++k)
			{
				double p3 = p2;
				p2 = p1;
				p1 = ((2 * k - 1) * z * p2 - (k - 1) * p3) / k;
			}

			dp = N * (z * p1 - p2) / (z * z - 1);
			zo = z;
			z = z - p1 / dp;
		} while (std::abs(z - zo) > 3e-14);

		mWA[i - 1] = xr / ((1 - z * z) * dp * dp);
		mWA[j - 1] = mWA[i - 1];
		mST[i - 1] = xm - xh * z;
		mST[j - 1] = xm + xh * z;
	}

	transform(mST.begin(), mST.end(), mSTS.begin(), [](double s)
		{ return s * s; });
}

// --------------------------------------------------------------------
// Calculate Radius integral over r of calculated density
// code inspired by radint.f in edstats

double DensityIntegration::integrateDensity(double r, int ks, const std::vector<double> &fst) const
{
	double y = 0;

	double rt = r;
	if (rt < 0)
		rt = 0;

	if (rt > 1e-10)
	{
		double t = 4 * cif::kPI * rt;
		y = 0;

		for (size_t i = 0; i < mST.size(); ++i)
			y += fst[i] * sineIntegration(t * mST[i]);

		if (r < 0)
			y = y - ks * r;
	}

	return ks * y;
}

#if HAVE_NEWUOA
template <class F>
NewuoaClosure make_closure(F &function)
{
	struct Wrap
	{
		static double call(void *data, long n, const double *values)
		{
			return reinterpret_cast<F *>(data)->operator()(n, values);
		}
	};
	return NewuoaClosure{ &function, &Wrap::call };
}
#endif

double DensityIntegration::integrateRadius(float perc, float occupancy, double yi, const std::vector<double> &fst) const
{
	double yt = perc * 0.25 * cif::kPI * occupancy * yi;

#if HAVE_NEWUOA
	double initialValue = 0.25;

	// code from newuoa-example
	const long variables_count = 1;
	const long number_of_interpolation_conditions = (variables_count + 1) * (variables_count + 2) / 2;
	double variables_values[] = { initialValue };
	const double initial_trust_region_radius = 1e-3;
	const double final_trust_region_radius = 1e3;
	const long max_function_calls_count = 100;
	const size_t working_space_size = NEWUOA_WORKING_SPACE_SIZE(variables_count,
		number_of_interpolation_conditions);
	double working_space[working_space_size];

	auto function = [&](long n, const double *x)
	{
		assert(n == 1);
		return this->integrateDensity(x[0], -1, fst);
	};
	auto closure = make_closure(function);

	double result = newuoa_closure(
		&closure,
		variables_count,
		number_of_interpolation_conditions,
		variables_values,
		initial_trust_region_radius,
		final_trust_region_radius,
		max_function_calls_count,
		working_space);

	double y1 = 0;
	double y2 = -result;
	double x1 = 0;
	double x2 = variables_values[0];
#else
	auto function = [&](const double x)
	{
		return this->integrateDensity(x, -1, fst);
	};

	auto r = dlib::find_min_global(function, { 1e-3 }, { 1e3 }, { false }, dlib::max_function_calls(10));

	double result = r.x(0);

	double y1 = 0;
	double y2 = -r.y;
	double x1 = 0;
	double x2 = result;
#endif

	const double kRE = 5e-5;

	if (y2 > yt)
	{
		for (int it = 0; it < 100; ++it)
		{
			double x = 0.5 * (x1 + x2);
			double y = integrateDensity(x, 1, fst);

			if (std::abs(y - yt) < kRE * std::abs(yt))
			{
				result = x;
				break;
			}

			if ((y1 < yt and y < yt) or (y1 > yt and y > yt))
			{
				x1 = x;
				y1 = y;
			}
			else
			{
				x2 = x;
				y2 = y;
			}
		}
	}
	else
		result = x2;

	return result;
}

// --------------------------------------------------------------------

struct AtomShapeImpl
{
	virtual ~AtomShapeImpl() = default;

	AtomShapeImpl(point location, atom_type symbol, int charge, float uIso, float occupancy, float resHigh, float resLow, bool electronScattering)
		: mSymbol(symbol)
		, mCharge(charge)
		, mUIso(uIso)
		, mOccupancy(occupancy)
		, mResHigh(resHigh)
		, mResLow(resLow)
		, mElectronScattering(electronScattering)
		, mLocation(location)
		, mIntegrator(DensityIntegration::instance(resLow, resHigh))
	{
		auto st = mIntegrator.st();
		auto sts = mIntegrator.sts();
		auto wa = mIntegrator.wa();

		mYi = 0;
		mFst = std::vector<double>(st.size(), 0);

		auto &D =
			mElectronScattering ? atom_type_traits(symbol).elsf() : atom_type_traits(symbol).wksf(charge);
		auto bIso = clipper::Util::u2b(uIso);

		if (bIso == 0 and cif::VERBOSE >= 0)
			std::cerr << "Zero b-factor?\n";

		float as = mIntegrator.a() * mIntegrator.a();
		float bs = mIntegrator.b() * mIntegrator.b();

		for (int i = 0; i < 6; ++i)
		{
			double bi = D.b[i] + bIso;
			mYi += D.a[i] * (std::exp(-bi * as) - std::exp(-bi * bs)) / bi;
		}

		for (size_t i = 0; i < st.size(); ++i)
		{
			double t = 0;
			for (int j = 0; j < 6; ++j)
			{
				double bj = D.b[j] + bIso;
				t += D.a[j] * std::exp(-bj * sts[i]);
			}

			mFst[i] = occupancy * wa[i] * t * st[i];
		}

		for (size_t i = 0; i < 6; ++i)
		{
			mBW[i] = static_cast<float>(-4 * kPI * kPI / (D.b[i] + bIso));
			mAW[i] = static_cast<float>(D.a[i] * std::pow(-mBW[i] / kPI, 1.5));
		}
	}

	atom_type mSymbol;
	int mCharge;
	float mUIso, mOccupancy;
	float mResHigh, mResLow;
	bool mElectronScattering;
	point mLocation;

	const DensityIntegration &mIntegrator;
	double mYi;
	std::vector<double> mFst;
	float mAW[6], mBW[6];

	virtual float integratedRadius(float perc) const
	{
		float result = static_cast<float>(mIntegrator.integrateRadius(perc, mOccupancy, mYi, mFst));

		assert(not std::isnan(result));

		return result;
	}

	virtual float calculatedDensity(float r) const
	{
		float rsq = r * r;
		return mOccupancy *
		       (mAW[0] * std::exp(mBW[0] * rsq) + mAW[1] * std::exp(mBW[1] * rsq) +
				   mAW[2] * std::exp(mBW[2] * rsq) + mAW[3] * std::exp(mBW[3] * rsq) +
				   mAW[4] * std::exp(mBW[4] * rsq) + mAW[5] * std::exp(mBW[5] * rsq));
	}

	virtual float calculatedDensity(point p) const
	{
		return calculatedDensity(distance(mLocation, p));
	}
};

struct AtomShapeAnisoImpl : public AtomShapeImpl
{
	AtomShapeAnisoImpl(point location, atom_type symbol, int charge, clipper::U_aniso_orth &anisou, float occupancy, float resHigh, float resLow, bool electronScattering)
		: AtomShapeImpl(location, symbol, charge, static_cast<float>(anisou.u_iso()), occupancy, resHigh, resLow, electronScattering)
		, mAnisoU(anisou)
	{
		auto &D =
			mElectronScattering ? atom_type_traits(symbol).elsf() : atom_type_traits(symbol).wksf(charge);

		const float fourpi2 = static_cast<float>(4 * kPI * kPI);
		const float pi3 = static_cast<float>(kPI * kPI * kPI);

		for (int i = 0; i < 6; ++i)
		{
			mAnisoInv[i] = clipper::Mat33sym<>(
				-2 * mAnisoU.mat00() - D.b[i] / fourpi2,
				-2 * mAnisoU.mat11() - D.b[i] / fourpi2,
				-2 * mAnisoU.mat22() - D.b[i] / fourpi2,
				-2 * mAnisoU.mat01(),
				-2 * mAnisoU.mat02(),
				-2 * mAnisoU.mat12())
			                   .inverse();

			double det = -mAnisoInv[i].det();

			mAW[i] = static_cast<float>(D.a[i] * std::sqrt(det / pi3));
			mBW[i] = static_cast<float>(-std::pow(det, 1.0f / 3));
		}
	}

	float calculatedDensity(point p) const override
	{
		const point l = p - mLocation;
		const clipper::Coord_orth dxyz(l.m_x, l.m_y, l.m_z);
		return mOccupancy *
		       static_cast<float>(
				   mAW[0] * std::exp(mAnisoInv[0].quad_form(dxyz)) +
				   mAW[1] * std::exp(mAnisoInv[1].quad_form(dxyz)) +
				   mAW[2] * std::exp(mAnisoInv[2].quad_form(dxyz)) +
				   mAW[3] * std::exp(mAnisoInv[3].quad_form(dxyz)) +
				   mAW[4] * std::exp(mAnisoInv[4].quad_form(dxyz)) +
				   mAW[5] * std::exp(mAnisoInv[5].quad_form(dxyz)));
	}

	clipper::U_aniso_orth mAnisoU;
	std::array<clipper::Mat33sym<>, 6> mAnisoInv;
};

// --------------------------------------------------------------------

AtomShape::AtomShape(cif::row_handle atom, cif::row_handle atom_aniso, float resHigh, float resLow, bool electronScattering, std::optional<float> bFactor)
	: mImpl(nullptr)
{
	const auto &[x, y, z, charge, type_symbol, occupancy, compound_id] =
		atom.get<float, float, float, std::optional<int>, std::string, float, std::string>(
			"Cartn_x", "Cartn_y", "Cartn_z", "pdbx_formal_charge", "type_symbol", "occupancy", "label_comp_id");

	atom_type type = type_symbol == "X" ? cif::Nn : atom_type_traits(type_symbol).type();

	int formal_charge = charge.value_or(0);
	if (not charge.has_value())
	{
		auto compound = cif::compound_factory::instance().create(compound_id);
		if (compound != nullptr and compound->atoms().size() == 1)
			formal_charge = compound->atoms().front().charge;
	}

	if (bFactor.has_value() and *bFactor != 0)
		mImpl = new AtomShapeImpl({ x, y, z }, type, formal_charge, static_cast<float>(clipper::Util::b2u(*bFactor)), 1.0, resHigh, resLow, electronScattering);
	else if (not atom_aniso.empty())
	{
		const auto &[u11, u12, u13, u22, u23, u33] =
			atom_aniso.get<float, float, float, float, float, float>("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");

		clipper::U_aniso_orth u(u11, u22, u33, u12, u13, u23);

		while (u.det() < 1.0e-20)
			u = u + clipper::U_aniso_orth(0.01, 0.01, 0.01, 0, 0, 0);

		mImpl = new AtomShapeAnisoImpl({ x, y, z }, type, formal_charge, u, occupancy, resHigh, resLow, electronScattering);
	}
	else
	{
		const auto &[u_iso, b_iso] = atom.get<std::optional<float>, std::optional<float>>("U_iso_or_equiv", "B_iso_or_equiv");

		float iso = 0;

		if (u_iso.has_value())
			iso = *u_iso;
		else if (b_iso.has_value())
			iso = *b_iso / static_cast<float>(8 * kPI * kPI);

		if (iso == 0)
			iso = 2.0f / static_cast<float>(8 * kPI * kPI);;

		mImpl = new AtomShapeImpl({ x, y, z }, type, formal_charge, iso, occupancy, resHigh, resLow, electronScattering);
	}
}

// AtomShape::AtomShape(const cif::mm::atom &atom, float resHigh, float resLow, bool electronScattering, float bFactor)
// 	: mImpl(new AtomShapeImpl(atom.get_location(), atom.type(), atom.charge(), static_cast<float>(clipper::Util::b2u(bFactor)),
// 		  1.0, resHigh, resLow, electronScattering))
// {
// }

AtomShape::~AtomShape()
{
	delete mImpl;
}

float AtomShape::radius() const
{
	return mImpl->integratedRadius(0.95f);
}

float AtomShape::calculatedDensity(float r) const
{
	return mImpl->calculatedDensity(r);
}

float AtomShape::calculatedDensity(point p) const
{
	return mImpl->calculatedDensity(p);
}

} // namespace pdb_redo
