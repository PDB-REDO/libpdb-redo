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

#include "pdb-redo/Restraints.hpp"

#include "pdb-redo/Minimizer.hpp"

#include <algorithm>
#include <cif++/cif++.hpp>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <numeric>

namespace pdb_redo
{

double BondRestraint::f(const AtomLocationProvider &atoms) const
{
	double d = mDist - distance(atoms[mA], atoms[mB]);
	double result = (d * d) / (mDistESD * mDistESD);

	return result;
}

void BondRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	auto a1 = atoms[mA], a2 = atoms[mB];

	auto bi = distance(a1, a2);
	if (bi < 0.1)
		bi = 0.1;

	auto c = 2 * (1 - mDist / bi) / (mDistESD * mDistESD);

	if (cif::VERBOSE > 2)
		std::cerr << "bond::df(): " << atoms.atom(mA) << " <> " << atoms.atom(mB) << ' '
				  << bi << ' ' << mDist << ' ' << mDistESD << '\n';

	df.add(mA, (a1 - a2) * c);
	df.add(mB, (a2 - a1) * c);
}

void BondRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "bond " << atoms.atom(mA) << " to " << atoms.atom(mB) << " => " << mDist << " / " << mDistESD << '\n';
}

// --------------------------------------------------------------------

double AngleRestraint::f(const AtomLocationProvider &atoms) const
{
	DPoint p[3] = { atoms[mA], atoms[mB], atoms[mC] };

	double c = cif::cosinus_angle(p[1], p[0], p[1], p[2]);
	double angle = std::atan2(std::sqrt(1 - c * c), c) * 180 / std::numbers::pi;

	double d = mAngle - angle;
	double result = (d * d) / (mESD * mESD);

	if (cif::VERBOSE > 2)
		std::cerr << "angle::f() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << ' ' << " = " << result << '\n';

	return result;
}

void AngleRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	const double kRadToDegree = 180.0 / std::numbers::pi, kDegreeToRad = 1 / kRadToDegree;

	if (cif::VERBOSE > 2)
		std::cerr << "angle::df() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << ' ' << ":\n";

	DPoint k = atoms[mA], l = atoms[mB], m = atoms[mC];

	auto aVec = (k - l);
	auto bVec = (m - l);

	auto a = distance(k, l);
	if (a < 0.01)
	{
		a = 0.01;
		aVec = DPoint{ 0.01, 0.01, 0.01 };
	};

	auto b = distance(m, l);
	if (b < 0.01)
	{
		b = 0.01;
		bVec = DPoint{ 0.01, 0.01, 0.01 };
	};

	auto cosTheta = glm::dot(aVec, bVec) / (a * b);
	if (cosTheta > 1.0)
		cosTheta = 1.0;
	if (cosTheta < -1.0)
		cosTheta = -1.0;
	auto theta = std::acos(cosTheta);
	if (theta < 0.001)
		theta = 0.001;

	auto target = mAngle * kDegreeToRad;

	auto wf = 2 * (theta - target) * kRadToDegree * kRadToDegree / (mESD * mESD);
	auto prem = -wf / std::sin(theta);

	df.add(mA, prem * (cosTheta * (l - k) / (a * a) + (m - l) / (a * b)));
	df.add(mC, prem * (cosTheta * (l - m) / (b * b) + (k - l) / (a * b)));

	auto term1 = (l - k) * -cosTheta / (a * a) + (l - m) * -cosTheta / (b * b);
	auto term2 = ((l - k) + (l - m)) / (a * b);

	df.add(mB, prem * (term1 + term2));
}

void AngleRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "angle " << atoms.atom(mA) << " ; " << atoms.atom(mB) << " ; " << atoms.atom(mC) << " => " << mAngle << " / " << mESD << '\n';
}

// --------------------------------------------------------------------

std::tuple<DPoint, DPoint, DPoint, DPoint>
TorsionRestraint::CalculateTorsionGradients(double theta, DPoint p[4]) const
{
	auto a = p[1] - p[0], b = p[2] - p[1], c = p[3] - p[2];

	auto blen = std::sqrt(glm::length(b));

	if (blen < 0.01)
		blen = 0.01;

	auto H = -glm::dot(a, c),
		 J = glm::dot(a, b),
		 K = glm::dot(b, c),
		 L = 1 / (blen * blen);

	auto E = glm::dot(a, glm::cross(b, c)) / blen;
	auto G = H + J * K * L;
	auto F = 1 / G;

	if (G == 0)
		F = 999999999.9;

	DPoint dH[4] = { c, -c, a, -a };
	DPoint dK[4] = { {}, -c, c - b, b };
	DPoint dJ[4] = { -b, b - a, a, {} };
	DPoint dL[4] = { {}, 2.0 * (p[2] - p[1]) * L * L, -2.0 * (p[2] - p[1]) * L * L, {} };
	DPoint dM[4] = {
		{ -(b.y * c.z - b.z * c.y),
			-(b.z * c.x - b.x * c.z),
			-(b.x * c.y - b.y * c.x) },
		{ (b.y * c.z - b.z * c.y) + (a.y * c.z - a.z * c.y),
			(b.z * c.x - b.x * c.z) + (a.z * c.x - a.x * c.z),
			(b.x * c.y - b.y * c.x) + (a.x * c.y - a.y * c.x) },
		{ (b.y * a.z - b.z * a.y) - (a.y * c.z - a.z * c.y),
			-(a.z * c.x - a.x * c.z) + (b.z * a.x - b.x * a.z),
			-(a.x * c.y - a.y * c.x) + (a.y * b.x - a.x * b.y) },
		{ -(b.y * a.z - b.z * a.y),
			-(b.z * a.x - b.x * a.z),
			-(a.y * b.x - a.x * b.y) }
	};

	DPoint dE[4]{
		dM[0] / blen,
		dM[1] / blen + E * (p[2] - p[1]) * L,
		dM[2] / blen - E * (p[2] - p[1]) * L,
		dM[3] / blen
	};

	auto eff = E * F * F;
	auto jl = J * L;
	auto kl = K * L;
	auto jk = J * K;

	return std::make_tuple(
		F * dE[0] - eff * (dH[0] + jl * dK[0] + kl * dJ[0] + jk * dL[0]),
		F * dE[1] - eff * (dH[1] + jl * dK[1] + kl * dJ[1] + jk * dL[1]),
		F * dE[2] - eff * (dH[2] + jl * dK[2] + kl * dJ[2] + jk * dL[2]),
		F * dE[3] - eff * (dH[3] + jl * dK[3] + kl * dJ[3] + jk * dL[3]));
}

double TorsionRestraint::f(const AtomLocationProvider &atoms) const
{
	double result = 0;

	double cos_a1 = cif::cosinus_angle(atoms[mB], atoms[mA], atoms[mC], atoms[mB]);
	double cos_a2 = cif::cosinus_angle(atoms[mC], atoms[mB], atoms[mD], atoms[mC]);

	if (cos_a1 <= 0.9 and cos_a2 <= 0.9)
	{
		double period = 360.0;
		if (mPeriodicity > 0)
			period /= mPeriodicity;

		double theta = cif::dihedral_angle(atoms[mA], atoms[mB], atoms[mC], atoms[mD]);
		double diff = std::fmod(std::abs(theta - mTarget) + period / 2, period) - period / 2;

		if (not std::isnan(diff))
			result = (diff * diff) / (mESD * mESD);

		if (cif::VERBOSE > 2)
			std::cerr << "torsion::f() = " << result << " for theta " << theta
					  << " diff: " << diff
					  << " target: " << mTarget << " sigma: " << mESD
					  << " atoms: " << atoms.atom(mA) << ", " << atoms.atom(mB) << ", " << atoms.atom(mC) << ", " << atoms.atom(mD)
					  << '\n';
	}

	return result;
}

void TorsionRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
		std::cerr << "torsion::df() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << "/" << atoms.atom(mD) << ' ' << ":\n";

	double cos_a1 = cif::cosinus_angle(atoms[mB], atoms[mA], atoms[mC], atoms[mB]);
	double cos_a2 = cif::cosinus_angle(atoms[mC], atoms[mB], atoms[mD], atoms[mC]);

	if (cos_a1 <= 0.9 and cos_a2 <= 0.9)
	{
		double period = 360.0;
		if (mPeriodicity > 0)
			period /= mPeriodicity;

		double theta = cif::dihedral_angle(atoms[mA], atoms[mB], atoms[mC], atoms[mD]);
		double diff = std::fmod(std::abs(theta - mTarget) + period / 2, period) - period / 2;

		if (not std::isnan(diff))
		{
			auto tt = std::tan(std::numbers::pi * theta / 180);
			double scale = 180.0 / ((1 + tt * tt) * std::numbers::pi);
			auto w = 1 / (mESD * mESD);

			DPoint p[4] = { atoms[mA], atoms[mB], atoms[mC], atoms[mD] };
			DPoint d[4];

			std::tie(d[0], d[1], d[2], d[3]) = CalculateTorsionGradients(theta, p);

			df.add(mA, /* 2.0 *  */ diff * d[0] * scale * w);
			df.add(mB, /* 2.0 *  */ diff * d[1] * scale * w);
			df.add(mC, /* 2.0 *  */ diff * d[2] * scale * w);
			df.add(mD, /* 2.0 *  */ diff * d[3] * scale * w);
		}
	}
}

void TorsionRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "torsion " << atoms.atom(mA) << " ; " << atoms.atom(mB) << " ; " << atoms.atom(mC) << " ; " << atoms.atom(mD) << " => " << mPeriodicity << " / " << mESD << '\n';
}

// --------------------------------------------------------------------

double ChiralVolumeRestraint::f(const AtomLocationProvider &atoms) const
{
	auto chiralVolume = glm::dot(atoms[mA1] - atoms[mCentre],
		glm::cross(atoms[mA2] - atoms[mCentre], atoms[mA3] - atoms[mCentre]));

	double d = mVolume - chiralVolume;
	double result = (d * d) / (mESD * mESD);

	if (cif::VERBOSE > 2)
		std::cerr << "chiral::f() = " << result << '\n';

	return result;
}

void ChiralVolumeRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
		std::cerr << "chiral::df():\n";

	DPoint centre = atoms[mCentre];
	DPoint a = atoms[mA1] - centre;
	DPoint b = atoms[mA2] - centre;
	DPoint c = atoms[mA3] - centre;

	auto chiralVolume = glm::dot(a, glm::cross(b, c));

	auto d = chiralVolume - mVolume;
	auto s = 2 * d / (mESD * mESD);

	df.add(mCentre, s * DPoint{
							-(b.y * c.z - b.z * c.y) - (a.z * c.y - a.y * c.z) - (a.y * b.z - a.z * b.y),
							-(b.z * c.x - b.x * c.z) - (a.x * c.z - a.z * c.x) - (a.z * b.x - a.x * b.z),
							-(b.x * c.y - b.y * c.x) - (a.y * c.x - a.x * c.y) - (a.x * b.y - a.y * b.x) });

	df.add(mA1, s * DPoint{ b.y * c.z - b.z * c.y, b.z * c.x - b.x * c.z, b.x * c.y - b.y * c.x });
	df.add(mA2, s * DPoint{ a.z * c.y - a.y * c.z, a.x * c.z - a.z * c.x, a.y * c.x - a.x * c.y });
	df.add(mA3, s * DPoint{ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x });
}

void ChiralVolumeRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "chiral volume " << atoms.atom(mA1) << " ; " << atoms.atom(mA2) << " ; " << atoms.atom(mA3) << " => " << mVolume << " / " << mESD << '\n';
}

// --------------------------------------------------------------------

void PlanarityRestraint::calculatePlaneFunction(const AtomLocationProvider &atoms, double abcd[4]) const
{
	DPoint center;

	for (auto &a : mAtoms)
		center += atoms[a];
	center /= mAtoms.size();

	double Cxx = 0, Cyy = 0, Czz = 0, Cxy = 0, Cxz = 0, Cyz = 0;

	for (auto &a : mAtoms)
	{
		Cxx += (atoms[a].x - center.x) * (atoms[a].x - center.x);
		Cyy += (atoms[a].y - center.y) * (atoms[a].y - center.y);
		Czz += (atoms[a].z - center.z) * (atoms[a].z - center.z);
		Cxy += (atoms[a].x - center.x) * (atoms[a].y - center.y);
		Cxz += (atoms[a].x - center.x) * (atoms[a].z - center.z);
		Cyz += (atoms[a].y - center.y) * (atoms[a].z - center.z);
	}

	double data[9] = {
		Cxx, Cxy, Cxz,
		Cxy, Cyy, Cyz,
		Cxz, Cyz, Czz
	};

	gsl_matrix_view m = gsl_matrix_view_array(data, 3, 3);

	gsl_vector *eval = gsl_vector_alloc(3);
	gsl_matrix *evec = gsl_matrix_alloc(3, 3);

	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);

	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	gsl_vector_view evec_i = gsl_matrix_column(evec, 0);

	abcd[0] = gsl_vector_get(&evec_i.vector, 0);
	abcd[1] = gsl_vector_get(&evec_i.vector, 1);
	abcd[2] = gsl_vector_get(&evec_i.vector, 2);

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	double sumSq = 1e-20 + abcd[0] * abcd[0] + abcd[1] * abcd[1] + abcd[2] * abcd[2];

	abcd[0] /= sumSq;
	abcd[1] /= sumSq;
	abcd[2] /= sumSq;

	abcd[3] = abcd[0] * center.x + abcd[1] * center.y + abcd[2] * center.z;
}

double PlanarityRestraint::f(const AtomLocationProvider &atoms) const
{
	double abcd[4];

	calculatePlaneFunction(atoms, abcd);

	double result = accumulate(mAtoms.begin(), mAtoms.end(), 0.,
		[&atoms, &abcd, esd = mESD](double sum, AtomRef a)
		{
			double v = abcd[0] * atoms[a].x +
		               abcd[1] * atoms[a].y +
		               abcd[2] * atoms[a].z -
		               abcd[3];

			double r = v / esd;

			return sum + r * r;
		});

	if (cif::VERBOSE > 2)
	{
		std::vector<std::string> as;
		std::ranges::transform(mAtoms, back_inserter(as),
			[](auto &a)
			{
				std::stringstream s;
				s << a;
				return s.str();
			});

		std::cerr << "plane::f() = " << result << " for " << mAtoms.size() << " atoms " << cif::join(as, ", ") << '\n';
	}

	return result;
}

void PlanarityRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
	{
		std::vector<std::string> as;
		std::ranges::transform(mAtoms, back_inserter(as),
			[](auto &a)
			{
				std::stringstream s;
				s << a;
				return s.str();
			});

		std::cerr << "plane::df() for " << mAtoms.size() << " atoms " << cif::join(as, ", ") << '\n';
	}

	double abcd[4];

	calculatePlaneFunction(atoms, abcd);

	for (auto &a : mAtoms)
	{
		auto l = atoms[a];
		auto deviLen = l.x * abcd[0] + l.y * abcd[1] + l.z * abcd[2] - abcd[3];

		df.add(a, 2 * deviLen * DPoint{ abcd[0], abcd[1], abcd[2] } / (mESD * mESD));
	}
}

void PlanarityRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "plane ";

	for (auto &a : mAtoms)
		std::cout << atoms.atom(a) << ' ';

	std::cout << "=> " << 0 << " / " << mESD << '\n';
}

// --------------------------------------------------------------------

double NonBondedContactRestraint::f(const AtomLocationProvider &atoms) const
{
	double result = 0;

	double distance = cif::distance_squared(atoms[mA], atoms[mB]);
	if (distance < mMinDistSq)
	{
		double d = mMinDist - std::sqrt(distance);
		result = (d * d) / (mDistESD * mDistESD);
	}

	return result;
}

void NonBondedContactRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	auto a1 = atoms[mA], a2 = atoms[mB];

	auto bi = cif::distance_squared(a1, a2);
	if (bi < mMinDistSq)
	{
		bi = std::sqrt(bi);
		if (bi < 0.1)
			bi = 0.1;

		if (cif::VERBOSE > 2)
			std::cerr << "non-bonded::df(): " << atoms.atom(mA) << " and " << atoms.atom(mB) << " "
					  << "distance: " << bi << " "
					  << "target: " << mMinDist << '\n';

		double c = 2 * (1 - mMinDist / bi) / (mDistESD * mDistESD);

		df.add(mA, (a1 - a2) * c);
		df.add(mB, (a2 - a1) * c);
	}
}

void NonBondedContactRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "nbc " << atoms.atom(mA) << " " << atoms.atom(mB)
			  << " => " << distance(atoms[mA], atoms[mB]) << ' ' << mMinDist << " / " << mDistESD << '\n';
}

// --------------------------------------------------------------------

DensityRestraint::DensityRestraint(std::vector<std::pair<AtomRef, double>> &&atoms,
	const Xmap &xMap, double mapWeight)
	: mAtoms(std::move(atoms))
	, mXMap(xMap)
	, mMapWeight(mapWeight)
{
}

double DensityRestraint::f(const AtomLocationProvider &atoms) const
{
	double result = 0;

	for (auto &a : mAtoms)
	{
		clipper::Coord_orth p{ atoms[a.first].x, atoms[a.first].y, atoms[a.first].z };
		clipper::Coord_frac pf = p.coord_frac(mXMap.cell());

		result += a.second * mXMap.interp<clipper::Interp_cubic>(pf);
	}

	if (cif::VERBOSE > 2)
		std::cerr << "density::f() = " << -result << '\n';

	return mMapWeight * -result;
}

void DensityRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
		std::cerr << "density::df():\n";

	for (auto &a : mAtoms)
	{
		clipper::Coord_orth p{ atoms[a.first].x, atoms[a.first].y, atoms[a.first].z };
		clipper::Coord_frac pf = p.coord_frac(mXMap.cell());
		auto pm = pf.coord_map(mXMap.grid_sampling());

		clipper::Grad_map<double> grad;
		double dv;

		clipper::Interp_cubic::interp_grad(mXMap, pm, dv, grad);
		auto gradFrac = grad.grad_frac(mXMap.grid_sampling());

		auto gradOrth = gradFrac.grad_orth(mXMap.cell());

		df.add(a.first, DPoint{ gradOrth.dx(), gradOrth.dy(), gradOrth.dz() } * mMapWeight * -a.second);
	}
}

void DensityRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "density\n";
}

} // namespace pdb_redo