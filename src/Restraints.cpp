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

#include "pdb-redo/Minimizer.hpp"
#include "pdb-redo/Restraints.hpp"

#include <cif++.hpp>

#include <Eigen/Eigenvalues>

#include <numeric>

namespace pdb_redo
{

double BondRestraint::f(const AtomLocationProvider &atoms) const
{
	double d = mDist - distance(atoms[mA], atoms[mB]);
	double result = (d * d) / (mDistESD * mDistESD);

	if (cif::VERBOSE > 2)
		std::cerr << "bond::f() = " << atoms.atom(mA) << " <> " << atoms.atom(mB)
				  << " " << atoms[mA] << " <> " << atoms[mB]
				  << " => " << result << std::endl;

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
				  << bi << ' ' << mDist << ' ' << mDistESD << std::endl;

	df.add(mA, (a1 - a2) * c);
	df.add(mB, (a2 - a1) * c);
}

void BondRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "bond " << atoms.atom(mA) << " to " << atoms.atom(mB) << " => " << mDist << " / " << mDistESD << std::endl;
}

// --------------------------------------------------------------------

double AngleRestraint::f(const AtomLocationProvider &atoms) const
{
	DPoint p[3] = {atoms[mA], atoms[mB], atoms[mC]};

	double c = cosinus_angle(p[1], p[0], p[1], p[2]);
	double angle = std::atan2(std::sqrt(1 - c * c), c) * 180 / cif::kPI;

	double d = mAngle - angle;
	double result = (d * d) / (mESD * mESD);

	if (cif::VERBOSE > 2)
		std::cerr << "angle::f() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << ' ' << " = " << result << std::endl;

	return result;
}

void AngleRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	const double kRadToDegree = 180.0 / cif::kPI, kDegreeToRad = 1 / kRadToDegree;

	if (cif::VERBOSE > 2)
		std::cerr << "angle::df() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << ' ' << ": " << std::endl;

	DPoint k = atoms[mA], l = atoms[mB], m = atoms[mC];

	auto aVec = (k - l);
	auto bVec = (m - l);

	auto a = distance(k, l);
	if (a < 0.01)
	{
		a = 0.01;
		aVec = DPoint{0.01, 0.01, 0.01};
	};

	auto b = distance(m, l);
	if (b < 0.01)
	{
		b = 0.01;
		bVec = DPoint{0.01, 0.01, 0.01};
	};

	auto cosTheta = dot_product(aVec, bVec) / (a * b);
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
	std::cout << "angle " << atoms.atom(mA) << " ; " << atoms.atom(mB) << " ; " << atoms.atom(mC) << " => " << mAngle << " / " << mESD << std::endl;
}

// --------------------------------------------------------------------

std::tuple<DPoint, DPoint, DPoint, DPoint>
TorsionRestraint::CalculateTorsionGradients(float theta, DPoint p[4]) const
{
	auto a = p[1] - p[0], b = p[2] - p[1], c = p[3] - p[2];

	auto blensq = b.length_sq();
	auto blen = std::sqrt(blensq);

	if (blen < 0.01)
	{
		blen = 0.01;
		blensq = 0.0001;
	}

	auto H = -dot_product(a, c),
		 J = dot_product(a, b),
		 K = dot_product(b, c),
		 L = 1 / blensq;

	auto E = dot_product(a, cross_product(b, c)) / blen;
	auto G = H + J * K * L;
	auto F = 1 / G;

	if (G == 0)
		F = 999999999.9;

	DPoint dH[4] = {c, -c, a, -a};
	DPoint dK[4] = {{}, -c, c - b, b};
	DPoint dJ[4] = {-b, b - a, a, {}};
	DPoint dL[4] = {{}, 2.0 * (p[2] - p[1]) * L * L, -2.0 * (p[2] - p[1]) * L * L, {}};
	DPoint dM[4] = {
		{-(b.m_y * c.m_z - b.m_z * c.m_y),
			-(b.m_z * c.m_x - b.m_x * c.m_z),
			-(b.m_x * c.m_y - b.m_y * c.m_x)},
		{(b.m_y * c.m_z - b.m_z * c.m_y) + (a.m_y * c.m_z - a.m_z * c.m_y),
			(b.m_z * c.m_x - b.m_x * c.m_z) + (a.m_z * c.m_x - a.m_x * c.m_z),
			(b.m_x * c.m_y - b.m_y * c.m_x) + (a.m_x * c.m_y - a.m_y * c.m_x)},
		{(b.m_y * a.m_z - b.m_z * a.m_y) - (a.m_y * c.m_z - a.m_z * c.m_y),
			-(a.m_z * c.m_x - a.m_x * c.m_z) + (b.m_z * a.m_x - b.m_x * a.m_z),
			-(a.m_x * c.m_y - a.m_y * c.m_x) + (a.m_y * b.m_x - a.m_x * b.m_y)},
		{-(b.m_y * a.m_z - b.m_z * a.m_y),
			-(b.m_z * a.m_x - b.m_x * a.m_z),
			-(a.m_y * b.m_x - a.m_x * b.m_y)}};

	DPoint dE[4]{
		dM[0] / blen,
		dM[1] / blen + E * (p[2] - p[1]) * L,
		dM[2] / blen - E * (p[2] - p[1]) * L,
		dM[3] / blen};

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

	double cos_a1 = cosinus_angle(atoms[mB], atoms[mA], atoms[mC], atoms[mB]);
	double cos_a2 = cosinus_angle(atoms[mC], atoms[mB], atoms[mD], atoms[mC]);

	if (cos_a1 <= 0.9 and cos_a2 <= 0.9)
	{
		double period = 360.0;
		if (mPeriodicity > 0)
			period /= mPeriodicity;

		double theta = dihedral_angle(atoms[mA], atoms[mB], atoms[mC], atoms[mD]);
		double diff = std::fmod(std::abs(theta - mTarget) + period / 2, period) - period / 2;

		if (not std::isnan(diff))
			result = (diff * diff) / (mESD * mESD);

		if (cif::VERBOSE > 2)
			std::cerr << "torsion::f() = " << result << " for theta " << theta
					  << " diff: " << diff
					  << " target: " << mTarget << " sigma: " << mESD
					  << " atoms: " << atoms.atom(mA) << ", " << atoms.atom(mB) << ", " << atoms.atom(mC) << ", " << atoms.atom(mD)
					  << std::endl;
	}

	return result;
}

void TorsionRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
		std::cerr << "torsion::df() " << atoms.atom(mA) << "/" << atoms.atom(mB) << "/" << atoms.atom(mC) << "/" << atoms.atom(mD) << ' ' << ": " << std::endl;

	double cos_a1 = cosinus_angle(atoms[mB], atoms[mA], atoms[mC], atoms[mB]);
	double cos_a2 = cosinus_angle(atoms[mC], atoms[mB], atoms[mD], atoms[mC]);

	if (cos_a1 <= 0.9 and cos_a2 <= 0.9)
	{
		double period = 360.0;
		if (mPeriodicity > 0)
			period /= mPeriodicity;

		double theta = dihedral_angle(atoms[mA], atoms[mB], atoms[mC], atoms[mD]);
		double diff = std::fmod(std::abs(theta - mTarget) + period / 2, period) - period / 2;

		if (not std::isnan(diff))
		{
			auto tt = std::tan(cif::kPI * theta / 180);
			double scale = 180.0 / ((1 + tt * tt) * cif::kPI);
			auto w = 1 / (mESD * mESD);

			DPoint p[4] = {atoms[mA], atoms[mB], atoms[mC], atoms[mD]};
			DPoint d[4];

			std::tie(d[0], d[1], d[2], d[3]) = CalculateTorsionGradients(theta, p);

			df.add(mA, 2.0 * diff * d[0] * scale * w);
			df.add(mB, 2.0 * diff * d[1] * scale * w);
			df.add(mC, 2.0 * diff * d[2] * scale * w);
			df.add(mD, 2.0 * diff * d[3] * scale * w);
		}
	}
}

void TorsionRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "torsion " << atoms.atom(mA) << " ; " << atoms.atom(mB) << " ; " << atoms.atom(mC) << " ; " << atoms.atom(mD) << " => " << mPeriodicity << " / " << mESD << std::endl;
}

// --------------------------------------------------------------------

double ChiralVolumeRestraint::f(const AtomLocationProvider &atoms) const
{
	auto chiralVolume = dot_product(atoms[mA1] - atoms[mCentre],
		cross_product(atoms[mA2] - atoms[mCentre], atoms[mA3] - atoms[mCentre]));

	double d = mVolume - chiralVolume;
	double result = (d * d) / (mESD * mESD);

	if (cif::VERBOSE > 2)
		std::cerr << "chiral::f() = " << result << std::endl;

	return result;
}

void ChiralVolumeRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
		std::cerr << "chiral::df(): " << std::endl;

	DPoint centre = atoms[mCentre];
	DPoint a = atoms[mA1] - centre;
	DPoint b = atoms[mA2] - centre;
	DPoint c = atoms[mA3] - centre;

	auto chiralVolume = dot_product(a, cross_product(b, c));

	auto d = chiralVolume - mVolume;
	auto s = 2 * d / (mESD * mESD);

	df.add(mCentre, s * DPoint{
							-(b.m_y * c.m_z - b.m_z * c.m_y) - (a.m_z * c.m_y - a.m_y * c.m_z) - (a.m_y * b.m_z - a.m_z * b.m_y),
							-(b.m_z * c.m_x - b.m_x * c.m_z) - (a.m_x * c.m_z - a.m_z * c.m_x) - (a.m_z * b.m_x - a.m_x * b.m_z),
							-(b.m_x * c.m_y - b.m_y * c.m_x) - (a.m_y * c.m_x - a.m_x * c.m_y) - (a.m_x * b.m_y - a.m_y * b.m_x)});

	df.add(mA1, s * DPoint{b.m_y * c.m_z - b.m_z * c.m_y, b.m_z * c.m_x - b.m_x * c.m_z, b.m_x * c.m_y - b.m_y * c.m_x});
	df.add(mA2, s * DPoint{a.m_z * c.m_y - a.m_y * c.m_z, a.m_x * c.m_z - a.m_z * c.m_x, a.m_y * c.m_x - a.m_x * c.m_y});
	df.add(mA3, s * DPoint{a.m_y * b.m_z - a.m_z * b.m_y, a.m_z * b.m_x - a.m_x * b.m_z, a.m_x * b.m_y - a.m_y * b.m_x});
}

void ChiralVolumeRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "chiral volume " << atoms.atom(mA1) << " ; " << atoms.atom(mA2) << " ; " << atoms.atom(mA3) << " => " << mVolume << " / " << mESD << std::endl;
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
		Cxx += (atoms[a].m_x - center.m_x) * (atoms[a].m_x - center.m_x);
		Cyy += (atoms[a].m_y - center.m_y) * (atoms[a].m_y - center.m_y);
		Czz += (atoms[a].m_z - center.m_z) * (atoms[a].m_z - center.m_z);
		Cxy += (atoms[a].m_x - center.m_x) * (atoms[a].m_y - center.m_y);
		Cxz += (atoms[a].m_x - center.m_x) * (atoms[a].m_z - center.m_z);
		Cyz += (atoms[a].m_y - center.m_y) * (atoms[a].m_z - center.m_z);
	}
	
	Eigen::Matrix3d mat;
	mat << Cxx, Cxy, Cxz,
		   Cxy, Cyy, Cyz,
		   Cxz, Cyz, Czz;

	Eigen::EigenSolver<Eigen::Matrix3d> es(mat);

	auto ev = es.eigenvalues();

	float b_ev = std::numeric_limits<float>::max();
	for (size_t i = 0; i < 3; ++i)
	{
		if (ev[i].real() > b_ev)
			continue;

		b_ev = ev[i].real();

		auto col = es.eigenvectors().col(i);

		abcd[0] = col(0).real();
		abcd[1] = col(1).real();
		abcd[2] = col(2).real();
	}

	double sumSq = 1e-20 + abcd[0] * abcd[0] + abcd[1] * abcd[1] + abcd[2] * abcd[2];

	abcd[0] /= sumSq;
	abcd[1] /= sumSq;
	abcd[2] /= sumSq;

	abcd[3] = abcd[0] * center.m_x + abcd[1] * center.m_y + abcd[2] * center.m_z;
}

double PlanarityRestraint::f(const AtomLocationProvider &atoms) const
{
	double abcd[4];

	calculatePlaneFunction(atoms, abcd);

	double result = accumulate(mAtoms.begin(), mAtoms.end(), 0.,
		[&atoms, &abcd, esd = mESD](double sum, AtomRef a)
		{
			double v = abcd[0] * atoms[a].m_x +
		               abcd[1] * atoms[a].m_y +
		               abcd[2] * atoms[a].m_z -
		               abcd[3];

			double r = v / esd;

			return sum + r * r;
		});

	if (cif::VERBOSE > 2)
	{
		std::vector<std::string> as;
		transform(mAtoms.begin(), mAtoms.end(), back_inserter(as),
			[](auto &a)
			{
				std::stringstream s;
				s << a;
				return s.str();
			});

		std::cerr << "plane::f() = " << result << " for " << mAtoms.size() << " atoms " << cif::join(as, ", ") << std::endl;
	}

	return result;
}

void PlanarityRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
	{
		std::vector<std::string> as;
		transform(mAtoms.begin(), mAtoms.end(), back_inserter(as),
			[](auto &a)
			{
				std::stringstream s;
				s << a;
				return s.str();
			});

		std::cerr << "plane::df() for " << mAtoms.size() << " atoms " << cif::join(as, ", ") << std::endl;
	}

	double abcd[4];

	calculatePlaneFunction(atoms, abcd);

	for (auto &a : mAtoms)
	{
		auto l = atoms[a];
		auto deviLen = l.m_x * abcd[0] + l.m_y * abcd[1] + l.m_z * abcd[2] - abcd[3];

		df.add(a, 2 * deviLen * DPoint{abcd[0], abcd[1], abcd[2]} / (mESD * mESD));
	}
}

void PlanarityRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "plane ";

	for (auto &a : mAtoms)
		std::cout << atoms.atom(a) << ' ';

	std::cout << "=> " << 0 << " / " << mESD << std::endl;
}

// --------------------------------------------------------------------

double NonBondedContactRestraint::f(const AtomLocationProvider &atoms) const
{
	double result = 0;

	double distance = distance_squared(atoms[mA], atoms[mB]);
	if (distance < mMinDistSq)
	{
		double d = mMinDist - std::sqrt(distance);
		result = (d * d) / (mDistESD * mDistESD);

		if (cif::VERBOSE > 2)
			std::cerr << "non-bonded-contact::f() = " << result << " min-dist is " << mMinDist << " and dist is " << std::sqrt(distance)
					<< " a1: " << atoms.atom(mA) << " a2: " << atoms.atom(mB) << std::endl
					<< " a1: " << atoms[mA] << " a2: " << atoms[mB] << std::endl;
	}


	return result;
}

void NonBondedContactRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	auto a1 = atoms[mA], a2 = atoms[mB];

	auto bi = distance_squared(a1, a2);
	if (bi < mMinDistSq)
	{
		bi = std::sqrt(bi);
		if (bi < 0.1)
			bi = 0.1;

		if (cif::VERBOSE > 2)
			std::cerr << "non-bonded::df(): " << atoms.atom(mA) << " and " << atoms.atom(mB) << " "
					  << "distance: " << bi << " "
					  << "target: " << mMinDist << std::endl;

		double c = 2 * (1 - mMinDist / bi) / (mDistESD * mDistESD);

		df.add(mA, (a1 - a2) * c);
		df.add(mB, (a2 - a1) * c);
	}
}

void NonBondedContactRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "nbc " << atoms.atom(mA) << " " << atoms.atom(mB)
			  << " => " << distance(atoms[mA], atoms[mB]) << ' ' << mMinDist << " / " << mDistESD << std::endl;
}

// --------------------------------------------------------------------

DensityRestraint::DensityRestraint(std::vector<std::pair<AtomRef, double>> &&atoms,
	const Xmap &xMap, double mapWeight)
	: mAtoms(move(atoms))
	, mXMap(xMap)
	, mMapWeight(mapWeight)
{
}

double DensityRestraint::f(const AtomLocationProvider &atoms) const
{
	double result = 0;

	for (auto &a : mAtoms)
	{
		clipper::Coord_orth p{atoms[a.first].m_x, atoms[a.first].m_y, atoms[a.first].m_z };
		clipper::Coord_frac pf = p.coord_frac(mXMap.cell());

		result += a.second * mXMap.interp<clipper::Interp_cubic>(pf);
	}

	if (cif::VERBOSE > 2)
		std::cerr << "density::f() = " << -result << std::endl;

	return mMapWeight * -result;
}

void DensityRestraint::df(const AtomLocationProvider &atoms, DFCollector &df) const
{
	if (cif::VERBOSE > 2)
		std::cerr << "density::df(): " << std::endl;

	for (auto &a : mAtoms)
	{
		clipper::Coord_orth p{atoms[a.first].m_x, atoms[a.first].m_y, atoms[a.first].m_z };
		clipper::Coord_frac pf = p.coord_frac(mXMap.cell());
		auto pm = pf.coord_map(mXMap.grid_sampling());

		clipper::Grad_map<double> grad;
		double dv;

		clipper::Interp_cubic::interp_grad(mXMap, pm, dv, grad);
		auto gradFrac = grad.grad_frac(mXMap.grid_sampling());

		auto gradOrth = gradFrac.grad_orth(mXMap.cell());

		df.add(a.first, DPoint{gradOrth.dx(), gradOrth.dy(), gradOrth.dz()} * mMapWeight * -a.second);
	}
}

void DensityRestraint::print(const AtomLocationProvider &atoms) const
{
	std::cout << "density " << std::endl;
}

} // namespace pdb_redo