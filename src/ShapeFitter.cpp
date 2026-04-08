/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "pdb-redo/ShapeFitter.hpp"

#include "pdb-redo/Minimizer.hpp"
#include "pdb-redo/Restraints.hpp"

#include <algorithm>
#include <cif++/atom_type.hpp>
#include <cif++/compound.hpp>
#include <cif++/matrix.hpp>
#include <cif++/model.hpp>
#include <cif++/point.hpp>
#include <cif++/symmetry.hpp>
#include <clipper/core/clipper_types.h>
#include <clipper/core/coords.h>
#include <cmath>
#include <glm/glm.hpp>
#include <gsl/gsl_blas.h> // for debugging norm of gradient
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector_double.h>
#include <memory>

namespace pdb_redo
{

auto createInertiaTensorForBlob(const std::vector<cif::point> pts, clipper::Xmap<float> &xmap)
{
	std::array<float, 6> If{};

	auto [c, r] = cif::smallest_sphere_around_points(pts);

	for (auto pt : pts)
	{
		clipper::Coord_orth cp{ pt.m_x, pt.m_y, pt.m_z };
		clipper::Coord_frac pf = cp.coord_frac(xmap.cell());
		auto dp = xmap.interp<clipper::Interp_cubic>(pf);

		pt -= c;

		If[0] += dp * (pt.m_y * pt.m_y + pt.m_z * pt.m_z); // 11
		If[1] += dp * (pt.m_x * pt.m_x + pt.m_z * pt.m_z); // 22
		If[2] += dp * (pt.m_x * pt.m_x + pt.m_y * pt.m_y); // 33
		If[3] -= dp * pt.m_x * pt.m_y;                     // 12
		If[4] -= dp * pt.m_x * pt.m_z;                     // 13
		If[5] -= dp * pt.m_y * pt.m_z;                     // 23
	}

	return glm::mat3{
		glm::normalize(glm::vec3{ If[0], If[3], If[4] }),
		glm::normalize(glm::vec3{ If[3], If[1], If[5] }),
		glm::normalize(glm::vec3{ If[4], If[5], If[2] })
	};
}

auto createInertiaTensorForLigand(const cif::mm::residue &res)
{
	std::array<float, 6> If{};

	auto [c, r] = cif::smallest_sphere_around_points(
		res.atoms() | std::views::transform(&cif::mm::atom::get_location) | std::ranges::to<std::vector>());

	for (auto atom : res.atoms())
	{
		auto pt = atom.get_location();
		pt -= c;

		cif::atom_type_traits t(atom.get_type());

		auto dp = t.weight();

		If[0] += dp * (pt.m_y * pt.m_y + pt.m_z * pt.m_z); // 11
		If[1] += dp * (pt.m_x * pt.m_x + pt.m_z * pt.m_z); // 22
		If[2] += dp * (pt.m_x * pt.m_x + pt.m_y * pt.m_y); // 33
		If[3] -= dp * pt.m_x * pt.m_y;                     // 12
		If[4] -= dp * pt.m_x * pt.m_z;                     // 13
		If[5] -= dp * pt.m_y * pt.m_z;                     // 23
	}

	return glm::mat3{
		glm::normalize(glm::vec3{ If[0], If[3], If[4] }),
		glm::normalize(glm::vec3{ If[3], If[1], If[5] }),
		glm::normalize(glm::vec3{ If[4], If[5], If[2] })
	};
}

auto principalAxis(const glm::mat3 &m)
{
	double data[9] = {
		m[0][0], m[0][1], m[0][2],
		m[0][1], m[1][1], m[1][2],
		m[0][2], m[1][2], m[2][2]
	};

	gsl_matrix_view g = gsl_matrix_view_array(data, 3, 3);

	gsl_vector *eval = gsl_vector_alloc(3);
	gsl_matrix *evec = gsl_matrix_alloc(3, 3);

	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv(&g.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);

	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	gsl_vector_view evec_i = gsl_matrix_column(evec, 0);

	cif::point result
	{
		static_cast<float>(gsl_vector_get(&evec_i.vector, 0)),
		static_cast<float>(gsl_vector_get(&evec_i.vector, 1)),
		static_cast<float>(gsl_vector_get(&evec_i.vector, 2))
	};

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	result.normalize();

	return result;
}

// --------------------------------------------------------------------

// Locate the single blob in the xmap
std::vector<clipper::Coord_grid> findSingleBlob(clipper::Xmap<float> &xmap, bool removeColinear)
{
	auto &sg = xmap.spacegroup();

	struct Vec3Less
	{
		bool operator()(const clipper::Vec3<int> &a, const clipper::Vec3<int> &b) const
		{
			int d = a[0] - b[0];
			if (d == 0)
				d = a[1] - b[1];
			if (d == 0)
				d = a[2] - b[2];
			return d < 0;
		}
	};

	// Minimal blob finding algo
	std::set<clipper::Coord_grid, Vec3Less> gridPoints;
	std::stack<clipper::Coord_grid> stack;

	for (clipper::Xmap<float>::Map_reference_coord i(xmap); not i.last(); i.next())
	{
		if (xmap[i] <= 0)
			continue;

		if (int symNr = i.sym(); symNr == 0)
			stack.push(i.coord());
		else
		{
			clipper::Coord_map cm(i.coord());
			cm = clipper::Coord_map(sg.symop(symNr).inverse() * cm);
			stack.push(cm.coord_grid());
		}

		break;
	}

	while (not stack.empty())
	{
		auto p = stack.top();
		stack.pop();

		clipper::Xmap<float>::Map_reference_coord iw(xmap, p);

		gridPoints.insert(p);

		for (int du : { -1, 0, 1 })
			for (int dv : { -1, 0, 1 })
				for (int dw : { -1, 0, 1 })
				{
					if (du == 0 and dv == 0 and dw == 0)
						continue;

					clipper::Coord_grid g(p.u() + du, p.v() + dv, p.w() + dw);
					clipper::Xmap<float>::Map_reference_coord iw(xmap, g);

					if (xmap[iw] == 0)
						continue;

					// if (iw.sym() != 0)
					// {
					// 	// Move point into the correct symmetry
					// 	clipper::Coord_map cm(g);
					// 	cm = clipper::Coord_map(sg.symop(iw.sym()).inverse() * cm);
					// 	g = cm.coord_grid();
					// }

					if (not gridPoints.contains(g))
						stack.push(g);
				}
	}

	std::vector<clipper::Coord_grid> result(gridPoints.begin(), gridPoints.end());

	// Very simplistic, only include the outer edges of the points that share an axis
	if (removeColinear)
	{
		std::set<float> vx, vy, vz;
		for (auto &p : result)
		{
			vx.insert(p[0]);
			vy.insert(p[1]);
			vz.insert(p[2]);
		}

		for (float x : vx)
		{
			for (float y : vy)
			{
				float min_z = 0, max_z = 0;
				for (bool first = true; auto &p : result)
				{
					if (p[0] != x or p[1] != y)
						continue;

					if (std::exchange(first, false))
						min_z = max_z = p[2];
					else
					{
						if (min_z > p[2])
							min_z = p[2];
						if (max_z < p[2])
							max_z = p[2];
					}
				}

				std::erase_if(result, [x, y, min_z, max_z](const clipper::Coord_grid &p)
					{ return p[0] == x and p[1] == y and p[2] > min_z and p[2] < max_z; });
			}
		}

		for (float x : vx)
		{
			for (float z : vz)
			{
				float min_y = 0, max_y = 0;
				for (bool first = true; auto &p : result)
				{
					if (p[0] != x or p[2] != z)
						continue;

					if (std::exchange(first, false))
						min_y = max_y = p[1];
					else
					{
						if (min_y > p[1])
							min_y = p[1];
						if (max_y < p[1])
							max_y = p[1];
					}
				}

				std::erase_if(result, [x, z, min_y, max_y](const clipper::Coord_grid &p)
					{ return p[0] == x and p[2] == z and p[1] > min_y and p[1] < max_y; });
			}
		}

		for (float y : vy)
		{
			for (float z : vz)
			{
				float min_x = 0, max_x = 0;
				for (bool first = true; auto &p : result)
				{
					if (p[1] != y or p[2] != z)
						continue;

					if (std::exchange(first, false))
						min_x = max_x = p[0];
					else
					{
						if (min_x > p[0])
							min_x = p[0];
						if (max_x < p[0])
							max_x = p[0];
					}
				}

				std::erase_if(result, [y, z, min_x, max_x](const clipper::Coord_grid &p)
					{ return p[1] == y and p[2] == z and p[0] > min_x and p[0] < max_x; });
			}
		}
	}

	return result;
}

// --------------------------------------------------------------------

class JiggleFitter
{
  public:
	JiggleFitter(cif::mm::residue &res, clipper::Xmap<float> &xmap, float mapWeight = 45.f);

	double refine();
	virtual void transform();

  protected:
	std::vector<cif::mm::atom> mAtoms;
	std::vector<cif::point> mLocations;

	cif::point mCenter;

	std::vector<double> mVariables;
	std::vector<float> mStepSizes;

	std::unique_ptr<pdb_redo::DensityRestraint> mDensityDestraint;

	static double F(const gsl_vector *v, void *params)
	{
		auto *self = reinterpret_cast<JiggleFitter *>(params);
		return self->F(v);
	}

	double F(const gsl_vector *v)
	{
		for (size_t i = 0; i < mVariables.size(); ++i)
			mVariables[i] = gsl_vector_get(v, i);

		transform();

		pdb_redo::AtomLocationProvider loc(mAtoms);

		auto result = mDensityDestraint->f(loc);

		return result;
	}
};

JiggleFitter::JiggleFitter(cif::mm::residue &res, clipper::Xmap<float> &xmap, float mapWeight)
	: mAtoms(res.atoms())
{
	for (auto &atom : mAtoms)
		mLocations.emplace_back(atom.get_location());

	std::tie(mCenter, std::ignore) = cif::smallest_sphere_around_points(mLocations);

	std::vector<std::pair<pdb_redo::AtomRef, double>> densityAtoms;
	densityAtoms.reserve(mAtoms.size());

	std::ranges::transform(mAtoms, std::back_inserter(densityAtoms),
		[&densityAtoms](const cif::mm::atom &a)
		{
			double z = static_cast<int>(a.get_type());
			double weight = 1;
			double occupancy = a.get_occupancy();

			if (occupancy > 1)
				occupancy = 1;

			// TODO: cryo_em support
			return std::make_pair(densityAtoms.size(), z * weight * occupancy);
		});

	mDensityDestraint = std::make_unique<DensityRestraint>(std::move(densityAtoms), xmap, mapWeight);

	// This is where to set step sizes that are used during the rigid body fit
	mVariables = { 0, 0, 0, 0, 0 };
	mStepSizes = { 10.f, 10.f, 0.25f, 0.25f, 0.25f };
}

void JiggleFitter::transform()
{
	// get the variables assigned
	// NOLINTNEXTLINE(bugprone-narrowing-conversions)
	const float alpha = mVariables[0], beta = mVariables[1], x = mVariables[2], y = mVariables[3], z = mVariables[4];

	// rotations
	auto q0 = cif::construct_from_angle_axis(alpha, { 1, 0, 0 }); // construct quaternion from float angle, point axis
	auto q1 = cif::construct_from_angle_axis(beta, { 0, 0, 1 });

	auto q = q0 * q1;

	// translations
	cif::point translation(x, y, z);

	// Move the molecule
	for (size_t i = 0; i < mLocations.size(); ++i)
	{
		auto a = mLocations[i];
		a.rotate(q, mCenter); // rotate a using quaternion q01, move it to the alpha_loc - rotate - move back
		mAtoms[i].set_location(a + translation);
	}
}

double JiggleFitter::refine()
{
	const int kMaxIterations = 500;

	gsl_multimin_function f = {
		.f = &JiggleFitter::F,
		.n = mVariables.size(),
		.params = this
	};

	auto T = gsl_multimin_fminimizer_nmsimplex;

	auto x = gsl_vector_alloc(mVariables.size());
	for (size_t i = 0; i < mVariables.size(); ++i)
		gsl_vector_set(x, i, mVariables[i]);

	auto ss = gsl_vector_alloc(mStepSizes.size());
	for (size_t i = 0; i < mStepSizes.size(); ++i)
		gsl_vector_set(ss, i, mStepSizes[i]);

	auto s = gsl_multimin_fminimizer_alloc(T, mVariables.size());

	gsl_multimin_fminimizer_set(s, &f, x, ss);

	int iter = 0, status;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		double size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-3);

		if (status == GSL_SUCCESS)
		{
			if (cif::VERBOSE > 1)
				std::cout << "Minimum reached after " << iter << " iterations" << std::endl;
		}
	} while (status == GSL_CONTINUE and iter < kMaxIterations);

	for (size_t i = 0; i < mVariables.size(); ++i)
		mVariables[i] = gsl_vector_get(s->x, i);

	transform();

	gsl_vector_free(x);
	gsl_vector_free(ss);

	auto result = s->fval;

	gsl_multimin_fminimizer_free(s);

	return result;
}

// --------------------------------------------------------------------

cif::point centerOfMassBlob(clipper::Xmap<float> &xmap, const std::vector<cif::point> &blob)
{
	double sumMass = 0;
	std::array<double, 3> c{ 0, 0, 0  };
	for (auto pt : blob)
	{
		clipper::Coord_orth cp{ pt.m_x, pt.m_y, pt.m_z };
		clipper::Coord_frac pf = cp.coord_frac(xmap.cell());
		auto dp = xmap.interp<clipper::Interp_cubic>(pf);

		sumMass += dp;
		c[0] += dp * pt.m_x;
		c[1] += dp * pt.m_y;
		c[2] += dp * pt.m_z;
	}

	return { static_cast<float>(c[0] / sumMass), static_cast<float>(c[1] / sumMass), static_cast<float>(c[2] / sumMass) };
}

cif::point centerOfMassLigand(const cif::mm::residue &lig)
{
	double sumMass = 0;
	std::array<double, 3> c{ 0, 0, 0  };
	for (auto a : lig.atoms())
	{
		auto pt = a.get_location();
		auto dp = cif::atom_type_traits(a.get_type()).weight();

		sumMass += dp;
		c[0] += dp * pt.m_x;
		c[1] += dp * pt.m_y;
		c[2] += dp * pt.m_z;
	}

	return { static_cast<float>(c[0] / sumMass), static_cast<float>(c[1] / sumMass), static_cast<float>(c[2] / sumMass) };
}

double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap,
	const std::vector<cif::point> &blob)
{
	const auto dots = cif::spherical_dots<30>::instance();

	auto &ligand = structure.get_residue(asym_id);

	std::vector<cif::point> atomLocations;
	for (auto a : ligand.atoms())
		atomLocations.emplace_back(a.get_location());

	// Locate the center of the blob
	auto blobCenterOfMass = centerOfMassBlob(xmap, blob);

	// Same for the ligand
	auto ligandCenterOfMass = centerOfMassLigand(ligand);

	// Move ligand to the correct center and store new positions
	auto d = blobCenterOfMass - ligandCenterOfMass;
	for (auto li = atomLocations.begin(); auto a : ligand.atoms())
		a.set_location(*li++ += d);

	// Calculate inertia tensors for both ligand and blob

	auto itb = createInertiaTensorForBlob(blob, xmap);
	auto itl = createInertiaTensorForLigand(ligand);

	// Take principal vector, using eigen values

	auto ivb = principalAxis(itb);
	auto ivl = principalAxis(itl);

	// rotate ligand to match blob
	auto q = cif::construct_from_angle_axis(cif::angle(ivb, { 0, 0, 0 }, ivl), cif::cross_product(ivl, ivb));

	for (auto &loc : atomLocations)
		loc.rotate(q, blobCenterOfMass);

	// Test to see of a 180° rotation is needed
	auto [blobCenterOfSphere, blobRadius] = cif::smallest_sphere_around_points(blob);
	auto [ligandCenterOfSphere, ligandRadius] = cif::smallest_sphere_around_points(atomLocations);
	if (auto a = cif::angle(ligandCenterOfSphere, blobCenterOfMass, blobCenterOfSphere); a > 90)
	{
		q = cif::construct_from_angle_axis(180, cif::cross_product(ligandCenterOfSphere - blobCenterOfMass, blobCenterOfSphere - blobCenterOfMass));

		for (auto &loc : atomLocations)
			loc.rotate(q, blobCenterOfMass);
	}

	cif::crystal crystal(structure.get_datablock());

	double bestScore = 0;
	std::vector<cif::point> bestLoc;

	for (size_t i = 0; i < dots.size(); ++i)
	{
		auto axis = cif::cross_product(dots[0], dots[i]);
		auto angle = cif::angle(dots[0], {}, dots[i]);

		auto q = cif::construct_from_angle_axis(angle, axis); // NOLINT(bugprone-narrowing-conversions)

		for (auto li = atomLocations.begin(); auto a : ligand.atoms())
		{
			auto loc = *li++;
			loc.rotate(q, blobCenterOfMass);
			a.set_location(loc);
		}

		JiggleFitter f(ligand, xmap);
		auto jScore = f.refine();

		if (cif::VERBOSE > 1)
			std::cout << "jigglefit score: " << jScore << " for iteration " << i << "\n";

		if (jScore >= bestScore)
			continue;

		bestLoc.clear();
		for (auto a : ligand.atoms())
			bestLoc.emplace_back(a.get_location());
		bestScore = jScore;
	}

	for (auto li = bestLoc.begin(); auto a : ligand.atoms())
		a.set_location(*li++);

	return bestScore;
}

} // namespace pdb_redo