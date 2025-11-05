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

#include "pdb-redo/Compound.hpp"
#include "pdb-redo/Minimizer.hpp"

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
#include <gsl/gsl_blas.h> // for debugging norm of gradient
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector_double.h>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>

namespace pdb_redo
{

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
				float min_z, max_z;
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

				result.erase(
					std::remove_if(result.begin(), result.end(), [x, y, min_z, max_z](const clipper::Coord_grid &p)
						{ return p[0] == x and p[1] == y and p[2] > min_z and p[2] < max_z; }),
					result.end());
			}
		}

		for (float x : vx)
		{
			for (float z : vz)
			{
				float min_y, max_y;
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

				result.erase(
					std::remove_if(result.begin(), result.end(), [x, z, min_y, max_y](const clipper::Coord_grid &p)
						{ return p[0] == x and p[2] == z and p[1] > min_y and p[1] < max_y; }),
					result.end());
			}
		}

		for (float y : vy)
		{
			for (float z : vz)
			{
				float min_x, max_x;
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

				result.erase(
					std::remove_if(result.begin(), result.end(), [y, z, min_x, max_x](const clipper::Coord_grid &p)
						{ return p[1] == y and p[2] == z and p[0] > min_x and p[0] < max_x; }),
					result.end());
			}
		}
	}

	return result;
}

// --------------------------------------------------------------------
// This code is highly inspired by the jiggle fit implemented in carbivore

class JiggleFitter
{
  public:
	JiggleFitter(pdb_redo::Minimizer *minimizer, cif::mm::residue &res);

	double refine();

	double score()
	{
		return m_minimizer->score();
	}

	virtual void transform();

  protected:
	pdb_redo::Minimizer *m_minimizer;

	std::vector<cif::mm::atom> m_atoms;
	std::vector<cif::point> m_locations;

	cif::point m_center;

	std::vector<double> m_variables;
	std::vector<float> m_stepsizes;
	const cif::mm::structure &m_structure;

	static double F(const gsl_vector *v, void *params)
	{
		JiggleFitter *self = reinterpret_cast<JiggleFitter *>(params);
		return self->F(v);
	}

	double F(const gsl_vector *v)
	{
		for (size_t i = 0; i < m_variables.size(); ++i)
			m_variables[i] = gsl_vector_get(v, i);

		transform();

		return score();
	}

	static double call(void *data, long n, const double *values)
	{
		return reinterpret_cast<JiggleFitter *>(data)->call(n, values);
	}

	double call(long n, const double *values)
	{
		// assert(n == m_variables.size());
		std::copy(values, values + n, m_variables.begin());
		transform();

		return score();
	}
};

JiggleFitter::JiggleFitter(pdb_redo::Minimizer *minimizer, cif::mm::residue &res)
	: m_minimizer(std::move(minimizer))
	, m_atoms(res.atoms())
	, m_structure(*res.get_structure())
{
	for (auto &atom : m_atoms)
		m_locations.emplace_back(atom.get_location());

	std::tie(m_center, std::ignore) = cif::smallest_sphere_around_points(m_locations);

	// This is where to set step sizes that are used during the rigid body fit
	m_variables = { 0, 0, 0, 0, 0 };
	m_stepsizes = { 1.f, 1.f, 0.2f, 0.2f, 0.2f };
}

void JiggleFitter::transform()
{
	// get the variables assigned
	const double alpha = m_variables[0], beta = m_variables[1], x = m_variables[2], y = m_variables[3], z = m_variables[4];

	// rotations
	auto q0 = cif::construct_from_angle_axis(alpha, { 1, 0, 0 }); // construct quaternion from float angle, point axis
	auto q1 = cif::construct_from_angle_axis(beta, { 0, 0, 1 });

	auto q = q0 * q1;

	// translations
	cif::point translation(x, y, z);

	// Move the molecule
	for (size_t i = 0; i < m_locations.size(); ++i)
	{
		auto a = m_locations[i];
		a.rotate(q, m_center); // rotate a using quaternion q01, move it to the alpha_loc - rotate - move back
		m_atoms[i].set_location(a + translation);
	}
}

double JiggleFitter::refine()
{
	const int kMaxIterations = 4000;

	gsl_multimin_function f = {
		.f = &JiggleFitter::F,
		.n = m_variables.size(),
		.params = this
	};

	auto T = gsl_multimin_fminimizer_nmsimplex2;
	auto x = gsl_vector_alloc(m_variables.size());

	for (size_t i = 0; i < m_variables.size(); ++i)
		gsl_vector_set(x, i, m_variables[i]);

	auto ss = gsl_vector_alloc(m_variables.size());
	for (size_t i = 0; i < m_stepsizes.size(); ++i)
		gsl_vector_set(ss, i, m_stepsizes[i]);

	auto s = gsl_multimin_fminimizer_alloc(T, m_variables.size());

	gsl_multimin_fminimizer_set(s, &f, x, ss);

	int iter = 0, status;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		double size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-2);

		if (status == GSL_SUCCESS)
		{
			if (cif::VERBOSE > 1)
				std::cout << "Minimum reached after " << iter << " iterations" << std::endl;
		}
	} while (status == GSL_CONTINUE and iter < kMaxIterations);

	for (size_t i = 0; i < m_variables.size(); ++i)
		m_variables[i] = gsl_vector_get(s->x, i);

	transform();

	gsl_vector_free(x);
	gsl_vector_free(ss);

	auto result = s->fval;

	gsl_multimin_fminimizer_free(s);

	// m_minimizer->printStats();

	return result;
}

// // --------------------------------------------------------------------

// double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap)
// {
// 	const auto dots = cif::spherical_dots<7>::instance();

// 	auto cellVolume = xmap.cell().volume();
// 	auto gridSize = xmap.grid_sampling().size();
// 	auto gridPointVolume = cellVolume / gridSize;
// 	auto gridPointRadius = std::pow((3 * gridPointVolume) / (4 * cif::kPI), 1 / 3.0);

// 	// Locate the center of the blob
// 	std::vector<cif::point> blob;
// 	for (auto p : findSingleBlob(xmap, false))
// 		blob.emplace_back(xmap.coord_orth(p.coord_map()));

// 	for (auto p : blob)
// 		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f}, r: {:.3f} }},\n", p.m_x, p.m_y, p.m_z, gridPointRadius);
// 	std::cout << "\n\n";

// 	auto [blobCenter, blobRadius] = cif::smallest_sphere_around_points(blob);

// 	for (auto p : blob)
// 	{
// 		auto d = cif::distance(p, blobCenter);
// 		assert(d < 1.01f * blobRadius);
// 	}

// 	// Same for the ligand
// 	auto &ligand = structure.get_residue(asym_id);

// 	std::vector<cif::point> atomLocations;
// 	for (auto a : ligand.atoms())
// 		atomLocations.emplace_back(a.get_location());
// 	auto [ligandCenter, ligandRadius] = cif::smallest_sphere_around_points(atomLocations);

// 	for (auto p : atomLocations)
// 	{
// 		auto d = cif::distance(p, ligandCenter);
// 		assert(d < 1.01f * ligandRadius);
// 	}

// 	// Move ligand to the correct center
// 	if (cif::distance(ligandCenter, blobCenter) > 0.1f)
// 	{
// 		auto d = blobCenter - ligandCenter;
// 		atomLocations.clear();

// 		for (auto a : ligand.atoms())
// 		{
// 			auto loc = a.get_location() + d;
// 			a.set_location(loc);
// 			atomLocations.emplace_back(loc);
// 		}
// 	}

// 	for (size_t ix = 0; auto p : atomLocations)
// 		std::cout << std::format("{{ x: {:.4f}, y: {:.4f}, z: {:.4f}, r: {:.3f} }},\n", p.m_x, p.m_y, p.m_z, cif::atom_type_traits(ligand.atoms()[ix++].get_type()).radius());
// 	std::cout << "\n\n";

// 	struct Score
// 	{
// 		cif::quaternion q;
// 		double v;
// 		std::vector<cif::point> loc;

// 		bool operator<(const Score &rhs) const
// 		{
// 			return v < rhs.v;
// 		}
// 	};

// 	std::vector<Score> best;
// 	cif::crystal crystal(structure.get_datablock());

// 	std::unique_ptr<pdb_redo::Minimizer> minimizer(pdb_redo::Minimizer::create(crystal, structure, ligand.atoms(), xmap));

// 	for (size_t i = 0; i < dots.size(); ++i)
// 	{
// 		auto axis = cif::cross_product(dots[0], dots[i]);
// 		auto angle = cif::angle(dots[0], {}, dots[i]);

// 		auto q = cif::construct_from_angle_axis(angle, axis);

// 		for (auto li = atomLocations.begin(); auto a : ligand.atoms())
// 		{
// 			auto loc = *li++;
// 			loc.rotate(q, blobCenter);
// 			a.set_location(loc);
// 		}

// 		JiggleFitter f(minimizer.get(), ligand);
// 		auto jScore = f.refine();

// 		std::cout << "jigglefit score: " << jScore << " for iteration " << i << "\n";

// 		if (jScore > 0)
// 			continue;

// 		auto score = minimizer->refine(false);
// 		std::cout << "score: " << score << " for iteration " << i << "\n";

// 		// minimizer->printStats();

// 		std::vector<cif::point> bestLoc;
// 		for (auto a : ligand.atoms())
// 			bestLoc.emplace_back(a.get_location());

// 		best.emplace_back(q, score, std::move(bestLoc));

// 		std::push_heap(best.begin(), best.end());
// 	}

// 	std::sort_heap(best.begin(), best.end());

// 	for (bool first = true; auto [q, v, loc] : best)
// 	{
// 		std::cout << "q: " << q << ", v: " << v << "\n";
// 		if (std::exchange(first, false))
// 		{
// 			for (size_t ix = 0; ix < loc.size(); ++ix)
// 				ligand.atoms()[ix].set_location(loc[ix]);
// 		}
// 	}

// 	return best.front().v;
// }

// --------------------------------------------------------------------

class ConformationIterator
{
  public:
	ConformationIterator(cif::mm::structure &structure, const std::string &asym_id);

	void next();
	bool last() const;

  private:
	struct TorsionData
	{
		cif::mm::atom dihedral_atoms[4];
		std::vector<cif::mm::atom> rotating_atoms;
	};

	struct Iteration
	{
		size_t index;
		float angle;
	};

	cif::mm::residue &m_residue;
	std::vector<TorsionData> m_torsions;
};

ConformationIterator::ConformationIterator(cif::mm::structure &structure, const std::string &asym_id)
	: m_residue(structure.get_residue(asym_id))
{
	auto compound = pdb_redo::CompoundFactory::instance().create(m_residue.get_compound_id());

	for (auto &torsion : compound->torsions())
	{
		if (torsion.period <= 1)
			continue;

		TorsionData td{
			{ m_residue.get_atom_by_atom_id(torsion.atomID[0]),
				m_residue.get_atom_by_atom_id(torsion.atomID[1]),
				m_residue.get_atom_by_atom_id(torsion.atomID[2]),
				m_residue.get_atom_by_atom_id(torsion.atomID[3]) }
		};

		bool ok = true;
		for (auto &a : td.dihedral_atoms)
			ok = ok and a and a.get_type() != cif::H;
		if (not ok)
			continue;
	}
}

void ConformationIterator::next()
{
}

bool ConformationIterator::last() const
{
}

// --------------------------------------------------------------------

double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap,
	const std::vector<cif::point> &blob)
{
	const auto dots = cif::spherical_dots<7>::instance();

	// Locate the center of the blob
	auto [blobCenter, blobRadius] = cif::smallest_sphere_around_points(blob);

	// Same for the ligand
	auto &ligand = structure.get_residue(asym_id);

	std::vector<cif::point> atomLocations;
	for (auto a : ligand.atoms())
		atomLocations.emplace_back(a.get_location());
	auto [ligandCenter, ligandRadius] = cif::smallest_sphere_around_points(atomLocations);

	// Move ligand to the correct center and store new positions
	auto d = blobCenter - ligandCenter;
	atomLocations.clear();

	for (auto a : ligand.atoms())
	{
		auto loc = a.get_location() + d;
		a.set_location(loc);
		atomLocations.emplace_back(loc);
	}

	struct Score
	{
		cif::quaternion q;
		double v;
		std::vector<cif::point> loc;

		bool operator<(const Score &rhs) const
		{
			return v < rhs.v;
		}
	};

	std::vector<Score> best;
	cif::crystal crystal(structure.get_datablock());
	std::unique_ptr<Minimizer> minimizer(Minimizer::create(crystal, structure, ligand.atoms(), xmap));

	for (size_t i = 0; i < dots.size(); ++i)
	{
		auto axis = cif::cross_product(dots[0], dots[i]);
		auto angle = cif::angle(dots[0], {}, dots[i]);

		auto q = cif::construct_from_angle_axis(angle, axis);

		for (auto li = atomLocations.begin(); auto a : ligand.atoms())
		{
			auto loc = *li++;
			loc.rotate(q, blobCenter);
			a.set_location(loc);
		}

		JiggleFitter f(minimizer.get(), ligand);
		auto jScore = f.refine();

		if (cif::VERBOSE > 1)
			std::cout << "jigglefit score: " << jScore << " for iteration " << i << "\n";

		if (jScore > 0)
			continue;

		auto score = minimizer->refine(false);
		if (cif::VERBOSE > 1)
			std::cout << "score: " << score << " for iteration " << i << "\n";

		// minimizer->printStats();

		std::vector<cif::point> bestLoc;
		for (auto a : ligand.atoms())
			bestLoc.emplace_back(a.get_location());

		best.emplace_back(q, score, std::move(bestLoc));

		std::push_heap(best.begin(), best.end());
	}

	std::sort_heap(best.begin(), best.end());

	for (bool first = true; auto [q, v, loc] : best)
	{
		if (cif::VERBOSE > 1)
			std::cout << "q: " << q << ", v: " << v << "\n";

		if (std::exchange(first, false))
		{
			for (size_t ix = 0; ix < loc.size(); ++ix)
				ligand.atoms()[ix].set_location(loc[ix]);
		}
	}

	return best.empty() ? 0 : best.front().v;
}

} // namespace pdb_redo