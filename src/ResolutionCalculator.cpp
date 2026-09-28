// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

#include <cmath>
#include <numbers>

#include "pdb-redo/ResolutionCalculator.hpp"

namespace pdb_redo
{

const double kPI = std::numbers::pi;

ResolutionCalculator::ResolutionCalculator(const clipper::Cell &cell)
	: ResolutionCalculator(cell.a(), cell.b(), cell.c(),
		  180 * cell.alpha() / kPI,
		  180 * cell.beta() / kPI,
		  180 * cell.gamma() / kPI)
{
}

ResolutionCalculator::ResolutionCalculator(double a, double b, double c,
	double alpha, double beta, double gamma)
{
	double deg2rad = std::atan(1.0) / 45.0;

	double ca = std::cos(deg2rad * alpha);
	double sa = std::sin(deg2rad * alpha);
	double cb = std::cos(deg2rad * beta);
	double sb = std::sin(deg2rad * beta);
	double cg = std::cos(deg2rad * gamma);
	double sg = std::sin(deg2rad * gamma);

	double cast = (cb * cg - ca) / (sb * sg);
	double cbst = (cg * ca - cb) / (sg * sa);
	double cgst = (ca * cb - cg) / (sa * sb);

	double sast = std::sqrt(1 - cast * cast);
	double sbst = std::sqrt(1 - cbst * cbst);
	double sgst = std::sqrt(1 - cgst * cgst);

	double ast = 1 / (a * sb * sgst);
	double bst = 1 / (b * sg * sast);
	double cst = 1 / (c * sa * sbst);

	mCoefs[0] = ast * ast;
	mCoefs[1] = 2 * ast * bst * cgst;
	mCoefs[2] = 2 * ast * cst * cbst;
	mCoefs[3] = bst * bst;
	mCoefs[4] = 2 * bst * cst * cast;
	mCoefs[5] = cst * cst;
}

} // namespace pdb_redo
