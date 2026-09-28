// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2025
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cif++/model.hpp>
#include <cif++/matrix.hpp>

#include <clipper/clipper.h>

namespace pdb_redo
{

double fitShape(cif::mm::structure &structure, const std::string &asym_id, clipper::Xmap<float> &xmap,
    const std::vector<cif::point> &blob);

// might be useful outside
cif::matrix3x3<float> createInertiaTensorForBlob(const std::vector<cif::point> pts, clipper::Xmap<float> &xmap);
cif::matrix3x3<float> createInertiaTensorForLigand(const cif::mm::residue &res);
cif::point principalAxis(const cif::matrix3x3<float> &m);

std::vector<cif::point> create_spherical_dots(int P);

}