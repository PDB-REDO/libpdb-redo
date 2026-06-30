// Copyright NKI/AVL 2026
//
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cif++/cif++.hpp>

namespace pdb_redo
{

std::vector<cif::file> createConformers(const std::string inCompoundID, int N);

}