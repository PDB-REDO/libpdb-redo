// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

#include "pdb-redo/exports.hpp"

#include <string>

namespace pdb_redo
{

// To force link the version code in pdb-redo, assign a value
// to the following global variable somewhere in your code.
extern PDB_REDO_EXPORT int force_link;

std::string get_version();

}