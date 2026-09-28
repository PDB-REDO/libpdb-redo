// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

// The only reason for this source file is to make sure the version info for
// libpdb-redo is available in client software.

#include "pdb-redo/Version.hpp"

#include "revision.hpp"

#include <sstream>

namespace pdb_redo
{

PDB_REDO_EXPORT int force_link = 1;

std::string get_version()
{
	std::ostringstream os;
	write_version_string(os, false);
	return os.str();
}

}
