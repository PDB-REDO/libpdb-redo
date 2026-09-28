// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cif++/datablock.hpp>

#include <string>
#include <tuple>
#include <vector>

namespace pdb_redo
{

struct tls_selection;
struct tls_residue;

struct tls_selection
{
	virtual ~tls_selection() = default;
	virtual void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel = 0) const = 0;
	std::vector<std::tuple<std::string, int, int>> get_ranges(cif::datablock &db) const;
};

// Low level: get the selections
std::unique_ptr<tls_selection> parse_tls_selection_details(const std::string &program, const std::string &selection);

} // namespace cif
