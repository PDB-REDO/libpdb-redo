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
    Date: maandag 07 januari, 2019

    Skip list for e.g. pepflip
*/

#pragma once

#include <filesystem>

#include <cif++.hpp>

#if __has_include(<experimental/optional>)
#include <experimental/optional>
using std::experimental::optional;
#else
#include <optional>
using std::optional;
#endif

namespace pdb_redo
{

struct ResidueSpec
{
	std::string auth_asym_id;
	std::string auth_comp_id;
	std::string auth_seq_id;
	optional<char> pdbx_PDB_ins_code;
	std::string label_asym_id;
	std::string label_comp_id;
	int label_seq_id;

	ResidueSpec() {}
	ResidueSpec(const std::string &auth_asym_id,
		const std::string &auth_comp_id,
		const std::string &auth_seq_id,
		const std::string &pdbx_PDB_ins_code,
		const std::string &label_asym_id,
		const std::string &label_comp_id,
		int label_seq_id)
		: auth_asym_id(auth_asym_id)
		, auth_comp_id(auth_comp_id)
		, auth_seq_id(auth_seq_id)
#if __has_include(<experimental/optional>)
		, pdbx_PDB_ins_code(pdbx_PDB_ins_code.empty() ? optional<char>{} : std::experimental::make_optional(pdbx_PDB_ins_code.c_str()[0]))
#else
		, pdbx_PDB_ins_code(pdbx_PDB_ins_code.empty() ? optional<char>{} : std::make_optional(pdbx_PDB_ins_code.c_str()[0]))
#endif
		, label_asym_id(label_asym_id)
		, label_comp_id(label_comp_id)
		, label_seq_id(label_seq_id)
	{
	}

	ResidueSpec(const ResidueSpec &rhs) = default;
	ResidueSpec &operator=(const ResidueSpec &rhs) = default;

	ResidueSpec(const cif::mm::residue &res)
		: auth_asym_id(res.get_auth_asym_id())
		, auth_comp_id(res.get_compound_id())
		, auth_seq_id(res.get_auth_seq_id())
		, label_asym_id(res.get_asym_id())
		, label_comp_id(res.get_compound_id())
		, label_seq_id(res.get_seq_id())
	{
		char ins_code = res.get_pdb_ins_code().c_str()[0];
		if (ins_code != 0 and ins_code != ' ')
			pdbx_PDB_ins_code = ins_code;
	}

	ResidueSpec(const cif::mm::atom &atom)
		: auth_asym_id(atom.get_auth_asym_id())
		, auth_comp_id(atom.get_label_comp_id())
		, auth_seq_id(atom.get_auth_seq_id())
		, label_asym_id(atom.get_label_asym_id())
		, label_comp_id(atom.get_label_comp_id())
		, label_seq_id(atom.get_label_seq_id())
	{
		char ins_code = atom.get_pdb_ins_code().c_str()[0];
		if (ins_code != 0 and ins_code != ' ')
			pdbx_PDB_ins_code = ins_code;
	}

	bool operator==(const ResidueSpec &rhs) const
	{
		return auth_asym_id == rhs.auth_asym_id and
		       auth_comp_id == rhs.auth_comp_id and
		       auth_seq_id == rhs.auth_seq_id and
		       pdbx_PDB_ins_code == rhs.pdbx_PDB_ins_code and
		       label_asym_id == rhs.label_asym_id and
		       label_comp_id == rhs.label_comp_id and
		       label_seq_id == rhs.label_seq_id;
	}
};

using SkipList = std::vector<ResidueSpec>;

enum class SkipListFormat
{
	OLD,
	CIF
};

// --------------------------------------------------------------------

SkipList readSkipList(const std::filesystem::path &file);
SkipList readSkipList(std::istream &is);

void writeSkipList(std::ostream &os, const SkipList &list, SkipListFormat format);

void writeSkipList(const std::filesystem::path &file, const SkipList &list, SkipListFormat format);

} // namespace pdb_redo