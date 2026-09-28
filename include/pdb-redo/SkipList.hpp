// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

/*
    Created by: Maarten L. Hekkelman
    Date: maandag 07 januari, 2019

    Skip list for e.g. pepflip
*/

#pragma once

#include <cif++/cif++.hpp>
#include <filesystem>
#include <optional>

namespace pdb_redo
{

struct ResidueSpec
{
	std::string pdb_asym_id;
	std::string pdb_comp_id;
	std::string pdb_seq_id;
	std::optional<char> pdbx_PDB_ins_code;
	std::string label_asym_id;
	std::string label_comp_id;
	int label_seq_id;

	ResidueSpec() = default;
	ResidueSpec(std::string auth_asym_id,
		std::string auth_comp_id,
		std::string auth_seq_id,
		const std::string &pdbx_PDB_ins_code,
		std::string label_asym_id,
		std::string label_comp_id,
		int label_seq_id)
		: pdb_asym_id(std::move(auth_asym_id))
		, pdb_comp_id(std::move(auth_comp_id))
		, pdb_seq_id(std::move(auth_seq_id))
		, pdbx_PDB_ins_code(pdbx_PDB_ins_code.empty() ? std::optional<char>{} : std::make_optional(pdbx_PDB_ins_code.c_str()[0]))
		, label_asym_id(std::move(label_asym_id))
		, label_comp_id(std::move(label_comp_id))
		, label_seq_id(label_seq_id)
	{
	}

	ResidueSpec(const ResidueSpec &rhs) = default;
	ResidueSpec &operator=(const ResidueSpec &rhs) = default;

	ResidueSpec(const cif::mm::residue &res)
		: pdb_asym_id(res.get_pdb_strand_id())
		, pdb_comp_id(res.get_compound_id())
		, pdb_seq_id(res.get_pdb_seq_num())
		, label_asym_id(res.get_asym_id())
		, label_comp_id(res.get_compound_id())
		, label_seq_id(res.get_seq_id())
	{
		char ins_code = res.get_pdb_ins_code().c_str()[0];
		if (ins_code != 0 and ins_code != ' ')
			pdbx_PDB_ins_code = ins_code;
	}

	ResidueSpec(const cif::mm::atom &atom)
		: pdb_asym_id(atom.get_auth_asym_id())
		, pdb_comp_id(atom.get_label_comp_id())
		, pdb_seq_id(atom.get_auth_seq_id())
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
		return pdb_asym_id == rhs.pdb_asym_id and
		       pdb_comp_id == rhs.pdb_comp_id and
		       pdb_seq_id == rhs.pdb_seq_id and
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