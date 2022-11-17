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
*/

#include <filesystem>
#include <fstream>

#include <zeep/xml/document.hpp>
#include <zeep/xml/serialize.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>
#include <zeep/json/serializer.hpp>

#include "pdb-redo/SkipList.hpp"

namespace fs = std::filesystem;
namespace zx = zeep::xml;

namespace pdb_redo
{

// --------------------------------------------------------------------

void writeOLDSkipList(std::ostream &os, const SkipList &list)
{
	os << ':';
	for (auto &res : list)
		os << res.auth_asym_id << res.auth_seq_id << (res.pdbx_PDB_ins_code and *res.pdbx_PDB_ins_code ? *res.pdbx_PDB_ins_code : ' ') << ':';
}

void writeJSONSkipList(std::ostream &os, const SkipList &list)
{
	zeep::json::element e;
	zeep::json::to_element(e, list);
	os << e;
}

void writeCIFSkipList(std::ostream &os, const SkipList &list)
{
	cif::file file;
	auto &db = file.emplace_back("skip");

	auto &&[cat, ignore] = db.emplace("skip_list");

	for (auto &res : list)
		cat->emplace({{"auth_asym_id", res.auth_asym_id},
			{"auth_comp_id", res.auth_comp_id},
			{"auth_seq_id", res.auth_seq_id},
			{"pdbx_PDB_ins_code", std::string{res.pdbx_PDB_ins_code and *res.pdbx_PDB_ins_code ? *res.pdbx_PDB_ins_code : '?'}},
			{"label_asym_id", res.label_asym_id},
			{"label_comp_id", res.label_comp_id},
			{"label_seq_id", res.label_seq_id}});

	file.save(os);
}

// --------------------------------------------------------------------

void writeSkipList(std::ostream &os, const SkipList &list, SkipListFormat format)
{
	switch (format)
	{
		case SkipListFormat::OLD:
			writeOLDSkipList(os, list);
			break;

		case SkipListFormat::JSON:
			writeJSONSkipList(os, list);
			break;

		case SkipListFormat::CIF:
			writeCIFSkipList(os, list);
			break;

		default:
			break;
	}
}

void writeSkipList(const fs::path &file, const SkipList &list, SkipListFormat format)
{
	std::ofstream os(file, std::ios::binary);
	writeSkipList(os, list, format);
}

// --------------------------------------------------------------------

SkipList readOLDSkipList(std::istream &is)
{
	SkipList result;

	char separator = 0;
	if (is.rdbuf()->in_avail() > 0)
		is.read(&separator, 1);
	if (separator != ':')
		throw std::runtime_error("Not an old format skip list");

	while (not is.eof())
	{
		ResidueSpec spec{};

		char chain = 0;

		is.read(&chain, 1);

		int seq_nr;
		is >> seq_nr;

		char ins_code = 0;
		is.read(&ins_code, 1);

		if (is.rdbuf()->in_avail() > 0)
			is.read(&separator, 1);

		if (chain == 0)
			break;

		if (separator != ':' and separator != 0)
			throw std::runtime_error("Invalid old format skiplist");

		spec.auth_asym_id.push_back(chain);
		spec.auth_seq_id = std::to_string(seq_nr);
		if (ins_code != ' ' and ins_code != 0)
			spec.pdbx_PDB_ins_code = ins_code;

		result.push_back(spec);
	}

	return result;
}

SkipList readJSONSkipList(std::istream &is)
{
	SkipList result;

	zeep::json::element e;
	zeep::json::parse_json(is, e);

	from_element(e, result);

	return result;
}

SkipList readCIFSkipList(std::istream &is)
{
	SkipList result;

	cif::file file;
	file.load(is);

	auto &db = file["skip"];
	auto &cat = db["skip_list"];

	for (const auto &[auth_asym_id, auth_comp_id, auth_seq_id, pdbx_PDB_ins_code, label_asym_id, label_comp_id, label_seq_id] :
		cat.rows<std::string, std::string, std::string, std::string, std::string, std::string, int>("auth_asym_id", "auth_comp_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_asym_id", "label_comp_id", "label_seq_id"))
	{
		result.emplace_back(auth_asym_id, auth_comp_id, auth_seq_id, pdbx_PDB_ins_code, label_asym_id, label_comp_id, label_seq_id);
	}

	return result;
}

SkipList readSkipList(std::istream &is)
{
	try
	{
		return readCIFSkipList(is);
	}
	catch (const std::exception &e)
	{
		if (cif::VERBOSE > 0)
			std::cerr << e.what() << std::endl;
		is.rdbuf()->pubseekpos(0);
	}

	try
	{
		return readJSONSkipList(is);
	}
	catch (const std::exception &e)
	{
		if (cif::VERBOSE > 0)
			std::cerr << e.what() << std::endl;
		is.rdbuf()->pubseekpos(0);
	}

	return readOLDSkipList(is);
}

SkipList readSkipList(std::filesystem::path &file)
{
	std::ifstream is(file);
	return readSkipList(is);
}

} // namespace pdb_redo