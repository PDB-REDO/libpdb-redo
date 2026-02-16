/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2017 NKI/AVL, Netherlands Cancer Institute
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

#include "pdb-redo/TLS.hpp"

#include <algorithm>
#include <cif++.hpp>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

namespace pdb_redo
{

const int
	kResidueNrWildcard = std::numeric_limits<int>::min(),
	kNoSeqNum = std::numeric_limits<int>::max() - 1;

// --------------------------------------------------------------------
// We parse selection statements and create a selection expression tree
// which is then interpreted by setting the selected flag for the
// residues. After that, the selected ranges are collected and printed.

struct tls_residue
{
	std::string chainID;
	int seqNr = 0;
	char iCode;
	std::string name;
	bool selected;

	// std::string asymID;
	// int seqID = 0;

	bool operator==(const tls_residue &rhs) const
	{
		return chainID == rhs.chainID and
		       seqNr == rhs.seqNr and
		       iCode == rhs.iCode and
		       cif::iequals(name, rhs.name) and
		       selected == rhs.selected;
	}
};

void dump_selection(const std::vector<tls_residue> &selected, std::size_t indentLevel)
{
	std::string indent(indentLevel * 2, ' ');

	auto i = selected.begin();
	bool first = true;

	// First print in PDB space
	while (i != selected.end())
	{
		auto b = find_if(i, selected.end(), [](auto s) -> bool
			{ return s.selected; });
		if (b == selected.end())
			break;

		if (first)
			std::cout << indent << "PDB:\n";
		first = false;

		auto e = find_if(b, selected.end(), [b](auto s) -> bool
			{ return s.chainID != b->chainID or not s.selected; });

		std::cout << indent << " >> " << b->chainID << ' ' << b->seqNr << ':' << (e - 1)->seqNr << '\n';
		i = e;
	}

	// Then in mmCIF space

	if (not first)
		std::cout << indent << "mmCIF:\n";

	i = selected.begin();
	while (i != selected.end())
	{
		auto b = find_if(i, selected.end(), [](auto s) -> bool
			{ return s.selected; });
		if (b == selected.end())
			break;

		// auto e = find_if(b, selected.end(), [b](auto s) -> bool
		// 	{ return s.asymID != b->asymID or not s.selected; });

		// std::string asymID = b->asymID;
		// int from = b->seqID, to = from;

		// for (auto j = b + 1; j != e; ++j)
		// {
		// 	if (j->seqID == to + 1)
		// 		to = j->seqID;
		// 	else if (j->seqID != to) // probably an insertion code
		// 	{
		// 		if (from == kNoSeqNum or to == kNoSeqNum)
		// 			std::cout << indent << " >> " << asymID << '\n';
		// 		else
		// 			std::cout << indent << " >> " << asymID << ' ' << from << ':' << to << '\n';
		// 		asymID = b->asymID;
		// 		from = to = b->seqID;
		// 	}
		// }

		auto e = find_if(b, selected.end(), [b](auto s) -> bool
			{ return s.chainID != b->chainID or not s.selected; });

		std::string chainID = b->chainID;
		int from = b->seqNr, to = from;

		for (auto j = b + 1; j != e; ++j)
		{
			if (j->seqNr == to + 1)
				to = j->seqNr;
			else if (j->seqNr != to) // probably an insertion code
			{
				if (from == kNoSeqNum or to == kNoSeqNum)
					std::cout << indent << " >> " << chainID << '\n';
				else
					std::cout << indent << " >> " << chainID << ' ' << from << ':' << to << '\n';
				chainID = b->chainID;
				from = to = b->seqNr;
			}
		}


		if (from == kNoSeqNum or to == kNoSeqNum)
			std::cout << indent << " >> " << chainID << '\n';
		else
			std::cout << indent << " >> " << chainID << ' ' << from << ':' << to << '\n';

		i = e;
	}

	if (first)
	{
		using namespace cif::colour;
		std::cout << indent << cif::coloured("Empty selection", white, red, bold) << '\n';
	}
}

std::vector<std::tuple<std::string, int, int>> tls_selection::get_ranges(cif::datablock &db) const
{
	std::vector<tls_residue> selected;

	// Collect the residues from poly seq scheme...
	for (const auto &[chain, seqNr, iCode, name] :
		db["pdbx_poly_seq_scheme"].rows<std::string,int,std::string,std::string>("pdb_strand_id", "pdb_seq_num", "pdb_ins_code", "pdb_mon_id"))
	{
		if (iCode.length() > 1)
			throw std::runtime_error("invalid iCode");

		selected.emplace_back(chain, seqNr, iCode[0], name);
	}

	// ... those from the nonpoly scheme
	for (const auto &[chain, iCode, name] :
		db["pdbx_nonpoly_scheme"].rows<std::string,std::string,std::string>("pdb_strand_id", "pdb_ins_code", "pdb_mon_id"))
	{
		if (cif::iequals(name, "HOH") or cif::iequals(name, "H2O"))
			continue;

		if (iCode.length() > 1)
			throw std::runtime_error("invalid iCode");

		selected.emplace_back(chain, 0, iCode[0], name);
	}

	// ... those from the nonpoly scheme
	for (const auto &[chain, iCode, name] :
		db["pdbx_branch_scheme"].rows<std::string,std::string,std::string>("pdb_strand_id", "pdb_ins_code", "pdb_mon_id"))
	{
		if (iCode.length() > 1)
			throw std::runtime_error("invalid iCode");

		selected.emplace_back(chain, 0, iCode[0], name);
	}

	// selected might consist of multiple ranges
	// output per chain

	std::ranges::stable_sort(selected, [](auto &a, auto &b) -> bool
		{
			int d = a.chainID.compare(b.chainID);
			if (d == 0)
				d = a.seqNr - b.seqNr;
			return d < 0; });

	collect_residues(db, selected);

	std::vector<std::tuple<std::string, int, int>> result;

	auto i = selected.begin();

	while (i != selected.end())
	{
		auto b = find_if(i, selected.end(), [](auto s) -> bool
			{ return s.selected; });
		if (b == selected.end())
			break;

		auto e = find_if(b, selected.end(), [b](auto s) -> bool
			{ return s.chainID != b->chainID or not s.selected; });

		// return ranges with strict increasing sequence numbers.
		// So when there's a gap in the sequence we split the range.
		// Beware of iCodes though
		result.emplace_back(b->chainID, b->seqNr, b->seqNr);
		for (auto j = b + 1; j != e; ++j)
		{
			if (j->seqNr == std::get<2>(result.back()) + 1)
				std::get<2>(result.back()) = j->seqNr;
			else if (j->seqNr != std::get<2>(result.back())) // probably an insertion code
				result.emplace_back(b->chainID, j->seqNr, j->seqNr);
		}

		i = e;
	}

	for (auto &&[name, i1, i2] : result)
	{
		if (i1 == kNoSeqNum) i1 = 0;
		if (i2 == kNoSeqNum) i2 = 0;
	}

	return result;
}

struct tls_selection_not : public tls_selection
{
	tls_selection_not(std::unique_ptr<tls_selection> selection)
		: selection(selection.release())
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		selection->collect_residues(db, residues, indentLevel + 1);

		for (auto &r : residues)
			r.selected = not r.selected;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "NOT\n";
			dump_selection(residues, indentLevel);
		}
	}

	std::unique_ptr<tls_selection> selection;
};

struct tls_selection_all : public tls_selection
{
	tls_selection_all() = default;

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		for (auto &r : residues)
			r.selected = true;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "ALL\n";
			dump_selection(residues, indentLevel);
		}
	}
};

struct tls_selection_chain : public tls_selection_all
{
	tls_selection_chain(std::string chainID)
		: m_chain(std::move(chainID))
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		bool allChains = m_chain == "*";

		for (auto &r : residues)
			r.selected = allChains or r.chainID == m_chain;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "CHAIN " << m_chain << '\n';
			dump_selection(residues, indentLevel);
		}
	}

	std::string m_chain;
};

struct tls_selection_res_id : public tls_selection_all
{
	tls_selection_res_id(int seqNr, char iCode)
		: m_seq_nr(seqNr)
		, m_icode(iCode)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		for (auto &r : residues)
			r.selected = r.seqNr == m_seq_nr and r.iCode == m_icode;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "ResID " << m_seq_nr << (m_icode ? std::string{ m_icode } : "") << '\n';
			dump_selection(residues, indentLevel);
		}
	}

	int m_seq_nr;
	char m_icode;
};

struct tls_selection_range_seq : public tls_selection_all
{
	tls_selection_range_seq(int first, int last)
		: m_first(first)
		, m_last(last)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		for (auto &r : residues)
		{
			r.selected = ((r.seqNr >= m_first or m_first == kResidueNrWildcard) and
						  (r.seqNr <= m_last or m_last == kResidueNrWildcard));
		}

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Range " << m_first << ':' << m_last << '\n';
			dump_selection(residues, indentLevel);
		}
	}

	int m_first, m_last;
};

struct tls_selection_range_id : public tls_selection_all
{
	tls_selection_range_id(int first, int last, char icodeFirst = 0, char icodeLast = 0)
		: m_first(first)
		, m_last(last)
		, m_icode_first(icodeFirst)
		, m_icode_last(icodeLast)
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		// need to do this per chain
		std::set<std::string> chains;
		for (auto &r : residues)
			chains.insert(r.chainID);

		for (std::string chain : chains)
		{
			auto f = std::ranges::find_if(residues,
				[this,chain](auto r) -> bool
				{
					return r.chainID == chain and r.seqNr == m_first and r.iCode == m_icode_first;
				});

			auto l = std::ranges::find_if(residues,
				[this,chain](auto r) -> bool
				{
					return r.chainID == chain and r.seqNr == m_last and r.iCode == m_icode_last;
				});

			if (f != residues.end() and l != residues.end() and f <= l)
			{
				++l;

				for (; f != l; ++f)
					f->selected = true;
			}
		}

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Through " << m_first << ':' << m_last << '\n';
			dump_selection(residues, indentLevel);
		}
	}

	int m_first, m_last;
	char m_icode_first, m_icode_last;
};

struct tls_selection_union : public tls_selection
{
	tls_selection_union(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	tls_selection_union(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &&rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		auto a = residues;
		std::ranges::for_each(a, [](auto &r)
			{ r.selected = false; });

		auto b = residues;
		std::ranges::for_each(b, [](auto &r)
			{ r.selected = false; });

		lhs->collect_residues(db, a, indentLevel + 1);
		rhs->collect_residues(db, b, indentLevel + 1);

		for (auto ai = a.begin(), bi = b.begin(), ri = residues.begin(); ri != residues.end(); ++ai, ++bi, ++ri)
			ri->selected = ai->selected or bi->selected;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Union\n";
			dump_selection(residues, indentLevel);
		}
	}

	std::unique_ptr<tls_selection> lhs;
	std::unique_ptr<tls_selection> rhs;
};

struct tls_selection_intersection : public tls_selection
{
	tls_selection_intersection(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	tls_selection_intersection(std::unique_ptr<tls_selection> &lhs, std::unique_ptr<tls_selection> &&rhs)
		: lhs(lhs.release())
		, rhs(rhs.release())
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		auto a = residues;
		std::ranges::for_each(a, [](auto &r)
			{ r.selected = false; });

		auto b = residues;
		std::ranges::for_each(b, [](auto &r)
			{ r.selected = false; });

		lhs->collect_residues(db, a, indentLevel + 1);
		rhs->collect_residues(db, b, indentLevel + 1);

		for (auto ai = a.begin(), bi = b.begin(), ri = residues.begin(); ri != residues.end(); ++ai, ++bi, ++ri)
			ri->selected = ai->selected and bi->selected;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Intersection\n";
			dump_selection(residues, indentLevel);
		}
	}

	std::unique_ptr<tls_selection> lhs;
	std::unique_ptr<tls_selection> rhs;
};

struct tls_selection_by_name : public tls_selection_all
{
  public:
	tls_selection_by_name(std::string resname)
		: m_name(std::move(resname))
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		for (auto &r : residues)
			r.selected = r.name == m_name;

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Name " << m_name << '\n';
			dump_selection(residues, indentLevel);
		}
	}

	std::string m_name;
};

struct tls_selection_by_element : public tls_selection_all
{
  public:
	tls_selection_by_element(std::string element)
		: m_element(std::move(element))
	{
	}

	void collect_residues(cif::datablock &db, std::vector<tls_residue> &residues, std::size_t indentLevel) const override
	{
		// rationale... We want to select residues only. So we select
		// residues that have just a single atom of type m_element.
		// And we assume these have as residue name... m_element.
		// ... Right?

		for (auto &r : residues)
			r.selected = cif::iequals(r.name, m_element);

		if (cif::VERBOSE > 0)
		{
			std::cout << std::string(indentLevel * 2, ' ') << "Element " << m_element << '\n';
			dump_selection(residues, indentLevel);
		}
	}

	std::string m_element;
};

// --------------------------------------------------------------------

class tls_selection_parser_impl
{
  public:
	tls_selection_parser_impl(std::string selection)
		: m_selection(std::move(selection))
		, m_p(m_selection.begin())
		, m_end(m_selection.end())
	{
	}

	virtual std::unique_ptr<tls_selection> Parse() = 0;

  protected:
	virtual int get_next_token() = 0;
	virtual void match(int token);
	virtual std::string to_string(int token) = 0;

	std::string m_selection;
	std::string::iterator m_p, m_end;
	int m_lookahead;
	std::string m_token;
};

void tls_selection_parser_impl::match(int token)
{
	if (m_lookahead == token)
		m_lookahead = get_next_token();
	else
	{
		std::string expected;
		if (token >= 256)
			expected = to_string(token);
		else
			expected = { static_cast<char>(token) };

		std::string found;
		if (m_lookahead >= 256)
			found = to_string(m_lookahead) + " (" + m_token + ')';
		else
			found = { static_cast<char>(m_lookahead) };

		throw std::runtime_error("Expected " + expected + " but found " + found);
	}
}

// --------------------------------------------------------------------

class TLSSelectionParserImplPhenix : public tls_selection_parser_impl
{
  public:
	TLSSelectionParserImplPhenix(const std::string &selection)
		: tls_selection_parser_impl(selection)
	{
		m_lookahead = get_next_token();
	}

	std::unique_ptr<tls_selection> Parse() override;

  private:
	std::unique_ptr<tls_selection> ParseAtomSelection();
	std::unique_ptr<tls_selection> ParseTerm();
	std::unique_ptr<tls_selection> ParseFactor();

	enum TOKEN
	{
		pt_NONE = 0,
		pt_IDENT = 256,
		pt_STRING = 257,
		pt_NUMBER = 258,
		pt_RESID = 259,
		pt_EOLN = 260,
		pt_KW_ALL = 261,
		pt_KW_CHAIN = 262,
		pt_KW_RESSEQ = 263,
		pt_KW_RESID = 264,
		pt_KW_ICODE = 265,
		pt_KW_RESNAME = 266,
		pt_KW_ELEMENT = 267,
		pt_KW_AND = 268,
		pt_KW_OR = 269,
		pt_KW_NOT = 270,
		pt_KW_PDB = 271,
		pt_KW_ENTRY = 272,
		pt_KW_THROUGH = 273
	};

	int get_next_token() override;
	std::string to_string(int token) override;

	int m_value_i;
	std::string m_value_s;
	char m_icode;
};

int TLSSelectionParserImplPhenix::get_next_token()
{
	int result = pt_NONE;
	enum STATE
	{
		st_START = 0,
		st_RESID = 200,
		st_NUM = 300,
		st_IDENT = 400,
		st_QUOTED = 500,
		st_DQUOTED = 550,
		st_OTHER = 600
	};
	int state = st_START;

	m_value_i = 0;
	m_icode = 0;
	m_value_s.clear();
	auto s = m_p;

	auto start = state;
	m_token.clear();

	auto restart = [&]()
	{
		switch (start)
		{
			case st_START: state = start = st_RESID; break;
			case st_RESID: state = start = st_NUM; break;
			case st_NUM: state = start = st_IDENT; break;
			case st_IDENT: state = start = st_QUOTED; break;
			case st_QUOTED: state = start = st_DQUOTED; break;
			case st_DQUOTED: state = start = st_OTHER; break;
			default:;
		}
		m_token.clear();
		m_p = s;
	};

	auto retract = [&]()
	{
		--m_p;
		m_token.pop_back();
	};

	while (result == pt_NONE)
	{
		char ch = *m_p++;
		if (m_p > m_end)
			ch = 0;
		else
			m_token += ch;

		switch (state)
		{
			// start block
			case st_START:
				if (ch == 0)
					result = pt_EOLN;
				else if (isspace(ch))
				{
					m_token.clear();
					++s;
				}
				else
					restart();
				break;

			// RESID block
			case st_RESID:
				if (ch == '-')
					state = st_RESID + 1;
				else if (isdigit(ch))
				{
					m_value_i = (ch - '0');
					state = st_RESID + 2;
				}
				else
					restart();
				break;

			case st_RESID + 1:
				if (isdigit(ch))
				{
					m_value_i = -(ch - '0');
					state = st_RESID + 2;
				}
				else
					restart();
				break;

			case st_RESID + 2:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (m_value_i < 0 ? -1 : 1) * (ch - '0');
				else if (isalpha(ch))
				{
					m_icode = ch;
					state = st_RESID + 3;
				}
				else
					restart();
				break;

			case st_RESID + 3:
				if (isalnum(ch))
					restart();
				else
				{
					retract();
					result = pt_RESID;
				}
				break;

				// NUM block

			case st_NUM:
				if (ch == '-')
					state = st_NUM + 1;
				else if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM + 2;
				}
				else
					restart();
				break;

			case st_NUM + 1:
				if (isdigit(ch))
				{
					m_value_i = -(ch - '0');
					state = st_NUM + 2;
				}
				else
					restart();
				break;

			case st_NUM + 2:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (m_value_i < 0 ? -1 : 1) * (ch - '0');
				else if (not isalpha(ch))
				{
					result = pt_NUMBER;
					retract();
				}
				else
					restart();
				break;

				// IDENT block

			case st_IDENT:
				if (isalnum(ch))
				{
					m_value_s = { ch };
					state = st_IDENT + 1;
				}
				else
					restart();
				break;

			case st_IDENT + 1:
				if (isalnum(ch) or ch == '\'')
					m_value_s += ch;
				else
				{
					--m_p;
					result = pt_IDENT;
				}
				break;

				// QUOTED block

			case st_QUOTED:
				if (ch == '\'')
				{
					m_value_s.clear();
					state = st_QUOTED + 1;
				}
				else
					restart();
				break;

			case st_QUOTED + 1:
				if (ch == '\'')
					result = pt_STRING;
				else if (ch == 0)
					throw std::runtime_error("Unexpected end of selection, missing quote character?");
				else
					m_value_s += ch;
				break;

				// QUOTED block

			case st_DQUOTED:
				if (ch == '\"')
				{
					m_value_s.clear();
					state = st_DQUOTED + 1;
				}
				else
					restart();
				break;

			case st_DQUOTED + 1:
				if (ch == '\"')
					result = pt_STRING;
				else if (ch == 0)
					throw std::runtime_error("Unexpected end of selection, missing quote character?");
				else
					m_value_s += ch;
				break;

			// OTHER block
			case st_OTHER:
				result = static_cast<unsigned char>(ch);
				break;

			default:
				throw std::runtime_error("invalid state in TLS parser");
		}
	}

	if (result == pt_IDENT)
	{
		if (cif::iequals(m_value_s, "CHAIN"))
			result = pt_KW_CHAIN;
		else if (cif::iequals(m_value_s, "ALL"))
			result = pt_KW_ALL;
		else if (cif::iequals(m_value_s, "AND"))
			result = pt_KW_AND;
		else if (cif::iequals(m_value_s, "OR"))
			result = pt_KW_OR;
		else if (cif::iequals(m_value_s, "NOT"))
			result = pt_KW_NOT;
		else if (cif::iequals(m_value_s, "RESSEQ"))
			result = pt_KW_RESSEQ;
		else if (cif::iequals(m_value_s, "RESID") or cif::iequals(m_value_s, "RESI"))
			result = pt_KW_RESID;
		else if (cif::iequals(m_value_s, "RESNAME"))
			result = pt_KW_RESNAME;
		else if (cif::iequals(m_value_s, "ELEMENT"))
			result = pt_KW_ELEMENT;
		else if (cif::iequals(m_value_s, "PDB"))
			result = pt_KW_PDB;
		else if (cif::iequals(m_value_s, "ENTRY"))
			result = pt_KW_ENTRY;
		else if (cif::iequals(m_value_s, "THROUGH"))
			result = pt_KW_THROUGH;
	}

	return result;
}

std::string TLSSelectionParserImplPhenix::to_string(int token)
{
	switch (token)
	{
		case pt_IDENT: return "identifier";
		case pt_STRING: return "std::string";
		case pt_NUMBER: return "number";
		case pt_RESID: return "resid";
		case pt_EOLN: return "end of line";

		case pt_KW_ALL: return "ALL";
		case pt_KW_CHAIN: return "CHAIN";
		case pt_KW_RESSEQ: return "RESSEQ";
		case pt_KW_RESID: return "RESID";
		case pt_KW_RESNAME: return "RESNAME";
		case pt_KW_ELEMENT: return "ELEMENT";
		case pt_KW_AND: return "AND";
		case pt_KW_OR: return "OR";
		case pt_KW_NOT: return "NOT";
		case pt_KW_PDB: return "PDB";
		case pt_KW_ENTRY: return "ENTRY";
		case pt_KW_THROUGH: return "THROUGH";

		default: return "character";
	}
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::Parse()
{
	if (m_lookahead == pt_KW_PDB)
	{
		match(pt_KW_PDB);
		//		Match(pt_KW_ENTRY);

		throw std::runtime_error("Unimplemented PDB ENTRY specification");
	}

	std::unique_ptr<tls_selection> result = ParseAtomSelection();

	bool extraParenthesis = false;

	if (m_lookahead == ')')
	{
		extraParenthesis = true;
		m_lookahead = get_next_token();
	}

	match(pt_EOLN);

	if (extraParenthesis)
		std::cerr << "WARNING: too many closing parenthesis in TLS selection statement\n";

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::ParseAtomSelection()
{
	std::unique_ptr<tls_selection> result = ParseTerm();

	while (m_lookahead == pt_KW_OR)
	{
		match(pt_KW_OR);
		result = std::make_unique<tls_selection_union>(result, ParseTerm());
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::ParseTerm()
{
	std::unique_ptr<tls_selection> result = ParseFactor();

	while (m_lookahead == pt_KW_AND)
	{
		match(pt_KW_AND);
		result = std::make_unique<tls_selection_intersection>(result, ParseFactor());
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplPhenix::ParseFactor()
{
	std::unique_ptr<tls_selection> result;

	switch (m_lookahead)
	{
		case '(':
			match('(');
			result = ParseAtomSelection();
			if (m_lookahead == pt_EOLN)
				std::cerr << "WARNING: missing closing parenthesis in TLS selection statement\n";
			else
				match(')');
			break;

		case pt_KW_NOT:
			match(pt_KW_NOT);
			result = std::make_unique<tls_selection_not>(ParseAtomSelection());
			break;

		case pt_KW_CHAIN:
		{
			match(pt_KW_CHAIN);

			std::string chainID = m_value_s;
			if (m_lookahead == pt_NUMBER) // sigh
			{
				chainID = to_string(m_value_i);
				match(pt_NUMBER);
			}
			else
				match(m_lookahead == pt_STRING ? pt_STRING : pt_IDENT);

			result = std::make_unique<tls_selection_chain>(chainID);
			break;
		}

		case pt_KW_RESNAME:
		{
			match(pt_KW_RESNAME);
			std::string name = m_value_s;
			match(pt_IDENT);
			result = std::make_unique<tls_selection_by_name>(name);
			break;
		}

		case pt_KW_ELEMENT:
		{
			match(pt_KW_ELEMENT);
			std::string element = m_value_s;
			match(pt_IDENT);
			result = std::make_unique<tls_selection_by_element>(element);
			break;
		}

		case pt_KW_RESSEQ:
		{
			match(pt_KW_RESSEQ);

			int from = m_value_i;
			match(pt_NUMBER);

			int to = from;
			if (m_lookahead == ':')
			{
				match(':');
				to = m_value_i;
				match(pt_NUMBER);
			}

			result = std::make_unique<tls_selection_range_seq>(from, to);
			break;
		}

		case pt_KW_RESID:
		{
			match(pt_KW_RESID);

			int from, to;
			char icode_from = 0, icode_to = 0;
			bool through = false;

			from = to = m_value_i;

			if (m_lookahead == pt_NUMBER)
				match(pt_NUMBER);
			else
			{
				icode_from = m_icode;
				match(pt_RESID);
			}

			if (m_lookahead == ':' or m_lookahead == pt_KW_THROUGH or m_lookahead == '-')
			{
				through = m_lookahead == pt_KW_THROUGH;

				match(m_lookahead);

				to = m_value_i;
				if (m_lookahead == pt_NUMBER)
					match(pt_NUMBER);
				else
				{
					icode_to = m_icode;
					match(pt_RESID);
				}

				if (through)
					result = std::make_unique<tls_selection_range_id>(from, to, icode_from, icode_to);
				else
				{
					if (cif::VERBOSE and (icode_from or icode_to))
						std::cerr << "Warning, ignoring insertion codes\n";

					result = std::make_unique<tls_selection_range_seq>(from, to);
				}
			}
			else
				result = std::make_unique<tls_selection_res_id>(from, icode_from);

			break;
		}

		case pt_KW_ALL:
			match(pt_KW_ALL);
			result = std::make_unique<tls_selection_all>();
			break;

		default:
			throw std::runtime_error("Unexpected token " + to_string(m_lookahead) + " (" + m_token + ')');
	}

	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserImplBuster : public tls_selection_parser_impl
{
  public:
	TLSSelectionParserImplBuster(const std::string &selection);

	std::unique_ptr<tls_selection> Parse() override;

  protected:
	enum TOKEN
	{
		bt_NONE = 0,
		bt_IDENT = 256,
		bt_NUMBER = 257,
		bt_EOLN = 258,
	};

	int get_next_token() override;
	std::string to_string(int token) override;

	std::unique_ptr<tls_selection> ParseGroup();
	std::tuple<std::string, int> ParseAtom();

	std::unique_ptr<tls_selection> ParseOldGroup();

	int m_value_i;
	std::string m_value_s;
	bool m_parsing_old_style = false;
};

TLSSelectionParserImplBuster::TLSSelectionParserImplBuster(const std::string &selection)
	: tls_selection_parser_impl(selection)
{
	m_lookahead = get_next_token();
}

int TLSSelectionParserImplBuster::get_next_token()
{
	int result = bt_NONE;
	enum STATE
	{
		st_START,
		st_NEGATE,
		st_NUM,
		st_IDENT
	} state = st_START;

	m_value_i = 0;
	m_value_s.clear();
	bool negative = false;

	while (result == bt_NONE)
	{
		char ch = *m_p++;
		if (m_p > m_end)
			ch = 0;

		switch (state)
		{
			case st_START:
				if (ch == 0)
					result = bt_EOLN;
				else if (isspace(ch))
					continue;
				else if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
				}
				else if (isalpha(ch))
				{
					m_value_s = { ch };
					state = st_IDENT;
				}
				else if (ch == '-')
				{
					state = st_NEGATE;
				}
				else
					result = static_cast<unsigned char>(ch);
				break;

			case st_NEGATE:
				if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
					negative = true;
				}
				else
				{
					--m_p;
					result = '-';
				}
				break;

			case st_NUM:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (ch - '0');
				else
				{
					if (negative)
						m_value_i = -m_value_i;

					result = bt_NUMBER;
					--m_p;
				}
				break;

			case st_IDENT:
				if (isalnum(ch))
					m_value_s += ch;
				else
				{
					--m_p;
					result = bt_IDENT;
				}
				break;
		}
	}

	return result;
}

std::string TLSSelectionParserImplBuster::to_string(int token)
{
	switch (token)
	{
		case bt_IDENT: return "identifier (" + m_value_s + ')';
		case bt_NUMBER: return "number (" + to_string(m_value_i) + ')';
		case bt_EOLN: return "end of line";

		default:
			assert(false);
			return "unknown token";
	}
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBuster::ParseGroup()
{
	std::unique_ptr<tls_selection> result;

	auto add = [&result](const std::string &chainID, int from, int to)
	{
		std::unique_ptr<tls_selection> sc(new tls_selection_chain(chainID));
		std::unique_ptr<tls_selection> sr(new tls_selection_range_seq(from, to));
		std::unique_ptr<tls_selection> s(new tls_selection_intersection(sc, sr));

		if (result == nullptr)
			result.reset(s.release());
		else
			result = std::make_unique<tls_selection_union>( result, s );
	};

	match('{');

	do
	{
		std::string chain1;
		int seqNr1;
		std::tie(chain1, seqNr1) = ParseAtom();

		if (m_lookahead == '-')
		{
			std::string chain2;
			int seqNr2 = seqNr1;

			match('-');

			if (m_lookahead == bt_NUMBER)
			{
				seqNr2 = m_value_i;
				match(bt_NUMBER);
			}
			else
			{
				std::tie(chain2, seqNr2) = ParseAtom();
				if (chain1 != chain2)
				{
					std::cerr << "Warning, ranges over multiple chains detected\n";

					std::unique_ptr<tls_selection> sc1(new tls_selection_chain(chain1));
					std::unique_ptr<tls_selection> sr1(new tls_selection_range_seq(seqNr1, kResidueNrWildcard));
					std::unique_ptr<tls_selection> s1(new tls_selection_intersection(sc1, sr1));

					std::unique_ptr<tls_selection> sc2(new tls_selection_chain(chain2));
					std::unique_ptr<tls_selection> sr2(new tls_selection_range_seq(kResidueNrWildcard, seqNr2));
					std::unique_ptr<tls_selection> s2(new tls_selection_intersection(sc2, sr2));

					std::unique_ptr<tls_selection> s(new tls_selection_union(s1, s2));

					if (result == nullptr)
						result.reset(s.release());
					else
						result = std::make_unique<tls_selection_union>( result, s );

					chain1.clear();
				}
			}

			if (not chain1.empty())
				add(chain1, seqNr1, seqNr2);
		}
		else
			add(chain1, seqNr1, seqNr1);
	} while (m_lookahead != '}');

	match('}');

	return result;
}

std::tuple<std::string, int> TLSSelectionParserImplBuster::ParseAtom()
{
	std::string chain = m_value_s;
	int seqNr = kResidueNrWildcard;

	if (m_lookahead == '*')
		match('*');
	else
		match(bt_IDENT);

	match('|');

	if (m_lookahead == '*')
		match('*');
	else
	{
		seqNr = m_value_i;
		match(bt_NUMBER);

		if (m_lookahead == ':')
		{
			match(':');
			std::string atom = m_value_s;

			if (cif::VERBOSE > 0)
				std::cerr << "Warning: ignoring atom ID '" << atom << "' in TLS selection\n";

			match(bt_IDENT);
		}
	}

	return std::make_tuple(chain, seqNr);
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBuster::Parse()
{
	std::unique_ptr<tls_selection> result = ParseGroup();
	match(bt_EOLN);
	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserImplBusterOld : public tls_selection_parser_impl
{
  public:
	TLSSelectionParserImplBusterOld(const std::string &selection)
		: tls_selection_parser_impl(selection)
	{
		m_lookahead = get_next_token();
	}

	std::unique_ptr<tls_selection> Parse() override;

  private:
	std::unique_ptr<tls_selection> ParseAtomSelection();
	std::unique_ptr<tls_selection> ParseTerm();
	std::unique_ptr<tls_selection> ParseFactor();

	std::unique_ptr<tls_selection> ParseResid();
	std::unique_ptr<tls_selection> ParseChainResid();

	enum TOKEN
	{
		pt_NONE = 0,
		pt_IDENT = 256,
		pt_CHAINRESID = 257,
		pt_STRING = 258,
		pt_NUMBER = 259,
		pt_RANGE = 260,
		pt_EOLN = 261,

		pt_KW_ALL = 262,
		pt_KW_CHAIN = 263,
		pt_KW_RESSEQ = 264,
		pt_KW_RESID = 265,
		pt_KW_RESNAME = 266,
		pt_KW_ELEMENT = 267,
		pt_KW_AND = 268,
		pt_KW_OR = 269,
		pt_KW_NOT = 270,
		pt_KW_PDB = 271,
		pt_KW_ENTRY = 272,
		pt_KW_THROUGH = 273
	};

	int get_next_token() override;
	std::string to_string(int token) override;

	int m_value_i;
	std::string m_value_s;
	int m_value_r[2];
};

int TLSSelectionParserImplBusterOld::get_next_token()
{
	int result = pt_NONE;
	enum STATE
	{
		st_START,
		st_NEGATE,
		st_NUM,
		st_RANGE,
		st_IDENT_1,
		st_IDENT,
		st_CHAINRESID,
		st_QUOTED_1,
		st_QUOTED_2
	} state = st_START;

	m_value_i = 0;
	m_value_s.clear();

	bool negative = false;

	while (result == pt_NONE)
	{
		char ch = *m_p++;
		if (m_p > m_end)
			ch = 0;

		switch (state)
		{
			case st_START:
				if (ch == 0)
					result = pt_EOLN;
				else if (isspace(ch))
					continue;
				else if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
				}
				else if (isalpha(ch))
				{
					m_value_s = { ch };
					state = st_IDENT_1;
				}
				else if (ch == '-')
				{
					state = st_NEGATE;
				}
				else if (ch == '\'')
				{
					state = st_QUOTED_1;
				}
				else
					result = static_cast<unsigned char>(ch);
				break;

			case st_NEGATE:
				if (isdigit(ch))
				{
					m_value_i = ch - '0';
					state = st_NUM;
					negative = true;
				}
				else
				{
					--m_p;
					result = '-';
				}
				break;

			case st_NUM:
				if (isdigit(ch))
					m_value_i = 10 * m_value_i + (ch - '0');
				else if (ch == '-' or ch == ':')
				{
					if (negative)
						m_value_i = -m_value_i;

					m_value_r[0] = m_value_i;
					m_value_r[1] = 0;
					state = st_RANGE;
				}
				else
				{
					if (negative)
						m_value_i = -m_value_i;

					result = pt_NUMBER;
					--m_p;
				}
				break;

			case st_RANGE: // TODO: question, is "-2--1" a valid range? We do not support that, yet
				if (isdigit(ch))
					m_value_r[1] = 10 * m_value_r[1] + (ch - '0');
				else if (m_value_r[1] != 0)
				{
					result = pt_RANGE;
					--m_p;
				}
				else
				{
					--m_p;
					--m_p;
					result = pt_NUMBER;
				}
				break;

			case st_IDENT_1:
				if (isalpha(ch))
				{
					m_value_s += ch;
					state = st_IDENT;
				}
				else if (isdigit(ch))
				{
					m_value_i = (ch - '0');
					state = st_CHAINRESID;
				}
				else
				{
					--m_p;
					result = pt_IDENT;
				}
				break;

			case st_CHAINRESID:
				if (isalpha(ch))
				{
					m_value_s += to_string(m_value_i);
					m_value_s += ch;
					state = st_IDENT;
				}
				else if (isdigit(ch))
					m_value_i = 10 * m_value_i + (ch - '0');
				else
				{
					--m_p;
					result = pt_CHAINRESID;
				}
				break;

			case st_IDENT:
				if (isalnum(ch))
					m_value_s += ch;
				else
				{
					--m_p;
					result = pt_IDENT;
				}
				break;

			case st_QUOTED_1:
				if (ch == '\'')
				{
					--m_p;
					result = '\'';
				}
				else
				{
					m_value_s = { ch };
					state = st_QUOTED_2;
				}
				break;

			case st_QUOTED_2:
				if (ch == '\'')
					result = pt_STRING;
				else if (ch == 0)
					throw std::runtime_error("Unexpected end of selection, missing quote character?");
				else
					m_value_s += ch;
				break;
		}
	}

	if (result == pt_IDENT)
	{
		if (cif::iequals(m_value_s, "CHAIN"))
			result = pt_KW_CHAIN;
		else if (cif::iequals(m_value_s, "ALL"))
			result = pt_KW_ALL;
		else if (cif::iequals(m_value_s, "AND"))
			result = pt_KW_AND;
		else if (cif::iequals(m_value_s, "OR"))
			result = pt_KW_OR;
		else if (cif::iequals(m_value_s, "NOT"))
			result = pt_KW_NOT;
		else if (cif::iequals(m_value_s, "RESSEQ"))
			result = pt_KW_RESSEQ;
		else if (cif::iequals(m_value_s, "RESID") or cif::iequals(m_value_s, "RESI") or cif::iequals(m_value_s, "RESIDUES"))
			result = pt_KW_RESID;
		else if (cif::iequals(m_value_s, "RESNAME"))
			result = pt_KW_RESNAME;
		else if (cif::iequals(m_value_s, "PDB"))
			result = pt_KW_PDB;
		else if (cif::iequals(m_value_s, "ENTRY"))
			result = pt_KW_ENTRY;
		else if (cif::iequals(m_value_s, "THROUGH"))
			result = pt_KW_THROUGH;
	}

	return result;
}

std::string TLSSelectionParserImplBusterOld::to_string(int token)
{
	switch (token)
	{
		case pt_IDENT: return "identifier (" + m_value_s + ')';
		case pt_STRING: return "std::string (" + m_value_s + ')';
		case pt_NUMBER: return "number (" + to_string(m_value_i) + ')';
		case pt_RANGE: return "range (" + to_string(m_value_r[0]) + ':' + to_string(m_value_r[1]) + ')';
		case pt_EOLN: return "end of line";

		case pt_KW_ALL: return "ALL";
		case pt_KW_CHAIN: return "CHAIN";
		case pt_KW_RESSEQ: return "RESSEQ";
		case pt_KW_RESID: return "RESID";
		case pt_KW_RESNAME: return "RESNAME";
		case pt_KW_ELEMENT: return "ELEMENT";
		case pt_KW_AND: return "AND";
		case pt_KW_OR: return "OR";
		case pt_KW_NOT: return "NOT";
		case pt_KW_PDB: return "PDB";
		case pt_KW_ENTRY: return "ENTRY";
		case pt_KW_THROUGH: return "THROUGH";
		default:
			assert(false);
			return "unknown token";
	}
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::Parse()
{
	if (m_lookahead == pt_KW_PDB)
	{
		match(pt_KW_PDB);
		//		Match(pt_KW_ENTRY);

		throw std::runtime_error("Unimplemented PDB ENTRY specification");
	}

	std::unique_ptr<tls_selection> result = ParseAtomSelection();

	match(pt_EOLN);

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseAtomSelection()
{
	std::unique_ptr<tls_selection> result = ParseTerm();

	while (m_lookahead == pt_KW_OR)
	{
		match(pt_KW_OR);
		result = std::make_unique<tls_selection_union>(result, ParseTerm());
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseTerm()
{
	std::unique_ptr<tls_selection> result = ParseFactor();

	while (m_lookahead == pt_KW_AND)
	{
		match(pt_KW_AND);
		result = std::make_unique<tls_selection_intersection>(result, ParseFactor());
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseFactor()
{
	std::unique_ptr<tls_selection> result;

	switch (m_lookahead)
	{
		case '(':
			match('(');
			result = ParseAtomSelection();
			match(')');
			break;

		case pt_KW_NOT:
			match(pt_KW_NOT);
			result = std::make_unique<tls_selection_not>(ParseAtomSelection());
			break;

		case pt_KW_CHAIN:
		{
			match(pt_KW_CHAIN);

			std::string chainID = m_value_s;
			if (m_lookahead == pt_NUMBER) // sigh
			{
				chainID = to_string(m_value_i);
				match(pt_NUMBER);
			}
			else
				match(m_lookahead == pt_STRING ? pt_STRING : pt_IDENT);

			result = std::make_unique<tls_selection_chain>(chainID);
			break;
		}

		case pt_KW_RESNAME:
		{
			match(pt_KW_RESNAME);
			std::string name = m_value_s;
			match(pt_IDENT);
			result = std::make_unique<tls_selection_by_name>(name);
			break;
		}

		case pt_KW_RESSEQ:
			match(pt_KW_RESSEQ);
			result = ParseResid();
			break;

		case pt_KW_RESID:
			match(pt_KW_RESID);
			result = ParseResid();
			break;

		case pt_KW_ALL:
			match(pt_KW_ALL);
			result = std::make_unique<tls_selection_all>();
			break;

		case pt_CHAINRESID:
			result = ParseChainResid();
			break;

		default:
			throw std::runtime_error("Unexpected token " + to_string(m_lookahead));
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseResid()
{
	std::unique_ptr<tls_selection> result;

	for (;;)
	{
		int from, to;

		if (m_lookahead == pt_RANGE)
		{
			from = m_value_r[0];
			to = m_value_r[1];
			match(pt_RANGE);
		}
		else
		{
			from = m_value_i;
			match(pt_NUMBER);

			to = from;
			if (m_lookahead == ':' or m_lookahead == '-' or m_lookahead == pt_KW_THROUGH)
			{
				match(m_lookahead);
				to = m_value_i;
				match(pt_NUMBER);
			}
		}

		std::unique_ptr<tls_selection> range(new tls_selection_range_seq(from, to));

		if (result)
			result = std::make_unique<tls_selection_union>(result, range);
		else
			result.reset(range.release());

		if (m_lookahead == ',')
		{
			match(',');
			continue;
		}

		break;
	}

	return result;
}

std::unique_ptr<tls_selection> TLSSelectionParserImplBusterOld::ParseChainResid()
{
	std::unique_ptr<tls_selection> result;

	for (;;)
	{
		int from, to;

		from = to = m_value_i;
		std::string chainID = m_value_s;

		match(pt_CHAINRESID);

		if (m_lookahead == '-')
		{
			match(m_lookahead);
			to = m_value_i;

			if (m_value_s != chainID)
				throw std::runtime_error("Cannot have two different chainIDs in a range selection");

			match(pt_CHAINRESID);
		}

		std::unique_ptr<tls_selection> sc(new tls_selection_chain(chainID));
		std::unique_ptr<tls_selection> sr(new tls_selection_range_seq(from, to));
		std::unique_ptr<tls_selection> range(new tls_selection_intersection(sc, sr));

		if (result)
			result = std::make_unique<tls_selection_union>(result, range);
		else
			result.reset(range.release());

		if (m_lookahead == ',')
		{
			match(',');
			continue;
		}

		break;
	}

	return result;
}

// --------------------------------------------------------------------

class TLSSelectionParserBase
{
  public:
	[[nodiscard]] virtual std::unique_ptr<tls_selection> Parse(const std::string &selection) const = 0;
	virtual ~TLSSelectionParserBase() = default;
};

template <typename IMPL>
class TLSSelectionParser
{
  public:
	[[nodiscard]] virtual std::unique_ptr<tls_selection> Parse(const std::string &selection) const
	{
		std::unique_ptr<tls_selection> result;

		try
		{
			IMPL p(selection);
			result = p.Parse();
		}
		catch (const std::exception &ex)
		{
			std::cerr << "ParseError: " << ex.what() << '\n';
		}

		return result;
	}
};

// --------------------------------------------------------------------

std::unique_ptr<tls_selection> parse_tls_selection_details(const std::string &program, const std::string &selection)
{
	TLSSelectionParser<TLSSelectionParserImplPhenix> phenix;
	TLSSelectionParser<TLSSelectionParserImplBuster> buster;
	TLSSelectionParser<TLSSelectionParserImplBusterOld> busterOld;

	std::unique_ptr<tls_selection> result;

	if (cif::icontains(program, "buster"))
	{
		result = buster.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to old BUSTER\n";
			result = busterOld.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to PHENIX\n";
			result = phenix.Parse(selection);
		}
	}
	else if (cif::icontains(program, "phenix"))
	{
		result = phenix.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to BUSTER\n";
			result = buster.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to old BUSTER\n";
			result = busterOld.Parse(selection);
		}
	}
	else
	{
		if (cif::VERBOSE > 0)
			std::cerr << "No known program specified, trying PHENIX\n";

		result = phenix.Parse(selection);

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to BUSTER\n";
			result = buster.Parse(selection);
		}

		if (not result)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Falling back to old BUSTER\n";
			result = busterOld.Parse(selection);
		}
	}

	return result;
}

} // namespace cif
