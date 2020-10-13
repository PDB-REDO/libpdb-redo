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

#pragma once

#include <set>

#include <zeep/xml/serialize.hpp>

#include "cif++/Structure.hpp"

class NotAHBondSet
{
  public:
	~NotAHBondSet() {}

	NotAHBondSet(const NotAHBondSet&) = delete;
	NotAHBondSet& operator=(const NotAHBondSet&) = delete;

	static NotAHBondSet* Create();

	void Save(zeep::xml::element& inNode);
	static NotAHBondSet* Load(zeep::xml::element& inNode);
	
	bool IsHBondDonorOrAcceptor(const std::string& inMonomer,
		const std::string& inAtomID) const
	{
		return mData.count({inMonomer, inAtomID}) > 0;
	}
  
	bool operator()(const mmcif::Atom& atom) const
	{
		return atom.type() == mmcif::N and
			IsHBondDonorOrAcceptor(atom.comp().id(), atom.labelAtomID());
	}
  
  private:

	struct NotAHBondAtom
	{
		std::string	monomer;
		std::string	atomId;
		
		bool operator<(const NotAHBondAtom& rhs) const
		{
			int d = monomer.compare(rhs.monomer);
			if (d == 0)
				d = atomId.compare(rhs.atomId);
			return d < 0;
		}
	};
	
	typedef std::set<NotAHBondAtom> NotAHBondAtomSet;

	NotAHBondSet(NotAHBondAtomSet&& data)
		: mData(std::move(data))
	{
	}
	
	NotAHBondAtomSet mData;
};
