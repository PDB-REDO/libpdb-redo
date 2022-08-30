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
   Date: dinsdag 19 juni, 2018
*/

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>

#include <map>
#include <mutex>

#include <clipper/clipper.h>

#include "pdb-redo/Ramachandran.hpp"

namespace pdb_redo
{

const double kPI = M_PI;

// --------------------------------------------------------------------

class RamachandranTables
{
  public:
	static RamachandranTables &instance()
	{
		std::lock_guard lock(sMutex);

		static RamachandranTables sInstance;
		return sInstance;
	}

	clipper::Ramachandran &table(const std::string &aa, bool prePro)
	{
		std::lock_guard lock(sMutex);

		auto i = mTables.find(std::make_tuple(aa, prePro));
		if (i == mTables.end())
		{
			clipper::Ramachandran::TYPE type;

			if (aa == "GLY")
				type = clipper::Ramachandran::Gly2;
			else if (aa == "PRO")
				type = clipper::Ramachandran::Pro2;
			else if (aa == "ILE" or aa == "VAL")
				type = clipper::Ramachandran::IleVal2;
			else if (prePro)
				type = clipper::Ramachandran::PrePro2;
			else
				type = clipper::Ramachandran::NoGPIVpreP2;

			i = mTables.emplace(make_pair(std::make_tuple(aa, prePro), clipper::Ramachandran(type))).first;
		}

		return i->second;
	}

  private:
	std::map<std::tuple<std::string, int>, clipper::Ramachandran> mTables;
	static std::mutex sMutex;
};

std::mutex RamachandranTables::sMutex;

float calculateRamachandranZScore(const std::string &aa, bool prePro, float phi, float psi)
{
	auto &table = RamachandranTables::instance().table(aa, prePro);
	return static_cast<float>(table.probability(phi * kPI / 180, psi * kPI / 180));
}

RamachandranScore calculateRamachandranScore(const std::string &aa, bool prePro, float phi, float psi)
{
	auto &table = RamachandranTables::instance().table(aa, prePro);

	phi *= static_cast<float>(kPI / 180);
	psi *= static_cast<float>(kPI / 180);

	RamachandranScore result;

	if (table.favored(phi, psi))
		result = rsFavoured;
	else if (table.allowed(phi, psi))
		result = rsAllowed;
	else
		result = rsNotAllowed;

	return result;
}

} // namespace pdb_redo