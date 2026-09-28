// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

/*
   Created by: Maarten L. Hekkelman
   Date: dinsdag 19 juni, 2018
*/

#include <cassert>
#include <cmath>

#include <map>
#include <mutex>

#include <clipper/clipper.h>

#include "pdb-redo/Ramachandran.hpp"

namespace pdb_redo
{

const float kPI = std::numbers::pi_v<float>;

// --------------------------------------------------------------------

class RamachandranTables
{
  public:
	static RamachandranTables &instance()
	{
		std::scoped_lock lock(sMutex);

		static RamachandranTables sInstance;
		return sInstance;
	}

	clipper::Ramachandran &table(const std::string &aa, bool prePro)
	{
		std::scoped_lock lock(sMutex);

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

			i = mTables.emplace(std::make_tuple(aa, prePro), clipper::Ramachandran(type)).first;
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

	phi *= kPI / 180;
	psi *= kPI / 180;

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