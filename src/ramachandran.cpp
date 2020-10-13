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

#include "cif++/Config.hpp"

#include <cmath>
#include <cassert>

#include <map>
#include <mutex>

#include <clipper/clipper.h>

#include "cif++/Point.hpp"

#include "ramachandran.h"

// --------------------------------------------------------------------

using namespace std;

// --------------------------------------------------------------------

class RamachandranTables
{
  public:
	static RamachandranTables& instance()
	{
		lock_guard<mutex> lock(sMutex);
		
		static RamachandranTables sInstance;
		return sInstance;
	}
	
	clipper::Ramachandran& table(const string& aa, bool prePro)
	{
		lock_guard<mutex> lock(sMutex);
		
		auto i = mTables.find(make_tuple(aa, prePro));
		if (i == mTables.end())
		{
			clipper::Ramachandran::TYPE type;
			
				 if (aa == "GLY")	type = clipper::Ramachandran::Gly2;
			else if (aa == "PRO")	type = clipper::Ramachandran::Pro2;
			else if (aa == "ILE" or aa == "VAL")
									type = clipper::Ramachandran::IleVal2;
			else if (prePro)
									type = clipper::Ramachandran::PrePro2;
			else					type = clipper::Ramachandran::NoGPIVpreP2;
			
			i = mTables.emplace(make_pair(make_tuple(aa, prePro), clipper::Ramachandran(type))).first;
		}
		
		return i->second;
	}
	
  private:
	map<tuple<string,int>,clipper::Ramachandran> mTables;
	static mutex sMutex;
};

mutex RamachandranTables::sMutex;

float calculateRamachandranZScore(const string& aa, bool prePro, float phi, float psi)
{
	auto& table = RamachandranTables::instance().table(aa, prePro);
	return table.probability(phi * mmcif::kPI / 180, psi * mmcif::kPI / 180);
}

RamachandranScore calculateRamachandranScore(const string& aa, bool prePro, float phi, float psi)
{
	auto& table = RamachandranTables::instance().table(aa, prePro);

	phi *= mmcif::kPI / 180;
	psi *= mmcif::kPI / 180;
	
	RamachandranScore result;

	if (table.favored(phi, psi))
		result = rsFavoured;
	else if (table.allowed(phi, psi))
		result = rsAllowed;
	else
		result = rsNotAllowed;
	
	return result;
}

//template<typename T, int N>
//struct RamachandranData
//{
//	const char									aa[4];
//	const RamachandranSecondaryStructureType	ss;
//	const float									average;
//	const float									sd;
//	const T										counts[N * N];
//	
//	T count(size_t bin1, size_t bin2) const
//	{
//		size_t ix = bin1 * N + bin2;
//		assert(ix < N * N);
//		return counts[ix];
//	}
//	
//	float operator()(float phi, float psi) const
//	{
//		size_t phiFloorIx = static_cast<size_t>(N * (phi + 180) / 360);
//		size_t psiFloorIx = static_cast<size_t>(N * (psi + 180) / 360);
//		
//		size_t phiCeilIx = phiFloorIx + 1 < N ? phiFloorIx + 1 : phiFloorIx;
//		size_t psiCeilIx = psiFloorIx + 1 < N ? psiFloorIx + 1 : psiFloorIx;
//
//		float phiFloorAngle = (phiFloorIx * 360.0f) / N - 180;
//		float psiFloorAngle = (psiFloorIx * 360.0f) / N - 180;
//
//		float phiCeilAngle = (phiCeilIx * 360.0f) / N - 180;
//		float psiCeilAngle = (psiCeilIx * 360.0f) / N - 180;
//
//		float phiFactor = phiCeilIx > phiFloorIx ? (phi - phiFloorAngle) / (phiCeilAngle - phiFloorAngle) : 1;
//		float psiFactor = psiCeilIx > psiFloorIx ? (psi - psiFloorAngle) / (psiCeilAngle - psiFloorAngle) : 1;
//
//		float c1 = count(phiFloorIx, psiFloorIx) + (count(phiCeilIx, psiFloorIx) - count(phiFloorIx, psiFloorIx)) * phiFactor;
//		float c2 = count(phiFloorIx, psiCeilIx) + (count(phiCeilIx, psiCeilIx) - count(phiFloorIx, psiCeilIx)) * phiFactor;
//		
//		float iCount = c1 + (c2 - c1) * psiFactor;
//		
//		return (iCount - average) / sd;
//	}
//};
//
//RamachandranData<uint16_t,120> kRamachandranData[] = {
//
//// The next file should be generated by the script generate-rama-datastructures.pl
//#include "ramachandran-data.inl"
//	
//};
//
//// --------------------------------------------------------------------
//
//float calculateRamachandranZScore(const string& aa, RamachandranSecondaryStructureType ss,
//	float phi, float psi)
//{
//	float result = -999;
//	
//	for (auto& d: kRamachandranData)
//	{
//		if (ss != d.ss or aa != d.aa)
//			continue;
//		
//		result = d(phi, psi);
//		break;
//	}
//	
//	return result;
//}
