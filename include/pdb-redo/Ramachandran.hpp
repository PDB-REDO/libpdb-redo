// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

/*
   Created by: Maarten L. Hekkelman
   Date: dinsdag 19 juni, 2018
*/

#pragma once

namespace pdb_redo
{

float calculateRamachandranZScore(const std::string &aa, bool prePro, float phi, float psi);

enum RamachandranScore
{
	rsNotAllowed,
	rsAllowed,
	rsFavoured
};

RamachandranScore calculateRamachandranScore(const std::string &aa, bool prePro, float phi, float psi);

} // namespace pdb_redo
