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

#include <cif++/datablock.hpp>
#include <vector>

#include "cif++.hpp"

namespace pdb_redo
{

// --------------------------------------------------------------------
// The chemical composition of the structure in an mmCIF file is
// defined in the class composition. A compositon consists of
// entities. Each Entity can be either a polymer, a non-polymer
// a macrolide or a water molecule.
// Entities themselves are made up of compounds. And compounds
// contain CompoundAtom records for each atom.

class Compound;
class Link;
struct CompoundAtom;

enum BondType
{
	singleBond,
	doubleBond,
	tripleBond,
	aromaticBond,
	delocalizedBond
};

// --------------------------------------------------------------------
// struct containing information about an atom in a chemical compound
// This information comes from the CCP4 monomer library.

struct CompoundAtom
{
	std::string id;
	cif::atom_type typeSymbol;
	std::string typeEnergy;
	float partialCharge;
	float x, y, z;
};

// --------------------------------------------------------------------
// struct containing information about the bonds
// This information comes from the CCP4 monomer library.

struct CompoundBond
{
	std::string atomID[2];
	BondType type;
	bool aromatic;
	float distance;
	float esd;
};

// --------------------------------------------------------------------
// struct containing information about the bond-angles
// This information comes from the CCP4 monomer library.

struct CompoundAngle
{
	std::string atomID[3];
	float angle;
	float esd;
};

// --------------------------------------------------------------------
// struct containing information about the bond-angles
// This information comes from the CCP4 monomer library.

struct CompoundTorsion
{
	std::string atomID[4];
	float angle;
	float esd;
	int period;
};

// --------------------------------------------------------------------
// struct containing information about the bond-angles
// This information comes from the CCP4 monomer library.

struct CompoundPlane
{
	std::string id;
	std::vector<std::string> atomID;
	float esd;
};

// --------------------------------------------------------------------
// struct containing information about a chiral centre
// This information comes from the CCP4 monomer library.

enum ChiralVolumeSign
{
	negativ,
	positiv,
	both
};

struct CompoundChiralCentre
{
	std::string id;
	std::string atomIDCentre;
	std::string atomID[3];
	ChiralVolumeSign volumeSign;
};

// --------------------------------------------------------------------
// a class that contains information about a chemical compound.
// This information is derived from the ccp4 monomer library by default.
// To create compounds, you'd best use the factory method.

class Compound
{
  public:
	Compound(const cif::datablock &db, const std::string &id, std::string name,
		std::string group);

	// accessors
	[[nodiscard]] std::string id() const { return mID; }
	[[nodiscard]] std::string name() const { return mName; }
	[[nodiscard]] std::string type() const;
	[[nodiscard]] std::string group() const { return mGroup; }
	[[nodiscard]] std::vector<CompoundAtom> atoms() const { return mAtoms; }
	[[nodiscard]] std::vector<CompoundBond> bonds() const { return mBonds; }
	[[nodiscard]] std::vector<CompoundAngle> angles() const { return mAngles; }
	[[nodiscard]] std::vector<CompoundChiralCentre> chiralCentres() const
	{
		return mChiralCentres;
	}
	[[nodiscard]] std::vector<CompoundPlane> planes() const { return mPlanes; }
	[[nodiscard]] std::vector<CompoundTorsion> torsions() const { return mTorsions; }

	[[nodiscard]] CompoundAtom get_atom_by_atom_id(const std::string &atomID) const;

	[[nodiscard]] bool atomsBonded(const std::string &atomId_1, const std::string &atomId_2) const;
	[[nodiscard]] float atomBondValue(const std::string &atomId_1, const std::string &atomId_2) const;
	[[nodiscard]] float bondAngle(const std::string &atomId_1, const std::string &atomId_2, const std::string &atomId_3) const;
	[[nodiscard]] float chiralVolume(const std::string &centreID) const;

	[[nodiscard]] std::string formula() const;
	[[nodiscard]] float formulaWeight() const;
	[[nodiscard]] int charge() const;
	[[nodiscard]] bool isWater() const;
	[[nodiscard]] bool isSugar() const;

	// std::vector<std::string> isomers() const;
	// bool isIsomerOf(const Compound &c) const;
	// std::vector<std::tuple<std::string, std::string>> mapToIsomer(const Compound &c) const;

	/// @brief Return the content of this restraint compound in a CCD format
	/// @return Datablock containing the CCD information for this compound
	[[nodiscard]] cif::datablock generateCCDCompound() const;

	[[nodiscard]] std::string getDescriptor(std::string_view type) const;

  private:
	cif::datablock mCF;
	std::string mID;
	std::string mName;
	std::string mGroup;
	std::vector<CompoundAtom> mAtoms;
	std::vector<CompoundBond> mBonds;
	std::vector<CompoundAngle> mAngles;
	std::vector<CompoundTorsion> mTorsions;
	std::vector<CompoundChiralCentre> mChiralCentres;
	std::vector<CompoundPlane> mPlanes;
};

// --------------------------------------------------------------------
// struct containing information about the bonds
// This information comes from the CCP4 monomer library.

struct LinkAtom
{
	int compID;
	std::string atomID;

	bool operator==(const LinkAtom &rhs) const { return compID == rhs.compID and atomID == rhs.atomID; }
};

struct LinkBond
{
	LinkAtom atom[2];
	BondType type;
	float distance;
	float esd;
};

// --------------------------------------------------------------------
// struct containing information about the bond-angles
// This information comes from the CCP4 monomer library.

struct LinkAngle
{
	LinkAtom atom[3];
	float angle;
	float esd;
};

// --------------------------------------------------------------------
// struct containing information about the bond-torsions
// This information comes from the CCP4 monomer library.

struct LinkTorsion
{
	LinkAtom atom[4];
	float angle;
	float esd;
	int period;
};

// --------------------------------------------------------------------
// struct containing information about the bond-angles
// This information comes from the CCP4 monomer library.

struct LinkPlane
{
	std::string id;
	std::vector<LinkAtom> atoms;
	float esd;
};

// --------------------------------------------------------------------
// struct containing information about a chiral centre
// This information comes from the CCP4 monomer library.

struct LinkChiralCentre
{
	std::string id;
	LinkAtom atomCentre;
	LinkAtom atom[3];
	ChiralVolumeSign volumeSign;
};

// --------------------------------------------------------------------
// a class that contains information about a chemical link between compounds.
// This information is derived from the ccp4 monomer library by default.

class Link
{
  public:
	Link(cif::datablock &db);

	// accessors
	[[nodiscard]] std::string id() const { return mID; }
	[[nodiscard]] std::vector<LinkBond> bonds() const { return mBonds; }
	[[nodiscard]] std::vector<LinkAngle> angles() const { return mAngles; }
	[[nodiscard]] std::vector<LinkChiralCentre> chiralCentres() const { return mChiralCentres; }
	[[nodiscard]] std::vector<LinkPlane> planes() const { return mPlanes; }
	[[nodiscard]] std::vector<LinkTorsion> torsions() const { return mTorsions; }

	[[nodiscard]] float atomBondValue(const LinkAtom &atomId_1, const LinkAtom &atomId_2) const;
	[[nodiscard]] float bondAngle(const LinkAtom &atomId_1, const LinkAtom &atomId_2, const LinkAtom &atomId_3) const;

	/// \brief Calculate the target chiral volume for \a id for the link between \a compound_id_1 and \a compound_id_2
	/// The compound id's are required to calculate standard bond lengths in case these are not recorded in the link record
	[[nodiscard]] float chiralVolume(const std::string &id, const std::string &compound_id_1, const std::string &compound_id_2) const;

  private:

	std::string mID;
	std::vector<LinkBond> mBonds;
	std::vector<LinkAngle> mAngles;
	std::vector<LinkTorsion> mTorsions;
	std::vector<LinkChiralCentre> mChiralCentres;
	std::vector<LinkPlane> mPlanes;
};

// --------------------------------------------------------------------
// Factory class for Compound and Link objects

class CompoundFactory
{
  public:
	static CompoundFactory &instance();

	const Compound *get(std::string id);
	const Compound *create(std::string id);

	const Link *getLink(std::string id);
	const Link *createLink(std::string id);

	~CompoundFactory();

	void pushDictionary(const std::filesystem::path &inDictFile);
	void pushDictionary(std::istream &inDictionary);
	void pushDictionary(cif::file cf);

	void popDictionary();

  private:
	CompoundFactory();

	class CompoundFactoryImpl *mImpl{};
};

} // namespace pdb_redo