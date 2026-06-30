// Copyright NKI/AVL 2026
//
// SPDX-License-Identifier: BSD-2-Clause

#include "pdb-redo/Conformers.hpp"

#include "pdb-redo/Compound.hpp"

#include <Geometry/point.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <random>

namespace pdb_redo
{

auto createConformers(const std::string &lig, int N)
{
	static std::random_device rd;
	static std::mt19937_64 rng(rd());
	std::uniform_int_distribution<int> rand(1, std::numeric_limits<int>::max());

	int seed = rand(rng);

	// First create an mmCIF file containing one copy of the compound name lig
	cif::file restraintFile;
	restraintFile.emplace(lig);
	restraintFile.front().load_dictionary("mmcif_pdbx.dic");

	auto savedVerbose = std::exchange(cif::VERBOSE, -1);
	cif::mm::structure s(restraintFile);
	cif::VERBOSE = savedVerbose;

	auto compound = pdb_redo::CompoundFactory::instance().create(lig);
	if (compound == nullptr)
		throw std::runtime_error("Failed to create compound for " + lig);

	// construct an entity from scratch based on the atoms
	std::vector<cif::row_initializer> atoms;
	std::vector<std::string> atomIDs;

	for (auto &&[id, typeSymbol, energy, charge, x, y, z] : compound->atoms())
	{
		cif::row_initializer atom({ { "label_atom_id", id },
			{ "auth_atom_id", id },
			{ "type_symbol", cif::atom_type_traits(typeSymbol).symbol() },
			{ "cartn_x", { x, 3 } },
			{ "cartn_y", { y, 3 } },
			{ "cartn_z", { z, 3 } } });

		atoms.emplace_back(std::move(atom));

		if (typeSymbol != cif::H)
			atomIDs.emplace_back(id);
	}

	auto entityID = s.create_non_poly_entity(lig);
	s.create_non_poly(entityID, std::move(atoms));

	// Convert this mmCIF to PDB format
	std::ostringstream os;
	cif::pdb::write(os, restraintFile);

	// Create an RDKit molecule for this PDB file
	auto mol = RDKit::v2::FileParsers::MolFromPDBBlock(os.str());

	RDKit::DGeomHelpers::EmbedMolecule(*mol);

	if (atomIDs.size() != mol->getNumAtoms())
		throw std::runtime_error(std::format("Number of atoms not the same in restraint and RDKit molecule: {} vs {}",
			atomIDs.size(), mol->getNumAtoms()));

	// for each conformer, write out a restraint

	cif::datablock db_1("comp_list");
	db_1["chem_comp"].emplace({ //
		{ "id", lig },
		{ "three_letter_code", lig },
		{ "name", cif::item_value_type::INAPPLICABLE },
		{ "group", "NON-POLYMER" },
		{ "number_atoms_all", mol->getNumAtoms() },
		{ "number_atoms_nh", mol->getNumHeavyAtoms() } });

	std::vector<cif::file> restraints;

	// write original as unmodified first
	{
		cif::datablock db_2("comp_" + lig);

		auto conformer = mol->getConformer(0);

		for (unsigned int i = 0; i < mol->getNumAtoms(); ++i)
		{
			auto atom = mol->getAtomWithIdx(i);
			auto pos = conformer.getAtomPos(i);

			auto symbol = atom->getSymbol();
			auto atomID = atomIDs[i];

			db_2["chem_comp_atom"].emplace({ { "comp_id", lig },
				{ "atom_id", atomID },
				{ "alt_atom_id", atomID },
				{ "type_symbol", symbol },
				{ "type_energy", cif::item_value_type::MISSING },
				{ "charge", atom->getFormalCharge() },
				{ "x", { pos.x, 3 } },
				{ "y", { pos.y, 3 } },
				{ "z", { pos.z, 3 } },
				{ "pdbx_model_Cartn_x_ideal", { pos.x, 3 } },
				{ "pdbx_model_Cartn_y_ideal", { pos.y, 3 } },
				{ "pdbx_model_Cartn_z_ideal", { pos.z, 3 } },
				{ "pdbx_leaving_atom_flag", "N" } });
		}

		for (unsigned int i = 0; i < mol->getNumBonds(); ++i)
		{
			auto bond = mol->getBondWithIdx(i);

			auto a1 = bond->getBeginAtom();
			auto a2 = bond->getEndAtom();

			auto na1 = atomIDs[a1->getIdx()];
			auto na2 = atomIDs[a2->getIdx()];

			auto p1 = conformer.getAtomPos(a1->getIdx());
			auto p2 = conformer.getAtomPos(a2->getIdx());

			auto dist = cif::distance(
				cif::point{ static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z) },
				cif::point{ static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z) });

			std::string type;
			if (bond->getBondType() == RDKit::Bond::BondType::SINGLE)
				type = "SINGLE";
			else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE)
				type = "DOUBLE";
			else if (bond->getBondType() == RDKit::Bond::BondType::TRIPLE)
				type = "TRIPLE";
			else if (bond->getBondType() == RDKit::Bond::BondType::AROMATIC)
				type = "AROMATIC";
			else
				type = "UNKNOWN";

			db_2["chem_comp_bond"].emplace({
				{ "comp_id", lig },
				{ "atom_id_1", na1 },
				{ "atom_id_2", na2 },
				{ "value_order", type },
				{ "pdbx_aromatic_flag", bond->getIsAromatic() ? "Y" : "N" },
				{ "value_dist_nucleus", cif::item_value_type::INAPPLICABLE },
				{ "value_dist_nucleus_esd", cif::item_value_type::INAPPLICABLE },
				{ "value_dist", { dist, 3 } },
				{ "value_dist_esd", { 0.01, 4 } },
			});
		}

		for (unsigned int ix1 = 0; ix1 + 1 < mol->getNumBonds(); ++ix1)
		{
			auto b1 = mol->getBondWithIdx(ix1);
			for (unsigned int ix2 = ix1 + 1; ix2 < mol->getNumBonds(); ++ix2)
			{
				auto b2 = mol->getBondWithIdx(ix2);

				unsigned int i, j, k;

				if (b1->getBeginAtomIdx() == b2->getBeginAtomIdx() or b1->getBeginAtomIdx() == b2->getEndAtomIdx())
					j = b1->getBeginAtomIdx();
				else if (b1->getEndAtomIdx() == b2->getBeginAtomIdx() or b1->getEndAtomIdx() == b2->getEndAtomIdx())
					j = b1->getEndAtomIdx();
				else
					continue;

				i = b1->getOtherAtomIdx(j);
				k = b2->getOtherAtomIdx(j);

				auto angle = MolTransforms::getAngleDeg(conformer, i, j, k);

				db_2["chem_comp_anlge"].emplace({
					{ "comp_id", lig },
					{ "atom_id_1", atomIDs[i] },
					{ "atom_id_2", atomIDs[j] },
					{ "atom_id_3", atomIDs[k] },
					{ "value_angle", { angle, 3 } },
					{ "value_angle_esd", { 1.5, 2 } },
				});
			}
		}

		cif::file file;
		file.emplace_back(db_1);
		file.emplace_back(std::move(db_2));

		if (cif::VERBOSE > 1)
			file.save(std::cout);

		restraints.emplace_back(std::move(file));
	}

	if (N > 1)
	{

		RDKit::DGeomHelpers::EmbedParameters params(RDKit::DGeomHelpers::ETKDG);
		params.randomSeed = seed;
		RDKit::DGeomHelpers::EmbedMolecule(*mol, params);

		RDKit::INT_VECT mol_cids;
		RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, mol_cids, N - 1, params);

		for (auto ci : mol_cids)
		{
			cif::datablock db_2("comp_" + lig);

			auto conformer = mol->getConformer(ci);

			RDKit::UFF::UFFOptimizeMolecule(*mol, 1000, 10.0, ci);

			for (unsigned int i = 0; i < mol->getNumAtoms(); ++i)
			{
				auto atom = mol->getAtomWithIdx(i);
				auto pos = conformer.getAtomPos(i);

				auto symbol = atom->getSymbol();
				auto atomID = atomIDs[i];

				db_2["chem_comp_atom"].emplace({ { "comp_id", lig },
					{ "atom_id", atomID },
					{ "alt_atom_id", atomID },
					{ "type_symbol", symbol },
					{ "type_energy", cif::item_value_type::MISSING },
					{ "charge", atom->getFormalCharge() },
					{ "x", { pos.x, 3 } },
					{ "y", { pos.y, 3 } },
					{ "z", { pos.z, 3 } },
					{ "pdbx_model_Cartn_x_ideal", { pos.x, 3 } },
					{ "pdbx_model_Cartn_y_ideal", { pos.y, 3 } },
					{ "pdbx_model_Cartn_z_ideal", { pos.z, 3 } },
					{ "pdbx_leaving_atom_flag", "N" } });
			}

			for (unsigned int i = 0; i < mol->getNumBonds(); ++i)
			{
				auto bond = mol->getBondWithIdx(i);

				auto a1 = bond->getBeginAtom();
				auto a2 = bond->getEndAtom();

				auto na1 = atomIDs[a1->getIdx()];
				auto na2 = atomIDs[a2->getIdx()];

				auto p1 = conformer.getAtomPos(a1->getIdx());
				auto p2 = conformer.getAtomPos(a2->getIdx());

				auto dist = cif::distance(
					cif::point{ static_cast<float>(p1.x), static_cast<float>(p1.y), static_cast<float>(p1.z) },
					cif::point{ static_cast<float>(p2.x), static_cast<float>(p2.y), static_cast<float>(p2.z) });

				std::string type;
				if (bond->getBondType() == RDKit::Bond::BondType::SINGLE)
					type = "SINGLE";
				else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE)
					type = "DOUBLE";
				else if (bond->getBondType() == RDKit::Bond::BondType::TRIPLE)
					type = "TRIPLE";
				else if (bond->getBondType() == RDKit::Bond::BondType::AROMATIC)
					type = "AROMATIC";
				else
					type = "UNKNOWN";

				db_2["chem_comp_bond"].emplace({
					{ "comp_id", lig },
					{ "atom_id_1", na1 },
					{ "atom_id_2", na2 },
					{ "value_order", type },
					{ "pdbx_aromatic_flag", bond->getIsAromatic() ? "Y" : "N" },
					{ "value_dist_nucleus", cif::item_value_type::INAPPLICABLE },
					{ "value_dist_nucleus_esd", cif::item_value_type::INAPPLICABLE },
					{ "value_dist", { dist, 3 } },
					{ "value_dist_esd", { 0.01, 4 } },
				});
			}

			for (unsigned int ix1 = 0; ix1 + 1 < mol->getNumBonds(); ++ix1)
			{
				auto b1 = mol->getBondWithIdx(ix1);
				for (unsigned int ix2 = ix1 + 1; ix2 < mol->getNumBonds(); ++ix2)
				{
					auto b2 = mol->getBondWithIdx(ix2);

					unsigned int i, j, k;

					if (b1->getBeginAtomIdx() == b2->getBeginAtomIdx() or b1->getBeginAtomIdx() == b2->getEndAtomIdx())
						j = b1->getBeginAtomIdx();
					else if (b1->getEndAtomIdx() == b2->getBeginAtomIdx() or b1->getEndAtomIdx() == b2->getEndAtomIdx())
						j = b1->getEndAtomIdx();
					else
						continue;

					i = b1->getOtherAtomIdx(j);
					k = b2->getOtherAtomIdx(j);

					auto angle = MolTransforms::getAngleDeg(conformer, i, j, k);

					db_2["chem_comp_anlge"].emplace({
						{ "comp_id", lig },
						{ "atom_id_1", atomIDs[i] },
						{ "atom_id_2", atomIDs[j] },
						{ "atom_id_3", atomIDs[k] },
						{ "value_angle", { angle, 3 } },
						{ "value_angle_esd", { 1.5, 2 } },
					});
				}
			}

			cif::file file;
			file.emplace_back(db_1);
			file.emplace_back(std::move(db_2));

			if (cif::VERBOSE > 1)
				file.save(std::cout);

			restraints.emplace_back(std::move(file));
		}
	}

	return restraints;
}

} // namespace pdb_redo