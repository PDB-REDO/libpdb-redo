// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

#include "pdb-redo/ClipperWrapper.hpp"

namespace pdb_redo
{

// --------------------------------------------------------------------

clipper::Atom toClipper(cif::const_row_handle atom, cif::const_row_handle aniso_row)
{
	const double kPI = std::numbers::pi;

	clipper::Atom result;

	cif::point location = atom.get<float, float, float>("Cartn_x", "Cartn_y", "Cartn_z");
	result.set_coord_orth({ location.m_x, location.m_y, location.m_z });

	if (atom["occupancy"].empty())
		result.set_occupancy(1.0);
	else
		result.set_occupancy(atom["occupancy"].get<float>());

	std::string element = atom["type_symbol"].get<std::string>();
	if (not atom["pdbx_formal_charge"].empty())
	{
		int charge = atom["pdbx_formal_charge"].get<int>();

		if (cif::atom_type_traits(element).has_sf(charge) and charge != 0)
		{
			element += std::to_string(std::abs(charge));

			if (charge < 0)
				element += '-';
			else
				element += '+';
		}
	}

	result.set_element(element);

	if (not atom["U_iso_or_equiv"].empty())
		result.set_u_iso(atom["U_iso_or_equiv"].get<float>());
	else if (not atom["B_iso_or_equiv"].empty())
		result.set_u_iso(atom["B_iso_or_equiv"].get<float>() / (8 * kPI * kPI));
	else
		throw std::runtime_error("Missing B_iso or U_iso");

	if (aniso_row.empty())
		result.set_u_aniso_orth(clipper::U_aniso_orth(nan("0"), 0, 0, 0, 0, 0));
	else
	{
		const auto &[u11, u12, u13, u22, u23, u33] =
			aniso_row.get<float, float, float, float, float, float>("U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]", "U[2][3]", "U[3][3]");
		result.set_u_aniso_orth(clipper::U_aniso_orth(u11, u22, u33, u12, u13, u23));
	}

	return result;
}

clipper::Atom toClipper(const cif::mm::atom &atom)
{
	return toClipper(atom.get_row(), atom.get_row_aniso());
}

// --------------------------------------------------------------------

clipper::Spacegroup getSpacegroup(const cif::datablock &db)
{
	auto spacegroup = db["symmetry"].find_first<std::string>(cif::key("entry_id") == db.name(), "space_group_name_H-M");

	if (spacegroup == "P 1-")
		spacegroup = "P -1";
	else if (spacegroup == "P 21 21 2 A")
		spacegroup = "P 21 21 2 (a)";
	else if (spacegroup.empty())
		throw std::runtime_error("No spacegroup, cannot continue");

	try
	{
		return clipper::Spacegroup{ clipper::Spgr_descr(cif::get_space_group_number(spacegroup)) };
	}
	catch (const clipper::Message_fatal &m)
	{
		try
		{
			return clipper::Spacegroup{ clipper::Spgr_descr(spacegroup) };
		}
		catch (const clipper::Message_fatal &e)
		{
			std::cerr << e.text() << '\n';
		}
	}

	throw std::runtime_error("Unsupported spacegroup: " + spacegroup);
}

clipper::Cell getCell(const cif::datablock &db)
{
	const auto &[a, b, c, alpha, beta, gamma] =
		db["cell"].find1<float, float, float, float, float, float>(cif::key("entry_id") == db.name(),
			"length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma");

	return clipper::Cell{ clipper::Cell_descr(a, b, c, alpha, beta, gamma) };
}

// --------------------------------------------------------------------

cif::symop_data GetSymOpDataForRTop_frac(const clipper::RTop_frac &rt)
{
	auto &rot = rt.rot();
	auto &trn = rt.trn();

	auto rte = [&rot](int i, int j)
	{ return static_cast<int8_t>(lrint(rot(i, j))); };

	std::array<int, 15> krt{
		rte(0, 0), rte(0, 1), rte(0, 2),
		rte(1, 0), rte(1, 1), rte(1, 2),
		rte(2, 0), rte(2, 1), rte(2, 2)};

	for (int i = 0; i < 3; ++i)
	{
		int n = lrint(trn[i] * 24);
		int d = 24;

		if (n == 0 or std::abs(n) == 24)
			continue; // is 0, 0 in our table

		for (int j = 5; j > 1; --j)
			if (n % j == 0 and d % j == 0)
			{
				n /= j;
				d /= j;
			}

		n = (n + d) % d;

		switch (i)
		{
			case 0:
				krt[9] = n;
				krt[10] = d;
				break;
			case 1:
				krt[11] = n;
				krt[12] = d;
				break;
			case 2:
				krt[13] = n;
				krt[14] = d;
				break;
			default:;
		}
	}

	return cif::symop_data{ krt };
}

// --------------------------------------------------------------------

// std::ostream &operator<<(std::ostream &os, const cif::symop_data &s)
// {
// 	os << '[';

// 	bool first = true;
// 	for (auto i : s.data())
// 	{
// 		if (not std::exchange(first, false))
// 			os << ", ";
// 		os << i;
// 	}
	
// 	os << ']';

// 	return os;
// }

int getSpacegroupNumber(const clipper::Spacegroup &sg)
{
	std::set<cif::symop_data> sg_ops;

	for (int i = 0; i < sg.num_symops(); ++i)
	{
		const auto &symop = sg.symop(i);

		for (int u : {-1, 0, 1})
			for (int v : {-1, 0, 1})
				for (int w : {-1, 0, 1})
				{
					if (i == 0 and u == 0 and v == 0 and w == 0)
						continue;

					auto rtop = clipper::RTop_frac(
						symop.rot(), symop.trn() + clipper::Vec3<>(u, v, w));
					
					sg_ops.insert(GetSymOpDataForRTop_frac(rtop));
				}
	}

	auto s = cif::kSymopNrTable, e = s + cif::kSymopNrTableSize;
	int result = 0;

	while (s != e)
	{
		auto t = s + 1;
		while (t->spacegroup() == s->spacegroup())
			++t;
		
		if (static_cast<std::size_t>(t - s) != sg_ops.size())
		{
			s = t;
			continue;
		}

		std::size_t seen = 0;

		for (auto &k : sg_ops)
		{
			if (std::find_if(s, t, [&k](const cif::symop_datablock &b) { return b.symop() == k; }) != t)
				++seen;
		}

		if (seen != sg_ops.size())
		{
			s = t;
			continue;
		}

		result = s->spacegroup();
		break;
	}

	return result;
}

} // namespace pdb_redo