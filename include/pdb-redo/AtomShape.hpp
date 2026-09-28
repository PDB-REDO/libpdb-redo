// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

// AtomShape, analogue to the similarly named code in clipper

#pragma once

#include <cif++/cif++.hpp>

namespace pdb_redo
{

// --------------------------------------------------------------------
// Class used in calculating radii

class AtomShape
{
  public:
	AtomShape(cif::const_row_handle atom, cif::const_row_handle atom_aniso, float resHigh, float resLow,
		bool electronScattering, std::optional<float> bFactor = {});

	AtomShape(const cif::mm::atom &atom, float resHigh, float resLow, bool electronScattering, std::optional<float> bFactor = {})
		: AtomShape(atom.get_row(), atom.get_row_aniso(), resHigh, resLow, electronScattering, bFactor)
	{
	}

	AtomShape(const AtomShape &) = default;
	AtomShape(AtomShape &&rhs)
	{
		swap(*this, rhs);
	}

	AtomShape &operator=(AtomShape rhs)
	{
		swap(*this, rhs);
		return *this;
	}

	~AtomShape();

	friend void swap(AtomShape &a, AtomShape &b) noexcept
	{
		std::swap(a.mImpl, b.mImpl);
	}


	[[nodiscard]] float radius() const;
	[[nodiscard]] float calculatedDensity(float r) const;
	[[nodiscard]] float calculatedDensity(cif::point p) const;

  private:
	std::shared_ptr<struct AtomShapeImpl> mImpl;
};

} // namespace pdb_redo
