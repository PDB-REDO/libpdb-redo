// SPDX-FileCopyrightText: NKI/AVL, Netherlands Cancer Institute, 2020
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <cif++/cif++.hpp>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace pdb_redo
{

class BondMapException : public std::runtime_error
{
  public:
	BondMapException(const std::string &msg)
		: runtime_error(msg)
	{
	}
};

class BondMap
{
  public:
	BondMap(const cif::datablock &db, std::optional<std::tuple<cif::point, float>> around = {}, int model_nr = 1);

	BondMap(const BondMap &) = delete;
	BondMap &operator=(const BondMap &) = delete;

	BondMap(BondMap &&);
	BondMap &operator=(BondMap &&);

	bool operator()(const std::string &atom_1, const std::string &atom_2) const
	{
		auto aix1 = index.find(atom_1);
		auto aix2 = index.find(atom_2);
		return aix1 != index.end() and aix2 != index.end() and isBonded(aix1->second, aix2->second);
	}

	bool operator()(const cif::mm::atom &atom_1, const cif::mm::atom &atom_2) const
	{
		return operator()(atom_1.id(), atom_2.id());
	}

	bool is1_4(const std::string &atom_1, const std::string &atom_2) const
	{
		auto aix1 = index.find(atom_1);
		auto aix2 = index.find(atom_2);
		return aix1 != index.end() and aix2 != index.end() and bond_1_4.count(key(aix1->second, aix2->second));
	}

	bool is1_4(const cif::mm::atom &atom_1, const cif::mm::atom &atom_2) const
	{
		return is1_4(atom_1.id(), atom_2.id());
	}

	// links coming from the struct_conn records:
	std::vector<std::string> linked(const std::string &atom) const;

	// This list of atomID's is comming from either CCD or the CCP4 dictionaries loaded
	static std::vector<std::string> atomIDsForCompound(const std::string &compoundID);

	// iterator

	class iterator
	{
	  public:
		using iterator_category = std::bidirectional_iterator_tag;
		using value_type = std::tuple<std::string, std::string>;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		iterator() = default;
		iterator(const iterator &) = default;
		iterator &operator=(const iterator &) = default;

		auto operator*()
		{
			updateValue();
			return mValue;
		}

		auto operator*() const
		{
			updateValue();
			return mValue;
		}

		auto operator->()
		{
			updateValue();
			return &mValue;
		}

		auto operator->() const
		{
			updateValue();
			return &mValue;
		}

		auto &operator++()
		{
			++mCurrent;
			return *this;
		}

		auto operator++(int)
		{
			iterator result(*this);
			this->operator++();
			return result;
		}

		auto &operator--()
		{
			--mCurrent;
			return *this;
		}

		auto operator--(int)
		{
			iterator result(*this);
			this->operator--();
			return result;
		}

		constexpr auto operator==(const iterator &i) const noexcept
		{
			return mCurrent == i.mCurrent;
		}

		constexpr auto operator!=(const iterator &i) const noexcept
		{
			return mCurrent != i.mCurrent;
		}

	  private:
		friend class BondMap;
		using base_iterator = std::set<std::tuple<uint32_t, uint32_t>>::iterator;

		iterator(const BondMap *bm, base_iterator i)
			: mBondmap(bm)
			, mCurrent(i)
		{
		}

		void updateValue() const
		{
			auto [one, two] = *mCurrent;
			mValue = { mBondmap->rIndex[one], mBondmap->rIndex[two] };
		}

		const BondMap *mBondmap;
		base_iterator mCurrent;
		mutable value_type mValue;
	};

	auto begin() const { return iterator(this, bond.begin()); }
	auto end() const { return iterator(this, bond.end()); }

	auto cbegin() const { return std::make_const_iterator(iterator(this, bond.begin())); }
	auto cend() const { return std::make_const_iterator(iterator(this, bond.end())); }

	auto rbegin() const { return std::make_reverse_iterator(iterator(this, bond.end())); }
	auto rend() const { return std::make_reverse_iterator(iterator(this, bond.begin())); }

	size_t size() const { return bond.size(); }

  private:
	friend class iterator;

	constexpr std::tuple<uint32_t, uint32_t> key(uint32_t a, uint32_t b) const
	{
		if (a > b)
			std::swap(a, b);
		return { a, b };
	}

	constexpr bool isBonded(uint32_t ai, uint32_t bi) const
	{
		if (ai > bi)
			std::swap(ai, bi);

		return bond.count({ ai, bi }) != 0;
	}

	uint32_t dim;
	std::unordered_map<std::string, uint32_t> index;
	std::vector<std::string> rIndex;
	std::set<std::tuple<uint32_t, uint32_t>> bond, bond_1_4;
	std::map<std::string, std::set<std::string>> link;
};

} // namespace pdb_redo
