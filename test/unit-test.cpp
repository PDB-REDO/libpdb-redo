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

#define BOOST_TEST_MODULE Libpdb_redo_Test
#include <boost/test/included/unit_test.hpp>

#include <stdexcept>
#include <filesystem>

#include <cif++/Structure.hpp>

#include "pdb-redo/AtomShape.hpp"


namespace fs = std::filesystem;
namespace tt = boost::test_tools;
namespace utf = boost::unit_test;

// --------------------------------------------------------------------

// atom radii calculated with AtomShape and NEWUOA
struct TestRadius {
	std::string		type;
	float			radius;
} kTestRadii[] = {
	{ "N", 1.07327 },
	{ "C", 1.07747 },
	{ "C", 1.06093 },
	{ "O", 1.03793 },
	{ "C", 1.08041 },
	{ "C", 1.0807 },
	{ "C", 1.09096 },
	{ "N", 1.02885 },
	{ "C", 1.01706 },
	{ "C", 0.98581 },
	{ "O", 0.949879 },
	{ "C", 1.03256 },
	{ "C", 1.04372 },
	{ "O", 1.00052 },
	{ "N", 1.05783 },
	{ "N", 0.957396 },
	{ "C", 0.95347 },
	{ "C", 0.952071 },
	{ "O", 0.926261 },
	{ "C", 0.938776 },
	{ "C", 0.947439 },
	{ "C", 0.940042 },
	{ "C", 0.958542 },
	{ "C", 0.92616 },
	{ "C", 0.946796 },
	{ "C", 0.935361 },
	{ "N", 0.930847 },
	{ "C", 0.943131 },
	{ "C", 0.93617 },
	{ "O", 0.911906 },
	{ "C", 0.96053 },
	{ "O", 0.961552 },
	{ "N", 0.92031 },
	{ "C", 0.939679 },
	{ "C", 0.942494 },
	{ "O", 0.909548 },
	{ "N", 0.924386 },
	{ "C", 0.941311 },
	{ "C", 0.93563 },
	{ "O", 0.919858 },
	{ "C", 0.955999 },
	{ "C", 0.996134 },
	{ "O", 0.982804 },
	{ "N", 1.00745 },
	{ "N", 0.911032 },
	{ "C", 0.923198 },
	{ "C", 0.930114 },
	{ "O", 0.90597 },
	{ "C", 0.917515 },
	{ "C", 0.917601 }
};

BOOST_AUTO_TEST_CASE(atom_shape_1, *utf::tolerance(0.001f))
{
	const fs::path example("/usr/share/doc/libcifpp-dev/examples/1cbs.cif.gz");

	mmcif::File file(example);
	mmcif::Structure structure(file);

	const float kResHi = 1.80009, kResLo = 7.99918; 

	const size_t N = sizeof(kTestRadii) / sizeof(TestRadius);
	size_t i = 0;

	for (auto& atom: structure.atoms())
	{
		if (i >= N)
			break;

		mmcif::AtomShape shape(atom, kResHi, kResLo, false);

		BOOST_CHECK(mmcif::AtomTypeTraits(atom.type()).symbol() == kTestRadii[i].type);

		float radius = shape.radius();
		float test = kTestRadii[i].radius;

		BOOST_TEST(radius == test);

		++i;
	}
}
