/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
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

#include "cif++/compound.hpp"
#include "cif++/model.hpp"
#include "pdb-redo/BlobFinder.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/ShapeFitter.hpp"
#include <catch2/matchers/catch_matchers.hpp>

#define CATCH_CONFIG_RUNNER

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cif++/cif++.hpp>
#include <cif++/pdb.hpp>
#include <clipper/core/xmap.h>
#include <filesystem>
#include <glm/glm.hpp>
#include <gsl/gsl_blas.h> // for debugging norm of gradient
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector_double.h>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, std::size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, std::size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

// --------------------------------------------------------------------

std::filesystem::path gTestDir = std::filesystem::current_path();

int main(int argc, char *argv[])
{
	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
#if CATCH22
	using namespace Catch::clara;
#else
	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;
#endif

	auto cli = session.cli()                                // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                  // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                   // the option names it will respond to
	           ("The directory containing the data files"); // description string for the help output

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	if (fs::exists(gTestDir / "minimal-components.cif"))
		cif::compound_factory::instance().push_dictionary(gTestDir / "minimal-components.cif");

	return session.run();
}

// --------------------------------------------------------------------

// TEST_CASE("sf-1")
// {
// 	pdb_redo::Map<float> map;
// 	map.read(gTestDir / "1cbs-REA-blob.map");
// 	clipper::Xmap<float> &xmap = map.get();

// 	auto cf = R"(
// data_1CBS
// #
// _entry.id   1CBS
// #
// _cell.entry_id           1CBS
// _cell.length_a           45.650
// _cell.length_b           47.560
// _cell.length_c           77.610
// _cell.angle_alpha        90.00
// _cell.angle_beta         90.00
// _cell.angle_gamma        90.00
// _cell.Z_PDB              4
// _cell.pdbx_unique_axis   ?
// #
// _symmetry.entry_id                         1CBS
// _symmetry.space_group_name_H-M             'P 21 21 21'
// _symmetry.pdbx_full_space_group_name_H-M   ?
// _symmetry.cell_setting                     ?
// _symmetry.Int_Tables_number                19
// #
// 	)"_cf;

// 	// Create a ligand
// 	cif::datablock &db = cf.front(); // almost empty
// 	db.set_validator(cif::validator_factory::instance().get("mmcif_pdbx.dic"));
// 	cif::mm::structure s(db);

// 	pdb_redo::BlobFinder bf(xmap, 0);

// 	auto blob = bf.next();

// 	CHECK(blob.size() == 831);

// 	auto ligand_asym_id = s.create_non_poly("REA", true);

// 	std::cout << "ligand: " << ligand_asym_id << " created\n";

// 	auto score = pdb_redo::fitShape(s, ligand_asym_id, xmap, blob);

// 	CHECK(score < 0);

// 	// std::ofstream file(std::filesystem::temp_directory_path() / "test.cif");
// 	// cf.save(file);
// }

// --------------------------------------------------------------------

TEST_CASE("sf-2")
{
	const fs::path example(gTestDir / ".." / "examples" / "1cbs.cif.gz");
	cif::file file(example.string());
	cif::mm::structure s(file);
	s.remove_residue(s.get_residue("B"));

	pdb_redo::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ(gTestDir / ".." / "examples" / "1cbs_map.mtz", samplingRate);

	auto &mm_fb = mm.fb();
	auto maskedmap = mm_fb.masked(s, s.atoms());

	pdb_redo::BlobFinder blobFinder(maskedmap, s);

	auto ligand_asym_id = s.create_non_poly("REA", true);

	for (;;)
	{
		auto blob = blobFinder.next();

		auto score = pdb_redo::fitShape(s, ligand_asym_id, mm_fb, blob);

		CHECK(score < 0);

		// std::ofstream of(std::filesystem::temp_directory_path() / "test-2.cif");
		// file.save(of);

		break;
	}
}

// // --------------------------------------------------------------------

// TEST_CASE("sf-3")
// {
// 	const fs::path example(gTestDir / "3aba_besttls.cif.gz");
// 	cif::file file = cif::pdb::read(example.string());

// 	cif::pdb::reconstruct_pdbx(file);

// 	cif::mm::structure s(file);
// 	// auto &db = s.get_datablock();

// 	for (std::string asymID : { "H", "I", "J", "K", "L"})
// 	{
// 		s.remove_residue(s.get_residue(asymID));

// 		pdb_redo::MapMaker<float> mm;
// 		float samplingRate = 0.75;
// 		mm.loadMTZ(gTestDir / "3aba_loopwhole.mtz", samplingRate);

// 		auto &mm_fb = mm.fb();
// 		auto maskedmap = mm_fb.masked(s, s.atoms());

// 		pdb_redo::BlobFinder blobFinder(maskedmap, s);

// 		auto ligand_asym_id = s.create_non_poly("GOL", true);

// 		for (int i = 0;; ++i)
// 		{
// 			auto blob = blobFinder.next();
// 			if (blob.empty())
// 				break;

// 			auto score = pdb_redo::fitShape(s, ligand_asym_id, mm_fb, blob);

// 			if (score < 0)
// 			{
// 				std::ofstream of(std::filesystem::temp_directory_path() / std::format("{}-{}-{}.cif", "3aba", asymID, i));
// 				file.save(of);
// 			}
// 		}
// 	}
// }

auto createInertiaTensorForBlob(const std::vector<cif::point> pts, clipper::Xmap<float> &xmap)
{
	std::array<float, 6> If{};

	auto [c, r] = cif::smallest_sphere_around_points(pts);

	for (auto pt : pts)
	{
		clipper::Coord_orth cp{ pt.m_x, pt.m_y, pt.m_z };
		clipper::Coord_frac pf = cp.coord_frac(xmap.cell());
		auto dp = xmap.interp<clipper::Interp_cubic>(pf);

		pt -= c;

		If[0] += dp * (pt.m_y * pt.m_y + pt.m_z * pt.m_z); // 11
		If[1] += dp * (pt.m_x * pt.m_x + pt.m_z * pt.m_z); // 22
		If[2] += dp * (pt.m_x * pt.m_x + pt.m_y * pt.m_y); // 33
		If[3] -= dp * pt.m_x * pt.m_y;                     // 12
		If[4] -= dp * pt.m_x * pt.m_z;                     // 13
		If[5] -= dp * pt.m_y * pt.m_z;                     // 23
	}

	return glm::mat3{
		glm::normalize(glm::vec3{ If[0], If[3], If[4] }),
		glm::normalize(glm::vec3{ If[3], If[1], If[5] }),
		glm::normalize(glm::vec3{ If[4], If[5], If[2] })
	};
}

auto createInertiaTensorForLigand(const cif::mm::residue &res)
{
	std::array<float, 6> If{};

	auto [c, r] = cif::smallest_sphere_around_points(
		res.atoms() | std::views::transform(&cif::mm::atom::get_location) | std::ranges::to<std::vector>());

	for (auto atom : res.atoms())
	{
		auto pt = atom.get_location();
		pt -= c;

		cif::atom_type_traits t(atom.get_type());

		auto dp = t.weight();

		If[0] += dp * (pt.m_y * pt.m_y + pt.m_z * pt.m_z); // 11
		If[1] += dp * (pt.m_x * pt.m_x + pt.m_z * pt.m_z); // 22
		If[2] += dp * (pt.m_x * pt.m_x + pt.m_y * pt.m_y); // 33
		If[3] -= dp * pt.m_x * pt.m_y;                     // 12
		If[4] -= dp * pt.m_x * pt.m_z;                     // 13
		If[5] -= dp * pt.m_y * pt.m_z;                     // 23
	}

	return glm::mat3{
		glm::normalize(glm::vec3{ If[0], If[3], If[4] }),
		glm::normalize(glm::vec3{ If[3], If[1], If[5] }),
		glm::normalize(glm::vec3{ If[4], If[5], If[2] })
	};
}

auto principalAxis(const glm::mat3 &m)
{
	// Eigen::Matrix3f M;

	// M(0, 0) = m[0][0];
	// M(1, 1) = m[1][1];
	// M(2, 2) = m[2][2];
	// M(0, 1) = M(1, 0) = m[0][1];
	// M(0, 2) = M(2, 0) = m[0][2];
	// M(1, 2) = M(2, 1) = m[1][2];

	// Eigen::EigenSolver<Eigen::Matrix3f> es(M);

	// auto v = es.eigenvectors()[0];

	// return glm::vec3 { v[0], v[1], v[2] };

	double data[9] = {
		m[0][0], m[0][1], m[0][2],
		m[0][1], m[1][1], m[1][2],
		m[0][2], m[1][2], m[2][2]
	};

	gsl_matrix_view g = gsl_matrix_view_array(data, 3, 3);

	gsl_vector *eval = gsl_vector_alloc(3);
	gsl_matrix *evec = gsl_matrix_alloc(3, 3);

	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv(&g.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);

	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	gsl_vector_view evec_i = gsl_matrix_column(evec, 0);

	cif::point result;

	result.m_x = gsl_vector_get(&evec_i.vector, 0);
	result.m_y = gsl_vector_get(&evec_i.vector, 1);
	result.m_z = gsl_vector_get(&evec_i.vector, 2);

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	result.normalize();

	return result;
}

// --------------------------------------------------------------------

TEST_CASE("i-1")
{
	std::array<float, 6> If{};


	// std::vector<cif::point> pts{
	// 	{ 0, 0, 0 },
	// 	{ 0.5, 0.5, 0.5 },
	// 	{ 1, 1, 1 }
	// };

	// std::vector<cif::point> pts{
	// 	{ 0, 0, 0 },
	// 	{ 0.5, 0.5, 0 },
	// 	{ 1, 1, 0 }
	// };

	std::vector<cif::point> pts{
		{ 0, 0.5, 0 },
		{ 0.5, 0.5, 0 },
		{ 1, 0.5, 0 }
	};



	auto [c, r] = cif::smallest_sphere_around_points(pts);

	for (auto pt : pts)
	{
		pt -= c;

		float dp = 1;

		If[0] += dp * (pt.m_y * pt.m_y + pt.m_z * pt.m_z); // 11
		If[1] += dp * (pt.m_x * pt.m_x + pt.m_z * pt.m_z); // 22
		If[2] += dp * (pt.m_x * pt.m_x + pt.m_y * pt.m_y); // 33
		If[3] -= dp * pt.m_x * pt.m_y;                     // 12 21
		If[4] -= dp * pt.m_x * pt.m_z;                     // 13 31
		If[5] -= dp * pt.m_y * pt.m_z;                     // 23 32
	}

	// glm::mat3 im{
	// 	glm::normalize(glm::vec3{ If[0], If[3], If[4] }),
	// 	glm::normalize(glm::vec3{ If[3], If[1], If[5] }),
	// 	glm::normalize(glm::vec3{ If[4], If[5], If[2] })
	// };

	glm::mat3 im{
		glm::vec3{ If[0], If[3], If[4] },
		glm::vec3{ If[3], If[1], If[5] },
		glm::vec3{ If[4], If[5], If[2] }
	};


	auto v = principalAxis(im);
	CHECK_THAT(v.m_x, Catch::Matchers::WithinAbs(1.0f, 0.1f));
	CHECK_THAT(v.m_y, Catch::Matchers::WithinAbs(1.0f, 0.1f));
	CHECK_THAT(v.m_z, Catch::Matchers::WithinAbs(1.0f, 0.1f));

}