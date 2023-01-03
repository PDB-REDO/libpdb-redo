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

#include <filesystem>
#include <fstream>
#include <iomanip>

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <cif++.hpp>

#include "pdb-redo/ClipperWrapper.hpp"
#include "pdb-redo/MapMaker.hpp"
#include "pdb-redo/ResolutionCalculator.hpp"
#include "pdb-redo/Statistics.hpp"
#include "pdb-redo/Symmetry-2.hpp"

#ifdef _MSC_VER
#include <io.h>
#define mkstemp _mktemp
#endif

namespace fs = std::filesystem;

namespace pdb_redo
{

// --------------------------------------------------------------------
// a private ccp4 map file implementation

//  1      NC              # of Columns    (fastest changing in map)
//  2      NR              # of Rows
//  3      NS              # of Sections   (slowest changing in map)
//  4      MODE            Data type
//                           0 = envelope stored as signed bytes (from
//                               -128 lowest to 127 highest)
//                           1 = Image     stored as Integer*2
//                           2 = Image     stored as Reals
//                           3 = Transform stored as Complex Integer*2
//                           4 = Transform stored as Complex Reals
//                           5 == 0

//                           Note: Mode 2 is the normal mode used in
//                                 the CCP4 programs. Other modes than 2 and 0
//                                 may NOT WORK

//  5      NCSTART         Number of first COLUMN  in map
//  6      NRSTART         Number of first ROW     in map
//  7      NSSTART         Number of first SECTION in map
//  8      NX              Number of intervals along X
//  9      NY              Number of intervals along Y
// 10      NZ              Number of intervals along Z
// 11      X length        Cell Dimensions (Angstroms)
// 12      Y length                     "
// 13      Z length                     "
// 14      Alpha           Cell Angles     (Degrees)
// 15      Beta                         "
// 16      Gamma                        "
// 17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
// 18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
// 19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
// 20      AMIN            Minimum density value
// 21      AMAX            Maximum density value
// 22      AMEAN           Mean    density value    (Average)
// 23      ISPG            Space group number
// 24      NSYMBT          Number of bytes used for storing symmetry operators
// 25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
// 26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
//                         LSKFLG .ne. 0.
// 35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
//                         Skew transformation is from standard orthogonal
//                         coordinate frame (as used for atoms) to orthogonal
//                         map frame, as

//                                 Xo(map) = S * (Xo(atoms) - t)

// 38      future use       (some of these are used by the MSUBSX routines
//  .          "              in MAPBRICK, MAPCONT and FRODO)
//  .          "   (all set to zero by default)
//  .          "
// 52          "

// 53	MAP	        Character string 'MAP ' to identify file type
// 54	MACHST		Machine stamp indicating the machine type
// 			which wrote file
// 55      ARMS            Rms deviation of map from mean density
// 56      NLABL           Number of labels being used
// 57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)

enum CCP4MapFileMode : uint32_t
{
	AS_REALS = 2 // do not support anything else for now...
};

struct CCP4MapFileHeader
{
	uint32_t NC, NR, NS;
	CCP4MapFileMode MODE;
	int32_t NCSTART, NRSTART, NSSTART;
	uint32_t NX, NY, NZ;
	float cellLengths[3];
	float cellAngles[3];
	uint32_t MAPC, MAPR, MAPS;
	float AMIN, AMAX, AMEAN;
	uint32_t ISPG;
	uint32_t NSYMBT;
	uint32_t LSKFLG;
	float SKWMAT[9];
	float SKWTRN[3];
	uint32_t UNUSED[15];
	char MAP[4] = {'M', 'A', 'P', ' '};
	uint32_t MACHST = 0x00004144;
	float ARMS;
	uint32_t NLABL = 1;
	char LABEL[200 * 4];
};

template <typename FTYPE>
std::tuple<FTYPE, FTYPE, FTYPE, FTYPE> CalculateMapStatistics(const clipper::Xmap<FTYPE> &xmap, clipper::Grid_range r)
{
	FTYPE
	amin = std::numeric_limits<FTYPE>::max(),
	amax = std::numeric_limits<FTYPE>::min();
	long double asum = 0, asum2 = 0;
	size_t n = 0;

	clipper::Xmap_base::Map_reference_coord c(xmap);
	for (int g0 = r.min()[0]; g0 <= r.max()[0]; ++g0)
		for (int g1 = r.min()[1]; g1 <= r.max()[1]; ++g1)
			for (int g2 = r.min()[2]; g2 <= r.max()[2]; ++g2)
			{
				c.set_coord({g0, g1, g2});
				FTYPE v = xmap[c];

				asum += v;
				asum2 += v * v;
				if (amin > v)
					amin = v;
				if (amax < v)
					amax = v;

				++n;
			}

	FTYPE mean = static_cast<FTYPE>(asum / n);
	FTYPE rmsd = static_cast<FTYPE>(std::sqrt((asum2 / n) - (mean * mean)));

	return std::make_tuple(amin, amax, mean, rmsd);
}

template <typename FTYPE>
void writeCCP4MapFile(std::ostream &os, clipper::Xmap<FTYPE> &xmap, clipper::Grid_range range)
{
	static_assert(sizeof(CCP4MapFileHeader) == 256 * 4, "Map header is of incorrect size");
	// static_assert(__BYTE_ORDER == __LITTLE_ENDIAN, "Code for big endian systems is not implemented yet");

	auto &spacegroup = xmap.spacegroup();
	int spaceGroupNumber = spacegroup.descr().spacegroup_number();
	int orderFMS[3] = {3, 1, 2};

	switch (spaceGroupNumber)
	{
		case 1:
		case 2:
		case 3:
		case 4:
		case 10:
		case 16:
		case 17:
		case 18:
		case 20:
		case 21:
		case 23:
			orderFMS[0] = 2;
			orderFMS[2] = 3;
			break;
	}

	int orderXYZ[3];
	for (int i = 0; i < 3; ++i)
		orderXYZ[orderFMS[i] - 1] = i;

	int grid[3], gridFMSMin[3], gridFMSMax[3], dim[3];
	for (int i = 0; i < 3; ++i)
	{
		grid[i] = xmap.grid_sampling()[i];

		gridFMSMin[orderXYZ[i]] = range.min()[i];
		gridFMSMax[orderXYZ[i]] = range.max()[i];
	}

	for (size_t i = 0; i < 3; ++i)
		dim[i] = gridFMSMax[i] - gridFMSMin[i] + 1;

	auto cellDescription = xmap.cell().descr();

	CCP4MapFileHeader h = {};

	int r = snprintf(h.LABEL, sizeof(h.LABEL), "%s", "Map created with map-maker from the PDB-REDO suite of applications");
	for (size_t i = r; i < sizeof(h.LABEL); ++i)
		h.LABEL[i] = ' ';

	h.NC = dim[0];
	h.NR = dim[1];
	h.NS = dim[2];
	h.MODE = AS_REALS;
	h.NCSTART = gridFMSMin[0];
	h.NRSTART = gridFMSMin[1];
	h.NSSTART = gridFMSMin[2];
	h.NX = grid[0];
	h.NY = grid[1];
	h.NZ = grid[2];
	h.cellLengths[0] = cellDescription.a();
	h.cellLengths[1] = cellDescription.b();
	h.cellLengths[2] = cellDescription.c();
	h.cellAngles[0] = cellDescription.alpha_deg();
	h.cellAngles[1] = cellDescription.beta_deg();
	h.cellAngles[2] = cellDescription.gamma_deg();
	h.MAPC = orderFMS[0];
	h.MAPR = orderFMS[1];
	h.MAPS = orderFMS[2];
	h.ISPG = spaceGroupNumber;
	h.NSYMBT = spacegroup.num_symops() * 80;

	std::tie(h.AMIN, h.AMAX, h.AMEAN, h.ARMS) = CalculateMapStatistics(xmap, range);

	os.write(reinterpret_cast<char *>(&h), sizeof(h));

	const std::string kSpaces(80, ' ');
	for (int si = 0; si < spacegroup.num_symops(); ++si)
	{
		std::string symop = spacegroup.symop(si).format();
		os.write(symop.c_str(), symop.length());
		os.write(kSpaces.c_str(), 80 - symop.length());
	}

	clipper::Xmap_base::Map_reference_coord c(xmap);
	const uint32_t kSectionLength = dim[0] * dim[1];
	std::vector<float> section(kSectionLength);

	int g[3];
	for (g[2] = gridFMSMin[2]; g[2] <= gridFMSMax[2]; ++g[2])
	{
		auto si = section.begin();
		for (g[1] = gridFMSMin[1]; g[1] <= gridFMSMax[1]; ++g[1])
			for (g[0] = gridFMSMin[0]; g[0] <= gridFMSMax[0]; ++g[0])
			{
				c.set_coord({g[orderXYZ[0]], g[orderXYZ[1]], g[orderXYZ[2]]});
				*si++ = static_cast<float>(xmap[c]);
			}

		assert(si == section.end());
		os.write(reinterpret_cast<char *>(section.data()), kSectionLength * sizeof(float));
	}
}

// --------------------------------------------------------------------

bool IsMTZFile(const fs::path &p)
{
	bool result = false;

	std::ifstream f(p);
	if (f.is_open())
	{
		char sig[5] = {};
		f.read(sig, 4);
		result = sig == std::string("MTZ ");
	}

	return result;
}

// --------------------------------------------------------------------

template <typename FTYPE>
Map<FTYPE>::Map()
{
}

template <typename FTYPE>
Map<FTYPE>::~Map()
{
}

template <typename FTYPE>
void Map<FTYPE>::calculateStats()
{
	double sum = 0, sum2 = 0;
	int count = 0;

	mMinDensity = std::numeric_limits<double>::max();
	mMaxDensity = std::numeric_limits<double>::min();

	for (auto ix = mMap.first(); not ix.last(); ix.next())
	{
		auto v = mMap[ix];

		if (std::isnan(v))
			throw std::runtime_error("map contains NaN values");

		if (mMinDensity > v)
			mMinDensity = v;
		if (mMaxDensity < v)
			mMaxDensity = v;

		++count;
		sum += v;
		sum2 += v * v;
	}

	mMeanDensity = sum / count;
	mRMSDensity = std::sqrt((sum2 / count) - (mMeanDensity * mMeanDensity));
}

template <typename FTYPE>
void Map<FTYPE>::read(const std::filesystem::path &f)
{
	fs::path mapFile(f);
	fs::path dataFile = mapFile;

	if (cif::VERBOSE > 0)
		std::cout << "Reading map from " << mapFile << std::endl;

	if (mapFile.extension() == ".gz")
	{
		// file is compressed

		fs::path p = mapFile.parent_path();
		std::string s = mapFile.filename().string();

		cif::gzio::ifstream in(mapFile);

		char tmpFileName[] = "/tmp/map-tmp-XXXXXX";
		if (mkstemp(tmpFileName) < 0)
			throw std::runtime_error(std::string("Could not create temp file for map: ") + strerror(errno));

		dataFile = fs::path(tmpFileName);
		std::ofstream out(dataFile);

		if (not in.is_open() or not out.is_open())
			throw std::runtime_error("Could not handle compressed map file");

		out << in.rdbuf();
	}

	if (not fs::exists(dataFile))
		throw std::runtime_error("Could not open map file " + mapFile.string());

	using namespace clipper;

	CCP4MAPfile mapin;
	mapin.open_read(dataFile.string());
	mapin.import_xmap(mMap);
	mapin.close_read();

	if (dataFile != mapFile)
		fs::remove(dataFile);

	calculateStats();
}

template <typename FTYPE>
void Map<FTYPE>::write(const std::filesystem::path &f)
{
	write_masked(f, mMap.grid_asu());
}

template <typename FTYPE>
void Map<FTYPE>::write_masked(std::ostream &os, clipper::Grid_range r)
{
	writeCCP4MapFile(os, mMap, r);
}

template <typename FTYPE>
void Map<FTYPE>::write_masked(const std::filesystem::path &f, clipper::Grid_range r)
{
	std::ofstream file(f, std::ios_base::binary);
	if (not file.is_open())
		throw std::runtime_error("Could not open map file for writing: " + f.string());

	write_masked(file, r);
}

// --------------------------------------------------------------------

template <typename FTYPE>
Map<FTYPE> Map<FTYPE>::masked(const cif::mm::structure &structure, const std::vector<cif::mm::atom> &atoms) const
{
	Map<FTYPE> result(*this);

	auto rtops = AlternativeSites(getSpacegroup(structure.get_datablock()), getCell(structure.get_datablock()));

	for (auto &atom : atoms)
	{
		float radius = cif::atom_type_traits(atom.get_type()).radius(cif::radius_type::van_der_waals);
		if (std::isnan(radius))
			radius = cif::atom_type_traits(atom.get_type()).radius(cif::radius_type::calculated);
		if (std::isnan(radius))	// TODO: now what?
			radius = 200;

		auto cloc = toClipper(atom.get_location());

		for (auto &rt : rtops)
		{
			auto rcloc = cloc.transform(rt);

			iterateGrid(toClipper(toPoint(rcloc)), radius, result.mMap,
				[&result, radiusSq = radius * radius, a = atom.get_location()](auto iw)
				{
					cif::point p = toPoint(iw.coord_orth());

					if (distance_squared(a, p) < radiusSq)
						result.mMap[iw] = -10;
				});
		}
	}

	return result;
}

template <typename FTYPE>
float Map<FTYPE>::z_weighted_density(const cif::mm::structure &structure, const std::vector<cif::mm::atom> &atoms) const
{
	FTYPE result = 0;

	for (auto &atom : atoms)
	{
		auto co = toClipper(atom.get_location());
		auto a_cf = co.coord_frac(mMap.cell());
		auto a_cm = a_cf.coord_map(mMap.grid_sampling());

		FTYPE dv;
		clipper::Interp_nearest::interp(mMap, a_cm, dv);

		result += dv * static_cast<int>(atom.get_type() == cif::atom_type::D ? cif::atom_type::H : atom.get_type());
	}

	return result;
}

// --------------------------------------------------------------------

template class Map<float>;
template class Map<double>;

// --------------------------------------------------------------------

template <typename FTYPE>
MapMaker<FTYPE>::MapMaker()
{
}

template <typename FTYPE>
MapMaker<FTYPE>::~MapMaker()
{
}

template <typename FTYPE>
void MapMaker<FTYPE>::loadMTZ(const fs::path &f, float samplingRate,
	std::initializer_list<std::string> fbLabels, std::initializer_list<std::string> fdLabels,
	std::initializer_list<std::string> foLabels, std::initializer_list<std::string> fcLabels,
	std::initializer_list<std::string> faLabels)
{
	fs::path hklin(f);

	if (cif::VERBOSE > 0)
		std::cerr << "Reading map from " << hklin << std::endl
				  << "  with labels: FB: " << cif::join(fbLabels, ",") << std::endl
				  << "  with labels: FD: " << cif::join(fdLabels, ",") << std::endl
				  << "  with labels: FA: " << cif::join(faLabels, ",") << std::endl
				  << "  with labels: FO: " << cif::join(foLabels, ",") << std::endl
				  << "  with labels: FC: " << cif::join(fcLabels, ",") << std::endl;

	fs::path dataFile = hklin;

	if (hklin.extension() == ".gz")
	{
		// file is compressed

		fs::path p = hklin.parent_path();
		std::string s = hklin.filename().string();

		cif::gzio::ifstream in(hklin);

		char tmpFileName[] = "/tmp/mtz-tmp-XXXXXX";
		if (mkstemp(tmpFileName) < 0)
			throw std::runtime_error(std::string("Could not create temp file for mtz: ") + strerror(errno));

		dataFile = fs::path(tmpFileName);
		std::ofstream out(dataFile);
		
		out << in.rdbuf();
	}

	if (not fs::exists(dataFile))
		throw std::runtime_error("Could not open mtz file " + hklin.string());

	const std::string kBasePath("/%s/%s/[%s]");

	using clipper::CCP4MTZfile;

	CCP4MTZfile mtzin;
	mtzin.open_read(dataFile.string());

	mtzin.import_hkl_info(mHKLInfo);

	bool hasFAN = false, hasFREE = false;
	const std::regex rx(R"(^/[^/]+/[^/]+/(.+) \S$)");

	for (auto &label : mtzin.column_labels())
	{
		std::smatch m;
		if (not std::regex_match(label, m, rx))
			continue;

		if (m[1] == "FAN")
		{
			hasFAN = true;
			continue;
		}

		if (m[1] == "FREE")
		{
			hasFREE = true;
			continue;
		}
	}

	mtzin.import_hkl_data(mFbData,
		cif::format(kBasePath, "*", "*", cif::join(fbLabels, ",")).str());
	mtzin.import_hkl_data(mFdData,
		cif::format(kBasePath, "*", "*", cif::join(fdLabels, ",")).str());
	if (hasFAN)
		mtzin.import_hkl_data(mFaData,
			cif::format(kBasePath, "*", "*", cif::join(faLabels, ",")).str());
	mtzin.import_hkl_data(mFoData,
		cif::format(kBasePath, "*", "*", cif::join(foLabels, ",")).str());
	mtzin.import_hkl_data(mFcData,
		cif::format(kBasePath, "*", "*", cif::join(fcLabels, ",")).str());

	if (hasFREE)
		mtzin.import_hkl_data(mFreeData,
			cif::format(kBasePath, "*", "*", "FREE").str());

	mtzin.import_hkl_data(mPhiFomData,
		cif::format(kBasePath, "*", "*", "PHWT,FOM").str());

	mtzin.close_read();

	if (dataFile != hklin)
		fs::remove(dataFile);

	Cell cell = mHKLInfo.cell();
	Spacegroup spacegroup = mHKLInfo.spacegroup();

	ResolutionCalculator rc(cell);
	mResHigh = 99;
	mResLow = 0;

	for (auto hi = mFoData.first_data(); not hi.last(); hi = mFoData.next_data(hi))
	{
		auto res = rc(hi.hkl().h(), hi.hkl().k(), hi.hkl().l());

		if (mResHigh > res)
			mResHigh = res;

		if (mResLow < res)
			mResLow = res;
	}

	if (mResLow == 0 and mResHigh == 99)
		throw std::runtime_error("Empty Fo map");

	//	fixMTZ();

	mGrid.init(spacegroup, cell,
		mHKLInfo.resolution(), samplingRate); // define grid

	clipper::Xmap<FTYPE> &fbMap = mFb;
	clipper::Xmap<FTYPE> &fdMap = mFd;
	clipper::Xmap<FTYPE> &faMap = mFa;

	fbMap.init(spacegroup, cell, mGrid); // define map
	fbMap.fft_from(mFbData);             // generate map

	fdMap.init(spacegroup, cell, mGrid); // define map
	fdMap.fft_from(mFdData);             // generate map

	if (not mFaData.is_null())
	{
		faMap.init(spacegroup, cell, mGrid);
		faMap.fft_from(mFaData);
	}

	if (cif::VERBOSE > 0)
	{
		std::cerr << "Read Xmaps with sampling rate: " << samplingRate << std::endl
				  << "  stored resolution: " << mHKLInfo.resolution().limit() << std::endl
				  << "  calculated reshi = " << mResHigh << " reslo = " << mResLow << std::endl
				  << "  spacegroup: " << spacegroup.symbol_hm() << std::endl
				  << "  cell: " << cell.format() << std::endl
				  << "  grid: " << mGrid.format() << std::endl;

		printStats();
	}

	mFb.calculateStats();
	mFd.calculateStats();
}

// --------------------------------------------------------------------

template <typename FTYPE>
void MapMaker<FTYPE>::loadMaps(const fs::path &fbMapFile, const fs::path &fdMapFile, float reshi, float reslo)
{
	mResHigh = reshi;
	mResLow = reslo;

	mFb.read(fbMapFile);
	mFd.read(fdMapFile);

	if (not mFb.cell().equals(mFd.cell()))
		throw std::runtime_error("Fb and Fd map do not contain the same cell");

	clipper::Resolution reso(reshi);

	mHKLInfo.init(mFb.spacegroup(), mFb.cell(), reso, true);
	mGrid = mFb.get().grid_sampling();
}

// --------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const clipper::HKL &hkl)
{
	os << "h: " << hkl.h() << ", "
	   << "k: " << hkl.k() << ", "
	   << "l: " << hkl.l();

	return os;
};

// --------------------------------------------------------------------

template <typename FTYPE>
void MapMaker<FTYPE>::calculate(const fs::path &hklin,
	const cif::mm::structure &structure, bool noBulk, AnisoScalingFlag anisoScaling,
	float samplingRate, bool electronScattering,
	std::initializer_list<std::string> foLabels, std::initializer_list<std::string> freeLabels)
{
	if (IsMTZFile(hklin))
		loadFoFreeFromMTZFile(hklin, foLabels, freeLabels);
	else
		loadFoFreeFromReflectionsFile(hklin);

	recalc(structure, noBulk, anisoScaling, samplingRate, electronScattering);
}

// --------------------------------------------------------------------

template <typename FTYPE>
void MapMaker<FTYPE>::loadFoFreeFromReflectionsFile(const fs::path &hklin)
{
	using clipper::HKL;

	cif::file reflnsFile(hklin);
	auto &reflns = reflnsFile.front();

	//	m_xname = reflns["exptl_crystal"].front()["id"].as<std::string>();
	//	m_pname = reflns["entry"].front()["id"].as<std::string>();

	float a, b, c, alpha, beta, gamma;
	cif::tie(a, b, c, alpha, beta, gamma) = reflns["cell"].front().get(
		"length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma");

	using clipper::Cell_descr;
	Cell cell = Cell(Cell_descr{a, b, c, alpha, beta, gamma});

	//	if (not cell2.equals(m_cell))
	//		throw std::runtime_error("Reflections file and coordinates file do not agree upon the cell parameters");

	// --------------------------------------------------------------------

	// Read reflections file to calculate resolution low and high
	ResolutionCalculator rc(a, b, c, alpha, beta, gamma);
	double hires = 99;

	for (auto r : reflns["refln"])
	{
		int h, k, l;

		cif::tie(h, k, l) = r.get("index_h", "index_k", "index_l");

		double res = rc(h, k, l);

		if (hires > res)
			hires = res;
	}

	std::string spacegroupDescr = reflns["symmetry"].front()["space_group_name_H-M"].as<std::string>();
	auto spacegroup = Spacegroup(clipper::Spgr_descr{spacegroupDescr});
	mHKLInfo = HKL_info(spacegroup, cell, clipper::Resolution{hires}, true);

	//	m_crystal = MTZcrystal(m_xname, m_pname, m_cell);

	mFoData.init(mHKLInfo, mHKLInfo.cell());
	mFreeData.init(mHKLInfo, mHKLInfo.cell());

	for (auto ih = mFreeData.first(); not ih.last(); ih.next())
		mFreeData[ih].set_null();

	// --------------------------------------------------------------------

	enum FreeRConvention
	{
		frXPLO,
		frCCP4
	} freeRConvention = frXPLO;
	int freeRefl = 1, workRefl = 0;

	if (false /*m_statusXPLO*/)
	{
		freeRConvention = frCCP4;
		freeRefl = 0;
		workRefl = 1;
	}

	bool first = false;
	for (auto r : reflns["refln"])
	{
		int h, k, l;
		char flag;
		float F, sigF;

		cif::tie(h, k, l, flag, F, sigF) = r.get("index_h", "index_k", "index_l", "status", "F_meas_au", "F_meas_sigma_au");

		int ix = mHKLInfo.index_of(HKL{h, k, l});

		if (ix < 0)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Ignoring hkl(" << h << ", " << k << ", " << l << ")" << std::endl;
			continue;
		}

		if (first and (flag == freeRefl or flag == workRefl))
		{
			std::cerr << "Non-standard _refln.status column detected" << std::endl
					  << "Assuming " << (freeRConvention == frXPLO ? "XPLOR" : "CCP4") << " convention for free R flag" << std::endl;
			first = false;
		}

		mFoData[ix] = F_sigF(F, sigF);

		switch (flag)
		{
			case 'o':
			case 'h':
			case 'l':
				mFreeData[ix] = Flag(1);
				break;

			case 'f':
				mFreeData[ix] = Flag(0);
				break;

			case '0':
			case '1':
				mFreeData[ix] = Flag(workRefl == flag ? 1 : 0);
				break;

			default:
				if (cif::VERBOSE > 1)
					std::cerr << "Unexpected value in status: '" << flag << "' for hkl(" << h << ", " << k << ", " << l << ")" << std::endl;
				break;
		}
	}
}

// --------------------------------------------------------------------

template <typename FTYPE>
void MapMaker<FTYPE>::loadFoFreeFromMTZFile(const fs::path &hklin,
	std::initializer_list<std::string> foLabels, std::initializer_list<std::string> freeLabels)
{
	if (cif::VERBOSE > 0)
		std::cerr << "Recalculating maps from " << hklin << std::endl;

	const std::string kBasePath("/%s/%s/[%s]");

	using clipper::CCP4MTZfile;

	CCP4MTZfile mtzin;
	mtzin.open_read(hklin.string());

	mtzin.import_hkl_info(mHKLInfo);
	mtzin.import_hkl_data(mFoData,
		cif::format(kBasePath, "*", "*", cif::join(foLabels, ",")).str());
	mtzin.import_hkl_data(mFreeData,
		cif::format(kBasePath, "*", "*", cif::join(freeLabels, ",")).str());

	mtzin.close_read();
}

// --------------------------------------------------------------------

template <typename FTYPE>
void MapMaker<FTYPE>::recalc(const cif::mm::structure &structure,
	bool noBulk, AnisoScalingFlag anisoScaling,
	float samplingRate, bool electronScattering)
{
	Cell cell = mHKLInfo.cell();
	Spacegroup spacegroup = mHKLInfo.spacegroup();

	// The calculation work
	std::vector<clipper::Atom> atoms;

	for (auto a : structure.atoms())
		atoms.push_back(toClipper(a));

	mFcData.init(mHKLInfo, cell);

	if (not electronScattering)
	{
		auto &exptl = structure.get_category("exptl");
		electronScattering = not exptl.empty() and exptl.front()["method"] == "ELECTRON CRYSTALLOGRAPHY";
	}

	clipper::ScatteringFactors::selectScattteringFactorsType(
		electronScattering ? clipper::SF_ELECTRON : clipper::SF_WAASMAIER_KIRFEL);

	if (noBulk)
	{
		clipper::SFcalc_aniso_fft<float> sfc;
		sfc(mFcData, atoms);
	}
	else
	{
		clipper::SFcalc_obs_bulk<float> sfcb;
		sfcb(mFcData, mFoData, atoms);

		if (cif::VERBOSE > 0)
			std::cerr << "Bulk correction volume: " << sfcb.bulk_frac() << std::endl
					  << "Bulk correction factor: " << sfcb.bulk_scale() << std::endl;
	}

	if (anisoScaling != as_None)
	{
		clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
		clipper::SFscale_aniso<float> sfscl;
		if (anisoScaling == as_Observed)
			sfscl(mFoData, mFcData); // scale Fobs
		else
			sfscl(mFcData, mFoData); // scale Fcal

		if (cif::VERBOSE > 0)
			std::cerr << "Anisotropic scaling:" << std::endl
					  << sfscl.u_aniso_orth(F).format() << std::endl;
	}

	// now do sigmaa calc
	mFbData.init(mHKLInfo, cell);
	mFdData.init(mHKLInfo, cell);
	mPhiFomData.init(mHKLInfo, cell);

	HKL_data<Flag> flag(mHKLInfo, cell);

	const int freeflag = 0;
	for (auto ih = mFreeData.first(); not ih.last(); ih.next())
	{
		if (not mFoData[ih].missing() and (mFreeData[ih].missing() or mFreeData[ih].flag() == freeflag))
			flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
		else
			flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	}

	// do sigmaa calc
	clipper::SFweight_spline<float> sfw(mNumRefln, mNumParam);
	sfw(mFbData, mFdData, mPhiFomData, mFoData, mFcData, flag);

	// mFbData now contains 2mFo - DFc
	// mFdData now contains  mFo - DFc

	fixMTZ();

	ResolutionCalculator rc(cell);
	mResHigh = 99;
	mResLow = 0;

	for (auto hi = mFoData.first_data(); not hi.last(); hi = mFoData.next_data(hi))
	{
		auto res = rc(hi.hkl().h(), hi.hkl().k(), hi.hkl().l());

		if (mResHigh > res)
			mResHigh = res;

		if (mResLow < res)
			mResLow = res;
	}

	if (cif::VERBOSE > 1)
		std::cerr << "calculated reshi = " << mResHigh << " reslo = " << mResLow << std::endl;

	// samplingRate /= 2;

	mGrid.init(spacegroup, cell,
		mHKLInfo.resolution(), samplingRate); // define grid

	clipper::Xmap<FTYPE> &fbMap = mFb;
	clipper::Xmap<FTYPE> &fdMap = mFd;

	fbMap.init(spacegroup, cell, mGrid); // define map
	fbMap.fft_from(mFbData);             // generate map

	fdMap.init(spacegroup, cell, mGrid); // define map
	fdMap.fft_from(mFdData);             // generate map

	if (cif::VERBOSE > 0)
	{
		std::cerr << "Read Xmaps with sampling rate: " << samplingRate << std::endl
				  << "  resolution: " << mResHigh << std::endl
				  << "  cell: " << cell.format() << std::endl
				  << "  grid: " << mGrid.format() << std::endl;

		printStats();
	}

	mFb.calculateStats();
	mFd.calculateStats();
}

template <typename FTYPE>
void MapMaker<FTYPE>::fixMTZ()
{
	Spacegroup spacegroup = mHKLInfo.spacegroup();

	enum
	{
		A1,  // A1:  FC = 2mFo - FM
		A2,  // A2:  FC >= 2mFo - FM
		A3,  // A3:  FD = FM - mFo
		A4,  // A4:  FD = 2(FM - mFo)
		C5,  // C5:  FC = 2mFo - FM
		C6,  // C6:  FM = mFo
		C7,  // C7:  FD = mFo - FC
		C8,  // C8:  FD = 2(mFo - FC)
		C9,  // C9:  FD <= mFo - FC
		T10, // 10:  FM = FC (unobserved only)
		T11, // 11:  FD = 0 (unobserved only)
		TestCount
	};

	std::vector<bool> tests(TestCount, true);

	// first run the tests to see if we need to fix anything

	if (cif::VERBOSE > 0)
		std::cerr << "Testing MTZ file" << std::endl;

	for (auto ih = mFbData.first(); not ih.last(); ih.next())
	{
		clipper::HKL_class cls(spacegroup, ih.hkl());

		auto W = mPhiFomData[ih].fom();

		auto FM = mFbData[ih].f();
		auto PM = mFbData[ih].phi() * 180 / kPI;
		auto FD = mFdData[ih].f();
		auto PD = mFdData[ih].phi() * 180 / kPI;
		auto FO = mFoData[ih].f();
		auto FC = mFcData[ih].f();
		auto PC = mFcData[ih].phi() * 180 / kPI;

		auto WFO = W * FO;

		if (std::abs(std::fmod(std::abs(PM - PC) + 180, 360) - 180) > 90)
			FM = -FM;

		if (std::abs(std::fmod(std::abs(PD - PC) + 180, 360) - 180) > 90)
			FD = -FD;

		if (mFoData[ih].missing() or W == 0)
		{
			if (tests[T10] and std::abs(FM - FC) > 0.05)
			{
				tests[T10] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test 10 failed at " << ih.hkl() << std::endl;
			}

			if (tests[T11] and std::abs(FD) > 0.05)
			{
				tests[T11] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test 11 failed at " << ih.hkl() << std::endl;
			}
		}
		else if (cls.centric())
		{
			if (tests[C5] and std::abs(FC + FM - 2 * WFO) > 0.05)
			{
				tests[C5] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test C5 failed at " << ih.hkl() << std::endl;
			}

			if (tests[C6] and std::abs(FM - WFO) > 0.05)
			{
				tests[C6] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test C6 failed at " << ih.hkl() << std::endl;
			}

			if (tests[C7] and std::abs(FC + FD - WFO) > 0.05)
			{
				tests[C7] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test C7 failed at " << ih.hkl() << std::endl;
			}

			if (tests[C8] and std::abs(FC + 0.5 * FD - WFO) > 0.05)
			{
				tests[C8] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test C8 failed at " << ih.hkl() << std::endl;
			}

			if (tests[C9] and (1.01 * FC + FD - WFO) < -0.05)
			{
				tests[C9] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test C9 failed at " << ih.hkl() << std::endl;
			}
		}
		else
		{
			if (tests[A1] and std::abs(FC + FM - 2 * WFO) > 0.05)
			{
				tests[A1] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test A1 failed at " << ih.hkl() << std::endl;
			}

			if (tests[A2] and 1.01 * FC + FM - 2 * WFO < -0.05)
			{
				tests[A2] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test A2 failed at " << ih.hkl() << std::endl;
			}

			if (tests[A3] and std::abs(FM - FD - WFO) > 0.05)
			{
				tests[A3] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test A3 failed at " << ih.hkl() << std::endl;
			}

			if (tests[A4] and std::abs(FM - 0.5 * FD - WFO) > 0.05)
			{
				tests[A4] = false;
				if (cif::VERBOSE > 0)
					std::cerr << "Test A4 failed at " << ih.hkl() << std::endl;
			}
		}
	}

	using clipper::HKL_class;
	using clipper::data32::F_phi;

	const F_phi fzero(0, 0);

	// mtzfix...
	for (auto ih = mFbData.first(); not ih.last(); ih.next())
	{
		if (mFbData[ih].missing() or mFdData[ih].missing())
			continue;

		auto PM = mFbData[ih].phi() * 180 / kPI;
		auto PD = mFdData[ih].phi() * 180 / kPI;
		auto PC = mFcData[ih].phi() * 180 / kPI;

		if (std::abs(std::fmod(std::abs(PM - PC) + 180, 360) - 180) > 90)
		{
			mFbData[ih].f() = -mFbData[ih].f();
			mFbData[ih].phi() = mFcData[ih].phi();
		}

		if (std::abs(std::fmod(std::abs(PD - PC) + 180, 360) - 180) > 90)
		{
			mFdData[ih].f() = -mFdData[ih].f();
			mFdData[ih].phi() = mFcData[ih].phi();
		}

		auto mFo = mFbData[ih] - mFdData[ih];

		HKL_class cls(spacegroup, ih.hkl());

		if (not mFoData[ih].missing() and mPhiFomData[ih].fom() > 0)
		{
			if (cls.centric())
			{
				if (not tests[C6])
					mFbData[ih] = mFo;
				if (not tests[C7] and tests[C8])
					mFdData[ih].f() = mFdData[ih].f() / 2;
			}
			else
			{
				if (tests[A3] and not tests[A4])
					mFdData[ih] = mFdData[ih] + mFdData[ih];
			}
		}
		else
		{
			if (not tests[T10])
			{
				if ((not cls.centric() and tests[A1]) or
					(cls.centric() and (tests[C5] or tests[C7] or tests[C8])))
				{
					mFbData[ih] = mFcData[ih];
				}
			}

			if (not tests[T11])
				mFdData[ih] = fzero;
		}
	}
}

template <typename FTYPE>
void MapMaker<FTYPE>::printStats()
{
	// calc R and R-free
	std::vector<double> params(mNumParam, 1.0);

	clipper::BasisFn_spline basisfn(mFoData, mNumParam, 1.0);
	clipper::TargetFn_scaleF1F2<clipper::data32::F_phi, clipper::data32::F_sigF> targetfn(mFcData, mFoData);
	clipper::ResolutionFn rfn(mHKLInfo, basisfn, targetfn, params);

	double r1w = 0, f1w = 0, r1f = 0, f1f = 0;
	const int freeflag = 0;

	for (auto ih = mFoData.first_data(); not ih.last(); ih = mFoData.next_data(ih))
	{
		if (mFcData[ih].missing())
			continue;
		//			throw std::runtime_error("missing Fc");

		double Fo = mFoData[ih].f();
		double Fc = std::sqrt(rfn.f(ih)) * mFcData[ih].f();

		if (mFreeData[ih].flag() == freeflag)
		{
			r1f += fabs(Fo - Fc);
			f1f += Fo;
		}
		else
		{
			r1w += fabs(Fo - Fc);
			f1w += Fo;
		}
	}

	if (f1f < 0.1)
		f1f = 0.1;
	r1f /= f1f;

	if (f1w < 0.1)
		f1w = 0.1;
	r1w /= f1w;

	std::cerr << "R-factor      : " << r1w << std::endl
			  << "Free R-factor : " << r1f << std::endl;
}

template <typename FTYPE>
void MapMaker<FTYPE>::writeMTZ(const fs::path &file, const std::string &pname, const std::string &cname)
{
	if (mHKLInfo.is_null())
		throw std::runtime_error("HKL info not initialized");

	clipper::CCP4MTZfile mtz;
	clipper::MTZdataset dataset(pname, 0);
	clipper::MTZcrystal crystal(cname, pname, mHKLInfo.cell());

	const std::string col = "/" + pname + "/" + cname + "/";

	mtz.open_write(file.string());
	mtz.export_hkl_info(mHKLInfo);
	mtz.export_crystal(crystal, col);
	mtz.export_dataset(dataset, col);
	if (not mFreeData.is_null())
		mtz.export_hkl_data(mFreeData, col + "[FREE]");
	if (not mFoData.is_null())
		mtz.export_hkl_data(mFoData, col + "[FP,SIGFP]");
	if (not mFcData.is_null())
		mtz.export_hkl_data(mFcData, col + "[FC_ALL,PHIC_ALL]");
	if (not mFbData.is_null())
		mtz.export_hkl_data(mFbData, col + "[FWT,PHWT]");
	if (not mFdData.is_null())
		mtz.export_hkl_data(mFdData, col + "[DELFWT,PHDELWT]");
	if (not mPhiFomData.is_null())
		mtz.export_hkl_data(mPhiFomData, col + "[PHI,FOM]");
	mtz.close_write();
}

template class MapMaker<float>;
template class MapMaker<double>;

} // namespace pdb_redo
