// Simple example showing the use of libpdb-redo to calculate density statistics

#include <filesystem>

#include <pdb-redo/Statistics.hpp>
#include "cif++/BondMap.hpp"

namespace fs = std::filesystem;


int main()
{
	namespace c = mmcif;

	// a sample structure file, lets use 1CBS
	const fs::path example("1cbs.cif.gz");

	// load the mmCIF
	mmcif::File file(example);

	// and load this into a structure (note, structure caches data from the file, so order is important)
	mmcif::Structure structure(file);

	// now create the maps based on the MTZ file
	c::MapMaker<float> mm;
	float samplingRate = 0.75;
	mm.loadMTZ("1cbs_map.mtz", samplingRate);

	// create a map for the bonds in this structure
	mmcif::BondMap bm(structure);

	// and finally collect the statistics
	mmcif::EDIAStatsCollector collector(mm, structure, false, bm);
	auto r = collector.collect();

	for (auto& ri: r)
	{
		// and do something with the data
		std::cout << ri.EDIAm << std::endl;
	}

	return 0;
}
