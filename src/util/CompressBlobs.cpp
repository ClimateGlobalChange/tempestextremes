///////////////////////////////////////////////////////////////////////////////
///
///	\file    CompressBlobs.cpp
///	\author  Paul Ullrich
///	\version October 16, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "SimpleGrid.h"
#include "NetCDFUtilities.h"

#include <set>
#include <queue>
#include <fstream>

///////////////////////////////////////////////////////////////////////////////

void Compress(
	const std::string & strInputData,
	const std::string & strConnectivity,
	const std::string & strBlobVar,
	const std::string & strOutputData,
	bool fRegional,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName
) {
	AnnounceBanner();

	// Load in the blob file
	NcFile infile(strInputData.c_str());
	if (!infile.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\"", strInputData.c_str());
	}

	// Check for blob variable
	NcVar * varBlob = infile.get_var(strBlobVar.c_str());
	if (varBlob == NULL) {
		_EXCEPTION2("No variable \"%s\" found in file \"%s\"",
			strBlobVar.c_str(), strInputData.c_str());
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		grid.GenerateLatitudeLongitude(
			&infile,
			strLatitudeName,
			strLongitudeName,
			fRegional,
			false);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Number of variable and grid dimensions
	long lVarDims = varBlob->num_dims();
	long lGridDims = grid.DimCount();
	long lAuxDims = lVarDims - lGridDims;
	if (lAuxDims < 0) {
		_EXCEPTION2("Variable \"%s\" must have at least %lu dimensions",
			strBlobVar.c_str(), lGridDims);
	}

	// Get aux dimensions
	long lAllAuxDimCount = 1;
	std::vector<long> lAuxDimSize(lAuxDims);
	for (long l = 0; l < lAuxDims; l++) {
		lAuxDimSize[l] = varBlob->get_dim(l)->size();
		lAllAuxDimCount *= lAuxDimSize[l];
	}

	// Open output file
	std::ofstream ofile(strOutputData.c_str(), std::ios::out | std::ios::binary);
	if (!ofile.is_open()) {
		_EXCEPTION1("Unable to open output file \"%s\"", strOutputData.c_str());
	}

	// Data storage
	DataArray1D<int> iBlobData(grid.GetSize());

	// Loop over all aux dimensions
	std::vector<long> lPos(lVarDims, 0);
	std::vector<long> lSize(lVarDims, 1);

	lSize[lVarDims-2] = varBlob->get_dim(lVarDims-2)->size();
	lSize[lVarDims-1] = varBlob->get_dim(lVarDims-1)->size();

	for (long l = 0; l < lAllAuxDimCount; l++) {

		// File position
		long lx = l;
		for (long d = lAuxDims-1; d >= 0; d--) {
			lPos[d] = lx % lAuxDimSize[d];
			lx = lx / lAuxDimSize[d];
		}

		// Read in this data
		varBlob->set_cur(&(lPos[0]));
		varBlob->get(&(iBlobData[0]), &(lSize[0]));

		// Set of points in blobs
		std::set<size_t> setAllBlobPoints;

		for (size_t s = 0; s < grid.GetSize(); s++) {

			// Find all points that are blobs but unvisited
			if (iBlobData[s] == 0) {
				continue;
			}
			if (setAllBlobPoints.find(s) != setAllBlobPoints.end()) {
				continue;
			}

			// Queue of points to visit
			std::queue<size_t> queueToVisit;
			queueToVisit.push(s);

			// Set of all perimeter points in blob
			std::set<size_t> setPerimeter;

			// Set of all interior points in blob
			std::set<size_t> setNonPerimeter;

			// Build up the sets of perimeter and interior points in blob
			size_t sCurrent = s;

			while (queueToVisit.size() != 0) {
				size_t sNext = queueToVisit.front();
				queueToVisit.pop();

				if (setPerimeter.find(sNext) != setPerimeter.end()) {
					continue;
				}
				if (setNonPerimeter.find(sNext) != setNonPerimeter.end()) {
					continue;
				}

				// Add neighbors to queue and identify if this point is a
				// perimeter point or not
				bool fPerimeterPoint = false;
				for (size_t n = 0; n < grid.m_vecConnectivity[sNext].size(); n++) {
					size_t sNeighbor = grid.m_vecConnectivity[sNext][n];
					if (iBlobData[sNeighbor] == 0) {
						fPerimeterPoint = true;
					} else {
						queueToVisit.push(sNeighbor);
					}
				}

				if (fPerimeterPoint) {
					setPerimeter.insert(sNext);
					setAllBlobPoints.insert(sNext);
				} else {
					setNonPerimeter.insert(sNext);
					setAllBlobPoints.insert(sNext);
				}
			}

			// Find number of interior points needed
			std::set<size_t> setNeededNonPerimeter;

			while (setNonPerimeter.size() != 0) {
				size_t sInterior = *(setNonPerimeter.begin());
				setNeededNonPerimeter.insert(sInterior);

				queueToVisit.push(sInterior);

				std::set<size_t> setVisited;

				while (queueToVisit.size() != 0) {
					size_t sNext = queueToVisit.front();
					queueToVisit.pop();

					setNonPerimeter.erase(sNext);

					if (setPerimeter.find(sNext) != setPerimeter.end()) {
						continue;
					}
					if (setVisited.find(sNext) != setVisited.end()) {
						continue;
					}

					setVisited.insert(sNext);

					for (size_t n = 0; n < grid.m_vecConnectivity[sNext].size(); n++) {
						queueToVisit.push(grid.m_vecConnectivity[sNext][n]);
					}
				}
			}

			// Write perimeter points
			std::cout << setPerimeter.size() << " " << setNonPerimeter.size() << " " << setNeededNonPerimeter.size() << std::endl;

			ofile << setPerimeter.size();
			for (auto it = setPerimeter.begin(); it != setPerimeter.end(); it++) {
				ofile << *it;
			}

			ofile << setNeededNonPerimeter.size();
			for (auto it = setNeededNonPerimeter.begin(); it != setNeededNonPerimeter.end(); it++) {
				ofile << *it;
			}

		}
	}

	AnnounceBanner();

}

///////////////////////////////////////////////////////////////////////////////

void Compress2(
	const std::string & strInputData,
	const std::string & strConnectivity,
	const std::string & strBlobVar,
	const std::string & strOutputData,
	bool fRegional,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName
) {
	AnnounceBanner();

	// Load in the blob file
	NcFile infile(strInputData.c_str());
	if (!infile.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\"", strInputData.c_str());
	}

	// Check for blob variable
	NcVar * varBlob = infile.get_var(strBlobVar.c_str());
	if (varBlob == NULL) {
		_EXCEPTION2("No variable \"%s\" found in file \"%s\"",
			strBlobVar.c_str(), strInputData.c_str());
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		grid.GenerateLatitudeLongitude(
			&infile,
			strLatitudeName,
			strLongitudeName,
			fRegional,
			false);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Number of variable and grid dimensions
	long lVarDims = varBlob->num_dims();
	long lGridDims = grid.DimCount();
	long lAuxDims = lVarDims - lGridDims;
	if (lAuxDims < 0) {
		_EXCEPTION2("Variable \"%s\" must have at least %lu dimensions",
			strBlobVar.c_str(), lGridDims);
	}

	// Get aux dimensions
	long lAllAuxDimCount = 1;
	std::vector<long> lAuxDimSize(lAuxDims);
	for (long l = 0; l < lAuxDims; l++) {
		lAuxDimSize[l] = varBlob->get_dim(l)->size();
		lAllAuxDimCount *= lAuxDimSize[l];
	}

	// Open output file
	std::ofstream ofile(strOutputData.c_str(), std::ios::out | std::ios::binary);
	if (!ofile.is_open()) {
		_EXCEPTION1("Unable to open output file \"%s\"", strOutputData.c_str());
	}

	// Data storage
	DataArray1D<int> iBlobData(grid.GetSize());

	// Loop over all aux dimensions
	std::vector<long> lPos(lVarDims, 0);
	std::vector<long> lSize(lVarDims, 1);

	lSize[lVarDims-2] = varBlob->get_dim(lVarDims-2)->size();
	lSize[lVarDims-1] = varBlob->get_dim(lVarDims-1)->size();

	for (long l = 0; l < lAllAuxDimCount; l++) {

		// File position
		long lx = l;
		for (long d = lAuxDims-1; d >= 0; d--) {
			lPos[d] = lx % lAuxDimSize[d];
			lx = lx / lAuxDimSize[d];
		}

		// Read in this data
		varBlob->set_cur(&(lPos[0]));
		varBlob->get(&(iBlobData[0]), &(lSize[0]));

		// Write counts and value
		long lCount = 1;
		int iCurrent = iBlobData[0];
		for (long l = 1; l <= grid.GetSize(); l++) {
			if (l == grid.GetSize()) {
				ofile << lCount;
			}
			if (iBlobData[l] != iCurrent) {
				ofile << lCount;
				lCount = 0;
				iCurrent = iBlobData[l];
			}
			lCount++;
		}
	}

	AnnounceBanner();
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input data file
	std::string strInputData;

	// Connectivity file
	std::string strConnectivity;

	// Blob variable
	std::string strBlobVar;

	// Output data file
	std::string strOutputData;

	// Decompress
	bool fDecompress;

	// Regional data
	bool fRegional;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strBlobVar, "var", "");
		CommandLineString(strOutputData, "out", "");
		CommandLineBool(fDecompress, "decompress");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Check arguments
	if (strInputData == "") {
		_EXCEPTIONT("No input file (--in) specified");
	}
	if (strOutputData == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}
	if (strBlobVar == "") {
		_EXCEPTIONT("No blob variable (--var) specified");
	}

	// Decompress
	if (fDecompress) {
		_EXCEPTION();
	} else {
		Compress2(
			strInputData,
			strConnectivity,
			strBlobVar,
			strOutputData,
			fRegional,
			strLatitudeName,
			strLongitudeName
		); 
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

