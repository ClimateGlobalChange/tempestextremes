///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchNodes.cpp
///	\author  Paul Ullrich
///	\version October 1st, 2014
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
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

#include "DataVector.h"
#include "DataMatrix.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A structure for storing a latitude / longitude index.
///	</summary>
struct LatLonPair {

	int lat;
	int lon;

	///	<summary>
	///		Constructor.
	///	</summary>
	LatLonPair(int a_lat, int a_lon) :
		lat(a_lat), lon(a_lon)
	{ }

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const LatLonPair & pr) const {
		if (lat < pr.lat) {
			return true;
		} else if (lat > pr.lat) {
			return false;
		}

		if (lon < pr.lon) {
			return true;
		} else {
			return false;
		}
	}

	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const LatLonPair & pr) const {
		if ((lat == pr.lat) && (lon == pr.lon)) {
			return true;
		} else {
			return false;
		}
	}
};

///	<summary>
///		A structure for storing a bounding box in latitude / longitude space.
///	</summary>
struct LatLonBox {
	int lat[2];
	int lon[2];
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {

	// Input file
	std::string strInputFile;

	// Output file
	std::string strOutputFile;

	// Variable name
	std::string strVariable;

	// Minimum patch size
	int nMinPatchSize;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariable, "var", "");
		CommandLineInt(nMinPatchSize, "minsize", 1);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if (strInputFile.length() == 0) {
		_EXCEPTIONT("No input file (--in) specified");
	}

	// Check output
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check variable
	if (strVariable.length() == 0) {
		_EXCEPTIONT("No variable name (--var) specified");
	}

	// Load the netcdf file
	NcFile ncInput(strInputFile.c_str());

	// Get latitude/longitude dimensions
	NcDim * dimLat = ncInput.get_dim("lat");
	if (dimLat == NULL) {
		_EXCEPTIONT("No dimension \"lat\" found in input file");
	}

	NcDim * dimLon = ncInput.get_dim("lon");
	if (dimLon == NULL) {
		_EXCEPTIONT("No dimension \"lon\" found in input file");
	}

	NcVar * varLat = ncInput.get_var("lat");
	if (varLat == NULL) {
		_EXCEPTIONT("No variable \"lat\" found in input file");
	}

	NcVar * varLon = ncInput.get_var("lon");
	if (varLon == NULL) {
		_EXCEPTIONT("No variable \"lon\" found in input file");
	}

	int nLat = dimLat->size();
	int nLon = dimLon->size();

	DataVector<double> dataLat(nLat);
	varLat->get(dataLat, nLat);

	for (int j = 0; j < nLat; j++) {
		dataLat[j] *= M_PI / 180.0;
	}

	DataVector<double> dataLon(nLon);
	varLon->get(dataLon, nLon);

	for (int i = 0; i < nLon; i++) {
		dataLon[i] *= M_PI / 180.0;
	}

	// Get time dimension
	NcDim * dimTime = ncInput.get_dim("time");
	if (dimTime == NULL) {
		_EXCEPTIONT("No dimension \"time\" found in input file");
	}

	NcVar * varTime = ncInput.get_var("time");
	if (varTime == NULL) {
		_EXCEPTIONT("No variable \"time\" found in input file");
	}

	int nTime = dimTime->size();

	DataVector<double> dTime;
	dTime.Initialize(nTime);

	varTime->get(dTime, nTime);

	// Load in indicator variable
	NcVar * varIndicator = ncInput.get_var(strVariable.c_str());

	if (varIndicator->num_dims() != 3) {
		_EXCEPTION1("Incorrect number of dimensions for \"%s\" (3 expected)",
			strVariable.c_str());
	}

	// Allocate indicator data
	DataMatrix<int> dataIndicator(nLat, nLon);

	// Set of indicator locations
	typedef std::set<LatLonPair> IndicatorSet;
	typedef IndicatorSet::iterator IndicatorSetIterator;

	IndicatorSet setIndicators;
	IndicatorSet setNeighbors;

	// Set of nodes contained in all patches
	std::vector< std::vector<IndicatorSet> > vecAllPatches;
	vecAllPatches.resize(nTime);

	// Bounding boxes for all patches
	std::vector< std::vector<LatLonBox> > vecAllPatchBoxes;
	vecAllPatchBoxes.resize(nTime);

	// Loop through all times
	for (int t = 0; t < nTime; t++) {

		// Get the current patch vector
		std::vector<IndicatorSet> & vecPatches = vecAllPatches[t];

		std::vector<LatLonBox> & vecPatchBoxes = vecAllPatchBoxes[t];

		// New announcement block for timestep
		char szStartBlock[128];
		sprintf(szStartBlock, "Time %i", t);
		AnnounceStartBlock(szStartBlock);

		// Load in the data at this time
		varIndicator->set_cur(t, 0, 0);
		varIndicator->get(&(dataIndicator[0][0]), 1, nLat, nLon);

		// Insert all detected locations into set
		for (int j = 0; j < nLat; j++) {
		for (int i = 0; i < nLon; i++) {
			if (dataIndicator[j][i] != 0) {
				setIndicators.insert( LatLonPair(j, i));
			}
		}
		}

		Announce("Tagged points: %i", setIndicators.size());

		// Rejections due to insufficient node count
		int nRejectedMinSize = 0;

		// Find all patches
		for (; setIndicators.size() != 0;) {

			// Next starting location
			LatLonPair pr = *(setIndicators.begin());

			// Current patch
			int ixPatch = vecPatches.size();
			vecPatches.resize(ixPatch+1);
			vecPatchBoxes.resize(ixPatch+1);

			// Initialize bounding box
			LatLonBox & box = vecPatchBoxes[ixPatch];

			box.lat[0] = pr.lat;
			box.lat[1] = pr.lat;
			box.lon[0] = pr.lon;
			box.lon[1] = pr.lon;

			// Find all connecting elements in patch
			setNeighbors.clear();
			setNeighbors.insert(pr);
			while (setNeighbors.size() != 0) {
				pr = *(setNeighbors.begin());
				setNeighbors.erase(setNeighbors.begin());

				// This node is already included in the patch
				if (vecPatches[ixPatch].find(pr) != vecPatches[ixPatch].end()) {
					continue;
				}

				// This node has not been tagged
				IndicatorSetIterator iterIndicator = setIndicators.find(pr);
				if (iterIndicator == setIndicators.end()) {
					continue;
				}

				// Remove this indicator from the set of available indicators
				setIndicators.erase(iterIndicator);

				// Insert the node into the patch
				vecPatches[ixPatch].insert(pr);

				// Update bounding box
				if (pr.lat > box.lat[1]) {
					box.lat[1] = pr.lat;
				}
				if (pr.lat < box.lat[0]) {
					box.lat[0] = pr.lat;
				}
				if (pr.lon == nLon-1) {
					if (box.lon[0] <= box.lon[1]) {
						if (nLon - 1 - box.lon[1] < box.lon[0] + 1) {
							box.lon[1] = nLon - 1;
						} else {
							box.lon[0] = nLon - 1;
						}
					}

				} else if (pr.lon == 0) {
					if (box.lon[0] <= box.lon[1]) {
						if (nLon - box.lon[1] < box.lon[0]) {
							box.lon[1] = 0;
						} else {
							box.lon[0] = 0;
						}
					}

				} else {
					if (box.lon[0] <= box.lon[1]) {
						if (pr.lon < box.lon[0]) {
							box.lon[0] = pr.lon;
						}
						if (pr.lon > box.lon[1]) {
							box.lon[1] = pr.lon;
						}

					} else {
						if ((pr.lon >= box.lon[0]) || (pr.lon <= box.lon[1])) {
						} else if (pr.lon - box.lon[1] < box.lon[0] - pr.lon) {
							box.lon[1] = pr.lon;
						} else {
							box.lon[0] = pr.lon;
						}
					}
				}

				// Insert neighbors
				if (pr.lat == 0) {
					for (int i = 0; i < nLon; i++) {
						setNeighbors.insert(LatLonPair(0, i));
					}
					setNeighbors.insert(LatLonPair(1, (pr.lon+nLon-1)%nLon));
					setNeighbors.insert(LatLonPair(1, pr.lon));
					setNeighbors.insert(LatLonPair(1, (pr.lon+1)%nLon));

				} else if (pr.lat == nLat - 1) {
					for (int i = 0; i < nLon; i++) {
						setNeighbors.insert(LatLonPair(nLat-1, i));
					}
					setNeighbors.insert(
						LatLonPair(nLat-2, (pr.lon+nLon-1)%nLon));
					setNeighbors.insert(
						LatLonPair(nLat-2, pr.lon));
					setNeighbors.insert(
						LatLonPair(nLat-2, (pr.lon+1)%nLon));
	
				} else {

					setNeighbors.insert(
						LatLonPair(pr.lat+1, (pr.lon+nLon-1)%nLon));
					setNeighbors.insert(
						LatLonPair(pr.lat  , (pr.lon+nLon-1)%nLon));
					setNeighbors.insert(
						LatLonPair(pr.lat-1, (pr.lon+nLon-1)%nLon));

					setNeighbors.insert(
						LatLonPair(pr.lat+1, pr.lon));
					setNeighbors.insert(
						LatLonPair(pr.lat-1, pr.lon));

					setNeighbors.insert(
						LatLonPair(pr.lat+1, (pr.lon+1)%nLon));
					setNeighbors.insert(
						LatLonPair(pr.lat,   (pr.lon+1)%nLon));
					setNeighbors.insert(
						LatLonPair(pr.lat-1, (pr.lon+1)%nLon));
				}
			}

			// Check patch size
			if (vecPatches[ixPatch].size() < nMinPatchSize) {
				nRejectedMinSize++;
				vecPatches.resize(ixPatch);
				vecPatchBoxes.resize(ixPatch);
			}
		}

		Announce("Patches detected: %i", vecPatches.size());
		Announce("Rejected (min size): %i", nRejectedMinSize);

		for (int p = 0; p < vecPatchBoxes.size(); p++) {
			printf("Patch %i [%i, %i] x [%i, %i]\n",
				p+1,
				vecPatchBoxes[p].lat[0],
				vecPatchBoxes[p].lat[1],
				vecPatchBoxes[p].lon[0],
				vecPatchBoxes[p].lon[1]);
		}

		AnnounceEndBlock("Done");
	}

	ncInput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

