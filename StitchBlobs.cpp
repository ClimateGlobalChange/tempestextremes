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
#include "NetCDFUtilities.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <map>

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

	///	<summary>
	///		Bounding latitudes (endpoints are included).
	///	</summar>
	int lat[2];

	///	<summary>
	///		Bounding longitudes (endpoints are included).
	///	</summary>
	int lon[2];

	///	<summary>
	///		Determine this LatLonBox overlaps with box.
	///	</summary>
	bool Overlaps(const LatLonBox & box) const {

		// Check latitudes
		if (lat[0] > box.lat[1]) {
			return false;
		}
		if (lat[1] < box.lat[0]) {
			return false;
		}

		// Both boxes cross lon 360
		if ((lon[0] > lon[1]) && (box.lon[0] > box.lon[1])) {
			return true;
		}

		// This box crosses lon 360
		if (lon[0] > lon[1]) {
			if (box.lon[1] >= lon[0]) {
				return true;
			}
			if (box.lon[0] <= lon[1]) {
				return true;
			}
			return false;
		}

		// That box crosses lon 360
		if (box.lon[0] > box.lon[1]) {
			if (lon[1] >= box.lon[0]) {
				return true;
			}
			if (lon[0] <= box.lon[1]) {
				return true;
			}
			return false;
		}

		// No boxes cross lon 360
		if (box.lon[1] < lon[0]) {
			return false;
		}
		if (box.lon[0] > lon[1]) {
			return false;
		}
		return true;
	}
};

///////////////////////////////////////////////////////////////////////////////

struct Tag {

	///	<summary>
	///		Identifier associated with this tag.
	///	</summary>
	int id;

	///	<summary>
	///		Time index associated with this tag.
	///	</summary>
	int time;

	///	<summary>
	///		Global id associated with this tag.
	///	</summary>
	int global_id;

	///	<summary>
	///		Default constructor.
	///	</summary>
	Tag() :
		id(0),
		time(0),
		global_id(0)
	{ }

	///	<summary>
	///		Value-based constructor.
	///	</summary>
	Tag(int a_id, int a_time) :
		id(a_id),
		time(a_time),
		global_id(0)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	Tag(const Tag & tag) :
		id(tag.id),
		time(tag.time),
		global_id(tag.global_id)
	{ }

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const Tag & tag) const {
		if (time < tag.time) {
			return true;
		}
		if (time > tag.time) {
			return false;
		}
		if (id < tag.id) {
			return true;
		}
		return false;
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load in the contents of a text file containing one filename per
///		line and store in a vector of strings.
///	</summary>
void GetInputFileList(
	const std::string & strInputFileList,
	std::vector<std::string> & vecInputFiles
) {
	FILE * fp = fopen(strInputFileList.c_str(), "r");

	char szBuffer[1024];
	for (;;) {
		fgets(szBuffer, 1024, fp);

		if (feof(fp)) {
			break;
		}

		// Remove end-of-line characters
		for (;;) {
			int nLen = strlen(szBuffer);
			if ((szBuffer[nLen-1] == '\n') ||
				(szBuffer[nLen-1] == '\r') ||
				(szBuffer[nLen-1] == ' ')
			) {
				szBuffer[nLen-1] = '\0';
				continue;
			}
			break;
		}

		vecInputFiles.push_back(szBuffer);
	}

	if (vecInputFiles.size() == 0) {
		_EXCEPTION1("No files found in file \"%s\"", strInputFileList.c_str());
	}

	fclose(fp);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get a DataVector containing the time variable across a list of
///		input files.
///	</summary>
void GetAllTimes(
	const std::vector<std::string> & vecInputFiles,
	DataVector<double> & dataTimes
) {
	std::vector<double> vecTimes;

	for (int f = 0; f < vecInputFiles.size(); f++) {
		NcFile ncFile(vecInputFiles[f].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		NcDim * dimTime = ncFile.get_dim("time");
		if (dimTime == NULL) {
			_EXCEPTIONT("No dimension \"time\" found in input file");
		}

		NcVar * varTime = ncFile.get_var("time");
		if (varTime == NULL) {
			_EXCEPTIONT("No variable \"time\" found in input file");
		}

		int nTime = dimTime->size();

		DataVector<double> dTime;
		dTime.Initialize(nTime);

		varTime->get(dTime, nTime);

		for (int t = 0; t < nTime; t++) {
			vecTimes.push_back(dTime[t]);
		}
	}

	dataTimes.Initialize(vecTimes.size());
	memcpy(&(dataTimes[0]), &(vecTimes[0]), vecTimes.size() * sizeof(double));
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {

	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Output file
	std::string strOutputFile;

	// Variable name
	std::string strVariable;

	// Output variable name
	std::string strOutputVariable;

	// Minimum patch size
	int nMinPatchSize;

	// Minimum number of timesteps
	int nMinTime;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "inlist", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariable, "var", "");
		CommandLineString(strOutputVariable, "outvar", "");
		CommandLineInt(nMinPatchSize, "minsize", 1);
		CommandLineInt(nMinTime, "mintime", 1);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if ((strInputFile == "") && (strInputFileList == "")) {
		_EXCEPTIONT("No input file (--in) or (--inlist) specified");
	}
	if ((strInputFile != "") && (strInputFileList != "")) {
		_EXCEPTIONT("Only one input file (--in) or (--inlist) allowed");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check variable
	if (strVariable == "") {
		_EXCEPTIONT("No variable name (--var) specified");
	}

	// Check output variable
	if (strOutputVariable.length() == 0) {
		strOutputVariable = strVariable + "tag";
	}

	// Input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile != "") {
		vecInputFiles.push_back(strInputFile);
	}
	if (strInputFileList != "") {
		GetInputFileList(strInputFileList, vecInputFiles);
	}

	int nFiles = vecInputFiles.size();

	// Load in spatial dimension data
	int nLat;
	int nLon;

	DataVector<double> dataLatDeg;
	DataVector<double> dataLat;

	DataVector<double> dataLonDeg;
	DataVector<double> dataLon;

	{
		// Load the first netcdf input file
		NcFile ncInput(vecInputFiles[0].c_str());

		if (!ncInput.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"",
				vecInputFiles[0].c_str());
		}

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

		nLat = dimLat->size();
		nLon = dimLon->size();

		dataLatDeg.Initialize(nLat);
		dataLat.Initialize(nLat);

		dataLonDeg.Initialize(nLon);
		dataLon.Initialize(nLon);

		varLat->get(dataLatDeg, nLat);
		for (int j = 0; j < nLat; j++) {
			dataLat[j] = dataLatDeg[j] * M_PI / 180.0;
		}

		varLon->get(dataLonDeg, nLon);
		for (int i = 0; i < nLon; i++) {
			dataLon[i] = dataLonDeg[i] * M_PI / 180.0;
		}

		// Close first netcdf file
		ncInput.close();
	}

	// Get time dimension over all files
	DataVector<double> dTime;
	GetAllTimes(vecInputFiles, dTime);

	int nTime = dTime.GetRows();

	// Allocate indicator data
	DataMatrix<int> dataIndicator(nLat, nLon);

	// Build blobs at each time level
	AnnounceStartBlock("Building blob set at each time level");

	// Set of indicator locations
	typedef std::set<LatLonPair> IndicatorSet;
	typedef IndicatorSet::iterator IndicatorSetIterator;
	typedef IndicatorSet::const_iterator IndicatorSetConstIterator;

	IndicatorSet setIndicators;
	IndicatorSet setNeighbors;

	// Set of nodes contained in all patches
	std::vector< std::vector<IndicatorSet> > vecAllPatches;
	vecAllPatches.resize(nTime);

	// Bounding boxes for all patches
	std::vector< std::vector<LatLonBox> > vecAllPatchBoxes;
	vecAllPatchBoxes.resize(nTime);

	// Time index across all files
	int iTime = 0;

	// Loop through all files
	for (int f = 0; f < nFiles; f++) {

		// Load in each file
		NcFile ncInput(vecInputFiles[f].c_str());
		if (!ncInput.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		// Get current time dimension
		NcDim * dimTime = ncInput.get_dim("time");

		int nLocalTimes = dimTime->size();

		// Load in indicator variable
		NcVar * varIndicator = ncInput.get_var(strVariable.c_str());

		if (varIndicator->num_dims() != 3) {
			_EXCEPTION1("Incorrect number of dimensions for \"%s\" (3 expected)",
				strVariable.c_str());
		}

		// Loop through all times
		for (int t = 0; t < nLocalTimes; t++, iTime++) {

			// Get the current patch vector
			std::vector<IndicatorSet> & vecPatches = vecAllPatches[iTime];

			std::vector<LatLonBox> & vecPatchBoxes = vecAllPatchBoxes[iTime];

			// New announcement block for timestep
			char szStartBlock[128];
			sprintf(szStartBlock, "Time %i (%i)", iTime, t);
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
					if (vecPatches[ixPatch].find(pr) !=
						vecPatches[ixPatch].end()
					) {
						continue;
					}

					// This node has not been tagged
					IndicatorSetIterator iterIndicator = setIndicators.find(pr);
					if (iterIndicator == setIndicators.end()) {
						continue;
					}

					// Remove this from the set of available indicators
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
							if ((pr.lon >= box.lon[0]) ||
								(pr.lon <= box.lon[1])
							) {
							} else if (
								pr.lon - box.lon[1] < box.lon[0] - pr.lon
							) {
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
				Announce("Blob %i [%i, %i] x [%i, %i]",
					p+1,
					vecPatchBoxes[p].lat[0],
					vecPatchBoxes[p].lat[1],
					vecPatchBoxes[p].lon[0],
					vecPatchBoxes[p].lon[1]);
			}

			AnnounceEndBlock("Done");
		}

		// Close NetCDF file
		ncInput.close();
	}

	AnnounceEndBlock("Done");

	// Stitch blobs together in time
	AnnounceStartBlock("Stitching Blobs");

	// Patch tags for each blob
	std::vector< std::vector<Tag> > vecAllPatchTags;
	vecAllPatchTags.resize(nTime);

	// Next available patch tag
	Tag tagNextPatch(1, 0);

	// Give patch tags to the initial set of patches
	std::set<Tag> setAllTags;

	for (int t = 0; t < nTime; t++) {
		vecAllPatchTags[t].resize(vecAllPatches[t].size());

		tagNextPatch.id = 0;
		for (int p = 0; p < vecAllPatchTags[t].size(); p++) {
			vecAllPatchTags[t][p] = tagNextPatch;
			setAllTags.insert(tagNextPatch);

			tagNextPatch.id++;
		}
		tagNextPatch.time++;
	}

	// Array of equivalent tags
	typedef std::multimap<Tag, Tag> MapGraph;
	typedef MapGraph::const_iterator MapGraphConstIterator;
	typedef MapGraph::iterator MapGraphIterator;

	MapGraph multimapTagGraph;

	// Loop through all remaining time steps
	for (int t = 1; t < nTime; t++) {

		// Get the current patch vector
		const std::vector<Tag> & vecPrevPatchTags = vecAllPatchTags[t-1];

		std::vector<Tag> & vecPatchTags = vecAllPatchTags[t];

		const std::vector<IndicatorSet> & vecPrevPatches = vecAllPatches[t-1];

		const std::vector<IndicatorSet> & vecPatches = vecAllPatches[t];

		const std::vector<LatLonBox> & vecPrevPatchBoxes = vecAllPatchBoxes[t-1];

		const std::vector<LatLonBox> & vecPatchBoxes = vecAllPatchBoxes[t];

		// Determine overlaps between these patches and previous patches
		vecPatchTags.resize(vecPatches.size());
		for (int p = 0; p < vecPatchTags.size(); p++) {

			const LatLonBox & boxP = vecPatchBoxes[p];

			// Find overlap with bounding boxes at previous time
			for (int q = 0; q < vecPrevPatchTags.size(); q++) {

				const LatLonBox & boxQ = vecPrevPatchBoxes[q];

				// Check if bounding boxes overlap
				if (!boxP.Overlaps(boxQ)) {
					continue;
				}
					
				// Verify that at least one node overlaps between patches
				bool fHasOverlapNode = false;
				IndicatorSetConstIterator iter = vecPatches[p].begin();
				for (; iter != vecPatches[p].end(); iter++) {
					if (vecPrevPatches[q].find(*iter) !=
						vecPrevPatches[q].end()
					) {
						fHasOverlapNode = true;
						break;
					}
				}

				if (!fHasOverlapNode) {
					continue;
				}

				// Insert bidirectional edge in graph
				multimapTagGraph.insert(
					std::pair<Tag, Tag>(
						vecPatchTags[p], vecPrevPatchTags[q]));

				multimapTagGraph.insert(
					std::pair<Tag, Tag>(
						vecPrevPatchTags[q], vecPatchTags[p]));
			}
		}
	}

	// Total number of blobs
	int nTotalBlobCount = 0;

	// Find all cliques
	std::map<Tag, Tag> mapEquivalentTags;

	std::set<Tag>::const_iterator iterTag = setAllTags.begin();

	for (; iterTag != setAllTags.end(); iterTag++) {

		std::set<Tag> setTagsVisited;

		std::set<Tag> setTagsToVisit;
		setTagsToVisit.insert(*iterTag);

		Tag tagMinimum = *iterTag;

		// Check if this tag is already part of an explored clique
		if (mapEquivalentTags.find(*iterTag) != mapEquivalentTags.end()) {
			continue;
		}

		// New blob
		nTotalBlobCount++;

		// All time values associated with this blob
		std::set<int> setBlobTimes;

		// Find clique containing this node
		for (;;) {

			// Out of tags to visit; done
			if (setTagsToVisit.size() == 0) {
				break;
			}

			// Get next tag to visit
			Tag tagNext = *(setTagsToVisit.begin());
			setTagsToVisit.erase(setTagsToVisit.begin());

			// Verify we haven't visited this tag already
			if (setTagsVisited.find(tagNext) != setTagsVisited.end()) {
				continue;
			}
			setTagsVisited.insert(tagNext);

			// Check minimum tag
			if (tagNext < tagMinimum) {
				tagMinimum = tagNext;
			}

			// Times
			setBlobTimes.insert(tagNext.time);

			// Get edges from this node
			std::pair<MapGraphIterator, MapGraphIterator> iterGraphEdges
				= multimapTagGraph.equal_range(tagNext);

			MapGraphIterator iter = iterGraphEdges.first;
			for (; iter != iterGraphEdges.second; iter++) {
				setTagsToVisit.insert(iter->second);
			}
		}

		// Set global id
		if (setBlobTimes.size() >= nMinTime) {
			tagMinimum.global_id = nTotalBlobCount;
		} else {
			tagMinimum.global_id = 0;
			nTotalBlobCount--;
		}

		// Refer all tags in clique to minimum tag
		std::set<Tag>::const_iterator iterTagsVisited = setTagsVisited.begin();
		for (; iterTagsVisited != setTagsVisited.end(); iterTagsVisited++) {
			mapEquivalentTags.insert(
				std::pair<Tag,Tag>(*iterTagsVisited, tagMinimum));
		}
	}

	// Apply equivalent tags
	for (int t = 0; t < nTime; t++) {

		std::vector<Tag> & vecPatchTags = vecAllPatchTags[t];

		for (int p = 0; p < vecPatchTags.size(); p++) {
			std::map<Tag, Tag>::const_iterator iterTagPair =
				mapEquivalentTags.find(vecPatchTags[p]);

			if (iterTagPair != mapEquivalentTags.end()) {
				vecPatchTags[p] = iterTagPair->second;
			}
		}
	}

	Announce("Blobs found: %i", nTotalBlobCount);

	AnnounceEndBlock("Done");

	// Output
	AnnounceStartBlock("Output Blobs");

	// Load the netcdf output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\"",
			strOutputFile.c_str());
	}

	// Create output
	NcDim * dimOutputTime = ncOutput.add_dim("time", nTime);
	NcDim * dimOutputLat = ncOutput.add_dim("lat", nLat);
	NcDim * dimOutputLon = ncOutput.add_dim("lon", nLon);

	if (dimOutputTime == NULL) {
		_EXCEPTIONT("Unable to create dimension \"time\" in output");
	}
	if (dimOutputLat == NULL) {
		_EXCEPTIONT("Unable to create dimension \"lat\" in output");
	}
	if (dimOutputLon == NULL) {
		_EXCEPTIONT("Unable to create dimension \"lon\" in output");
	}

	NcVar * varOutputTime = ncOutput.add_var("time", ncDouble, dimOutputTime);
	NcVar * varOutputLat  = ncOutput.add_var("lat", ncDouble, dimOutputLat);
	NcVar * varOutputLon  = ncOutput.add_var("lon", ncDouble, dimOutputLon);

	if (varOutputTime == NULL) {
		_EXCEPTIONT("Unable to create variable \"time\" in output");
	}
	if (varOutputLat == NULL) {
		_EXCEPTIONT("Unable to create variable \"lat\" in output");
	}
	if (varOutputLon == NULL) {
		_EXCEPTIONT("Unable to create variable \"lon\" in output");
	}

	varOutputTime->set_cur((long)0);
	varOutputTime->put(dTime, nTime);

	varOutputLat->set_cur((long)0);
	varOutputLat->put(dataLatDeg, nLat);

	varOutputLon->set_cur((long)0);
	varOutputLon->put(dataLonDeg, nLon);

	// Copy variable attributes from first input file
	{
		NcFile ncInput(vecInputFiles[0].c_str());

		NcVar * varTime = ncInput.get_var("time");
		NcVar * varLat = ncInput.get_var("lat");
		NcVar * varLon = ncInput.get_var("lon");

		CopyNcVarAttributes(varTime, varOutputTime);
		CopyNcVarAttributes(varLat, varOutputLat);
		CopyNcVarAttributes(varLon, varOutputLon);
	}

	NcVar * varData =
		ncOutput.add_var(
			strOutputVariable.c_str(),
			ncInt,
			dimOutputTime,
			dimOutputLat,
			dimOutputLon);

	// Loop through all time steps
	DataMatrix<int> dataPatchTag;
	dataPatchTag.Initialize(nLat, nLon);

	for (int t = 0; t < nTime; t++) {

		dataPatchTag.Zero();

		// Get the current patch vectors
		const std::vector<Tag> & vecPatchTags = vecAllPatchTags[t];

		const std::vector<IndicatorSet> & vecPatches = vecAllPatches[t];

		// Put patch information into matrix
		for (int p = 0; p < vecPatchTags.size(); p++) {

			if (vecPatchTags[p].global_id == 0) {
				continue;
			}

			IndicatorSetConstIterator iter = vecPatches[p].begin();
			for (; iter != vecPatches[p].end(); iter++) {
				dataPatchTag[iter->lat][iter->lon] =
					vecPatchTags[p].global_id;
			}
		}

		// Write to file
		varData->set_cur(t, 0, 0);
		varData->put(&(dataPatchTag[0][0]), 1, nLat, nLon);
	}

	ncOutput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

