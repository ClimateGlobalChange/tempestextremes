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

#include "BlobUtilities.h"

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "TimeObj.h"

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
///		Storage structure for data associated with Blobs.
///	</summary>
struct BlobQuantities {

	///	<summary>
	///		Output quantities.
	///	</summary>
	enum OutputQuantity {
		MinLat,
		MaxLat,
		MinLon,
		MaxLon,
		CentroidLat,
		CentroidLon,
		Area
	};

	///	<summary>
	///		Latitude-longitude bounding box for blob.
	///	</summary>
	LatLonBox box;

	///	<summary>
	///		Total area of blob.
	///	</summary>
	double dArea;

	///	<summary>
	///		Array of output variables associated with blob.
	///	</summary>
	std::vector<double> dOutputVars;
};

///	<summary>
///		A map between time and quantities associated with the blob.
///	</summary>
typedef std::map<int, BlobQuantities> TimedBlobQuantitiesMap;

///	<summary>
///		A map between blob index and TimedBlobQuantitiesMap;
///	</summary>
typedef std::map<int, TimedBlobQuantitiesMap> AllTimedBlobQuantitiesMap;

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {

	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Output file
	std::string strOutputFile;

	// Input variable
	std::string strInputVariable;

	// Summary quantities
	std::string strOutputQuantities;

	// Output headers
	bool fOutputHeaders;

	// Output the full times rather than time indexes
	bool fOutputFullTimes;

	// Display help message
	bool fHelp;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "infile", "");
		CommandLineString(strInputFileList, "inlist", "");
		CommandLineString(strOutputFile, "outfile", "");
		CommandLineString(strInputVariable, "invar", "");
		CommandLineString(strOutputQuantities, "out", "");
		CommandLineBool(fOutputHeaders, "out_headers");
		CommandLineBool(fOutputFullTimes, "out_fulltime");
		CommandLineBool(fHelp, "help");

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
		_EXCEPTIONT("No output file (--outfile) specified");
	}

	// Check variable
	if (strInputVariable == "") {
		_EXCEPTIONT("No variable name (--invar) specified");
	}

	// Check output variable
	if (strOutputQuantities == "") {
		_EXCEPTIONT("No output quantities (--out) specified");
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

	DataArray1D<double> dataLatDeg;
	DataArray1D<double> dataLat;

	DataArray1D<double> dataLonDeg;
	DataArray1D<double> dataLon;

	bool fFlippedLat = false;

	double dAreaElement;

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

		dataLatDeg.Allocate(nLat);
		dataLat.Allocate(nLat);

		dataLonDeg.Allocate(nLon);
		dataLon.Allocate(nLon);

		varLat->get(dataLatDeg, nLat);
		for (int j = 0; j < nLat; j++) {
			dataLat[j] = dataLatDeg[j] * M_PI / 180.0;
		}

		varLon->get(dataLonDeg, nLon);
		for (int i = 0; i < nLon; i++) {
			dataLon[i] = dataLonDeg[i] * M_PI / 180.0;
		}

		dAreaElement =
			M_PI / static_cast<double>(nLat)
			* 2.0 * M_PI / static_cast<double>(nLon);

		// Check for flipped latitude
		if (dataLatDeg.GetRows() >= 2) {
			if (dataLatDeg[1] < dataLatDeg[0]) {
				fFlippedLat = true;
			}
		}

		// Close first netcdf file
		ncInput.close();
	}

	// Parse the list of output quantities
	std::vector<BlobQuantities::OutputQuantity> vecOutputVars;

	{
		// Parse output quantities and store in vecOutputVars
		AnnounceStartBlock("Parsing output quantities");

		int iLast = 0;
		for (int i = 0; i <= strOutputQuantities.length(); i++) {

			if ((i == strOutputQuantities.length()) ||
				(strOutputQuantities[i] == ',') ||
				(strOutputQuantities[i] == ';')
			) {
				std::string strSubStr =
					strOutputQuantities.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecOutputVars.size());
				vecOutputVars.resize(iNextOp + 1);

				if (strSubStr == "minlat") {
					vecOutputVars[iNextOp] = BlobQuantities::MinLat;
				} else if (strSubStr == "maxlat") {
					vecOutputVars[iNextOp] = BlobQuantities::MaxLat;
				} else if (strSubStr == "minlon") {
					vecOutputVars[iNextOp] = BlobQuantities::MinLon;
				} else if (strSubStr == "maxlon") {
					vecOutputVars[iNextOp] = BlobQuantities::MaxLon;
				} else if (strSubStr == "centlat") {
					vecOutputVars[iNextOp] = BlobQuantities::CentroidLat;
				} else if (strSubStr == "centlon") {
					vecOutputVars[iNextOp] = BlobQuantities::CentroidLon;
				} else if (strSubStr == "area") {
					vecOutputVars[iNextOp] = BlobQuantities::Area;
				} else {
					_EXCEPTIONT("Invalid output quantity:  Expected\n"
						"[minlat, maxlat, minlon, maxlon, "
						"centlat, centlon, area]");
				}

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Computed quantities associated with each blob
	AllTimedBlobQuantitiesMap mapAllQuantities;

	// Time index across all files
	int iTime = 0;

	// Open output file
	FILE * fpout = fopen(strOutputFile.c_str(), "w");
	if (fpout == NULL) {
		_EXCEPTIONT("Error opening output file");
	}

	// Output header
	if (fOutputHeaders) {
		if (fOutputFullTimes) {
			fprintf(fpout, "time,%s\n", strOutputQuantities.c_str());
		} else {
			fprintf(fpout, "time_ix,%s\n", strOutputQuantities.c_str());
		}
	}

	// Loop through all files
	for (int f = 0; f < nFiles; f++) {

		// Time objects for each time in this file
		std::vector<Time> vecFileTimes;

		// Load in each file
		NcFile ncInput(vecInputFiles[f].c_str());
		if (!ncInput.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		// Blob index data
		DataArray2D<int> dataIndex(nLat, nLon);

		// Get current time dimension
		NcDim * dimTime = ncInput.get_dim("time");
		if (dimTime == NULL) {
			_EXCEPTIONT("Dimension \"time\" missing from input file");
		}

		int nLocalTimes = dimTime->size();

		// Load in file times
		if (fOutputFullTimes) {
			ReadCFTimeDataFromNcFile(&ncInput, vecInputFiles[f], vecFileTimes, true);
		}

		// Load in indicator variable
		NcVar * varIndicator = ncInput.get_var(strInputVariable.c_str());

		if (varIndicator == NULL) {
			_EXCEPTION1("Unable to load variable \"%s\"",
				strInputVariable.c_str());
		}

		if (varIndicator->num_dims() != 3) {
			_EXCEPTION1("Incorrect number of dimensions for \"%s\""
				" (3 expected)", strInputVariable.c_str());
		}

		// Loop through all times
		for (int t = 0; t < nLocalTimes; t++, iTime++) {

			// Load in the data at this time
			varIndicator->set_cur(t, 0, 0);
			varIndicator->get(&(dataIndex[0][0]), 1, nLat, nLon);

			// Last data index
			int iLastDataIndex = 0;

			// Iterator to BlobQuantities with iLastDataIndex and iTime
			TimedBlobQuantitiesMap::iterator iterBlobQuantities;

			// Loop over all locations
			for (int j = 0; j < nLat; j++) {
			for (int i = 0; i < nLon; i++) {

				// Ignore non-blob data
				if (dataIndex[j][i] == 0) {
					continue;
				}

				// Check if iterator already points to correct data
				if (dataIndex[j][i] != iLastDataIndex) {

					// Iterator to TimedBlobQuantities with iLastDataIndex
					AllTimedBlobQuantitiesMap::iterator
						iterTimedBlobQuantities;

					// Get new iterator for this dataIndex
					iterTimedBlobQuantities =
						mapAllQuantities.find(
							dataIndex[j][i]);

					if (iterTimedBlobQuantities == mapAllQuantities.end()) {

						std::pair<AllTimedBlobQuantitiesMap::iterator,bool> pr =
							mapAllQuantities.insert(
								AllTimedBlobQuantitiesMap::value_type(
									dataIndex[j][i],
									TimedBlobQuantitiesMap()));

						if (pr.second == false) {
							_EXCEPTIONT("Map insertion failed");
						}
						iterTimedBlobQuantities = pr.first;
					}

					// Get the BlobQuantities at the current time
					TimedBlobQuantitiesMap & mapTimedBlobQuantities =
						iterTimedBlobQuantities->second;

					iterBlobQuantities =
						mapTimedBlobQuantities.find(iTime);

					if (iterBlobQuantities == mapTimedBlobQuantities.end()) {
						std::pair<TimedBlobQuantitiesMap::iterator,bool> pr =
							mapTimedBlobQuantities.insert(
								TimedBlobQuantitiesMap::value_type(
									iTime,
									BlobQuantities()));

						if (pr.second == false) {
							_EXCEPTIONT("Map insertion failed");
						}
						iterBlobQuantities = pr.first;
					}

					// Update last data index
					iLastDataIndex = dataIndex[j][i];
				}

				//std::cout << iterBlobQuantities->first << std::endl;
/*
				if ((dataIndex[j][i] == 1) && (iTime == 1)) {
					printf("%i %i : %i %i\n",
						j, i,
						iterBlobQuantities->second.box.lon[0],
						iterBlobQuantities->second.box.lon[1]);
				}
*/
				// Insert point into array
				iterBlobQuantities->second.box.InsertPoint(j, i, nLat, nLon);

				// Add blob area
				iterBlobQuantities->second.dArea +=
					cos(dataLat[j]) * dAreaElement;

				if ((t == 3) && (dataIndex[j][i] == 2)) {
					printf ("%1.5e %1.5e %1.5e %1.5e %1.5e\n",
						dataLat[j] * 180.0 / M_PI, dataLon[i] * 180.0 / M_PI, cos(dataLat[j]), dAreaElement, iterBlobQuantities->second.dArea);
				}

			}
			}
		}

		// Output all BlobQuantities
		{
			AllTimedBlobQuantitiesMap::iterator iterBlobs =
				mapAllQuantities.begin();

			for (; iterBlobs != mapAllQuantities.end(); iterBlobs++) {

				fprintf(fpout, "Blob %i (%lu)\n",
					iterBlobs->first,
					iterBlobs->second.size());

				TimedBlobQuantitiesMap::iterator iterTimes =
					iterBlobs->second.begin();

				for (; iterTimes != iterBlobs->second.end(); iterTimes++) {

					if (fOutputFullTimes) {
						fprintf(fpout, "%s", vecFileTimes[iterTimes->first].ToShortString().c_str());
					} else {
						fprintf(fpout, "%i", iterTimes->first);
					}

					BlobQuantities & quants = iterTimes->second;
					for (int i = 0; i < vecOutputVars.size(); i++) {

						// Bounding box coordinates
						double dLat0 = dataLatDeg[quants.box.lat[0]];
						double dLat1 = dataLatDeg[quants.box.lat[1]];
						double dLon0 = dataLonDeg[quants.box.lon[0]];
						double dLon1 = dataLonDeg[quants.box.lon[1]];

						// Minimum latitude
						if (vecOutputVars[i] == BlobQuantities::MinLat) {
							if (fFlippedLat) {
								fprintf(fpout, "\t%1.5f", dLat1);
							} else {
								fprintf(fpout, "\t%1.5f", dLat0);
							}
						}

						// Maximum latitude
						if (vecOutputVars[i] == BlobQuantities::MaxLat) {
							if (fFlippedLat) {
								fprintf(fpout, "\t%1.5f", dLat0);
							} else {
								fprintf(fpout, "\t%1.5f", dLat1);
							}
						}

						// Minimum longitude
						if (vecOutputVars[i] == BlobQuantities::MinLon) {
							fprintf(fpout, "\t%1.5f", dLon0);
						}

						// Maximum longitude
						if (vecOutputVars[i] == BlobQuantities::MaxLon) {
							fprintf(fpout, "\t%1.5f", dLon1);
						}

						// Centroid latitude
						if (vecOutputVars[i] == BlobQuantities::CentroidLat) {
							double dMidLat = 0.5 * (dLat0 + dLat1);

							fprintf(fpout, "\t%1.5f", dMidLat);
						}

						// Centroid longitude
						if (vecOutputVars[i] == BlobQuantities::CentroidLon) {
							double dMidLon;
							if (dLon0 <= dLon1) {
								dMidLon = 0.5 * (dLon0 + dLon1);
							} else {
								dMidLon = 0.5 * (dLon0 + dLon1 + 360.0);
								if (dMidLon > 360.0) {
									dMidLon -= 360.0;
								}
							}

							fprintf(fpout, "\t%1.5f", dMidLon);
						}

						// Area
						if (vecOutputVars[i] == BlobQuantities::Area) {
							fprintf(fpout, "\t%1.5f", quants.dArea);
						}
					}

					// Endline
					fprintf(fpout, "\n");
				}
			}
		}
	}

	fclose(fpout);

	// Done
	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

