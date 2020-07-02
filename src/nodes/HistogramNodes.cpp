///////////////////////////////////////////////////////////////////////////////
///
///	\file    DensityNodes.cpp
///	\author  Paul Ullrich
///	\version January 20, 2015
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

#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <set>

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

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {

	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Input file format
	std::string strInputFormat;

	// NetCDF file containing latitude and longitude arrays
	std::string strLatLonFile;

	// Output file (NetCDF)
	std::string strOutputFile;

	// Output variable name
	std::string strOutputVariable;

	// Column in which the longitude index appears
	int iLonIxCol;

	// Column in which the latitude index appears
	int iLatIxCol;

	// Begin latitude
	double dLatBegin;

	// End latitude
	double dLatEnd;

	// Begin longitude
	double dLonBegin;

	// End longitude
	double dLonEnd;

	// Number of latitudes in output
	int nLat;

	// Number of longitudes in output
	int nLon;

	// Only include genesis points
	bool fGenesis;

	// Online include termination points
	bool fTermination;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "inlist", "");
		CommandLineStringD(strInputFormat, "in_format", "std", "(std|visit)");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputVariable, "outvar", "density");
		CommandLineInt(iLonIxCol, "iloncol", 8);
		CommandLineInt(iLatIxCol, "ilatcol", 9);
		CommandLineDouble(dLatBegin, "lat_begin", -90.0);
		CommandLineDouble(dLatEnd, "lat_end", 90.0);
		CommandLineDouble(dLonBegin, "lon_begin", 0.0);
		CommandLineDouble(dLonEnd, "lon_end", 360.0);
		CommandLineInt(nLat, "nlat", 180);
		CommandLineInt(nLon, "nlon", 360);
		CommandLineBool(fGenesis, "genesis_only");
		CommandLineBool(fTermination, "termination_only");

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
	if (strInputFormat != "std") {
		_EXCEPTIONT("UNIMPLEMENTED:  Only \"--in_format std\" supported");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check output variable
	if (strOutputVariable == "") {
		_EXCEPTIONT("No output variable name (--outvar) specified");
	}

	// Number of latitudes and longitudes
	if (nLat == 0) {
		_EXCEPTIONT("UNIMPLEMENTED: --nlat must be specified currently");
	}
	if (nLon == 0) {
		_EXCEPTIONT("UNIMPLEMENTED: --nlon must be specified currently");
	}

	// Check genesis / termination flags
	if (fGenesis && fTermination) {
		_EXCEPTIONT("Only one of --genesis_only and --termination_only allowed");
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

	// Density
	DataArray2D<int> nCounts(nLat, nLon);

	// Loop through all files in list
	AnnounceStartBlock("Processing files");

	std::string strBuffer;
	strBuffer.reserve(1024);

	for (int f = 0; f < nFiles; f++) {
		Announce("File \"%s\"", vecInputFiles[f].c_str());

		FILE * fp = fopen(vecInputFiles[f].c_str(), "r");
		if (fp == NULL) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		int iLine = 0;
		int nNodesTotal = 0;
		int nNodesRemaining = 0;
		for (;;) {
			iLine++;

			// Read in the next line
			fgets(&(strBuffer[0]), 1024, fp);

			int nLength = strlen(&(strBuffer[0]));

			// Check for end of file
			if (feof(fp)) {
				break;
			}

			// Check for comment line
			if (strBuffer[0] == '#') {
				continue;
			}

			// Check for new storm
			if (strncmp(&(strBuffer[0]), "start", 5) == 0) {

				for (int i = 6; i <= nLength; i++) {
					if ((strBuffer[i] < '0') || (strBuffer[i] > '9')) {
						if (i == 6) {
							_EXCEPTION1("Missing track length on line %i", iLine);
						}
						nNodesRemaining = atoi(&(strBuffer[6]));
						nNodesTotal = nNodesRemaining;
						break;
					}
				}

				continue;
			}

			// Count down number of tracks
			nNodesRemaining--;
			if (nNodesRemaining < 0) {
				_EXCEPTION1("Malformed nodefile on line %i: too many nodes in track", iLine);
			}

			// Parse line
			double dLon;
			double dLat;

			int iCol = 0;
			int iLast = 0;

			bool fWhitespace = true;

			for (int i = 0; i <= nLength; i++) {
				if ((strBuffer[i] == ' ') ||
					(strBuffer[i] == ',') ||
					(strBuffer[i] == '\t') ||
					(strBuffer[i] == '\0')
				) {
					if (!fWhitespace) {
						if (iCol == iLonIxCol) {
							strBuffer[i] = '\0';
							dLon = atof(&(strBuffer[iLast]));
						}
						if (iCol == iLatIxCol) {
							strBuffer[i] = '\0';
							dLat = atof(&(strBuffer[iLast]));
						}
					}

					fWhitespace = true;

				} else {
					if (fWhitespace) {
						iLast = i;
						iCol++;
					}
					fWhitespace = false;
				}
			}

			// Latitude and longitude index
			int iLon =
				static_cast<int>(static_cast<double>(nLon)
					* (dLon - dLonBegin) / (dLonEnd - dLonBegin));
			int iLat =
				static_cast<int>(static_cast<double>(nLat)
					* (dLat - dLatBegin) / (dLatEnd - dLatBegin));

			if (iLon == (-1)) {
				iLon = 0;
			}
			if (iLon == nLon) {
				iLon = nLon - 1;
			}
			if (iLat == (-1)) {
				iLat = 0;
			}
			if (iLat == nLat) {
				iLat = nLat - 1;
			}

			if ((iLat < 0) || (iLat >= nLat)) {
				_EXCEPTION1("Latitude index (%i) out of range", iLat);
			}
			if ((iLon < 0) || (iLon >= nLon)) {
				_EXCEPTION1("Longitude index (%i) out of range", iLon);
			}

			if (fGenesis) {
				if (nNodesRemaining == nNodesTotal - 1) {
					nCounts[iLat][iLon]++;
				}

			} else if (fTermination) {
				if (nNodesRemaining == 0) {
					nCounts[iLat][iLon]++;
				}

			} else {
				nCounts[iLat][iLon]++;
			}
		}

		fclose(fp);
	}

	AnnounceEndBlock("Done");

	// Output results
	AnnounceStartBlock("Output results");

	// Load the netcdf output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\"",
			strOutputFile.c_str());
	}

	// Create output
	NcDim * dimLat = ncOutput.add_dim("lat", nLat);
	NcDim * dimLon = ncOutput.add_dim("lon", nLon);

	NcVar * varLat = ncOutput.add_var("lat", ncDouble, dimLat);
	NcVar * varLon = ncOutput.add_var("lon", ncDouble, dimLon);

	varLat->add_att("units", "degrees_north");
	varLon->add_att("units", "degrees_east");

	DataArray1D<double> dLat(nLat);
	DataArray1D<double> dLon(nLon);

	for (int j = 0; j < nLat; j++) {
		dLat[j] = dLatBegin
			+ (dLatEnd - dLatBegin)
				* (static_cast<double>(j) + 0.5)
				/ static_cast<double>(nLat);
	}
	for (int i = 0; i < nLon; i++) {
		dLon[i] = dLonBegin
			+ (dLonEnd - dLonBegin)
			* (static_cast<double>(i) + 0.5)
			/ static_cast<double>(nLon);
	}

	varLat->put(&(dLat[0]), nLat);
	varLon->put(&(dLon[0]), nLon);

	// Output counts
	NcVar * varCount =
		ncOutput.add_var(
			strOutputVariable.c_str(),
			ncInt,
			dimLat,
			dimLon);

	varCount->put(&(nCounts[0][0]), nLat, nLon);

	ncOutput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

