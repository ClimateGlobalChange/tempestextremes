///////////////////////////////////////////////////////////////////////////////
///
///	\file    JiangARs.cpp
///	\author  Paul Ullrich
///	\version September 11th, 2016
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

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input file
	std::string strInputFile;

	// Output file
	std::string strOutputFile;

	// Input integrated water vapor variable name
	std::string strIWVVariable;

	// Output variable name
	std::string strOutputVariable;

	// Minimum absolute latitude
	double dMinAbsLat;

	// Minimum pointwise integrated water vapor
	double dMinIWV;

	// Minimum area
	double dMinArea;

	// Maximal areal fraction
	double dMaxArealFraction;

	// Zonal mean weight
	double dZonalMeanWeight;

	// Zonal max weight
	double dZonalMaxWeight;

	// Meridional mean weight
	double dMeridMeanWeight;

	// Meridional max weight
	double dMeridMaxWeight;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		//CommandLineString(strInputFileList, "inlist", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strIWVVariable, "iwvvar", "");
		CommandLineString(strOutputVariable, "outvar", "");
		CommandLineDouble(dMinAbsLat, "minabslat", 15.0);
		CommandLineDouble(dMinIWV, "miniwv", 20.0);
		CommandLineDouble(dZonalMeanWeight, "zonalmeanwt", 0.7);
		CommandLineDouble(dZonalMaxWeight, "zonalmaxwt", 0.3);
		CommandLineDouble(dMeridMeanWeight, "meridmeanwt", 0.9);
		CommandLineDouble(dMeridMaxWeight, "meridmaxwt", 0.1);
		//CommandLineBool(fRegional, "regional");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	AnnounceStartBlock("Loading data");

	// Check input
	if (strInputFile == "") {
		_EXCEPTIONT("No input file (--in) specified");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check variable
	if (strIWVVariable == "") {
		_EXCEPTIONT("No IWV variable name (--iwvvar) specified");
	}

	// Check output variable
	if (strOutputVariable.length() == 0) {
		strOutputVariable = strIWVVariable + "tag";
	}

	// Open the NetCDF input file
	NcFile ncInput(strInputFile.c_str());

	if (!ncInput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for reading",
			strInputFile.c_str());
	}

#pragma message "Change to file create mode"
	// Open the NetCDF output file
	NcFile ncOutput(strOutputFile.c_str());

	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Get the longitude dimension
	NcDim * dimLon = ncInput.get_dim("lon");
	if (dimLon == NULL) {
		_EXCEPTIONT("Error accessing dimension \"lon\"");
	}

	// Get the longitude variable
	NcVar * varLon = ncInput.get_var("lon");
	if (varLon == NULL) {
		_EXCEPTIONT("Error accessing variable \"lon\"");
	}

	DataVector<double> dLonDeg(dimLon->size());
	varLon->get(&(dLonDeg[0]), dimLon->size());

	// Get the latitude dimension
	NcDim * dimLat = ncInput.get_dim("lat");
	if (dimLat == NULL) {
		_EXCEPTIONT("Error accessing dimension \"lat\"");
	}

	// Get the latitude variable
	NcVar * varLat = ncInput.get_var("lat");
	if (varLat == NULL) {
		_EXCEPTIONT("Error accessing variable \"lat\"");
	}

	DataVector<double> dLatDeg(dimLat->size());
	varLat->get(&(dLatDeg[0]), dimLat->size());

	// Get the integrated water vapor variable
	NcVar * varIWV = ncInput.get_var(strIWVVariable.c_str());
	if (varIWV == NULL) {
		_EXCEPTION1("Error accessing variable \"%s\"",
			strIWVVariable.c_str());
	}

	DataMatrix<float> dIWV(dimLon->size(), dimLat->size());
	varIWV->get(&(dIWV[0][0]), dimLon->size(), dimLat->size());

	// Announce
	AnnounceEndBlock("Done");

	AnnounceStartBlock("Compute zonal/meridional thresholds");

	// Compute zonal threshold
	DataVector<float> dZonalThreshold(dimLat->size());
	for (int j = 0; j < dimLat->size(); j++) {
		float dMaxZonalIWV = dIWV[0][j];
		for (int i = 0; i < dimLon->size(); i++) {
			dZonalThreshold[j] += dIWV[i][j];
			if (dIWV[i][j] > dMaxZonalIWV) {
				dMaxZonalIWV = dIWV[i][j];
			}
		}
		dZonalThreshold[j] /= static_cast<float>(dimLon->size());

		dZonalThreshold[j] =
			dZonalMeanWeight * dZonalThreshold[j]
			+ dZonalMaxWeight * dMaxZonalIWV;
	}

	// Compute meridional threshold
	DataVector<float> dMeridThreshold(dimLon->size());
	for (int i = 0; i < dimLon->size(); i++) {
		float dMaxMeridIWV = dIWV[i][0];
		for (int j = 0; j < dimLat->size(); j++) {
			dMeridThreshold[i] += dIWV[i][j];
			if (dIWV[i][j] > dMaxMeridIWV) {
				dMaxMeridIWV = dIWV[i][j];
			}
		}
		dMeridThreshold[i] /= static_cast<float>(dimLon->size());

		dMeridThreshold[i] =
			dMeridMeanWeight * dMeridThreshold[i]
			+ dMeridMaxWeight * dMaxMeridIWV;
	}

	// Announce
	AnnounceEndBlock("Done");

	AnnounceStartBlock("Build tagged cell array");

	// Build tagged cell array
	DataMatrix<int> dIWVtag(dimLon->size(), dimLat->size());

	for (int i = 0; i < dimLon->size(); i++) {
	for (int j = 0; j < dimLat->size(); j++) {
		if (fabs(dLatDeg[j]) < dMinAbsLat) {
			continue;
		}
		if (dIWV[i][j] < dMinIWV) {
			continue;
		}
		if (dIWV[i][j] < dZonalThreshold[j]) {
			continue;
		}
		if (dIWV[i][j] < dMeridThreshold[i]) {
			continue;
		}

		dIWVtag[i][j] = 1;
	}
	}

	// Announce
	AnnounceEndBlock("Done");

	AnnounceStartBlock("Writing results");

	// Output tagged cell array
	CopyNcVar(ncInput, ncOutput, "lat", true);
	CopyNcVar(ncInput, ncOutput, "lon", true);

	NcVar * varIWVtag =
		ncOutput.add_var(strOutputVariable.c_str(), ncInt, dimLon, dimLat);

	varIWVtag->put(&(dIWVtag[0][0]), dimLon->size(), dimLat->size());

	AnnounceEndBlock("Done");

/*
	1) Restrict domain to >15N/S
2) Grid point IWV > 20mm*
3) IWV ≥ IWVzonal mean  +  0.3**(IWVzonal max - IWVzonal mean)  		and
IWV ≥ IWVmeridional mean  +  0.1**(IWVmeridional max - IWVmeridional mean)
4) A connected component labeling algorithm to discard small structures (< 2.5x10^5 km2)***
5) Shape analysis: confine AR region to rectagular box "of minimal size" (?) and calculate areal fraction AR occupies in box. If the fraction <0.75, the AR is "elongated" and it passes.
6) Check orientation by regressing lats onto lons (make sure it is not northeast)****
*/

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

