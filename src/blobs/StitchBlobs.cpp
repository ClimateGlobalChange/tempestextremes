///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchBlobs.cpp
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

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

#include "Constants.h"
#include "CoordTransforms.h"
#include "BlobUtilities.h"

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"

#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"
#include "SimpleGrid.h"
#include "Variable.h"

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
	///		Global id associated with this tag (minimum value 1).
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

// Set of indicator locations stored as latitude-longitude pairs
typedef std::set<int> IndicatorSet;
typedef IndicatorSet::iterator IndicatorSetIterator;
typedef IndicatorSet::const_iterator IndicatorSetConstIterator;

///////////////////////////////////////////////////////////////////////////////

class BlobThresholdOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum ThresholdQuantity {
		MinArea,
		MaxArea,
		MinArealFraction,
		MaxArealFraction,
		//EastwardOrientation,
		//WestwardOrientation
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Quantity,
			ReadMode_Value,
			//ReadMode_MinCount,
			ReadMode_Invalid
		} eReadMode = ReadMode_Quantity;

		// Loop through string
		int iLast = 0;
		for (int i = 0; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in threshold quantity
				if (eReadMode == ReadMode_Quantity) {
					if (strSubStr == "minarea") {
						m_eQuantity = MinArea;
						eReadMode = ReadMode_Value;

					} else if (strSubStr == "maxarea") {
						m_eQuantity = MaxArea;
						eReadMode = ReadMode_Value;

					} else if (strSubStr == "minarealfraction") {
						m_eQuantity = MinArealFraction;
						eReadMode = ReadMode_Value;

					} else if (strSubStr == "maxarealfraction") {
						m_eQuantity = MaxArealFraction;
						eReadMode = ReadMode_Value;
/*
					} else if (strSubStr == "eastwardorientation") {
						m_eQuantity = EastwardOrientation;
						eReadMode = ReadMode_Invalid;

					} else if (strSubStr == "westwardorientation") {
						m_eQuantity = WestwardOrientation;
						eReadMode = ReadMode_Invalid;
*/
					} else {
						_EXCEPTION1("Threshold invalid quantity \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					m_dValue = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;
/*
				// Read in minimum count
				} else if (eReadMode == ReadMode_MinCount) {
					if (strSubStr == "all") {
						m_nMinimumCount = (-1);
					} else {
						m_nMinimumCount = atoi(strSubStr.c_str());
					}

					if (m_nMinimumCount < -1) {
						_EXCEPTION1("Invalid minimum count \"%i\"",
							m_nMinimumCount);
					}

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;
*/
				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("Too many entries in threshold string \"%s\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("Insufficient entries in threshold string \"%s\"",
					strOp.c_str());
		}

		// Output announcement
		std::string strDescription;

		char szValue[128];
		sprintf(szValue, "%f", m_dValue);

		if (m_eQuantity == MinArea) {
			strDescription += "Minimum area ";
			strDescription += szValue;
		} else if (m_eQuantity == MaxArea) {
			strDescription += "Maximum area ";
			strDescription += szValue;
		} else if (m_eQuantity == MinArealFraction) {
			strDescription += "Minimum areal fraction ";
			strDescription += szValue;
		} else if (m_eQuantity == MaxArealFraction) {
			strDescription += "Maximum areal fraction ";
			strDescription += szValue;
		//} else if (m_eQuantity == EastwardOrientation) {
		//	strDescription += "Eastward orientation ";
		//} else if (m_eQuantity == WestwardOrientation) {
		//	strDescription += "Westward orientation ";
		}

		Announce("%s", strDescription.c_str());
	}

	///	<summary>
	///		Verify that the specified path satisfies the threshold op.
	///	</summary>
	bool Apply(
		const SimpleGrid & grid,
		const IndicatorSet & setBlobPoints,
		const LatLonBox<double> & boxBlobDeg
	) {
		// Thresholds related to area
		if ((m_eQuantity == MinArea) ||
		    (m_eQuantity == MaxArea) ||
		    (m_eQuantity == MinArealFraction) ||
		    (m_eQuantity == MaxArealFraction)
		) {
			// Calculate the area of each blob box
			double dBoxArea = DegToRad(boxBlobDeg.lon[1]) - DegToRad(boxBlobDeg.lon[0]);
			if (boxBlobDeg.lon[0] > boxBlobDeg.lon[1]) {
				dBoxArea += 2.0 * M_PI;
			}
			dBoxArea *=
				fabs(sin(DegToRad(boxBlobDeg.lat[1])) - sin(DegToRad(boxBlobDeg.lat[0])));

			// Calculate the area and mean lat/lon of each blob
			double dBlobArea = 0.0;
			IndicatorSetIterator iterBlob = setBlobPoints.begin();
			for (; iterBlob != setBlobPoints.end(); iterBlob++) {
				dBlobArea += grid.m_dArea[*iterBlob];
			}

			// Minimum area
			if (m_eQuantity == MinArea) {
				if (dBlobArea < m_dValue) {
					return false;
				}

			// Maximum area
			} else if (m_eQuantity == MaxArea) {
				if (dBlobArea > m_dValue) {
					return false;
				}

			// Minimum areal fraction
			} else if (m_eQuantity == MinArealFraction) {
				if (dBlobArea < m_dValue * dBoxArea) {
					return false;
				}

			// Maximum areal fraction
			} else if (m_eQuantity == MaxArealFraction) {
				if (dBlobArea > m_dValue * dBoxArea) {
					return false;
				}
			}
/*
		// Thresholds related to orientation
		} else if ((m_eQuantity == EastwardOrientation) ||
		           (m_eQuantity == WestwardOrientation)
		) {
			double dNorthHemiMeanLat = 0.0;
			double dNorthHemiMeanLon = 0.0;
			double dNorthHemiMeanLon2 = 0.0;
			double dNorthHemiCoLatLon = 0.0;

			double dSouthHemiMeanLat = 0.0;
			double dSouthHemiMeanLon = 0.0;
			double dSouthHemiMeanLon2 = 0.0;
			double dSouthHemiCoLatLon = 0.0;

			// Calculate regression coefficients for this blob
			IndicatorSetIterator iterBlob = setBlobPoints.begin();
			for (; iterBlob != setBlobPoints.end(); iterBlob++) {

				double dAltLon = 0.0;
				if (dLatDeg[iterBlob->lat] > 0.0) {
					if (iterBlob->lon < boxBlob.lon[0]) {
						dAltLon = dLonDeg[iterBlob->lon] + 360.0;
					} else {
						dAltLon = dLonDeg[iterBlob->lon];
					}

					dNorthHemiMeanLat += dLatDeg[iterBlob->lat];
					dNorthHemiMeanLon += dAltLon;
					dNorthHemiMeanLon2 += dAltLon * dAltLon;
					dNorthHemiCoLatLon += dLatDeg[iterBlob->lat] * dAltLon;

				} else if (dLatDeg[iterBlob->lat] < 0.0) {
					if (iterBlob->lon < boxBlob.lon[0]) {
						dAltLon = dLonDeg[iterBlob->lon] + 360.0;
					} else {
						dAltLon = dLonDeg[iterBlob->lon];
					}

					dSouthHemiMeanLat += dLatDeg[iterBlob->lat];
					dSouthHemiMeanLon += dAltLon;
					dSouthHemiMeanLon2 += dAltLon * dAltLon;
					dSouthHemiCoLatLon += dLatDeg[iterBlob->lat] * dAltLon;
				}
			}

			double dBlobCount = static_cast<double>(setBlobPoints.size());

			dNorthHemiMeanLat /= dBlobCount;
			dNorthHemiMeanLon /= dBlobCount;
			dNorthHemiMeanLon2 /= dBlobCount;
			dNorthHemiCoLatLon /= dBlobCount;

			dSouthHemiMeanLat /= dBlobCount;
			dSouthHemiMeanLon /= dBlobCount;
			dSouthHemiMeanLon2 /= dBlobCount;
			dSouthHemiCoLatLon /= dBlobCount;

			// Calculate the slope of the regression line
			double dNorthSlopeNum =
				dNorthHemiCoLatLon
					- dNorthHemiMeanLat * dNorthHemiMeanLon;

			double dNorthSlopeDen =
				dNorthHemiMeanLon2
					- dNorthHemiMeanLon * dNorthHemiMeanLon;

			double dSouthSlopeNum =
				dSouthHemiCoLatLon
					- dSouthHemiMeanLat * dSouthHemiMeanLon;

			double dSouthSlopeDen =
				dSouthHemiMeanLon2
					- dSouthHemiMeanLon * dSouthHemiMeanLon;

			// Check orientation
			if (m_eQuantity == EastwardOrientation) {

				if (dNorthSlopeNum * dNorthSlopeDen < 0.0) {
					return false;
				}
				if (dSouthSlopeNum * dSouthSlopeDen > 0.0) {
					return false;
				}

			} else if (m_eQuantity == WestwardOrientation) {
				if (dNorthSlopeNum * dNorthSlopeDen > 0.0) {
					return false;
				}
				if (dSouthSlopeNum * dSouthSlopeDen < 0.0) {
					return false;
				}
			}
*/
		// Invalid quantity
		} else {
			_EXCEPTIONT("Invalid quantity");
		}

		return true;
	}

protected:
	///	<summary>
	///		Threshold quantity.
	///	</summary>
	ThresholdQuantity m_eQuantity;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

class RestrictRegion {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RestrictRegion() :
		m_fActive(false),
		m_dBox(),
		m_nCount(0)
	{ }

	///	<summary>
	///		Check if this RestrictRegion is active.
	///	</summary>
	bool IsActive() {
		return m_fActive;
	}

	///	<summary>
	///		Get the minimum count associated with this threshold.
	///	</summary>
	int GetMinimumCount() const {
		return m_nCount;
	}

	///	<summary>
	///		Parse a restrict region argument string, consisting of a list of
	///		5 elements.
	///	</summary>
	void Parse(
		const std::string & strArg
	) {
		int iMode = 0;

		int iLast = 0;
		for (int i = 0; i <= strArg.length(); i++) {

			if ((i == strArg.length()) || (strArg[i] == ',')) {

				std::string strSubStr =
					strArg.substr(iLast, i - iLast);

				if (iMode == 0) {
					m_dBox.lat[0] = atof(strSubStr.c_str());

				} else if (iMode == 1) {
					m_dBox.lat[1] = atof(strSubStr.c_str());

				} else if (iMode == 2) {
					m_dBox.lon[0] = atof(strSubStr.c_str());

				} else if (iMode == 3) {
					m_dBox.lon[1] = atof(strSubStr.c_str());

				} else if (iMode == 4) {
					m_nCount = atoi(strSubStr.c_str());
				}

				iMode++;
				iLast = i + 1;
			}
		}

		if (iMode != 5) {
			_EXCEPTION1("Incorrect number of arguments in --restrict_region: "
				"5 arguments expected, %i found", iMode);
		}

		if (m_nCount < 1) {
			_EXCEPTION1("Invalid value of count in --restrict_region (%i)", m_nCount);
		}

		m_fActive = true;
		m_dBox.is_null = false;

		Announce("Blob must enter region [%1.5f %1.5f] x [%1.5f %1.5f] for at least %i times",
			m_dBox.lat[0], m_dBox.lat[1], m_dBox.lon[0], m_dBox.lon[1], m_nCount);
	}

	///	<summary>
	///		Apply the operator.
	///	</summary>
	bool ContainsPoint(
		double dLatDeg,
		double dLonDeg
	) {
		return m_dBox.contains(dLatDeg, dLonDeg);
	}

protected:
	///	<summary>
	///		A flag indicating this operator is active
	///	</summary>
	bool m_fActive;

	///	<summary>
	///		A latitude-longitude box indicating the region of interest.
	///	</summary>
	LatLonBox<double> m_dBox;

	///	<summary>
	///		The count of the number of time points the blob must be in
	///		this region.
	///	</summary>
	int m_nCount;
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	NcError error(NcError::silent_nonfatal);

try {

#if defined(TEMPEST_MPIOMP)
	// Throw an error if this executable is being run in parallel
	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

	if (nMPISize != 1) {
		_EXCEPTIONT("StitchBlobs does not support parallel MPI operation: "
			"Rerun with one thread.");
	}
#endif

	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Output file
	std::string strOutputFile;

	// Output file list
	std::string strOutputFileList;

	// Variable name
	std::string strVariable;

	// Output variable name
	std::string strOutputVariable;

	// Minimum blob size (in grid points)
	int nMinBlobSize;

	// Minimum duration of blob
	std::string strMinTime;

	// Minimum percentage of the earlier blob that needs to be covered by the later blob
	double dMinPercentOverlapPrev;

	// Maximum percentage of the earlier blob that needs to be covered by the later blob
	double dMaxPercentOverlapPrev;

	// Minimum percentage of the later blob that needs to be covered by the earlier blob
	double dMinPercentOverlapNext;

	// Maximum percentage of the later blob that needs to be covered by the earlier blob
	double dMaxPercentOverlapNext;

	// Restrict blobs to those passing through a given region
	std::string strRestrictRegion;

	// Minimum latitude for detections
	double dMinLatDeg;

	// Maximum latitude for detections
	double dMaxLatDeg;

	// Minimum longitude for detections
	double dMinLonDeg;

	// Maximum longitude for detections
	double dMaxLonDeg;

	// Operate on a regional area (no periodic boundaries)
	bool fRegional;

	// Threshold commands
	std::string strThresholdCmd;

	// Latitude variable dimension name
	std::string strLatitudeName;

	// Longitude variable dimension name
	std::string strLongitudeName;

	// Time variable dimension name
	//std::string strTimeName;

	// Time variable units
	std::string strOutTimeUnits;

	// Verbose output
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "in_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputFileList, "out_list", "");
		CommandLineString(strVariable, "var", "");
		CommandLineString(strOutputVariable, "outvar", "");
		CommandLineInt(nMinBlobSize, "minsize", 1);
		CommandLineString(strMinTime, "mintime", "1");
		CommandLineDoubleD(dMinPercentOverlapPrev, "min_overlap_prev", 0.0, "(%)")
		CommandLineDoubleD(dMaxPercentOverlapPrev, "max_overlap_prev", 100.0, "(%)")
		CommandLineDoubleD(dMinPercentOverlapNext, "min_overlap_next", 0.0, "(%)")
		CommandLineDoubleD(dMaxPercentOverlapNext, "max_overlap_next", 100.0, "(%)")
		CommandLineStringD(strRestrictRegion, "restrict_region", "", "(lat0,lat1,lon0,lon1,count)");
		CommandLineBool(fRegional, "regional");
		CommandLineDouble(dMinLatDeg, "minlat", -90.0);
		CommandLineDouble(dMaxLatDeg, "maxlat", 90.0);
		CommandLineDouble(dMinLonDeg, "minlon", 0.0);
		CommandLineDouble(dMaxLonDeg, "maxlon", 360.0);
		CommandLineString(strLatitudeName, "latname","lat");
		CommandLineString(strLongitudeName, "lonname","lon");
		//CommandLineString(strTimeName, "timename", "time");
		CommandLineString(strOutTimeUnits,"outtimeunits","");
		CommandLineString(strThresholdCmd, "thresholdcmd", "");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Check input
	if ((strInputFile == "") && (strInputFileList == "")) {
		_EXCEPTIONT("No input file (--in) or (--in_list) specified");
	}
	if ((strInputFile != "") && (strInputFileList != "")) {
		_EXCEPTIONT("Only one of input file (--in) or (--in_list) allowed");
	}

	// Check output
	if ((strOutputFile == "") && (strOutputFileList == "")) {
		_EXCEPTIONT("No output file (--out) or (--out_list) specified");
	}
	if ((strOutputFile != "") && (strOutputFileList != "")) {
		_EXCEPTIONT("Only one of input file (--in) or (--in_list) allowed");
	}

	// Check variable
	if (strVariable == "") {
		_EXCEPTIONT("No variable name (--var) specified");
	}

	// Register variable
	int varix = varreg.FindOrRegister(strVariable);

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

	// Output file list
	std::vector<std::string> vecOutputFiles;

	if (strOutputFile != "") {
		vecOutputFiles.push_back(strOutputFile);
	}
	if (strOutputFileList != "") {
		GetInputFileList(strOutputFileList, vecOutputFiles);

		if (vecOutputFiles.size() != vecInputFiles.size()) {
			_EXCEPTION2("Mismatch in number of rows of --in_list (%lu) and --out_list (%lu)",
				vecInputFiles.size(), vecOutputFiles.size());
		}
	}

	// Parse --mintime
	int nMinTime = 1;
	double dMinTimeSeconds = 0.0;
	if (STLStringHelper::IsIntegerIndex(strMinTime)) {
		nMinTime = std::stoi(strMinTime);
		if (nMinTime < 1) {
			_EXCEPTIONT("Invalid value of --mintime; expected positive integer or time delta");
		}

	} else {
		Time timeMinTime;
		timeMinTime.FromFormattedString(strMinTime);
		dMinTimeSeconds = timeMinTime.AsSeconds();
		if (dMinTimeSeconds <= 0.0) {
			_EXCEPTIONT("Invalid value for --mintime; expected positive integer or positive time delta");
		}
	}

	// Convert percent overlap into a decimal value
	if ((dMinPercentOverlapPrev < 0.0) || (dMinPercentOverlapPrev > 100.0)) {
		_EXCEPTIONT("--min_overlap_prev must take on values between 0 and 100");
	}
	if ((dMaxPercentOverlapPrev < 0.0) || (dMaxPercentOverlapPrev > 100.0)) {
		_EXCEPTIONT("--max_overlap_prev must take on values between 0 and 100");
	}
	if ((dMinPercentOverlapNext < 0.0) || (dMinPercentOverlapNext > 100.0)) {
		_EXCEPTIONT("--min_overlap_next must take on values between 0 and 100");
	}
	if ((dMaxPercentOverlapNext < 0.0) || (dMaxPercentOverlapNext > 100.0)) {
		_EXCEPTIONT("--max_overlap_next must take on values between 0 and 100");
	}

	dMinPercentOverlapPrev /= 100.0;
	dMaxPercentOverlapPrev /= 100.0;
	dMinPercentOverlapNext /= 100.0;
	dMaxPercentOverlapNext /= 100.0;

	// Ensure --minlat and --maxlat values are in range
	if ((dMinLatDeg < -90.0) || (dMinLatDeg > 90.0)) {
		_EXCEPTION1("--minlat must be in range [-90,90] (given %1.5f)", dMinLatDeg);
	}
	if ((dMaxLatDeg < -90.0) || (dMaxLatDeg > 90.0)) {
		_EXCEPTION1("--maxlat must be in range [-90,90] (given %1.5f)", dMaxLatDeg);
	}

	// Convert longitude values to standard range
	if (dMinLonDeg != 0.0) {
		dMinLonDeg = LonDegToStandardRange(dMinLonDeg);
	}
	if (dMaxLonDeg != 360.0) {
		dMaxLonDeg = LonDegToStandardRange(dMaxLonDeg);
	}
	if (dMinLonDeg == dMaxLonDeg) {
		_EXCEPTION2("--minlon (%1.5f) and --maxlon (%1.5f) must not be equal",
			dMinLonDeg, dMaxLonDeg);
	}
	if (dMinLatDeg >= dMaxLonDeg) {
		_EXCEPTION2("--minlat (%1.5f) must be strictly less than --maxlat (%1.5f)",
			dMinLatDeg, dMaxLatDeg);
	}

	// Parse the restrict region string
	RestrictRegion opRestrictRegion;
	if (strRestrictRegion != "") {
		opRestrictRegion.Parse(strRestrictRegion);
	}

	// Parse the threshold string
	std::vector<BlobThresholdOp> vecThresholdOp;

	if (strThresholdCmd != "") {
		AnnounceStartBlock("Parsing thresholds");

		int iLast = 0;
		for (int i = 0; i <= strThresholdCmd.length(); i++) {

			if ((i == strThresholdCmd.length()) ||
			    (strThresholdCmd[i] == ';')
			) {
				std::string strSubStr =
					strThresholdCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecThresholdOp.size());
				vecThresholdOp.resize(iNextOp + 1);
				vecThresholdOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		// Load in the benchmark file
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[0]);
		_ASSERT(vecNcFiles.size() > 0);

		if (vecNcFiles.size() < 1) {
			_EXCEPTIONT("No data files specified; unable to generate grid");
		}

		grid.GenerateLatitudeLongitude(
			vecNcFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		_ASSERT(grid.m_nGridDim.size() == 2);

		AnnounceEndBlock("Done");
	}

	// Check for area
	if (!grid.HasAreas()) {
		_EXCEPTIONT("SimpleGrid has no area information (needed for StitchBlobs)");
	}

	// Get time dimension over all files
	// strOutTimeUnits is either predetermined or set at the command line
	NcType nctypeTime;
	std::vector< std::pair<int, int> > vecGlobalTimeIxToFileTimeIx;

	std::vector< std::vector<Time> > vecGlobalTimes;
	vecGlobalTimes.resize(vecOutputFiles.size());

	for (int f = 0; f < vecInputFiles.size(); f++){

		// Load in the benchmark file
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[f]);
		_ASSERT(vecNcFiles.size() > 0);

		// Get the time variable
		NcVar * varTime = vecNcFiles[0]->get_var("time");
		if (varTime == NULL) {
			_EXCEPTION1("File \"%s\" does not contain \"time\" variable",
				vecNcFiles.GetFilename(0).c_str());
		}

		nctypeTime = varTime->type();

		// Get time units (if not specified on command line)
		if (strOutTimeUnits == "") {
			NcAtt * attTime = varTime->get_att("units");
			if (attTime == NULL) {
				_EXCEPTION1("File \"%s\" missing \"time:units\" attribute",
					vecNcFiles.GetFilename(0).c_str());
			}

			strOutTimeUnits = attTime->as_string(0);
		}

		// Load in CF-compliant time data
		std::vector<Time> vecTimes;
		ReadCFTimeDataFromNcFile(
			vecNcFiles[0],
			vecInputFiles[f].c_str(),
			vecTimes,
			true);

		std::cout << vecTimes[0].ToString() << " " << vecTimes[0].GetCalendarType() << std::endl;

		if (vecOutputFiles.size() == 1) {
			for (int t = 0; t < vecTimes.size(); t++) {
				int iGlobalTime = vecGlobalTimes[0].size();
				vecGlobalTimes[0].push_back(vecTimes[t]);
				vecGlobalTimeIxToFileTimeIx.push_back( std::pair<int,int>(0,iGlobalTime) );
			}
		} else {
			for (int t = 0; t < vecTimes.size(); t++) {
				vecGlobalTimes[f].push_back(vecTimes[t]);
				vecGlobalTimeIxToFileTimeIx.push_back( std::pair<int,int>(f,t) );
			}
		}
	}

	int nGlobalTimes = 0;
	for (int f = 0; f < vecGlobalTimes.size(); f++) {
		nGlobalTimes += vecGlobalTimes[f].size();
	}
	_ASSERT(nGlobalTimes > 0);
	_ASSERT(nGlobalTimes == vecGlobalTimeIxToFileTimeIx.size());

	///////////////////////////////////////////////////////////////////////////
	// Build the set of nodes at each time contained in each blob
	///////////////////////////////////////////////////////////////////////////

	// Build blobs at each time level
	AnnounceStartBlock("Building blob set at each time level");

	// Set of nodes at each time contained in each blob
	std::vector< std::vector<IndicatorSet> > vecAllBlobs;
	vecAllBlobs.resize(nGlobalTimes);

	// Bounding boxes at each time for each blob
	std::vector< std::vector< LatLonBox<double> > > vecAllBlobBoxesDeg;
	vecAllBlobBoxesDeg.resize(nGlobalTimes);

	// Time index across all files
	int iTime = 0;

	// Loop through all files
	for (int f = 0; f < nFiles; f++) {

		// Clear existing data in the register
		varreg.UnloadAllGridData();

		// Load in the benchmark file
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[f]);
		_ASSERT(vecNcFiles.size() > 0);

		// Number of times in this input file
		NcDim * dimTimeInput = vecNcFiles[0]->get_dim("time");
		if (dimTimeInput == NULL) {
			_EXCEPTION1("No dimension \"time\" in file \"%s\"",
				vecNcFiles.GetFilename(0).c_str());
		}
		int nLocalTimes = dimTimeInput->size();

		// Loop through all times
		for (int t = 0; t < nLocalTimes; t++, iTime++) {

			// Get the current patch vector
			std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[iTime];

			std::vector< LatLonBox<double> > & vecBlobBoxesDeg = vecAllBlobBoxesDeg[iTime];

			// New announcement block for timestep
			if (vecGlobalTimes.size() == 1) {
				_ASSERT((iTime >= 0) && (iTime < vecGlobalTimes[0].size()));
				AnnounceStartBlock("Time %i (%s)", iTime,
					vecGlobalTimes[0][iTime].ToString().c_str());
			} else {
				_ASSERT((t >= 0) && (t < vecGlobalTimes[f].size()));
				AnnounceStartBlock("Time %i (%s)", iTime,
					vecGlobalTimes[f][t].ToString().c_str());
			}

			// Load the search variable data
			Variable & var = varreg.Get(varix);
			vecNcFiles.SetConstantTimeIx(t);
			var.LoadGridData(varreg, vecNcFiles, grid);
			const DataArray1D<float> & dataIndicator = var.GetData();
/*
			float dChecksum = 0.0;
			for (int i = 0; i < dataState.GetRows(); i++) {
				dChecksum += dataState[i];
			}
			std::cout << dChecksum << std::endl;
*/

			// Buffer variables used for looping through indicators
			IndicatorSet setIndicators;

			// Insert all detected locations into set
			// (accounting for points out of range)
			if ((dMinLatDeg != -90.0) ||
				(dMaxLatDeg != 90.0) ||
			    (dMinLonDeg != 0.0) ||
				(dMaxLonDeg != 360.0)
			) {
				std::cout << dMinLatDeg << std::endl;
				std::cout << dMaxLatDeg << std::endl;
				std::cout << dMinLonDeg << std::endl;
				std::cout << dMaxLonDeg << std::endl;

				LatLonBox<double> boxBoundsDeg(fRegional);
				boxBoundsDeg.lon[0] = dMinLonDeg;
				boxBoundsDeg.lon[1] = dMaxLonDeg;
				boxBoundsDeg.lat[0] = dMinLatDeg;
				boxBoundsDeg.lat[1] = dMaxLatDeg;

				for (int i = 0; i < grid.GetSize(); i++) {
					if (dataIndicator[i] == 0.0f) {
						continue;
					}

					double dLonDeg = RadToDeg(grid.m_dLon[i]);
					double dLatDeg = RadToDeg(grid.m_dLat[i]);

					dLonDeg = LonDegToStandardRange(dLonDeg);

					if (!boxBoundsDeg.contains(dLatDeg, dLonDeg)) {
						continue;
					}

					setIndicators.insert(i);
				}

			// Insert all detected locations into set
			// (no bounds checking)
			} else {
				for (int i = 0; i < grid.GetSize(); i++) {
					if (dataIndicator[i] != 0.f) {
						setIndicators.insert(i);
					}
				}
			}

			Announce("Tagged points: %i", setIndicators.size());

			// Rejections due to insufficient node count
			int nRejectedMinSize = 0;

			DataArray1D<int> nRejectedThreshold(vecThresholdOp.size());

			// Find all patches
			for (; setIndicators.size() != 0;) {

				// Next starting location
				int ixNode = *(setIndicators.begin());

				// Current patch
				int ixBlob = vecBlobs.size();
				vecBlobs.resize(ixBlob+1);
				vecBlobBoxesDeg.resize(ixBlob+1, LatLonBox<double>(fRegional));

				// Initialize bounding box
				LatLonBox<double> & boxDeg = vecBlobBoxesDeg[ixBlob];
				boxDeg.lat[0] = RadToDeg(grid.m_dLat[ixNode]);
				boxDeg.lat[1] = RadToDeg(grid.m_dLat[ixNode]);
				boxDeg.lon[0] = LonDegToStandardRange(RadToDeg(grid.m_dLon[ixNode]));
				boxDeg.lon[1] = LonDegToStandardRange(RadToDeg(grid.m_dLon[ixNode]));

				// Find all connecting nodes in patch
				IndicatorSet setNeighbors;
				setNeighbors.insert(ixNode);
				while (setNeighbors.size() != 0) {
					ixNode = *(setNeighbors.begin());
					setNeighbors.erase(setNeighbors.begin());

					// This node is already included in the blob
					if (vecBlobs[ixBlob].find(ixNode) != vecBlobs[ixBlob].end()) {
						continue;
					}

					// This node has not been tagged
					IndicatorSetIterator iterIndicator = setIndicators.find(ixNode);
					if (iterIndicator == setIndicators.end()) {
						continue;
					}

					// Remove this from the set of available indicators
					setIndicators.erase(iterIndicator);

					// Insert the node into the blob
					vecBlobs[ixBlob].insert(ixNode);

					// Update bounding box
					boxDeg.insert(
						RadToDeg(grid.m_dLat[ixNode]),
						LonDegToStandardRange(RadToDeg(grid.m_dLon[ixNode])));

					// Insert neighbors
					for (int i = 0; i < grid.m_vecConnectivity[ixNode].size(); i++) {
						setNeighbors.insert(
							grid.m_vecConnectivity[ixNode][i]);
					}
				}

				// Check blob size
				if (vecBlobs[ixBlob].size() < nMinBlobSize) {
					nRejectedMinSize++;
					vecBlobs.resize(ixBlob);
					vecBlobBoxesDeg.resize(ixBlob);

				// Check other thresholds
				} else {
					for (int x = 0; x < vecThresholdOp.size(); x++) {

						bool fSatisfies =
							vecThresholdOp[x].Apply(
								grid,
								vecBlobs[ixBlob],
								boxDeg);

						if (!fSatisfies) {
							nRejectedThreshold[x]++;
							vecBlobs.resize(ixBlob);
							vecBlobBoxesDeg.resize(ixBlob);
							break;
						}
					}
				}
			}

			Announce("Blobs detected: %i", vecBlobs.size());
			Announce("Rejected (min size): %i", nRejectedMinSize);
			for (int x = 0; x < vecThresholdOp.size(); x++) {
				Announce("Rejected (threshold %i): %i",
					x, nRejectedThreshold[x]);
			}

			if (fVerbose) {
				for (int p = 0; p < vecBlobBoxesDeg.size(); p++) {
					Announce("Blob %i [%1.5f, %1.5f] x [%1.5f, %1.5f]",
						p+1,
						vecBlobBoxesDeg[p].lat[0],
						vecBlobBoxesDeg[p].lat[1],
						vecBlobBoxesDeg[p].lon[0],
						vecBlobBoxesDeg[p].lon[1]);
				}
			}

			AnnounceEndBlock("Done");
		}
	}

	AnnounceEndBlock("Done");

	///////////////////////////////////////////////////////////////////////////
	// Stitch blobs together in time using graph search
	///////////////////////////////////////////////////////////////////////////

	AnnounceStartBlock("Assign local tags to each blob");

	// Tags for each blob at each time slice
	std::vector< std::vector<Tag> > vecAllBlobTags;
	vecAllBlobTags.resize(nGlobalTimes);

	// Next available patch tag
	Tag tagNextBlob(1, 0);

	// Give blob tags to the initial set of blobs
	std::set<Tag> setAllTags;

	for (int t = 0; t < nGlobalTimes; t++) {
		vecAllBlobTags[t].resize(vecAllBlobs[t].size());

		tagNextBlob.id = 0;
		for (int p = 0; p < vecAllBlobTags[t].size(); p++) {
			vecAllBlobTags[t][p] = tagNextBlob;
			setAllTags.insert(tagNextBlob);

			tagNextBlob.id++;
		}
		tagNextBlob.time++;
	}

	AnnounceEndBlock("Done");

	// Array of equivalent tags
	typedef std::multimap<Tag, Tag> MapGraph;
	typedef MapGraph::const_iterator MapGraphConstIterator;
	typedef MapGraph::iterator MapGraphIterator;

	MapGraph multimapTagGraph;

	// Set of Tags that are within the restrict region
	std::set<Tag> setRestrictRegion;

	AnnounceStartBlock("Building connectivity graph");

	// Loop through all remaining time steps
	int iFileLocal = 0;
	int iTimeLocal = 0;
	for (int t = 1; t < nGlobalTimes; t++) {

		// New announcement block for timestep
		if (vecGlobalTimes.size() == 1) {
			_ASSERT((t >= 0) && (t < vecGlobalTimes[0].size()));
			Announce("Time %i (%s)", t,
				vecGlobalTimes[0][t].ToString().c_str());
		} else {
			iTimeLocal++;
			_ASSERT(iFileLocal < vecGlobalTimes.size());
			if (iTimeLocal >= vecGlobalTimes[iFileLocal].size()) {
				iTimeLocal = 0;
				iFileLocal++;
				_ASSERT(iFileLocal < vecGlobalTimes.size());
			}
			_ASSERT(iTimeLocal < vecGlobalTimes[iFileLocal].size());
			Announce("Time %i (%s)", t,
				vecGlobalTimes[iFileLocal][iTimeLocal].ToString().c_str());
		}

		// Get the current blob vector
		const std::vector<Tag> & vecPrevBlobTags = vecAllBlobTags[t-1];

		std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		const std::vector<IndicatorSet> & vecPrevBlobs = vecAllBlobs[t-1];

		const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

		const std::vector< LatLonBox<double> > & vecPrevBlobBoxesDeg
			= vecAllBlobBoxesDeg[t-1];

		const std::vector< LatLonBox<double> > & vecBlobBoxesDeg
			= vecAllBlobBoxesDeg[t];

		// Determine overlaps between these blobs and previous blobs
		vecBlobTags.resize(vecBlobs.size());
		int nCountRemove = 0;
		for (int p = 0; p < vecBlobTags.size(); p++) {

			const LatLonBox<double> & boxP = vecBlobBoxesDeg[p];

			// Find overlap with bounding boxes at previous time
			for (int q = 0; q < vecPrevBlobTags.size(); q++) {

				const LatLonBox<double> & boxQ = vecPrevBlobBoxesDeg[q];

				// Check if bounding boxes overlap
				if (!boxP.overlaps(boxQ)) {
					continue;
				}

				// Verify that at least one node overlaps between blobs
				{
					bool fHasOverlapNode = false;
					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						if (vecPrevBlobs[q].find(*iter) !=
							vecPrevBlobs[q].end()
						) {
							fHasOverlapNode = true;
							break;
						}
					}

					if (!fHasOverlapNode) {
						continue;
					}
				}

				// Verify that blobs meet percentage overlap criteria
				// (as a percentage of the current blob)
				if ((dMinPercentOverlapNext != 0.0) ||
				    (dMaxPercentOverlapNext != 1.0)
				) {
					double dCurrentArea = 0.0;
					double dOverlapArea = 0.0;

					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						dCurrentArea += grid.m_dArea[*iter];
						if (vecPrevBlobs[q].find(*iter) != vecPrevBlobs[q].end()) {
							dOverlapArea += grid.m_dArea[*iter];
						}
					}
					if (dCurrentArea == 0.0) {
						_EXCEPTIONT("Logic error (zero area blob)");
					}
					if (dOverlapArea == 0.0) {
						_EXCEPTIONT("Logic error (zero overlap area)");
					}
					if (dOverlapArea < dCurrentArea * dMinPercentOverlapNext) {
						continue;
					}
					if (dOverlapArea > dCurrentArea * dMaxPercentOverlapNext) {
						continue;
					}
				}

				// Verify the blobs meet percentage overlap criteria
				// (as a percentage of the earlier blob)
				if ((dMinPercentOverlapPrev != 0.0) ||
				    (dMaxPercentOverlapPrev != 1.0)
				) {
					double dCurrentArea = 0.0;
					double dOverlapArea = 0.0;

					IndicatorSetConstIterator iter = vecPrevBlobs[q].begin();
					for (; iter != vecPrevBlobs[q].end(); iter++) {
						dCurrentArea += grid.m_dArea[*iter];
						if (vecBlobs[p].find(*iter) != vecBlobs[p].end()) {
							dOverlapArea += grid.m_dArea[*iter];
						}
					}
					if (dCurrentArea == 0.0) {
						_EXCEPTIONT("Logic error (zero area blob)");
					}
					if (dOverlapArea == 0.0) {
						_EXCEPTIONT("Logic error (zero overlap area)");
					}
					if (dOverlapArea < dCurrentArea * dMinPercentOverlapPrev) {
						continue;
					}
					if (dOverlapArea > dCurrentArea * dMaxPercentOverlapPrev) {
						continue;
					}
				}

				// Check restrict_region criteria
				if (opRestrictRegion.IsActive()) {
					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						bool fInRestrictRegion =
							opRestrictRegion.ContainsPoint(
								grid.m_dLat[*iter],
								grid.m_dLon[*iter]);

						if (fInRestrictRegion) {
							setRestrictRegion.insert(vecBlobTags[p]);
							break;
						}
					}
				}

				// Insert bidirectional edge in graph
				multimapTagGraph.insert(
					std::pair<Tag, Tag>(
						vecBlobTags[p], vecPrevBlobTags[q]));

				multimapTagGraph.insert(
					std::pair<Tag, Tag>(
						vecPrevBlobTags[q], vecBlobTags[p]));
			}
		}
	}

	AnnounceEndBlock("Done");

	AnnounceStartBlock("Identify components in connectivity graph");

	// Total number of blobs
	int nTotalBlobCount = 0;

	// Find all components using a bidirectional graph search
	std::map<Tag, Tag> mapEquivalentTags;

	std::set<Tag>::const_iterator iterTag = setAllTags.begin();

	for (; iterTag != setAllTags.end(); iterTag++) {

		std::set<Tag> setTagsVisited;

		std::set<Tag> setTagsToVisit;
		setTagsToVisit.insert(*iterTag);

		Tag tagMinimum = *iterTag;

		// Check if this tag is already part of an explored component
		if (mapEquivalentTags.find(*iterTag) != mapEquivalentTags.end()) {
			continue;
		}

		// New blob
		nTotalBlobCount++;

		// All time indices associated with this blob
		std::set<int> setBlobTimes;

		// All time indices where this blob is in the restrict region
		std::set<int> setBlobTimesInRestrictRegion;

		// Find component containing this node
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

			setBlobTimes.insert(tagNext.time);

			if (setRestrictRegion.find(tagNext) != setRestrictRegion.end()) {
				setBlobTimesInRestrictRegion.insert(tagNext.time);
			}

			// Get edges from this node
			std::pair<MapGraphIterator, MapGraphIterator> iterGraphEdges
				= multimapTagGraph.equal_range(tagNext);

			MapGraphIterator iter = iterGraphEdges.first;
			for (; iter != iterGraphEdges.second; iter++) {
				setTagsToVisit.insert(iter->second);
			}
		}

		// Apply filters to the blob
		bool fAcceptBlob = true;

		// Filter on RestrictRegion count for this global_id
		if (opRestrictRegion.IsActive()) {
			if (setBlobTimesInRestrictRegion.size() < opRestrictRegion.GetMinimumCount()) {
				fAcceptBlob = false;
			}
		}

		// Filter on min_time
		if (setBlobTimes.size() < nMinTime) {
			fAcceptBlob = false;
		}
		if (dMinTimeSeconds > 0.0) {
			int iFirstTime = *(setBlobTimes.begin());
			int iLastTime = *(setBlobTimes.rbegin());

			_ASSERT(iFirstTime >= 0);
			_ASSERT(iLastTime < vecGlobalTimeIxToFileTimeIx.size());

			std::pair<int,int> iFileTimeIxFirst = vecGlobalTimeIxToFileTimeIx[iFirstTime];
			std::pair<int,int> iFileTimeIxLast = vecGlobalTimeIxToFileTimeIx[iLastTime];

			const Time & timeFirst = vecGlobalTimes[iFileTimeIxFirst.first][iFileTimeIxFirst.second];
			const Time & timeLast = vecGlobalTimes[iFileTimeIxLast.first][iFileTimeIxLast.second];

			double dBlockDurationSeconds = timeFirst.DeltaSeconds(timeLast);
			if (dBlockDurationSeconds < dMinTimeSeconds) {
				fAcceptBlob = false;
			}
		}

		// Set global id
		if (fAcceptBlob) {
			tagMinimum.global_id = nTotalBlobCount;
		} else {
			tagMinimum.global_id = 0;
			nTotalBlobCount--;
		}

		// Refer all tags in component to minimum tag
		std::set<Tag>::const_iterator iterTagsVisited = setTagsVisited.begin();
		for (; iterTagsVisited != setTagsVisited.end(); iterTagsVisited++) {
			mapEquivalentTags.insert(
				std::pair<Tag,Tag>(*iterTagsVisited, tagMinimum));
		}
	}

	AnnounceEndBlock("Done");

	// Merge blobs at each time step with equivalent tags
	AnnounceStartBlock("Reassign blob tags");

	for (int t = 0; t < nGlobalTimes; t++) {

		std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		for (int p = 0; p < vecBlobTags.size(); p++) {
			std::map<Tag, Tag>::const_iterator iterTagPair =
				mapEquivalentTags.find(vecBlobTags[p]);

			if (iterTagPair != mapEquivalentTags.end()) {
				vecBlobTags[p] = iterTagPair->second;
			}
		}
	}

	Announce("Unique tags found: %i", nTotalBlobCount);

	AnnounceEndBlock("Done");
/*
	// Apply post-hoc threshold operators
	std::vector<bool> fRejectedBlob;
	fRejectedBlob.resize(nTotalBlobCount+1, false);

	if (vecThresholdOp.size() != 0) {
		AnnounceStartBlock("Applying threshold commands");

		// Whether or not each blob (global_id) satisfies each threshold
		// at each time
		std::vector< std::vector< std::vector<bool> > >
			vecBlobSatisfiesThreshold;

		vecBlobSatisfiesThreshold.resize(nTotalBlobCount+1);
		for (int gid = 1; gid <= nTotalBlobCount; gid++) {
			vecBlobSatisfiesThreshold[gid].resize(
				vecThresholdOp.size());
		}

		// Determine if blobs satisfy the threshold at each time
		for (int t = 0; t < nTime; t++) {

			std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

			std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

			std::vector<LatLonBox> & vecBlobBoxesDeg = vecAllBlobBoxesDeg[t];

			for (int x = 0; x < vecThresholdOp.size(); x++) {

				std::map<int, bool> mapSatisfiesThreshold;

				vecThresholdOp[x].Apply(
					dCellArea,
					dataLatDeg,
					dataLonDeg,
					vecBlobTags,
					vecBlobs,
					vecBlobBoxesDeg,
					mapSatisfiesThreshold);

				std::map<int,bool>::const_iterator iter =
					mapSatisfiesThreshold.begin();

				for (; iter != mapSatisfiesThreshold.end(); iter++) {
					int iGlobalBlobIndex = iter->first;

					if ((iGlobalBlobIndex < 1) ||
					    (iGlobalBlobIndex > nTotalBlobCount)
					) {
						_EXCEPTION1("global_id out of range (%i)",
							iGlobalBlobIndex);
					}
					vecBlobSatisfiesThreshold[iGlobalBlobIndex][x].
						push_back(iter->second);
				}
			}
		}

		// Number of blobs rejected for each reason
		std::vector<int> nBlobsRejected;
		nBlobsRejected.resize(vecThresholdOp.size());

		for (int gid = 1; gid <= nTotalBlobCount; gid++) {
		for (int x = 0; x < vecThresholdOp.size(); x++) {
			int nCount = 0;
			for (int t = 0; t < vecBlobSatisfiesThreshold[gid][x].size(); t++) {
				if (vecBlobSatisfiesThreshold[gid][x][t]) {
					nCount++;
				}
			}

			// Criteria must be satisfied at all times
			if (vecThresholdOp[x].GetMinimumCount() == (-1)) {
				if (nCount != vecBlobSatisfiesThreshold[gid][x].size()) {
					fRejectedBlob[gid] = true;
					nBlobsRejected[x]++;
					break;
				}

			// Criteria must be satisfied at a minimum number of times
			} else if (nCount < vecThresholdOp[x].GetMinimumCount()) {
				fRejectedBlob[gid] = true;
				nBlobsRejected[x]++;
				break;
			}
		}
		}

		// Announce rejections
		for (int x = 0; x < vecThresholdOp.size(); x++) {
			Announce("Rejected (threshold %i): %i", x, nBlobsRejected[x]);
		}

		AnnounceEndBlock("Done");
	}
*/


	///////////////////////////////////////////////////////////////////////////
	// Output results
	///////////////////////////////////////////////////////////////////////////

	AnnounceStartBlock("Output blobs");

	{
		// Load in the benchmark file
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[0]);
		_ASSERT(vecNcFiles.size() > 0);

		// Loop through all output files
		_ASSERT(vecOutputFiles.size() == vecGlobalTimes.size());

		int iGlobalTimeIx = 0;

		for (int f = 0; f < vecOutputFiles.size(); f++) {

			Announce("Writing file \"%s\"", vecOutputFiles[f].c_str());

			// Open output file
			NcFile ncOutput(vecOutputFiles[f].c_str(), NcFile::Replace);
			if (!ncOutput.is_valid()) {
				_EXCEPTION1("Unable to open output file \"%s\"",
					vecOutputFiles[f].c_str());
			}

			// Output time dimension
			int nLocalTimes = vecGlobalTimes[f].size();

			NcDim * dimOutputTime = ncOutput.add_dim("time", nLocalTimes);
			if (dimOutputTime == NULL) {
				_EXCEPTIONT("Unable to create dimension \"time\" in output");
			}
			NcVar * varOutputTime =
				ncOutput.add_var("time", ncDouble, dimOutputTime);

			DataArray1D<double> dOutputTimes(nLocalTimes);
			for (int t = 0; t < vecGlobalTimes[f].size(); t++) {
				dOutputTimes[t] =
					vecGlobalTimes[f][t].GetCFCompliantUnitsOffsetDouble(strOutTimeUnits);
			}

			varOutputTime->add_att("long_name","time");
			varOutputTime->add_att("units",strOutTimeUnits.c_str());
			varOutputTime->add_att("calendar",vecGlobalTimes[f][0].GetCalendarName().c_str());

			varOutputTime->put(&(dOutputTimes[0]), nLocalTimes);

			// Create output variable
			NcDim * dimOut0 = NULL;
			NcDim * dimOut1 = NULL;
			NcVar * varTagOut = NULL;

			PrepareBlobOutputVar(
				*(vecNcFiles[0]),
				ncOutput,
				vecOutputFiles[f],
				grid,
				strOutputVariable,
				strLatitudeName,
				strLongitudeName,
				ncInt,
				dimOutputTime,
				&dimOut0,
				&dimOut1,
				&varTagOut);

			_ASSERT(varTagOut != NULL);

			int nDimOutSize0 = 0;
			int nDimOutSize1 = 0;

			if (dimOut0 != NULL) {
				nDimOutSize0 = dimOut0->size();
			}
			if (dimOut1 != NULL) {
				nDimOutSize1 = dimOut1->size();
			}

			// Write all time steps
			DataArray1D<int> dataBlobTag(grid.GetSize());

			for (int t = 0; t < nLocalTimes; t++) {

				dataBlobTag.Zero();

				_ASSERT(iGlobalTimeIx + t < vecAllBlobTags.size());
				_ASSERT(iGlobalTimeIx + t < vecAllBlobs.size());
	
				// Get the current blob vectors
				const std::vector<Tag> & vecBlobTags = vecAllBlobTags[iGlobalTimeIx + t];
				const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[iGlobalTimeIx + t];

				_ASSERT(vecBlobTags.size() == vecBlobs.size());

				// Put blob information into dataBlobTag
				for (int p = 0; p < vecBlobTags.size(); p++) {

					if (vecBlobTags[p].global_id == 0) {
						continue;
					}
		
					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						dataBlobTag[*iter] =
							vecBlobTags[p].global_id;
					}
				}

				// Write to file
				if (grid.DimCount() == 1) {
					varTagOut->set_cur(t, 0);
					varTagOut->put(&(dataBlobTag[0]), 1, nDimOutSize0);
				} else {
					varTagOut->set_cur(t, 0);
					varTagOut->put(&(dataBlobTag[0]), 1, nDimOutSize0, nDimOutSize1);
				}
			}

			// Update global time index
			iGlobalTimeIx += nLocalTimes;

			// Close the output file
			ncOutput.close();

		}

		AnnounceEndBlock("Done");
	}
/*
	// Copy variable attributes from first input file
	{
		NcFile ncInput(vecInputFiles[0].c_str());

		NcVar * varLat = ncInput.get_var(strLatitudeName.c_str());
		NcVar * varLon = ncInput.get_var(strLongitudeName.c_str());

		CopyNcVarAttributes(varLat, varOutputLat);
		CopyNcVarAttributes(varLon, varOutputLon);

		NcVar * varTime = ncInput.get_var(strTimeName.c_str());
		if (varTime != NULL) {
			CopyNcVarAttributes(varTime, varOutputTime);
		}
	}
*/
/*
	// Change the time units to the preset
	NcAtt * oldUnits = varOutputTime->get_att("units");
	oldUnits->remove();
	varOutputTime->add_att("units",strOutTimeUnits.c_str());
	NcVar * varData =
		ncOutput.add_var(
			strOutputVariable.c_str(),
			ncInt,
			dimOutputTime,
			dimOutputLat,
			dimOutputLon);

	// Loop through all time steps
	DataArray2D<int> dataBlobTag(nLat, nLon);

	for (int t = 0; t < nTime; t++) {

		dataBlobTag.Zero();

		// Get the current blob vectors
		const std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

		// Put blob information into matrix
		for (int p = 0; p < vecBlobTags.size(); p++) {

			if (vecBlobTags[p].global_id == 0) {
				continue;
			}

			//if (fRejectedBlob[vecBlobTags[p].global_id]) {
			//	continue;
			//}

			IndicatorSetConstIterator iter = vecBlobs[p].begin();
			for (; iter != vecBlobs[p].end(); iter++) {
				dataBlobTag[iter->lat][iter->lon] =
					vecBlobTags[p].global_id;
			}
		}

		// Write to file
		varData->set_cur(t, 0, 0);
		varData->put(&(dataBlobTag[0][0]), 1, nLat, nLon);
	}
*/

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

