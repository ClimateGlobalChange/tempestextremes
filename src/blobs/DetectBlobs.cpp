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

#include "BlobUtilities.h"

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

static const double EarthRadius = 6.37122e6;

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
typedef std::set<LatLonPair> IndicatorSet;
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
		EastwardOrientation,
		WestwardOrientation
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

					} else if (strSubStr == "eastwardorientation") {
						m_eQuantity = EastwardOrientation;
						eReadMode = ReadMode_Invalid;

					} else if (strSubStr == "westwardorientation") {
						m_eQuantity = WestwardOrientation;
						eReadMode = ReadMode_Invalid;

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
		} else if (m_eQuantity == EastwardOrientation) {
			strDescription += "Eastward orientation ";
		} else if (m_eQuantity == WestwardOrientation) {
			strDescription += "Westward orientation ";
		}
/*
		char szMinCount[160];
		if (m_nMinimumCount == -1) {
			strDescription += " at all times";
		} else {
			sprintf(szMinCount, " at least %i time(s)", m_nMinimumCount);
			strDescription += szMinCount;
		}
*/
		Announce("%s", strDescription.c_str());
	}

	///	<summary>
	///		Verify that the specified path satisfies the threshold op.
	///	</summary>
	bool Apply(
		const DataVector<double> & dCellArea,
		const DataVector<double> & dLatDeg,
		const DataVector<double> & dLonDeg,
		const IndicatorSet & setBlobPoints,
		const LatLonBox & boxBlob
	) {
		// Number of longitudes
		int nLonCount = dLonDeg.GetRows();

		// Thresholds related to area
		if ((m_eQuantity == MinArea) ||
		    (m_eQuantity == MaxArea) ||
		    (m_eQuantity == MinArealFraction) ||
		    (m_eQuantity == MaxArealFraction)
		) {
			double dBoxArea = 0.0;
			double dBlobArea = 0.0;

			// Calculate the area of each blob box
			for (int i = boxBlob.lat[0]; i <= boxBlob.lat[1]; i++) {
				dBoxArea += dCellArea[i]
					* static_cast<double>(boxBlob.Width(nLonCount));
			}

			// Calculate the area and mean lat/lon of each blob
			IndicatorSetIterator iterBlob = setBlobPoints.begin();
			for (; iterBlob != setBlobPoints.end(); iterBlob++) {
				dBlobArea += dCellArea[iterBlob->lat];
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

		// Invalid quantity
		} else {
			_EXCEPTIONT("Invalid quantity");
		}

		return true;
	}
/*
	///	<summary>
	///		Verify that the specified path satisfies the threshold op.
	///	</summary>
	void Apply(
		const DataVector<double> & dCellArea,
		const DataVector<double> & dLatDeg,
		const DataVector<double> & dLonDeg,
		const std::vector<Tag> & vecBlobTags,
		const std::vector<IndicatorSet> & vecBlobs,
		const std::vector<LatLonBox> & vecBlobBoxes,
		std::map<int, bool> & mapSatisfiesThreshold
	) {
		if (vecBlobTags.size() != vecBlobs.size()) {
			_EXCEPTIONT("ASSERT: vecBlobTags.size() != vecBlobs.size()");
		}
		if (vecBlobTags.size() != vecBlobBoxes.size()) {
			_EXCEPTIONT("ASSERT: vecBlobTags.size() != vecBlobBoxes.size()");
		}

		mapSatisfiesThreshold.clear();

		// Number of longitudes
		int nLonCount = dLonDeg.GetRows();

		// Thresholds related to area
		if ((m_eQuantity == MinArea) ||
		    (m_eQuantity == MaxArea) ||
		    (m_eQuantity == MinArealFraction) ||
		    (m_eQuantity == MaxArealFraction)
		) {
			// Loop through all blobs
			for (int p = 0; p < vecBlobTags.size(); p++) {

				double dBoxArea = 0.0;
				double dBlobArea = 0.0;

				// Obtain an iterator to this blob in the satisfies
				// threshold map.
				std::map<int, bool>::iterator iterSatisfies =
					mapSatisfiesThreshold.find(vecBlobTags[p].id);

				if (iterSatisfies == mapSatisfiesThreshold.end()) {
					iterSatisfies = mapSatisfiesThreshold.insert(
						std::pair<int,bool>(
							vecBlobTags[p].id, true)).first;

				} else if (!iterSatisfies->second) {
					continue;
				}

				// Calculate the area of each blob box
				for (int i = vecBlobBoxes[p].lat[0];
				     i <= vecBlobBoxes[p].lat[1]; i++
				) {
					dBoxArea += dCellArea[i]
						* static_cast<double>(
							vecBlobBoxes[p].Width(nLonCount));
				}

				// Calculate the area and mean lat/lon of each blob
				IndicatorSetIterator iterBlob = vecBlobs[p].begin();
				for (; iterBlob != vecBlobs[p].end(); iterBlob++) {
					dBlobArea += dCellArea[iterBlob->lat];
				}

				// Minimum area
				if (m_eQuantity == MinArea) {
					if (dBlobArea < m_dValue) {
						iterSatisfies->second = false;
					}

				// Maximum area
				} else if (m_eQuantity == MaxArea) {
					if (dBlobArea > m_dValue) {
						iterSatisfies->second = false;
					}

				// Minimum areal fraction
				} else if (m_eQuantity == MinArealFraction) {
					if (dBlobArea < m_dValue * dBoxArea) {
						iterSatisfies->second = false;
					}

				// Maximum areal fraction
				} else if (m_eQuantity == MaxArealFraction) {
					if (dBlobArea > m_dValue * dBoxArea) {
						iterSatisfies->second = false;
					}
				}
			}

		// Thresholds related to orientation
		} else if ((m_eQuantity == EastwardOrientation) ||
		           (m_eQuantity == WestwardOrientation)
		) {
			// Loop through all blobs
			for (int p = 0; p < vecBlobTags.size(); p++) {

				double dNorthHemiMeanLat = 0.0;
				double dNorthHemiMeanLon = 0.0;
				double dNorthHemiMeanLon2 = 0.0;
				double dNorthHemiCoLatLon = 0.0;

				double dSouthHemiMeanLat = 0.0;
				double dSouthHemiMeanLon = 0.0;
				double dSouthHemiMeanLon2 = 0.0;
				double dSouthHemiCoLatLon = 0.0;

				// Obtain an iterator to this blob in the satisfies
				// threshold map.
				std::map<int, bool>::iterator iterSatisfies =
					mapSatisfiesThreshold.find(vecBlobTags[p].id);

				if (iterSatisfies == mapSatisfiesThreshold.end()) {
					iterSatisfies = mapSatisfiesThreshold.insert(
						std::pair<int,bool>(
							vecBlobTags[p].id, true)).first;

				} else if (!iterSatisfies->second) {
					continue;
				}

				// Calculate regression coefficients for this blob
				IndicatorSetIterator iterBlob = vecBlobs[p].begin();
				for (; iterBlob != vecBlobs[p].end(); iterBlob++) {

					double dAltLon = 0.0;
					if (dLatDeg[iterBlob->lat] > 0.0) {
						if (iterBlob->lon < vecBlobBoxes[p].lon[0]) {
							dAltLon = dLonDeg[iterBlob->lon] + 360.0;
						} else {
							dAltLon = dLonDeg[iterBlob->lon];
						}

						dNorthHemiMeanLat += dLatDeg[iterBlob->lat];
						dNorthHemiMeanLon += dAltLon;
						dNorthHemiMeanLon2 += dAltLon * dAltLon;
						dNorthHemiCoLatLon += dLatDeg[iterBlob->lat] * dAltLon;

					} else if (dLatDeg[iterBlob->lat] < 0.0) {
						if (iterBlob->lon < vecBlobBoxes[p].lon[0]) {
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

				double dBlobCount = static_cast<double>(vecBlobs[p].size());

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
						iterSatisfies->second = false;
					}
					if (dSouthSlopeNum * dSouthSlopeDen > 0.0) {
						iterSatisfies->second = false;
					}

				} else if (m_eQuantity == WestwardOrientation) {
					if (dNorthSlopeNum * dNorthSlopeDen > 0.0) {
						iterSatisfies->second = false;
					}
					if (dSouthSlopeNum * dSouthSlopeDen < 0.0) {
						iterSatisfies->second = false;
					}
				}
			}

		// Invalid quantity
		} else {
			_EXCEPTIONT("Invalid quantity");
		}
	}
*/
/*
	///	<summary>
	///		Get the minimum count associated with this threshold.
	///	</summary>
	int GetMinimumCount() const {
		return m_nMinimumCount;
	}
*/
protected:
	///	<summary>
	///		Threshold quantity.
	///	</summary>
	ThresholdQuantity m_eQuantity;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;
/*
	///	<summary>
	///		Minimum number of segments that need to satisfy the op.
	///	</summary>
	int m_nMinimumCount;
*/
};

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

	// Variable name
	std::string strVariable;

	// Output variable name
	std::string strOutputVariable;

	// Minimum blob size (in grid points)
	int nMinBlobSize;

	// Minimum number of timesteps for blob
	int nMinTime;

	// Minimum latitude for detections
	double dMinLat;

	// Maximum latitude for detections
	double dMaxLat;

	// Minimum longitude for detections
	double dMinLon;

	// Maximum longitude for detections
	double dMaxLon;

	// Operate on a regional area (no periodic boundaries)
	bool fRegional;

	// Threshold commands
	std::string strThresholdCmd;

        std::string latname;
	std::string lonname;
	std::string timename;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "inlist", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariable, "var", "");
		CommandLineString(strOutputVariable, "outvar", "");
		CommandLineInt(nMinBlobSize, "minsize", 1);
		CommandLineInt(nMinTime, "mintime", 1);
		CommandLineBool(fRegional, "regional");
		CommandLineDouble(dMinLat, "minlat", -90.0);
		CommandLineDouble(dMaxLat, "maxlat", 90.0);
		CommandLineDouble(dMinLon, "minlon", 0.0);
		CommandLineDouble(dMaxLon, "maxlon", 360.0);
		CommandLineString(latname, "latname", "lat");
		CommandLineString(lonname, "lonname", "lon");
		CommandLineString(timename,"timename", "time");
		CommandLineString(strThresholdCmd, "thresholdcmd", "");

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

	// Load in spatial dimension data
	int nLat;
	int nLon;

	DataVector<double> dataLatDeg;
	DataVector<double> dataLatRad;

	DataVector<double> dataLonDeg;
	DataVector<double> dataLonRad;

	// Cell areas as a function of latitude index
	DataVector<double> dCellArea;

	{
		// Load the first netcdf input file
		NcFile ncInput(vecInputFiles[0].c_str());

		if (!ncInput.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"",
				vecInputFiles[0].c_str());
		}

		// Get latitude/longitude dimensions
		NcDim * dimLat = ncInput.get_dim(latname.c_str());
		if (dimLat == NULL) {
			_EXCEPTIONT("No dimension \"lat\" found in input file");
		}

		NcDim * dimLon = ncInput.get_dim(lonname.c_str());
		if (dimLon == NULL) {
			_EXCEPTIONT("No dimension \"lon\" found in input file");
		}

		NcVar * varLat = ncInput.get_var(latname.c_str());
		if (varLat == NULL) {
			_EXCEPTIONT("No variable \"lat\" found in input file");
		}

		NcVar * varLon = ncInput.get_var(lonname.c_str());
		if (varLon == NULL) {
			_EXCEPTIONT("No variable \"lon\" found in input file");
		}

		nLat = dimLat->size();
		nLon = dimLon->size();

		dataLatDeg.Initialize(nLat);
		dataLatRad.Initialize(nLat);

		dataLonDeg.Initialize(nLon);
		dataLonRad.Initialize(nLon);

		varLat->get(dataLatDeg, nLat);
		for (int j = 0; j < nLat; j++) {
			dataLatRad[j] = dataLatDeg[j] * M_PI / 180.0;
		}

		varLon->get(dataLonDeg, nLon);
		for (int i = 0; i < nLon; i++) {
			dataLonRad[i] = dataLonDeg[i] * M_PI / 180.0;
		}

		// Calculate cell areas in m^2
		double dDeltaLon = 2.0 * M_PI / static_cast<double>(nLon);
		double dDeltaLat = M_PI / static_cast<double>(nLat);

		dCellArea.Initialize(nLat);
		for (int j = 0; j < nLat; j++) {
			dCellArea[j] =
				EarthRadius
				* EarthRadius
				* cos(dataLatRad[j])
				* dDeltaLon
				* dDeltaLat;
		}

		// Close first netcdf file
		ncInput.close();
	}

	// Get time dimension over all files
	DataVector<double> dTime;
	GetAllTimes(vecInputFiles, dTime,timename);

	int nTime = dTime.GetRows();

	// Allocate indicator data
	DataMatrix<int> dataIndicator(nLat, nLon);

	// Build blobs at each time level
	AnnounceStartBlock("Building blob set at each time level");

	// Buffer variables used for looping through indicators
	IndicatorSet setIndicators;
	IndicatorSet setNeighbors;

	///////////////////////////////////////////////////////////////////////////
	// Build the set of nodes at each time contained in each blob
	///////////////////////////////////////////////////////////////////////////

	// Set of nodes at each time contained in each blob
	std::vector< std::vector<IndicatorSet> > vecAllBlobs;
	vecAllBlobs.resize(nTime);

	// Bounding boxes at each time for each blob
	std::vector< std::vector<LatLonBox> > vecAllBlobBoxes;
	vecAllBlobBoxes.resize(nTime);

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
		NcDim * dimTime = ncInput.get_dim(timename.c_str());

		int nLocalTimes = dimTime->size();

		// Load in indicator variable
		NcVar * varIndicator = ncInput.get_var(strVariable.c_str());

		if (varIndicator == NULL) {
			_EXCEPTION1("Unable to load variable \"%s\"",
				strVariable.c_str());
		}

		if (varIndicator->num_dims() != 3) {
			_EXCEPTION1("Incorrect number of dimensions for \"%s\" (3 expected)",
				strVariable.c_str());
		}

		// Loop through all times
		for (int t = 0; t < nLocalTimes; t++, iTime++) {

			// Get the current patch vector
			std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[iTime];

			std::vector<LatLonBox> & vecBlobBoxes = vecAllBlobBoxes[iTime];

			// New announcement block for timestep
			char szStartBlock[128];
			sprintf(szStartBlock, "Time %i (%i)", iTime, t);
			AnnounceStartBlock(szStartBlock);

			// Load in the data at this time
			varIndicator->set_cur(t, 0, 0);
			varIndicator->get(&(dataIndicator[0][0]), 1, nLat, nLon);

			// Elminate detections out of range
			if ((dMinLat != -90.0) || (dMaxLat != 90.0) ||
			    (dMinLon != 0.0) || (dMaxLon != 360.0)
			) {
				for (int j = 0; j < nLat; j++) {
				for (int i = 0; i < nLon; i++) {
					if (dataIndicator[j][i] != 0) {
						if ((dMinLat != -90.0) || (dMaxLat != 90.0)) {
							if (dataLatDeg[j] < dMinLat) {
								dataIndicator[j][i] = 0;
							}
							if (dataLatDeg[j] > dMaxLat) {
								dataIndicator[j][i] = 0;
							}
						}
						if ((dMinLon != 0.0) || (dMaxLon != 360.0)) {
							if (dMinLon < dMaxLon) {
								if (dataLonDeg[i] < dMinLon) {
									dataIndicator[j][i] = 0;
								}
								if (dataLonDeg[i] > dMaxLon) {
									dataIndicator[j][i] = 0;
								}

							} else {
								if ((dataLonDeg[i] < dMinLon) &&
								    (dataLonDeg[i] > dMaxLon)
								) {
									dataIndicator[j][i] = 0;
								}
							}
						}
					}
				}
				}
			}

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

			DataVector<int> nRejectedThreshold;
			nRejectedThreshold.Initialize(vecThresholdOp.size());

			// Find all patches
			for (; setIndicators.size() != 0;) {

				// Next starting location
				LatLonPair pr = *(setIndicators.begin());

				// Current patch
				int ixBlob = vecBlobs.size();
				vecBlobs.resize(ixBlob+1);
				vecBlobBoxes.resize(ixBlob+1);

				// Initialize bounding box
				LatLonBox & box = vecBlobBoxes[ixBlob];

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
					if (vecBlobs[ixBlob].find(pr) !=
						vecBlobs[ixBlob].end()
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
					vecBlobs[ixBlob].insert(pr);

					// Update bounding box
					box.InsertPoint(pr.lat, pr.lon, nLat, nLon);

					// Insert neighbors (regional case)
					if (fRegional) {
						setNeighbors.insert(
							LatLonPair(pr.lat, pr.lon));

						if (pr.lon != 0) {
							setNeighbors.insert(
								LatLonPair(pr.lat, pr.lon-1));
						}
						if (pr.lon != nLon-1) {
							setNeighbors.insert(
								LatLonPair(pr.lat, pr.lon+1));
						}

						if (pr.lat != 0) {
							setNeighbors.insert(
								LatLonPair(pr.lat-1, pr.lon));

							if (pr.lon != 0) {
								setNeighbors.insert(
									LatLonPair(pr.lat-1, pr.lon-1));
							}
							if (pr.lon != nLon-1) {
								setNeighbors.insert(
									LatLonPair(pr.lat-1, pr.lon+1));
							}
						}

						if (pr.lat != nLat-1) {
							setNeighbors.insert(
								LatLonPair(pr.lat+1, pr.lon));

							if (pr.lon != 0) {
								setNeighbors.insert(
									LatLonPair(pr.lat+1, pr.lon-1));
							}
							if (pr.lon != nLon-1) {
								setNeighbors.insert(
									LatLonPair(pr.lat+1, pr.lon+1));
							}
						}

					// Insert neighbors (global case)
					} else if (pr.lat == 0) {
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
				if (vecBlobs[ixBlob].size() < nMinBlobSize) {
					nRejectedMinSize++;
					vecBlobs.resize(ixBlob);
					vecBlobBoxes.resize(ixBlob);

				// Check other thresholds
				} else {
					for (int x = 0; x < vecThresholdOp.size(); x++) {

						bool fSatisfies =
							vecThresholdOp[x].Apply(
								dCellArea,
								dataLatDeg,
								dataLonDeg,
								vecBlobs[ixBlob],
								box);

						if (!fSatisfies) {
							nRejectedThreshold[x]++;
							vecBlobs.resize(ixBlob);
							vecBlobBoxes.resize(ixBlob);
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

			for (int p = 0; p < vecBlobBoxes.size(); p++) {
				Announce("Blob %i [%i, %i] x [%i, %i]",
					p+1,
					vecBlobBoxes[p].lat[0],
					vecBlobBoxes[p].lat[1],
					vecBlobBoxes[p].lon[0],
					vecBlobBoxes[p].lon[1]);
			}

			AnnounceEndBlock("Done");
		}

		// Close NetCDF file
		ncInput.close();
	}

	AnnounceEndBlock("Done");

	///////////////////////////////////////////////////////////////////////////
	// Stitch blobs together in time using graph search
	///////////////////////////////////////////////////////////////////////////

	AnnounceStartBlock("Stitching Blobs");

	// Tags for each blob at each time slice
	std::vector< std::vector<Tag> > vecAllBlobTags;
	vecAllBlobTags.resize(nTime);

	// Next available patch tag
	Tag tagNextBlob(1, 0);

	// Give blob tags to the initial set of blobs
	std::set<Tag> setAllTags;

	for (int t = 0; t < nTime; t++) {
		vecAllBlobTags[t].resize(vecAllBlobs[t].size());

		tagNextBlob.id = 0;
		for (int p = 0; p < vecAllBlobTags[t].size(); p++) {
			vecAllBlobTags[t][p] = tagNextBlob;
			setAllTags.insert(tagNextBlob);

			tagNextBlob.id++;
		}
		tagNextBlob.time++;
	}

	// Array of equivalent tags
	typedef std::multimap<Tag, Tag> MapGraph;
	typedef MapGraph::const_iterator MapGraphConstIterator;
	typedef MapGraph::iterator MapGraphIterator;

	MapGraph multimapTagGraph;

	// Loop through all remaining time steps
	for (int t = 1; t < nTime; t++) {

		// Get the current blob vector
		const std::vector<Tag> & vecPrevBlobTags = vecAllBlobTags[t-1];

		std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		const std::vector<IndicatorSet> & vecPrevBlobs = vecAllBlobs[t-1];

		const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

		const std::vector<LatLonBox> & vecPrevBlobBoxes = vecAllBlobBoxes[t-1];

		const std::vector<LatLonBox> & vecBlobBoxes = vecAllBlobBoxes[t];

		// Determine overlaps between these blobs and previous blobs
		vecBlobTags.resize(vecBlobs.size());
		for (int p = 0; p < vecBlobTags.size(); p++) {

			const LatLonBox & boxP = vecBlobBoxes[p];

			// Find overlap with bounding boxes at previous time
			for (int q = 0; q < vecPrevBlobTags.size(); q++) {

				const LatLonBox & boxQ = vecPrevBlobBoxes[q];

				// Check if bounding boxes overlap
				if (!boxP.Overlaps(boxQ)) {
					continue;
				}

				// Verify that at least one node overlaps between blobs
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
		//if (mapEquivalentTags.find(*iterTag) != mapEquivalentTags.end()) {
		//	continue;
		//}

		// New blob
		nTotalBlobCount++;
/*
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
*/
		// Set global id
		//if (setBlobTimes.size() >= nMinTime) {
			tagMinimum.global_id = nTotalBlobCount;
		//} else {
		//	tagMinimum.global_id = 0;
		//	nTotalBlobCount--;
		//}

		// Refer all tags in clique to minimum tag
		std::set<Tag>::const_iterator iterTagsVisited = setTagsVisited.begin();
		for (; iterTagsVisited != setTagsVisited.end(); iterTagsVisited++) {
			//mapEquivalentTags.insert(
			//	std::pair<Tag,Tag>(*iterTagsVisited, tagMinimum));
		}
	}
/*
	// Merge blobs at each time step with equivalent tags
	for (int t = 0; t < nTime; t++) {

		std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		for (int p = 0; p < vecBlobTags.size(); p++) {
			std::map<Tag, Tag>::const_iterator iterTagPair =
				mapEquivalentTags.find(vecBlobTags[p]);

			if (iterTagPair != mapEquivalentTags.end()) {
				vecBlobTags[p] = iterTagPair->second;
			}
		}
	}
*/
	Announce("Blobs found: %i", nTotalBlobCount);
/*
	// Apply threshold operators
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

			std::vector<LatLonBox> & vecBlobBoxes = vecAllBlobBoxes[t];

			for (int x = 0; x < vecThresholdOp.size(); x++) {

				std::map<int, bool> mapSatisfiesThreshold;

				vecThresholdOp[x].Apply(
					dCellArea,
					dataLatDeg,
					dataLonDeg,
					vecBlobTags,
					vecBlobs,
					vecBlobBoxes,
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

		NcVar * varLat = ncInput.get_var(latname.c_str());
		NcVar * varLon = ncInput.get_var(lonname.c_str());

		CopyNcVarAttributes(varLat, varOutputLat);
		CopyNcVarAttributes(varLon, varOutputLon);

		NcVar * varTime = ncInput.get_var(timename.c_str());
		if (varTime != NULL) {
			CopyNcVarAttributes(varTime, varOutputTime);
		}
	}

	NcVar * varData =
		ncOutput.add_var(
			strOutputVariable.c_str(),
			ncInt,
			dimOutputTime,
			dimOutputLat,
			dimOutputLon);

	// Loop through all time steps
	DataMatrix<int> dataBlobTag;
	dataBlobTag.Initialize(nLat, nLon);

	int b = 0;
	for (int t = 0; t < nTime; t++) {

		dataBlobTag.Zero();

		// Get the current blob vectors
		const std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

		// Put blob information into matrix
		for (int p = 0; p < vecBlobTags.size(); p++) {

			//if (vecBlobTags[p].global_id == 0) {
			//	continue;
			//}

			b++;
/*
			if (fRejectedBlob[vecBlobTags[p].global_id]) {
				continue;
			}
*/
			IndicatorSetConstIterator iter = vecBlobs[p].begin();
			for (; iter != vecBlobs[p].end(); iter++) {
				dataBlobTag[iter->lat][iter->lon] = b;
					//vecBlobTags[p].global_id;
			}
		}

		// Write to file
		varData->set_cur(t, 0, 0);
		varData->put(&(dataBlobTag[0][0]), 1, nLat, nLon);
	}

	ncOutput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

