///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectCyclones.cpp
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
#include "TimeObj.h"

#include "kdtree.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all minima in the given DataMatrix.
///	</summary>
void FindAllLocalMinima(
	const DataMatrix<float> & data,
	std::set< std::pair<int, int> > & setMinima,
	bool fRegional
) {
	const int nLon = data.GetColumns();
	const int nLat = data.GetRows();

	const int nLonBegin = (fRegional)?(1):(0);
	const int nLonEnd   = (fRegional)?(nLon-1):(nLon);

	// Check interior nodes
	for (int j = 1; j < nLat-1; j++) {
	for (int i = nLonBegin; i < nLonEnd; i++) {

		bool fMinimum = true;
		for (int q = -1; q <= 1; q++) {
		for (int p = -1; p <= 1; p++) {
			int ix = (i + nLon + p) % nLon;
			int jx = (j + q);

			if (data[jx][ix] < data[j][i]) {
				fMinimum = false;
				goto DoneCandidates;
			}
		}
		}

DoneCandidates:
		if (fMinimum) {
			setMinima.insert(std::pair<int,int>(j,i));
		}
	}
	}
/*
	// Check south pole (one node at south pole)
	{
		bool fMinimum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[1][i] < data[0][i]) {
				fMinimum = false;
				break;
			}
		}
		if (fMinimum) {
			_EXCEPTION();
			setMinima.insert(std::pair<int,int>(0,0));
		}
	}

	// Check north pole (one node at north pole)
	{
		bool fMinimum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[nLat-2][i] < data[nLat-1][i]) {
				fMinimum = false;
				break;
			}
		}
		if (fMinimum) {
			_EXCEPTION();
			setMinima.insert(std::pair<int,int>(nLat-1,0));
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all maxima in the given DataMatrix.
///	</summary>
void FindAllLocalMaxima(
	const DataMatrix<float> & data,
	std::set< std::pair<int, int> > & setMaxima,
	bool fRegional
) {
	const int nLon = data.GetColumns();
	const int nLat = data.GetRows();

	const int nLonBegin = (fRegional)?(1):(0);
	const int nLonEnd   = (fRegional)?(nLon-1):(nLon);

	for (int j = 1; j < nLat-1; j++) {
	for (int i = nLonBegin; i < nLonEnd; i++) {

		bool fMaximum = true;
		for (int q = -1; q <= 1; q++) {
		for (int p = -1; p <= 1; p++) {
			int ix = (i + nLon + p) % nLon;
			int jx = (j + q);

			if (data[jx][ix] > data[j][i]) {
				fMaximum = false;
				goto DonePressureMaxima;
			}
		}
		}

DonePressureMaxima:
		if (fMaximum) {
			setMaxima.insert(std::pair<int,int>(j,i));
		}
	}
	}
/*
	// Check south pole (one node at south pole)
	{
		bool fMaximum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[1][i] > data[0][i]) {
				fMaximum = false;
				break;
			}
		}
		if (fMaximum) {
			_EXCEPTION();
			setMaxima.insert(std::pair<int,int>(0,0));
		}
	}

	// Check north pole (one node at north pole)
	{
		bool fMaximum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[nLat-2][i] > data[nLat-1][i]) {
				fMaximum = false;
				break;
			}
		}
		if (fMaximum) {
			_EXCEPTION();
			setMaxima.insert(std::pair<int,int>(nLat-1,0));
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the minimum/maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
void FindLocalMinMax(
	bool fMinimum,
	const DataMatrix<float> & data,
	const DataVector<double> & dataLon,
	const DataVector<double> & dataLat,
	int iLat,
	int iLon,
	double dMaxDist,
	std::pair<int, int> & prMaximum,
	float & dMaxValue,
	float & dRMax,
	bool fRegional
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Number of latitudes and longitudes
	const int nLat = dataLat.GetRows();
	const int nLon = dataLon.GetRows();

	// Initialize the maximum to the central location
	prMaximum = std::pair<int, int>(iLat, iLon);
	dMaxValue = data[iLat][iLon];
	dRMax = 0.0;

	// Queue of nodes that remain to be visited
	std::queue< std::pair<int, int> > queueNodes;
	queueNodes.push(prMaximum);

	// Set of nodes that have already been visited
	std::set< std::pair<int, int> > setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = dataLat[iLat];
	double dLon0 = dataLon[iLon];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		std::pair<int, int> pr = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(pr) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(pr);

		double dLatThis = dataLat[pr.first];
		double dLonThis = dataLon[pr.second];

		// Great circle distance to this element
		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local maximum
		if (fMinimum) {
			if (data[pr.first][pr.second] < dMaxValue) {
				prMaximum = pr;
				dMaxValue = data[pr.first][pr.second];
				dRMax = dR;
			}

		} else {
			if (data[pr.first][pr.second] > dMaxValue) {
				prMaximum = pr;
				dMaxValue = data[pr.first][pr.second];
				dRMax = dR;
			}
		}

		// Add all neighbors of this point
		std::pair<int,int> prWest(pr.first, (pr.second + nLon - 1) % nLon);
		if (setNodesVisited.find(prWest) == setNodesVisited.end()) {
			if (!fRegional) {
				queueNodes.push(prWest);
			} else if (pr.second >= 0) {
				queueNodes.push(prWest);
			}
		}

		std::pair<int,int> prEast(pr.first, (pr.second + 1) % nLon);
		if (setNodesVisited.find(prEast) == setNodesVisited.end()) {
			if (!fRegional) {
				queueNodes.push(prEast);
			} else if (pr.second <= nLon-1) {
				queueNodes.push(prEast);
			}
		}

		std::pair<int,int> prNorth(pr.first + 1, pr.second);
		if ((prNorth.first < nLat) &&
			(setNodesVisited.find(prNorth) == setNodesVisited.end())
		) {
			queueNodes.push(prNorth);
		}

		std::pair<int,int> prSouth(pr.first - 1, pr.second);
		if ((prSouth.first >= 0) &&
			(setNodesVisited.find(prSouth) == setNodesVisited.end())
		) {
			queueNodes.push(prSouth);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse a pair of Date values.
///	</summary>
void ParseDate(
	int nDate,
	int nDateSec,
	int & nDateYear,
	int & nDateMonth,
	int & nDateDay,
	int & nDateHour
) {
	nDateYear  = nDate / 10000;
	nDateMonth = (nDate % 10000) / 100;
	nDateDay   = (nDate % 100);
	nDateHour  = nDateSec / 3600;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse a time value.
///	</summary>
void ParseTimeDouble(
	const std::string & strTimeUnits,
	const std::string & strTimeCalendar,
	double dTime,
	int & nDateYear,
	int & nDateMonth,
	int & nDateDay,
	int & nDateHour
) {
	// Get calendar type
	Time::CalendarType cal;
	if ((strTimeCalendar.length() >= 6) &&
		(strncmp(strTimeCalendar.c_str(), "noleap", 6) == 0)
	) {
		cal = Time::CalendarNoLeap;

	} else if (
		(strTimeCalendar.length() >= 8) &&
		(strncmp(strTimeCalendar.c_str(), "standard", 8) == 0)
	) {
		cal = Time::CalendarStandard;

	} else {
		_EXCEPTION1("Unknown calendar type \"%s\"", strTimeCalendar.c_str());
	}

	// Time format is "days since ..."
	if ((strTimeUnits.length() >= 11) &&
	    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(11);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		int nDays = static_cast<int>(dTime);
		time.AddDays(nDays);

		int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 86400.0);
		time.AddSeconds(nSeconds);

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	// Time format is "hours since ..."
	} else if (
	    (strTimeUnits.length() >= 12) &&
	    (strncmp(strTimeUnits.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(12);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 3600.0);
		time.AddSeconds(nSeconds);

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
	//_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field has a closed contour about this point.
///	</summary>
bool HasClosedContour(
	const DataVector<double> & dataLat,
	const DataVector<double> & dataLon,
	const DataMatrix<float> & dataState,
	const int iLat0,
	const int iLon0,
	double dDeltaAmt,
	double dDeltaDist,
	int nDeltaCount,
	double dMinMaxDist,
	bool fRegional
) {
	// Verify arguments
	if (dDeltaAmt == 0.0) {
		_EXCEPTIONT("Closed contour amount must be non-zero");
	}
	if (dDeltaDist <= 0.0) {
		_EXCEPTIONT("Closed contour distance must be positive");
	}
	if (nDeltaCount <= 0) {
		_EXCEPTIONT("Number of sample points for closed contour "
			"must be positive");
	}

	// Grid spacing
	const int nLat = dataLat.GetRows();
	const int nLon = dataLon.GetRows();

	const double dDeltaLat = (dataLat[1] - dataLat[0]);
	const double dDeltaLon = (dataLon[1] - dataLon[0]);

	// Find min/max near point
	int iLat;
	int iLon;

	if (dMinMaxDist == 0.0) {
		iLat = iLat0;
		iLon = iLon0;

	// Find a local minimum / maximum
	} else {
		std::pair<int, int> prExtrema;
		float dValue;
		float dR;

		FindLocalMinMax(
			(dDeltaAmt > 0.0),
			dataState,
			dataLon,
			dataLat,
			iLat0,
			iLon0,
			dMinMaxDist,
			prExtrema,
			dValue,
			dR,
			fRegional);

		iLat = prExtrema.first;
		iLon = prExtrema.second;
	}

	// Coordinates of detection point	
	double dLat = dataLat[iLat];
	double dLon = dataLon[iLon];

	// This location in Cartesian geometry
	double dX0 = cos(dLon) * cos(dLat);
	double dY0 = sin(dLon) * cos(dLat);
	double dZ0 = sin(dLat);

	// Pick a quasi-arbitrary reference direction
	double dX1;
	double dY1;
	double dZ1;
	if ((fabs(dX0) >= fabs(dY0)) && (fabs(dX0) >= fabs(dZ0))) {
		dX1 = dX0;
		dY1 = dY0 + 1.0;
		dZ1 = dZ0;
	} else if ((fabs(dY0) >= fabs(dX0)) && (fabs(dY0) >= fabs(dZ0))) {
		dX1 = dX0;
		dY1 = dY0;
		dZ1 = dZ0 + 1.0;
	} else {
		dX1 = dX0 + 1.0;
		dY1 = dY0;
		dZ1 = dZ0;
	}

	// Project perpendicular to detection location
	double dDot = dX1 * dX0 + dY1 * dY0 + dZ1 * dZ0;

	dX1 -= dDot * dX0;
	dY1 -= dDot * dY0;
	dZ1 -= dDot * dZ0;

	// Normalize
	double dMag1 = sqrt(dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1);

	if (dMag1 < 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	double dScale1 = tan(dDeltaDist * M_PI / 180.0) / dMag1;

	dX1 *= dScale1;
	dY1 *= dScale1;
	dZ1 *= dScale1;

	// Verify dot product is zero
	dDot = dX0 * dX1 + dY0 * dY1 + dZ0 * dZ1;
	if (fabs(dDot) > 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	// Cross product (magnitude automatically 
	double dCrossX = dY0 * dZ1 - dZ0 * dY1;
	double dCrossY = dZ0 * dX1 - dX0 * dZ1;
	double dCrossZ = dX0 * dY1 - dY0 * dX1;

	// Obtain reference value
	float dRefValue = dataState[iLat][iLon];

	// Output closed contour details
	AnnounceStartBlock(2, "Closed contour:");
	Announce(2, "Storm at (%i, %i) : (%1.5f, %1.5f)",
		iLat0, iLon0, dataLat[iLat0] * 180.0 / M_PI, dataLon[iLon0] * 180.0 / M_PI);
	Announce(2, "Extrema at (%i, %i) : (%1.5f, %1.5f) : %1.5e",
		iLat, iLon, dLat * 180.0 / M_PI, dLon * 180.0 / M_PI, dRefValue);

	// Loop through all sample points
	for (int a = 0; a < nDeltaCount; a++) {

		// Angle of rotation
		double dAngle = 2.0 * M_PI
			* static_cast<double>(a)
			/ static_cast<double>(nDeltaCount);

		// Calculate new rotated vector
		double dX2 = dX0 + dX1 * cos(dAngle) + dCrossX * sin(dAngle);
		double dY2 = dY0 + dY1 * cos(dAngle) + dCrossY * sin(dAngle);
		double dZ2 = dZ0 + dZ1 * cos(dAngle) + dCrossZ * sin(dAngle);

		double dMag2 = sqrt(dX2 * dX2 + dY2 * dY2 + dZ2 * dZ2);

		dX2 /= dMag2;
		dY2 /= dMag2;
		dZ2 /= dMag2;

		// Calculate new lat/lon location
		double dLat2 = asin(dZ2);
		double dLon2 = atan2(dY2, dX2);

		if (dLon2 < dataLon[0]) {
			dLon2 += 2.0 * M_PI;
		}
/*
		printf("%i : %3.2f %3.2f : %3.2f %3.2f\n",
			a,
			dLat * 180.0 / M_PI,
			dLon * 180.0 / M_PI,
			dLat2 * 180.0 / M_PI,
			dLon2 * 180.0 / M_PI);
*/
		int j = static_cast<int>((dLat2 - dataLat[0]) / dDeltaLat + 0.5);
		int i = static_cast<int>((dLon2 - dataLon[0]) / dDeltaLon + 0.5);

		if (!fRegional) {
			if (i == nLon) {
				i = nLon-1;
			}
			if (j == nLat) {
				j = nLat-1;
			}
		}

		//printf("%i %i : %i %i\n", iLat, iLon, j, i);

		// Check for insufficient distance
		if ((i == iLon) && (j == iLat)) {
			_EXCEPTIONT("Closed contour distance insufficient; increase value");
		}

		// Check for out of bounds
		if ((i < 0) || (j < 0) || (i >= nLon) || (j >= nLat)) {
			if (fRegional) {
				Announce(2, "Point (%i, %i) out of bounds, returning", j, i);
				AnnounceEndBlock(2, NULL);
				return false;
			}
			_EXCEPTION4("Logic error %i/%i, %i/%i", i, nLon, j, nLat);
		}

		Announce(2, "Checking point (%i, %i) : (%1.5f %1.5f) : %1.5e",
			j, i, dLat2 * 180.0 / M_PI, dLon2 * 180.0 / M_PI, dataState[j][i]);

		// Verify sufficient increase in value
		if (dDeltaAmt > 0.0) {
			if (dataState[j][i] - dRefValue < dDeltaAmt) {
				Announce(2, "Failed criteria; returning");
				AnnounceEndBlock(2, NULL);
				return false;
			}

		// Verify sufficient decrease in value
		} else {
			if (dRefValue - dataState[j][i] < -dDeltaAmt) {
				Announce(2, "Failed criteria; returning");
				AnnounceEndBlock(2, NULL);
				return false;
			}
		}
	}

	Announce(2, "Passed criteria; returning");
	AnnounceEndBlock(2, NULL);
	return true;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class describing a general closed contour operation.
///	</summary>
class ClosedContourOp {

public:
	///	<summary>
	///		Parse a closed contour operation string.
	///	</summary>
	void Parse(
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Variable,
			ReadMode_Distance,
			ReadMode_Amount,
			ReadMode_Count,
			ReadMode_MinMaxDist,
			ReadMode_Invalid
		} eReadMode = ReadMode_Variable;

		// Loop through string
		int iLast = 0;
		for (int i = 0; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in variable name
				if (eReadMode == ReadMode_Variable) {
					m_strVariable = strSubStr;

					iLast = i + 1;
					eReadMode = ReadMode_Amount;

				// Read in amount
				} else if (eReadMode == ReadMode_Amount) {
					m_dDeltaAmount = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in distance
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Count;

				// Read in count
				} else if (eReadMode == ReadMode_Count) {
					m_nCount = atoi(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_MinMaxDist;

				// Read in min/max distance
				} else if (eReadMode == ReadMode_MinMaxDist) {
					m_dMinMaxDist = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nToo many entries in closed contour op \"%s\""
						"\nRequired: \"<name>,<amount>,<distance>,"
						"<count>,<minmaxdist>\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in closed contour op \"%s\""
					"\nRequired: \"<name>,<amount>,<distance>,"
					"<count>,<minmaxdist>\"",
					strOp.c_str());
		}

		// Output announcement
		char szBuffer[128];

		std::string strDescription;
		strDescription += m_strVariable;

		if (m_dDeltaAmount == 0.0) {
			_EXCEPTIONT("For closed contour op, delta amount must be non-zero");
		}
		if (m_dDistance <= 0.0) {
			_EXCEPTIONT("For closed contour op, distance must be positive");
		}
		if (m_nCount <= 0) {
			_EXCEPTIONT("For closed contour op, count must be positive");
		}
		if (m_dMinMaxDist < 0.0) {
			_EXCEPTIONT("For closed contour op, min/max dist must be nonnegative");
		}

		if (m_dDeltaAmount < 0.0) {
			Announce("%s decreases by %f over %f deg (%i samples)"
				   " (max search %f deg)",
				m_strVariable.c_str(),
				-m_dDeltaAmount,
				m_dDistance,
				m_nCount,
				m_dMinMaxDist);

		} else {
			Announce("%s increases by %f over %f degrees (%i samples)"
					" (min search %f deg)",
				m_strVariable.c_str(),
				m_dDeltaAmount,
				m_dDistance,
				m_nCount,
				m_dMinMaxDist);
		}
	}

public:
	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strVariable;

	///	<summary>
	///		Threshold amount.  If positive this represents a minimum
	///		increase.  If negative this represents a minimum decrease.
	///	</summary>
	double m_dDeltaAmount;

	///	<summary>
	///		Threshold distance.
	///	</summary>
	double m_dDistance;

	///	<summary>
	///		Number of nodes used to evaluate closed contour.
	///	</summary>
	int m_nCount;

	///	<summary>
	///		Distance to search for min or max.
	///	</summary>
	double m_dMinMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {
	// Gravitational constant
	const float ParamGravity = 9.80616;

	// Radius of the Earth
	const float ParamEarthRadius = 6.37122e6;

	// Input file
	std::string strInputFile;

	// Output file
	std::string strOutputFile;

	// Variable to search for the minimum
	std::string strSearchByMin;

	// Variable to search for the maximum
	std::string strSearchByMax;

	// Zonal velocity name
	std::string strUName;

	// Meridional velocity name
	std::string strVName;

	// Data is regional
	bool fRegional;

	// Maximum latitude for detection
	double dMaxLatitude;

	// Minimum latitude for detection
	double dMinLatitude;

	// File containing information on topography
	std::string strTopoFile;

	// Maximum topographic height for a detection
	double dMaxTopoHeight;

	// Merge distance
	double dMergeDist;

	// Require temperature maxima at T200 and T500 within this distance
	double dWarmCoreDist;

	// No temperature maxima at T200 and T500 within this distance
	double dNoWarmCoreDist;

	// Require 850hPa vorticity maxima within this distance
	double dVortDist;

	// Closed contour commands
	std::string strClosedContourCmd;

	// Minimum Laplacian value
	double dMinLaplacian;

	// Distance to search for maximum wind speed
	double dWindSpDist;

	// Append to output file
	bool fAppend;

	// Output header
	bool fOutputHeader;

	// Verbosity level
	int iVerbosityLevel;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineStringD(strSearchByMin, "searchbymin", "", "(default PSL)");
		CommandLineString(strSearchByMax, "searchbymax", "");
		CommandLineString(strUName, "uname", "U850");
		CommandLineString(strVName, "vname", "V850");
		CommandLineBool(fRegional, "regional");
		CommandLineDoubleD(dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineString(strTopoFile, "topofile", "");
		CommandLineDoubleD(dMaxTopoHeight, "maxtopoht", 0.0, "(m)");
		CommandLineDoubleD(dMergeDist, "mergedist", 0.0, "(degrees)");
		CommandLineDoubleD(dWarmCoreDist, "warmcoredist", 0.0, "(degrees)");
		CommandLineDoubleD(dNoWarmCoreDist, "nowarmcoredist", 0.0, "(degrees)");
		CommandLineDoubleD(dVortDist, "vortdist", 0.0, "(degrees)");
		CommandLineString(strClosedContourCmd, "closedcontourcmd", "");
		//CommandLineDoubleD(dMinLaplacian, "minlaplacian", 0.0, "(Pa / degree^2)");
		CommandLineDoubleD(dWindSpDist, "windspdist", 0.0, "(degrees)");
		//CommandLineBool(fAppend, "append");
		CommandLineBool(fOutputHeader, "out_header");
		CommandLineInt(iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Set verbosity level
	AnnounceSetVerbosityLevel(iVerbosityLevel);

	// Check input
	if (strInputFile.length() == 0) {
		_EXCEPTIONT("No input file (--in) specified");
	}

	// Check output
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check warm core distance
	if ((dWarmCoreDist != 0.0) && (dNoWarmCoreDist != 0.0)) {
		_EXCEPTIONT("Only one of --warmcoredist and --nowarmcoredist"
			   " may be active");
	}

	// Only one of search by min or search by max should be specified
	if ((strSearchByMin == "") && (strSearchByMax == "")) {
		strSearchByMin = "PSL";
	}
	if ((strSearchByMin != "") && (strSearchByMax != "")) {
		_EXCEPTIONT("Only one of --searchbymin or --searchbymax can"
			" be specified");
	}

	// Parse the closed contour command string
	std::vector<ClosedContourOp> vecClosedContourOp;

	if (strClosedContourCmd != "") {
		AnnounceStartBlock("Parsing closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strClosedContourCmd.length(); i++) {

			if ((i == strClosedContourCmd.length()) ||
				(strClosedContourCmd[i] == ';')
			) {
				std::string strSubStr =
					strClosedContourCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecClosedContourOp.size());
				vecClosedContourOp.resize(iNextOp + 1);
				vecClosedContourOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Check minimum latitude and maximum latitude
	if (dMaxLatitude < 0.0) {
		_EXCEPTIONT("--maxlat must be nonnegative");
	}
	if (dMinLatitude < 0.0) {
		_EXCEPTIONT("--minlat must be nonnegative");
	}

	dMaxLatitude *= M_PI / 180.0;
	dMinLatitude *= M_PI / 180.0;

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

	const int nLat = dimLat->size();
	const int nLon = dimLon->size();

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

	// Get auxiliary variables
	bool fSearchByMinima = false;
	NcVar * varSearch;
	if (strSearchByMin != "") {
		fSearchByMinima = true;
		varSearch = ncInput.get_var(strSearchByMin.c_str());
		if (varSearch == NULL) {
			_EXCEPTION1("Unable to find variable \"%s\"",
				strSearchByMin.c_str());
		}

	} else {
		fSearchByMinima = false;
		varSearch = ncInput.get_var(strSearchByMax.c_str());
		if (varSearch == NULL) {
			_EXCEPTION1("Unable to find variable \"%s\"",
				strSearchByMax.c_str());
		}
	}

	NcVar * varPSL = ncInput.get_var("PSL");
	NcVar * varUvel = ncInput.get_var(strUName.c_str());
	NcVar * varVvel = ncInput.get_var(strVName.c_str());
	NcVar * varT200 = ncInput.get_var("T200");
	NcVar * varT500 = ncInput.get_var("T500");

	if (varUvel == NULL) {
		_EXCEPTION1("No variable \"%s\" found in input file", strUName.c_str());
	}
	if (varVvel == NULL) {
		_EXCEPTION1("No variable \"%s\" found in input file", strVName.c_str());
	}
	if (varT200 == NULL) {
		_EXCEPTIONT("No variable \"T200\" found in input file");
	}
	if (varT500 == NULL) {
		_EXCEPTIONT("No variable \"T500\" found in input file");
	}

	// Storage for auxiliary variables
	DataMatrix<float> dataSearch(nLat, nLon);
	DataMatrix<float> dataPSL(nLat, nLon);
	DataMatrix<float> dataUvel(nLat, nLon);
	DataMatrix<float> dataVvel(nLat, nLon);
	DataMatrix<float> dataT200(nLat, nLon);
	DataMatrix<float> dataT500(nLat, nLon);

	DataMatrix<float> dataUMag850(nLat, nLon);

	DataMatrix<float> dataDel2PSL(nLat, nLon);

	DataMatrix<float> dataPHIS(nLat, nLon);

	// Topography variable
	NcVar * varPHIS = NULL;

	if (strTopoFile != "") {
		NcFile ncTopo(strTopoFile.c_str());
		NcVar * varPHIS = ncTopo.get_var("PHIS");

		if (varPHIS == NULL) {
			_EXCEPTIONT("No variable \"PHIS\" found in topography file");
		}

		if (varPHIS->num_dims() == 3) {
			varPHIS->get(&(dataPHIS[0][0]), 1, nLat, nLon);
		} else if (varPHIS->num_dims() == 2) {
			varPHIS->get(&(dataPHIS[0][0]), nLat, nLon);
		} else {
			_EXCEPTIONT("Invalid number of dimensions on variable \"PHIS\"");
		}
	}

	if ((strTopoFile == "") && (dMaxTopoHeight != 0.0)) {
		_EXCEPTIONT("No topography file specified; required for --maxtopoht");
	}

	// Open output file
	FILE * fpOutput = fopen(strOutputFile.c_str(), "w");
	if (fpOutput == NULL) {
		_EXCEPTION1("Could not open output file \"%s\"",
			strOutputFile.c_str());
	}

	if (fOutputHeader) {
		fprintf(fpOutput, "#day\tmonth\tyear\tcount\thour\n");
		fprintf(fpOutput, "#\t#\ti\tj\tpsl_lon\tpsl_lat\twind_max\tr_wind_max\tpsl_min\n");
	}

	// Loop through all times
	for (int t = 0; t < nTime; t++) {
	//for (int t = 0; t < 1; t++) {

		char szStartBlock[128];
		sprintf(szStartBlock, "Time %i", t);
		AnnounceStartBlock(szStartBlock);

		// Get the auxiliary variables
		varSearch->set_cur(t,0,0);
		varSearch->get(&(dataSearch[0][0]), 1, nLat, nLon);

		varPSL->set_cur(t,0,0);
		varPSL->get(&(dataPSL[0][0]), 1, nLat, nLon);

		varUvel->set_cur(t,0,0);
		varUvel->get(&(dataUvel[0][0]), 1, nLat, nLon);

		varVvel->set_cur(t,0,0);
		varVvel->get(&(dataVvel[0][0]), 1, nLat, nLon);

		varT200->set_cur(t,0,0);
		varT200->get(&(dataT200[0][0]), 1, nLat, nLon);

		varT500->set_cur(t,0,0);
		varT500->get(&(dataT500[0][0]), 1, nLat, nLon);

		// Compute wind magnitude
		for (int j = 0; j < nLat; j++) {
		for (int i = 0; i < nLon; i++) {
			dataUMag850[j][i] = sqrt(
				  dataUvel[j][i] * dataUvel[j][i]
				+ dataVvel[j][i] * dataVvel[j][i]);
		}
		}

		// Tag all pressure minima
		std::set< std::pair<int, int> > setCandidates;
		if (strSearchByMin != "") {
			FindAllLocalMinima(dataSearch, setCandidates, fRegional);
		} else {
			FindAllLocalMaxima(dataSearch, setCandidates, fRegional);
		}

		// Total number of pressure minima
		int nTotalCandidates = setCandidates.size();
		int nRejectedLatitude = 0;
		int nRejectedTopography = 0;
		int nRejectedMerge = 0;
		int nRejectedVortMax = 0;
		int nRejectedWarmCore = 0;
		int nRejectedNoWarmCore = 0;
		int nRejectedLaplacian = 0;

		DataVector<int> vecRejectedClosedContour(vecClosedContourOp.size());

		// Eliminate based on maximum latitude
		if (dMaxLatitude != 0.0) {
			std::set< std::pair<int, int> > setNewCandidates;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = dataLat[iterCandidate->first];

				if (fabs(dLat) <= dMaxLatitude) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedLatitude++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on minimum latitude
		if (dMinLatitude != 0.0) {
			std::set< std::pair<int, int> > setNewCandidates;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = dataLat[iterCandidate->first];

				if (fabs(dLat) >= dMinLatitude) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedLatitude++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on maximum topographic height
		if (dMaxTopoHeight != 0.0) {
			std::set< std::pair<int, int> > setNewCandidates;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				float dPHIS = dataPHIS[iterCandidate->first][iterCandidate->second];

				double dTopoHeight =
					static_cast<double>(dPHIS / ParamGravity);

				if (dTopoHeight <= dMaxTopoHeight) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedTopography++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on merge distance
		if (dMergeDist != 0.0) {
			std::set< std::pair<int, int> > setNewCandidates;

			// Calculate spherical distance
			double dSphDist = 2.0 * sin(0.5 * dMergeDist / 180.0 * M_PI);

			// Create a new KD Tree containing all nodes
			kdtree * kdMerge = kd_create(3);

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = dataLat[iterCandidate->first];
				double dLon = dataLon[iterCandidate->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdMerge, dX, dY, dZ, (void*)(&(*iterCandidate)));
			}

			// Loop through all PSL find set of nearest neighbors
			iterCandidate = setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = dataLat[iterCandidate->first];
				double dLon = dataLon[iterCandidate->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				// Find all neighbors within dSphDist
				kdres * kdresMerge =
					kd_nearest_range3(kdMerge, dX, dY, dZ, dSphDist);

				// Number of neighbors
				int nNeighbors = kd_res_size(kdresMerge);
				if (nNeighbors == 0) {
					setNewCandidates.insert(*iterCandidate);

				} else {
					double dValue =
						dataSearch[iterCandidate->first][iterCandidate->second];

					bool fExtrema = true;
					for (;;) {
						std::pair<int,int> * ppr =
							(std::pair<int,int> *)(kd_res_item_data(kdresMerge));

						if (fSearchByMinima) {
							if (dataSearch[ppr->first][ppr->second] < dValue) {
								fExtrema = false;
								break;
							}

						} else {
							if (dataSearch[ppr->first][ppr->second] > dValue) {
								fExtrema = false;
								break;
							}
						}

						int iHasMore = kd_res_next(kdresMerge);
						if (!iHasMore) {
							break;
						}
					}

					if (fExtrema) {
						setNewCandidates.insert(*iterCandidate);
					} else {
						nRejectedMerge++;
					}
				}

				kd_res_free(kdresMerge);
			}

			// Destroy the KD Tree
			kd_free(kdMerge);

			// Update set of pressure minima
			setCandidates = setNewCandidates;
		}

		// Eliminate based on presence of vorticity maximum
		if (dVortDist != 0.0) {

			float dDeltaLat = static_cast<float>(dataLat[1] - dataLat[0]);
			float dDeltaLon = static_cast<float>(dataLon[1] - dataLon[0]);

			int iLonBegin = (fRegional)?(1):(0);
			int iLonEnd   = (fRegional)?(nLon-1):(nLon);

			DataMatrix<float> dataZETA850(nLat, nLon);

			// Compute vorticity
			for (int j = 1; j < nLat-1; j++) {
			for (int i = iLonBegin; i < iLonEnd; i++) {

				int inext = (i + 1) % nLon;
				int jnext = (j + 1);

				int ilast = (i + nLon - 1) % nLon;
				int jlast = (j - 1);

				float dDyUvel = (dataUvel[jnext][i] - dataUvel[jlast][i])
					/ (2.0 * dDeltaLat);
				float dDxVvel = (dataVvel[j][inext] - dataVvel[j][ilast])
					/ (2.0 * dDeltaLon);

				dataZETA850[j][i] = dDxVvel - dDyUvel
					+ dataUvel[j][i] / ParamEarthRadius * tan(dataLat[j]);

				if (dataLat[j] < 0.0) {
					dataZETA850[j][i] *= -1.0;
				}
			}
			}

			// Find all vorticity maxima
			std::set< std::pair<int, int> > setZETA850Maxima;
			FindAllLocalMaxima(dataZETA850, setZETA850Maxima, fRegional);

			// Construct KD tree for ZETA850
			kdtree * kdZETA850Maxima = kd_create(3);
			std::set< std::pair<int, int> >::const_iterator iterZETA850
				= setZETA850Maxima.begin();

			for (; iterZETA850 != setZETA850Maxima.end(); iterZETA850++) {
				double dLat = dataLat[iterZETA850->first];
				double dLon = dataLon[iterZETA850->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdZETA850Maxima, dX, dY, dZ, NULL);
			}

			// Remove pressure minima that are near temperature maxima
			std::set< std::pair<int, int> > setNewCandidates;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = dataLat[iterCandidate->first];
				double dLon = dataLon[iterCandidate->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kdres * kdresZETA850 = kd_nearest3(kdZETA850Maxima, dX, dY, dZ);

				double dZETA850pos[3];

				kd_res_item(kdresZETA850, dZETA850pos);

				kd_res_free(kdresZETA850);

				double dZETA850dist = sqrt(
					  (dX - dZETA850pos[0]) * (dX - dZETA850pos[0])
					+ (dY - dZETA850pos[1]) * (dY - dZETA850pos[1])
					+ (dZ - dZETA850pos[2]) * (dZ - dZETA850pos[2]));

				dZETA850dist = 2.0 * asin(0.5 * dZETA850dist) * 180.0 / M_PI;

				// Reject storms with no warm core
				if (dZETA850dist <= dVortDist) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedVortMax++;
				}
			}

			kd_free(kdZETA850Maxima);

			setCandidates = setNewCandidates;
		}

		// Detect presence of warm core near PSL min
		if ((dWarmCoreDist != 0.0) || (dNoWarmCoreDist != 0.0)) {

			std::set< std::pair<int, int> > setT200Maxima;
			FindAllLocalMaxima(dataT200, setT200Maxima, fRegional);

			std::set< std::pair<int, int> > setT500Maxima;
			FindAllLocalMaxima(dataT500, setT500Maxima, fRegional);

			// Construct KD tree for T200
			kdtree * kdT200Maxima = kd_create(3);
			std::set< std::pair<int, int> >::const_iterator iterT200
				= setT200Maxima.begin();

			for (; iterT200 != setT200Maxima.end(); iterT200++) {
				double dLat = dataLat[iterT200->first];
				double dLon = dataLon[iterT200->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdT200Maxima, dX, dY, dZ, NULL);
			}

			// Construct KD tree for T500
			kdtree * kdT500Maxima = kd_create(3);
			std::set< std::pair<int, int> >::const_iterator iterT500
				= setT500Maxima.begin();

			for (; iterT500 != setT500Maxima.end(); iterT500++) {
				double dLat = dataLat[iterT500->first];
				double dLon = dataLon[iterT500->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdT500Maxima, dX, dY, dZ, NULL);
			}

			// Remove pressure minima that are near temperature maxima
			std::set< std::pair<int, int> > setNewCandidates;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = dataLat[iterCandidate->first];
				double dLon = dataLon[iterCandidate->second];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kdres * kdresT200 = kd_nearest3(kdT200Maxima, dX, dY, dZ);
				kdres * kdresT500 = kd_nearest3(kdT500Maxima, dX, dY, dZ);

				double dT200pos[3];
				double dT500pos[3];

				kd_res_item(kdresT200, dT200pos);
				kd_res_item(kdresT500, dT500pos);

				kd_res_free(kdresT200);
				kd_res_free(kdresT500);

				double dT200dist = sqrt(
					  (dX - dT200pos[0]) * (dX - dT200pos[0])
					+ (dY - dT200pos[1]) * (dY - dT200pos[1])
					+ (dZ - dT200pos[2]) * (dZ - dT200pos[2]));

				dT200dist = 2.0 * asin(0.5 * dT200dist) * 180.0 / M_PI;

				double dT500dist = sqrt(
					  (dX - dT500pos[0]) * (dX - dT500pos[0])
					+ (dY - dT500pos[1]) * (dY - dT500pos[1])
					+ (dZ - dT500pos[2]) * (dZ - dT500pos[2]));

				dT500dist = 2.0 * asin(0.5 * dT500dist) * 180.0 / M_PI;

				// Reject storms with warm core
				if (dNoWarmCoreDist != 0.0) {
					if ((dT200dist >= dNoWarmCoreDist) ||
						(dT500dist >= dNoWarmCoreDist)
					) {
						setNewCandidates.insert(*iterCandidate);
					} else {
						nRejectedWarmCore++;
					}
				}

				// Reject storms with no warm core
				if (dWarmCoreDist != 0.0) {
					if ((dT200dist <= dWarmCoreDist) &&
						(dT500dist <= dWarmCoreDist)
					) {
						setNewCandidates.insert(*iterCandidate);
					} else {
						nRejectedNoWarmCore++;
					}
				}
			}

			kd_free(kdT200Maxima);
			kd_free(kdT500Maxima);

			setCandidates = setNewCandidates;
		}

		// Eliminate based on closed contours
		for (int ccc = 0; ccc < vecClosedContourOp.size(); ccc++) {
			std::set< std::pair<int, int> > setNewCandidates;

			// Load relevant data
			DataMatrix<float> dataState(nLat, nLon);
			
			NcVar * varDataState =
				ncInput.get_var(vecClosedContourOp[ccc].m_strVariable.c_str());

			varDataState->set_cur(t,0,0);
			varDataState->get(&(dataState[0][0]), 1, nLat, nLon);

			// Loop through all pressure minima
			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if pressure minima have a closed contour
				bool fHasClosedContour =
					HasClosedContour(
						dataLat,
						dataLon,
						dataState,
						iterCandidate->first,
						iterCandidate->second,
						vecClosedContourOp[ccc].m_dDeltaAmount,
						vecClosedContourOp[ccc].m_dDistance,
						vecClosedContourOp[ccc].m_nCount,
						vecClosedContourOp[ccc].m_dMinMaxDist,
						fRegional
					);

				// If not rejected, add to new pressure minima array
				if (fHasClosedContour) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					vecRejectedClosedContour[ccc]++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on minimum the Laplacian of pressure at the PSL min
		if (dMinLaplacian != 0.0) {
			float dDeltaLat = static_cast<float>(dataLat[1] - dataLat[0]);
			float dDeltaLon = static_cast<float>(dataLon[1] - dataLon[0]);

			std::set< std::pair<int, int> > setNewCandidates;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				int i = iterCandidate->second;
				int j = iterCandidate->first;

				int inext = (i + 1) % nLon;
				int jnext = (j + 1);

				int ilast = (i + nLon - 1) % nLon;
				int jlast = (j - 1);

				float dDphiPSL =
					(dataPSL[jnext][i] - dataPSL[jlast][i]) / (2.0 * dDeltaLat);
				float dDlambdaPSL =
					(dataPSL[j][inext] - dataPSL[j][ilast]) / (2.0 * dDeltaLon);
				float dD2phiPSL =
					(dataPSL[jnext][i] - 2.0 * dataPSL[j][i] + dataPSL[jlast][i])
					/ (dDeltaLat * dDeltaLat);
				float dD2lambdaPSL =
					(dataPSL[j][inext] - 2.0 * dataPSL[j][i] + dataPSL[j][ilast])
					/ (dDeltaLon * dDeltaLon);

				float dSecLat = 1.0 / cos(dataLat[j]);
			
				float dLaplacian =
					dD2phiPSL - tan(dataLat[j]) * dDphiPSL
					+ dSecLat * dSecLat * dD2lambdaPSL;

				// Convert to Pa / degree^2
				dLaplacian *= (M_PI / 180.0) * (M_PI / 180.0);

				if (dLaplacian >= dMinLaplacian) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedLaplacian++;
				}
			}

			setCandidates = setNewCandidates;
		}

		Announce("Total candidates: %i", setCandidates.size());
		Announce("Rejected (    latitude): %i", nRejectedLatitude);
		Announce("Rejected (  topography): %i", nRejectedTopography);
		Announce("Rejected (      merged): %i", nRejectedMerge);
		Announce("Rejected ( no vort max): %i", nRejectedVortMax);
		Announce("Rejected (   warm core): %i", nRejectedWarmCore);
		Announce("Rejected (no warm core): %i", nRejectedNoWarmCore);
		//Announce("Rejected ( slp contour): %i", nRejectedSLPClosedContour);
		//Announce("Rejected (temp contour): %i", nRejectedTempClosedContour);

		for (int ccc = 0; ccc < vecRejectedClosedContour.GetRows(); ccc++) {
			Announce("Rejected: (contour %s): %i",
					vecClosedContourOp[ccc].m_strVariable.c_str(),
					vecRejectedClosedContour[ccc]);
		}

		// Determine wind maximum at all pressure minima
		AnnounceStartBlock("Searching for maximum winds");
		std::vector<float> vecMaxWindSp;
		std::vector<float> vecRMaxWindSp;
		{
			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				std::pair<int, int> prMaximum;

				float dRMaxWind;
				float dMaxWindSp;

				FindLocalMinMax(
					false,
					dataUMag850,
					dataLon,
					dataLat,
					iterCandidate->first,
					iterCandidate->second,
					dWindSpDist,
					prMaximum,
					dMaxWindSp,
					dRMaxWind,
					fRegional);

				vecMaxWindSp.push_back(dMaxWindSp);
				vecRMaxWindSp.push_back(dRMaxWind);
			}
		}

		AnnounceEndBlock("Done");

		// Write results to file
		{
			// Parse time information
			//NcVar * varDate = ncInput.get_var("date");
			//NcVar * varDateSec = ncInput.get_var("datesec");

			int nDateYear;
			int nDateMonth;
			int nDateDay;
			int nDateHour;
/*
			if ((varDate != NULL) && (varDateSec != NULL)) {
				int nDate;
				int nDateSec;

				varDate->set_cur(t);
				varDate->get(&nDate, 1);

				varDateSec->set_cur(t);
				varDateSec->get(&nDateSec, 1);

				ParseDate(
					nDate,
					nDateSec,
					nDateYear,
					nDateMonth,
					nDateDay,
					nDateHour);

			} else {
*/
				NcAtt * attTimeUnits = varTime->get_att("units");
				if (attTimeUnits == NULL) {
					_EXCEPTIONT("Variable \"time\" has no \"units\" attribute");
				}

				std::string strTimeUnits = attTimeUnits->as_string(0);


				std::string strTimeCalendar = "noleap";
				NcAtt * attTimeCalendar = varTime->get_att("calendar");
				if (attTimeUnits != NULL) {
					strTimeCalendar = attTimeCalendar->as_string(0);
				}


				ParseTimeDouble(
					strTimeUnits,
					strTimeCalendar,
					dTime[t],
					nDateYear,
					nDateMonth,
					nDateDay,
					nDateHour);
//			}

			// Write time information
			fprintf(fpOutput, "%i\t%i\t%i\t%i\t%i\n",
				nDateYear,
				nDateMonth,
				nDateDay,
				static_cast<int>(setCandidates.size()),
				nDateHour);

			// Write candidate information
			int iCandidateCount = 0;

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				fprintf(fpOutput, "%i\t%i\t%i\t%3.6f\t%3.6f\t%2.6f\t%3.6f\t%6.6f\n",
					iCandidateCount,
					iterCandidate->second,
					iterCandidate->first,
					dataLon[iterCandidate->second] * 180.0 / M_PI,
					dataLat[iterCandidate->first]  * 180.0 / M_PI,
					vecMaxWindSp[iCandidateCount],
					vecRMaxWindSp[iCandidateCount],
					dataPSL[iterCandidate->first][iterCandidate->second]);

				iCandidateCount++;
			}
		}

		AnnounceEndBlock("Done");
	}

	fclose(fpOutput);

	ncInput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

