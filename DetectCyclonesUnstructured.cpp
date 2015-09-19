///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectCyclonesUnstructured.cpp
///	\author  Paul Ullrich
///	\version September 18, 2015
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
///		A data structure describing the grid, including coordinates of
///		each data point and graph connectivity of elements.
///	</summary>
class SimpleGrid {

public:
	///	<summary>
	///		Generate the unstructured grid information for a
	///		longitude-latitude grid.
	///	</summary>
	void GenerateLatitudeLongitude(
		DataVector<double> vecLat,
		DataVector<double> vecLon
	) {
		int nLat = vecLat.GetRows();
		int nLon = vecLon.GetRows();

		m_dLat.Initialize(nLon * nLat);
		m_dLon.Initialize(nLon * nLat);
		m_vecConnectivity.resize(nLon * nLat);

		m_nGridDim.resize(2);
		m_nGridDim[0] = nLat;
		m_nGridDim[1] = nLon;

		int ixs = 0;
		for (int j = 0; j < nLat; j++) {
		for (int i = 0; i < nLon; i++) {

			// Vectorize coordinates
			m_dLat[ixs] = vecLat[j];
			m_dLon[ixs] = vecLon[i];

			// Connectivity in each compass direction
			if (j != 0) {
				m_vecConnectivity[ixs].push_back((j-1) * nLon + i);
			}
			if (j != nLat-1) {
				m_vecConnectivity[ixs].push_back((j+1) * nLon + i);
			}

			m_vecConnectivity[ixs].push_back(j * nLon + ((i + 1) % nLon));
			m_vecConnectivity[ixs].push_back(j * nLon + ((i + nLon - 1) % nLon));

			ixs++;
		}
		}

	}

	///	<summary>
	///		Load the grid information from a file.
	///	</summary>
	void FromFile(
		const std::string & strGridInfoFile
	) {
		FILE * fp = fopen(strGridInfoFile.c_str(), "r");
		if (fp == NULL) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strGridInfoFile.c_str());
		}

		size_t nFaces;
		fscanf(fp, "%lu", &nFaces);

		m_nGridDim.resize(1);
		m_nGridDim[0] = nFaces;

		m_dLon.Initialize(nFaces);
		m_dLat.Initialize(nFaces);
		m_vecConnectivity.resize(nFaces);

		for (size_t f = 0; f < nFaces; f++) {
			size_t sNeighbors;
			fscanf(fp, "%lf", &(m_dLon[f]));
			fscanf(fp, "%lf", &(m_dLat[f]));
			fscanf(fp, "%lu", &sNeighbors);

			m_vecConnectivity[f].resize(sNeighbors);
			for (size_t n = 0; n < sNeighbors; n++) {
				fscanf(fp, "%i", &(m_vecConnectivity[f][n]));
				m_vecConnectivity[f][n]--;
			}
			if (feof(fp)) {
				_EXCEPTIONT("Premature end of file");
			}
		}
	}

	///	<summary>
	///		Get the size of the SimpleGrid (number of points).
	///	</summary>
	size_t GetSize() const {
		return (m_vecConnectivity.size());
	}

public:
	///	<summary>
	///		Longitude of each grid point.
	///	</summary>
	DataVector<double> m_dLon;

	///	<summary>
	///		Latitude of each grid point.
	///	</summary>
	DataVector<double> m_dLat;

	///	<summary>
	///		Connectivity of each grid point.
	///	</summary>
	std::vector< std::vector<int> > m_vecConnectivity;

	///	<summary>
	///		Grid dimensions.
	///	</summary>
	std::vector<size_t> m_nGridDim;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a thresholding operator.
///	</summary>
class ThresholdOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		GreaterThan,
		LessThan,
		GreaterThanEqualTo,
		LessThanEqualTo,
		EqualTo,
		NotEqualTo
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
			ReadMode_Variable,
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Distance,
			ReadMode_Invalid
		} eReadMode = ReadMode_Variable;

		// Loop through string
		int iLast = 0;
		for (int i = 0; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in column name
				if (eReadMode == ReadMode_Variable) {

					m_strVariable = strSubStr;

					iLast = i + 1;
					eReadMode = ReadMode_Op;

				// Read in operation
				} else if (eReadMode == ReadMode_Op) {
					if (strSubStr == ">") {
						m_eOp = GreaterThan;
					} else if (strSubStr == "<") {
						m_eOp = LessThan;
					} else if (strSubStr == ">=") {
						m_eOp = GreaterThanEqualTo;
					} else if (strSubStr == "<=") {
						m_eOp = LessThanEqualTo;
					} else if (strSubStr == "=") {
						m_eOp = EqualTo;
					} else if (strSubStr == "!=") {
						m_eOp = NotEqualTo;
					} else {
						_EXCEPTION1("Threshold invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Value;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					m_dValue = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in minimum count
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
							"\nRequired: \"<name>,<operation>"
							",<value>,<distance>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
					"\nRequired: \"<name>,<operation>,<value>,<distance>\"",
					strOp.c_str());
		}

		if (m_dDistance < 0.0) {
			_EXCEPTIONT("For threshold op, distance must be nonnegative");
		}

		// Output announcement
		char szBuffer[128];
		sprintf(szBuffer, "%s", m_strVariable.c_str());

		std::string strDescription = szBuffer;
		if (m_eOp == GreaterThan) {
			strDescription += " is greater than ";
		} else if (m_eOp == LessThan) {
			strDescription += " is less than ";
		} else if (m_eOp == GreaterThanEqualTo) {
			strDescription += " is greater than or equal to ";
		} else if (m_eOp == LessThanEqualTo) {
			strDescription += " is less than or equal to ";
		} else if (m_eOp == EqualTo) {
			strDescription += " is equal to ";
		} else if (m_eOp == NotEqualTo) {
			strDescription += " is not equal to ";
		}

		sprintf(szBuffer, "%f within %f degrees",
			m_dValue, m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for thresholding.
	///	</summary>
	std::string m_strVariable;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Distance to search for threshold value
	///	</summary>
	double m_dDistance;
};

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

///	<summary>
///		Find the locations of all minima in the given DataMatrix.
///	</summary>
template <typename real>
void FindAllLocalMinima(
	const SimpleGrid & grid,
	const DataVector<real> & data,
	std::set<int> & setMinima
) {
	int sFaces = grid.m_vecConnectivity.size();
	for (int f = 0; f < sFaces; f++) {
		
		bool fMinimum = true;

		real dValue = data[f];
		int sNeighbors = grid.m_vecConnectivity[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (data[grid.m_vecConnectivity[f][n]] < dValue) {
				fMinimum = false;
				break;
			}
		}

		if (fMinimum) {
			setMinima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all maxima in the given DataMatrix.
///	</summary>
template <typename real>
void FindAllLocalMaxima(
	const SimpleGrid & grid,
	const DataVector<real> & data,
	std::set<int> & setMaxima
) {
	int sFaces = grid.m_vecConnectivity.size();
	for (int f = 0; f < sFaces; f++) {
		
		bool fMaximum = true;

		real dValue = data[f];
		int sNeighbors = grid.m_vecConnectivity[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (data[grid.m_vecConnectivity[f][n]] > dValue) {
				fMaximum = false;
				break;
			}
		}

		if (fMaximum) {
			setMaxima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the minimum/maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
template <typename real>
void FindLocalMinMax(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataVector<real> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	real & dMaxValue,
	float & dRMax
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Initialize the maximum to the central location
	ixExtremum = ix0;
	dMaxValue = data[ix0];
	dRMax = 0.0;

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ixExtremum);

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = grid.m_dLat[ix0];
	double dLon0 = grid.m_dLon[ix0];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		// Great circle distance to this element
		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local extremum
		if (fMinimum) {
			if (data[ix] < dMaxValue) {
				ixExtremum = ix;
				dMaxValue = data[ix];
				dRMax = dR;
			}

		} else {
			if (data[ix] > dMaxValue) {
				ixExtremum = ix;
				dMaxValue = data[ix];
				dRMax = dR;
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
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
/*
	Time time(Time::CalendarStandard);
	time.FromFormattedString("1800-01-01 00:00:00");
	printf("%1.15e %i\n", 3600.0 * 1577832.0, (int)(3600.0 * 1577832.0));
	time.AddHours(1577832);

	Announce("Time (YMDS): %i %i %i %i",
			time.GetYear(),
			time.GetMonth(),
			time.GetDay(),
			time.GetSecond());

	_EXCEPTION();
*/
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

		time.AddHours(static_cast<int>(dTime));

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
template <typename real>
bool HasClosedContour(
	const SimpleGrid & grid,
	const DataVector<real> & dataState,
	const int ix0,
	double dDeltaAmt,
	double dDeltaDist,
	int nDeltaCount,
	double dMinMaxDist
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

	// Find min/max near point
	int ixOrigin;

	if (dMinMaxDist == 0.0) {
		ixOrigin = ix0;

	// Find a local minimum / maximum
	} else {
		real dValue;
		float dR;

		FindLocalMinMax<real>(
			grid,
			(dDeltaAmt > 0.0),
			dataState,
			ix0,
			dMinMaxDist,
			ixOrigin,
			dValue,
			dR);
	}

	// Set of visited nodes
	std::set<int> setNodesVisited;

	// Set of nodes to visit
	std::queue<int> queueToVisit;
	queueToVisit.push(ixOrigin);

	// Reference value
	real dRefValue = dataState[ixOrigin];

	const double dLat0 = grid.m_dLat[ixOrigin];
	const double dLon0 = grid.m_dLon[ixOrigin];

	// Build up nodes
	while (queueToVisit.size() != 0) {
		int ix = queueToVisit.front();
		queueToVisit.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element
		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		// Check great circle distance
		if (dR > dDeltaDist) {
			Announce(2, "Failed criteria; returning");
			AnnounceEndBlock(2, NULL);
			return false;
		}

		// Verify sufficient increase in value
		if (dDeltaAmt > 0.0) {
			if (dataState[ix] - dRefValue >= dDeltaAmt) {
				continue;
			}

		// Verify sufficient decrease in value
		} else {
			if (dRefValue - dataState[ix] >= -dDeltaAmt) {
				continue;
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueToVisit.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Report success with criteria
	Announce(2, "Passed criteria; returning");
	AnnounceEndBlock(2, NULL);
	return true;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field satisfies the threshold.
///	</summary>
template <typename real>
bool SatisfiesThreshold(
	const SimpleGrid & grid,
	const DataVector<real> & dataState,
	const int ix0,
	const ThresholdOp::Operation op,
	const double dTargetValue,
	const double dMaxDist
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ix0);

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = grid.m_dLat[ix0];
	double dLon0 = grid.m_dLon[ix0];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element
		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		if ((ix != ix0) && (dR > dMaxDist)) {
			continue;
		}

		// Value at this location
		double dValue = dataState[ix];

		// Apply operator
		if (op == ThresholdOp::GreaterThan) {
			if (dValue > dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThan) {
			if (dValue < dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::GreaterThanEqualTo) {
			if (dValue >= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThanEqualTo) {
			if (dValue <= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::EqualTo) {
			if (dValue == dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::NotEqualTo) {
			if (dValue != dTargetValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		// Special case: zero distance
		if (dMaxDist == 0.0) {
			return false;
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
void FindLocalAverage(
	const DataMatrix<float> & data,
	const DataVector<double> & dataLat,
	const DataVector<double> & dataLon,
	int iLat,
	int iLon,
	double dMaxDist,
	float & dAverage,
	float & dMaxValue
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Number of latitudes and longitudes
	const int nLat = dataLat.GetRows();
	const int nLon = dataLon.GetRows();

	// Queue of nodes that remain to be visited
	std::queue< std::pair<int, int> > queueNodes;
	queueNodes.push( std::pair<int, int>(iLat, iLon) );

	// Set of nodes that have already been visited
	std::set< std::pair<int, int> > setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = dataLat[iLat];
	double dLon0 = dataLon[iLon];

	// Sum
	float dSum = 0.0f;
	float dArea = 0.0f;

	// Reset average
	dAverage = 0.0f;

	// Reset maximum value
	dMaxValue = data[iLat][iLon];

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

		// Add value to sum
		float dLocalArea = cos(dLatThis);

		dArea += dLocalArea;
		dSum += data[pr.first][pr.second] * dLocalArea;

		if (data[pr.first][pr.second] > dMaxValue) {
			dMaxValue = data[pr.first][pr.second];
		}

		// Add all neighbors of this point
		std::pair<int,int> prWest(pr.first, (pr.second + nLon - 1) % nLon);
		if (setNodesVisited.find(prWest) == setNodesVisited.end()) {
			queueNodes.push(prWest);
		}

		std::pair<int,int> prEast(pr.first, (pr.second + 1) % nLon);
		if (setNodesVisited.find(prEast) == setNodesVisited.end()) {
			queueNodes.push(prEast);
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

	// Set average
	dAverage = dSum / dArea;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load a data block from the NcFile.
///	</summary>
template <typename real>
void LoadGridData(
	NcFile & ncFile,
	const std::string & strVariable,
	const SimpleGrid & grid,
	DataVector<real> & data,
	int iDim0 = 0
) {
	// Get pointer to variable
	NcVar * var = ncFile.get_var(strVariable.c_str());
	if (var == NULL) {
		_EXCEPTION1("No variable \"%s\" found in topography file",
			strVariable.c_str());
	}

	// Grid dimensions
	int nSize = 0;
	int nLat = 0;
	int nLon = 0;

	if (grid.m_nGridDim.size() == 2) {
		nLat = grid.m_nGridDim[0];
		nLon = grid.m_nGridDim[1];
	} else if (grid.m_nGridDim.size() == 1) {
		nSize = grid.m_nGridDim[0];
	}

	// Load data
	if (var->num_dims() == 3) {
		if (nSize != 0) {
			_EXCEPTION1("%s field inconsistent with grid",
				strVariable.c_str());
		}
		if (var->get_dim(1)->size() != nLat) {
			_EXCEPTION1("%s field inconsistent with grid",
				strVariable.c_str());
		}
		if (var->get_dim(2)->size() != nLon) {
			_EXCEPTION1("%s field inconsistent with grid",
				strVariable.c_str());
		}
		var->set_cur(iDim0, 0, 0);
		var->get(&(data[0]), 1, nLat, nLon);

	} else if (var->num_dims() == 2) {
		if (nSize != 0) {
			if (var->get_dim(1)->size() != nSize) {
				_EXCEPTION1("%s field inconsistent with grid",
					strVariable.c_str());
			}
			var->set_cur(iDim0, 0);
			var->get(&(data[0]), 1, nSize);

		} else {
			if (var->get_dim(0)->size() != nLat) {
				_EXCEPTION1("%s field inconsistent with grid",
					strVariable.c_str());
			}
			if (var->get_dim(1)->size() != nLon) {
				_EXCEPTION1("%s field inconsistent with grid",
					strVariable.c_str());
			}
			var->get(&(data[0]), nLat, nLon);
		}

	} else if (var->num_dims() == 1) {
		if (var->get_dim(0)->size() != nSize) {
			_EXCEPTION1("%s field inconsistent with grid",
				strVariable.c_str());
		}
		var->get(&(data[0]), nSize);

	} else {
		_EXCEPTIONT("Invalid number of dimensions on variable \"PHIS\"");
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {
	// Gravitational constant
	const float ParamGravity = 9.80616;

	// Radius of the Earth
	const float ParamEarthRadius = 6.37122e6;

	// Input dat file
	std::string strInputFile;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutputFile;

	// Variable to search for the minimum
	std::string strSearchByMin;

	// Variable to search for the maximum
	std::string strSearchByMax;

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

	// Closed contour commands
	std::string strClosedContourCmd;

	// Threshold commands
	std::string strThresholdCmd;

	// Minimum Laplacian value
	double dMinLaplacian;

	// Distance to search for maximum wind speed
	double dWindSpDist;

	// Distance to search for precipitation statistics (avg and max)
	double dPrectDist;

	// Output header
	bool fOutputHeader;

	// Verbosity level
	int iVerbosityLevel;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineStringD(strSearchByMin, "searchbymin", "", "(default PSL)");
		CommandLineString(strSearchByMax, "searchbymax", "");
		CommandLineDoubleD(dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineString(strTopoFile, "topofile", "");
		CommandLineDoubleD(dMaxTopoHeight, "maxtopoht", 0.0, "(m)");
		CommandLineDoubleD(dMergeDist, "mergedist", 0.0, "(degrees)");
		CommandLineString(strClosedContourCmd, "closedcontourcmd", "");
		CommandLineString(strThresholdCmd, "thresholdcmd", "");
		CommandLineDoubleD(dWindSpDist, "windspdist", 0.0, "(degrees)");
		CommandLineDoubleD(dPrectDist, "prectdist", -1.0, "(degrees)");
		CommandLineBool(fOutputHeader, "out_header");
		CommandLineInt(iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Set verbosity level
	AnnounceSetVerbosityLevel(iVerbosityLevel);

	// Check input
	if (strInputFile.length() == 0) {
		_EXCEPTIONT("No input data file (--in_data) specified");
	}

	// Check output
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Only one of search by min or search by max should be specified
	if ((strSearchByMin == "") && (strSearchByMax == "")) {
		strSearchByMin = "PSL";
	}
	if ((strSearchByMin != "") && (strSearchByMax != "")) {
		_EXCEPTIONT("Only one of --searchbymin or --searchbymax can"
			" be specified");
	}

	bool fSearchByMinima = false;
	std::string strSearchBy;
	if (strSearchByMin != "") {
		strSearchBy = strSearchByMin;
		fSearchByMinima = true;
	}
	if (strSearchByMax != "") {
		strSearchBy = strSearchByMax;
		fSearchByMinima = false;
	}

	// Parse the closed contour command string
	std::vector<ClosedContourOp> vecClosedContourOp;

	if (strClosedContourCmd != "") {
		AnnounceStartBlock("Parsing closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strClosedContourCmd.length(); i++) {

			if ((i == strClosedContourCmd.length()) ||
				(strClosedContourCmd[i] == ';') ||
				(strClosedContourCmd[i] == ':')
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

	// Parse the threshold operator command string
	std::vector<ThresholdOp> vecThresholdOp;

	if (strThresholdCmd != "") {
		AnnounceStartBlock("Parsing threshold operations");

		int iLast = 0;
		for (int i = 0; i <= strThresholdCmd.length(); i++) {

			if ((i == strThresholdCmd.length()) ||
				(strThresholdCmd[i] == ';') ||
				(strThresholdCmd[i] == ':')
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
	if (!ncInput.is_valid()) {
		_EXCEPTION1("Cannot open input file \"%s\"", strInputFile.c_str());
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Dimensions
	int nSize = 0;
	int nLon = 0;
	int nLat = 0;

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

		nSize = grid.GetSize();

	// No connectivity file; check for latitude/longitude dimension
	} else {

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

		DataVector<double> vecLat(nLat);
		varLat->get(vecLat, nLat);

		for (int j = 0; j < nLat; j++) {
			vecLat[j] *= M_PI / 180.0;
		}

		DataVector<double> vecLon(nLon);
		varLon->get(vecLon, nLon);

		for (int i = 0; i < nLon; i++) {
			vecLon[i] *= M_PI / 180.0;
		}

		// Generate the SimpleGrid
		grid.GenerateLatitudeLongitude(vecLat, vecLon);
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

	// Search variable data
	DataVector<float> dataSearch(grid.GetSize());

	// Load topography data (if requested)
	NcVar * varPHIS = NULL;

	DataVector<float> dataPHIS;

	if (strTopoFile != "") {
		NcFile ncTopo(strTopoFile.c_str());
		if (!ncTopo.is_valid()) {
			_EXCEPTION1("Unable to open file \"%s\"", strTopoFile.c_str());
		}

		LoadGridData<float>(ncTopo, "PHIS", grid, dataPHIS);
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

		// Load the data for the search variable
		LoadGridData<float>(ncInput, strSearchBy.c_str(), grid, dataSearch, t);

		// Tag all minima
		std::set<int> setCandidates;

		if (strSearchByMin != "") {
			FindAllLocalMinima<float>(grid, dataSearch, setCandidates);
		} else {
			FindAllLocalMaxima<float>(grid, dataSearch, setCandidates);
		}

		// Total number of candidates
		int nTotalCandidates = setCandidates.size();

		int nRejectedLatitude = 0;
		int nRejectedTopography = 0;
		int nRejectedMerge = 0;

		DataVector<int> vecRejectedClosedContour(vecClosedContourOp.size());
		DataVector<int> vecRejectedThreshold(vecThresholdOp.size());

		// Eliminate based on maximum latitude
		if (dMaxLatitude != 0.0) {
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];

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
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];

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
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				float dPHIS = dataPHIS[*iterCandidate];

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
			std::set<int> setNewCandidates;

			// Calculate spherical distance
			double dSphDist = 2.0 * sin(0.5 * dMergeDist / 180.0 * M_PI);

			// Create a new KD Tree containing all nodes
			kdtree * kdMerge = kd_create(3);

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];
				double dLon = grid.m_dLon[*iterCandidate];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdMerge, dX, dY, dZ, (void*)(&(*iterCandidate)));
			}

			// Loop through all candidates find set of nearest neighbors
			iterCandidate = setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];
				double dLon = grid.m_dLon[*iterCandidate];

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
						static_cast<double>(dataSearch[*iterCandidate]);

					bool fExtrema = true;
					for (;;) {
						int * ppr = (int *)(kd_res_item_data(kdresMerge));

						if (fSearchByMinima) {
							if (static_cast<double>(dataSearch[*ppr]) < dValue) {
								fExtrema = false;
								break;
							}

						} else {
							if (static_cast<double>(dataSearch[*ppr]) > dValue) {
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

		// Eliminate based on closed contours
		for (int ccc = 0; ccc < vecClosedContourOp.size(); ccc++) {
			std::set<int> setNewCandidates;

			// Load the search variable data
			DataVector<float> dataState(grid.GetSize());
			
			LoadGridData<float>(
				ncInput, vecClosedContourOp[ccc].m_strVariable.c_str(),
				grid, dataState, t);

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if pressure minima have a closed contour
				bool fHasClosedContour =
					HasClosedContour<float>(
						grid,
						dataState,
						*iterCandidate,
						vecClosedContourOp[ccc].m_dDeltaAmount,
						vecClosedContourOp[ccc].m_dDistance,
						vecClosedContourOp[ccc].m_nCount,
						vecClosedContourOp[ccc].m_dMinMaxDist
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

		// Eliminate based on thresholds
		for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

			std::set<int> setNewCandidates;

			// Load the search variable data
			DataVector<float> dataState(grid.GetSize());
			
			LoadGridData<float>(
				ncInput, vecThresholdOp[tc].m_strVariable.c_str(),
				grid, dataState, t);

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if the threshold is satisfied
				bool fSatisfiesThreshold =
					SatisfiesThreshold<float>(
						grid,
						dataState,
						*iterCandidate,
						vecThresholdOp[tc].m_eOp,
						vecThresholdOp[tc].m_dValue,
						vecThresholdOp[tc].m_dDistance
					);

				// If not rejected, add to new pressure minima array
				if (fSatisfiesThreshold) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					vecRejectedThreshold[tc]++;
				}
			}

			setCandidates = setNewCandidates;
		}

		Announce("Total candidates: %i", setCandidates.size());
		Announce("Rejected (  latitude): %i", nRejectedLatitude);
		Announce("Rejected (topography): %i", nRejectedTopography);
		Announce("Rejected (    merged): %i", nRejectedMerge);

		for (int ccc = 0; ccc < vecRejectedClosedContour.GetRows(); ccc++) {
			Announce("Rejected: (contour %s): %i",
					vecClosedContourOp[ccc].m_strVariable.c_str(),
					vecRejectedClosedContour[ccc]);
		}

		for (int tc = 0; tc < vecRejectedThreshold.GetRows(); tc++) {
			Announce("Rejected: (thresh. %s): %i",
					vecThresholdOp[tc].m_strVariable.c_str(),
					vecRejectedThreshold[tc]);
		}

/*
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
					dataLat,
					dataLon,
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
*/
/*
		// Determine average and maximum PRECT
		std::vector<float> vecPrectAvg;
		std::vector<float> vecPrectMax;

		if (dPrectDist >= 0.0) {
			AnnounceStartBlock("Finding PRECT average and maximum");

			std::set< std::pair<int, int> >::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				float dAverage;
				float dMaxValue;

				FindLocalAverage(
					dataPrect,
					dataLat,
					dataLon,
					iterCandidate->first,
					iterCandidate->second,
					dPrectDist,
					dAverage,
					dMaxValue);

				vecPrectAvg.push_back(dAverage);
				vecPrectMax.push_back(dMaxValue);
			}
			AnnounceEndBlock("Done");
		}
*/
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

			std::set<int>::const_iterator iterCandidate = setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
/*
				fprintf(fpOutput, "%i\t%i\t%i\t%3.6f\t%3.6f\t%2.6f\t%3.6f\t%6.6f",
					iCandidateCount,
					iterCandidate->second,
					iterCandidate->first,
					grid.m_dLon[*iterCandidate] * 180.0 / M_PI,
					grid.m_dLat[*iterCandidate]  * 180.0 / M_PI,
					vecMaxWindSp[iCandidateCount],
					vecRMaxWindSp[iCandidateCount]);
					//dataPSL[*iterCandidate]);
*/
/*
				if (dPrectDist >= 0.0) {
					fprintf(fpOutput, "\t%1.5e\t%1.5e\n",
						vecPrectAvg[iCandidateCount],
						vecPrectMax[iCandidateCount]);
				} else {
					fprintf(fpOutput, "\n");
				}
*/
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

