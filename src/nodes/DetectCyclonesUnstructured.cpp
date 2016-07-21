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

#include "Variable.h"
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
#include <string>
#include <set>
#include <queue>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

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
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Distance,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Parse variable
		Variable var;
		int iLast = var.ParseFromString(varreg, strOp) + 1;
		m_varix = varreg.FindOrRegister(var);

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
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
		std::string strDescription = var.ToString(varreg);
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

		char szBuffer[128];
		sprintf(szBuffer, "%f within %f degrees",
			m_dValue, m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for thresholding.
	///	</summary>
	VariableIndex m_varix;

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
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Amount,
			ReadMode_Distance,
			ReadMode_MinMaxDist,
			ReadMode_Invalid
		} eReadMode = ReadMode_Amount;

		// Get variable information
		Variable var;
		int iLast = var.ParseFromString(varreg, strOp) + 1;
		m_varix = varreg.FindOrRegister(var);

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in amount
				if (eReadMode == ReadMode_Amount) {
					m_dDeltaAmount = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in distance
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

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
						"<minmaxdist>\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in closed contour op \"%s\""
					"\nRequired: \"<name>,<amount>,<distance>,"
					"<minmaxdist>\"",
					strOp.c_str());
		}

		// Output announcement
		char szBuffer[128];

		std::string strDescription;
		strDescription += var.ToString(varreg);

		if (m_dDeltaAmount == 0.0) {
			_EXCEPTIONT("For closed contour op, delta amount must be non-zero");
		}
		if (m_dDistance <= 0.0) {
			_EXCEPTIONT("For closed contour op, distance must be positive");
		}
		if (m_dMinMaxDist < 0.0) {
			_EXCEPTIONT("For closed contour op, min/max dist must be nonnegative");
		}

		if (m_dDeltaAmount < 0.0) {
			Announce("%s decreases by %f over %f degrees"
				   " (max search %f deg)",
				var.ToString(varreg).c_str(),
				-m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);

		} else {
			Announce("%s increases by %f over %f degrees"
					" (min search %f deg)",
				var.ToString(varreg).c_str(),
				m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);
		}
	}

public:
	///	<summary>
	///		Variable to use for closed contour op.
	///	</summary>
	VariableIndex m_varix;

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
	///		Distance to search for min or max.
	///	</summary>
	double m_dMinMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing an output operator.
///	</summary>
class OutputOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		Max,
		Min,
		Avg,
		MaxDist,
		MinDist
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Distance,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Get variable information
		Variable var;
		int iLast = var.ParseFromString(varreg, strOp) + 1;
		m_varix = varreg.FindOrRegister(var);

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
					if (strSubStr == "max") {
						m_eOp = Max;
					} else if (strSubStr == "min") {
						m_eOp = Min;
					} else if (strSubStr == "avg") {
						m_eOp = Avg;
					} else if (strSubStr == "maxdist") {
						m_eOp = MaxDist;
					} else if (strSubStr == "mindist") {
						m_eOp = MinDist;
					} else {
						_EXCEPTION1("Output invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in minimum count
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in output op \"%s\""
							"\nRequired: \"<name>,<operation>,<distance>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in output op \"%s\""
					"\nRequired: \"<name>,<operation>,<distance>\"",
					strOp.c_str());
		}

		if (m_dDistance < 0.0) {
			_EXCEPTIONT("For output op, distance must be nonnegative");
		}

		// Output announcement
		std::string strDescription;

		if (m_eOp == Max) {
			strDescription += "Maximum of ";
		} else if (m_eOp == Min) {
			strDescription += "Minimum of ";
		} else if (m_eOp == Avg) {
			strDescription += "Average of ";
		} else if (m_eOp == MaxDist) {
			strDescription += "Distance to maximum of ";
		} else if (m_eOp == MinDist) {
			strDescription += "Distance to minimum of ";
		}

		char szBuffer[128];

		sprintf(szBuffer, "%s", var.ToString(varreg).c_str());
		strDescription += szBuffer;

		sprintf(szBuffer, " within %f degrees", m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for output.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Distance to use when applying operation.
	///	</summary>
	double m_dDistance;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse the list of input files.
///	</summary>
void ParseInputFiles(
	const std::string & strInputFile,
	std::vector<NcFile *> & vecFiles
) {
	int iLast = 0;
	for (int i = 0; i <= strInputFile.length(); i++) {
		if ((i == strInputFile.length()) ||
		    (strInputFile[i] == ';')
		) {
			std::string strFile =
				strInputFile.substr(iLast, i - iLast);

			NcFile * pNewFile = new NcFile(strFile.c_str());

			if (!pNewFile->is_valid()) {
				_EXCEPTION1("Cannot open input file \"%s\"",
					strFile.c_str());
			}

			vecFiles.push_back(pNewFile);
			iLast = i+1;
		}
	}

	if (vecFiles.size() == 0) {
		_EXCEPTION1("No input files found in \"%s\"",
			strInputFile.c_str());
	}
}

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
		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

		if (dR >= 1.0) {
			dR = 0.0;
		} else if (dR <= -1.0) {
			dR = 180.0;
		} else {
			dR = 180.0 / M_PI * acos(dR);
		}
		if (dR != dR) {
			_EXCEPTIONT("NaN value detected");
		}

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

	// Time format is "minutes since ..."
	} else if (
	    (strTimeUnits.length() >= 14) &&
	    (strncmp(strTimeUnits.c_str(), "minutes since ", 14) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(14);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		time.AddMinutes(static_cast<int>(dTime));

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
	double dMinMaxDist
) {
	// Verify arguments
	if (dDeltaAmt == 0.0) {
		_EXCEPTIONT("Closed contour amount must be non-zero");
	}
	if (dDeltaDist <= 0.0) {
		_EXCEPTIONT("Closed contour distance must be positive");
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

	//printf("%lu %lu : %lu %lu : %1.5f %1.5f\n", ix0 % grid.m_nGridDim[1], ix0 / grid.m_nGridDim[1], ixOrigin % grid.m_nGridDim[1], ixOrigin / grid.m_nGridDim[1], dataState[ix0], dataState[ixOrigin]);

	// Set of visited nodes
	std::set<int> setNodesVisited;

	// Set of nodes to visit
	std::queue<int> queueToVisit;
	queueToVisit.push(ixOrigin);

	// Reference value
	real dRefValue = dataState[ixOrigin];

	const double dLat0 = grid.m_dLat[ixOrigin];
	const double dLon0 = grid.m_dLon[ixOrigin];

	Announce(2, "Checking (%lu) : (%1.5f %1.5f)",
		ixOrigin, dLat0, dLon0);


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

		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

		if (dR >= 1.0) {
			dR = 0.0;
		} else if (dR <= -1.0) {
			dR = 180.0;
		} else {
			dR = 180.0 / M_PI * acos(dR);
		}
		if (dR != dR) {
			_EXCEPTIONT("NaN value detected");
		}

		Announce(2, "-- (%lu) : (%1.5f %1.5f) : dx %1.5f",
			ix, dLatThis, dLonThis, dR);

		//printf("-- (%lu %lu) %1.5f %1.5f\n", ix % grid.m_nGridDim[1], ix / grid.m_nGridDim[1], dR, dataState[ix] - dRefValue);

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

		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

		if (dR >= 1.0) {
			dR = 0.0;
		} else if (dR <= -1.0) {
			dR = 180.0;
		} else {
			dR = 180.0 / M_PI * acos(dR);
		}
		if (dR != dR) {
			_EXCEPTIONT("NaN value detected");
		}

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
template <typename real>
void FindLocalAverage(
	const SimpleGrid & grid,
	const DataVector<real> & data,
	int ix0,
	double dMaxDist,
	float & dAverage
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

	// Number of points
	float dSum = 0.0;
	int nCount = 0;

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
		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

		if (dR >= 1.0) {
			dR = 0.0;
		} else if (dR <= -1.0) {
			dR = 180.0;
		} else {
			dR = 180.0 / M_PI * acos(dR);
		}
		if (dR != dR) {
			_EXCEPTIONT("NaN value detected");
		}

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local extremum
		dSum += data[ix];
		nCount++;

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	dAverage = dSum / static_cast<float>(nCount);
}

///////////////////////////////////////////////////////////////////////////////

class DetectCyclonesParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DetectCyclonesParam() :
		fpLog(NULL),
		ixSearchBy(0),
		fSearchByMinima(false),
		dMaxLatitude(0.0),
		dMinLatitude(0.0),
		dMaxLongitude(0.0),
		dMinLongitude(0.0),
		dMergeDist(0.0),
		pvecClosedContourOp(NULL),
		pvecNoClosedContourOp(NULL),
		pvecThresholdOp(NULL),
		pvecOutputOp(NULL),
		nTimeStride(1),
		fRegional(false),
		fOutputHeader(false),
		iVerbosityLevel(0)
	{ }

public:
	// Log
	FILE * fpLog;

	// Variable index to search on
	VariableIndex ixSearchBy;

	// Serach on minima
	bool fSearchByMinima;

	// Maximum latitude for detection
	double dMaxLatitude;

	// Minimum latitude for detection
	double dMinLatitude;

	// Maximum longitude for detection
	double dMaxLongitude;

	// Minimum longitude for detection
	double dMinLongitude;

	// Merge distance
	double dMergeDist;

	// Vector of closed contour operators
	std::vector<ClosedContourOp> * pvecClosedContourOp;

	// Vector of no closed contour operators
	std::vector<ClosedContourOp> * pvecNoClosedContourOp;

	// Vector of threshold operators
	std::vector<ThresholdOp> * pvecThresholdOp;

	// Vector of output operators
	std::vector<OutputOp> * pvecOutputOp;

	// Time stride
	int nTimeStride;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Output header
	bool fOutputHeader;

	// Verbosity level
	int iVerbosityLevel;

};

///////////////////////////////////////////////////////////////////////////////

void DetectCyclonesUnstructured(
	int iFile,
	const std::string & strInputFiles,
	const std::string & strOutputFile,
	const std::string & strConnectivity,
	VariableRegistry & varreg,
	const DetectCyclonesParam & param
) {
	// Set the Announce buffer
	if (param.fpLog == NULL) {
		_EXCEPTIONT("Invalid log buffer");
	}

	AnnounceSetOutputBuffer(param.fpLog);
	AnnounceOutputOnAllRanks();

	// Check minimum longitude / latitude
	if ((param.dMinLongitude < 0.0) || (param.dMinLongitude >= 360.0)) {
		_EXCEPTIONT("Invalid MinLongitude");
	}
	if ((param.dMaxLongitude < 0.0) || (param.dMaxLongitude >= 360.0)) {
		_EXCEPTIONT("Invalid MaxLongitude");
	}
	if ((param.dMaxLatitude < -90.0) || (param.dMaxLatitude > 90.0)) {
		_EXCEPTIONT("--maxlat must in the range [-90,90]");
	}
	if ((param.dMinLatitude < -90.0) || (param.dMinLatitude > 90.0)) {
		_EXCEPTIONT("--minlat must in the range [-90,90]");
	}

	// Dereference pointers to operators
	std::vector<ClosedContourOp> & vecClosedContourOp =
		*(param.pvecClosedContourOp);

	std::vector<ClosedContourOp> & vecNoClosedContourOp =
		*(param.pvecNoClosedContourOp);

	std::vector<ThresholdOp> & vecThresholdOp =
		*(param.pvecThresholdOp);

	std::vector<OutputOp> & vecOutputOp =
		*(param.pvecOutputOp);

	// Unload data from the VariableRegistry
	varreg.UnloadAllGridData();

	// Define the SimpleGrid
	SimpleGrid grid;

	// Dimensions
	int nSize = 0;
	int nLon = 0;
	int nLat = 0;

	// Load in the benchmark file
	NcFileVector vecFiles;

	ParseInputFiles(strInputFiles, vecFiles);

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

		nSize = grid.GetSize();

	// No connectivity file; check for latitude/longitude dimension
	} else {

		NcDim * dimLat = vecFiles[0]->get_dim("lat");
		if (dimLat == NULL) {
			_EXCEPTIONT("No dimension \"lat\" found in input file");
		}

		NcDim * dimLon = vecFiles[0]->get_dim("lon");
		if (dimLon == NULL) {
			_EXCEPTIONT("No dimension \"lon\" found in input file");
		}

		NcVar * varLat = vecFiles[0]->get_var("lat");
		if (varLat == NULL) {
			_EXCEPTIONT("No variable \"lat\" found in input file");
		}

		NcVar * varLon = vecFiles[0]->get_var("lon");
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
		grid.GenerateLatitudeLongitude(vecLat, vecLon, param.fRegional);
	}

	// Get time dimension
	NcDim * dimTime = vecFiles[0]->get_dim("time");
	if (dimTime == NULL) {
		_EXCEPTIONT("No dimension \"time\" found in first input file");
	}

	NcVar * varTime = vecFiles[0]->get_var("time");
	if (varTime == NULL) {
		_EXCEPTIONT("No variable \"time\" found in input file");
	}

	int nTime = dimTime->size();

	DataVector<double> dTime;
	dTime.Initialize(nTime);

	if (varTime->type() == ncDouble) {
		varTime->get(dTime, nTime);

	} else if (varTime->type() == ncFloat) {
		DataVector<float> dTimeFloat;
		dTimeFloat.Initialize(nTime);

		varTime->get(dTimeFloat, nTime);
		for (int t = 0; t < nTime; t++) {
			dTime[t] = static_cast<double>(dTimeFloat[t]);
		}

	} else if (varTime->type() == ncInt) {
		DataVector<int> dTimeInt;
		dTimeInt.Initialize(nTime);

		varTime->get(dTimeInt, nTime);
		for (int t = 0; t < nTime; t++) {
			dTime[t] = static_cast<double>(dTimeInt[t]);
		}

	} else {
		_EXCEPTIONT("Variable \"time\" has an invalid type:\n"
			"Expected \"float\", \"double\" or \"int\"");
	}

	// Open output file
	FILE * fpOutput = fopen(strOutputFile.c_str(), "w");
	if (fpOutput == NULL) {
		_EXCEPTION1("Could not open output file \"%s\"",
			strOutputFile.c_str());
	}

	if (param.fOutputHeader) {
		fprintf(fpOutput, "#year\tmonth\tday\tcount\thour\n");

		if (grid.m_nGridDim.size() == 1) {
			fprintf(fpOutput, "#\ti\tlon\tlat");
		} else {
			fprintf(fpOutput, "#\ti\tj\tlon\tlat");
		}

		for (int i = 0; i < vecOutputOp.size(); i++) {
			Variable & varOp = varreg.Get(vecOutputOp[i].m_varix);
			fprintf(fpOutput, "\t%s", varOp.ToString(varreg).c_str());
		}
		fprintf(fpOutput, "\n");
	}

	// Loop through all times
	for (int t = 0; t < nTime; t += param.nTimeStride) {
	//for (int t = 0; t < 1; t++) {

		char szStartBlock[128];
		sprintf(szStartBlock, "Time %i", t);
		AnnounceStartBlock(szStartBlock);

		// Load the data for the search variable
		Variable & varSearchBy = varreg.Get(param.ixSearchBy);
		varSearchBy.LoadGridData(varreg, vecFiles, grid, t);

		const DataVector<float> & dataSearch = varSearchBy.GetData();

		// Tag all minima
		std::set<int> setCandidates;

		if (param.fSearchByMinima) {
			FindAllLocalMinima<float>(grid, dataSearch, setCandidates);
		} else {
			FindAllLocalMaxima<float>(grid, dataSearch, setCandidates);
		}

		// Total number of candidates
		int nTotalCandidates = setCandidates.size();

		int nRejectedLocation = 0;
		int nRejectedTopography = 0;
		int nRejectedMerge = 0;

		DataVector<int> vecRejectedClosedContour(vecClosedContourOp.size());
		DataVector<int> vecRejectedNoClosedContour(vecNoClosedContourOp.size());
		DataVector<int> vecRejectedThreshold(vecThresholdOp.size());

		// Eliminate based on interval
		if ((param.dMinLatitude != param.dMaxLatitude) ||
		    (param.dMinLongitude != param.dMaxLongitude)
		) {
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];
				double dLon = grid.m_dLon[*iterCandidate];

				if (param.dMinLatitude != param.dMaxLatitude) {
					if (dLat < param.dMinLatitude) {
						nRejectedLocation++;
						continue;
					}
					if (dLat > param.dMaxLatitude) {
						nRejectedLocation++;
						continue;
					}
				}
				if (param.dMinLongitude != param.dMaxLongitude) {
					if (dLon < 0.0) {
						int iLonShift = static_cast<int>(dLon / (2.0 * M_PI));
						dLon += static_cast<double>(iLonShift + 1) * 2.0 * M_PI;
					}
					if (dLon >= 2.0 * M_PI) {
						int iLonShift = static_cast<int>(dLon / (2.0 * M_PI));
						dLon -= static_cast<double>(iLonShift - 1) * 2.0 * M_PI;
					}
					if (param.dMinLongitude < param.dMaxLongitude) {
						if (dLon < param.dMinLongitude) {
							nRejectedLocation++;
							continue;
						}
						if (dLon > param.dMaxLongitude) {
							nRejectedLocation++;
							continue;
						}

					} else {
						if ((dLon > param.dMaxLongitude) &&
						    (dLon < param.dMinLongitude)
						) {
							nRejectedLocation++;
							continue;
						}
					}
				}
				setNewCandidates.insert(*iterCandidate);
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on merge distance
		if (param.dMergeDist != 0.0) {
			std::set<int> setNewCandidates;

			// Calculate chord distance
			double dSphDist =
				2.0 * sin(0.5 * param.dMergeDist / 180.0 * M_PI);

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

						if (param.fSearchByMinima) {
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

		// Eliminate based on thresholds
		for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

			std::set<int> setNewCandidates;

			// Load the search variable data
			Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);
			var.LoadGridData(varreg, vecFiles, grid, t);
			const DataVector<float> & dataState = var.GetData();

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

		// Eliminate based on closed contours
		for (int ccc = 0; ccc < vecClosedContourOp.size(); ccc++) {
			std::set<int> setNewCandidates;

			// Load the search variable data
			Variable & var = varreg.Get(vecClosedContourOp[ccc].m_varix);
			var.LoadGridData(varreg, vecFiles, grid, t);
			const DataVector<float> & dataState = var.GetData();

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if a closed contour is present
				bool fHasClosedContour =
					HasClosedContour<float>(
						grid,
						dataState,
						*iterCandidate,
						vecClosedContourOp[ccc].m_dDeltaAmount,
						vecClosedContourOp[ccc].m_dDistance,
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

		// Eliminate based on no closed contours
		for (int ccc = 0; ccc < vecNoClosedContourOp.size(); ccc++) {
			std::set<int> setNewCandidates;

			// Load the search variable data
			Variable & var = varreg.Get(vecNoClosedContourOp[ccc].m_varix);
			var.LoadGridData(varreg, vecFiles, grid, t);
			const DataVector<float> & dataState = var.GetData();

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if a closed contour is present
				bool fHasClosedContour =
					HasClosedContour<float>(
						grid,
						dataState,
						*iterCandidate,
						vecNoClosedContourOp[ccc].m_dDeltaAmount,
						vecNoClosedContourOp[ccc].m_dDistance,
						vecNoClosedContourOp[ccc].m_dMinMaxDist
					);

				// If a closed contour is present, reject this candidate
				if (fHasClosedContour) {
					vecRejectedNoClosedContour[ccc]++;
				} else {
					setNewCandidates.insert(*iterCandidate);
				}
			}

			setCandidates = setNewCandidates;
		}

		Announce("Total candidates: %i", setCandidates.size());
		Announce("Rejected (  location): %i", nRejectedLocation);
		Announce("Rejected (topography): %i", nRejectedTopography);
		Announce("Rejected (    merged): %i", nRejectedMerge);

		for (int tc = 0; tc < vecRejectedThreshold.GetRows(); tc++) {
			Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);

			Announce("Rejected (thresh. %s): %i",
					var.m_strName.c_str(),
					vecRejectedThreshold[tc]);
		}

		for (int ccc = 0; ccc < vecRejectedClosedContour.GetRows(); ccc++) {
			Variable & var = varreg.Get(vecClosedContourOp[ccc].m_varix);

			Announce("Rejected (contour %s): %i",
					var.m_strName.c_str(),
					vecRejectedClosedContour[ccc]);
		}

		for (int ccc = 0; ccc < vecRejectedNoClosedContour.GetRows(); ccc++) {
			Variable & var = varreg.Get(vecNoClosedContourOp[ccc].m_varix);

			Announce("Rejected (nocontour %s): %i",
					var.m_strName.c_str(),
					vecRejectedNoClosedContour[ccc]);
		}

		// Write results to file
		{
			// Parse time information
			//NcVar * varDate = ncInput.get_var("date");
			//NcVar * varDateSec = ncInput.get_var("datesec");

			int nDateYear;
			int nDateMonth;
			int nDateDay;
			int nDateHour;

			NcAtt * attTimeUnits = varTime->get_att("units");
			if (attTimeUnits == NULL) {
				_EXCEPTIONT("Variable \"time\" has no \"units\" attribute");
			}

			std::string strTimeUnits = attTimeUnits->as_string(0);

			std::string strTimeCalendar = "standard";
			NcAtt * attTimeCalendar = varTime->get_att("calendar");
			if (attTimeCalendar != NULL) {
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

			// Write time information
			fprintf(fpOutput, "%i\t%i\t%i\t%i\t%i\n",
				nDateYear,
				nDateMonth,
				nDateDay,
				static_cast<int>(setCandidates.size()),
				nDateHour);

			// Write candidate information
			int iCandidateCount = 0;

			// Apply output operators
			DataMatrix<float> dOutput(setCandidates.size(), vecOutputOp.size());
			for (int outc = 0; outc < vecOutputOp.size(); outc++) {

				// Load the search variable data
				Variable & var = varreg.Get(vecOutputOp[outc].m_varix);
				var.LoadGridData(varreg, vecFiles, grid, t);
				const DataVector<float> & dataState = var.GetData();

				// Loop through all pressure minima
				std::set<int>::const_iterator iterCandidate
					= setCandidates.begin();

				iCandidateCount = 0;
				for (; iterCandidate != setCandidates.end(); iterCandidate++) {

					int ixExtremum;
					float dValue;
					float dRMax;

					if (vecOutputOp[outc].m_eOp == OutputOp::Max) {
						FindLocalMinMax<float>(
							grid,
							false,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dValue;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::MaxDist) {
						FindLocalMinMax<float>(
							grid,
							false,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dRMax;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::Min) {
						FindLocalMinMax<float>(
							grid,
							true,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dValue;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::MinDist) {
						FindLocalMinMax<float>(
							grid,
							true,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dRMax;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::Avg) {
						FindLocalAverage<float>(
							grid,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							dValue);

						dOutput[iCandidateCount][outc] = dValue;

					} else {
						_EXCEPTIONT("Invalid Output operator");
					}

					iCandidateCount++;
				}
			}

			// Output all candidates
			iCandidateCount = 0;

			std::set<int>::const_iterator iterCandidate = setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				if (grid.m_nGridDim.size() == 1) {
					fprintf(fpOutput, "\t%i", *iterCandidate);

				} else if (grid.m_nGridDim.size() == 2) {
					fprintf(fpOutput, "\t%i\t%i",
						(*iterCandidate) % static_cast<int>(grid.m_nGridDim[1]),
						(*iterCandidate) / static_cast<int>(grid.m_nGridDim[1]));
				}

				fprintf(fpOutput, "\t%3.6f\t%3.6f",
					grid.m_dLon[*iterCandidate] * 180.0 / M_PI,
					grid.m_dLat[*iterCandidate] * 180.0 / M_PI);

				for (int outc = 0; outc < vecOutputOp.size(); outc++) {
					fprintf(fpOutput, "\t%3.6e",
						dOutput[iCandidateCount][outc]);
				}

				fprintf(fpOutput, "\n");

				iCandidateCount++;
			}
		}

		AnnounceEndBlock("Done");
	}

	fclose(fpOutput);

	for (int i = 0; i < vecFiles.size(); i++) {
		vecFiles[i]->close();
	}

	// Reset the Announce buffer
	AnnounceSetOutputBuffer(stdout);
	AnnounceOnlyOutputOnRankZero();
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {
	// Parameters for DetectCycloneUnstructured
	DetectCyclonesParam dcuparam;

	// Input dat file
	std::string strInputFile;

	// Input list file
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutput;

	// Output file list
	std::string strOutputFileList;

	// Variable to search for the minimum
	std::string strSearchByMin;

	// Variable to search for the maximum
	std::string strSearchByMax;
/*
	// Maximum latitude for detection
	double dMaxLatitude;

	// Minimum latitude for detection
	double dMinLatitude;

	// Maximum longitude for detection
	double dMaxLongitude;

	// Minimum longitude for detection
	double dMinLongitude;

	// Merge distance
	double dMergeDist;
*/
	// Closed contour commands
	std::string strClosedContourCmd;

	// Closed contour commands
	std::string strNoClosedContourCmd;

	// Threshold commands
	std::string strThresholdCmd;

	// Output commands
	std::string strOutputCmd;
/*
	// Time stride
	int nTimeStride;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Output header
	bool fOutputHeader;

	// Verbosity level
	int iVerbosityLevel;
*/
	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutput, "out", "");
		CommandLineString(strOutputFileList, "out_file_list", "");
		CommandLineStringD(strSearchByMin, "searchbymin", "", "(default PSL)");
		CommandLineString(strSearchByMax, "searchbymax", "");
		CommandLineDoubleD(dcuparam.dMinLongitude, "minlon", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMaxLongitude, "maxlon", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMergeDist, "mergedist", 0.0, "(degrees)");
		CommandLineStringD(strClosedContourCmd, "closedcontourcmd", "", "[var,delta,dist,minmaxdist;...]");
		CommandLineStringD(strNoClosedContourCmd, "noclosedcontourcmd", "", "[var,delta,dist,minmaxdist;...]");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		CommandLineStringD(strOutputCmd, "outputcmd", "", "[var,op,dist;...]");
		CommandLineInt(dcuparam.nTimeStride, "timestride", 1);
		CommandLineBool(dcuparam.fRegional, "regional");
		CommandLineBool(dcuparam.fOutputHeader, "out_header");
		CommandLineInt(dcuparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(dcuparam.iVerbosityLevel);

	// Check input
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}

	// Check output
	if ((strOutput.length() != 0) && (strOutputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out) or (--out_data_list)"
			" may be specified");
	}

	// Load input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		std::ifstream ifInputFileList(strInputFileList);
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputFiles.push_back(strFileLine);
		}
	}

	// Load output file list
	std::vector<std::string> vecOutputFiles;

	if (strOutputFileList.length() != 0) {

		std::ifstream ifOutputFileList(strOutputFileList);
		if (!ifOutputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strOutputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecOutputFiles.push_back(strFileLine);
		}

		if (vecOutputFiles.size() != vecInputFiles.size()) {
			_EXCEPTIONT("File --in_file_list must match --out_file_list");
		}
	}

	// Only one of search by min or search by max should be specified
	if ((strSearchByMin == "") && (strSearchByMax == "")) {
		strSearchByMin = "PSL";
	}
	if ((strSearchByMin != "") && (strSearchByMax != "")) {
		_EXCEPTIONT("Only one of --searchbymin or --searchbymax can"
			" be specified");
	}

	dcuparam.fSearchByMinima = false;
	{
		Variable varSearchByArg;
		if (strSearchByMin != "") {
			varSearchByArg.ParseFromString(varreg, strSearchByMin);
			dcuparam.fSearchByMinima = true;
		}
		if (strSearchByMax != "") {
			varSearchByArg.ParseFromString(varreg, strSearchByMax);
			dcuparam.fSearchByMinima = false;
		}

		dcuparam.ixSearchBy = varreg.FindOrRegister(varSearchByArg);
	}

	// Parse the closed contour command string
	std::vector<ClosedContourOp> vecClosedContourOp;
	dcuparam.pvecClosedContourOp = &vecClosedContourOp;

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
				vecClosedContourOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the no closed contour command string
	std::vector<ClosedContourOp> vecNoClosedContourOp;
	dcuparam.pvecNoClosedContourOp = &vecNoClosedContourOp;

	if (strNoClosedContourCmd != "") {
		AnnounceStartBlock("Parsing no closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strNoClosedContourCmd.length(); i++) {

			if ((i == strNoClosedContourCmd.length()) ||
				(strNoClosedContourCmd[i] == ';') ||
				(strNoClosedContourCmd[i] == ':')
			) {
				std::string strSubStr =
					strNoClosedContourCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecNoClosedContourOp.size());
				vecNoClosedContourOp.resize(iNextOp + 1);
				vecNoClosedContourOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the threshold operator command string
	std::vector<ThresholdOp> vecThresholdOp;
	dcuparam.pvecThresholdOp = &vecThresholdOp;

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
				vecThresholdOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the output operator command string
	std::vector<OutputOp> vecOutputOp;
	dcuparam.pvecOutputOp = &vecOutputOp;

	if (strOutputCmd != "") {
		AnnounceStartBlock("Parsing output operations");

		int iLast = 0;
		for (int i = 0; i <= strOutputCmd.length(); i++) {

			if ((i == strOutputCmd.length()) ||
				(strOutputCmd[i] == ';') ||
				(strOutputCmd[i] == ':')
			) {
				std::string strSubStr =
					strOutputCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecOutputOp.size());
				vecOutputOp.resize(iNextOp + 1);
				vecOutputOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Check minimum/maximum latitude/longitude
	if ((dcuparam.dMaxLatitude < -90.0) || (dcuparam.dMaxLatitude > 90.0)) {
		_EXCEPTIONT("--maxlat must in the range [-90,90]");
	}
	if ((dcuparam.dMinLatitude < -90.0) || (dcuparam.dMinLatitude > 90.0)) {
		_EXCEPTIONT("--minlat must in the range [-90,90]");
	}
	if (dcuparam.dMinLatitude > dcuparam.dMaxLatitude) {
		_EXCEPTIONT("--minlat must be less than --maxlat");
	}

	dcuparam.dMaxLatitude *= M_PI / 180.0;
	dcuparam.dMinLatitude *= M_PI / 180.0;

	if (dcuparam.dMinLongitude < 0.0) {
		int iMinLongitude =
			static_cast<int>(-dcuparam.dMinLongitude / 360.0);
		dcuparam.dMinLongitude +=
			static_cast<double>(iMinLongitude + 1) * 360.0;
	}
	if (dcuparam.dMinLongitude >= 360.0) {
		int iMinLongitude =
			static_cast<int>(dcuparam.dMinLongitude / 360.0);
		dcuparam.dMinLongitude -=
			static_cast<double>(iMinLongitude - 1) * 360.0;
	}
	if (dcuparam.dMaxLongitude < 0.0) {
		int iMaxLongitude =
			static_cast<int>(-dcuparam.dMaxLatitude / 360.0);
		dcuparam.dMaxLongitude +=
			static_cast<double>(iMaxLongitude + 1) * 360.0;
	}
	if (dcuparam.dMaxLongitude >= 360.0) {
		int iMaxLongitude =
			static_cast<int>(dcuparam.dMaxLongitude / 360.0);
		dcuparam.dMaxLongitude -=
			static_cast<double>(iMaxLongitude - 1) * 360.0;
	}

	dcuparam.dMaxLongitude *= M_PI / 180.0;
	dcuparam.dMinLongitude *= M_PI / 180.0;

#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
#endif

	AnnounceStartBlock("Begin search operation");
	if (vecInputFiles.size() != 1) {
		if (vecOutputFiles.size() != 0) {
			Announce("Output will be written following --out_file_list");
		} else if (strOutput == "") {
			Announce("Output will be written to outXXXXXX.dat");
		} else {
			Announce("Output will be written to %sXXXXXX.dat",
				strOutput.c_str());
		}
		Announce("Logs will be written to logXXXXXX.txt");
	}

	// Loop over all files to be processed
	for (int f = 0; f < vecInputFiles.size(); f++) {
#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}
#endif
		// Generate output file name
		std::string strOutputFile;
		if (vecInputFiles.size() == 1) {
			dcuparam.fpLog = stdout;

			if (strOutput == "") {
				strOutputFile = "out.dat";
			} else {
				strOutputFile = strOutput;
			}

		} else {
			char szFileIndex[32];
			sprintf(szFileIndex, "%06i", f);

			if (vecOutputFiles.size() != 0) {
				strOutputFile = vecOutputFiles[f];
			} else {
				if (strOutput == "") {
					strOutputFile =
						"out" + std::string(szFileIndex) + ".dat";
				} else {
					strOutputFile =
						strOutput + std::string(szFileIndex) + ".dat";
				}
			}

			std::string strLogFile = "log" + std::string(szFileIndex) + ".txt";
			dcuparam.fpLog = fopen(strLogFile.c_str(), "w");
		}

		// Perform DetectCyclonesUnstructured
		DetectCyclonesUnstructured(
			f,
			vecInputFiles[f],
			strOutputFile,
			strConnectivity,
			varreg,
			dcuparam);

		// Close the log file
		if (vecInputFiles.size() != 1) {
			fclose(dcuparam.fpLog);
		}
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}

