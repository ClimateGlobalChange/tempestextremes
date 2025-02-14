///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeOutputOp.h
///	\author  Paul Ullrich
///	\version August 14, 2018
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _NODEOUTPUTOP_H_
#define _NODEOUTPUTOP_H_

#include "CoordTransforms.h"
#include "SimpleGridUtilities.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing an output operator.
///	</summary>
class NodeOutputOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		InvalidOperation,
		Max,
		Min,
		Avg,
		MaxDist,
		MinDist,
		MaxCoordinate,
		MinCoordinate,
		MaxIndex,
		MinIndex,
		PosClosedContour,
		NegClosedContour,
		PosMinusNegWtArea,
		MaxPoleward
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	NodeOutputOp() :
		m_varix(InvalidVariableIndex),
		m_strVariableString(""),
		m_eOp(InvalidOperation),
		m_dDistance(0.0),
		m_dMinMaxDist(0.0)
	{ }

public:
	///	<summary>
	///		Parse a output operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp,
		bool fVerbose = true
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Distance,
			ReadMode_MinMaxDist,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		std::string strOpCommand;

		// Get variable information
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
					strOpCommand = strSubStr;

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
					} else if (strSubStr == "maxcoord") {
						m_eOp = MaxCoordinate;
					} else if (strSubStr == "mincoord") {
						m_eOp = MinCoordinate;
					} else if (strSubStr == "maxix") {
						m_eOp = MaxIndex;
					} else if (strSubStr == "minix") {
						m_eOp = MinIndex;
					} else if (strSubStr == "posclosedcontour") {
						m_eOp = PosClosedContour;
					} else if (strSubStr == "negclosedcontour") {
						m_eOp = NegClosedContour;
					} else if (strSubStr == "posminusnegwtarea") {
						m_eOp = PosMinusNegWtArea;
					} else if (strSubStr == "maxpoleward") {
						m_eOp = MaxPoleward;

					} else {
						_EXCEPTION1("Output invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in distance
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					if ((m_eOp == PosClosedContour) || (m_eOp == NegClosedContour)) {
						eReadMode = ReadMode_MinMaxDist;
					} else {
						eReadMode = ReadMode_Invalid;
					}

				// Read in min-max distance
				} else if (eReadMode == ReadMode_MinMaxDist) {
					m_dMinMaxDist = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					if ((m_eOp == PosClosedContour) || (m_eOp == NegClosedContour)) {
						_EXCEPTION1("\nToo many arguments in output op \"%s\""
								"\nRequired: \"<name>,<operation>,<distance>[,<minmaxdist>]\"",
								strOp.c_str());
					} else {
						_EXCEPTION1("\nToo many arguments in output op \"%s\""
								"\nRequired: \"<name>,<operation>,<distance>\"",
								strOp.c_str());
					}
				}
			}
		}

		// Min-max distance is an optional argument
		if (eReadMode == ReadMode_MinMaxDist) {
			eReadMode = ReadMode_Invalid;
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient arguments in output op \"%s\""
					"\nRequired: \"<name>,<operation>,<distance>\"",
					strOp.c_str());
		}

		if ((m_eOp == PosClosedContour) || (m_eOp == NegClosedContour)) {
			if (m_dDistance <= 0.0) {
				_EXCEPTION1("For output op \"%s\", distance must be nonnegative", strOpCommand.c_str());
			}
			if (m_dMinMaxDist < 0.0) {
				_EXCEPTION1("For output op \"%s\", min-max distance must be nonnegative", strOpCommand.c_str());
			}

		} else{
			if (m_dDistance < 0.0) {
				_EXCEPTION1("For output op \"%s\", distance must be nonnegative", strOpCommand.c_str());
			}
		}

		// Store variable string
		m_strVariableString = varreg.GetVariableString(m_varix);

		// Output announcement
		if (fVerbose) {
			Announce("%s", ToString().c_str());
		}
	}

	///	<summary>
	///		Get the string representation of this operator.
	///	</summary>
	std::string ToString() const {
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
		} else if (m_eOp == MaxCoordinate) {
			strDescription += "Coordinates of maximum of ";
		} else if (m_eOp == MinCoordinate) {
			strDescription += "Coordinates of minimum of ";
		} else if (m_eOp == PosClosedContour) {
			strDescription += "Greatest positive closed contour delta of ";
		} else if (m_eOp == NegClosedContour) {
			strDescription += "Greatest negative closed contour delta of ";
		} else if (m_eOp == PosMinusNegWtArea) {
			strDescription += "Positive minus negative weighted area of ";
		} else if (m_eOp == MaxPoleward) {
			strDescription += "Maximum poleward value of ";
		}

		strDescription += m_strVariableString;

		char szBuffer[128];
		snprintf(szBuffer, 128, " within %f degrees", m_dDistance);
		strDescription += szBuffer;

		if (m_eOp == MaxPoleward) {
			strDescription += " longitude";
		}

		return strDescription;
	}

public:
	///	<summary>
	///		Variable to use for output.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Variable string.
	///	</summary>
	std::string m_strVariableString;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Distance to use when applying operation.
	///	</summary>
	double m_dDistance;

	///	<summary>
	///		MinMax distance.
	///	</summary>
	double m_dMinMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void ApplyNodeOutputOp(
	const NodeOutputOp & op,
	const SimpleGrid & grid,
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const Time & time,
	int ixCandidate,
	std::string & strResult
) {
	static const char * szFormat = "%3.6e";
	char buf[100];

	// Load the search variable data
	Variable & var = varreg.Get(op.m_varix);
	vecFiles.SetTime(time);
	var.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataState = var.GetData();

	// Return values from the output operators
	int ixExtremum;
	float dValue;
	float dRMax;

	// Value of the minimum or maximum value within given range
	if ((op.m_eOp == NodeOutputOp::Max) ||
	    (op.m_eOp == NodeOutputOp::Min)
	) {
		FindLocalMinMax<real>(
			grid,
			(op.m_eOp == NodeOutputOp::Min),
			dataState,
			ixCandidate,
			op.m_dDistance,
			ixExtremum,
			dValue,
			dRMax);

		snprintf(buf, 100, szFormat, dValue);
		strResult = buf;

	// Distance to the minimum or maximum value within given range
	} else if (
		(op.m_eOp == NodeOutputOp::MaxDist) ||
		(op.m_eOp == NodeOutputOp::MinDist)
	) {
		FindLocalMinMax<float>(
			grid,
			(op.m_eOp == NodeOutputOp::MinDist),
			dataState,
			ixCandidate,
			op.m_dDistance,
			ixExtremum,
			dValue,
			dRMax);

		snprintf(buf, 100, szFormat, dRMax);
		strResult = buf;

	// Coordinates (lon,lat) of the minimum or maximum value within given range
	} else if (
		(op.m_eOp == NodeOutputOp::MaxCoordinate) ||
		(op.m_eOp == NodeOutputOp::MinCoordinate)
	) {
		FindLocalMinMax<float>(
			grid,
			(op.m_eOp == NodeOutputOp::MinCoordinate),
			dataState,
			ixCandidate,
			op.m_dDistance,
			ixExtremum,
			dValue,
			dRMax);

		if ((ixExtremum <= 0) || (ixExtremum >= grid.m_dLon.GetRows())) {
			_EXCEPTIONT("Extremum position out of range");
		}
		if (ixExtremum >= grid.m_dLat.GetRows()) {
			_EXCEPTIONT("Longitude/latitude array size inconsistency in SimpleGrid");
		}
		snprintf(buf, 100, "%3.6f", RadToDeg(grid.m_dLon[ixExtremum]));
		strResult = buf;

		strResult += "\t";

		snprintf(buf, 100, "%3.6f", RadToDeg(grid.m_dLat[ixExtremum]));
		strResult += buf;

	// Coordinates (lon,lat) of the minimum or maximum value within given range
	} else if (
		(op.m_eOp == NodeOutputOp::MaxIndex) ||
		(op.m_eOp == NodeOutputOp::MinIndex)
	) {
		FindLocalMinMax<float>(
			grid,
			(op.m_eOp == NodeOutputOp::MinIndex),
			dataState,
			ixCandidate,
			op.m_dDistance,
			ixExtremum,
			dValue,
			dRMax);

		snprintf(buf, 100, "%i", ixExtremum);
		strResult = buf;

	// Average of the field over a given distance
	} else if (op.m_eOp == NodeOutputOp::Avg) {
		FindLocalAverage<float>(
			grid,
			dataState,
			ixCandidate,
			op.m_dDistance,
			dValue);

		snprintf(buf, 100, szFormat, dValue);
		strResult = buf;

	// Positive closed contour deltas (valleys)
	} else if (op.m_eOp == NodeOutputOp::PosClosedContour) {
		FindMaxClosedContourDelta<float>(
			grid,
			dataState,
			ixCandidate,
			op.m_dDistance,
			op.m_dMinMaxDist,
			true,
			dValue);

		snprintf(buf, 100, szFormat, dValue);
		strResult = buf;

	// Negative closed contour deltas (hills)
	} else if (op.m_eOp == NodeOutputOp::NegClosedContour) {
		FindMaxClosedContourDelta<float>(
			grid,
			dataState,
			ixCandidate,
			op.m_dDistance,
			op.m_dMinMaxDist,
			false,
			dValue);

		snprintf(buf, 100, szFormat, dValue);
		strResult = buf;

	// Positive minus negative weighted area
	} else if (op.m_eOp == NodeOutputOp::PosMinusNegWtArea) {
		PositiveMinusNegativeWeightedArea<float>(
			grid,
			dataState,
			ixCandidate,
			op.m_dDistance,
			dValue);

		snprintf(buf, 100, szFormat, dValue);
		strResult = buf;

	// Positive minus negative weighted area
	} else if (op.m_eOp == NodeOutputOp::MaxPoleward) {
		MaxPolewardValue<float>(
			grid,
			dataState,
			ixCandidate,
			op.m_dDistance,
			dValue);

		snprintf(buf, 100, szFormat, dValue);
		strResult = buf;

	} else {
		_EXCEPTIONT("Invalid Output operator");
	}
}

///////////////////////////////////////////////////////////////////////////////

#endif // _NODEOUTPUTOP_H_

