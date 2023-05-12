///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectBlobs.cpp
///	\author  Paul Ullrich
///	\version July 17th, 2018
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
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

#include "Variable.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "SimpleGrid.h"
#include "BlobUtilities.h"
#include "CoordTransforms.h"
#include "Units.h"
#include "TimeMatch.h"

#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include "ThresholdOp.h"

#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, int> Point;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field satisfies the threshold.
///	</summary>
template <typename real>
bool SatisfiesThreshold(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	const int ix0,
	const ThresholdOp::Operation op,
	const double dTargetValue,
	const double dMaxDist
) {
	// Special case if dMaxDist is zero
	if (dMaxDist < 1.0e-12) {
		double dValue = dataState[ix0];

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

	}

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

// Set of indicator locations stored as grid indices
typedef std::set<int> IndicatorSet;
typedef IndicatorSet::iterator IndicatorSetIterator;
typedef IndicatorSet::const_iterator IndicatorSetConstIterator;

///////////////////////////////////////////////////////////////////////////////

class GeoFilterOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum GeoFilterProperty {
		Area,
		ArealFraction
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		const std::string & strOp
	) {
		m_strThis = strOp;

		// Read mode
		enum {
			ReadMode_Property,
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Invalid
		} eReadMode = ReadMode_Property;

		std::string strAreaString;

		// Loop through string
		int iLast = 0;
		for (int i = 0; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in property
				if (eReadMode == ReadMode_Property) {
					if (strSubStr == "area") {
						m_eProperty = Area;
						eReadMode = ReadMode_Op;

					} else if (strSubStr == "areafrac") {
						m_eProperty = ArealFraction;
						eReadMode = ReadMode_Op;

					} else {
						_EXCEPTION1("--geofiltercmd invalid property \"%s\", expected \"area\" or \"areafrac\".",
							strSubStr.c_str());
					}

					iLast = i + 1;

				// Read in operation
				} else if (eReadMode == ReadMode_Op) {
					if (strSubStr == ">") {
						m_eOp = ThresholdOp::GreaterThan;
					} else if (strSubStr == "<") {
						m_eOp = ThresholdOp::LessThan;
					} else if (strSubStr == ">=") {
						m_eOp = ThresholdOp::GreaterThanEqualTo;
					} else if (strSubStr == "<=") {
						m_eOp = ThresholdOp::LessThanEqualTo;
					} else if (strSubStr == "=") {
						m_eOp = ThresholdOp::EqualTo;
					} else if (strSubStr == "!=") {
						m_eOp = ThresholdOp::NotEqualTo;
					} else {
						_EXCEPTION1("--geofiltercmd invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Value;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					strAreaString = strSubStr;

					std::string strValue;
					std::string strUnits;

					SplitIntoValueAndUnits(strAreaString, strValue, strUnits);
					if (!STLStringHelper::IsFloat(strValue)) {
						_EXCEPTIONT("Value entry in --geofiltercmd must be of floating point type");
					}

					m_dValue = std::atof(strValue.c_str());

					if (m_eProperty == Area) {
						if (strUnits == "") {
							strUnits = "sr";
						}

						bool fIsArea = ConvertUnits<double>(m_dValue, strUnits, std::string("sr"), false);
						if (!fIsArea) {
							_EXCEPTION1("Value entry in --geofiltercmd \"%s\" unknown units; expected \"sr\", \"m2\", \"km2\"",
								strUnits.c_str());
						}

					} else if (m_eProperty == ArealFraction) {
						if (strUnits == "") {
							strUnits = "%";
						}
						if (strUnits != "%") {
							_EXCEPTION1("Value entry in --geofiltercmd \"%s\" unknown units; expected \"%\"",
								strUnits.c_str());
						}
						strAreaString = strValue + "%";
					}

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("Too many entries in --geofiltercmd string \"%s\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("Insufficient entries in --geofiltercmd string \"%s\"",
					strOp.c_str());
		}

		// Check values
		if (m_eProperty == Area) {
			if (m_dValue < 0.0) {
				_EXCEPTIONT("Area threshold must be nonnegative in --geofiltercmd");
			}
		}
		if (m_eProperty == ArealFraction) {
			if ((m_dValue < 0.0) || (m_dValue > 100.0)) {
				_EXCEPTIONT("Areal fraction threshold (%%) must be between 0 and 100 in --geofiltercmd");
			}
		}

		// Output announcement
		std::string strDescription;
		if (m_eProperty == Area) {
			strDescription = "Area";
		} else {
			strDescription = "Areal fraction";
		}

		if (m_eOp == ThresholdOp::GreaterThan) {
			strDescription += " is greater than ";
		} else if (m_eOp == ThresholdOp::LessThan) {
			strDescription += " is less than ";
		} else if (m_eOp == ThresholdOp::GreaterThanEqualTo) {
			strDescription += " is greater than or equal to ";
		} else if (m_eOp == ThresholdOp::LessThanEqualTo) {
			strDescription += " is less than or equal to ";
		} else if (m_eOp == ThresholdOp::EqualTo) {
			strDescription += " is equal to ";
		} else if (m_eOp == ThresholdOp::NotEqualTo) {
			strDescription += " is not equal to ";
		}

		strDescription += strAreaString;

		Announce("%s", strDescription.c_str());
	}

	///	<summary>
	///		Verify that the specified blob satisfies the geofilter op.
	///	</summary>
	bool Apply(
		const SimpleGrid & grid,
		const IndicatorSet & setBlobPoints,
		const LatLonBox<double> & boxBlobRad
	) const {
		_ASSERT(
			(m_eOp == ThresholdOp::GreaterThan) ||
			(m_eOp == ThresholdOp::LessThan) ||
			(m_eOp == ThresholdOp::GreaterThanEqualTo) ||
			(m_eOp == ThresholdOp::LessThanEqualTo) ||
			(m_eOp == ThresholdOp::EqualTo) ||
			(m_eOp == ThresholdOp::NotEqualTo));

		// Thresholds related to area
		if ((m_eProperty == Area) ||
		    (m_eProperty == ArealFraction)
		) {
			if (grid.m_dArea.GetRows() == 0) {
				_EXCEPTIONT("Face area, which is needed for GeoFilter operation, is not defined for this grid.");
			}

			// Calculate the area of each blob
			double dBlobArea = 0.0;
			auto iterBlob = setBlobPoints.begin();
			for (; iterBlob != setBlobPoints.end(); iterBlob++) {
				dBlobArea += grid.m_dArea[*iterBlob];
			}

			// Minimum area
			if (m_eProperty == Area) {
				if ((m_eOp == ThresholdOp::GreaterThan) && (dBlobArea > m_dValue)) {
					return true;

				} else if ((m_eOp == ThresholdOp::LessThan) && (dBlobArea < m_dValue)) {
					return true;

				} else if ((m_eOp == ThresholdOp::GreaterThanEqualTo) && (dBlobArea >= m_dValue)) {
					return true;

				} else if ((m_eOp == ThresholdOp::LessThanEqualTo) && (dBlobArea <= m_dValue)) {
					return true;

				} else if ((m_eOp == ThresholdOp::EqualTo) && (dBlobArea == m_dValue)) {
					return true;

				} else if ((m_eOp == ThresholdOp::NotEqualTo) && (dBlobArea != m_dValue)) {
					return true;
				}

			} else if (m_eProperty == ArealFraction) {

				// Calculate the area of the blob box
				double dBoxArea = fabs(sin(boxBlobRad.lat[1]) - sin(boxBlobRad.lat[0]));
				if (boxBlobRad.lon[1] > boxBlobRad.lon[0]) {
					dBoxArea *= (boxBlobRad.lon[1] - boxBlobRad.lon[0]);
				} else {
					dBoxArea *= (2.0 * M_PI - boxBlobRad.lon[1] + boxBlobRad.lon[0]);
				}

				double dBoxAreaThresh = dBoxArea * m_dValue / 100.0;

				if ((m_eOp == ThresholdOp::GreaterThan) && (dBlobArea > dBoxAreaThresh)) {
					return true;

				} else if ((m_eOp == ThresholdOp::LessThan) && (dBlobArea < dBoxAreaThresh)) {
					return true;

				} else if ((m_eOp == ThresholdOp::GreaterThanEqualTo) && (dBlobArea >= dBoxAreaThresh)) {
					return true;

				} else if ((m_eOp == ThresholdOp::LessThanEqualTo) && (dBlobArea <= dBoxAreaThresh)) {
					return true;

				} else if ((m_eOp == ThresholdOp::EqualTo) && (dBlobArea == dBoxAreaThresh)) {
					return true;

				} else if ((m_eOp == ThresholdOp::NotEqualTo) && (dBlobArea != dBoxAreaThresh)) {
					return true;
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
		// Invalid property
		} else {
			_EXCEPTIONT("Invalid property");
		}

		return false;
	}

	///	<summary>
	///		Get a string representation of this operator.
	///	</summary>
	const std::string & ToOperatorString() const {
		return m_strThis;
	}

protected:
	///	<summary>
	///		A string representation of this operator.
	///	</summary>
	std::string m_strThis;

	///	<summary>
	///		Threshold quantity.
	///	</summary>
	GeoFilterProperty m_eProperty;

	///	<summary>
	///		Operation that must be satisfied.
	///	</summary>
	ThresholdOp::Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a filtering operator.
///	</summary>
class FilterOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FilterOp() :
		m_varix(InvalidVariableIndex),
		m_eOp(ThresholdOp::GreaterThan),
		m_dValue(0.0),
		m_nCount(0)
	{ }

public:
	///	<summary>
	///		Parse a filter operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Count,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
					if (strSubStr == ">") {
						m_eOp = ThresholdOp::GreaterThan;
					} else if (strSubStr == "<") {
						m_eOp = ThresholdOp::LessThan;
					} else if (strSubStr == ">=") {
						m_eOp = ThresholdOp::GreaterThanEqualTo;
					} else if (strSubStr == "<=") {
						m_eOp = ThresholdOp::LessThanEqualTo;
					} else if (strSubStr == "=") {
						m_eOp = ThresholdOp::EqualTo;
					} else if (strSubStr == "!=") {
						m_eOp = ThresholdOp::NotEqualTo;
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
					eReadMode = ReadMode_Count;

				// Read in minimum count
				} else if (eReadMode == ReadMode_Count) {
					m_nCount = atoi(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in filter op \"%s\""
							"\nRequired: \"<name>,<operation>"
							",<value>,<count>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in filter op \"%s\""
					"\nRequired: \"<name>,<operation>,<value>,<count>\"",
					strOp.c_str());
		}

		if (m_nCount < 1) {
			_EXCEPTIONT("For filter op, count must be positive");
		}

		// Output announcement
		std::string strDescription = varreg.GetVariableString(m_varix);

		if (m_eOp == ThresholdOp::GreaterThan) {
			strDescription += " is greater than ";
		} else if (m_eOp == ThresholdOp::LessThan) {
			strDescription += " is less than ";
		} else if (m_eOp == ThresholdOp::GreaterThanEqualTo) {
			strDescription += " is greater than or equal to ";
		} else if (m_eOp == ThresholdOp::LessThanEqualTo) {
			strDescription += " is less than or equal to ";
		} else if (m_eOp == ThresholdOp::EqualTo) {
			strDescription += " is equal to ";
		} else if (m_eOp == ThresholdOp::NotEqualTo) {
			strDescription += " is not equal to ";
		}

		Announce("%s %f at at least %i points",
			strDescription.c_str(),
			m_dValue,
			m_nCount);
	}

public:
	///	<summary>
	///		Variable to use for filtering.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Operation that must be satisfied.
	///	</summary>
	ThresholdOp::Operation m_eOp;

	///	<summary>
	///		Filter threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Number of points that exceed the threshold.
	///	</summary>
	int m_nCount;
};

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void ApplyFilters(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	const ThresholdOp::Operation op,
	const double dTargetValue,
	const int nCount,
	DataArray1D<int> & bTag
) {
	_ASSERT(bTag.GetRows() == grid.GetSize());
	_ASSERT(bTag.GetRows() == dataState.GetRows());

	// Number of blobs
	int nBlobs = 0;
	int nBlobsFiltered = 0;

	// Visit all tagged points
	std::set<int> setNodesVisited;
	for (int i = 0; i < grid.GetSize(); i++) {

		// Verify point it tagged and we haven't visited before
		if (bTag[i] == 0) {
			continue;
		}
		if (setNodesVisited.find(i) != setNodesVisited.end()) {
			continue;
		}

		// Number of blobs
		nBlobs++;

		// New blob
		IndicatorSet setCurrentBlob;

		// Number of points within blob that satisfy threshold
		int nThresholdPoints = 0;

		// Build connectivity
		std::queue<int> queueToVisit;
		queueToVisit.push(i);

		while (queueToVisit.size() != 0) {
			int iNext = queueToVisit.front();
			queueToVisit.pop();

			// Verify point is tagged and we haven't visited before
			if (bTag[iNext] == 0) {
				continue;
			}
			if (setNodesVisited.find(iNext) != setNodesVisited.end()) {
				continue;
			}
			setNodesVisited.insert(iNext);
			setCurrentBlob.insert(iNext);

			// Check if this point satisfies threshold
			bool fSatisfiesThreshold =
				SatisfiesThreshold<real>(
					grid,
					dataState,
					iNext,
					op,
					dTargetValue,
					0.0);

			if (fSatisfiesThreshold) {
				nThresholdPoints++;
			}

			// Insert all connected neighbors into "to visit" queue
			for (int n = 0; n < grid.m_vecConnectivity[iNext].size(); n++) {
				queueToVisit.push(grid.m_vecConnectivity[iNext][n]);
			}
		}

		// If not enough points satisfy the filter then eliminate this blob
		if (nThresholdPoints < nCount) {
			nBlobsFiltered++;
			for (auto it = setCurrentBlob.begin(); it != setCurrentBlob.end(); it++) {
				bTag[*it] = 0;
			}
		}
	}

	// Announce results
	Announce("Filter removed %i of %i blobs", nBlobsFiltered, nBlobs);
}

///////////////////////////////////////////////////////////////////////////////

void ApplyGeoFilter(
	const SimpleGrid & grid,
	bool fRegional,
	const GeoFilterOp & op,
	DataArray1D<int> & bTag
) {
	_ASSERT(bTag.GetRows() == grid.GetSize());

	// Number of blobs
	int nBlobs = 0;
	int nBlobsFiltered = 0;

	// Visit all tagged points
	std::set<int> setNodesVisited;
	for (int i = 0; i < grid.GetSize(); i++) {

		// Verify point it tagged and we haven't visited before
		if (bTag[i] == 0) {
			continue;
		}
		if (setNodesVisited.find(i) != setNodesVisited.end()) {
			continue;
		}

		// Number of blobs
		nBlobs++;

		// New blob
		IndicatorSet setCurrentBlob;
		LatLonBox<double> boxCurrentBlob(!fRegional, 2.0 * M_PI);

		// Build connectivity
		std::queue<int> queueToVisit;
		queueToVisit.push(i);

		while (queueToVisit.size() != 0) {
			int iNext = queueToVisit.front();
			queueToVisit.pop();

			// Verify point is tagged and we haven't visited before
			if (bTag[iNext] == 0) {
				continue;
			}
			if (setNodesVisited.find(iNext) != setNodesVisited.end()) {
				continue;
			}
			setNodesVisited.insert(iNext);
			setCurrentBlob.insert(iNext);

			//std::cout << LonRadToStandardRange(grid.m_dLon[iNext]) << std::endl;
			boxCurrentBlob.insert(
				grid.m_dLat[iNext],
				LonRadToStandardRange(grid.m_dLon[iNext]));

			// Insert all connected neighbors into "to visit" queue
			for (int n = 0; n < grid.m_vecConnectivity[iNext].size(); n++) {
				queueToVisit.push(grid.m_vecConnectivity[iNext][n]);
			}
		}

		bool fSuccess = op.Apply(grid, setCurrentBlob, boxCurrentBlob);

		// If not enough points satisfy the filter then eliminate this blob
		if (!fSuccess) {
			nBlobsFiltered++;
			for (auto it = setCurrentBlob.begin(); it != setCurrentBlob.end(); it++) {
				bTag[*it] = 0;
			}
		}
	}

	// Announce results
	Announce("Geofilter \"%s\" removed %i of %i blobs", op.ToOperatorString().c_str(), nBlobsFiltered, nBlobs);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a output operator.
///	</summary>
class BlobOutputOp {

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
			ReadMode_Name,
			ReadMode_Invalid
		} eReadMode = ReadMode_Name;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in name
				if (eReadMode == ReadMode_Name) {
					m_strName = strSubStr;
					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in output op \"%s\""
							"\nRequired: \"<variable>,<name>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in output op \"%s\""
					"\nRequired: \"<variable>,<name>\"",
					strOp.c_str());
		}

		// Output announcement
		std::string strDescription = varreg.GetVariableString(m_varix);
		strDescription += " with name \"";
		strDescription += m_strName;
		strDescription += "\"";

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for output.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Name of variable in NetCDF output file.
	///	</summary>
	std::string m_strName;
};

///////////////////////////////////////////////////////////////////////////////

class DetectBlobsParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DetectBlobsParam() :
		strGridFile(),
		dMinAbsLat(0.0),
		dMaxAbsLat(90.0),
		dMinLat(-90.0),
		dMaxLat(90.0),
		dMinLon(0.0),
		dMaxLon(0.0),
		fOutFloat(false),
		fRegional(false),
		fDiagonalConnectivity(false),
		iVerbosityLevel(0),
		strTagVar("binary_tag"),
		strLongitudeName("lon"),
		strLatitudeName("lat"),
		pvecThresholdOp(NULL),
		pvecFilterOp(NULL),
		pvecGeoFilterOp(NULL),
		pvecOutputOp(NULL),
		strTimeFilter("")
	{ }

public:
	// Grid file
	std::string strGridFile;

	// Minimum absolute latitude (in degrees)
	double dMinAbsLat;

	// Minimum absolute latitude (in degrees)
	double dMaxAbsLat;

	// Minimum latitude (in degrees)
	double dMinLat;

	// Maximum latitude (in degrees)
	double dMaxLat;

	// Minimum longitude (in degrees)
	double dMinLon;

	// Maximum latitude (in degrees)
	double dMaxLon;

	// Write output as float
	bool fOutFloat;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Verbosity level
	int iVerbosityLevel;

	// Name of output variable for tag
	std::string strTagVar;

	// Name of longitude variabe
	std::string strLongitudeName;

	// Name of latitude variable
	std::string strLatitudeName;

	// Vector of threshold operators
	std::vector<ThresholdOp> * pvecThresholdOp;

	// Vector of filter operators
	std::vector<FilterOp> * pvecFilterOp;

	// Vector of geometric filter operators
	std::vector<GeoFilterOp> * pvecGeoFilterOp;

	// Vector of output operators
	std::vector<BlobOutputOp> * pvecOutputOp;

	// Time filter
	std::string strTimeFilter;
};

///////////////////////////////////////////////////////////////////////////////

void DetectBlobs(
	int iFile,
	const std::string & strInputFiles,
	const std::string & strOutputFile,
	const std::string & strConnectivity,
	VariableRegistry & varreg,
	const DetectBlobsParam & param
) {
	// Dereference pointers to operators
	_ASSERT(param.pvecThresholdOp != NULL);
	std::vector<ThresholdOp> & vecThresholdOp =
		*(param.pvecThresholdOp);

	// Dereference pointers to operators
	_ASSERT(param.pvecFilterOp != NULL);
	std::vector<FilterOp> & vecFilterOp =
		*(param.pvecFilterOp);

	// Dereference pointers to operators
	_ASSERT(param.pvecGeoFilterOp != NULL);
	std::vector<GeoFilterOp> & vecGeoFilterOp =
		*(param.pvecGeoFilterOp);

#ifdef TEMPEST_NOREGEX
	if (param.strTimeFilter != "") {
		_EXCEPTIONT("Cannot use --timefilter with -DTEMPEST_NOREGEX compiler flag");
	}
#endif
#ifndef TEMPEST_NOREGEX
	// Parse --timefilter
	std::regex reTimeSubset;
	if (param.strTimeFilter != "") {
		
		// Test regex support
		TestRegex();

		std::string strTimeFilter = param.strTimeFilter;
		if (strTimeFilter == "3hr") {
			strTimeFilter = "....-..-.. (00|03|06|09|12|15|18|21):00:00";
		}
		if (strTimeFilter == "6hr") {
			strTimeFilter = "....-..-.. (00|06|12|18):00:00";
		}
		if (strTimeFilter == "daily") {
			strTimeFilter = "....-..-.. 00:00:00";
		}

		try {
			reTimeSubset.assign(strTimeFilter);
		} catch(std::regex_error & reerr) {
			_EXCEPTION2("Parse error in --timefilter regular expression \"%s\" (code %i)",
				strTimeFilter.c_str(), reerr.code());
		}
	}
#endif

	// Unload data from the VariableRegistry
	varreg.UnloadAllGridData();

	// Define the SimpleGrid
	SimpleGrid grid;

	// Load in the benchmark file
	NcFileVector vecFiles;
	vecFiles.ParseFromString(strInputFiles);

	// Latitude/longitude names
	std::string strLatitudeName(param.strLatitudeName);
	std::string strLongitudeName(param.strLongitudeName);

	NcFile * pncGridFile = NULL;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for grid file
	} else if (param.strGridFile != "") {
		AnnounceStartBlock("Generating grid information from grid file");
		pncGridFile = new NcFile(param.strGridFile.c_str());
		if ((pncGridFile == NULL) || (!pncGridFile->is_valid())) {
			_EXCEPTION1("Unable to open grid file \"%s\"", param.strGridFile.c_str());
		}

		grid.GenerateLatitudeLongitude(
			pncGridFile,
			strLatitudeName,
			strLongitudeName,
			param.fRegional,
			param.fDiagonalConnectivity);
		AnnounceEndBlock("Done");

	// Try generating grid information from data file
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		grid.GenerateLatitudeLongitude(
			vecFiles[0],
			strLatitudeName,
			strLongitudeName,
			param.fRegional,
			param.fDiagonalConnectivity);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Get time dimension
	NcVar * varTime = vecFiles[0]->get_var("time");
	if (varTime == NULL) {
		_EXCEPTION1("File \"%s\" missing \"time\" variable",
			vecFiles.GetFilename(0).c_str());
	}
	NcDim * dimTime = vecFiles[0]->get_dim("time");
	if (dimTime == NULL) {
		if (varTime->num_dims() != 0) {
			_EXCEPTION1("File \"%s\" missing \"time\" dimension",
				vecFiles.GetFilename(0).c_str());
		}
	}

	// Read the time data
	const NcTimeDimension & vecTimes = vecFiles.GetNcTimeDimension(0);

	std::vector<bool> vecTimeRetained;
	std::vector<Time> vecOutputTimes;
#ifndef TEMPEST_NOREGEX
	if (param.strTimeFilter != "") {
		vecTimeRetained.resize(vecTimes.size(), false);
		for (int t = 0; t < vecTimes.size(); t++) {
			std::string strTime = vecTimes[t].ToString();
			std::smatch match;
			if (std::regex_search(strTime, match, reTimeSubset)) {
				vecOutputTimes.push_back(vecTimes[t]);
				vecTimeRetained[t] = true;
			}
		}

	} else {
		vecOutputTimes = vecTimes;
		vecTimeRetained.resize(vecTimes.size(), true);
	}
#else
	vecOutputTimes = vecTimes;
	vecTimeRetained.resize(vecTimes.size(), true);
#endif

	// Create reference to NetCDF input file
	NcFile & ncInput = *(vecFiles[0]);

	// Open the NetCDF output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Copy over time variables to output file
	NcDim * dimTimeOut = NULL;
	if ((dimTime != NULL) && (varTime != NULL)) {
		if (param.strTimeFilter != "") {
			CopyNcVarTimeSubset(ncInput, ncOutput, "time", vecOutputTimes);

		} else {
			CopyNcVar(ncInput, ncOutput, "time", true);
		}

		dimTimeOut = ncOutput.get_dim("time");
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

	} /*else if (nAddTimeDim != -1) {
		dimTimeOut = ncOutput.add_dim("time", 0);
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error creating dimension \"time\" in output file");
		}

		NcVar * varTimeOut = ncOutput.add_var("time", ncDouble, dimTimeOut);
		if (varTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

		double dTime = static_cast<double>(nAddTimeDim);
		varTimeOut->set_cur((long)0);
		varTimeOut->put(&dTime, 1);

		if (strAddTimeDimUnits != "") {
			varTimeOut->add_att("units", strAddTimeDimUnits.c_str());
		}
		varTimeOut->add_att("long_name", "time");
		varTimeOut->add_att("calendar", "standard");
		varTimeOut->add_att("standard_name", "time");
	}*/

	// Create output variable
	NcDim * dim0 = NULL;
	NcDim * dim1 = NULL;
	NcVar * varTag = NULL;

	PrepareBlobOutputVar(
		(pncGridFile == NULL)?(ncInput):(*pncGridFile),
		ncOutput,
		strOutputFile,
		grid,
		param.strTagVar,
		strLatitudeName,
		strLongitudeName,
		(param.fOutFloat)?(ncFloat):(ncByte),
		dimTimeOut,
		&dim0,
		&dim1,
		&varTag);

	_ASSERT(varTag != NULL);

	if (pncGridFile != NULL) {
		delete pncGridFile;
	}

/*
	CopyNcVarIfExists(ncInput, ncOutput, param.strLatitudeName, true);
	CopyNcVarIfExists(ncInput, ncOutput, param.strLongitudeName, true);


	if (grid.m_nGridDim.size() == 1) {
		dim0 = ncOutput.get_dim("ncol");
		if (dim0 == NULL) {
			dim0 = ncOutput.add_dim("ncol", grid.GetSize());
		}
		if (dim0 == NULL) {
			_EXCEPTION1("Error creating dim \"ncol\" in file \"%s\"",
				strOutputFile.c_str());
		}
		if (dim0->size() != grid.GetSize()) {
			_EXCEPTION4("Dimension \"%s\" in file \"%s\" has inconsistent length (%i vs %i)",
				dim0->name(),
				strOutputFile.c_str(),
				dim0->size(),
				grid.GetSize());
		}

		if (dimTimeOut != NULL) {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dimTimeOut,
				dim0);

		} else {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dim0);
		}

	} else if (grid.m_nGridDim.size() == 2) {

		// Copy over latitude and longitude variable
		dim0 = ncOutput.get_dim(param.strLatitudeName.c_str());
		if (dim0 == NULL) {
			_EXCEPTION1("Error copying variable \"%s\" to output file",
				param.strLatitudeName.c_str());
		}
		dim1 = ncOutput.get_dim(param.strLongitudeName.c_str());
		if (dim1 == NULL) {
			_EXCEPTION1("Error copying variable \"%s\" to output file",
				param.strLongitudeName.c_str());
		}

		// Create output tag
		if (dimTimeOut != NULL) {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dimTimeOut,
				dim0,
				dim1);

		} else {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dim0,
				dim1);
		}

	} else {
		_EXCEPTIONT("Invalid grid dimension -- value must be 1 or 2");
	}

	_ASSERT(varTag != NULL);
*/
	AnnounceEndBlock("Done");

	// Tagged cell array
	DataArray1D<int> bTag(grid.GetSize());

	// Loop through all times
	int to = (-1);
	for (int t = 0; t < vecTimes.size(); t++) {

		// Announce
		std::string strTime = vecTimes[t].ToString();
		AnnounceStartBlock("Time %s", strTime.c_str());

		if (!vecTimeRetained[t]) {
			AnnounceEndBlock("Skipping (timefilter)");
			continue;
		}
		to++;

		AnnounceStartBlock("Build tagged cell array");
		bTag.Zero();

		// Set all points within the specified latitude/longitude bounds to 1
		for (int i = 0; i < grid.GetSize(); i++) {
			if (fabs(grid.m_dLat[i]) < DegToRad(param.dMinAbsLat)) {
				continue;
			}
			if (fabs(grid.m_dLat[i]) > DegToRad(param.dMaxAbsLat)) {
				continue;
			}
			if (grid.m_dLat[i] < DegToRad(param.dMinLat)) {
				continue;
			}
			if (grid.m_dLat[i] > DegToRad(param.dMaxLat)) {
				continue;
			}
			if (param.dMinLon != param.dMaxLon) {
				double dLon = LonDegToStandardRange(RadToDeg(grid.m_dLon[i]));
				double dMinLon = LonDegToStandardRange(param.dMinLon);
				double dMaxLon = LonDegToStandardRange(param.dMaxLon);

				if (dMinLon > dMaxLon) {
					if ((dLon < dMinLon) && (dLon > dMaxLon)) {
						continue;
					}

				} else {
					if ((dLon < dMinLon) || (dLon > dMaxLon)) {
						continue;
					}
				}
			}

			bTag[i] = 1;
		}
		AnnounceEndBlock("Done");

		// Eliminate based on threshold commands
		AnnounceStartBlock("Apply threshold commands");
		for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

			// Load the search variable data
			Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);
			vecFiles.SetTime(vecTimes[t]);
			var.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataState = var.GetData();

			// Loop through data
			for (int i = 0; i < grid.GetSize(); i++) {
				if (bTag[i] == 0) {
					continue;
				}

				// Determine if the threshold is satisfied
				bool fSatisfiesThreshold =
					SatisfiesThreshold<float>(
						grid,
						dataState,
						i,
						vecThresholdOp[tc].m_eOp,
						vecThresholdOp[tc].m_dValue,
						vecThresholdOp[tc].m_dDistance
					);

				if (!fSatisfiesThreshold) {
					bTag[i] = 0;
				}
			}
		}
		AnnounceEndBlock("Done");

		// Eliminate based on filter commands
		if (vecFilterOp.size() != 0) {
			AnnounceStartBlock("Apply filters");
			for (int fc = 0; fc < vecFilterOp.size(); fc++) {

				// Load the search variable data
				Variable & var = varreg.Get(vecFilterOp[fc].m_varix);
				vecFiles.SetTime(vecTimes[t]);
				var.LoadGridData(varreg, vecFiles, grid);
				const DataArray1D<float> & dataState = var.GetData();

				// Check the blobs against the filter
				ApplyFilters<float>(
					grid,
					dataState,
					vecFilterOp[fc].m_eOp,
					vecFilterOp[fc].m_dValue,
					vecFilterOp[fc].m_nCount,
					bTag);
			}
			AnnounceEndBlock("Done");
		}

		// Eliminate based on geofilter commands
		if (vecGeoFilterOp.size() != 0) {
			AnnounceStartBlock("Apply geometric filters");
			for (int gc = 0; gc < vecGeoFilterOp.size(); gc++) {
				ApplyGeoFilter(
					grid,
					param.fRegional,
					vecGeoFilterOp[gc],
					bTag);
			}
			AnnounceEndBlock("Done");
		}

		// Output tagged cell array
		AnnounceStartBlock("Writing results");
		if (dimTimeOut != NULL) {
			if (grid.m_nGridDim.size() == 1) {
				varTag->set_cur(to, 0);
				varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0]);

			} else if (grid.m_nGridDim.size() == 2) {
				varTag->set_cur(to, 0, 0);
				varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);

			} else {
				_EXCEPTION();
			}

			for (int oc = 0; oc < param.pvecOutputOp->size(); oc++) {
				if (grid.m_nGridDim.size() == 1) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dimTimeOut,
							dim0);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					vecFiles.SetTime(vecTimes[t]);
					var.LoadGridData(varreg, vecFiles, grid);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur(to, 0);
					ncvar->put(&(dataState[0]), 1, grid.m_nGridDim[0]);

				} else if (grid.m_nGridDim.size() == 2) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dimTimeOut,
							dim0,
							dim1);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					vecFiles.SetTime(vecTimes[t]);
					var.LoadGridData(varreg, vecFiles, grid);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur(to, 0, 0);
					ncvar->put(&(dataState[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);

				} else {
					_EXCEPTION();
				}

			}

		} else {
			if (grid.m_nGridDim.size() == 1) {
				varTag->set_cur((long)0);
				varTag->put(&(bTag[0]), grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varTag->set_cur(0, 0);
				varTag->put(&(bTag[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}

			for (int oc = 0; oc < param.pvecOutputOp->size(); oc++) {
				if (grid.m_nGridDim.size() == 1) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dim0);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					vecFiles.SetTime(vecTimes[t]);
					var.LoadGridData(varreg, vecFiles, grid);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur((long)0);
					ncvar->put(&(dataState[0]), grid.m_nGridDim[0]);

				} else if (grid.m_nGridDim.size() == 2) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dim0,
							dim1);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					vecFiles.SetTime(vecTimes[t]);
					var.LoadGridData(varreg, vecFiles, grid);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur(0, 0);
					ncvar->put(&(dataState[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);

				} else {
					_EXCEPTION();
				}
			}
		}

		AnnounceEndBlock("Done");

		AnnounceEndBlock(NULL);
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {
	// Parameters for DetectBlobs
	DetectBlobsParam dbparam;

	// Input dat file
	std::string strInputFile;

	// Input list file
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutputFile;

	// Output file list
	std::string strOutputFileList;

	// Threshold commands
	std::string strThresholdCmd;

	// Filter commands
	std::string strFilterCmd;

	// Filter commands
	std::string strGeoFilterCmd;

	// Output commands
	std::string strOutputCmd;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(dbparam.fDiagonalConnectivity, "diag_connect");
		CommandLineString(dbparam.strGridFile, "grid_file", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputFileList, "out_list", "");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		CommandLineStringD(strFilterCmd, "filtercmd", "", "[var,op,value,count]");
		CommandLineStringD(strGeoFilterCmd, "geofiltercmd", "", "[prop,op,value]");
		CommandLineStringD(strOutputCmd, "outputcmd", "", "[var,name;...]");
		CommandLineString(dbparam.strTimeFilter, "timefilter", "");
		CommandLineDouble(dbparam.dMinAbsLat, "minabslat", 0.0);
		CommandLineDouble(dbparam.dMaxAbsLat, "maxabslat", 90.0);
		CommandLineDouble(dbparam.dMinLat, "minlat", -90.0);
		CommandLineDouble(dbparam.dMaxLat, "maxlat", 90.0);
		CommandLineDouble(dbparam.dMinLon, "minlon", 0.0);
		CommandLineDouble(dbparam.dMaxLon, "maxlon", 0.0);
		CommandLineBool(dbparam.fRegional, "regional");
		CommandLineBool(dbparam.fOutFloat, "out_float");
		CommandLineString(dbparam.strTagVar, "tagvar", "binary_tag");
		CommandLineString(dbparam.strLongitudeName, "lonname", "lon");
		CommandLineString(dbparam.strLatitudeName, "latname", "lat");
		CommandLineInt(dbparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(dbparam.iVerbosityLevel);

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
	if ((strOutputFile.length() != 0) && (strOutputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out) or (--out_list)"
			" may be specified");
	}

	// Load input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		std::ifstream ifInputFileList(strInputFileList.c_str());
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

		std::ifstream ifOutputFileList(strOutputFileList.c_str());
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

	// Parse the threshold operator command string
	std::vector<ThresholdOp> vecThresholdOp;
	dbparam.pvecThresholdOp = &vecThresholdOp;

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

	// Parse the filter operator command string
	std::vector<FilterOp> vecFilterOp;
	dbparam.pvecFilterOp = &vecFilterOp;

	if (strFilterCmd != "") {
		AnnounceStartBlock("Parsing filter operations");

		int iLast = 0;
		for (int i = 0; i <= strFilterCmd.length(); i++) {

			if ((i == strFilterCmd.length()) ||
				(strFilterCmd[i] == ';') ||
				(strFilterCmd[i] == ':')
			) {
				std::string strSubStr =
					strFilterCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecFilterOp.size());
				vecFilterOp.resize(iNextOp + 1);
				vecFilterOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the geofilter operator command string
	std::vector<GeoFilterOp> vecGeoFilterOp;
	dbparam.pvecGeoFilterOp = &vecGeoFilterOp;

	if (strGeoFilterCmd != "") {
		AnnounceStartBlock("Parsing geofilter operations");

		int iLast = 0;
		for (int i = 0; i <= strGeoFilterCmd.length(); i++) {

			if ((i == strGeoFilterCmd.length()) ||
				(strGeoFilterCmd[i] == ';') ||
				(strGeoFilterCmd[i] == ':')
			) {
				std::string strSubStr =
					strGeoFilterCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecGeoFilterOp.size());
				vecGeoFilterOp.resize(iNextOp + 1);
				vecGeoFilterOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the output operator command string
	std::vector<BlobOutputOp> vecOutputOp;
	dbparam.pvecOutputOp = &vecOutputOp;

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
			Announce("Output will be written following --out_list");
		} else if (strOutputFile == "") {
			Announce("Output will be written to outXXXXXX.dat");
		} else {
			Announce("Output will be written to %sXXXXXX.dat",
				strOutputFile.c_str());
		}
#if defined(TEMPEST_MPIOMP)
		if (vecInputFiles.size() != 1) {
			Announce("Logs will be written to logXXXXXX.txt");
		}
#endif
	}

	// Loop over all files to be processed
	for (int f = 0; f < vecInputFiles.size(); f++) {
#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}

		FILE * fpLog = NULL;
#endif
		// Output file for this file
		std::string strOutputFileCurrent;

		// Log file needs closing
		bool fCloseLogFile = false;

		// Generate output file name
		if (vecInputFiles.size() == 1) {
			if (vecOutputFiles.size() == 1) {
				strOutputFileCurrent = vecOutputFiles[0];
			} else if (strOutputFile == "") {
				strOutputFileCurrent = "out.nc";
			} else {
				strOutputFileCurrent = strOutputFile;
			}

		} else {
			char szFileIndex[32];
			sprintf(szFileIndex, "%06i", f);

			if (vecOutputFiles.size() != 0) {
				strOutputFileCurrent = vecOutputFiles[f];
			} else {
				if (strOutputFile == "") {
					strOutputFileCurrent =
						"out" + std::string(szFileIndex) + ".nc";
				} else {
					strOutputFileCurrent =
						strOutputFile + std::string(szFileIndex) + ".nc";
				}
			}

			std::string strLogFile = "log" + std::string(szFileIndex) + ".txt";
#if defined(TEMPEST_MPIOMP)
			fpLog = fopen(strLogFile.c_str(), "w");
			if (fpLog == NULL) {
				_EXCEPTION1("Unable to open log file \"%s\" for writing",
					strLogFile.c_str());
			}

			AnnounceSetOutputBuffer(fpLog);
			AnnounceOutputOnAllRanks();
#endif
		}

		// Perform DetectBlobs
		DetectBlobs(
			f,
			vecInputFiles[f],
			strOutputFileCurrent,
			strConnectivity,
			varreg,
			dbparam);

#if defined(TEMPEST_MPIOMP)
		// Close the log file
		if (fpLog != NULL) {
			AnnounceSetOutputBuffer(stdout);
			AnnounceOnlyOutputOnRankZero();
			fclose(fpLog);
		}
#endif
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

///////////////////////////////////////////////////////////////////////////////

