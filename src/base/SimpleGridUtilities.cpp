///////////////////////////////////////////////////////////////////////////////
///
///	\file    SimpleGridUtilities.cpp
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

#include "SimpleGridUtilities.h"
#include "CoordTransforms.h"
#include "ThresholdOp.h"

#include <queue>

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void FindLocalMinMax(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataArray1D<real> & data,
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

template <typename real>
void FindAllLocalMinima(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
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

template <typename real>
void FindAllLocalMaxima(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
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

template <typename real>
void FindAllLocalMinMaxWithThreshold(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	bool fMinima,
	const std::string & strThreshold,
	std::set<int> & setMinima
) {
	ThresholdOp::Operation opThreshold = ThresholdOp::NoThreshold;
	real dThresholdValue;

	if (strThreshold.length() != 0) {
		if (strThreshold.length() < 2) {
			_EXCEPTION1("Threshold operator \"%s\" must be of form \"<op><value>\"", strThreshold.c_str());
		}

		std::string strThreshold1 = strThreshold.substr(0,1);
		std::string strThreshold2 = strThreshold.substr(0,2);
		std::string strThresholdValue;

		if (strThreshold2 == ">=") {
			opThreshold = ThresholdOp::GreaterThanEqualTo;
			strThresholdValue = strThreshold.substr(2);

		} else if (strThreshold2 == "<=") {
			opThreshold = ThresholdOp::LessThanEqualTo;
			strThresholdValue = strThreshold.substr(2);

		} else if (strThreshold2 == "!=") {
			opThreshold = ThresholdOp::NotEqualTo;
			strThresholdValue = strThreshold.substr(2);

		} else if (strThreshold1 == ">") {
			opThreshold = ThresholdOp::GreaterThan;
			strThresholdValue = strThreshold.substr(1);

		} else if (strThreshold1 == "<") {
			opThreshold = ThresholdOp::LessThan;
			strThresholdValue = strThreshold.substr(1);

		} else if (strThreshold1 == "=") {
			opThreshold = ThresholdOp::EqualTo;
			strThresholdValue = strThreshold.substr(1);

		} else {
			_EXCEPTION1("Threshold operator \"%s\" must be of form \"<op><value>\" where <op> is one of >=,<=,!=,>,<,=",
				strThreshold.c_str());
		}

		// Extract value
		if (!STLStringHelper::IsFloat(strThresholdValue)) {
			_EXCEPTION1("Threshold operator \"%s\" must be of form \"<op><value>\" where <value> is a floating point number",
				strThreshold.c_str());
		}
		dThresholdValue = static_cast<real>(std::stod(strThresholdValue));
	}

	real dSign = (fMinima)?(-1.0):(1.0);

	int sFaces = grid.m_vecConnectivity.size();
	for (int f = 0; f < sFaces; f++) {
		
		real dValue = data[f];

		// Check thresholds
		if (opThreshold != ThresholdOp::NoThreshold) {
			if (opThreshold == ThresholdOp::GreaterThan) {
				if (dValue <= static_cast<real>(dThresholdValue)) {
					continue;
				}

			} else if (opThreshold == ThresholdOp::LessThan) {
				if (dValue >= static_cast<real>(dThresholdValue)) {
					continue;
				}

			} else if (opThreshold == ThresholdOp::GreaterThanEqualTo) {
				if (dValue < static_cast<real>(dThresholdValue)) {
					continue;
				}

			} else if (opThreshold == ThresholdOp::LessThanEqualTo) {
				if (dValue > static_cast<real>(dThresholdValue)) {
					continue;
				}

			} else if (opThreshold == ThresholdOp::EqualTo) {
				if (dValue != static_cast<real>(dThresholdValue)) {
					continue;
				}

			} else if (opThreshold == ThresholdOp::NotEqualTo) {
				if (dValue == static_cast<real>(dThresholdValue)) {
					continue;
				}

			} else {
				_EXCEPTIONT("Invalid operator");
			}
		}

		// Check neighbors
		bool fExtrema = true;
		int sNeighbors = grid.m_vecConnectivity[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (dSign * data[grid.m_vecConnectivity[f][n]] > dSign * dValue) {
				fExtrema = false;
				break;
			}
		}

		if (fExtrema) {
			setMinima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void FindAllLocalMinMaxWithGraphDistance(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	bool fMinima,
	int nMaxGraphDistance,
	std::set<int> & setMinMax
) {
	const real dSign = (fMinima)?(-1.0):(1.0);

	// If max graph distance is 1 then this is just FindAllLocalMinima/Maxima
	if (nMaxGraphDistance <= 1) {
		if (fMinima) {
			FindAllLocalMinima<real>(grid, data, setMinMax);
		} else {
			FindAllLocalMaxima<real>(grid, data, setMinMax);
		}
	}

	// Get first pass set
	std::set<int> setFirstPass;
	if (fMinima) {
		FindAllLocalMinima<real>(grid, data, setFirstPass);
	} else {
		FindAllLocalMaxima<real>(grid, data, setFirstPass);
	}

	// Use graph search to eliminate candidates within specified graph distance
	for (auto it = setFirstPass.begin(); it != setFirstPass.end(); it++) {

		real dValue = data[*it];

		std::queue<int> queueToVisit;
		std::set<int> setVisited;

		// Insert this node and its neighbors into the "visited" set
		setVisited.insert(*it);
		size_t sNeighbors = grid.m_vecConnectivity[*it].size();
		for (size_t n = 0; n < sNeighbors; n++) {
			int ixNeighbor = grid.m_vecConnectivity[*it][n];
			setVisited.insert(ixNeighbor);
		}

		// Insert this node's neighbor's neighbors into the "tovisit" set
		for (size_t n = 0; n < sNeighbors; n++) {
			int ixNeighbor = grid.m_vecConnectivity[*it][n];
			size_t sSubNeighbors = grid.m_vecConnectivity[ixNeighbor].size();
			for (size_t m = 0; m < sSubNeighbors; m++) {
				int ixSubNeighbor = grid.m_vecConnectivity[*it][n];
				if (setVisited.find(ixSubNeighbor) != setVisited.end()) {
					queueToVisit.push(ixSubNeighbor);
				}
			}
		}

		// Perform breadth-first search
		bool fExtrema = true;

		int nCurrentDistance = 2;
		int nNodesVisitedAtDist = 0;
		int nCurrentDistanceSize = queueToVisit.size();
		if (nCurrentDistanceSize == 0) {
			setMinMax.insert(*it);
			continue;
		}

		for (;;) {
			if (nNodesVisitedAtDist == nCurrentDistanceSize) {
				nNodesVisitedAtDist = 0;
				nCurrentDistanceSize = queueToVisit.size();
				nCurrentDistance++;

				if (nCurrentDistanceSize == 0) {
					break;
				}
			}

			int ix = queueToVisit.front();
			queueToVisit.pop();
			nNodesVisitedAtDist++;

			setVisited.insert(ix);

			if (dSign * data[ix] > dSign * dValue) {
				fExtrema = false;
				break;
			}

			if (nCurrentDistance < nMaxGraphDistance) {
				for (size_t n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
					int ixNeighbor = grid.m_vecConnectivity[ix][n];
					if (setVisited.find(ixNeighbor) == setVisited.end()) {
						queueToVisit.push(grid.m_vecConnectivity[ix][n]);
					}
				}
			}
		}
		if (fExtrema) {
			setMinMax.insert(*it);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void FindLocalAverage(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	int ix0,
	double dMaxDist,
	real & dAverage
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
	real dSum = 0.0;
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

template <typename real>
void FindMaxClosedContourDelta(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	int ix0,
	double dDistDeg,
	double dMinMaxDistDeg,
	bool fMaxClosedContourDeltaSign,
	real & dMaxClosedContourDelta
) {
	_ASSERT((ix0 >= 0) && (ix0 < grid.GetSize()));

	// Verify that dDistDeg is less than 180.0
	if ((dDistDeg <= 0.0) || (dDistDeg > 180.0)) {
		_EXCEPTIONT("DistDeg must be positive, between 0.0 and 180.0");
	}
	if ((dMinMaxDistDeg < 0.0) || (dMinMaxDistDeg > 180.0)) {
		_EXCEPTIONT("MinMaxDistDeg must be nonnegative, between 0.0 and 180.0");
	}

	if ((grid.m_vecConnectivity.size() != grid.m_dLon.GetRows()) ||
		(grid.m_vecConnectivity.size() != grid.m_dLat.GetRows())
	) {
		_EXCEPTION3("Inconsistent SimpleGrid (%i %i %i); connectivity information mismatch",
			grid.m_vecConnectivity.size(),
			grid.m_dLon.GetRows(),
			grid.m_dLat.GetRows());
	}

	// Convert fMaxClosedContourDeltaSign to double
	double dMCCDSign = 1.0;
	if (!fMaxClosedContourDeltaSign) {
		dMCCDSign = -1.0;
	}

	// Find new minimum (if using positive closed contour delta)
	// or maximum (if using negative closed contour delta)
	if (dMinMaxDistDeg != 0.0) {
		int ixExtremum;
		real dMaxValue;
		float dRMax;

		FindLocalMinMax<real>(
			grid,
			fMaxClosedContourDeltaSign,
			dataState,
			ix0,
			dMinMaxDistDeg,
			ixExtremum,
			dMaxValue,
			dRMax);

		ix0 = ixExtremum;
	}

	// Value of the field at the index point
	double dValue0 = dataState[ix0];

	// Central lat/lon and Cartesian coord
	double dLon0 = grid.m_dLon[ix0];
	double dLat0 = grid.m_dLat[ix0];

	// Initial max closed contour delta
	dMaxClosedContourDelta = 0.0;

	// Priority queue mapping deltas to indices
	std::map<double, int> mapPriorityQueue;
	mapPriorityQueue.insert(
		std::pair<double, int>(0.0, ix0));

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
		queueNodes.push(grid.m_vecConnectivity[ix0][n]);
	}

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Loop through all elements
	while (mapPriorityQueue.size() != 0) {
		std::map<double, int>::iterator iter = mapPriorityQueue.begin();

		const double dDelta = iter->first;
		const int ix = iter->second;

		if ((ix < 0) || (ix >= grid.m_vecConnectivity.size())) {
			_EXCEPTION2("Out of range index in connectivity matrix (%i/%i)",
				ix, grid.m_vecConnectivity.size());
		}

		mapPriorityQueue.erase(iter);

		setNodesVisited.insert(ix);

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			int ixNeighbor = grid.m_vecConnectivity[ix][n];
			if (setNodesVisited.find(ixNeighbor) != setNodesVisited.end()) {
				continue;
			}

			double dValue = dataState[ixNeighbor];
			mapPriorityQueue.insert(
				std::pair<double, int>(
					dMCCDSign * (dValue - dValue0), ixNeighbor));
		}

		// Don't perform calculation on central node
		if (ix == ix0) {
			continue;
		}

		// lat/lon and Cartesian coords of this point
		double dLat = grid.m_dLat[ix];
		double dLon = grid.m_dLon[ix];

		// Great circle distance to this element (in degrees)
		double dR = GreatCircleDistance_Deg(dLon0, dLat0, dLon, dLat);

		if (dR >= dDistDeg) {
			break;
		}

		// Store new maximum delta
		if (dDelta > dMaxClosedContourDelta) {
			dMaxClosedContourDelta = dDelta;
		}
	}

	dMaxClosedContourDelta *= dMCCDSign;
}

///////////////////////////////////////////////////////////////////////////////
// Explicit template instantiation

template void FindLocalMinMax<float>(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataArray1D<float> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	float & dMaxValue,
	float & dRMax
);

template void FindLocalMinMax<double>(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataArray1D<double> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	double & dMaxValue,
	float & dRMax
);

template void FindAllLocalMinima<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & data,
	std::set<int> & setMinima
);

template void FindAllLocalMinima<double>(
	const SimpleGrid & grid,
	const DataArray1D<double> & data,
	std::set<int> & setMinima
);

template void FindAllLocalMaxima<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & data,
	std::set<int> & setMaxima
);

template void FindAllLocalMaxima<double>(
	const SimpleGrid & grid,
	const DataArray1D<double> & data,
	std::set<int> & setMaxima
);

template void FindAllLocalMinMaxWithThreshold<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & data,
	bool fMinima,
	const std::string & strThreshold,
	std::set<int> & setMinima
);

template void FindAllLocalMinMaxWithThreshold<double>(
	const SimpleGrid & grid,
	const DataArray1D<double> & data,
	bool fMinima,
	const std::string & strThreshold,
	std::set<int> & setMinima
);

template void FindAllLocalMinMaxWithGraphDistance<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & data,
	bool fMinima,
	int nMaxGraphDistance,
	std::set<int> & setMinMax
);

template void FindAllLocalMinMaxWithGraphDistance<double>(
	const SimpleGrid & grid,
	const DataArray1D<double> & data,
	bool fMinima,
	int nMaxGraphDistance,
	std::set<int> & setMinMax
);

template void FindLocalAverage<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & data,
	int ix0,
	double dMaxDist,
	float & dAverage
);

template void FindLocalAverage<double>(
	const SimpleGrid & grid,
	const DataArray1D<double> & data,
	int ix0,
	double dMaxDist,
	double & dAverage
);

template void FindMaxClosedContourDelta<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & data,
	int ix0,
	double dDist,
	double dMinMaxDist,
	bool fMaxClosedContourDeltaSign,
	float & dMaxClosedContourDelta
);

template void FindMaxClosedContourDelta<double>(
	const SimpleGrid & grid,
	const DataArray1D<double> & data,
	int ix0,
	double dDist,
	double dMinMaxDist,
	bool fMaxClosedContourDeltaSign,
	double & dMaxClosedContourDelta
);

///////////////////////////////////////////////////////////////////////////////
