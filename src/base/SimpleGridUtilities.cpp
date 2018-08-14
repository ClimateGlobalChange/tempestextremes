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

#include <queue>

///////////////////////////////////////////////////////////////////////////////

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

template <typename real>
void FindLocalAverage(
	const SimpleGrid & grid,
	const DataVector<real> & data,
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
// Explicit template instantiation

template void FindLocalMinMax<float>(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataVector<float> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	float & dMaxValue,
	float & dRMax
);

template void FindLocalMinMax<double>(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataVector<double> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	double & dMaxValue,
	float & dRMax
);

template void FindAllLocalMinima<float>(
	const SimpleGrid & grid,
	const DataVector<float> & data,
	std::set<int> & setMinima
);

template void FindAllLocalMinima<double>(
	const SimpleGrid & grid,
	const DataVector<double> & data,
	std::set<int> & setMinima
);

template void FindAllLocalMaxima<float>(
	const SimpleGrid & grid,
	const DataVector<float> & data,
	std::set<int> & setMaxima
);

template void FindAllLocalMaxima<double>(
	const SimpleGrid & grid,
	const DataVector<double> & data,
	std::set<int> & setMaxima
);

template void FindLocalAverage<float>(
	const SimpleGrid & grid,
	const DataVector<float> & data,
	int ix0,
	double dMaxDist,
	float & dAverage
);

template void FindLocalAverage<double>(
	const SimpleGrid & grid,
	const DataVector<double> & data,
	int ix0,
	double dMaxDist,
	double & dAverage
);

///////////////////////////////////////////////////////////////////////////////
