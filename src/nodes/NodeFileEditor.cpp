///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileEditor.cpp
///	\author  Paul Ullrich
///	\version August 30, 2018
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

/*
#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif
*/

#include "Constants.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "Variable.h"
#include "AutoCurator.h"
#include "DataArray2D.h"
#include "CalculationList.h"
#include "STLStringHelper.h"
#include "NodeFileUtilities.h"
#include "GridElements.h"
#include "RLLPolygonArray.h"
#include "TimeMatch.h"
#include "NodeOutputOp.h"

#include "netcdfcpp.h"

#include <fstream>
#include <queue>
#include <set>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a filtering operator.
///	</summary>
class FilterOp {

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
	///		Constructor.
	///	</summary>
	FilterOp() :
		m_iColumn(-1),
		m_eOp(GreaterThan),
		m_dValue(0.0)
	{ }

public:
	///	<summary>
	///		Parse a filter operator string.
	///	</summary>
	void Parse(
		const ColumnDataHeader & cdhInput,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Column,
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Invalid
		} eReadMode = ReadMode_Column;

		// Loop through string
		int iLast = 0;
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in column name
				if (eReadMode == ReadMode_Column) {
					m_iColumn = cdhInput.GetIndexFromString(strSubStr);
					if (m_iColumn == (-1)) {
						_EXCEPTION1("Unknown column header \"%s\"", strSubStr.c_str());
					}

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
						_EXCEPTION1("Filter invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Value;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					m_dValue = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in filter op \"%s\""
							"\nRequired: \"<name>,<operation>,<value>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in filter op \"%s\""
					"\nRequired: \"<name>,<operation>,<value>\"",
					strOp.c_str());
		}

		// Output announcement
		std::string strDescription = cdhInput[m_iColumn];
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

		Announce("%s %f", strDescription.c_str(), m_dValue);
	}

	///	<summary>
	///		Check if the given value satisfies the filter.
	///	</summary>
	bool Satisfies(
		double dValue
	) const {
		if (m_eOp == FilterOp::GreaterThan) {
			if (dValue > m_dValue) {
				return true;
			}

		} else if (m_eOp == FilterOp::LessThan) {
			if (dValue < m_dValue) {
				return true;
			}

		} else if (m_eOp == FilterOp::GreaterThanEqualTo) {
			if (dValue >= m_dValue) {
				return true;
			}

		} else if (m_eOp == FilterOp::LessThanEqualTo) {
			if (dValue <= m_dValue) {
				return true;
			}

		} else if (m_eOp == FilterOp::EqualTo) {
			if (dValue == m_dValue) {
				return true;
			}

		} else if (m_eOp == FilterOp::NotEqualTo) {
			if (dValue != m_dValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		return false;
	}

public:
	///	<summary>
	///		Column index to use for filter.
	///	</summary>
	int m_iColumn;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

void CalculateRadialProfile(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varix,
	std::string strBins,
	std::string strBinWidth
) {
	// Get number of bins
	int nBins = pathnode.GetColumnDataAsInteger(cdh, strBins);

	// Get bin width
	double dBinWidth = pathnode.GetColumnDataAsDouble(cdh, strBinWidth);

	// Check arguments
	if (nBins <= 0) {
		_EXCEPTIONT("\nNonpositive value of <bins> argument given");
	}
	if (dBinWidth <= 0.0) {
		_EXCEPTIONT("\nNonpositive value of <bin_width> argument given");
	}
	if (static_cast<double>(nBins) * dBinWidth > 180.0) {
		_EXCEPTIONT("\n<bins> x <bin_width> must be no larger than 180 (degrees)");
	}

	// Get the center grid index
	const int ix0 = static_cast<int>(pathnode.m_gridix);
	if (ix0 < 0) {
		_EXCEPTION1("Invalid grid index (%i) in node file", ix0);
	}

	// Load the data
	Variable & var = varreg.Get(varix);
	var.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataState = var.GetData();

	// Verify that dRadius is less than 180.0
	double dRadius = dBinWidth * static_cast<double>(nBins);

	if ((dRadius < 0.0) || (dRadius > 180.0)) {
		_EXCEPTIONT("Radius must be in the range [0.0, 180.0]");
	}

	// Check grid index
	if (ix0 >= grid.m_vecConnectivity.size()) {
		_EXCEPTION2("Grid index (%i) out of range (< %i)",
			ix0, static_cast<int>(grid.m_vecConnectivity.size()));
	}

	// Central lat/lon and Cartesian coord
	double dLon0 = grid.m_dLon[ix0];
	double dLat0 = grid.m_dLat[ix0];

	// Allocate bins
	std::vector< std::vector<double> > dValues;
	dValues.resize(nBins);

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
		queueNodes.push(grid.m_vecConnectivity[ix0][n]);
	}

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Don't perform calculation on central node
		if (ix == ix0) {
			continue;
		}

		// lat/lon and Cartesian coords of this point
		double dLat = grid.m_dLat[ix];
		double dLon = grid.m_dLon[ix];

		// Great circle distance to this element (in degrees)
		double dR = GreatCircleDistance_Deg(dLon0, dLat0, dLon, dLat);

		if (dR >= dRadius) {
			continue;
		}

		// Determine bin
		int iBin = static_cast<int>(dR / dBinWidth);
		if (iBin >= nBins) {
			_EXCEPTIONT("Logic error");
		}

		dValues[iBin].push_back(dataState[ix]);

		if (iBin < nBins-1) {
			dValues[iBin+1].push_back(dataState[ix]);
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Construct radial profile
	ColumnDataRadialProfile * pdat =
		new ColumnDataRadialProfile;

	pathnode.PushColumnData(pdat);

	pdat->m_dR.resize(nBins);
	pdat->m_dValues.resize(nBins);

	pdat->m_dR[0] = 0.0;
	pdat->m_dValues[0] = 0.0;
	for (int i = 1; i < nBins; i++) {
		double dAvg = 0.0;
		if (dValues[i].size() != 0) {
			for (int j = 0; j < dValues[i].size(); j++) {
				dAvg += dValues[i][j];
			}
			dAvg /= static_cast<double>(dValues[i].size());
		}
		pdat->m_dR[i] = static_cast<double>(i) * dBinWidth;
		pdat->m_dValues[i] = dAvg;
	}

}

///////////////////////////////////////////////////////////////////////////////

void CalculateRadialWindProfile(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varixU,
	VariableIndex varixV,
	std::string strBins,
	std::string strBinWidth
) {
	// Get number of bins
	int nBins = pathnode.GetColumnDataAsInteger(cdh, strBins);

	// Get bin width
	double dBinWidth = pathnode.GetColumnDataAsDouble(cdh, strBinWidth);

	// Check arguments
	if (nBins <= 0) {
		_EXCEPTIONT("\nNonpositive value of <bins> argument given");
	}
	if (dBinWidth <= 0.0) {
		_EXCEPTIONT("\nNonpositive value of <bin_width> argument given");
	}
	if (static_cast<double>(nBins) * dBinWidth > 180.0) {
		_EXCEPTIONT("\n<bins> x <bin_width> must be no larger than 180 (degrees)");
	}

	// Get the center grid index
	const int ix0 = static_cast<int>(pathnode.m_gridix);
	if (ix0 < 0) {
		_EXCEPTION1("Invalid grid index (%i) in node file", ix0);
	}

	// Load the zonal wind data
	Variable & varU = varreg.Get(varixU);
	varU.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataStateU = varU.GetData();

	// Load the meridional wind data
	Variable & varV = varreg.Get(varixV);
	varV.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataStateV = varV.GetData();

	// Verify that dRadius is less than 180.0
	double dRadius = dBinWidth * static_cast<double>(nBins);

	if ((dRadius < 0.0) || (dRadius > 180.0)) {
		_EXCEPTIONT("Radius must be in the range [0.0, 180.0]");
	}

	// Check grid index
	if (ix0 >= grid.m_vecConnectivity.size()) {
		_EXCEPTION2("Grid index (%i) out of range (< %i)",
			ix0, static_cast<int>(grid.m_vecConnectivity.size()));
	}

	// Central lat/lon and Cartesian coord
	double dLon0 = grid.m_dLon[ix0];
	double dLat0 = grid.m_dLat[ix0];

	double dX0 = cos(dLon0) * cos(dLat0);
	double dY0 = sin(dLon0) * cos(dLat0);
	double dZ0 = sin(dLat0);

	// Allocate bins
	std::vector< std::vector<double> > dVelocities;
	dVelocities.resize(nBins);

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
		queueNodes.push(grid.m_vecConnectivity[ix0][n]);
	}

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Don't perform calculation on central node
		if (ix == ix0) {
			continue;
		}

		// lat/lon and Cartesian coords of this point
		double dLat = grid.m_dLat[ix];
		double dLon = grid.m_dLon[ix];

		double dX = cos(dLon) * cos(dLat);
		double dY = sin(dLon) * cos(dLat);
		double dZ = sin(dLat);

		// Great circle distance to this element (in degrees)
		double dR = GreatCircleDistance_Deg(dLon0, dLat0, dLon, dLat);

		if (dR >= dRadius) {
			continue;
		}

		// Velocities at this location
		double dUlon = dataStateU[ix];
		double dUlat = dataStateV[ix];

		// Cartesian velocities at this location
		double dUx = - sin(dLat) * cos(dLon) * dUlat - sin(dLon) * dUlon;
		double dUy = - sin(dLat) * sin(dLon) * dUlat + cos(dLon) * dUlon;
		double dUz = cos(dLat) * dUlat;

		double dUtot = sqrt(dUlon * dUlon + dUlat * dUlat);

		// Calculate local radial vector from central lat/lon
		// i.e. project \vec{X} - \vec{X}_0 to the surface of the
		//      sphere and normalize to unit length.
		double dRx = dX - dX0;
		double dRy = dY - dY0;
		double dRz = dZ - dZ0;

		double dDot = dRx * dX + dRy * dY + dRz * dZ;

		dRx -= dDot * dX;
		dRy -= dDot * dY;
		dRz -= dDot * dZ;

		double dMag = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

		dRx /= dMag;
		dRy /= dMag;
		dRz /= dMag;

		// Calculate local azimuthal velocity vector
		double dAx = dY * dRz - dZ * dRy;
		double dAy = dZ * dRx - dX * dRz;
		double dAz = dX * dRy - dY * dRx;
/*
		dDot = dAx * dAx + dAy * dAy + dAz * dAz;
		if (fabs(dDot - 1.0) > 1.0e-10) {
			std:: cout << dDot << std::endl;
			_EXCEPTIONT("Logic error");
		}
*/
		// Calculate radial velocity
		//double dUr = dUx * dRx + dUy * dRy + dUz * dRz;

		// Calculate azimuthal velocity
		double dUa = dUx * dAx + dUy * dAy + dUz * dAz;
		
		// Azimuthal convention positive if cyclonic, flip in SH
		if (dLat0 < 0.0) {
			dUa = -dUa;
		}

		//printf("%1.5e %1.5e :: %1.5e %1.5e\n", dUlon, dUlat, dUr, dUa);

		// Determine bin
		int iBin = static_cast<int>(dR / dBinWidth);
		if (iBin >= nBins) {
			_EXCEPTIONT("Logic error");
		}

		dVelocities[iBin].push_back(dUa);

		if (iBin < nBins-1) {
			dVelocities[iBin+1].push_back(dUa);
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Construct radial profile
	ColumnDataRadialVelocityProfile * pdat =
		new ColumnDataRadialVelocityProfile;

	pathnode.PushColumnData(pdat);

	pdat->m_dR.resize(nBins);
	pdat->m_dUa.resize(nBins);
	pdat->m_dUr.resize(nBins);

	pdat->m_dR[0] = 0.0;
	pdat->m_dUa[0] = 0.0;
	pdat->m_dUr[0] = 0.0;
	for (int i = 1; i < nBins; i++) {
		double dAvg = 0.0;
		if (dVelocities[i].size() != 0) {
			for (int j = 0; j < dVelocities[i].size(); j++) {
				dAvg += dVelocities[i][j];
			}
			dAvg /= static_cast<double>(dVelocities[i].size());
		}
		pdat->m_dR[i] = static_cast<double>(i) * dBinWidth;
		pdat->m_dUa[i] = dAvg;
		pdat->m_dUr[i] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void CalculateQuadrantProfile(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varix,
	std::string strBins,
	std::string strBinWidth
) {
	// Get number of bins
	int nBins = pathnode.GetColumnDataAsInteger(cdh, strBins);

	// Get bin width
	double dBinWidth = pathnode.GetColumnDataAsDouble(cdh, strBinWidth);

	// Check arguments
	if (nBins <= 0) {
		_EXCEPTIONT("\nNonpositive value of <bins> argument given");
	}
	if (dBinWidth <= 0.0) {
		_EXCEPTIONT("\nNonpositive value of <bin_width> argument given");
	}
	if (static_cast<double>(nBins) * dBinWidth > 180.0) {
		_EXCEPTIONT("\n<bins> x <bin_width> must be no larger than 180 (degrees)");
	}

	// Get the center grid index
	const int ix0 = static_cast<int>(pathnode.m_gridix);
	if (ix0 < 0) {
		_EXCEPTION1("Invalid grid index (%i) in node file", ix0);
	}

	// Load the zonal wind data
	Variable & varData = varreg.Get(varix);
	varData.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataState = varData.GetData();

	// Verify that dRadius is less than 180.0
	double dRadius = dBinWidth * static_cast<double>(nBins);

	if ((dRadius < 0.0) || (dRadius > 180.0)) {
		_EXCEPTIONT("Radius must be in the range [0.0, 180.0]");
	}

	// Check grid index
	if (ix0 >= grid.m_vecConnectivity.size()) {
		_EXCEPTION2("Grid index (%i) out of range (< %i)",
			ix0, static_cast<int>(grid.m_vecConnectivity.size()));
	}

	// Central lat/lon and Cartesian coord
	double dLonRad0 = grid.m_dLon[ix0];
	double dLatRad0 = grid.m_dLat[ix0];

	double dX0 = cos(dLonRad0) * cos(dLatRad0);
	double dY0 = sin(dLonRad0) * cos(dLatRad0);
	double dZ0 = sin(dLatRad0);

	// Allocate bins
	std::vector< std::vector<double> > dBinTotal;
	dBinTotal.resize(4);
	dBinTotal[0].resize(nBins, 0.0);
	dBinTotal[1].resize(nBins, 0.0);
	dBinTotal[2].resize(nBins, 0.0);
	dBinTotal[3].resize(nBins, 0.0);

	std::vector< std::vector<int> > dBinCount;
	dBinCount.resize(4);
	dBinCount[0].resize(nBins, 0.0);
	dBinCount[1].resize(nBins, 0.0);
	dBinCount[2].resize(nBins, 0.0);
	dBinCount[3].resize(nBins, 0.0);

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
		queueNodes.push(grid.m_vecConnectivity[ix0][n]);
	}

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Don't perform calculation on central node
		if (ix == ix0) {
			continue;
		}

		// lat/lon and Cartesian coords of this point
		double dLatRad = grid.m_dLat[ix];
		double dLonRad = grid.m_dLon[ix];

		double dX = cos(dLonRad) * cos(dLatRad);
		double dY = sin(dLonRad) * cos(dLatRad);
		double dZ = sin(dLatRad);

		// Great circle distance to this element (in degrees)
		double dR = GreatCircleDistance_Deg(dLonRad0, dLatRad0, dLonRad, dLatRad);

		if (dR >= dRadius) {
			continue;
		}

		double dDataValue = dataState[ix];

		// Determine quadrant
		double dXs;
		double dYs;
		StereographicProjection(dLonRad0, dLatRad0, dLonRad, dLatRad, dXs, dYs);
		
		int iQuadrant = 0;
		if (dXs > 0.0) {
			if (dYs > 0.0) {
				iQuadrant = 0;
			} else {
				iQuadrant = 1;
			}
		} else {
			if (dYs > 0.0) {
				iQuadrant = 3;
			} else {
				iQuadrant = 2;
			}
		}

		// Determine bin and update total velocity in bin
		int iBin = static_cast<int>(dR / dBinWidth);
		if (iBin >= nBins) {
			_EXCEPTIONT("Logic error");
		}

		dBinTotal[iQuadrant][iBin] += dataState[ix];
		dBinCount[iQuadrant][iBin]++;

		if (iBin < nBins-1) {
			dBinTotal[iQuadrant][iBin+1] += dataState[ix];
			dBinCount[iQuadrant][iBin+1]++;
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Construct radial profile
	ColumnData2DArray * pdat =
		new ColumnData2DArray;

	pathnode.PushColumnData(pdat);

	pdat->m_dValues.resize(5);
	pdat->m_dValues[0].resize(nBins);
	pdat->m_dValues[1].resize(nBins);
	pdat->m_dValues[2].resize(nBins);
	pdat->m_dValues[3].resize(nBins);
	pdat->m_dValues[4].resize(nBins);

	for (int j = 1; j < nBins; j++) {
		pdat->m_dValues[0][j] = static_cast<double>(j) * dBinWidth;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 1; j < nBins; j++) {
			if (dBinCount[i][j] != 0) {
				pdat->m_dValues[i+1][j] = dBinTotal[i][j] / static_cast<double>(dBinCount[i][j]);
			} else {
				pdat->m_dValues[i+1][j] = 0.0;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void CalculateQuadrantWindProfile(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varixU,
	VariableIndex varixV,
	std::string strBins,
	std::string strBinWidth
) {
	// Get number of bins
	int nBins = pathnode.GetColumnDataAsInteger(cdh, strBins);

	// Get bin width
	double dBinWidth = pathnode.GetColumnDataAsDouble(cdh, strBinWidth);

	// Check arguments
	if (nBins <= 0) {
		_EXCEPTIONT("\nNonpositive value of <bins> argument given");
	}
	if (dBinWidth <= 0.0) {
		_EXCEPTIONT("\nNonpositive value of <bin_width> argument given");
	}
	if (static_cast<double>(nBins) * dBinWidth > 180.0) {
		_EXCEPTIONT("\n<bins> x <bin_width> must be no larger than 180 (degrees)");
	}

	// Get the center grid index
	const int ix0 = static_cast<int>(pathnode.m_gridix);
	if (ix0 < 0) {
		_EXCEPTION1("Invalid grid index (%i) in node file", ix0);
	}

	// Load the zonal wind data
	Variable & varU = varreg.Get(varixU);
	varU.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataStateU = varU.GetData();

	// Load the meridional wind data
	Variable & varV = varreg.Get(varixV);
	varV.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataStateV = varV.GetData();

	// Verify that dRadius is less than 180.0
	double dRadius = dBinWidth * static_cast<double>(nBins);

	if ((dRadius < 0.0) || (dRadius > 180.0)) {
		_EXCEPTIONT("Radius must be in the range [0.0, 180.0]");
	}

	// Check grid index
	if (ix0 >= grid.m_vecConnectivity.size()) {
		_EXCEPTION2("Grid index (%i) out of range (< %i)",
			ix0, static_cast<int>(grid.m_vecConnectivity.size()));
	}

	// Central lat/lon and Cartesian coord
	double dLonRad0 = grid.m_dLon[ix0];
	double dLatRad0 = grid.m_dLat[ix0];

	double dX0 = cos(dLonRad0) * cos(dLatRad0);
	double dY0 = sin(dLonRad0) * cos(dLatRad0);
	double dZ0 = sin(dLatRad0);

	// Allocate bins
	std::vector< std::vector<double> > dVelocities;
	dVelocities.resize(4);
	dVelocities[0].resize(nBins, 0.0);
	dVelocities[1].resize(nBins, 0.0);
	dVelocities[2].resize(nBins, 0.0);
	dVelocities[3].resize(nBins, 0.0);

	std::vector< std::vector<int> > dVelocitiesCount;
	dVelocitiesCount.resize(4);
	dVelocitiesCount[0].resize(nBins, 0.0);
	dVelocitiesCount[1].resize(nBins, 0.0);
	dVelocitiesCount[2].resize(nBins, 0.0);
	dVelocitiesCount[3].resize(nBins, 0.0);

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
		queueNodes.push(grid.m_vecConnectivity[ix0][n]);
	}

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Don't perform calculation on central node
		if (ix == ix0) {
			continue;
		}

		// lat/lon and Cartesian coords of this point
		double dLatRad = grid.m_dLat[ix];
		double dLonRad = grid.m_dLon[ix];

		double dX = cos(dLonRad) * cos(dLatRad);
		double dY = sin(dLonRad) * cos(dLatRad);
		double dZ = sin(dLatRad);

		// Great circle distance to this element (in degrees)
		double dR = GreatCircleDistance_Deg(dLonRad0, dLatRad0, dLonRad, dLatRad);

		if (dR >= dRadius) {
			continue;
		}

		// Velocities at this location
		double dUlon = dataStateU[ix];
		double dUlat = dataStateV[ix];

		double dUtot = sqrt(dUlon * dUlon + dUlat * dUlat);

		// Determine quadrant
		double dXs;
		double dYs;
		StereographicProjection(dLonRad0, dLatRad0, dLonRad, dLatRad, dXs, dYs);
		
		int iQuadrant = 0;
		if (dXs > 0.0) {
			if (dYs > 0.0) {
				iQuadrant = 0;
			} else {
				iQuadrant = 1;
			}
		} else {
			if (dYs > 0.0) {
				iQuadrant = 3;
			} else {
				iQuadrant = 2;
			}
		}

		// Determine bin and update total velocity in bin
		int iBin = static_cast<int>(dR / dBinWidth);
		if (iBin >= nBins) {
			_EXCEPTIONT("Logic error");
		}

		dVelocities[iQuadrant][iBin] += dUtot;
		dVelocitiesCount[iQuadrant][iBin]++;

		if (iBin < nBins-1) {
			dVelocities[iQuadrant][iBin+1] += dUtot;
			dVelocitiesCount[iQuadrant][iBin+1]++;
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Construct radial profile
	ColumnData2DArray * pdat =
		new ColumnData2DArray;

	pathnode.PushColumnData(pdat);

	pdat->m_dValues.resize(5);
	pdat->m_dValues[0].resize(nBins);
	pdat->m_dValues[1].resize(nBins);
	pdat->m_dValues[2].resize(nBins);
	pdat->m_dValues[3].resize(nBins);
	pdat->m_dValues[4].resize(nBins);

	for (int j = 1; j < nBins; j++) {
		pdat->m_dValues[0][j] = static_cast<double>(j) * dBinWidth;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 1; j < nBins; j++) {
			if (dVelocitiesCount[i][j] != 0) {
				pdat->m_dValues[i+1][j] = dVelocities[i][j] / static_cast<double>(dVelocitiesCount[i][j]);
			} else {
				pdat->m_dValues[i+1][j] = 0.0;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void SumRadius(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varix,
	std::string strRadius
) {
	// Get the radius 
	double dRadius = pathnode.GetColumnDataAsDouble(cdh, strRadius);

	if (dRadius == 0.0) {
		pathnode.PushColumnData(
			new ColumnDataDouble(0.0));
	} else {
		pathnode.PushColumnData(
			new ColumnDataDouble(1.0));
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Polynomial cubic Hermite interpolation function using Fritsch-Carlson
///		method.
///	</summary>
///	<reference>
///		https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
///	</reference>
///	<param name="n">
///		Number of points to use in interpolant.
///	</param>
///	<param name="x">
///		Array of size n with independent coordinate values of interpolant.
///		The values of x must be monotone increasing or behavior is undefined.
///	</param>
///	<param name="y">
///		Array of size n with dependent coordinate values of interpolant.
///	</param>
///	<param name="c">
///		Output array of cubic polynomial coefficients of size (n-1) x 4.
///		The interpolation function is then defined by
///		  xdiff = x - x[i]
///		  f_i(x) = c[i][0]
///		         + xdiff * (c[i][1] + xdiff * (c[i][2] + xdiff * c[i][3]))
///	</param>
///	<param name="work">
///		Work array of size 3*n.
///	</param>
///	<returns>
///		(-1) if the value of n is less than 2.
///		(0) if the calculation completed successfully.
///	</returns>
int pchip_fc_coeff(
	int n,
	double * const x,
	double * const y,
	double * c,
	double * work
) {
	if (n < 2) {
		return (-1);
	}

	double * dx = work;
	double * dy = work+n;
	double * m = work+2*n;

	for (int i = 0; i < n-1; i++) {
		dx[i] = x[i+1] - x[i];
		dy[i] = y[i+1] - y[i];
		m[i] = dy[i] / dx[i];
		c[i*4] = y[i];
	}
	m[n-1] = dy[n-2] / dx[n-2];

	c[1] = m[0];
	for (int i = 0; i < n-2; i++) {
		double mi = m[i];
		double mn = m[i+1];
		if (mi * mn <= 0.0) {
			c[4*(i+1)+1] = 0.0;
		} else {
			double dc = dx[i] + dx[i+1];
			c[4*(i+1)+1] = 3.0 * dc / ((dc + dx[i+1]) / mi + (dc + dx[i]) / mn);
		}
	}

	for (int i = 0; i < n-2; i++) {
		double c1 = c[4*i+1];
		double invdx = 1.0/dx[i];
		double dc = c1 + c[4*(i+1)+1] - 2.0 * m[i];

		c[4*i+2] = (m[i] - c1 - dc) * invdx;
		c[4*i+3] = dc * invdx * invdx;
	}
	{
		int i = n-2;
		double c1 = c[4*i+1];
		double invdx = 1.0/dx[i];
		double dc = c1 + m[n-1] - 2.0 * m[i];

		c[4*i+2] = (m[i] - c1 - dc) * invdx;
		c[4*i+3] = dc * invdx * invdx;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

void CalculateStormVelocity(
	const SimpleGrid & grid,
	Path & path
) {
	const double MetersPerRadian = 111.325 * 1000.0 * 180.0 / M_PI;

	int nPathNodes = path.size();
	if (nPathNodes < 2) {
		_EXCEPTIONT("Path must contain at least two nodes");
	}

	DataArray1D<double> dX(nPathNodes);
	DataArray1D<double> dY(nPathNodes);
	DataArray2D<double> dC(nPathNodes-1, 4);
	DataArray1D<double> dWork(nPathNodes*3);

	// Velocity column data
	std::vector<ColumnDataRLLVelocity *> vecPathNodeVelocity;
	vecPathNodeVelocity.resize(nPathNodes);

	// Calculate meridional velocity
	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path[i];

		dX[i] = pathnode.m_time - path.m_timeStart;

		if ((pathnode.m_gridix < 0) || (pathnode.m_gridix >= grid.m_dLat.GetRows())) {
			_EXCEPTION2("Grid index of path (%li) out of range "
				"(only %i nodes in grid)",
				pathnode.m_gridix,
				grid.m_dLat.GetRows());
		}

		dY[i] = grid.m_dLat[pathnode.m_gridix];

		ColumnDataRLLVelocity * pcd = new ColumnDataRLLVelocity();
		vecPathNodeVelocity.push_back(pcd);
		pathnode.m_vecColumnData.push_back(pcd);
	}

	if (vecPathNodeVelocity.size() != nPathNodes) {
		_EXCEPTIONT("Logic error");
	}

	const double dFinalDeltaT = dX[nPathNodes-1] - dX[nPathNodes-2];

	pchip_fc_coeff(nPathNodes, &(dX[0]), &(dY[0]), &(dC[0][0]), &(dWork[0]));

	// Array of velocity column data
	for (int i = 0; i < nPathNodes-1; i++) {
		PathNode & pathnode = path[i];
		vecPathNodeVelocity[i]->m_dV = dC[i][1];
	}
	vecPathNodeVelocity[nPathNodes-1]->m_dV =
		dC[nPathNodes-2][1]
		+ 2.0 * dC[nPathNodes-2][2] * dFinalDeltaT
		+ 3.0 * dC[nPathNodes-2][3] * dFinalDeltaT * dFinalDeltaT;

	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path[i];
		vecPathNodeVelocity[i]->m_dV *= MetersPerRadian;
	}

	// Calculate zonal velocity
	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path[i];
		if ((pathnode.m_gridix < 0) || (pathnode.m_gridix >= grid.m_dLon.GetRows())) {
			_EXCEPTION1("Grid index (%li) out of range", pathnode.m_gridix);
		}
		dY[i] = grid.m_dLon[pathnode.m_gridix];
	}

	pchip_fc_coeff(nPathNodes, &(dX[0]), &(dY[0]), &(dC[0][0]), &(dWork[0]));

	for (int i = 0; i < nPathNodes-1; i++) {
		vecPathNodeVelocity[i]->m_dU = dC[i][1];
	}
	vecPathNodeVelocity[nPathNodes-1]->m_dU =
		dC[nPathNodes-2][1]
		+ 2.0 * dC[nPathNodes-2][2] * dFinalDeltaT
		+ 3.0 * dC[nPathNodes-2][3] * dFinalDeltaT * dFinalDeltaT;

	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path[i];
		vecPathNodeVelocity[i]->m_dU *=
			MetersPerRadian * cos(grid.m_dLat[pathnode.m_gridix]);
	}

/*
	for (int i = 0; i < nPathNodes-1; i++) {
		printf("%1.5e\n", path[i].m_dVelocityLat);
		printf("%1.5e %1.5e :: %1.5e %1.5e %1.5e %1.5e\n", dX[i], dY[i], dC[i][0], dC[i][1], dC[i][2], dC[i][3]); //path[i].m_dVelocityLat);
	}
	printf("%1.5e\n", path[nPathNodes-1].m_dVelocityLat);
	printf("%1.5e %1.5e\n", dX[nPathNodes-1], dY[nPathNodes-1]);
*/
/*
	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path[i];
		printf("%1.5e %1.5e\n", pathnode.m_dVelocityLat, pathnode.m_dVelocityLon);
	}
*/
	//_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

void MaxClosedContourDelta(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varix,
	std::string strRadius,
	std::string strIndex,
	bool fAddRootValue
) {
	// Get radius
	double dRadius = pathnode.GetColumnDataAsDouble(cdh, strRadius);

	// Check arguments
	if (dRadius <= 0.0) {
		_EXCEPTIONT("\nNonpositive value of <radius> argument given");
	}
	if (dRadius > 180.0) {
		_EXCEPTIONT("\n<radius> must be no larger than 180 (degrees)");
	}

	// Get index of the centerpoint
	int ix0;
	if (strIndex == "") {
		ix0 = pathnode.m_gridix;
	} else {
		ix0 = pathnode.GetColumnDataAsInteger(cdh, strIndex);
	}

	if ((ix0 < 0) || (ix0 >= grid.GetSize())) {
		_EXCEPTION2("Initial index point out of range (%i/%i)",
			ix0, grid.GetSize());
	}

	// Load the variable data
	Variable & var = varreg.Get(varix);
	var.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataState = var.GetData();

	if (dataState.GetRows() != grid.GetSize()) {
		_EXCEPTION2("Inconsistent data array (size %i) / grid (size %i)",
			dataState.GetRows(),
			grid.GetSize());
	}
	if ((grid.m_vecConnectivity.size() != grid.m_dLon.GetRows()) ||
		(grid.m_vecConnectivity.size() != grid.m_dLat.GetRows())
	) {
		_EXCEPTION3("Inconsistent SimpleGrid (%i %i %i)",
			grid.m_vecConnectivity.size(),
			grid.m_dLon.GetRows(),
			grid.m_dLat.GetRows());
	}

	// Maximum closed contour delta
	double dMaxDelta = 0.0;

	// Value of the field at the index point
	double dValue0 = dataState[ix0];

	// Central lat/lon and Cartesian coord
	double dLon0 = grid.m_dLon[ix0];
	double dLat0 = grid.m_dLat[ix0];

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

	// Loop through all latlon elements
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
					dValue - dValue0, ixNeighbor));
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

		if (dR >= dRadius) {
			break;
		}

		// Store new maximum delta
		if (dDelta > dMaxDelta) {
			dMaxDelta = dDelta;
		}
	}

	// Add the root value in
	if (fAddRootValue) {
		dMaxDelta += dValue0;
	}

	// Store the maximum closed contour delta as new column data
	ColumnDataDouble * pdat =
		new ColumnDataDouble(dMaxDelta);
	_ASSERT(pdat != NULL);

	pathnode.PushColumnData(pdat);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Types of standard cyclone metrics.
///	</summary>
enum CycloneMetric {
	CycloneMetric_ACEPSL,
	CycloneMetric_ACE,
	CycloneMetric_IKE,
	CycloneMetric_PDI,
	CycloneMetric_MIN,
	CycloneMetric_MAX
};

///	<summary>
///		Calculate popular cyclone metrics, including:
///		Accumulated Cyclone Energy (ACE),
///		Integrated Kinetic Energy (IKE),
///		Potential Dissipation Index (PDI)
///	</summary>
void CalculateCycloneMetrics(
	CycloneMetric eCycloneMetric,
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	VariableIndex varixU,
	VariableIndex varixV,
	std::string strRadius
) {
	// Get the radius of the calculation (in great circle degrees)
	double dRadius = pathnode.GetColumnDataAsDouble(cdh, strRadius);

	// Check arguments
	if ((dRadius < 0.0) || (dRadius > 180.0)) {
		_EXCEPTIONT("Radius must be in the range [0.0, 180.0]");
	}

	// Get the center grid index
	const int ix0 = pathnode.m_gridix;

	// Load the zonal wind data
	Variable & varU = varreg.Get(varixU);
	varU.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataStateU = varU.GetData();

	// Load the meridional wind data
	Variable & varV = varreg.Get(varixV);
	varV.LoadGridData(varreg, vecFiles, grid);
	const DataArray1D<float> & dataStateV = varV.GetData();

	// Check grid index
	if ((ix0 < 0) || (ix0 >= grid.m_vecConnectivity.size())) {
		_EXCEPTION2("Grid index (%i) out of range (< %i)",
			ix0, static_cast<int>(grid.m_vecConnectivity.size()));
	}

	// Central lat/lon and Cartesian coord
	double dLon0 = grid.m_dLon[ix0];
	double dLat0 = grid.m_dLat[ix0];

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
		queueNodes.push(grid.m_vecConnectivity[ix0][n]);
	}

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Value
	double dValue = 0.0;
	if ((eCycloneMetric == CycloneMetric_MIN) || (eCycloneMetric == CycloneMetric_MAX)) {
        dValue = dataStateU[ix0];
        }

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Don't perform calculation on central node
		if (ix == ix0) {
			continue;
		}

		// lat/lon and Cartesian coords of this point
		double dLat = grid.m_dLat[ix];
		double dLon = grid.m_dLon[ix];

		// Great circle distance to this element (in degrees)
		double dR = GreatCircleDistance_Deg(dLon0, dLat0, dLon, dLat);

		if (dR >= dRadius) {
			continue;
		}

		// Accumulated Cyclone Energy from PSL (ACEPSL)
		// Formula:  Holland (2008) Revised Hurricane Pressure-Wind Model
		if (eCycloneMetric == CycloneMetric_ACEPSL) {
			double dPSL_hPa = dataStateU[ix] / 100.0;

			if (dPSL_hPa > 1016.0) {
				continue;
			}

			// Since Holland gives their formula in m/s we need to convert to kts
			double dTempValue = KnotsPerMetersPerSecond * 3.92 * pow(1016.0 - dPSL_hPa, 0.644);

			dTempValue *= 1.0e-4 * dTempValue;

			if (dTempValue > dValue) {
				dValue = dTempValue;
			}
		
		} else if (eCycloneMetric == CycloneMetric_MIN) {
				double dTempValue = dataStateU[ix];
				if (dTempValue < dValue) {
					dValue = dTempValue;
				}	
		
		} else if (eCycloneMetric == CycloneMetric_MAX) {
				double dTempValue = dataStateU[ix];
				if (dTempValue > dValue) {
					dValue = dTempValue;
				}

		} else {
			// Velocities at this location
			double dUlon = dataStateU[ix];
			double dUlat = dataStateV[ix];

			// Velocity magnitude
			double dUmag = sqrt(dUlon * dUlon + dUlat * dUlat);

			// Accumulated Cyclone Energy (ACE)
			if (eCycloneMetric == CycloneMetric_ACE) {
				double dUmag_kt = KnotsPerMetersPerSecond * dUmag;
				double dTempValue = dUmag_kt * dUmag_kt * 1.0e-4;
				if (dTempValue > dValue) {
					dValue = dTempValue;
				}

			// Integrated Kinetic Energy (IKE)
			} else if (eCycloneMetric == CycloneMetric_IKE) {
				dValue +=
					0.5
					* dUmag * dUmag
					* grid.m_dArea[ix]
					* EarthRadius * EarthRadius;

			// Potential Dissipation Index (PDI)
			} else if (eCycloneMetric == CycloneMetric_PDI) {
				double dTempValue = dUmag * dUmag * dUmag;
				if (dTempValue > dValue) {
					dValue = dTempValue;
				}

			} else {
				_EXCEPTIONT("Invalid eCycloneMetric value");
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Store the maximum closed contour delta as new column data
	ColumnDataDouble * pdat =
		new ColumnDataDouble(dValue);

	pathnode.PushColumnData(pdat);
}

///////////////////////////////////////////////////////////////////////////////

double FirstWhere(
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	const std::string & strOp,
	const std::string & strThreshold,
	std::string strRadius,
	const std::vector<double> & dIndices,
	const std::vector<double> & dArray
) {
	if (dArray.size() == 0) {
		_EXCEPTIONT("PathNode RadialProfile has zero size");
	}
	if (dIndices.size() != dArray.size()) {
		_EXCEPTIONT("PathNode R array size different from Ua size");
	}

	// Get the radius
	double dRadius = pathnode.GetColumnDataAsDouble(cdh, strRadius);

	if (dRadius > dIndices[dIndices.size() - 1]) {
		_EXCEPTIONT("Radius should be smaller than the largest radius of the radial profile");
	}

    // Get the bin width
	double dBinWidth = dIndices[dIndices.size() - 1] / (dIndices.size() - 1);

	// Get the threshold
	double dThreshold;
	if (strThreshold == "max") {
		dThreshold = dArray[0];
		for (int k = 0; k < dRadius / dBinWidth; k++) {
			if (dArray[k] > dThreshold) {
				dThreshold = dArray[k];
			}
		}

	} else if (strThreshold == "min") {
		dThreshold = dArray[0];
		for (int k = 0; k < dRadius / dBinWidth; k++) {
			if (dArray[k] < dThreshold) {
				dThreshold = dArray[k];
			}
		}

	} else {
		dThreshold =
			pathnode.GetColumnDataAsDouble(cdh, strThreshold);
	}

	// Initial Value
	double dValue = dArray[0];
	int dIndex = 0;

	// Find array index
	int j = 0;

	if (strOp == "fallsbelow") {
		for (; j <= dRadius / dBinWidth; j++) {
			double dTempValue = dArray[j];
			if (dTempValue > dValue) {
				dValue = dTempValue;
				dIndex = j;
			}
		}
		if (dValue < dThreshold) {
			j = 0;
		} else {
			// Get the index of the maximum value
			j = dIndex + 1;
			for (; j < dArray.size(); j++) {
				if (dArray[j] < dThreshold) {
					break;
				}
			}
		}
		
	} else if (strOp == "risesabove") {
		for (; j <= dRadius / dBinWidth ; j++) {
			double dTempValue = dArray[j];
			if (dTempValue < dValue) {
				dValue = dTempValue;
				dIndex = j;
			}
		}
		if (dValue > dThreshold) {
			j = 0;
		} else {
			// Get the index of the minimum value
			j = dIndex + 1;
			for (; j < dArray.size(); j++) {
				if (dArray[j] > dThreshold) {
					break;
				}
			}
		}
        
    } else if (strOp == "=") {
		for (; j <= dRadius / dBinWidth ; j++) {
            if (dArray[j] == dThreshold) {
				break;
			}
		}
        
	} else {
		_EXCEPTION1("Invalid operator \"%s\" in function firstwhere()",
			strOp.c_str());
	}

	if (j < dArray.size()) {
		return dIndices[j];
	} else {
		return dIndices[dArray.size()-1];
	}
}

///////////////////////////////////////////////////////////////////////////////

double LastWhere(
	const ColumnDataHeader & cdh,
	PathNode & pathnode,
	const std::string & strOp,
	const std::string & strThreshold,
	const std::vector<double> & dIndices,
	const std::vector<double> & dArray
) {
	if (dArray.size() == 0) {
		_EXCEPTIONT("PathNode RadialProfile has zero size");
	}
	if (dIndices.size() != dArray.size()) {
		_EXCEPTIONT("PathNode R array size different from Ua size");
	}

	// Get the threshold
	double dThreshold;
	if (strThreshold == "max") {
		dThreshold = dArray[0];
		for (int k = 0; k < dArray.size(); k++) {
			if (dArray[k] > dThreshold) {
				dThreshold = dArray[k];
			}
		}

	} else if (strThreshold == "min") {
		dThreshold = dArray[0];
		for (int k = 0; k < dArray.size(); k++) {
			if (dArray[k] < dThreshold) {
				dThreshold = dArray[k];
			}
		}

	} else {
		dThreshold =
			pathnode.GetColumnDataAsDouble(cdh, strThreshold);
	}

	// Find array index
	int j = dArray.size()-1;
	if (strOp == ">=") {
		for (; j > 0; j--) {
			if (dArray[j] >= dThreshold) {
				break;
			}
		}

	} else if (strOp == ">") {
		for (; j > 0; j--) {
			if (dArray[j] > dThreshold) {
				break;
			}
		}

	} else if (strOp == "<=") {
		for (; j > 0; j--) {
			if (dArray[j] <= dThreshold) {
				break;
			}
		}

	} else if (strOp == "<") {
		for (; j > 0; j--) {
			if (dArray[j] < dThreshold) {
				break;
			}
		}

	} else if (strOp == "=") {
		for (; j > 0; j--) {
			if (dArray[j] == dThreshold) {
				break;
			}
		}

	} else {
		_EXCEPTION1("Invalid operator \"%s\" in function lastwhere()",
			strOp.c_str());
	}

	return dIndices[j];
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
/*
#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Not yet implemented
	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
	if (nMPISize > 1) {
		std::cout << "Sorry!  Parallel processing with MPI is not yet implemented" << std::endl;
		return (-1);
	}
#endif
*/
	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);
/*
	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();
*/
try {

	// Input node file
	std::string strInputNodeFile;

	// Input list of node files
	std::string strInputNodeFileList;

	// Input node file format
	//std::string strInputFileFormat;

	// Input file type
	std::string strPathType;

	// Input data file
	std::string strInputData;

	// Input list of data files
	std::string strInputDataList;

	// Input data index
	std::string strInputDataIndex;

	// Connectivity file
	std::string strConnectivity;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Data is regional
	bool fRegional;

	// Input format (columns of in_file)
	std::string strInputFormat;

	// Output format (columns of out_file)
	std::string strOutputFormat;

	// Output file
	std::string strOutputFile;

	// Output file format
	std::string strOutputFileFormat;

	// Time subset
	std::string strTimeFilter;

	// Filter commands
	std::string strColumnFilter;

	// Calculation commands
	std::string strCalculate;

	// Time delta
	std::string strTimeDelta;

	// List of variables to append
	std::string strAppend;

	// Append output to input file
	bool fOutputAppend;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_nodefile", "");
		//CommandLineString(strInputNodeFileList, "in_file_list", "");
		//CommandLineStringD(strInputFileFormat, "in_nodefile_format", "gfdl", "[gfdl|track]");
		CommandLineStringD(strPathType, "in_nodefile_type", "SN", "[DN|SN]");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strInputDataIndex, "in_data_index", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strInputFormat, "in_fmt", "(auto)");
		CommandLineString(strOutputFormat, "out_fmt", "");

		CommandLineString(strOutputFile, "out_nodefile", "");
		CommandLineStringD(strOutputFileFormat, "out_nodefile_format", "gfdl", "[gfdl|csv|csvnohead]");
		//CommandLineBool(fOutputAppend, "out_append");

		CommandLineString(strTimeFilter, "timefilter", "");
		CommandLineStringD(strColumnFilter, "colfilter", "", "[col,op,value;...]");
		CommandLineString(strCalculate, "calculate", "");
		CommandLineString(strTimeDelta, "apply_time_delta", "");
		//CommandLineString(strAppend, "append", "");

		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Create autocurator
	AutoCurator autocurator;

	// Check arguments
	if ((strInputNodeFile.length() == 0) && (strInputNodeFileList.length() == 0)) {
		_EXCEPTIONT("No input file (--in_nodefile) or (--in_nodefile_list)"
			" specified");
	}
	if ((strInputNodeFile.length() != 0) && (strInputNodeFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_nodefile) or (--in_nodefile_list)"
			" may be specified");
	}
/*
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
*/
	int nInputArguments =
		  ((strInputData.length() != 0)?(1):(0))
		+ ((strInputDataList.length() != 0)?(1):(0))
		+ ((strInputDataIndex.length() != 0)?(1):(0));
	if (nInputArguments > 1) {
		_EXCEPTIONT("Only one of (--in_data), (--in_data_list) or (--in_data_index)"
			" may be specified");
	}
	if ((strOutputFileFormat != "gfdl") &&
		(strOutputFileFormat != "csv") &&
		(strOutputFileFormat != "csvnohead")
	) {
		_EXCEPTIONT("Output format must be either \"gfdl\", \"csv\", or \"csvnohead\"");
	}

	// Input file type
	NodeFile::PathType iftype;
	if (strPathType == "DN") {
		iftype = NodeFile::PathTypeDN;
	} else if (strPathType == "SN") {
		iftype = NodeFile::PathTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_nodefile_type, expected \"SN\" or \"DN\"");
	}

	// NodeFile
	NodeFile nodefile;

	// Parse --in_fmt string
	ColumnDataHeader cdhInput;
	cdhInput.Parse(strInputFormat);

	// Parse --out_fmt string
	ColumnDataHeader cdhOutput;
	cdhOutput.Parse(strOutputFormat);

	// Parse --timefilter
	if (strTimeFilter == "3hr") {
		strTimeFilter = "....-..-.. (00|03|06|09|12|15|18|21):00:00";
	}
	if (strTimeFilter == "6hr") {
		strTimeFilter = "....-..-.. (00|06|12|18):00:00";
	}
	if (strTimeFilter == "daily") {
		strTimeFilter = "....-..-.. 00:00:00";
	}

#ifdef TEMPEST_NOREGEX
	if (strTimeFilter != "") {
		_EXCEPTIONT("Cannot use --timefilter with -DTEMPEST_NOREGEX compiler flag");
	}
#endif
#ifndef TEMPEST_NOREGEX
	std::regex reTimeSubset;
	if (strTimeFilter != "") {
		// Test regex support
		TestRegex();

		try {
			reTimeSubset.assign(strTimeFilter);
		} catch(std::regex_error & reerr) {
			_EXCEPTION2("Parse error in --timefilter regular expression \"%s\" (code %i)",
				strTimeFilter.c_str(), reerr.code());
		}
	}
#endif

	// Parse --calculate
	CalculationList calclist;
	if (strCalculate != "") {
		calclist.Parse(strCalculate);
	}

	// Parse --apply_time_delta
	bool fAddTimeDelta = true;
	Time timeDelta;
	if (strTimeDelta != "") {
		if (strTimeDelta[0] == '-') {
			fAddTimeDelta = false;
			strTimeDelta = strTimeDelta.substr(1);
		}
		timeDelta.FromFormattedString(strTimeDelta);
	}

	// Curate input data
	if (strInputData.length() != 0) {
		AnnounceStartBlock("Autocurating in_data");
		autocurator.IndexFiles(strInputData);

	} else if (strInputDataList.length() != 0) {
		AnnounceStartBlock("Autocurating in_data_list");
		std::ifstream ifInputDataList(strInputDataList.c_str());
		if (!ifInputDataList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputDataList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputDataList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			Announce(strFileLine.c_str());
			autocurator.IndexFiles(strFileLine);
		}

	} else if (strInputDataIndex.length() != 0) {
		AnnounceStartBlock("Reading autocurator index from \"%s\"", strInputDataIndex.c_str());
		autocurator.FromYAMLFile(strInputDataIndex);
	}
	AnnounceEndBlock("Done");

	// Define the SimpleGrid
	SimpleGrid grid;

	const std::vector<std::string> & vecFiles = autocurator.GetFilenames();

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No data files
	//} else if (vecFiles.size() < 1) {
	//	Announce("No data files specified; operations requiring a grid not supported");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");
		if (vecFiles.size() < 1) {
			_EXCEPTIONT("No data files specified -- unable to proceed without being able to determine grid dimensionality");
		}

		NcFile ncFile(vecFiles[0].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"", vecFiles[0].c_str());
		}

		grid.GenerateLatitudeLongitude(
			&ncFile,
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Load input file list
	std::vector<std::string> vecInputNodeFiles;

	if (strInputNodeFile.length() != 0) {
		vecInputNodeFiles.push_back(strInputNodeFile);

	} else {
		std::ifstream ifInputFileList(strInputNodeFileList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputNodeFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputNodeFiles.push_back(strFileLine);
		}
	}

	// A map from Times to file lines, used for StitchNodes formatted output
	TimeToPathNodeMap & mapTimeToPathNode = nodefile.GetTimeToPathNodeMap();

	// Vector of Path information
	PathVector & pathvec = nodefile.GetPathVector();

	// Column filters
	std::vector<FilterOp> vecFilterOp;

	// Loop over all files
	for (int f = 0; f < vecInputNodeFiles.size(); f++) {

		AnnounceStartBlock("Processing input (%s)", vecInputNodeFiles[f].c_str());

		// Read contents of NodeFile into PathVector
		AnnounceStartBlock("Reading file");
		nodefile.Read(
			vecInputNodeFiles[f],
			iftype,
			cdhInput,
			grid,
			autocurator.GetCalendarType());
		AnnounceEndBlock("Done");

		// If the nodefile contains negative grid indices, build the kd-tree
		int ixCDHLon = (-1);
		int ixCDHLat = (-1);
		if (nodefile.ContainsNegativeGridIx()) {
			AnnounceStartBlock("Negative grid indices found in nodefile");
			ixCDHLon = cdhInput.GetIndexFromString("lon");
			if (ixCDHLon == (-1)) {
				_EXCEPTIONT("Missing column header \"lon\", needed if negative grid indices present");
			}
			ixCDHLat = cdhInput.GetIndexFromString("lat");
			if (ixCDHLat == (-1)) {
				_EXCEPTIONT("Missing column header \"lat\", needed if negative grid indices present");
			}
			AnnounceStartBlock("Building kd-tree on grid");
			grid.BuildKDTree();
			AnnounceEndBlock("Done");
			AnnounceEndBlock(NULL);
		}

		// Time filter the nodefile
		if (strTimeFilter != "") {
#ifndef TEMPEST_NOREGEX
			for (int p = 0; p < nodefile.m_pathvec.size(); p++) {
				Path & path = nodefile.m_pathvec[p];
				for (int n = 0; n < path.size(); n++) {
					PathNode & pathnode = path[n];
					std::string strTime = pathnode.m_time.ToString();
					std::smatch match;
					if (!std::regex_search(strTime, match, reTimeSubset)) {
						path.erase(path.begin()+n);
						n--;
					}
				}
				if (path.size() == 0) {
					nodefile.m_pathvec.erase(nodefile.m_pathvec.begin()+p);
					p--;
				}
			}
#endif
		}

		// Parse --col_filter (need to do it here since cdh information
		// could be part of Read)
		if ((f == 0) && (strColumnFilter != "")) {
			AnnounceStartBlock("Parsing filter operations");

			int iLast = 0;
			for (int i = 0; i <= strColumnFilter.length(); i++) {

				if ((i == strColumnFilter.length()) ||
					(strColumnFilter[i] == ';') ||
					(strColumnFilter[i] == ':')
				) {
					std::string strSubStr =
						strColumnFilter.substr(iLast, i - iLast);
			
					int iNextOp = (int)(vecFilterOp.size());
					vecFilterOp.resize(iNextOp + 1);
					vecFilterOp[iNextOp].Parse(cdhInput, strSubStr);

					iLast = i + 1;
				}
			}
			AnnounceEndBlock("Done");
		}

		// Column filter the nodefile
		for (int op = 0; op < vecFilterOp.size(); op++) {
			for (int p = 0; p < nodefile.m_pathvec.size(); p++) {
				Path & path = nodefile.m_pathvec[p];
				for (int n = 0; n < path.size(); n++) {
					PathNode & pathnode = path[n];
					double dValue =
						pathnode.GetColumnDataAsDouble(
							vecFilterOp[op].m_iColumn);

					if (!vecFilterOp[op].Satisfies(dValue)) {
						path.erase(path.begin()+n);
						n--;
					}
				}
				if (path.size() == 0) {
					nodefile.m_pathvec.erase(nodefile.m_pathvec.begin()+p);
					p--;
				}
			}
		}

		// Generate the TimeToPathNodeMap
		nodefile.GenerateTimeToPathNodeMap();

		// Verify that all data is available for the nodefile
		{
			AnnounceStartBlock("Verifying that all data is available");
			TimeToPathNodeMap::iterator iterPathNode =
				mapTimeToPathNode.begin();
			for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
				const Time & time = iterPathNode->first;
				NcFileVector vecncDataFiles;
				autocurator.FindFilesAtTime(time, vecncDataFiles);
				if (vecncDataFiles.size() == 0) {
					_EXCEPTION1("Time (%s) does not exist in input data fileset",
						time.ToString().c_str());
				}
			}
			AnnounceEndBlock("Done");
		}

		// Working ColumnDataHeader
		ColumnDataHeader & cdhWorking = nodefile.m_cdh;

		// Perform calculations
		for (int i = 0; i < calclist.size(); i++) {
			const CalculationCommand & calccomm = calclist[i];
			AnnounceStartBlock("Calculating \"%s\"",
				calccomm.ToString().c_str());

			// Assignment operation
			if (calccomm.arg.size() == 0) {
				int ix = cdhWorking.GetIndexFromString(calccomm.rhs);
				if (ix == (-1)) {
					_EXCEPTION1("Unknown column header \"%s\"", calccomm.rhs.c_str());
				}
				nodefile.m_cdh.push_back(calccomm.lhs);
				pathvec.Duplicate(ix);
				AnnounceEndBlock("Done");
				continue;
			}

			// outputcmd
			if (calccomm.rhs == "outputcmd") {

				std::string strOutputCmd;
				if (calccomm.arg.size() == 3) {
					strOutputCmd = calccomm.arg[0] + "," + calccomm.arg[1] + "," + calccomm.arg[2];
				} else if (calccomm.arg.size() == 4) {
					strOutputCmd = calccomm.arg[0] + "(" + calccomm.arg[1] + ")," + calccomm.arg[2] + "," + calccomm.arg[3];
				} else {
					_EXCEPTIONT("Syntax error: Function \"outputcmd\" "
						"requires three arguments:\n"
						"outputcmd(<variable>, <op>, <dist>)");
				}

				NodeOutputOp opOutputCmd;
				opOutputCmd.Parse(varreg, strOutputCmd, false);

				std::string strResult;

				// Loop through all Times
				TimeToPathNodeMap::iterator iterPathNode =
					mapTimeToPathNode.begin();
				for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
					const Time & time = iterPathNode->first;

					Announce("%s", time.ToString().c_str());

					// Unload data from the VariableRegistry
					varreg.UnloadAllGridData();

					// Open NetCDF files with data at this time
					NcFileVector vecncDataFiles;
					autocurator.FindFilesAtTime(time, vecncDataFiles);
					if (vecncDataFiles.size() == 0) {
						_EXCEPTION1("Time (%s) does not exist in input data fileset",
							time.ToString().c_str());
					}

					// Loop through all PathNodes which need calculating at this Time
					for (int i = 0; i < iterPathNode->second.size(); i++) {
						int iPath = iterPathNode->second[i].first;
						int iPathNode = iterPathNode->second[i].second;

						PathNode & pathnode =
							pathvec[iPath][iPathNode];

						int ix0 = pathnode.m_gridix;
						if (ix0 < 0) {
							double dLonNodeDeg = pathnode.GetColumnDataAsDouble(ixCDHLon);
							double dLatNodeDeg = pathnode.GetColumnDataAsDouble(ixCDHLat);
							ix0 = static_cast<int>(
								grid.NearestNode(
									DegToRad(dLonNodeDeg),
									DegToRad(dLatNodeDeg)));
/*
							Announce("[%f,%f]->[%f,%f]",
								dLonNodeDeg,
								dLatNodeDeg,
								RadToDeg(grid.m_dLon[ix0]),
								RadToDeg(grid.m_dLat[ix0]));
*/
						}

						std::string strResult;

						ApplyNodeOutputOp<float>(
							opOutputCmd,
							grid,
							varreg,
							vecncDataFiles,
							time,
							ix0,
							strResult);

						// Add this data to the pathnode
						pathnode.PushColumnData(
							new ColumnDataString(strResult));
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;

			}

			// eval_ace, eval_ike, eval_pdi
			if ((calccomm.rhs == "eval_ace") ||
			    (calccomm.rhs == "eval_ike") ||
			    (calccomm.rhs == "eval_pdi") ||
				(calccomm.rhs == "eval_acepsl") ||
				(calccomm.rhs == "eval_min") ||
				(calccomm.rhs == "eval_max")
			) {

				// Cyclone metric
				CycloneMetric eCycloneMetric;
				if (calccomm.rhs == "eval_ace") {
					eCycloneMetric = CycloneMetric_ACE;
				} else if (calccomm.rhs == "eval_acepsl") {
					eCycloneMetric = CycloneMetric_ACEPSL;
				} else if (calccomm.rhs == "eval_ike") {
					eCycloneMetric = CycloneMetric_IKE;
				} else if (calccomm.rhs == "eval_pdi") {
					eCycloneMetric = CycloneMetric_PDI;
				} else if (calccomm.rhs == "eval_min") {
					eCycloneMetric = CycloneMetric_MIN;
				} else if (calccomm.rhs == "eval_max") {
					eCycloneMetric = CycloneMetric_MAX;
				} else {
					_EXCEPTION();
				}

				std::string strRadiusArg;

				if (eCycloneMetric == CycloneMetric_ACEPSL || eCycloneMetric == CycloneMetric_MIN || eCycloneMetric == CycloneMetric_MAX) {
					if (calccomm.arg.size() != 2) {
						_EXCEPTION2("Syntax error: Function \"%s\" "
							"requires two arguments:\n"
							"%s(<psl variable>, <radius>)",
							calccomm.rhs.c_str(),
							calccomm.rhs.c_str());
					}

					strRadiusArg = calccomm.arg[1];

				} else {
					if (calccomm.arg.size() != 3) {
						_EXCEPTION2("Syntax error: Function \"%s\" "
							"requires three arguments:\n"
							"%s(<u variable>, <v variable>, <radius>)",
							calccomm.rhs.c_str(),
							calccomm.rhs.c_str());
					}

					strRadiusArg = calccomm.arg[2];
				}

				// Parse zonal wind variable
				VariableIndex varixU = varreg.FindOrRegister(calccomm.arg[0]);

				// Parse meridional wind variable (if present)
				VariableIndex varixV = varixU;
				if (eCycloneMetric == CycloneMetric_ACE || eCycloneMetric == CycloneMetric_PDI || eCycloneMetric == CycloneMetric_IKE) {
					varixV = varreg.FindOrRegister(calccomm.arg[1]);
				}

				// Loop through all Times
				TimeToPathNodeMap::iterator iterPathNode =
					mapTimeToPathNode.begin();
				for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
					const Time & time = iterPathNode->first;

					// Unload data from the VariableRegistry
					varreg.UnloadAllGridData();

					// Open NetCDF files with data at this time
					NcFileVector vecncDataFiles;
					autocurator.FindFilesAtTime(time, vecncDataFiles);
					if (vecncDataFiles.size() == 0) {
						_EXCEPTION1("Time (%s) does not exist in input data fileset",
							time.ToString().c_str());
					}

					// Loop through all PathNodes which need calculating at this Time
					for (int i = 0; i < iterPathNode->second.size(); i++) {
						int iPath = iterPathNode->second[i].first;
						int iPathNode = iterPathNode->second[i].second;

						PathNode & pathnode =
							pathvec[iPath][iPathNode];

						CalculateCycloneMetrics(
							eCycloneMetric,
							varreg,
							vecncDataFiles,
							grid,
							cdhWorking,
							pathnode,
							varixU,
							varixV,
							strRadiusArg);
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}

			// radial_profile
			if ((calccomm.rhs == "radial_profile") || (calccomm.rhs == "radial_quadrant_profile")) {
				if (calccomm.arg.size() != 3) {
					_EXCEPTION2("Syntax error: Function \"%s\" "
						"requires three arguments:\n"
						"%s(<variable>, <# bins>, <bin width>)",
						calccomm.rhs.c_str(), calccomm.rhs.c_str());
				}

				bool fRadialProfile = (calccomm.rhs == "radial_profile");

				// Parse zonal wind variable
				VariableIndex varix = varreg.FindOrRegister(calccomm.arg[0]);

				// Loop through all Times
				TimeToPathNodeMap::iterator iterPathNode =
					mapTimeToPathNode.begin();
				for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
					const Time & time = iterPathNode->first;

					// Unload data from the VariableRegistry
					varreg.UnloadAllGridData();

					// Open NetCDF files with data at this time
					NcFileVector vecncDataFiles;
					autocurator.FindFilesAtTime(time, vecncDataFiles);
					if (vecncDataFiles.size() == 0) {
						_EXCEPTION1("Time (%s) does not exist in input data fileset",
							time.ToString().c_str());
					}

					// Loop through all PathNodes which need calculating at this Time
					for (int i = 0; i < iterPathNode->second.size(); i++) {
						int iPath = iterPathNode->second[i].first;
						int iPathNode = iterPathNode->second[i].second;

						PathNode & pathnode =
							pathvec[iPath][iPathNode];

						if (fRadialProfile) {
							CalculateRadialProfile(
								varreg,
								vecncDataFiles,
								grid,
								cdhWorking,
								pathnode,
								varix,
								calccomm.arg[1],
								calccomm.arg[2]);
	
						} else {
							CalculateQuadrantProfile(
								varreg,
								vecncDataFiles,
								grid,
								cdhWorking,
								pathnode,
								varix,
								calccomm.arg[1],
								calccomm.arg[2]);
						}
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}

			// radial_wind_profile
			if ((calccomm.rhs == "radial_wind_profile") || (calccomm.rhs == "quadrant_wind_profile")) {
				if (calccomm.arg.size() != 4) {
					_EXCEPTION2("Syntax error: Function \"%s\" "
						"requires four arguments:\n"
						"%s(<u variable>, <v variable>, <# bins>, <bin width>)",
						calccomm.rhs.c_str(), calccomm.rhs.c_str());
				}

				bool fRadialWindProfile = (calccomm.rhs == "radial_wind_profile");

				// Parse zonal wind variable
				VariableIndex varixU = varreg.FindOrRegister(calccomm.arg[0]);

				// Parse meridional wind variable
				VariableIndex varixV = varreg.FindOrRegister(calccomm.arg[1]);

				// Loop through all Times
				TimeToPathNodeMap::iterator iterPathNode =
					mapTimeToPathNode.begin();
				for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
					const Time & time = iterPathNode->first;

					// Unload data from the VariableRegistry
					varreg.UnloadAllGridData();

					// Open NetCDF files with data at this time
					NcFileVector vecncDataFiles;
					autocurator.FindFilesAtTime(time, vecncDataFiles);
					if (vecncDataFiles.size() == 0) {
						_EXCEPTION1("Time (%s) does not exist in input data fileset",
							time.ToString().c_str());
					}

					// Loop through all PathNodes which need calculating at this Time
					for (int i = 0; i < iterPathNode->second.size(); i++) {
						int iPath = iterPathNode->second[i].first;
						int iPathNode = iterPathNode->second[i].second;

						PathNode & pathnode =
							pathvec[iPath][iPathNode];

						if (fRadialWindProfile) {
							CalculateRadialWindProfile(
								varreg,
								vecncDataFiles,
								grid,
								cdhWorking,
								pathnode,
								varixU,
								varixV,
								calccomm.arg[2],
								calccomm.arg[3]);

						} else {
							CalculateQuadrantWindProfile(
								varreg,
								vecncDataFiles,
								grid,
								cdhWorking,
								pathnode,
								varixU,
								varixV,
								calccomm.arg[2],
								calccomm.arg[3]);
						}
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}

			// lastwhere
			if ((calccomm.rhs == "lastwhere") || (calccomm.rhs == "firstwhere")) {

				bool fFirstWhere = (calccomm.rhs == "firstwhere");

				if ((fFirstWhere) && (calccomm.arg.size() != 4)) {
					_EXCEPTIONT("Syntax error: Function \"firstwhere\" "
						"requires four arguments:\n"
						"firstwhere(<column name>, <op>, <value>, <radius>)");
				}

				if ((!fFirstWhere) && (calccomm.arg.size() != 3)) {
					_EXCEPTION2("Syntax error: Function \"%s\" "
						"requires three arguments:\n"
						"%s(<column name>, <op>, <value>)",
						calccomm.rhs.c_str(), calccomm.rhs.c_str());
				}

				// Get arguments
				int ixCol = cdhWorking.GetIndexFromString(calccomm.arg[0]);
				if (ixCol == (-1)) {
					_EXCEPTION1("Invalid column header \"%s\"", calccomm.arg[0].c_str());
				}

				const std::string & strOp = calccomm.arg[1];

				const std::string & strThreshold = calccomm.arg[2];

				std::string strRadiusArg = "";
				if (fFirstWhere) {
					strRadiusArg = calccomm.arg[3];
				} 

				// Loop through all PathNodes
				for (int p = 0; p < pathvec.size(); p++) {
					Path & path = pathvec[p];

					for (int i = 0; i < pathvec[p].size(); i++) {
						PathNode & pathnode = path[i];

						// Process a ColumnDataDoubleArrayTemplate
						ColumnDataDoubleArrayTemplate * pdatDA =
							dynamic_cast<ColumnDataDoubleArrayTemplate *>(
								pathnode.m_vecColumnData[ixCol]);

						if (pdatDA != NULL) {
							const std::vector<double> & dIndices = pdatDA->GetIndices();
							const std::vector<double> & dArray = pdatDA->GetValues();

							double dIndex;
							if (fFirstWhere) {
								dIndex = FirstWhere(cdhWorking, pathnode, strOp, strThreshold, strRadiusArg, dIndices, dArray);
							} else {
								dIndex = LastWhere(cdhWorking, pathnode, strOp, strThreshold, dIndices, dArray);
							}

							// Add this data to the pathnode
							pathnode.PushColumnData(new ColumnDataDouble(dIndex));

							continue;
						}

						// Process a ColumnData2DArray
						ColumnData2DArray * pdat2DA =
							dynamic_cast<ColumnData2DArray *>(
								pathnode.m_vecColumnData[ixCol]);

						if (pdat2DA != NULL) {
							if (pdat2DA->m_dValues.size() < 2) {
								_EXCEPTION1("2DArray \"%s\" must have at least two entries",
									(*pargfunc)[0].c_str());
							}
							ColumnData1DArray * pdatOut = new ColumnData1DArray;

							const std::vector<double> & dIndices = pdat2DA->m_dValues[0];
							for (int i = 1; i < pdat2DA->m_dValues.size(); i++) {
								const std::vector<double> & dArray = pdat2DA->m_dValues[i];

								double dIndex;
								if (fFirstWhere) {
									dIndex = FirstWhere(cdhWorking, pathnode, strOp, strThreshold, strRadiusArg, dIndices, dArray);
								} else {
									dIndex = LastWhere(cdhWorking, pathnode, strOp, strThreshold, dIndices, dArray);
								}

								pdatOut->m_dValues.push_back(dIndex);
							}
							pathnode.PushColumnData(pdatOut);
							continue;
						}

						_EXCEPTION1("Cannot cast \"%s\" to DoubleArray or 2DArray type",
							calccomm.arg[0].c_str());
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}

			// Extract the value of an array at a given radius
			if (calccomm.rhs == "value") {
				if (calccomm.arg.size() != 2) {
					_EXCEPTIONT("Syntax error: Function \"value\" "
						"requires two arguments:\n"
						"value(<column name>, <index>)");
				}

				// Get arguments
				int ix = cdhWorking.GetIndexFromString(calccomm.arg[0]);
				if (ix == (-1)) {
					_EXCEPTION1("Invalid column header \"%s\"", calccomm.arg[0].c_str());
				}

				const std::string & strIndex = calccomm.arg[1];

				// Loop through all PathNodes
				for (int p = 0; p < pathvec.size(); p++) {
					Path & path = pathvec[p];

					for (int i = 0; i < pathvec[p].size(); i++) {
						PathNode & pathnode = path[i];

						ColumnDataDoubleArrayTemplate * pdat =
							dynamic_cast<ColumnDataDoubleArrayTemplate *>(
								pathnode.m_vecColumnData[ix]);

						if (pdat == NULL) {
							_EXCEPTION1("Cannot cast \"%s\" to DoubleArray type",
								calccomm.arg[0].c_str());
						}

						const std::vector<double> & dR = pdat->GetIndices();
						const std::vector<double> & dUa = pdat->GetValues();

						if (dR.size() == 0) {
							_EXCEPTIONT("PathNode radial profile has zero size");
						}
						if (dR.size() != dUa.size()) {
							_EXCEPTIONT("PathNode R array size different from Ua size");
						}

						double dIndex =
							pathnode.GetColumnDataAsDouble(cdhWorking, strIndex);

						if (dIndex < 0.0) {
							_EXCEPTION1("Negative index value (%3.6f) found in call to value()",
								dIndex);
						}

						// Extract the value using linear interpolation
						double dValue;
						bool fIndexOutOfRange = true;
						for (int j = 1; j < dR.size(); j++) {
							if (dR[j] > dIndex) {
								if (dR[j] == dR[j-1]) {
									dValue = 0.5 * (dUa[j] + dUa[j-1]);
								} else {
									dValue = (
										dUa[j-1] * (dR[j] - dIndex)
										+ dUa[j] * (dIndex - dR[j-1])
										) / (dR[j] - dR[j-1]);
								}
								fIndexOutOfRange = false;
								break;
							}
						}
						if (fIndexOutOfRange) {
							dValue = dUa[dUa.size()-1];
						}

						// Add this data to the pathnode
						pathnode.PushColumnData(
							new ColumnDataDouble(dValue));
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}

			// max_closed_contour_delta
			if ((calccomm.rhs == "max_closed_contour_delta") || (calccomm.rhs == "max_closed_contour_value")) {
				if ((calccomm.arg.size() < 2) && (calccomm.arg.size() > 3)) {
					_EXCEPTION3("Syntax error: Function \"%s\" "
						"requires two or three arguments:\n"
						"%s(<variable>, <radius>)\n"
						"%s(<variable>, <radius>, <index>)",
						calccomm.rhs.c_str(), calccomm.rhs.c_str(), calccomm.rhs.c_str());
				}

				// Add the root value back into the max_closed_contour_delta
				bool fAddRootValue = (calccomm.rhs == "max_closed_contour_value");

				// Parse variable
				VariableIndex varix = varreg.FindOrRegister(calccomm.arg[0]);

				// Radius
				std::string strRadius(calccomm.arg[1]);
				
				// Index
				std::string strIndex("");
				if (calccomm.arg.size() == 3) {
					strIndex = calccomm.arg[2];
				}

				// Loop through all Times
				TimeToPathNodeMap::iterator iterPathNode =
					mapTimeToPathNode.begin();
				for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
					const Time & time = iterPathNode->first;

					// Unload data from the VariableRegistry
					varreg.UnloadAllGridData();

					// Open NetCDF files with data at this time
					NcFileVector vecncDataFiles;
					autocurator.FindFilesAtTime(time, vecncDataFiles);
					if (vecncDataFiles.size() == 0) {
						_EXCEPTION1("Time (%s) does not exist in input data fileset",
							time.ToString().c_str());
					}

					// Loop through all PathNodes which need calculating at this Time
					for (int i = 0; i < iterPathNode->second.size(); i++) {
						int iPath = iterPathNode->second[i].first;
						int iPathNode = iterPathNode->second[i].second;

						PathNode & pathnode =
							pathvec[iPath][iPathNode];

						MaxClosedContourDelta(
							varreg,
							vecncDataFiles,
							grid,
							cdhWorking,
							pathnode,
							varix,
							strRadius,
							strIndex,
							fAddRootValue);
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}

			// region_name
			if (calccomm.rhs == "region_name") {
				if (calccomm.arg.size() != 1) {
					_EXCEPTIONT("Syntax error: Function \"region_name\" "
						"requires one argument:\n"
						"region_name(<filename>)");
				}

				std::string strFilename = calccomm.arg[0];
				if (strFilename[0] == '\"') {
					if ((strFilename.length() == 1) ||
						(strFilename[strFilename.length()-1] != '\"')
					) {
						_EXCEPTION1("Unterminated quotation mark in filename \"%s\"",
							strFilename.c_str());
					}
					strFilename = strFilename.substr(1,strFilename.length()-2);
				}

				RLLPolygonArray rllpolyarray;

				rllpolyarray.FromFile(strFilename);

				// Loop through all PathNodes
				for (int p = 0; p < pathvec.size(); p++) {
					Path & path = pathvec[p];

					for (int i = 0; i < pathvec[p].size(); i++) {
						PathNode & pathnode = path[i];

						int ix0 = pathnode.m_gridix;

						RLLPoint pt;
						pt.lon = grid.m_dLon[ix0] * 180.0 / M_PI;
						pt.lat = grid.m_dLat[ix0] * 180.0 / M_PI;

						pathnode.PushColumnData(
							new ColumnDataString(
								rllpolyarray.NameOfRegionContainingPoint(pt)));

					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				continue;
			}
/*
			// min,max
			if ((calccomm.rhs == "min") ||
			    (calccomm.rhs == "max")
			) {
				if (calccomm.arg.size() < 2) {
					_EXCEPTION1("Syntax error: Function \"%s\" requires at least two arguments",
						(*pargtree)[2].c_str());
				}

				for (int i = 0; i < pargtree->size(); i++) {
					int ix = cdhWorking.GetIndexFromString((*pargfunc)[0]);
					if (ix == (-1)) {
						_EXCEPTION1("Invalid column header \"%s\"", (*pargfunc)[0].c_str());
					}
				}

				continue;
			}
*/
/*
			// sum_radius
			if (calccomm.rhs == "sum_radius") {
				if (calccomm.arg.size() != 2) {
					_EXCEPTIONT("Syntax error: Function \"sum_radius\" "
						"requires two arguments:\n"
						"lastwhere(<field>, <radius>)");
				}

				// Parse variable
				VariableIndex varix = varreg.FindOrRegister(calccomm.arg[0]);

				// Loop through all Times
				TimeToPathNodeMap::iterator iterPathNode =
					mapTimeToPathNode.begin();
				for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
					const Time & time = iterPathNode->first;

					// Unload data from the VariableRegistry
					varreg.UnloadAllGridData();

					// Open NetCDF files with data at this time
					NcFileVector vecncDataFiles;
					autocurator.Find(time, vecncDataFiles);
					if (vecncDataFiles.size() == 0) {
						_EXCEPTION1("Time (%s) does not exist in input data fileset",
							time.ToString().c_str());
					}

					// Loop through all PathNodes which need calculating at this Time
					for (int i = 0; i < iterPathNode->second.size(); i++) {
						int iPath = iterPathNode->second[i].first;
						int iPathNode = iterPathNode->second[i].second;

						PathNode & pathnode =
							pathvec[iPath][iPathNode];

						SumRadius(
							varreg,
							vecncDataFiles,
							grid,
							cdhWorking,
							pathnode,
							varix,
							calccomm.arg[1]);
					}
				}

				// Add new variable to ColumnDataHeader
				cdhWorking.push_back(calccomm.lhs);

				AnnounceEndBlock("Done");
				continue;
			}
*/
			// Unknown operation
			Announce(NULL);
			AnnounceEndBlock("WARNING: Unknown function \"%s\" no operation performed",
				calccomm.rhs.c_str());

		}

		// Apply time delta
		if (strTimeDelta != "") {
			nodefile.ApplyTimeDelta(timeDelta, fAddTimeDelta);
		}

		// Output
		if (strOutputFile != "") {
			std::vector<int> vecColumnDataOutIx;
			cdhWorking.GetIndicesFromColumnDataHeader(cdhOutput, vecColumnDataOutIx);

			AnnounceStartBlock("Writing output");

			NodeFile::FileFormat nodefileformat;
			if (strOutputFileFormat == "gfdl") {
				nodefile.Write(
					strOutputFile,
					&grid,
					&vecColumnDataOutIx,
					NodeFile::FileFormatGFDL);

			} else if (strOutputFileFormat == "csv") {
				nodefile.Write(
					strOutputFile,
					&grid,
					&vecColumnDataOutIx,
					NodeFile::FileFormatCSV,
					true);

			} else if (strOutputFileFormat == "csvnohead") {
				nodefile.Write(
					strOutputFile,
					&grid,
					&vecColumnDataOutIx,
					NodeFile::FileFormatCSV,
					false);

			} else {
				_EXCEPTION();
			}

			AnnounceEndBlock("Done");
		}
		AnnounceEndBlock("Done");
	}

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
/*
#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
*/
}

