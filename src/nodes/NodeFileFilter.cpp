///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileFilter.cpp
///	\author  Paul Ullrich
///	\version October 4, 2018
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

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "Variable.h"
#include "AutoCurator.h"
#include "DataArray2D.h"
#include "ArgumentTree.h"
#include "STLStringHelper.h"
#include "NodeFileUtilities.h"
#include "NetCDFUtilities.h"
#include "ClosedContourOp.h"
#include "SimpleGridUtilities.h"
#include "GridElements.h"

#include "netcdfcpp.h"

#include <fstream>
#include <queue>
#include <set>
#include <cmath>
/*
#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif
*/

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a nearbyblob operator.
///	</summary>
class NearbyBlobsOp {

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
	NearbyBlobsOp() :
		m_varix(InvalidVariableIndex),
		m_dDistance(0.0),
		m_eOp(GreaterThan),
		m_dValue(0.0),
		m_dMaxDist(180.0)
	{ }

public:
	///	<summary>
	///		Parse a nearbyblob operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Distance,
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_MaxDist,
			ReadMode_Invalid
		} eReadMode = ReadMode_Distance;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delimited
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in distance
				if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

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
					eReadMode = ReadMode_MaxDist;

				// Read in maximum distance
				} else if (eReadMode == ReadMode_MaxDist) {
					m_dMaxDist = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
							"\nRequired: \"<variable>,<distance>,<operation>,<value>[,<max distance>]\"",
							strOp.c_str());
				}
			}
		}

		if ((eReadMode != ReadMode_Invalid) && (eReadMode != ReadMode_MaxDist)) {
			_EXCEPTION1("\nInsufficient entries in --nearbyblobs op \"%s\""
					"\nRequired: \"<name>,<distance>,<operation>,<value>[,<max distance>]\"",
					strOp.c_str());
		}

		if (m_dDistance < 0.0) {
			_EXCEPTION1("For --nearbyblobs, distance (%2.6f) must be nonnegative", m_dDistance);
		}
		if (m_dMaxDist < 0.0) {
			_EXCEPTION1("For --nearbyblobs, max distance (%2.6f) must be nonnegative", m_dMaxDist);
		}
		if (m_dDistance > m_dMaxDist) {
			_EXCEPTION2("For --nearbyblobs, distance (%2.6f) must be less than or equal to max distance (%2.6f)", m_dDistance, m_dMaxDist);
		}

		// Output announcement
		char szBuffer[128];

		sprintf(szBuffer, "%f", m_dDistance);
		std::string strDescription =
			std::string("Include blobs with at least 1 point within ") + szBuffer
			+ std::string(" degrees where ") + varreg.GetVariableString(m_varix);
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

		if (fabs(m_dValue) < 1.0e-4) {
			sprintf(szBuffer, "%e", m_dValue);
		} else {
			sprintf(szBuffer, "%f", m_dValue);
		}
		strDescription += szBuffer;

		sprintf(szBuffer, "%f", m_dMaxDist);
		strDescription += std::string(" (max dist ") + szBuffer + " degrees)";

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Determine if a particular value satisfies this threshold.
	///	</summary>
	bool SatisfiedBy(
		double dValue
	) const {
		if (m_eOp == GreaterThan) {
			if (dValue > m_dValue) {
				return true;
			}

		} else if (m_eOp == LessThan) {
			if (dValue < m_dValue) {
				return true;
			}

		} else if (m_eOp == GreaterThanEqualTo) {
			if (dValue >= m_dValue) {
				return true;
			}

		} else if (m_eOp == LessThanEqualTo) {
			if (dValue <= m_dValue) {
				return true;
			}

		} else if (m_eOp == EqualTo) {
			if (dValue == m_dValue) {
				return true;	
			}

		} else if (m_eOp == NotEqualTo) {
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
	///		Variable to use for nearbyblobs.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Distance to search for blob.
	///	</summary>
	double m_dDistance;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Maximum distance.
	///	</summary>
	double m_dMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

void BuildMask_ByDist(
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	const PathVector & pathvec,
	const PathNodeIndexVector & vecPathNodes,
	const std::string & strDist,
	DataArray1D<double> & dataMask
) {
	// Get filter width (either fixed value or column data header)
	bool fFixedFilterWidth = false;
	double dFilterWidth = 0.0;
	int iFilterWidthIx = 0;
	if (STLStringHelper::IsFloat(strDist)) {
		fFixedFilterWidth = true;
		dFilterWidth = atof(strDist.c_str());
		if (dFilterWidth == 0.0) {
			return;
		}

	} else {
		iFilterWidthIx = cdh.GetIndexFromString(strDist);
		if (iFilterWidthIx == (-1)) {
			_EXCEPTION1("Unknown column header \"%s\"", strDist.c_str());
		}
	}

	// Loop through all PathNodes
	for (int j = 0; j < vecPathNodes.size(); j++) {
		const Path & path = pathvec[vecPathNodes[j].first];
		const PathNode & pathnode = path[vecPathNodes[j].second];

		int ixOrigin = static_cast<int>(pathnode.m_gridix);

		// Extract the filter width for this PathNode from ColumnData
		if (!fFixedFilterWidth) {
			if ((iFilterWidthIx < 0) ||
			    (iFilterWidthIx > pathnode.m_vecColumnData.size())
			) {
				_EXCEPTIONT("Logic error");
			}

			std::string str = pathnode.m_vecColumnData[iFilterWidthIx]->ToString();
			if (!STLStringHelper::IsFloat(str)) {
				_EXCEPTION1("Column header \"%s\" cannot be cast to type double",
					strDist.c_str());
			}
			dFilterWidth = atof(str.c_str());

			if (dFilterWidth == 0.0) {
				continue;
			}
		}

		// Set of visited nodes
		std::set<int> setNodesVisited;

		// Set of nodes to visit
		std::queue<int> queueToVisit;
		queueToVisit.push(ixOrigin);

		const double dLat0 = grid.m_dLat[ixOrigin];
		const double dLon0 = grid.m_dLon[ixOrigin];

		// Using the grid connectivity find nodes within a specified
		// distance of each PathNode.
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

			// Check great circle distance
			if (dR > dFilterWidth) {
				continue;
			}

			// Add this point to the filter
			dataMask[ix] = 1.0;

			// Add all neighbors of this point
			for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
				queueToVisit.push(grid.m_vecConnectivity[ix][n]);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void BuildMask_ByContour(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	const ColumnDataHeader & cdh,
	const PathVector & pathvec,
	const PathNodeIndexVector & vecPathNodes,
	const ClosedContourOp & op,
	DataArray1D<double> & dataMask
) {
	// Get the variable
	double dDeltaAmt = op.m_dDeltaAmount;
	double dDeltaDist = op.m_dDistance;
	double dMinMaxDist = op.m_dMinMaxDist;

	_ASSERT(dDeltaAmt != 0.0);
	_ASSERT(dDeltaDist > 0.0);

	// Loop through all PathNodes
	for (int j = 0; j < vecPathNodes.size(); j++) {
		const Path & path = pathvec[vecPathNodes[j].first];
		const PathNode & pathnode = path[vecPathNodes[j].second];

		int ix0 = static_cast<int>(pathnode.m_gridix);

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

			// Check great circle distance
			if (dR > dDeltaDist) {
				continue;
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

			// Tag this point
			dataMask[ix] = 1.0;

			// Add all neighbors of this point
			for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
				queueToVisit.push(grid.m_vecConnectivity[ix][n]);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void BuildMask_NearbyBlobs(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	const ColumnDataHeader & cdh,
	const PathVector & pathvec,
	const PathNodeIndexVector & vecPathNodes,
	const NearbyBlobsOp & nearbyblobsop,
	DataArray1D<double> & dataMask
) {
	// Get the variable
	const double dDist = nearbyblobsop.m_dDistance;
	const double dMaxDist = nearbyblobsop.m_dMaxDist;

	_ASSERT(dataState.GetRows() == grid.GetSize());
	_ASSERT(dDist >= 0.0);
	_ASSERT(dDist <= 180.0);
	_ASSERT(dMaxDist >= dDist);
	_ASSERT(dMaxDist <= 180.0);

	// Loop through all PathNodes
	for (int j = 0; j < vecPathNodes.size(); j++) {
		const Path & path = pathvec[vecPathNodes[j].first];
		const PathNode & pathnode = path[vecPathNodes[j].second];

		int ix0 = static_cast<int>(pathnode.m_gridix);
		_ASSERT((ix0 >= 0) && (ix0 < grid.GetSize()));

		// Queue of nodes that remain to be visited
		std::queue<int> queueNodes;
		queueNodes.push(ix0);

		// Set of nodes that have already been visited
		std::set<int> setNodesVisited;

		// Latitude and longitude at the origin
		const double dLat0 = grid.m_dLat[ix0];
		const double dLon0 = grid.m_dLon[ix0];

		// Loop through all elements
		while (queueNodes.size() != 0) {
			int ix = queueNodes.front();
			queueNodes.pop();

			if (setNodesVisited.find(ix) != setNodesVisited.end()) {
				continue;
			}

			setNodesVisited.insert(ix);

			// Great circle distance to this dof
			_ASSERT((ix >= 0) && (ix < grid.GetSize()));

			double dR =
				GreatCircleDistance_Deg(
					dLon0, dLat0,
					grid.m_dLon[ix], grid.m_dLat[ix]);

			// Check if we have exceeded the maximum distance to find blobs
			if ((ix != ix0) && (dR > dDist)) {
				continue;
			}

			// Add all neighbors of this point
			for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
				queueNodes.push(grid.m_vecConnectivity[ix][n]);
			}

			// Check if this point satisfies the nearbyblobs criteria
			if (!nearbyblobsop.SatisfiedBy(static_cast<double>(dataState[ix]))) {
				continue;
			}

			// Tag this point
			dataMask[ix] = 1.0;

			// Operator satisfied; add all points in this blob to mask up to maxdist
			std::queue<int> queueThresholdedNodes;
			queueThresholdedNodes.push(ix);

			std::set<int> setThresholdNodesVisited;

			while (queueThresholdedNodes.size() != 0) {
				int ixblob = queueThresholdedNodes.front();
				queueThresholdedNodes.pop();

				if (setThresholdNodesVisited.find(ixblob) !=
				    setThresholdNodesVisited.end()
				) {
					continue;
				}

				setThresholdNodesVisited.insert(ixblob);

				// Make sure great circle distance to this dof is closer than dMaxDist
				_ASSERT((ixblob >= 0) && (ixblob < grid.GetSize()));

				double dRblob =
					GreatCircleDistance_Deg(
						dLon0, dLat0,
						grid.m_dLon[ixblob], grid.m_dLat[ixblob]);

				if ((ix != ixblob) && (dRblob > dMaxDist)) {
					continue;
				}

				// Verify this point satisfies the condition
				if (!nearbyblobsop.SatisfiedBy(static_cast<double>(dataState[ixblob]))) {

					// Isn't part of the blob, but add it to the list of
					// nodes to visit.
					if (dRblob <= dDist) {
						queueNodes.push(ixblob);
					}
					continue;
				}

				// Add this point to the set of visited points to avoid it
				// again triggering a blob search.
				setNodesVisited.insert(ixblob);

				// Add all neighbors of this point to search
				for (int n = 0; n < grid.m_vecConnectivity[ixblob].size(); n++) {
					queueThresholdedNodes.push(grid.m_vecConnectivity[ixblob][n]);
				}

				dataMask[ixblob] = 1.0;
			}
		}
	}
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

	// Input nodefile
	std::string strInputNodeFile;

	// Input nodefile type
	std::string strInputNodeFileType;

	// Input format (columns of in_file)
	std::string strInputFormat;

	// Input data file
	std::string strInputData;

	// Input list of data files
	std::string strInputDataList;

	// Connectivity file
	std::string strConnectivity;

	// Data is regional
	bool fRegional;

	// Output data file
	std::string strOutputData;

	// Output list of data files
	std::string strOutputDataList;

	// List of variables to output
	std::string strVariables;

	// Output the mask to the given variable
	std::string strMaskVariable;

	// List of variables to preserve
	std::string strPreserve;

	// Fill value for areas outside of mask
	std::string strFillValue;

	// Filter variables by distance
	std::string strFilterByDist;

	// Filter variables by contour
	std::string strFilterByContour;

	// Detect nearby blobs
	std::string strNearbyBlobs;

	// Apply the inverted filter
	bool fInvert;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_nodefile", "");
		CommandLineStringD(strInputNodeFileType, "in_nodefile_type", "SN", "[DCU|SN]");
		CommandLineString(strInputFormat, "in_fmt", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputDataList, "out_data_list", "");

		CommandLineString(strVariables, "var", "");
		CommandLineString(strMaskVariable, "maskvar", "");
		CommandLineString(strPreserve, "preserve", "");
		CommandLineStringD(strFillValue, "fillvalue", "", "[<value>|nan|att]");

		CommandLineStringD(strFilterByDist, "bydist", "", "[dist]");
		CommandLineStringD(strFilterByContour, "bycontour", "", "[var,delta,dist,minmaxdist]");
		CommandLineStringD(strNearbyBlobs, "nearbyblobs", "", "[var,dist,op,value[,maxdist]]");
		CommandLineBool(fInvert, "invert");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Check arguments
	if (strInputNodeFile.length() == 0) {
		_EXCEPTIONT("No input file (--in_nodefile) specified");
	}
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if ((strOutputData.length() == 0) && (strOutputDataList.length() == 0)) {
		_EXCEPTIONT("No output file (--out_data) or (--out_data_list)"
			" specified");
	}
	if ((strOutputData.length() != 0) && (strOutputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out_data) or (--out_data_list)"
			" may be specified");
	}
	if ((strFilterByDist.length() == 0) &&
		(strFilterByContour.length() == 0) &&
		(strNearbyBlobs.length() == 0)
	) {
		_EXCEPTIONT("No command (--bydist, --bycontour, or --nearbyblobs) specified");
	}
	if ((strFilterByDist.length() != 0) && (strFilterByContour.length() != 0)) {
		_EXCEPTIONT("Only one filter command (--bydist or --bycontour) may be specified");
	}
	if ((strVariables.length() == 0) && (strMaskVariable.length() == 0)) {
		_EXCEPTIONT("One of (--var) or (--maskvar) must be specified");
	}

	bool fHasFillValue = false;
	float dFillValue = 0.0;
	if (strFillValue != "") {
		fHasFillValue = true;
		STLStringHelper::ToLower(strFillValue);
		if (strFillValue == "att") {
		} else if (strFillValue == "nan") {
			dFillValue = nanf(NULL);
		} else if (STLStringHelper::IsFloat(strFillValue)) {
			dFillValue = std::stof(strFillValue);
		} else {
			_EXCEPTIONT("Invalid value specified for --fillvalue");
		}
	}

	// Input file type
	NodeFile::PathType iftype;
	if (strInputNodeFileType == "DCU") {
		iftype = NodeFile::PathTypeDCU;
	} else if (strInputNodeFileType == "SN") {
		iftype = NodeFile::PathTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_nodefile_type, expected \"SN\" or \"DCU\"");
	}

	// NodeFile
	NodeFile nodefile;

	// Parse --in_fmt string
	ColumnDataHeader cdhInput;
	cdhInput.Parse(strInputFormat);

	// Parse --var argument
	std::vector<std::string> vecVarNames;
	std::vector<VariableIndex> vecVarIx;
	if (strVariables != "") {
		std::string strVariablesTemp = strVariables;
		for (;;) {
			VariableIndex varix;
			int iLast = varreg.FindOrRegisterSubStr(strVariablesTemp, &varix) + 1;

			vecVarIx.push_back(varix);
			vecVarNames.push_back( strVariablesTemp.substr(0,iLast-1) );
			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}
	}
	_ASSERT(vecVarNames.size() == vecVarIx.size());

	// Parse --bydist
	if (strFilterByDist != "") {
		if (!STLStringHelper::IsFloat(strFilterByDist)) {
			bool fFound = false;
			for (int i = 0; i < cdhInput.size(); i++) {
				if (cdhInput[i] == strFilterByDist) {
					fFound = true;
					break;
				}
			}
			if (!fFound) {
				_EXCEPTION1("Invalid format of --bydist;"
					" \"%s\" does not appear in --in_fmt",
					strFilterByDist.c_str());
			}
		}
	}

	// Parse --bycontour
	std::vector<ClosedContourOp> vecClosedContourOp;

	if (strFilterByContour != "") {
		AnnounceStartBlock("Parsing --byclosedcontour");
		ClosedContourOp op;
		op.Parse(varreg, strFilterByContour);
		vecClosedContourOp.push_back(op);
		AnnounceEndBlock(NULL);
	}

	// Parse --nearbyblobs
	std::vector<NearbyBlobsOp> vecNearbyBlobsOp;

	if (strNearbyBlobs != "") {
		AnnounceStartBlock("Parsing --nearbyblobs");
		NearbyBlobsOp op;
		op.Parse(varreg, strNearbyBlobs);
		vecNearbyBlobsOp.push_back(op);
		AnnounceEndBlock(NULL);
	}

	// Parse variable preservation list
	std::vector< std::string > vecPreserveVariables;
	STLStringHelper::ParseVariableList(strPreserve, vecPreserveVariables);

	// Define the SimpleGrid
	SimpleGrid grid;

	// Store input data
	std::vector<std::string> vecInputFileList;

	if (strInputData.length() != 0) {
		vecInputFileList.push_back(strInputData);

	} else {
		AnnounceStartBlock("Building input data list");
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
			for (int i = 0; i < strFileLine.length(); i++) {
				if (strFileLine[i] == ';') {
					_EXCEPTIONT("Only one filename allowed per line in --in_data_list");
				}
			}
			vecInputFileList.push_back(strFileLine);
		}
		AnnounceEndBlock("Done");
	}

	// Store output data
	std::vector<std::string> vecOutputFileList;

	if (strOutputData.length() != 0) {
		vecOutputFileList.push_back(strOutputData);

	} else {
		AnnounceStartBlock("Building output data list");
		std::ifstream ifOutputDataList(strOutputDataList.c_str());
		if (!ifOutputDataList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strOutputDataList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputDataList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			for (int i = 0; i < strFileLine.length(); i++) {
				if (strFileLine[i] == ';') {
					_EXCEPTIONT("Only one filename allowed per line in --out_data_list");
				}
			}
			vecOutputFileList.push_back(strFileLine);
		}
		AnnounceEndBlock("Done");
	}

	// Verify both input data and output data have the same length
	if (vecInputFileList.size() != vecOutputFileList.size()) {
		_EXCEPTION2("Mismatch in --in_data_list (%lu files) and --out_data_list length (%lu files)",
			vecInputFileList.size(), vecOutputFileList.size());
	}

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		if (vecInputFileList.size() < 1) {
			_EXCEPTIONT("No data files specified; unable to generate grid");
		}

		NcFile ncFile(vecInputFileList[0].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"", vecInputFileList[0].c_str());
		}

		grid.GenerateLatitudeLongitude(&ncFile, fRegional);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Get the CalendarType
	Time::CalendarType caltype;
	{
		if (vecInputFileList.size() == 0) {
			_EXCEPTION();
		}
		NcFile ncfile(vecInputFileList[0].c_str(), NcFile::ReadOnly);
		if (!ncfile.is_valid()) {
			_EXCEPTION1("Unable to open input datafile \"%s\"",
				vecInputFileList[0].c_str());
		}

		NcVar * varTime = ncfile.get_var("time");
		if (varTime == NULL) {
			_EXCEPTION1("Missing \"time\" variable in datafile \"%s\"",
				vecInputFileList[0].c_str());
		}

		NcAtt * attCalendar = varTime->get_att("calendar");
		if (attCalendar == NULL) {
			_EXCEPTION1("Missing \"calendar\" in datafile \"%s\"",
				vecInputFileList[0].c_str());
		}

		caltype = Time::CalendarTypeFromString(attCalendar->as_string(0));
		if (caltype == Time::CalendarUnknown) {
			_EXCEPTION1("Invalid \"calendar\" in datafile \"%s\"",
				vecInputFileList[0].c_str());
		}
	}

	// Parse NodeFile
	PathVector & pathvec = nodefile.GetPathVector();
	TimeToPathNodeMap & mapTimeToPathNode = nodefile.GetTimeToPathNodeMap();

	AnnounceStartBlock("Parsing nodefile");
	nodefile.Read(
		strInputNodeFile,
		iftype,
		cdhInput,
		grid,
		caltype);
	AnnounceEndBlock("Done");

	// Create data array
	DataArray1D<double> data(grid.GetSize());

	// Create mask
	DataArray1D<double> dataMask(grid.GetSize());

	// Loop over all files
	_ASSERT(vecInputFileList.size() == vecOutputFileList.size());

	for (int f = 0; f < vecInputFileList.size(); f++) {

		AnnounceStartBlock("Processing input (%s)", vecInputFileList[f].c_str());

		// Open input file
		NcFileVector vecInFiles;
		vecInFiles.ParseFromString(vecInputFileList[f]);
		_ASSERT(vecInFiles.size() > 0);

		NcFile & ncinfile = *(vecInFiles[0]);
		//if (!ncinfile.is_valid()) {
		//	_EXCEPTIONT("Unable to open input datafile");
		//}

		// Get time variable and verify calendar is compatible
		NcVar * varTime = ncinfile.get_var("time");
		if (varTime == NULL) {
			_EXCEPTION1("Missing \"time\" variable in datafile \"%s\"",
				vecInputFileList[f].c_str());
		}

		NcAtt * attCalendar = varTime->get_att("calendar");
		if (attCalendar == NULL) {
			_EXCEPTION1("Missing \"calendar\" metadata in datafile \"%s\"",
				vecInputFileList[f].c_str());
		}

		Time::CalendarType caltypeThis =
			Time::CalendarTypeFromString(attCalendar->as_string(0));
		if (caltype != caltypeThis) {
			_EXCEPTION1("\"calendar\" mismatch in datafile \"%s\" ",
				vecInputFileList[f].c_str());
		}

		// Load in vector of Times
		std::vector<Time> vecTimes;
		ReadCFTimeDataFromNcFile(&ncinfile, vecInputFileList[f], vecTimes, true);

		// Open output file
		NcFile ncoutfile(vecOutputFileList[f].c_str(), NcFile::Replace);
		if (!ncoutfile.is_valid()) {
			_EXCEPTION1("Unable to open output datafile \"%s\"",
				vecOutputFileList[f].c_str());
		}

		// Copy metadata
		CopyNcFileAttributes(&ncinfile, &ncoutfile);

		// Copy time
		CopyNcVar(ncinfile, ncoutfile, "time");

		// Get time dimension
		NcDim * dimTimeOut = ncoutfile.get_dim("time");
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Unable to find \"time\" dimension");
		}

		// Copy latitude and/or longitude
		CopyNcVarIfExists(ncinfile, ncoutfile, "lat");
		CopyNcVarIfExists(ncinfile, ncoutfile, "latitude");
		CopyNcVarIfExists(ncinfile, ncoutfile, "LAT");
		CopyNcVarIfExists(ncinfile, ncoutfile, "lat_0");

		CopyNcVarIfExists(ncinfile, ncoutfile, "lon");
		CopyNcVarIfExists(ncinfile, ncoutfile, "longitude");
		CopyNcVarIfExists(ncinfile, ncoutfile, "LON");
		CopyNcVarIfExists(ncinfile, ncoutfile, "lon_0");

		// Copy preserve variables
		for (int p = 0; p < vecPreserveVariables.size(); p++) {
			CopyNcVar(ncinfile, ncoutfile, vecPreserveVariables[p]);

			if (vecPreserveVariables[p] == strMaskVariable) {
				_EXCEPTION1("Variable \"%s\" specified as --maskvar also "
					"appears in --preserve",
					strMaskVariable.c_str());
			}
		}

		// Copy masked variables
		std::vector<NcVar *> vecNcVarIn;
		std::vector<NcVar *> vecNcVarOut;

		std::vector<std::string> strGridDimNames;

		NcDim * dimNcVarOut0 = NULL;
		NcDim * dimNcVarOut1 = NULL;

		for (int v = 0; v < vecVarNames.size(); v++) {
			const std::string & strVariable = vecVarNames[v];

			if (vecVarNames[v] == strMaskVariable) {
				_EXCEPTION1("Variable \"%s\" specified as --maskvar also appears in --var",
					strMaskVariable.c_str());
			}

			// Copy the variable from input file to output file
			CopyNcVar(ncinfile, ncoutfile, strVariable, true, false);

			// Get pointers to the NcVars associated with these variables
			NcVar * varIn = ncinfile.get_var(strVariable.c_str());
			if (varIn == NULL) {
				_EXCEPTION();
			}
			vecNcVarIn.push_back(varIn);

			NcVar * varOut = ncoutfile.get_var(strVariable.c_str());
			if (varOut == NULL) {
				_EXCEPTION();
			}
			vecNcVarOut.push_back(varOut);

			// Count variable dimensions
			long nVarDims = varIn->num_dims();
			if (nVarDims < 1 + grid.m_nGridDim.size()) {
				_EXCEPTION2("Insufficient dimensions in variable \"%s\" in file \"%s\"",
					strVariable.c_str(), vecOutputFileList[f].c_str()); 
			}
			if (strcmp(varIn->get_dim(0)->name(), "time") != 0) {
				_EXCEPTION2("First dimension of variable \"%s\" in file \"%s\" must be \"time\"",
					strVariable.c_str(), vecOutputFileList[f].c_str());
			}

			// Get output dimensions
			if (grid.m_nGridDim.size() == 1) {
				if (varIn->get_dim(nVarDims-1)->size() != grid.m_nGridDim[0]) {
					_EXCEPTION3("Last dimension of variable \"%s\" in file \"%s\" must have size equal to grid size (%i)",
						strVariable.c_str(), vecOutputFileList[f].c_str(), grid.m_nGridDim[0]);
				}
				dimNcVarOut0 = ncoutfile.get_dim(varIn->get_dim(nVarDims-1)->name());
				if (dimNcVarOut0 == NULL) {
					_EXCEPTION1("Unable to match dimension \"%s\" in output file",
						varIn->get_dim(nVarDims-1)->name());
				}
			}
			if (grid.m_nGridDim.size() == 2) {
				if ((varIn->get_dim(nVarDims-2)->size() != grid.m_nGridDim[0]) &&
					(varIn->get_dim(nVarDims-1)->size() != grid.m_nGridDim[1])
				) {
					_EXCEPTION4("Last dimensions of variable \"%s\" in file \"%s\" must have size equal to grid size (%i,%i)",
						strVariable.c_str(), vecOutputFileList[f].c_str(), grid.m_nGridDim[0], grid.m_nGridDim[1]);
				}
				dimNcVarOut0 = ncoutfile.get_dim(varIn->get_dim(nVarDims-2)->name());
				if (dimNcVarOut0 == NULL) {
					_EXCEPTION1("Unable to match dimension \"%s\" in output file",
						varIn->get_dim(nVarDims-2)->name());
				}
				dimNcVarOut1 = ncoutfile.get_dim(varIn->get_dim(nVarDims-1)->name());
				if (dimNcVarOut1 == NULL) {
					_EXCEPTION1("Unable to match dimension \"%s\" in output file",
						varIn->get_dim(nVarDims-1)->name());
				}

			}
		}
		_ASSERT(vecVarNames.size() == vecNcVarIn.size());
		_ASSERT(vecVarNames.size() == vecNcVarOut.size());

		// No variables specified, try to guess mask variable dimensions
		if (dimNcVarOut0 == NULL) {
			if (grid.m_nGridDim.size() == 1) {
				dimNcVarOut0 = ncoutfile.get_dim("ncol");
				if (dimNcVarOut0 == NULL) {
					dimNcVarOut0 = ncoutfile.add_dim("ncol", grid.m_nGridDim[0]);
				} else if (dimNcVarOut0->size() != grid.m_nGridDim[0]) {
					_EXCEPTIONT("Unable to determine grid dimension in output file");
				}
			}
			if (grid.m_nGridDim.size() == 2) {
				dimNcVarOut0 = ncoutfile.get_dim("lat");
				dimNcVarOut1 = ncoutfile.get_dim("lon");
				if (dimNcVarOut0 == NULL) {
					dimNcVarOut0 = ncoutfile.add_dim("lat", grid.m_nGridDim[0]);
				} else if (dimNcVarOut0->size() != grid.m_nGridDim[0]) {
					_EXCEPTIONT("Unable to determine grid dimension in output file");
				}
				if (dimNcVarOut1 == NULL) {
					dimNcVarOut1 = ncoutfile.add_dim("lon", grid.m_nGridDim[1]);
				} else if (dimNcVarOut1->size() != grid.m_nGridDim[1]) {
					_EXCEPTIONT("Unable to determine grid dimension in output file");
				}
			}
		}
		_ASSERT(dimNcVarOut0 != NULL);
		if (grid.m_nGridDim.size() == 2) {
			_ASSERT(dimNcVarOut1 != NULL);
		}

		// Create masked variable
		NcVar * varMask = NULL;
		if (strMaskVariable != "") {
			if (grid.m_nGridDim.size() == 1) {
				varMask = ncoutfile.add_var(
					strMaskVariable.c_str(),
					ncDouble,
					dimTimeOut,
					dimNcVarOut0
				);

			} else {
				varMask = ncoutfile.add_var(
					strMaskVariable.c_str(),
					ncDouble,
					dimTimeOut,
					dimNcVarOut0,
					dimNcVarOut1
				);
			}
		}

		// Loop through all times
		for (int t = 0; t < vecTimes.size(); t++) {

			AnnounceStartBlock("Processing time \"%s\"",
				vecTimes[t].ToString().c_str());
	
			// Build mask
			TimeToPathNodeMap::const_iterator iter =
				mapTimeToPathNode.find(vecTimes[t]);

			dataMask.Zero();

			if (iter != mapTimeToPathNode.end()) {
				if (strFilterByDist != "") {
					BuildMask_ByDist(
						grid,
						cdhInput,
						pathvec,
						iter->second,
						strFilterByDist,
						dataMask);
				}
				if (strFilterByContour != "") {
					Variable & varOp = varreg.Get(vecClosedContourOp[0].m_varix);
					varOp.LoadGridData(varreg, vecInFiles, grid, t);
					const DataArray1D<float> & dataState = varOp.GetData();
					_ASSERT(dataState.GetRows() == grid.GetSize());

					BuildMask_ByContour<float>(
						grid,
						dataState,
						cdhInput,
						pathvec,
						iter->second,
						vecClosedContourOp[0],
						dataMask);
				}
				if (vecNearbyBlobsOp.size() != 0) {
					Variable & varOp = varreg.Get(vecNearbyBlobsOp[0].m_varix);
					varOp.LoadGridData(varreg, vecInFiles, grid, t);
					const DataArray1D<float> & dataState = varOp.GetData();
					_ASSERT(dataState.GetRows() == grid.GetSize());

					BuildMask_NearbyBlobs<float>(
						grid,
						dataState,
						cdhInput,
						pathvec,
						iter->second,
						vecNearbyBlobsOp[0],
						dataMask);
				}
			}
			if (fInvert) {
				for (int i = 0; i < dataMask.GetRows(); i++) {
					dataMask[i] = 1.0 - dataMask[i];
				}
			}

			// Write mask to file
			if (varMask != NULL) {
				if (grid.m_nGridDim.size() == 1) {
					varMask->set_cur(t,0);
					varMask->put(&(dataMask[0]), 1, grid.m_nGridDim[0]);
				} else {
					varMask->set_cur(t,0,0);
					varMask->put(&(dataMask[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
				}
			}

			// Loop through all variables
			for (int v = 0; v < vecVarIx.size(); v++) {

				const std::string & strVariable = vecVarNames[v];

				Announce("Processing variable \"%s\"", strVariable.c_str());

				// Load input and output variables
				NcVar * varIn = vecNcVarIn[v];
				NcVar * varOut = vecNcVarOut[v];

				// Load in data
				int nGridDims = grid.m_nGridDim.size();

				std::vector<NcDim *> vecDim;
				for (int d = 0; d < varIn->num_dims(); d++) {
					vecDim.push_back(varIn->get_dim(d));
				}
/*
				if (vecDim.size() < 1 + nGridDims) {
					_EXCEPTION2("Insufficient dimensions in variable \"%s\" in file \"%s\"",
						strVariable.c_str(), vecOutputFileList[f].c_str());
				}
				if (strcmp(vecDim[0]->name(), "time") != 0) {
					_EXCEPTION2("First dimension of variable \"%s\" in file \"%s\" must be \"time\"",
						strVariable.c_str(), vecOutputFileList[f].c_str());
				}
*/
				int nVarSize = 1;
				for (int d = 0; d < nGridDims; d++) {
					nVarSize *= vecDim[vecDim.size()-d-1]->size();
				}
				if (grid.GetSize() != nVarSize) {
					_EXCEPTION2("Variable \"%s\" in file \"%s\" dimensionality inconsistent with grid:"
						" Verify final variable dimensions match grid",
						strVariable.c_str(), vecOutputFileList[f].c_str());
				}

				// Number of auxiliary dimensions and array of sizes for each slice
				int nAuxDims = 1;
				for (int d = 1; d < vecDim.size() - nGridDims; d++) {
					nAuxDims *= vecDim[d]->size();
				}

				DataArray1D<long> vecDataSize(vecDim.size());
				for (int d = 0; d < vecDim.size(); d++) {
					if (d >= vecDim.size() - nGridDims) {
						vecDataSize[d] = vecDim[d]->size();
					} else {
						vecDataSize[d] = 1;
					}
				}

				// Check for _FillValue
				if (strFillValue == "att") {
					NcAtt * attFillValue = varIn->get_att("_FillValue");
					if (attFillValue == NULL) {
						fHasFillValue = false;
						Announce("WARNING: Variable \"%s\" in file \"%s\" does not have a _FillValue attribute", strVariable.c_str(), vecInputFileList[f].c_str());

					} else {
						fHasFillValue = true;
						dFillValue = attFillValue->as_float(0);
					}
				}

				// Add _FillValue to output variable
				if ((fHasFillValue) && (strFillValue != "nan")) {
					NcAtt * attFillValueOut = varOut->get_att("_FillValue");
					if (attFillValueOut != NULL) {
						attFillValueOut->remove();
					}
					if (varOut->type() == ncFloat) {
						varOut->add_att("_FillValue", static_cast<float>(dFillValue));
					} else if (varOut->type() == ncDouble) {
						varOut->add_att("_FillValue", static_cast<double>(dFillValue));
					} else {
						_EXCEPTION1("Invalid type for variable \"%s\": Expected \"float\" or \"double\"", strVariable.c_str());
					}
				}

				// Loop through all auxiliary dimensions (holding time fixed)
				DataArray1D<long> vecDataPos(vecDim.size());
				vecDataPos[0] = t;

				for (int i = 0; i < nAuxDims; i++) {

					// Load data
					int ixDim = i;
					for (int d = vecDim.size() - nGridDims - 1; d >= 1; d--) {
						vecDataPos[d] = ixDim % vecDim[d]->size();
						ixDim /= vecDim[d]->size();
					}
					if (ixDim != 0) {
						_EXCEPTIONT("Logic error");
					}

					varIn->set_cur(&(vecDataPos[0]));
					varIn->get(&(data[0]), &(vecDataSize[0]));

					// Apply mask
					if (!fHasFillValue) {
						for (int i = 0; i < data.GetRows(); i++) {
							data[i] *= dataMask[i];
						}

					} else {
						for (int i = 0; i < data.GetRows(); i++) {
							if (dataMask[i] == 0.0) {
								data[i] = dFillValue;
							}
						}
					}

					// Write data
					varOut->set_cur(&(vecDataPos[0]));
					varOut->put(&(data[0]), &(vecDataSize[0]));
				}
			}

			AnnounceEndBlock("Done");
		}

		AnnounceEndBlock("Done");
	}

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

///////////////////////////////////////////////////////////////////////////////


