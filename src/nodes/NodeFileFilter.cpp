///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileFilter.cpp
///	\author  Paul Ullrich
///	\version July 2, 2020
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
#include "FilenameList.h"
#include "TimeMatch.h"

#include "netcdfcpp.h"

#include <fstream>
#include <queue>
#include <set>
#include <cmath>

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

		snprintf(szBuffer, 128, "%f", m_dDistance);
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
			snprintf(szBuffer, 128, "%e", m_dValue);
		} else {
			snprintf(szBuffer, 128, "%f", m_dValue);
		}
		strDescription += szBuffer;

		snprintf(szBuffer, 128, "%f", m_dMaxDist);
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
	const NodeFile & nodefile,
	const std::vector<double> & vecLonRad,
	const std::vector<double> & vecLatRad,
	const std::string & strDist,
	DataArray1D<double> & dataMask
) {
	_ASSERT(grid.GetSize() == dataMask.GetRows());

	// Get path ids from interpolation
	const std::vector<size_t> & vecPathId = nodefile.GetInterpolatedPathIds();

	_ASSERT(vecPathId.size() == vecLonRad.size());
	_ASSERT(vecPathId.size() == vecLatRad.size());

	// Get filter width (either fixed value or column data header)
	if (STLStringHelper::IsFloat(strDist)) {
		double dFilterWidthDeg = std::stod(strDist);
		if (dFilterWidthDeg == 0.0) {
			return;
		}

		std::vector<size_t> vecNodeIxs;
		for (int i = 0; i < vecPathId.size(); i++) {
			grid.NearestNodes(vecLonRad[i], vecLatRad[i], dFilterWidthDeg, vecNodeIxs);
			for (int n = 0; n < vecNodeIxs.size(); n++) {
				dataMask[vecNodeIxs[n]] = 1.0;
			}
		}

	} else {
		std::vector<double> vecFilterWidthDeg;
		nodefile.InterpolatedColumnDouble(
			strDist, vecFilterWidthDeg);

		std::vector<size_t> vecNodeIxs;
		for (int i = 0; i < vecPathId.size(); i++) {
			grid.NearestNodes(vecLonRad[i], vecLatRad[i], vecFilterWidthDeg[i], vecNodeIxs);
			for (int n = 0; n < vecNodeIxs.size(); n++) {
				dataMask[vecNodeIxs[n]] = 1.0;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void BuildMask_ByContour(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	std::vector<double> & vecLonRad,
	std::vector<double> & vecLatRad,
	const ClosedContourOp & op,
	DataArray1D<double> & dataMask
) {
	// Get the variable
	double dDeltaAmt = op.m_dDeltaAmount;
	double dDeltaDistDeg = op.m_dDistance;
	double dMinMaxDistDeg = op.m_dMinMaxDist;

	double dDeltaDistRad = DegToRad(dDeltaDistDeg);

	_ASSERT(vecLonRad.size() == vecLatRad.size());
	_ASSERT(dDeltaAmt != 0.0);
	_ASSERT(dDeltaDistDeg > 0.0);

	// Loop through all PathNodes
	for (int j = 0; j < vecLonRad.size(); j++) {

		double dLon0Rad = vecLonRad[j];
		double dLat0Rad = vecLatRad[j];

		int ix0 = static_cast<int>(
			grid.NearestNode(dLon0Rad, dLat0Rad));
		_ASSERT((ix0 >= 0) && (ix0 < grid.GetSize()));

		// Find min/max near point
		int ixOrigin;

		if (dMinMaxDistDeg == 0.0) {
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
				dMinMaxDistDeg,
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

		// Build up nodes
		while (queueToVisit.size() != 0) {
			int ix = queueToVisit.front();
			queueToVisit.pop();

			if (setNodesVisited.find(ix) != setNodesVisited.end()) {
				continue;
			}

			setNodesVisited.insert(ix);

			// Great circle distance to this element
			double dLatThisRad = grid.m_dLat[ix];
			double dLonThisRad = grid.m_dLon[ix];

			double dRRad =
				GreatCircleDistance_Rad(
					dLon0Rad, dLat0Rad,
					dLonThisRad, dLatThisRad);

			// Check great circle distance
			if (dRRad > dDeltaDistRad) {
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
	std::vector<double> & vecLonRad,
	std::vector<double> & vecLatRad,
	const NearbyBlobsOp & nearbyblobsop,
	DataArray1D<double> & dataMask
) {
	// Get the variable
	const double dDist = nearbyblobsop.m_dDistance;
	const double dMaxDist = nearbyblobsop.m_dMaxDist;

	_ASSERT(vecLonRad.size() == vecLatRad.size());
	_ASSERT(dataState.GetRows() == grid.GetSize());
	_ASSERT(dDist >= 0.0);
	_ASSERT(dDist <= 180.0);
	_ASSERT(dMaxDist >= dDist);
	_ASSERT(dMaxDist <= 180.0);

	// Loop through all PathNodes
	for (int j = 0; j < vecLonRad.size(); j++) {

		double dLon0 = vecLonRad[j];
		double dLat0 = vecLatRad[j];

		int ix0 = static_cast<int>(
			grid.NearestNode(dLon0, dLat0));
		_ASSERT((ix0 >= 0) && (ix0 < grid.GetSize()));

		// Queue of nodes that remain to be visited
		std::queue<int> queueNodes;
		queueNodes.push(ix0);

		// Set of nodes that have already beenvisited
		std::set<int> setNodesVisited;

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
				GreatCircleDistance_Rad(
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

void NodeFileFilter(
	const std::string & strInputFile,
	const std::string & strOutputFile,
	VariableRegistry & varreg,
	SimpleGrid & grid,
	NodeFile & nodefile,
	const Time::CalendarType & caltype,
	ColumnDataHeader & cdhInput,
	const std::vector<std::string> & vecVarNames,
	const std::vector<std::string> & vecVarOutNames,
	const std::vector<VariableIndex> & vecVarIx,
	const std::string & strMaskVariable,
	const std::string & strFilterByDist,
	const std::vector<ClosedContourOp> & vecClosedContourOp,
	const std::vector<NearbyBlobsOp> & vecNearbyBlobsOp,
	bool fInvert,
	const std::vector<std::string> & vecPreserveVariables,
	const std::string & strFillValueIn,
	const std::string & strTimeFilterIn,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName
) {
	_ASSERT(vecVarNames.size() == vecVarOutNames.size());

	// Begin processing
	AnnounceStartBlock("Processing input (%s)", strInputFile.c_str());

	// Parse fillvalue
	std::string strFillValue = strFillValueIn;
	bool fHasFillValue = false;
	float dFillValue = 0.0;
	if (strFillValue != "") {
		fHasFillValue = true;
		STLStringHelper::ToLower(strFillValue);
		if (strFillValue == "att") {
		} else if (strFillValue == "nan") {
			dFillValue = nanf(NULL);
		} else if (STLStringHelper::IsFloat(strFillValue)) {
			dFillValue = std::stod(strFillValue);
		} else {
			_EXCEPTIONT("Invalid value specified for --fillvalue");
		}
	}

	// Parse --timefilter
	std::string strTimeFilter = strTimeFilterIn;
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

	// Open input file
	NcFileVector vecInFiles;
	vecInFiles.ParseFromString(strInputFile);
	_ASSERT(vecInFiles.size() > 0);

	NcFile & ncinfile = *(vecInFiles[0]);
	//if (!ncinfile.is_valid()) {
	//	_EXCEPTIONT("Unable to open input datafile");
	//}

	// Get time variable and verify calendar is compatible
	NcVar * varTime = NcGetTimeVariable(ncinfile);
	if (varTime == NULL) {
		_EXCEPTION1("Missing \"time\" variable in datafile \"%s\"",
			strInputFile.c_str());
	}

	NcAtt * attCalendar = varTime->get_att("calendar");
	if (attCalendar == NULL) {
		_EXCEPTION1("Missing \"calendar\" metadata in datafile \"%s\"",
			strInputFile.c_str());
	}

	Time::CalendarType caltypeThis =
		Time::CalendarTypeFromString(attCalendar->as_string(0));
	if (caltype != caltypeThis) {
		_EXCEPTION1("\"calendar\" mismatch in datafile \"%s\" ",
			strInputFile.c_str());
	}

	// Load in vector of Times
	NcTimeDimension vecTimes;
	ReadCFTimeDataFromNcFile(&ncinfile, strInputFile, vecTimes, true);

	// Open output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputFile.c_str());
	}

	// Copy metadata
	CopyNcFileAttributes(&ncinfile, &ncoutfile);

	// Construct time vector and output
	NcTimeDimension vecTimesOut;

#ifndef TEMPEST_NOREGEX
	if (strTimeFilter != "") {
		vecTimesOut.m_nctype = vecTimes.m_nctype;
		vecTimesOut.m_units = vecTimes.m_units;
		vecTimesOut.m_dimtype = vecTimes.m_dimtype;
		for (long t = 0; t < vecTimes.size(); t++) {
			std::string strTime = vecTimes[t].ToString();
			std::smatch match;
			if (!std::regex_search(strTime, match, reTimeSubset)) {
				continue;
			}
			vecTimesOut.push_back(vecTimes[t]);
		}
	} else {
		vecTimesOut = vecTimes;
	}
#endif
#ifdef TEMPEST_NOREGEX
	vecTimesOut = vecTimes;
#endif

	WriteCFTimeDataToNcFile(
		&ncoutfile,
		strOutputFile,
		vecTimesOut);

	NcDim * dimTimeOut = NcGetTimeDimension(ncoutfile);
	if (dimTimeOut == NULL) {
		_EXCEPTION1("Error writing \"time\" dimension to file \"%s\"",
			strOutputFile.c_str());
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

	// Create data array
	DataArray1D<double> data(grid.GetSize());

	// Create mask
	DataArray1D<double> dataMask(grid.GetSize());

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

		// Get pointers to the NcVars associated with these variables
		NcVar * varIn = ncinfile.get_var(strVariable.c_str());
		if (varIn == NULL) {
			_EXCEPTION();
		}
		vecNcVarIn.push_back(varIn);

		// Copy the variable from input file to output file
		//CopyNcVar(ncinfile, ncoutfile, strVariable, true, false);

		{
			std::vector<NcDim *> vecDimOut;
			for (long d = 0; d < varIn->num_dims(); d++) {
				if (NcIsTimeDimension(varIn->get_dim(d))) {
					vecDimOut.push_back(dimTimeOut);
					continue;
				}

				NcDim * dimOut = ncoutfile.get_dim(varIn->get_dim(d)->name());
				if (dimOut == NULL) {
					dimOut = ncoutfile.add_dim(varIn->get_dim(d)->name(), varIn->get_dim(d)->size());
					if (dimOut == NULL) {
						_EXCEPTION2("Unable to create dimension \"%s\" in file \"%s\"",
							varIn->get_dim(d)->name(), strOutputFile.c_str());
					}
				}
				vecDimOut.push_back(dimOut);
			}

			NcVar * varOut = ncoutfile.add_var(
				vecVarOutNames[v].c_str(),
				varIn->type(),
				vecDimOut.size(),
				(const NcDim**)&(vecDimOut[0]));

			if (varOut == NULL) {
				_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
					varIn->name(), strOutputFile.c_str());
			}

			CopyNcVarAttributes(varIn, varOut);

			vecNcVarOut.push_back(varOut);
		}

		// Count variable dimensions
		long nVarDims = varIn->num_dims();
		if (nVarDims < 1 + grid.m_nGridDim.size()) {
			_EXCEPTION2("Insufficient dimensions in variable \"%s\" in file \"%s\"",
				strVariable.c_str(), strInputFile.c_str()); 
		}
		if (!NcIsTimeDimension(varIn->get_dim(0))) {
			_EXCEPTION2("First dimension of variable \"%s\" in file \"%s\" must be \"time\"",
				strVariable.c_str(), strInputFile.c_str());
		}

		// Get output dimensions
		if (grid.m_nGridDim.size() == 1) {
			if (varIn->get_dim(nVarDims-1)->size() != grid.m_nGridDim[0]) {
				_EXCEPTION3("Last dimension of variable \"%s\" in file \"%s\" must have size equal to grid size (%i)",
					strVariable.c_str(), strInputFile.c_str(), grid.m_nGridDim[0]);
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
					strVariable.c_str(), strOutputFile.c_str(), grid.m_nGridDim[0], grid.m_nGridDim[1]);
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
			dimNcVarOut0 = ncoutfile.get_dim(strLatitudeName.c_str());
			dimNcVarOut1 = ncoutfile.get_dim(strLongitudeName.c_str());
			if (dimNcVarOut0 == NULL) {
				dimNcVarOut0 = ncoutfile.add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]);
			} else if (dimNcVarOut0->size() != grid.m_nGridDim[0]) {
				_EXCEPTIONT("Unable to determine grid dimension in output file");
			}
			if (dimNcVarOut1 == NULL) {
				dimNcVarOut1 = ncoutfile.add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]);
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

		if (varMask == NULL) {
			_EXCEPTION1("Unable to add variable \"%s\" to output file", strMaskVariable.c_str());
		}
	}

	// Loop through all times
	for (int t = 0, to = -1; t < vecTimes.size(); t++) {

		AnnounceStartBlock("Processing time \"%s\"",
			vecTimes[t].ToString().c_str());

#ifndef TEMPEST_NOREGEX
		if (strTimeFilter != "") {
			std::string strTime = vecTimes[t].ToString();
			std::smatch match;
			if (!std::regex_search(strTime, match, reTimeSubset)) {
				AnnounceEndBlock("(skipping)");
				continue;
			}
		}
#endif
		to++;

		// Find all paths at this time
		nodefile.Interpolate(vecTimes[t]);

		// Get coordinates of nodes at this time
		std::vector<double> vecLonRad;
		std::vector<double> vecLatRad;
		nodefile.InterpolatedNodeCoordinatesRad("lon", "lat", vecLonRad, vecLatRad);

		// Build mask
		dataMask.Zero();

		if (strFilterByDist != "") {
			AnnounceStartBlock("Building mask (bydist)");
			BuildMask_ByDist(
				grid,
				nodefile,
				vecLonRad,
				vecLatRad,
				strFilterByDist,
				dataMask);
			AnnounceEndBlock("Done");
		}
		if (vecClosedContourOp.size() != 0) {
			AnnounceStartBlock("Building mask (bycontour)");
			Variable & varOp = varreg.Get(vecClosedContourOp[0].m_varix);
			vecInFiles.SetConstantTimeIx(t);
			varOp.LoadGridData(varreg, vecInFiles, grid);
			const DataArray1D<float> & dataState = varOp.GetData();
			_ASSERT(dataState.GetRows() == grid.GetSize());

			BuildMask_ByContour<float>(
				grid,
				dataState,
				vecLonRad,
				vecLatRad,
				vecClosedContourOp[0],
				dataMask);

			AnnounceEndBlock("Done");
		}
		if (vecNearbyBlobsOp.size() != 0) {
			AnnounceStartBlock("Building mask (nearbyblobs)");
			Variable & varOp = varreg.Get(vecNearbyBlobsOp[0].m_varix);
			vecInFiles.SetConstantTimeIx(t);
			varOp.LoadGridData(varreg, vecInFiles, grid);
			const DataArray1D<float> & dataState = varOp.GetData();
			_ASSERT(dataState.GetRows() == grid.GetSize());

			BuildMask_NearbyBlobs<float>(
				grid,
				dataState,
				vecLonRad,
				vecLatRad,
				vecNearbyBlobsOp[0],
				dataMask);

			AnnounceEndBlock("Done");
		}

		if (fInvert) {
			for (int i = 0; i < dataMask.GetRows(); i++) {
				dataMask[i] = 1.0 - dataMask[i];
			}
		}

		// Write mask to file
		if (varMask != NULL) {
			AnnounceStartBlock("Writing mask variable to file");
			if (grid.m_nGridDim.size() == 1) {
				varMask->set_cur(to,0);
				varMask->put(&(dataMask[0]), 1, grid.m_nGridDim[0]);
			} else {
				varMask->set_cur(to,0,0);
				varMask->put(&(dataMask[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
			}
			AnnounceEndBlock("Done");
		}

		// Loop through all variables
		for (int v = 0; v < vecVarNames.size(); v++) {

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

			int nVarSize = 1;
			for (int d = 0; d < nGridDims; d++) {
				nVarSize *= vecDim[vecDim.size()-d-1]->size();
			}
			if (grid.GetSize() != nVarSize) {
				_EXCEPTION2("Variable \"%s\" in file \"%s\" dimensionality inconsistent with grid:"
					" Verify final variable dimensions match grid",
					strVariable.c_str(), strInputFile.c_str());
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
					Announce("WARNING: Variable \"%s\" in file \"%s\" does not have a _FillValue attribute",
						strVariable.c_str(), strInputFile.c_str());

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
				} else if (varOut->type() == ncShort) {
					varOut->add_att("_FillValue", static_cast<short>(dFillValue));
				} else if (varOut->type() == ncInt) {
					varOut->add_att("_FillValue", static_cast<int>(dFillValue));
				} else {
					_EXCEPTION1("Invalid type for variable \"%s\"", strVariable.c_str());
				}
			}

			// Loop through all auxiliary dimensions (holding time fixed)
			DataArray1D<long> vecDataPos(vecDim.size());

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

				vecDataPos[0] = t;
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
				vecDataPos[0] = to;
				varOut->set_cur(&(vecDataPos[0]));
				varOut->put(&(data[0]), &(vecDataSize[0]));
			}
		}

		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");
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

	// Grid file
	std::string strGridFile;

	// Data is regional
	bool fRegional;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Output data file
	std::string strOutputData;

	// Output list of data files
	std::string strOutputDataList;

	// List of variables to filter
	std::string strVariables;

	// List of output names for the variables to filter
	std::string strVariableOutputs;

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

	// Time filter
	std::string strTimeFilter;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Log dir
	std::string strLogDir;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_nodefile", "");
		CommandLineStringD(strInputNodeFileType, "in_nodefile_type", "SN", "[DN|SN]");
		CommandLineString(strInputFormat, "in_fmt", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strGridFile, "in_gridfile", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputDataList, "out_data_list", "");

		CommandLineString(strVariables, "var", "");
		CommandLineString(strVariableOutputs, "varout", "");
		CommandLineString(strMaskVariable, "maskvar", "");
		CommandLineString(strPreserve, "preserve", "");
		CommandLineStringD(strFillValue, "fillvalue", "", "[<value>|nan|att]");

		CommandLineStringD(strFilterByDist, "bydist", "", "[dist]");
		CommandLineStringD(strFilterByContour, "bycontour", "", "[var,delta,dist,minmaxdist]");
		CommandLineStringD(strNearbyBlobs, "nearbyblobs", "", "[var,dist,op,value[,maxdist]]");
		CommandLineBool(fInvert, "invert");

		CommandLineString(strTimeFilter, "timefilter", "");
		CommandLineString(strLongitudeName, "lonname", "[auto]");
		CommandLineString(strLatitudeName, "latname", "[auto]");

		CommandLineString(strLogDir, "logdir", ".");

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
	//if ((strFilterByDist.length() != 0) && (strFilterByContour.length() != 0)) {
	//	_EXCEPTIONT("Only one filter command (--bydist or --bycontour) may be specified");
	//}
	if ((strVariables.length() == 0) && (strMaskVariable.length() == 0)) {
		_EXCEPTIONT("One of (--var) or (--maskvar) must be specified");
	}

	// Input file type
	NodeFile::PathType iftype;
	if (strInputNodeFileType == "DN") {
		iftype = NodeFile::PathTypeDN;
	} else if (strInputNodeFileType == "SN") {
		iftype = NodeFile::PathTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_nodefile_type, expected \"SN\" or \"DN\"");
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

	// Parse --varout argument
	std::vector<std::string> vecVarOutNames;
	if (strVariableOutputs == "") {
		vecVarOutNames = vecVarNames;
	} else {
		STLStringHelper::ParseVariableList(strVariableOutputs, vecVarOutNames);

		if (vecVarOutNames.size() != vecVarNames.size()) {
			_EXCEPTIONT("--var and --varout must have the same length");
		}
	}

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
		AnnounceStartBlock("Parsing --bycontour");
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

	// Store input data
	FilenameList vecInputFileList;

	if (strInputData.length() != 0) {
		vecInputFileList.push_back(strInputData);

	} else {
		AnnounceStartBlock("Building input data list");
		vecInputFileList.FromFile(strInputDataList, false);
		AnnounceEndBlock("Done");
	}

	// Store output data
	FilenameList vecOutputFileList;

	if (strOutputData.length() != 0) {
		vecOutputFileList.push_back(strOutputData);

	} else {
		AnnounceStartBlock("Building output data list");
		vecOutputFileList.FromFile(strOutputDataList, false);
		AnnounceEndBlock("Done");
	}

	// Verify both input data and output data have the same length
	if (vecInputFileList.size() != vecOutputFileList.size()) {
		_EXCEPTION2("Mismatch in --in_data_list (%lu files) and --out_data_list length (%lu files)",
			vecInputFileList.size(), vecOutputFileList.size());
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Check for grid file
	if (strGridFile != "") {
		AnnounceStartBlock("Generating grid information from grid file");

		NcFile ncFile(strGridFile.c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"", strGridFile.c_str());
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

	// Check for connectivity file
	} else if (strConnectivity != "") {
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

	// Build the KD tree for the grid
	AnnounceStartBlock("Generating KD tree on grid");
	grid.BuildKDTree();
	AnnounceEndBlock("Done");

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

		NcVar * varTime = NcGetTimeVariable(ncfile);
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
	AnnounceStartBlock("Parsing nodefile");
	nodefile.Read(
		strInputNodeFile,
		iftype,
		cdhInput,
		grid,
		caltype);

	nodefile.GenerateTimeToPathNodeMap();
	AnnounceEndBlock("Done");

#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

	// Set up logging
	if ((vecInputFileList.size() > 1) && (nMPISize > 1)) {
		Announce("Logs will be written to %s/logXXXXXX.txt", strLogDir.c_str());
	}
#endif

	// Loop over all files
	_ASSERT(vecInputFileList.size() == vecOutputFileList.size());

	for (size_t f = 0; f < vecInputFileList.size(); f++) {
#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}

		FILE * fpLog = NULL;
		if ((vecInputFileList.size() > 1) && (nMPISize > 1)) {
			char szFileIndex[32];
			snprintf(szFileIndex, 32, "%06lu", f);

			std::string strLogFile = strLogDir + std::string("/log") + szFileIndex + ".txt";

			fpLog = fopen(strLogFile.c_str(), "w");
			if (fpLog == NULL) {
				_EXCEPTION1("Unable to open log file \"%s\"", strLogFile.c_str());
			}
			AnnounceSetOutputBuffer(fpLog);
			AnnounceOutputOnAllRanks();
		}
#endif
		NodeFileFilter(
			vecInputFileList[f],
			vecOutputFileList[f],
			varreg,
			grid,
			nodefile,
			caltype,
			cdhInput,
			vecVarNames,
			vecVarOutNames,
			vecVarIx,
			strMaskVariable,
			strFilterByDist,
			vecClosedContourOp,
			vecNearbyBlobsOp,
			fInvert,
			vecPreserveVariables,
			strFillValue,
			strTimeFilter,
			strLatitudeName,
			strLongitudeName);

#if defined(TEMPEST_MPIOMP)
		AnnounceSetOutputBuffer(stdout);
		if (fpLog != NULL) {
			fclose(fpLog);
		}
		AnnounceOnlyOutputOnRankZero();
#endif

	}

#if defined(TEMPEST_MPIOMP)
	MPI_Barrier(MPI_COMM_WORLD);
#endif

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


