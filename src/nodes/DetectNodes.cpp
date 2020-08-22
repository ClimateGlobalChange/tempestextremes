///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectNodes.cpp
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

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

#include "Variable.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "FilenameList.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "TimeObj.h"
#include "NodeOutputOp.h"
#include "ClosedContourOp.h"
#include "ThresholdOp.h"
#include "NetCDFUtilities.h"
#include "SimpleGridUtilities.h"

#include "kdtree.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <set>
#include <queue>

#ifndef TEMPEST_NOREGEX
#include <regex>
#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field has a closed contour about this point.
///	</summary>
template <typename real>
bool HasClosedContour(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
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
	const DataArray1D<real> & dataState,
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
		dMinAbsLatitude(0.0),
		dMaxLongitude(0.0),
		dMinLongitude(0.0),
		dMergeDist(0.0),
		pvecClosedContourOp(NULL),
		pvecNoClosedContourOp(NULL),
		pvecThresholdOp(NULL),
		pvecOutputOp(NULL),
		nTimeStride(1),
		strLatitudeName("lat"),
		strLongitudeName("lon"),
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

	// Minimum absolute value of latitude for detection
	double dMinAbsLatitude;

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
	std::vector<NodeOutputOp> * pvecOutputOp;

	// Time stride
	int nTimeStride;

	// Time filter
	std::string strTimeFilter;

	// Name of the latitude dimension
	std::string strLatitudeName;

	// Name of the longitude dimension
	std::string strLongitudeName;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

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
	if ((param.dMinAbsLatitude < 0.0) || (param.dMinAbsLatitude > 90.0)) {
		_EXCEPTIONT("--minabslat must in the range [0,90]");
	}

	// Dereference pointers to operators
	_ASSERT(param.pvecClosedContourOp != NULL);
	std::vector<ClosedContourOp> & vecClosedContourOp =
		*(param.pvecClosedContourOp);

	_ASSERT(param.pvecNoClosedContourOp != NULL);
	std::vector<ClosedContourOp> & vecNoClosedContourOp =
		*(param.pvecNoClosedContourOp);

	_ASSERT(param.pvecThresholdOp != NULL);
	std::vector<ThresholdOp> & vecThresholdOp =
		*(param.pvecThresholdOp);

	_ASSERT(param.pvecOutputOp != NULL);
	std::vector<NodeOutputOp> & vecOutputOp =
		*(param.pvecOutputOp);

#ifdef TEMPEST_NOREGEX
	if (param.strTimeFilter != "") {
		_EXCEPTIONT("Cannot use --timefilter with -DTEMPEST_NOREGEX compiler flag");
	}
#endif
#ifndef TEMPEST_NOREGEX
	// Parse --timefilter
	std::regex reTimeSubset;
	{
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

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Loading grid data from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		_ASSERT(vecFiles.size() > 0);

		AnnounceStartBlock("Generating RLL grid data");
		grid.GenerateLatitudeLongitude(
			vecFiles[0],
			param.strLatitudeName,
			param.strLongitudeName,
			param.fRegional,
			param.fDiagonalConnectivity);
		AnnounceEndBlock("Done");
	}

	// Get time information
	std::vector<Time> vecTimes;
	ReadCFTimeDataFromNcFile(
		vecFiles[0],
		vecFiles.GetFilename(0),
		vecTimes,
		false);

#ifndef TEMPEST_NOREGEX
	{
		std::vector<Time> vecOutputTimes;
		for (int t = 0; t < vecTimes.size(); t++) {
			std::string strTime = vecTimes[t].ToString();
			std::smatch match;
			if (std::regex_search(strTime, match, reTimeSubset)) {
				vecOutputTimes.push_back(vecTimes[t]);
			}
		}
		vecTimes = vecOutputTimes;
	}
#endif

/*
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

	DataArray1D<double> dTime(nTime);

	if (varTime->type() == ncDouble) {
		varTime->get(dTime, nTime);

	} else if (varTime->type() == ncFloat) {
		DataArray1D<float> dTimeFloat(nTime);

		varTime->get(dTimeFloat, nTime);
		for (int t = 0; t < nTime; t++) {
			dTime[t] = static_cast<double>(dTimeFloat[t]);
		}

	} else if (varTime->type() == ncInt) {
		DataArray1D<int> dTimeInt(nTime);

		varTime->get(dTimeInt, nTime);
		for (int t = 0; t < nTime; t++) {
			dTime[t] = static_cast<double>(dTimeInt[t]);
		}

	} else if (varTime->type() == ncInt64) {
		DataArray1D<ncint64> dTimeInt(nTime);

		varTime->get(dTimeInt, nTime);
		for (int t = 0; t < nTime; t++) {
			dTime[t] = static_cast<double>(dTimeInt[t]);
		}

	} else {
		_EXCEPTIONT("Variable \"time\" has an invalid type:\n"
			"Expected \"float\", \"double\", \"int\", or \"int64\"");
	}
*/
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
	for (int t = 0; t < vecTimes.size(); t += param.nTimeStride) {

		// Announce
		AnnounceStartBlock("Time %s", vecTimes[t].ToString().c_str());

		// Load the data for the search variable
		Variable & varSearchBy = varreg.Get(param.ixSearchBy);
		vecFiles.SetTime(vecTimes[t]);
		varSearchBy.LoadGridData(varreg, vecFiles, grid);

		const DataArray1D<float> & dataSearch = varSearchBy.GetData();
/*
		// Parse time information
		NcAtt * attTimeUnits = varTime->get_att("units");
		if (attTimeUnits == NULL) {
			_EXCEPTIONT("Variable \"time\" has no \"units\" attribute");
		}

		std::string strTimeUnits = attTimeUnits->as_string(0);

		Time::CalendarType eCalendarType = Time::CalendarStandard;
		NcAtt * attTimeCalendar = varTime->get_att("calendar");
		if (attTimeCalendar != NULL) {
			eCalendarType = Time::CalendarTypeFromString(attTimeCalendar->as_string(0));
			if (eCalendarType == Time::CalendarUnknown) {
				_EXCEPTIONT("Unknown calendar type associated with variable \"time\"");
			}
		}

		Time time(eCalendarType);
		time.FromCFCompliantUnitsOffsetDouble(strTimeUnits, dTime[t]);
*/
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

		DataArray1D<int> vecRejectedClosedContour(vecClosedContourOp.size());
		DataArray1D<int> vecRejectedNoClosedContour(vecNoClosedContourOp.size());
		DataArray1D<int> vecRejectedThreshold(vecThresholdOp.size());

		// Eliminate based on interval
		if ((param.dMinLatitude != param.dMaxLatitude) ||
		    (param.dMinLongitude != param.dMaxLongitude) ||
			(param.dMinAbsLatitude != 0.0)
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
				if (param.dMinAbsLatitude != 0.0) {
					if (fabs(dLat) < param.dMinAbsLatitude) {
						nRejectedLocation++;
						continue;
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
			if (kdMerge == NULL) {
				_EXCEPTIONT("kd_create(3) failed");
			}

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
			vecFiles.SetTime(vecTimes[t]);
			var.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataState = var.GetData();

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
			vecFiles.SetTime(vecTimes[t]);
			var.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataState = var.GetData();

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
			vecFiles.SetTime(vecTimes[t]);
			var.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataState = var.GetData();

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
			Announce("Rejected (thresh. %s): %i",
					varreg.GetVariableString(vecThresholdOp[tc].m_varix).c_str(),
					vecRejectedThreshold[tc]);
		}

		for (int ccc = 0; ccc < vecRejectedClosedContour.GetRows(); ccc++) {
			Announce("Rejected (contour %s): %i",
					varreg.GetVariableString(vecClosedContourOp[ccc].m_varix).c_str(),
					vecRejectedClosedContour[ccc]);
		}

		for (int ccc = 0; ccc < vecRejectedNoClosedContour.GetRows(); ccc++) {
			Announce("Rejected (nocontour %s): %i",
					varreg.GetVariableString(vecNoClosedContourOp[ccc].m_varix).c_str(),
					vecRejectedNoClosedContour[ccc]);
		}

		// Write results to file
		{
			// Write time information
			fprintf(fpOutput, "%i\t%i\t%i\t%i\t%i\n",
				vecTimes[t].GetYear(),
				vecTimes[t].GetMonth(),
				vecTimes[t].GetDay(),
				static_cast<int>(setCandidates.size()),
				vecTimes[t].GetSecond() / 3600);
/*
			if (param.fOutputInfileInfo) {
				fprintf(fpOutput, "\t\"%s\"\t%i\n", strInputFiles.c_str(), t);
			} else {
				fprintf(fpOutput, "\n");
			}
*/
			// Write candidate information
			int iCandidateIx = 0;

			// Apply output operators
			std::vector< std::vector<std::string> > vecOutputValue;
			vecOutputValue.resize(setCandidates.size());
			for (int i = 0; i < setCandidates.size(); i++) {
				vecOutputValue[i].resize(vecOutputOp.size());
			}

			//DataArray2D<float> dOutput(setCandidates.size(), vecOutputOp.size());
			for (int outc = 0; outc < vecOutputOp.size(); outc++) {

				// Loop through all pressure minima
				std::set<int>::const_iterator iterCandidate
					= setCandidates.begin();

				iCandidateIx = 0;
				for (; iterCandidate != setCandidates.end(); iterCandidate++) {
					ApplyNodeOutputOp<float>(
						vecOutputOp[outc],
						grid,
						varreg,
						vecFiles,
						vecTimes[t],
						*iterCandidate,
						vecOutputValue[iCandidateIx][outc]);

					iCandidateIx++;
				}
			}

			// Output all candidates
			iCandidateIx = 0;

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
					fprintf(fpOutput, "\t%s", vecOutputValue[iCandidateIx][outc].c_str());
				}

				fprintf(fpOutput, "\n");

				iCandidateIx++;
			}
		}

		AnnounceEndBlock("Done");
	}

	fclose(fpOutput);

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

	// Input data file
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

	// Closed contour commands
	std::string strClosedContourCmd;

	// Closed contour commands
	std::string strNoClosedContourCmd;

	// Threshold commands
	std::string strThresholdCmd;

	// Output commands
	std::string strOutputCmd;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(dcuparam.fDiagonalConnectivity, "diag_connect");
		CommandLineString(strOutput, "out", "");
		CommandLineString(strOutputFileList, "out_file_list", "");
		CommandLineStringD(strSearchByMin, "searchbymin", "", "(default PSL)");
		CommandLineString(strSearchByMax, "searchbymax", "");
		CommandLineDoubleD(dcuparam.dMinLongitude, "minlon", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMaxLongitude, "maxlon", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMinAbsLatitude, "minabslat", 0.0, "(degrees)");
		CommandLineDoubleD(dcuparam.dMergeDist, "mergedist", 0.0, "(degrees)");
		CommandLineStringD(strClosedContourCmd, "closedcontourcmd", "", "[var,delta,dist,minmaxdist;...]");
		CommandLineStringD(strNoClosedContourCmd, "noclosedcontourcmd", "", "[var,delta,dist,minmaxdist;...]");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		CommandLineStringD(strOutputCmd, "outputcmd", "", "[var,op,dist;...]");
		CommandLineInt(dcuparam.nTimeStride, "timestride", 1);
		CommandLineString(dcuparam.strTimeFilter, "timefilter", "");
		CommandLineString(dcuparam.strLatitudeName, "latname", "lat");
		CommandLineString(dcuparam.strLongitudeName, "lonname", "lon");
		CommandLineBool(dcuparam.fRegional, "regional");
		CommandLineBool(dcuparam.fOutputHeader, "out_header");
		CommandLineInt(dcuparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Note timestride is deprecated
	if (dcuparam.nTimeStride != 1) {
		Announce("WARNING: --timestride is deprecated.  Consider using --timefilter instead.");
		if (dcuparam.strTimeFilter != "") {
			_EXCEPTIONT("Only one of --timestride and --timefilter can be used.");
		}
	}

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
	FilenameList vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		vecInputFiles.FromFile(strInputFileList, true);
	}

	// Load output file list
	FilenameList vecOutputFiles;

	if (strOutput != "") {
		vecOutputFiles.push_back(strOutput);

	} else {
		vecOutputFiles.FromFile(strOutputFileList);

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
		if (strSearchByMin != "") {
			dcuparam.ixSearchBy = varreg.FindOrRegister(strSearchByMin);
			dcuparam.fSearchByMinima = true;
		}
		if (strSearchByMax != "") {
			dcuparam.ixSearchBy = varreg.FindOrRegister(strSearchByMax);
			dcuparam.fSearchByMinima = false;
		}
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
	std::vector<NodeOutputOp> vecOutputOp;
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
	dcuparam.dMinAbsLatitude *= M_PI / 180.0;

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
			static_cast<int>(-dcuparam.dMaxLongitude / 360.0);
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

		if (vecOutputFiles.size() == 0) {
			Announce("Output will be written to outXXXXXX.dat");
		} else if (vecOutputFiles.size() == 1) {
			Announce("Output will be written to %sXXXXXX.dat",
				vecOutputFiles[0].c_str());
		} else {
			Announce("Output will be written following --out_file_list");
		}
		Announce("Logs will be written to logXXXXXX.txt");
	}

	// Loop over all files to be processed
	_ASSERT(vecInputFiles.size() == vecOutputFiles.size());

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

			if (vecOutputFiles.size() == 0) {
				strOutputFile = "out.dat";
			} else {
				strOutputFile = strOutput;
			}

		} else {
			char szFileIndex[32];
			sprintf(szFileIndex, "%06i", f);

			if (vecOutputFiles.size() == 0) {
				strOutputFile =
					"out" + std::string(szFileIndex) + ".dat";
			} else if (vecOutputFiles.size() == 1) {
				strOutputFile =
					vecOutputFiles[0] + std::string(szFileIndex) + ".dat";
			} else {
				strOutputFile =
					vecOutputFiles[f];
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

///////////////////////////////////////////////////////////////////////////////

