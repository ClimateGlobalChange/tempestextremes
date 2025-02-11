///////////////////////////////////////////////////////////////////////////////
///
///	\file    BlobStats.cpp
///	\author  Paul Ullrich
///	\version May 14th, 2021
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "BlobUtilities.h"
#include "Constants.h"
#include "CoordTransforms.h"
#include "Variable.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "FilenameList.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "TimeObj.h"
#include "TimeMatch.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Storage structure for data associated with Blobs.
///	</summary>
class BlobQuantities {

public:
	///	<summary>
	///		Output quantities.
	///	</summary>
	enum OutputQuantity {
		MinLat,
		MaxLat,
		MinLon,
		MaxLon,
		MeanLat,
		MeanLon,
		CentroidLat,
		CentroidLon,
		Area,
		WtCentroidLat,
		WtCentroidLon,
		WtArea,
		WtPCALength,
		WtPCAWidth
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	BlobQuantities(
		size_t sSumVarCount
	) :
		dAreaX(0.0),
		dAreaY(0.0),
		dAreaZ(0.0),
		dArea(0.0),
		dWtAreaX(0.0),
		dWtAreaY(0.0),
		dWtAreaZ(0.0),
		dWtArea(0.0),
		dWtCentLonDeg(0.0),
		dWtCentLatDeg(0.0),
		dSumVars(sSumVarCount)
	{
		dCorrComp[0] = 0.0;
		dCorrComp[1] = 0.0;
		dCorrComp[2] = 0.0;
	}

public:
	///	<summary>
	///		Latitude-longitude bounding box for blob.
	///	</summary>
	LatLonBox<double> box;

	///	<summary>
	///		Area-weighted 3D X coordinate.
	///	</summary>
	double dAreaX;

	///	<summary>
	///		Area-weighted 3D Y coordinate.
	///	</summary>
	double dAreaY;

	///	<summary>
	///		Area-weighted 3D Z coordinate.
	///	</summary>
	double dAreaZ;

	///	<summary>
	///		Total area of blob.
	///	</summary>
	double dArea;

	///	<summary>
	///		Area and variable-weighted 3D X coordinate.
	///	</summary>
	double dWtAreaX;

	///	<summary>
	///		Area and variable-weighted 3D Y coordinate.
	///	</summary>
	double dWtAreaY;

	///	<summary>
	///		Area and variable-weighted 3D Z coordinate.
	///	</summary>
	double dWtAreaZ;

	///	<summary>
	///		Variable-weighted area of blob.
	///	</summary>
	double dWtArea;

	///	<summary>
	///		Blob centroid longitude.
	///	</summary>
	double dWtCentLonDeg;

	///	<summary>
	///		Blob centroid latitude.
	///	</summary>
	double dWtCentLatDeg;

	///	<summary>
	///		Correlation matrix components.
	///	</summary>
	double dCorrComp[3];

	///	<summary>
	///		Array of sums associated with blob.
	///	</summary>
	std::vector<double> dSumVars;
};

///	<summary>
///		A map between time and quantities associated with the blob.
///	</summary>
typedef std::map<int, BlobQuantities> TimedBlobQuantitiesMap;

///	<summary>
///		A map between blob index and TimedBlobQuantitiesMap;
///	</summary>
typedef std::map<int, TimedBlobQuantitiesMap> AllTimedBlobQuantitiesMap;

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {

	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Find blobs
	bool fFindBlobs;

	// Connectivity file
	std::string strConnectivity;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Data is regional
	bool fRegional;

	// Output file
	std::string strOutputFile;

	// Output nodefile
	std::string strOutputNodeFile;

	// Input variable
	std::string strInputVariable;

	// Variables to sum over
	std::string strSumVariables;

	// Variables to use in weighting
	std::string strWeightVariable;

	// Summary quantities
	std::string strOutputQuantities;

	// Format BlobStats output as a nodefile
	bool fOutputNodefile;

	// When outputting BlobStats as a nodefile write out seconds instead of hours
	bool fOutputSeconds;

	// Output headers
	bool fOutputHeaders;

	// Output the full times rather than time indexes
	bool fOutputFullTimes;

	// Time filter
	std::string strTimeFilter;

	// Start time
	std::string strStartTime;

	// End time
	std::string strEndTime;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Display help message
	bool fHelp;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_file", "");
		CommandLineString(strInputFileList, "in_list", "");
		CommandLineBool(fFindBlobs, "findblobs");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");
		CommandLineString(strOutputFile, "out_file", "");
		CommandLineString(strInputVariable, "var", "object_id");
		CommandLineString(strSumVariables, "sumvar", "");
		CommandLineString(strWeightVariable, "wtvar", "");
		CommandLineString(strOutputQuantities, "out", "");
		CommandLineBool(fOutputNodefile, "out_nodefile");
		CommandLineBool(fOutputSeconds, "out_seconds");
		CommandLineBool(fOutputHeaders, "out_headers");
		CommandLineBool(fOutputFullTimes, "out_fulltime");
		CommandLineString(strTimeFilter, "timefilter", "");
		CommandLineString(strStartTime, "time_start", "");
		CommandLineString(strEndTime, "time_end", "");

		CommandLineString(strLatitudeName, "latname", "lat");
		CommandLineString(strLongitudeName, "lonname", "lon");

		CommandLineBool(fHelp, "help");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if ((strInputFile == "") && (strInputFileList == "")) {
		_EXCEPTIONT("No input file (--in) or (--in_list) specified");
	}
	if ((strInputFile != "") && (strInputFileList != "")) {
		_EXCEPTIONT("Only one input file (--in) or (--in_list) allowed");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out_file) specified");
	}

	// Check variable
	if (strInputVariable == "") {
		_EXCEPTIONT("No variable name (--var) specified");
	}

	// Check output variable
	if ((strOutputQuantities == "") && (!fOutputNodefile)) {
		_EXCEPTIONT("No output quantities (--out) specified");
	}

	// Nodefile output
	if (fOutputSeconds && !fOutputNodefile) {
		_EXCEPTIONT("--out_seconds must be in combination with --out_nodefile");
	}

	// Input file list
	FilenameList vecInputFiles;

	if (strInputFile != "") {
		vecInputFiles.push_back(strInputFile);
	} else {
		vecInputFiles.FromFile(strInputFileList);
	}

	int nFiles = vecInputFiles.size();

	// Define the SimpleGrid for the input
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[0]);

		grid.GenerateLatitudeLongitude(
			vecNcFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}
	if (!grid.HasAreas()) {
		_EXCEPTIONT("Grid is missing area information needed by BlobStats");
	}

	// Parse the list of output quantities
	bool fTwoPassCalculation = false;
	std::vector<BlobQuantities::OutputQuantity> vecOutputVars;
	if (strOutputQuantities != "") {
		// Parse output quantities and store in vecOutputVars
		AnnounceStartBlock("Parsing output quantities");

		int iLast = 0;
		for (int i = 0; i <= strOutputQuantities.length(); i++) {

			if ((i == strOutputQuantities.length()) ||
				(strOutputQuantities[i] == ',') ||
				(strOutputQuantities[i] == ';')
			) {
				std::string strSubStr =
					strOutputQuantities.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecOutputVars.size());
				vecOutputVars.resize(iNextOp + 1);

				if (strSubStr == "minlat") {
					vecOutputVars[iNextOp] = BlobQuantities::MinLat;
				} else if (strSubStr == "maxlat") {
					vecOutputVars[iNextOp] = BlobQuantities::MaxLat;
				} else if (strSubStr == "minlon") {
					vecOutputVars[iNextOp] = BlobQuantities::MinLon;
				} else if (strSubStr == "maxlon") {
					vecOutputVars[iNextOp] = BlobQuantities::MaxLon;
				} else if (strSubStr == "meanlon") {
					vecOutputVars[iNextOp] = BlobQuantities::MeanLon;
				} else if (strSubStr == "meanlat") {
					vecOutputVars[iNextOp] = BlobQuantities::MeanLat;
				} else if (strSubStr == "centlat") {
					vecOutputVars[iNextOp] = BlobQuantities::CentroidLat;
				} else if (strSubStr == "centlon") {
					vecOutputVars[iNextOp] = BlobQuantities::CentroidLon;
				} else if (strSubStr == "area") {
					vecOutputVars[iNextOp] = BlobQuantities::Area;
				} else if (strSubStr == "wtcentlat") {
					vecOutputVars[iNextOp] = BlobQuantities::WtCentroidLat;
				} else if (strSubStr == "wtcentlon") {
					vecOutputVars[iNextOp] = BlobQuantities::WtCentroidLon;
				} else if (strSubStr == "wtarea") {
					vecOutputVars[iNextOp] = BlobQuantities::WtArea;
				} else if (strSubStr == "wtpcawidth") {
					vecOutputVars[iNextOp] = BlobQuantities::WtPCAWidth;
					fTwoPassCalculation = true;
				} else if (strSubStr == "wtpcalength") {
					vecOutputVars[iNextOp] = BlobQuantities::WtPCALength;
					fTwoPassCalculation = true;
				} else {
					_EXCEPTIONT("Invalid output quantity:  Expected\n"
						"[minlat, maxlat, minlon, maxlon, meanlat, meanlon,"
						" centlat, centlon, wtarea, wtcentlat, wtcentlon,"
						" wtpcawidth, wtpcalength, area]");
				}

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Check output variables when fOutputNodefile is active
	if (fOutputNodefile) {
		if (vecOutputVars.size() == 0) {
			vecOutputVars.push_back(BlobQuantities::CentroidLon);
			vecOutputVars.push_back(BlobQuantities::CentroidLat);

		} else if (vecOutputVars.size() == 1) {
			_EXCEPTIONT("When using --out_nodefile, the first two arguments of --out must be a longitude and latitude");

		} else {
			if ((vecOutputVars[0] == BlobQuantities::MeanLon) &&
			    (vecOutputVars[1] == BlobQuantities::MeanLat)
			) {
			} else if (
			    (vecOutputVars[0] == BlobQuantities::CentroidLon) &&
			    (vecOutputVars[1] == BlobQuantities::CentroidLat)
			) {
			} else if (
			    (vecOutputVars[0] == BlobQuantities::WtCentroidLon) &&
			    (vecOutputVars[1] == BlobQuantities::WtCentroidLat)
			) {
			} else {
				_EXCEPTIONT("When using --out_nodefile, the first two arguments of --out must be a longitude and latitude");
			}
		}
	}

	// VariableRegistry
	VariableRegistry varreg;

	// Parse --var argument
	VariableIndex varixBlobTag = varreg.FindOrRegister(strInputVariable);

	// Parse --sumvar argument
	std::vector<std::string> vecSumVarNames;
	std::vector<VariableIndex> vecSumVarIx;
	if (strSumVariables != "") {
		std::string strSumVariablesTemp = strSumVariables;
		for (;;) {
			VariableIndex varixSum;
			int iLast = varreg.FindOrRegisterSubStr(strSumVariablesTemp, &varixSum) + 1;

			if (varixSum == varixBlobTag) {
				_EXCEPTIONT("--var and --sumvar cannot contain the same variables");
			}

			vecSumVarIx.push_back(varixSum);
			vecSumVarNames.push_back( strSumVariablesTemp.substr(0,iLast-1) );
			if (iLast >= strSumVariablesTemp.length()) {
				break;
			}
			strSumVariablesTemp = strSumVariablesTemp.substr(iLast);
		}
	}
	_ASSERT(vecSumVarNames.size() == vecSumVarIx.size());

	// Parse --wtvar argument
	VariableIndex varixWeight = InvalidVariableIndex;
	if (strWeightVariable != "") {
		varixWeight = varreg.FindOrRegister(strWeightVariable);
	}

	// Parse --timefilter
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

	// Start time and end time for filtering
	Time timeStartTime;
	Time timeEndTime;

	// Time index across all files
	int iTime = 0;

	// Open output file
	FILE * fpout = fopen(strOutputFile.c_str(), "w");
	if (fpout == NULL) {
		_EXCEPTIONT("Error opening output file");
	}

	// Output header
	if (fOutputHeaders) {
		if (fOutputFullTimes) {
			fprintf(fpout, "time,%s\n", strOutputQuantities.c_str());
		} else {
			fprintf(fpout, "time_ix,%s\n", strOutputQuantities.c_str());
		}
	}

	// Loop through all files
	int iBlobIndex = 0;

	int iBlobOutputIndex = 0;
	int iBlobOutputTime = 0;

	for (int f = 0; f < nFiles; f++) {
		NcFileVector vecInFiles;
		vecInFiles.ParseFromString(vecInputFiles[f].c_str());

		if (vecInFiles.size() == 0) {
			_EXCEPTION1("Unable to open files in \"%s\"", vecInputFiles[f].c_str());
		}

		AnnounceStartBlock("File \"%s\"", vecInputFiles[f].c_str());

		// Parse --time_start and --time_end
		if ((f == 0) && ((strStartTime != "") || (strEndTime != ""))) {
			NcVar * varInTime = NcGetTimeVariable(*(vecInFiles[0]));
			if (varInTime == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"time\"",
					vecInFiles.GetFilename(0).c_str());
			}

			std::string strTimeCalendar;
			NcAtt * attCalendar = varInTime->get_att("calendar");
			if (attCalendar == NULL) {
				Announce("Dataset \"time\" does not specify \"calendar\" attribute.  Assuming \"standard\".");
				strTimeCalendar = "standard";
			} else {
				strTimeCalendar = attCalendar->as_string(0);
			}

			Time::CalendarType caltype = Time::CalendarTypeFromString(strTimeCalendar);
			if (caltype == Time::CalendarUnknown) {
				_EXCEPTION1("Unknown calendar name \"%s\"; cannot determine number of days per year", strTimeCalendar.c_str());
			}

			if (strStartTime != "") {
				timeStartTime = Time(caltype);
				timeStartTime.FromFormattedString(strStartTime);
			}
			if (strEndTime != "") {
				timeEndTime = Time(caltype);
				timeEndTime.FromFormattedString(strEndTime);
			}
		}

		// Computed quantities associated with each blob
		AllTimedBlobQuantitiesMap mapAllQuantities;

		// Time objects for each time in this file
		NcTimeDimension vecFileTimes;

		// Load in file times
		ReadCFTimeDataFromNcFile(vecInFiles[0], vecInputFiles[f], vecFileTimes, true);
		int nLocalTimes = vecFileTimes.size();

		// Load in indicator variable and validate
		Variable & varBlobTag = varreg.Get(varixBlobTag);

		// Load in indicator variable and validate
		Variable * pvarWeight = NULL;
		if (varixWeight != InvalidVariableIndex) {
			pvarWeight = &(varreg.Get(varixWeight));
		}

		// Blob index data
		DataArray1D<int> dataIndex(grid.GetSize());

/*
		if (varIndicator == NULL) {
			_EXCEPTION1("Unable to load variable \"%s\"",
				strInputVariable.c_str());
		}
		if (varIndicator->num_dims() != 1 + grid.DimCount()) {
			_EXCEPTION2("Incorrect number of dimensions for \"%s\""
				" (%i expected)", strInputVariable.c_str(), 1 + grid.DimCount());
		}
		if (grid.DimCount() == 1) {
			if (varIndicator->get_dim(1)->size() != grid.GetSize()) {
				_EXCEPTION2("Variable size mismatch:  Grid size %lu, variable size %lu",
					grid.GetSize(), varIndicator->get_dim(1)->size());
			}
		} else if (grid.DimCount() == 2) {
			if ((varIndicator->get_dim(1)->size() != grid.m_nGridDim[0]) ||
			    (varIndicator->get_dim(2)->size() != grid.m_nGridDim[1])
			) {
				_EXCEPTION4("Variable size mismatch:  Grid size (%lu,%lu), variable size (%lu,%lu)",
					grid.m_nGridDim[0],
					grid.m_nGridDim[1],
					varIndicator->get_dim(1)->size(),
					varIndicator->get_dim(2)->size());
			}

		} else {
			_EXCEPTION();
		}
*/
		// Loop through all times
		int iFirstFileTime = iTime;
		for (int t = 0; t < nLocalTimes; t++, iTime++) {

#ifndef TEMPEST_NOREGEX
			if (strTimeFilter != "") {
				std::string strTime = vecFileTimes[t].ToString();
				std::smatch match;
				if (!std::regex_search(strTime, match, reTimeSubset)) {
					Announce("Time %s (skipping)", vecFileTimes[t].ToString().c_str());
					continue;
				}
			}
#endif
			if (strStartTime != "") {
				double dDeltaSeconds = timeStartTime - vecFileTimes[t];
				if (dDeltaSeconds > 0.0) {
					Announce("Time %s (before start time; skipping)",
						vecFileTimes[t].ToString().c_str());
					continue;
				}
			}
			if (strEndTime != "") {
				double dDeltaSeconds = vecFileTimes[t] - timeEndTime;
				if (dDeltaSeconds > 0.0) {
					Announce("Time %s (after end time; skipping)",
						vecFileTimes[t].ToString().c_str());
					continue;
				}
			}

			Announce("Time %s", vecFileTimes[t].ToString().c_str());

			// Set the time index
			vecInFiles.SetConstantTimeIx(t);

			// Load in the data
			varBlobTag.LoadGridData(varreg, vecInFiles, grid);
			DataArray1D<float> & dataIndex = varBlobTag.GetData();
			_ASSERT(dataIndex.GetRows() == grid.GetSize());

			// Load in the weight data (if specified)
			DataArray1D<float> * pdataWeight = NULL;
			if (pvarWeight != NULL) {
				pvarWeight->LoadGridData(varreg, vecInFiles, grid);
				pdataWeight = &(pvarWeight->GetData());
				_ASSERT(pdataWeight->GetRows() == grid.GetSize());
			}

			// If --findblobs then assign blob indices
			if (fFindBlobs) {
				std::set<int> setVisited;
				std::queue<int> queueToVisit;
				for (int i = 0; i < grid.GetSize(); i++) {
					if (dataIndex[i] == 0) {
						continue;
					}
					if (setVisited.find(i) != setVisited.end()) {
						continue;
					}

					iBlobIndex++;

					queueToVisit.push(i);
					while (!queueToVisit.empty()) {
						int ixNode = queueToVisit.front();
						queueToVisit.pop();
						if (setVisited.find(ixNode) != setVisited.end()) {
							continue;
						}
						setVisited.insert(ixNode);

						dataIndex[ixNode] = iBlobIndex;

						for (int k = 0; k < grid.m_vecConnectivity[ixNode].size(); k++) {
							int ixNeighbor = grid.m_vecConnectivity[ixNode][k];
							if (dataIndex[ixNeighbor] != 0) {
								queueToVisit.push(ixNeighbor);
							}
						}
					}
				}
			}

			// Load in all --sumvar data
			std::vector< DataArray1D<float> * > vecSumVarData;
			for (int v = 0; v < vecSumVarIx.size(); v++) {
				Variable & varSumVar = varreg.Get(vecSumVarIx[v]);
				varSumVar.LoadGridData(varreg, vecInFiles, grid);
				DataArray1D<float> & dataSumVar = varSumVar.GetData();
				vecSumVarData.push_back(&dataSumVar);
				_ASSERT(dataIndex.GetRows() == grid.GetSize());
			}

			// Last data index
			int iLastDataIndex = 0;

			// Iterator to BlobQuantities with iLastDataIndex and iTime
			TimedBlobQuantitiesMap::iterator iterBlobQuantities;

			// Loop over all locations
			for (int i = 0; i < grid.GetSize(); i++) {

				// Ignore non-blob data
				if (dataIndex[i] == 0) {
					continue;
				}

				// Check if iterator already points to correct data
				if (dataIndex[i] != iLastDataIndex) {

					// Iterator to TimedBlobQuantities with iLastDataIndex
					AllTimedBlobQuantitiesMap::iterator
						iterTimedBlobQuantities;

					// Get new iterator for this dataIndex
					iterTimedBlobQuantities =
						mapAllQuantities.find(
							dataIndex[i]);

					if (iterTimedBlobQuantities == mapAllQuantities.end()) {

						std::pair<AllTimedBlobQuantitiesMap::iterator,bool> pr =
							mapAllQuantities.insert(
								AllTimedBlobQuantitiesMap::value_type(
									dataIndex[i],
									TimedBlobQuantitiesMap()));

						if (pr.second == false) {
							_EXCEPTIONT("Map insertion failed");
						}
						iterTimedBlobQuantities = pr.first;
					}

					// Get the BlobQuantities at the current time
					TimedBlobQuantitiesMap & mapTimedBlobQuantities =
						iterTimedBlobQuantities->second;

					iterBlobQuantities =
						mapTimedBlobQuantities.find(iTime);

					if (iterBlobQuantities == mapTimedBlobQuantities.end()) {
						std::pair<TimedBlobQuantitiesMap::iterator,bool> pr =
							mapTimedBlobQuantities.insert(
								TimedBlobQuantitiesMap::value_type(
									iTime,
									BlobQuantities(vecSumVarIx.size())));

						if (pr.second == false) {
							_EXCEPTIONT("Map insertion failed");
						}
						iterBlobQuantities = pr.first;
					}

					// Update last data index
					iLastDataIndex = dataIndex[i];
				}

				// Associated BlobQuantities
				BlobQuantities & bq = iterBlobQuantities->second;

				// Insert point into array
				bq.box.insert(
					grid.m_dLat[i],
					LonRadToStandardRange(grid.m_dLon[i]));

				// Add blob area
				bq.dArea +=
					grid.m_dArea[i];

				// Add area-weighted 3D coordinates
				double dX, dY, dZ;
				RLLtoXYZ_Rad(grid.m_dLon[i], grid.m_dLat[i], dX, dY, dZ);
				bq.dAreaX += dX * grid.m_dArea[i];
				bq.dAreaY += dY * grid.m_dArea[i];
				bq.dAreaZ += dZ * grid.m_dArea[i];

				// Add area and variable-weighted 3D coordinates
				if (pdataWeight != NULL) {
					double dWt = grid.m_dArea[i] * (*pdataWeight)[i];
					bq.dWtArea += dWt;
					bq.dWtAreaX += dX * dWt;
					bq.dWtAreaY += dY * dWt;
					bq.dWtAreaZ += dZ * dWt;
				}

				// Add sum variables
				for (int v = 0; v < vecSumVarIx.size(); v++) {
					bq.dSumVars[v] += static_cast<double>((*(vecSumVarData[v]))[i]) * grid.m_dArea[i] * EarthRadius * EarthRadius;
				}
			}

			// Calculate weighted centroid
			if (pdataWeight != NULL) {
				for (auto & iterTimedBlobQuantities : mapAllQuantities) {
					auto iterBlobQuantities = iterTimedBlobQuantities.second.find(iTime);

					if (iterBlobQuantities != iterTimedBlobQuantities.second.end()) {
						BlobQuantities & bq = iterBlobQuantities->second;

						bq.dWtAreaX /= bq.dWtArea;
						bq.dWtAreaY /= bq.dWtArea;
						bq.dWtAreaZ /= bq.dWtArea;

						double dMag = sqrt(
							  bq.dWtAreaX * bq.dWtAreaX
							+ bq.dWtAreaY * bq.dWtAreaY
							+ bq.dWtAreaZ * bq.dWtAreaZ);

						bq.dWtAreaX /= dMag;
						bq.dWtAreaY /= dMag;
						bq.dWtAreaZ /= dMag;

						XYZtoRLL_Deg(bq.dWtAreaX, bq.dWtAreaY, bq.dWtAreaZ, bq.dWtCentLonDeg, bq.dWtCentLatDeg);
						//printf("%i %i %1.5f %1.5f\n", iterTimedBlobQuantities.first, iTime, bq.dWtCentLonDeg, bq.dWtCentLatDeg);
					}
				}
			}

			// If a two-pass calculation is needed then loop through all locations again
			if (fTwoPassCalculation) {
				for (int i = 0; i < grid.GetSize(); i++) {

					// Ignore non-blob data
					if (dataIndex[i] == 0) {
						continue;
					}

					// Check if iterator already points to correct data
					if (dataIndex[i] != iLastDataIndex) {

						// Iterator to TimedBlobQuantities with iLastDataIndex
						AllTimedBlobQuantitiesMap::iterator
							iterTimedBlobQuantities;

						// Get new iterator for this dataIndex
						iterTimedBlobQuantities =
							mapAllQuantities.find(
								dataIndex[i]);

						_ASSERT(iterTimedBlobQuantities != mapAllQuantities.end());

						// Get the BlobQuantities at the current time
						TimedBlobQuantitiesMap & mapTimedBlobQuantities =
							iterTimedBlobQuantities->second;

						iterBlobQuantities =
							mapTimedBlobQuantities.find(iTime);

						_ASSERT(iterBlobQuantities != mapTimedBlobQuantities.end());

						// Update last data index
						iLastDataIndex = dataIndex[i];
					}

					// Associated BlobQuantities
					BlobQuantities & bq = iterBlobQuantities->second;

					// Correlation matrix components
					//double dXs = DegToRad(bq.dWtCentLonDeg) - grid.m_dLon[i];
					//double dYs = DegToRad(bq.dWtCentLatDeg) - grid.m_dLat[i];

					double dXs;
					double dYs;
					StereographicProjection(
						DegToRad(bq.dWtCentLonDeg), DegToRad(bq.dWtCentLatDeg),
						grid.m_dLon[i], grid.m_dLat[i],
						dXs, dYs);

					//printf("%1.5f %1.5f\n", dXs, dYs);

					double dWt = grid.m_dArea[i] * (*pdataWeight)[i];

					bq.dCorrComp[0] += dWt * dXs * dXs;
					bq.dCorrComp[1] += dWt * dXs * dYs;
					bq.dCorrComp[2] += dWt * dYs * dYs;
				}
			}
		}

		// Output all BlobQuantities
		{
			Announce("Writing blob data");

			AllTimedBlobQuantitiesMap::iterator iterBlobs =
				mapAllQuantities.begin();

			for (; iterBlobs != mapAllQuantities.end(); iterBlobs++) {

				if (iterBlobs->second.size() == 0) {
					continue;
				}

				TimedBlobQuantitiesMap::iterator iterTimes =
					iterBlobs->second.begin();

				if (fOutputNodefile) {
					const Time & timeStart = vecFileTimes[iterTimes->first - iFirstFileTime];
					if (fOutputSeconds) {
						fprintf(fpout, "start\t%lu\t%i\t%i\t%i\t%05i\n",
							iterBlobs->second.size(),
							timeStart.GetYear(),
							timeStart.GetMonth(),
							timeStart.GetDay(),
							timeStart.GetSecond());
					} else {
						fprintf(fpout, "start\t%lu\t%i\t%i\t%i\t%i\n",
							iterBlobs->second.size(),
							timeStart.GetYear(),
							timeStart.GetMonth(),
							timeStart.GetDay(),
							timeStart.GetSecond() / 3600);
					}
				}

				if (iterBlobs->first != iBlobOutputIndex) {
					iBlobOutputTime = 0;
					iBlobOutputIndex = iterBlobs->first;
				}

				for (; iterTimes != iterBlobs->second.end(); iterTimes++) {

					if (!fOutputNodefile) {
						fprintf(fpout, "%i\t%i\t", iterBlobs->first, iBlobOutputTime);
					}
					iBlobOutputTime++;

					_ASSERT(iterTimes->first-iFirstFileTime < vecFileTimes.size());
					Time & timeCurrent = vecFileTimes[iterTimes->first-iFirstFileTime];

					if (fOutputNodefile) {
						fprintf(fpout, "\t0\t0");
					} else {
						if (fOutputFullTimes) {
							fprintf(fpout, "\"%s\"", vecFileTimes[iterTimes->first - iFirstFileTime].ToString().c_str());
						} else {
							fprintf(fpout, "%i", iterTimes->first);
						}
					}

					// BlobQuantities for this blob at this timestep
					BlobQuantities & bq = iterTimes->second;

					// Bounding box coordinates
					double dLatDeg0 = RadToDeg(bq.box.lat[0]);
					double dLatDeg1 = RadToDeg(bq.box.lat[1]);
					double dLonDeg0 = RadToDeg(bq.box.lon[0]);
					double dLonDeg1 = RadToDeg(bq.box.lon[1]);

					// Calculate centroid
					double dCentX = bq.dAreaX / bq.dArea;
					double dCentY = bq.dAreaY / bq.dArea;
					double dCentZ = bq.dAreaZ / bq.dArea;

					double dCentMag =
						sqrt(dCentX * dCentX + dCentY * dCentY + dCentZ * dCentZ);

					if (dCentMag == 0.0) {
						dCentX = 0.0;
						dCentY = 0.0;
						dCentZ = 1.0;
					} else {
						dCentX /= dCentMag;
						dCentY /= dCentMag;
						dCentZ /= dCentMag;
					}

					double dCentLonDeg;
					double dCentLatDeg;
					XYZtoRLL_Deg(dCentX, dCentY, dCentZ, dCentLonDeg, dCentLatDeg);

					// Length and width using PCA method
					double dWtPCAWidthDeg = 0.0;
					double dWtPCALengthDeg = 0.0;
					if (bq.dCorrComp[0] != 0.0) {
						double dCorrComp0 = bq.dCorrComp[0] / bq.dWtArea;
						double dCorrComp1 = bq.dCorrComp[1] / bq.dWtArea;
						double dCorrComp2 = bq.dCorrComp[2] / bq.dWtArea;

						double dDisc = dCorrComp0 - dCorrComp2;
						dDisc = sqrt(4.0 * dCorrComp1 * dCorrComp1 + dDisc * dDisc);
						dWtPCALengthDeg = sqrt(0.5 * (dCorrComp0 + dCorrComp2 + dDisc));
						dWtPCAWidthDeg = sqrt(0.5 * (dCorrComp0 + dCorrComp2 - dDisc));

						dWtPCALengthDeg = RadToDeg(2.0 * atan(dWtPCALengthDeg));
						dWtPCAWidthDeg = RadToDeg(2.0 * atan(dWtPCAWidthDeg));
					}

					// Write standard outputs
					for (int i = 0; i < vecOutputVars.size(); i++) {

						// Minimum latitude
						if (vecOutputVars[i] == BlobQuantities::MinLat) {
							fprintf(fpout, "\t%1.6f", dLatDeg0);
						}

						// Maximum latitude
						if (vecOutputVars[i] == BlobQuantities::MaxLat) {
							fprintf(fpout, "\t%1.6f", dLatDeg1);
						}

						// Minimum longitude
						if (vecOutputVars[i] == BlobQuantities::MinLon) {
							fprintf(fpout, "\t%1.6f", dLonDeg0);
						}

						// Maximum longitude
						if (vecOutputVars[i] == BlobQuantities::MaxLon) {
							fprintf(fpout, "\t%1.6f", dLonDeg1);
						}

						// Mean latitude
						if (vecOutputVars[i] == BlobQuantities::MeanLat) {
							double dMidLatDeg = 0.5 * (dLatDeg0 + dLatDeg1);

							fprintf(fpout, "\t%1.6f", dMidLatDeg);
						}

						// Mean longitude
						if (vecOutputVars[i] == BlobQuantities::MeanLon) {
							double dMidLonDeg;
							if (dLonDeg0 <= dLonDeg1) {
								dMidLonDeg = 0.5 * (dLonDeg0 + dLonDeg1);
							} else {
								dMidLonDeg = 0.5 * (dLonDeg0 + dLonDeg1 + 360.0);
								if (dMidLonDeg > 360.0) {
									dMidLonDeg -= 360.0;
								}
							}

							fprintf(fpout, "\t%1.6f", dMidLonDeg);
						}

						// Centroid latitude
						if (vecOutputVars[i] == BlobQuantities::CentroidLat) {
							fprintf(fpout, "\t%1.6f", dCentLatDeg);
						}

						// Centroid longitude
						if (vecOutputVars[i] == BlobQuantities::CentroidLon) {
							fprintf(fpout, "\t%1.6f", dCentLonDeg);
						}

						// Area
						if (vecOutputVars[i] == BlobQuantities::Area) {
							fprintf(fpout, "\t%1.6e", bq.dArea * EarthRadius * EarthRadius);
						}

						// Weighted area
						if (vecOutputVars[i] == BlobQuantities::WtArea) {
							fprintf(fpout, "\t%1.6e", bq.dWtArea * EarthRadius * EarthRadius);
						}

						// Weighted centroid latitude
						if (vecOutputVars[i] == BlobQuantities::WtCentroidLat) {
							fprintf(fpout, "\t%1.6f", bq.dWtCentLatDeg);
						}

						// Weighted centroid longitude
						if (vecOutputVars[i] == BlobQuantities::WtCentroidLon) {
							fprintf(fpout, "\t%1.6f", bq.dWtCentLonDeg);
						}

						// Weighted Width
						if (vecOutputVars[i] == BlobQuantities::WtPCAWidth) {
							fprintf(fpout, "\t%1.6f", dWtPCAWidthDeg);
						}

						// Weighted Length
						if (vecOutputVars[i] == BlobQuantities::WtPCALength) {
							fprintf(fpout, "\t%1.6f", dWtPCALengthDeg);
						}
					}

					// Sum variables
					for (int v = 0; v < bq.dSumVars.size(); v++) {
						fprintf(fpout,"\t%1.6e", bq.dSumVars[v]);
					}

					// Times
					if (fOutputNodefile) {
						if (fOutputSeconds) {
							fprintf(fpout, "\t%i\t%i\t%i\t%05i",
								timeCurrent.GetYear(),
								timeCurrent.GetMonth(),
								timeCurrent.GetDay(),
								timeCurrent.GetSecond());
						} else {
							fprintf(fpout, "\t%i\t%i\t%i\t%i",
								timeCurrent.GetYear(),
								timeCurrent.GetMonth(),
								timeCurrent.GetDay(),
								timeCurrent.GetSecond() / 3600);
						}
					}

					// Endline
					fprintf(fpout, "\n");
				}
			}
		}
		AnnounceEndBlock("Done");
	}

	fclose(fpout);

	// Done
	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

