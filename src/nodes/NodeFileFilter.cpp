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
#include "DataMatrix.h"
#include "ArgumentTree.h"
#include "STLStringHelper.h"
#include "NodeFileUtilities.h"
#include "NetCDFUtilities.h"

#include "netcdfcpp.h"

#include <fstream>
#include <queue>
#include <set>
/*
#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif
*/
///////////////////////////////////////////////////////////////////////////////

void BuildFilter(
	const SimpleGrid & grid,
	const ColumnDataHeader & cdh,
	const PathVector & pathvec,
	const PathNodeIndexVector & vecPathNodes,
	const std::string & strDist,
	DataVector<double> & dataFilter
) {
	dataFilter.Zero();

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
			dataFilter[ix] = 1.0;

			// Add all neighbors of this point
			for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
				queueToVisit.push(grid.m_vecConnectivity[ix][n]);
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

	// List of variables to preserve
	std::string strPreserve;

	// Output data file
	std::string strOutputData;

	// Output list of data files
	std::string strOutputDataList;

	// Filter variables by distance
	std::string strFilterByDist;

	// Only output the mask
	bool fMaskOnly;

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

		CommandLineString(strPreserve, "preserve", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputDataList, "out_data_list", "");

		CommandLineString(strFilterByDist, "bydist", "");
		CommandLineBool(fMaskOnly, "mask_only");
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
	if (strFilterByDist.length() == 0) {
		_EXCEPTIONT("No filter command (--bydist) specified");
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

	// Parse --by_dist and check arguments
	ArgumentTree atFilterByDist(true);
	atFilterByDist.Parse(strFilterByDist);

	for (int v = 0; v < atFilterByDist.size(); v++) {
		ArgumentTree * patThisFilterByDist =
			atFilterByDist.GetSubTree(v);

		if (patThisFilterByDist == NULL) {
			_EXCEPTION1("Invalid format of --bydist \"%s\"; expected <variable>,<dist>",
				atFilterByDist.GetArgumentString(v).c_str());
		}
		if (patThisFilterByDist->size() != 2) {
			_EXCEPTION1("Invalid format of --bydist \"%s\"; expected <variable>,<dist>",
				atFilterByDist.GetArgumentString(v).c_str());
		}
		if (!STLStringHelper::IsFloat((*patThisFilterByDist)[1])) {
			bool fFound = false;
			for (int i = 0; i < cdhInput.size(); i++) {
				if (cdhInput[i] == (*patThisFilterByDist)[1]) {
					fFound = true;
					break;
				}
			}
			if (!fFound) {
				_EXCEPTION1("Invalid format of --bydist;"
					" second argument <dist> \"%s\" does not appear in --in_fmt",
					(*patThisFilterByDist)[1].c_str());
			}
		}
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
		AnnounceStartBlock("Building input data list");
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
	DataVector<double> data(grid.GetSize());

	// Create filter
	DataVector<double> dataFilter(grid.GetSize());

	// Loop over all files
	if (vecInputFileList.size() != vecOutputFileList.size()) {
		_EXCEPTION();
	}

	for (int f = 0; f < vecInputFileList.size(); f++) {

		AnnounceStartBlock("Processing input (%s)", vecInputFileList[f].c_str());

		// Open input file
		NcFile ncinfile(vecInputFileList[f].c_str(), NcFile::ReadOnly);
		if (!ncinfile.is_valid()) {
			_EXCEPTIONT("Unable to open input datafile");
		}

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
		}

		// Loop through all variables
		for (int v = 0; v < atFilterByDist.size(); v++) {
			ArgumentTree * patThisFilterByDist =
				atFilterByDist.GetSubTree(v);

			const std::string & strVariable = (*patThisFilterByDist)[0];
			const std::string & strDist = (*patThisFilterByDist)[1];

			AnnounceStartBlock("Processing variable \"%s\"", strVariable.c_str());

			// Create an output variable
			CopyNcVar(ncinfile, ncoutfile, strVariable, true, false);

			NcVar * varIn = ncinfile.get_var(strVariable.c_str());
			if (varIn == NULL) {
				_EXCEPTION();
			}

			NcVar * varOut = ncoutfile.get_var(strVariable.c_str());
			if (varOut == NULL) {
				_EXCEPTION();
			}

			// Load in data
			int nGridDims = grid.m_nGridDim.size();

			std::vector<NcDim *> vecDim;
			for (int d = 0; d < varIn->num_dims(); d++) {
				vecDim.push_back(varIn->get_dim(d));
			}
			if (vecDim.size() < 1 + nGridDims) {
				_EXCEPTION2("Insufficient dimensions in variable \"%s\" in file \"%s\"",
					strVariable.c_str(), vecOutputFileList[f].c_str());
			}
			if (strcmp(vecDim[0]->name(), "time") != 0) {
				_EXCEPTION2("First dimension of variable \"%s\" in file \"%s\" must be \"time\"",
					strVariable.c_str(), vecOutputFileList[f].c_str());
			}
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
			for (int d = 0; d < vecDim.size() - nGridDims; d++) {
				nAuxDims *= vecDim[d]->size();
			}

			DataVector<long> vecDataSize(vecDim.size());
			for (int d = 0; d < vecDim.size(); d++) {
				if (d >= vecDim.size() - nGridDims) {
					vecDataSize[d] = vecDim[d]->size();
				} else {
					vecDataSize[d] = 1;
				}
			}

			// Current time index
			int iCurrentTime = (-1);

			// Loop through all auxiliary dimensions
			DataVector<long> vecDataPos(vecDim.size());
			for (int i = 0; i < nAuxDims; i++) {

				// Load data
				int ixDim = i;
				for (int d = vecDim.size() - nGridDims - 1; d >= 0; d--) {
					vecDataPos[d] = ixDim % vecDim[d]->size();
					ixDim /= vecDim[d]->size();
				}
				if (ixDim != 0) {
					_EXCEPTIONT("Logic error");
				}

				varIn->set_cur(&(vecDataPos[0]));
				varIn->get(&(data[0]), &(vecDataSize[0]));

				// Check if the time has changed
				if (iCurrentTime != vecDataPos[0]) {
					iCurrentTime = vecDataPos[0];

					if (iCurrentTime >= vecTimes.size()) {
						_EXCEPTION2("Logic error %i %i",
							iCurrentTime, vecTimes.size());
					}

					Announce("Processing time \"%s\"",
						vecTimes[iCurrentTime].ToString().c_str());

					// Build filter
					TimeToPathNodeMap::const_iterator iter =
						mapTimeToPathNode.find(vecTimes[iCurrentTime]);

					if (iter == mapTimeToPathNode.end()) {
						dataFilter.Zero();
					} else {
						BuildFilter(
							grid,
							cdhInput,
							pathvec,
							iter->second,
							strDist,
							dataFilter);
					}
				}

				// Apply filter
				if (!fMaskOnly) {
					if (!fInvert) {
						for (int i = 0; i < data.GetRows(); i++) {
							data[i] *= dataFilter[i];
						}
					} else {
						for (int i = 0; i < data.GetRows(); i++) {
							data[i] *= (1.0 - dataFilter[i]);
						}
					}

				// Only output the filter mask
				} else {
					if (!fInvert) {
						for (int i = 0; i < data.GetRows(); i++) {
							data[i] = dataFilter[i];
						}
					} else {
						for (int i = 0; i < data.GetRows(); i++) {
							data[i] = 1.0 - dataFilter[i];
						}
					}
				}

				// Write data
				varOut->set_cur(&(vecDataPos[0]));
				varOut->put(&(data[0]), &(vecDataSize[0]));
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


