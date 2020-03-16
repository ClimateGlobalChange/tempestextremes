///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileCompose.cpp
///	\author  Paul Ullrich
///	\version March 13, 2020
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

	// Output grid type
	std::string strOutputGrid;

	// Output data file
	std::string strOutputData;

	// Start of time range to use
	std::string strTimeBegin;

	// End of time range to use
	std::string strTimeEnd;

	// List of variables to output
	std::string strVariables;

	// Grid spacing of output (Cartesian great-circle or radial distance)
	double dDeltaX;

	// Resolution of output (Cartesian grid size or number of radial points)
	int nResolutionX;

	// Resolution of the output (azimuthal points)
	int nResolutionA;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_nodefile", "");
		CommandLineStringD(strInputNodeFileType, "in_nodefile_type", "SN", "[DCU|SN]");
		CommandLineString(strInputFormat, "in_fmt", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");

		CommandLineStringD(strOutputGrid, "out_grid", "XY", "[XY|RAD]");
		CommandLineString(strOutputData, "out_data", "");

		//CommandLineString(strTimeBegin, "time_begin", "");
		//CommandLineString(strTimeEnd, "time_end", "");

		CommandLineString(strVariables, "var", "");

		CommandLineDouble(dDeltaX, "dx", 0.5);
		CommandLineInt(nResolutionX, "resx", 11);
		CommandLineInt(nResolutionA, "resa", 10);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varregIn;

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
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}
	if (strVariables.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}

	// Output grid options
	STLStringHelper::ToLower(strOutputGrid);
	if ((strOutputGrid != "xy") && (strOutputGrid != "rad")) {
		_EXCEPTIONT("Grid type of output (--out_grid) must be \"xy\" or \"rad\"");
	}
	if (dDeltaX <= 0.0) {
		_EXCEPTIONT("Grid spacing of output (--dx) must be nonnegative");
	}
	if (nResolutionX < 1) {
		_EXCEPTIONT("Resolution of output (--resx) must be nonnegative");
	}
	if ((strOutputGrid == "rad") && (nResolutionA < 8)) {
		_EXCEPTIONT("Resolution of output (--resa) must be >= 8");
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
			int iLast = varregIn.FindOrRegisterSubStr(strVariablesTemp, &varix) + 1;

			vecVarIx.push_back(varix);
			vecVarNames.push_back( strVariablesTemp.substr(0,iLast-1) );
			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}
	}
	_ASSERT(vecVarNames.size() == vecVarIx.size());

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

	// Open output file
	NcFile ncoutfile(strOutputData.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputData.c_str());
	}

	std::vector< std::vector<int> > vecAuxDimSize;
	std::vector< NcVar * > vecNcVarOut;

	// Loop over all files
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

		// Load in vector of Times from the file
		std::vector<Time> vecTimes;
		ReadCFTimeDataFromNcFile(&ncinfile, vecInputFileList[f], vecTimes, true);

		// Initialize the output file
		if (f == 0) {
			if (strOutputGrid == "xy") {
				NcDim * dimX = ncoutfile.add_dim("x", nResolutionX);
				NcDim * dimY = ncoutfile.add_dim("y", nResolutionX);

			} else if (strOutputGrid == "rad") {
				NcDim * dimR = ncoutfile.add_dim("r", nResolutionX);
				NcDim * dimA = ncoutfile.add_dim("az", nResolutionA);

			} else {
				_EXCEPTIONT("Invalid grid");
			}
		}

		// Loop through all times
		for (int t = 0; t < vecTimes.size(); t++) {

			AnnounceStartBlock("Processing time \"%s\"",
				vecTimes[t].ToString().c_str());

			// Find any PathNodes at this time
			TimeToPathNodeMap::const_iterator iter =
				mapTimeToPathNode.find(vecTimes[t]);

			// No PathNodes at this time; continue to next file
			if (iter == mapTimeToPathNode.end()) {
				continue;
			}

			// Loop through all variables
			for (int v = 0; v < vecVarIx.size(); v++) {

				const std::string & strVariable = vecVarNames[v];

				Announce("Processing variable \"%s\"", strVariable.c_str());
/*
				// Load input and output variables
				NcVar * varIn = vecNcVarIn[v];
				NcVar * varOut = vecNcVarOut[v];

				// Load in data
				int nGridDims = grid.m_nGridDim.size();

				std::vector<NcDim *> vecDim;
				for (int d = 0; d < varIn->num_dims(); d++) {
					vecDim.push_back(varIn->get_dim(d));
				}
*/
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
/*
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

					// Do stuff

					// Write data
					varOut->set_cur(&(vecDataPos[0]));
					varOut->put(&(data[0]), &(vecDataSize[0]));
				}
*/
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


