///////////////////////////////////////////////////////////////////////////////
///
///	\file    VariableProcessor.cpp
///	\author  Paul Ullrich
///	\version August 23, 2020
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
#include "FilenameList.h"
#include "NetCDFUtilities.h"
#include "TimeMatch.h"

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

	// List of variables to preserve
	std::string strPreserveVarName;

	// List of variables to process
	std::string strVariables;

	// List of variables to output
	std::string strOutputVariables;

	// Time filter
	std::string strTimeFilter;

	// Override _FillValue
	std::string strFillValue;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Log dir
	std::string strLogDir;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strGridFile, "in_grid", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputDataList, "out_data_list", "");

		CommandLineString(strPreserveVarName, "preserve", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputVariables, "varout", "");
		CommandLineString(strTimeFilter, "timefilter", "");

		CommandLineString(strFillValue, "fillvalue", "");

		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");
		CommandLineString(strLogDir, "logdir", ".");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Check arguments
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

	// Parse --fillvalue argument
	float dFillValue = 0.0;
	if (strFillValue != "") {
		if (!STLStringHelper::IsFloat(strFillValue)) {
			_EXCEPTIONT("Argument --fillvalue must be a floating point value");
		}
		dFillValue = std::stof(strFillValue);
	}

	// Parse --var argument
	std::vector<std::string> vecVariableNamesIn;
	std::vector<VariableIndex> vecVarIxIn;
	if (strVariables != "") {
		std::string strVariablesTemp = strVariables;
		for (;;) {
			VariableIndex varix;
			int iLast = varreg.FindOrRegisterSubStr(strVariablesTemp, &varix) + 1;

			vecVarIxIn.push_back(varix);
			vecVariableNamesIn.push_back( strVariablesTemp.substr(0,iLast-1) );
			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}
	}
	_ASSERT(vecVariableNamesIn.size() == vecVarIxIn.size());
 
	// Parse --varout argument
	std::vector<std::string> vecVariableNamesOut;
	if (strOutputVariables != "") {
		STLStringHelper::ParseVariableList(strOutputVariables, vecVariableNamesOut);
	} else {
		vecVariableNamesOut = vecVariableNamesIn;
	}
	if (vecVariableNamesOut.size() != vecVariableNamesIn.size()) {
		_EXCEPTION2("Inconsistent number of variables in --var (%lu) and --varout (%lu)",
			vecVarIxIn.size(), vecVariableNamesOut.size());
	}

	// Parse preserve list (--preserve)
	std::vector<std::string> vecPreserveVariableStrings;
	STLStringHelper::ParseVariableList(strPreserveVarName, vecPreserveVariableStrings);

	// Store input data
	FilenameList vecInputFileList;

	if (strInputData.length() != 0) {
		vecInputFileList.push_back(strInputData);

	} else {
		AnnounceStartBlock("Building input data list");
		vecInputFileList.FromFile(strInputDataList);
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

		// Get the list of data files
		NcFileVector vecFiles;
		if (strGridFile == "") {
			vecFiles.ParseFromString(vecInputFileList[0]);
		} else {
			vecFiles.ParseFromString(strGridFile);
		}

		grid.GenerateLatitudeLongitude(
			vecFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
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

		// Test Regex
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

	/////////////////////////////////////////////
	// 	Actual variable processing

	AnnounceStartBlock("Begin processing");

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

	// Loop over all files to be processed
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
		if (vecInputFileList.size() != 1) {
			AnnounceStartBlock("DataList line %lu/%lu", f+1, vecInputFileList.size());
		}

		// Get the list of data files
		NcFileVector vecFiles;
		vecFiles.ParseFromString(vecInputFileList[f]);

		_ASSERT(vecFiles.size() > 0);

		// Get time information
		bool fFileHasTime =
			(vecFiles.GetFileType(0) != NcFileVector::FileType_NoTimeDim) &&
			(vecFiles.GetFileType(0) != NcFileVector::FileType_NoTimeVar);

		const NcTimeDimension & vecTimes =
			vecFiles.GetNcTimeDimension(0);
/*
		// No need to read in CFTimeData again, since it's already done in NcFileVector now
		bool fFileHasTime = (NcGetTimeVariable(*(vecFiles[0])) != NULL);
		NcTimeDimension vecTimes;
		if (fFileHasTime) {
			ReadCFTimeDataFromNcFile(
				vecFiles[0],
				vecFiles.GetFilename(0),
				vecTimes,
				false);
		}
*/
		NcTimeDimension vecOutputTimes;

#ifndef TEMPEST_NOREGEX
		if (strTimeFilter == "") {
			vecOutputTimes = vecTimes;
		} else {
			vecOutputTimes.m_units = vecTimes.m_units;
			vecOutputTimes.m_nctype = vecTimes.m_nctype;
			for (int t = 0; t < vecTimes.size(); t++) {
				std::string strTime = vecTimes[t].ToString();
				std::smatch match;
				if (std::regex_search(strTime, match, reTimeSubset)) {
					vecOutputTimes.push_back(vecTimes[t]);
				}
			}
		}
#else
		vecOutputTimes = vecTimes;
#endif

		if (!fFileHasTime) {
			vecOutputTimes.push_back(Time(Time::CalendarNone));
		}

		// Open output file
		NcFile ncout(vecOutputFileList[f].c_str(), NcFile::Replace);
		if (!ncout.is_valid()) {
			_EXCEPTION1("Unable to open output datafile \"%s\"",
				vecOutputFileList[f].c_str());
		}

		// Write new time dimension
		NcDim * dimTime = NULL;
		if (fFileHasTime) {
			WriteCFTimeDataToNcFile(
				&ncout,
				vecOutputFileList[f],
				vecOutputTimes);

			dimTime = NcGetTimeDimension(ncout);
			_ASSERT(dimTime != NULL);
		}

		// Write grid dimensions
		NcDim * dim0 = NULL;
		NcDim * dim1 = NULL;
		if (grid.m_nGridDim.size() == 1) {
			dim0 = ncout.add_dim("ncol", grid.m_nGridDim[0]);
		} else if (grid.m_nGridDim.size() == 2) {
			dim0 = ncout.add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]);
			dim1 = ncout.add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]);

			CopyNcVarIfExists(
				*(vecFiles[0]),
				ncout,
				strLatitudeName,
				true,
				true);

			CopyNcVarIfExists(
				*(vecFiles[0]),
				ncout,
				strLongitudeName,
				true,
				true);

		} else {
			_EXCEPTIONT("Only 1D or 2D spatial data supported");
		}

		// Copy preserve variables
		if (strPreserveVarName != "") {
			AnnounceStartBlock("Preserving variables");

			for (int v = 0; v < vecPreserveVariableStrings.size(); v++) {
				Announce("Variable %s", vecPreserveVariableStrings[v].c_str());

				NcVar * var;
				size_t sFileIx = vecFiles.FindContainingVariable(vecPreserveVariableStrings[v], &var);
				if (var == NULL) {
					_EXCEPTION1("Unable to find variable \"%s\" in specified files",
						vecPreserveVariableStrings[v].c_str());
				}

				NcFile * ncinfile = vecFiles[sFileIx];
				_ASSERT(ncinfile != NULL);

				CopyNcVar(
					(*ncinfile),
					ncout,
					vecPreserveVariableStrings[v]);
			}

			AnnounceEndBlock("Done");
		}

		// Build up the processing queue
		varreg.ClearProcessingQueue();

		std::vector<NcVar *> vecNcVarOut(vecVarIxIn.size());
		std::vector<bool> vecNcVarOutWroteUnits(vecVarIxIn.size(), false);
		for (int v = 0; v < vecVarIxIn.size(); v++) {
			NcVar * varOut = NULL;

			DimInfoVector vecAuxDimInfo;
			varreg.AppendVariableToProcessingQueue(
				vecFiles,
				grid,
				vecVarIxIn[v],
				&vecAuxDimInfo);

			std::vector<NcDim *> vecDimOut;

			if (fFileHasTime) {
				vecDimOut.push_back(dimTime);
			}

			for (long d = 0; d < vecAuxDimInfo.size(); d++) {
				NcDim * dim =
					AddNcDimOrUseExisting(
						ncout,
						vecAuxDimInfo[d].name,
						vecAuxDimInfo[d].size);

				vecDimOut.push_back(dim);
			}

			if (grid.m_nGridDim.size() == 1) {
				vecDimOut.push_back(dim0);
			} else if (grid.m_nGridDim.size() == 2) {
				vecDimOut.push_back(dim0);
				vecDimOut.push_back(dim1);
			}

			vecNcVarOut[v] = ncout.add_var(
				vecVariableNamesOut[v].c_str(),
				ncFloat,
				vecDimOut.size(),
				const_cast<const NcDim**>(&(vecDimOut[0])));

			if (vecNcVarOut[v] == NULL) {
				_EXCEPTION2("Error creating variable \"%s\" in file \"%s\"",
				vecVariableNamesOut[v].c_str(),
				vecOutputFileList[f].c_str());
			}

			// Add _FillValue
			if (strFillValue != "") {
				vecNcVarOut[v]->add_att("_FillValue", dFillValue);
			}

			// Copy attributes
			Variable & var = varreg.Get(vecVarIxIn[v]);
			if (!var.IsOp()) {
				NcVar * ncvarIn = NULL;
				size_t sFileIx = vecFiles.FindContainingVariable(var.GetName(), &ncvarIn);
				if (ncvarIn != NULL) {
					std::vector<std::string> vecDoNotCopyNames;
					if (strFillValue != "") {
						vecDoNotCopyNames.push_back("_FillValue");
						vecDoNotCopyNames.push_back("missing_value");
					}
					CopyNcVarAttributes(ncvarIn, vecNcVarOut[v], vecDoNotCopyNames);
				}
				vecNcVarOutWroteUnits[v] = true;
			} else {
				vecNcVarOut[v]->add_att("_FillValue", DataOp::DefaultFillValue);
			}
		}

		// Loop through all times in this file
		for (int t = 0; t < vecOutputTimes.size(); t++) {
			if (fFileHasTime) {
				AnnounceStartBlock("Time %s (%lu/%lu)",
					vecOutputTimes[t].ToString().c_str(),
					t+1, vecOutputTimes.size());
			}

			vecFiles.SetTime(vecOutputTimes[t]);

			// Loop through all variables
			varreg.ResetProcessingQueue();
			while(varreg.AdvanceProcessingQueue()) {
				size_t v = varreg.GetProcessingQueueVarPos();

				_ASSERT(v < vecVariableNamesIn.size());
				_ASSERT(v < vecVariableNamesOut.size());

				// Load the data for the search variable
				Variable & var = varreg.GetProcessingQueueVariable();
				Announce("Processing \"%s\" -> \"%s\"",
					var.ToString(varreg).c_str(),
					vecVariableNamesOut[v].c_str());

				var.LoadGridData(varreg, vecFiles, grid);
				const DataArray1D<float> & dataState = var.GetData();

				// Write units attribute, if present
				if (!vecNcVarOutWroteUnits[v]) {
					vecNcVarOutWroteUnits[v] = true;
					if (var.GetUnits() != "") {
						vecNcVarOut[v]->add_att("units", var.GetUnits().c_str());
					}
				}

				// Get the output position and size
				VariableAuxIndex lPos = varreg.GetProcessingQueueAuxIx();
				if (fFileHasTime) {
					lPos.insert(lPos.begin(), t);
				}
				VariableAuxIndex lSize(lPos.size(), 1);

				if (grid.m_nGridDim.size() == 1) {
					_ASSERT(dataState.GetRows() == dim0->size());
					lPos.push_back(0);
					lSize.push_back(dim0->size());

				} else {
					_ASSERT(dataState.GetRows() == dim0->size() * dim1->size());
					lPos.push_back(0);
					lPos.push_back(0);
					lSize.push_back(dim0->size());
					lSize.push_back(dim1->size());
				}

				// Write output
				vecNcVarOut[v]->set_cur(&(lPos[0]));
				vecNcVarOut[v]->put(&(dataState[0]), &(lSize[0]));

				NcError err;
				if (err.get_err() != NC_NOERR) {
					_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
				}
			}

			if (fFileHasTime) {
				AnnounceEndBlock(NULL);
			}
		}
		if (vecInputFileList.size() != 1) {
			AnnounceEndBlock("Done");
		}

#if defined(TEMPEST_MPIOMP)
		AnnounceSetOutputBuffer(stdout);
		if (fpLog != NULL) {
			fclose(fpLog);
		}
		AnnounceOnlyOutputOnRankZero();
#endif
	}

	AnnounceEndBlock(NULL);

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}

