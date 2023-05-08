///////////////////////////////////////////////////////////////////////////////
///
///	\file    PersistentBlobs.cpp
///	\author  Paul Ullrich
///	\version December 9th, 2022
///
///	<remarks>
///		Copyright 2022 Paul Ullrich
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
#include "SimpleGrid.h"
#include "BlobUtilities.h"
#include "CoordTransforms.h"
#include "Units.h"
#include "TimeMatch.h"

#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include "ThresholdOp.h"

#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a persistence thresholding operator.
///	</summary>
class PersistenceThresholdOp {

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
		NotEqualTo,
		NoThreshold
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PersistenceThresholdOp() :
		m_varix(InvalidVariableIndex),
		m_eOp(GreaterThan),
		m_dValue(0.0),
		m_nTimesteps(0)
	{ }

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Timesteps,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
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
					eReadMode = ReadMode_Timesteps;

				// Read in distance to point that satisfies threshold
				} else if (eReadMode == ReadMode_Timesteps) {
					m_nTimesteps = atoi(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
							"\nRequired: \"<name>,<operation>"
							",<value>,<timesteps>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
					"\nRequired: \"<name>,<operation>,<value>,<timesteps>\"",
					strOp.c_str());
		}

		if (m_nTimesteps < 1) {
			_EXCEPTIONT("For persistence threshold op, timesteps must be positive");
		}

		// Output announcement
		std::string strDescription = varreg.GetVariableString(m_varix);
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

		char szBuffer[128];
		sprintf(szBuffer, "%f for %i timesteps",
			m_dValue, m_nTimesteps);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
	}

	///	<summary>
	///		Return true if the given value satisfied the threshold.
	///	</summary>
	inline bool IsSatisfiedBy(
		double dGivenValue
	) {
		if (m_eOp == GreaterThan) {
			if (dGivenValue > m_dValue) {
				return true;
			}

		} else if (m_eOp == LessThan) {
			if (dGivenValue < m_dValue) {
				return true;
			}

		} else if (m_eOp == GreaterThanEqualTo) {
			if (dGivenValue >= m_dValue) {
				return true;
			}

		} else if (m_eOp == LessThanEqualTo) {
			if (dGivenValue <= m_dValue) {
				return true;
			}

		} else if (m_eOp == EqualTo) {
			if (dGivenValue == m_dValue) {
				return true;
			}

		} else if (m_eOp == NotEqualTo) {
			if (dGivenValue != m_dValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		return false;
	}

public:
	///	<summary>
	///		Variable to use for thresholding.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Number of timesteps where criteria must hold.
	///	</summary>
	int m_nTimesteps;
};

///////////////////////////////////////////////////////////////////////////////

class PersistentBlobsParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PersistentBlobsParam() :
		iVerbosityLevel(0),
		fOutFloat(false),
		fRegional(false),
		strTagVar("binary_tag"),
		strStartTime(""),
		strEndTime(""),
		strTimeFilter(""),
		strLongitudeName("lon"),
		strLatitudeName("lat"),
		pvecPersistenceThresholdOp(NULL)
	{ }

public:
	// Verbosity level
	int iVerbosityLevel;

	// Write output as float
	bool fOutFloat;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Variable name for tags
	std::string strTagVar;

	// Start time
	std::string strStartTime;

	// End time
	std::string strEndTime;

	// Time filter
	std::string strTimeFilter;

	// Name of longitude variabe
	std::string strLongitudeName;

	// Name of latitude variable
	std::string strLatitudeName;

	// Vector of threshold operators
	std::vector<PersistenceThresholdOp> * pvecPersistenceThresholdOp;
};

///////////////////////////////////////////////////////////////////////////////

void PersistentBlobs(
	const std::vector<std::string> & vecInputFiles,
	const std::vector<std::string> & vecOutputFiles,
	const std::string & strConnectivity,
	VariableRegistry & varreg,
	const PersistentBlobsParam & param
) {
	_ASSERT(vecInputFiles.size() > 0);
	_ASSERT(vecInputFiles.size() == vecOutputFiles.size());

	// Dereference pointers to operators
	_ASSERT(param.pvecPersistenceThresholdOp != NULL);
	std::vector<PersistenceThresholdOp> & vecPersistenceThresholdOp =
		*(param.pvecPersistenceThresholdOp);

	_ASSERT(vecPersistenceThresholdOp.size() == 1);

#ifdef TEMPEST_NOREGEX
	if (param.strTimeFilter != "") {
		_EXCEPTIONT("Cannot use --timefilter with -DTEMPEST_NOREGEX compiler flag");
	}
#endif
#ifndef TEMPEST_NOREGEX
	// Parse --timefilter
	std::regex reTimeSubset;
	if (param.strTimeFilter != "") {
		
		// Test regex support
		TestRegex();

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

	// Number of timesteps in threshold op
	int nPersistenceThresholdTimesteps = vecPersistenceThresholdOp[0].m_nTimesteps;

	// Unload data from the VariableRegistry
	varreg.UnloadAllGridData();

	// Define the SimpleGrid
	SimpleGrid grid;

	std::string strLatitudeName(param.strLatitudeName);
	std::string strLongitudeName(param.strLongitudeName);

	{
		// Load in the benchmark file
		NcFileVector vecFiles;
		vecFiles.ParseFromString(vecInputFiles[0]);

		// Check for connectivity file
		if (strConnectivity != "") {
			AnnounceStartBlock("Generating grid information from connectivity file");
			grid.FromFile(strConnectivity);
			AnnounceEndBlock("Done");

		// No connectivity file; check for latitude/longitude dimension
		} else {
			AnnounceStartBlock("No connectivity file specified");
			Announce("Attempting to generate latitude-longitude grid from data file");

			grid.GenerateLatitudeLongitude(
				vecFiles[0],
				strLatitudeName,
				strLongitudeName,
				param.fRegional,
				false);

			if (grid.m_nGridDim.size() != 2) {
				_EXCEPTIONT("Logic error when generating connectivity");
			}
			AnnounceEndBlock("Done");
		}
	}

	// Current tiemstep index for writing tags
	int iOutputTimestep = - nPersistenceThresholdTimesteps + 1;

	// Accumulated count
	DataArray1D<int> dCount(grid.GetSize());

	// Countdown
	DataArray1D<int> dCountdown(grid.GetSize());

	// Tag variable
	DataArray1D<int> bTag(grid.GetSize());

	// Array of pointers to output files
	std::vector<NcFile *> vecpncOutputFiles(vecOutputFiles.size());

	// Current output file pointer
	int iCurrentOutputFile = 0;

	// Loop through all files
	AnnounceStartBlock("Performing PersistentBlobs");
	for (size_t f = 0; f < vecInputFiles.size(); f++) {

		// Open the input file(s)
		NcFileVector vecFiles;
		vecFiles.ParseFromString(vecInputFiles[f]);

		// Get time dimension
		NcVar * varTime = vecFiles[0]->get_var("time");
		if (varTime == NULL) {
			_EXCEPTION1("File \"%s\" missing \"time\" variable",
				vecFiles.GetFilename(0).c_str());
		}
		NcDim * dimTime = vecFiles[0]->get_dim("time");
		if (dimTime == NULL) {
			if (varTime->num_dims() != 0) {
				_EXCEPTION1("File \"%s\" missing \"time\" dimension",
					vecFiles.GetFilename(0).c_str());
			}
		}

		// Read the time data and filter
		const NcTimeDimension & vecTimes = vecFiles.GetNcTimeDimension(0);

		std::vector<bool> vecTimeRetained(vecTimes.size(), false);
		std::vector<Time> vecOutputTimes;

		Time timeStartTime(Time::CalendarUnknown);
		Time timeEndTime(Time::CalendarUnknown);

		if (param.strStartTime != "") {
			timeStartTime = Time(vecTimes[0].GetCalendarType());
			timeStartTime.FromFormattedString(param.strStartTime);
		}
		if (param.strEndTime != "") {
			timeEndTime = Time(vecTimes[0].GetCalendarType());
			timeEndTime.FromFormattedString(param.strEndTime);
		}

		for (int t = 0; t < vecTimes.size(); t++) {
			if (timeStartTime.GetCalendarType() != Time::CalendarUnknown) {
				double dDeltaSeconds = timeStartTime - vecTimes[t];
				if (dDeltaSeconds > 0.0) {
					continue;
				}
			}
			if (timeEndTime.GetCalendarType() != Time::CalendarUnknown) {
				double dDeltaSeconds = vecTimes[t] - timeEndTime;
				if (dDeltaSeconds > 0.0) {
					continue;
				}
			}

#ifndef TEMPEST_NOREGEX
			if (param.strTimeFilter != "") {
				std::string strTime = vecTimes[t].ToString();
				std::smatch match;
				if (!std::regex_search(strTime, match, reTimeSubset)) {
					continue;
				}
			}
#endif
			vecOutputTimes.push_back(vecTimes[t]);
			vecTimeRetained[t] = true;	
		}

		// Create reference to NetCDF input file
		NcFile & ncInput = *(vecFiles[0]);

		// Open the NetCDF output file
		vecpncOutputFiles[f] = new NcFile(vecOutputFiles[f].c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
		if (!vecpncOutputFiles[f]->is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
				vecOutputFiles[f].c_str());
		}

		// Copy over time variable to output file
		NcDim * dimTimeOut = NULL;
		if ((dimTime != NULL) && (varTime != NULL)) {
			if (vecOutputTimes.size() != vecTimes.size()) {
				CopyNcVarTimeSubset(ncInput, *(vecpncOutputFiles[f]), "time", vecOutputTimes);

			} else {
				CopyNcVar(ncInput, *(vecpncOutputFiles[f]), "time", true);
			}

			dimTimeOut = vecpncOutputFiles[f]->get_dim("time");
			if (dimTimeOut == NULL) {
				_EXCEPTIONT("Error copying variable \"time\" to output file");
			}
		}

		// Create output variable
		NcDim * dim0 = NULL;
		NcDim * dim1 = NULL;
		NcVar * varTag = NULL;

		PrepareBlobOutputVar(
			ncInput,
			*(vecpncOutputFiles[f]),
			vecOutputFiles[f],
			grid,
			param.strTagVar,
			strLatitudeName,
			strLongitudeName,
			(param.fOutFloat)?(ncFloat):(ncByte),
			dimTimeOut,
			&dim0,
			&dim1,
			&varTag);

		_ASSERT(varTag != NULL);

		// Loop through all times
		int to = (-1);
		for (int t = 0; t < vecTimes.size(); t++) {

			std::string strTime = vecTimes[t].ToString();

			// Skip if not retained
			if (!vecTimeRetained[t]) {
				Announce("Input time %s.. Skipping (timefilter)", strTime.c_str());
				continue;
			}

			Announce("Input time %s", strTime.c_str());

			to++;

			// Load the search variable data
			Variable & var = varreg.Get(vecPersistenceThresholdOp[0].m_varix);
			vecFiles.SetTime(vecTimes[t]);
			var.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataState = var.GetData();

			float dFillValue = var.GetFillValueFloat();

			// Copy over count and update
			for (size_t i = 0; i < dCount.GetRows(); i++) {
				if ((dataState[i] == dFillValue) || (dataState[i] != dataState[i])) {
					dCount[i] = 0;
					continue;
				}
				if (!vecPersistenceThresholdOp[0].IsSatisfiedBy(dataState[i])) {
					dCount[i] = 0;
					continue;
				}

				dCount[i]++;
				if (dCount[i] >= nPersistenceThresholdTimesteps) {
					dCountdown[i] = nPersistenceThresholdTimesteps;
				}
			}

			// Perform output
			if (iOutputTimestep >= 0) {
				AnnounceStartBlock("Build tagged cell array at time %i in file %i",
					iOutputTimestep, iCurrentOutputFile);

				// Build tag array
				bTag.Zero();
				for (size_t i = 0; i < dCount.GetRows(); i++) {
					if (dCountdown[i] > 0) {
						bTag[i] = 1;
						dCountdown[i]--;
					}
				}

				// Write tag array
				_ASSERT(iCurrentOutputFile < vecpncOutputFiles.size());
				NcVar * varTag = vecpncOutputFiles[iCurrentOutputFile]->get_var(param.strTagVar.c_str());
				_ASSERT(varTag != NULL);

				if (grid.m_nGridDim.size() == 1) {
					varTag->set_cur(iOutputTimestep, 0);
					varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0]);

				} else if (grid.m_nGridDim.size() == 2) {
					varTag->set_cur(iOutputTimestep, 0, 0);
					varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);

				} else {
					_EXCEPTION();
				}

				// Advance to next timestep or file
				iOutputTimestep++;
				if (iOutputTimestep == varTag->get_dim(0)->size()) {
					iCurrentOutputFile++;
					iOutputTimestep = 0;
				}

				AnnounceEndBlock("Done");

			} else {
				iOutputTimestep++;
			}
		}

		// Finish off
		if (f == vecInputFiles.size()-1) {
			for (int t = 0; t < nPersistenceThresholdTimesteps - 1; t++) {
				AnnounceStartBlock("Build tagged cell array at time %i in file %i",
					iOutputTimestep, iCurrentOutputFile);

				// Build tag array
				bTag.Zero();
				for (size_t i = 0; i < dCount.GetRows(); i++) {
					if (dCountdown[i] > 0) {
						bTag[i] = 1;
						dCountdown[i]--;
					}
				}

				// Write tag array
				_ASSERT(iCurrentOutputFile < vecpncOutputFiles.size());
				NcVar * varTag = vecpncOutputFiles[iCurrentOutputFile]->get_var(param.strTagVar.c_str());
				_ASSERT(varTag != NULL);

				if (grid.m_nGridDim.size() == 1) {
					varTag->set_cur(iOutputTimestep, 0);
					varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0]);

				} else if (grid.m_nGridDim.size() == 2) {
					varTag->set_cur(iOutputTimestep, 0, 0);
					varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);

				} else {
					_EXCEPTION();
				}

				// Advance to next timestep or file
				iOutputTimestep++;
				if (iOutputTimestep == varTag->get_dim(0)->size()) {
					iCurrentOutputFile++;
					iOutputTimestep = 0;
				}

				AnnounceEndBlock("Done");
			}
		}
	}
	AnnounceEndBlock("Done");

	// Close all output files
	AnnounceStartBlock("Cleanup");
	for (int f = 0; f < vecpncOutputFiles.size(); f++) {
		delete vecpncOutputFiles[f];
	}
	AnnounceEndBlock("Done");

}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
/*
#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif
*/
	NcError error(NcError::silent_nonfatal);

try {
	// Parameters for PersistentBlobs
	PersistentBlobsParam pbparam;

	// Input dat file
	std::string strInputFile;

	// Input list file
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutputFile;

	// Output file list
	std::string strOutputFileList;

	// Threshold commands
	std::string strThresholdCmd;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputFileList, "out_list", "");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,timesteps;...]");
		CommandLineBool(pbparam.fRegional, "regional");
		CommandLineString(pbparam.strStartTime, "time_start", "");
		CommandLineString(pbparam.strEndTime, "time_end", "");
		CommandLineString(pbparam.strTimeFilter, "timefilter", "");
		CommandLineBool(pbparam.fOutFloat, "out_float");
		CommandLineString(pbparam.strTagVar, "tagvar", "binary_tag");
		CommandLineString(pbparam.strLongitudeName, "lonname", "lon");
		CommandLineString(pbparam.strLatitudeName, "latname", "lat");
		CommandLineInt(pbparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(pbparam.iVerbosityLevel);

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
	if ((strOutputFile.length() != 0) && (strOutputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out) or (--out_list)"
			" may be specified");
	}

	// Load input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		std::ifstream ifInputFileList(strInputFileList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputFiles.push_back(strFileLine);
		}
	}

	// Load output file list
	std::vector<std::string> vecOutputFiles;

	if (strOutputFileList.length() == 0) {
		vecOutputFiles.push_back(strOutputFile);

	} else {
		std::ifstream ifOutputFileList(strOutputFileList.c_str());
		if (!ifOutputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strOutputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecOutputFiles.push_back(strFileLine);
		}

		if (vecOutputFiles.size() != vecInputFiles.size()) {
			_EXCEPTIONT("File --in_file_list must match --out_file_list");
		}
	}

	// Ensure number of input files equals number of output files
	if (vecInputFiles.size() != vecOutputFiles.size()) {
		_EXCEPTION2("Number of input files (%lu) must match number of output files (%lu)",
			vecInputFiles.size(), vecOutputFiles.size());
	}

	// Parse the persistence threshold operator command string
	std::vector<PersistenceThresholdOp> vecPersistenceThresholdOp;
	pbparam.pvecPersistenceThresholdOp = &vecPersistenceThresholdOp;

	if (strThresholdCmd != "") {
		AnnounceStartBlock("Parsing persistence threshold operations");

		int iLast = 0;
		for (int i = 0; i <= strThresholdCmd.length(); i++) {

			if ((i == strThresholdCmd.length()) ||
				(strThresholdCmd[i] == ';') ||
				(strThresholdCmd[i] == ':')
			) {
				std::string strSubStr =
					strThresholdCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecPersistenceThresholdOp.size());
				vecPersistenceThresholdOp.resize(iNextOp + 1);
				vecPersistenceThresholdOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Currently only one threshold command supported
	if (vecPersistenceThresholdOp.size() != 1) {
		_EXCEPTIONT("At present exactly one --thresholdcmd must be specified");
	}

	// Perform PersistentBlobs
	PersistentBlobs(
		vecInputFiles,
		vecOutputFiles,
		strConnectivity,
		varreg,
		pbparam);

	AnnounceEndBlock("Done");

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

///////////////////////////////////////////////////////////////////////////////

