///////////////////////////////////////////////////////////////////////////////
///
///	\file    HeatWaveMetrics.cpp
///	\author  Paul Ullrich
///	\version December 12, 2022
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

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "Variable.h"

///////////////////////////////////////////////////////////////////////////////

class HeatWaveMetricsParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HeatWaveMetricsParam() :
		iVerbosityLevel(0),
		fRegional(false),
		fSplitYears(false),
		strTagVar("binary_tag"),
		strTempVar("t2m"),
		strLongitudeName("lon"),
		strLatitudeName("lat")
	{ }

public:
	// Verbosity level
	int iVerbosityLevel;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Split years
	bool fSplitYears;

	// Variable name for tags
	std::string strTagVar;

	// Variable name for temperature
	std::string strTempVar;

	// Name of longitude variabe
	std::string strLongitudeName;

	// Name of latitude variable
	std::string strLatitudeName;
};

///////////////////////////////////////////////////////////////////////////////

void HeatWaveMetrics(
	const std::vector<std::string> & vecInputFiles,
	const std::string & strOutputFile,
	const std::string & strConnectivity,
	VariableRegistry & varreg,
	const HeatWaveMetricsParam & param
) {
	_ASSERT(vecInputFiles.size() > 0);

	// Unload data from the VariableRegistry
	varreg.UnloadAllGridData();

	// Parse variable
	VariableIndex varixTag = varreg.FindOrRegister(param.strTagVar);
	if (varixTag == InvalidVariableIndex) {
		_EXCEPTION1("Unable to register variable \"%s\"", param.strTagVar.c_str());
	}

	VariableIndex varixTemp = varreg.FindOrRegister(param.strTempVar);
	if (varixTemp == InvalidVariableIndex) {
		_EXCEPTION1("Unable to register variable \"%s\"", param.strTempVar.c_str());
	}

	// Define the SimpleGrid
	SimpleGrid grid;

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
				param.strLatitudeName,
				param.strLongitudeName,
				param.fRegional,
				false);

			if (grid.m_nGridDim.size() != 2) {
				_EXCEPTIONT("Logic error when generating connectivity");
			}
			AnnounceEndBlock("Done");
		}
	}

	// Output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Copy over latitude, longitude to output file
	{
		NcFile ncInput(vecInputFiles[0].c_str());

		CopyNcVarIfExists(ncInput, ncOutput, param.strLatitudeName, true);
		CopyNcVarIfExists(ncInput, ncOutput, param.strLongitudeName, true);
	}

	// Create count array
	DataArray1D<double> dCount(grid.GetSize());

	// Create duration array
	DataArray1D<double> dDuration(grid.GetSize());

	// Create accumulated duration array
	DataArray1D<double> dAccumDuration(grid.GetSize());

	// Create intensity array
	DataArray1D<double> dIntensity(grid.GetSize());

	// Create accumulated intensity array
	DataArray1D<double> dAccumIntensity(grid.GetSize());

	// Current year
	int nCurrentYear = (-1);

	// Loop through all files
	AnnounceStartBlock("Performing HeatWaveMetrics");
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

		// Read the time data
		const NcTimeDimension & vecTimes = vecFiles.GetNcTimeDimension(0);

		// Stop accumulating statistics if there is a change of year
		if ((nCurrentYear != vecTimes[0].GetYear()) && (param.fSplitYears)) {
			if (nCurrentYear == (-1)) {
				nCurrentYear = vecTimes[0].GetYear();
			} else {
				for (size_t i = 0; i < dCount.GetRows(); i++) {
					if (dAccumDuration[i] != 0) {
						dCount[i] += 1.0;

						dDuration[i] += dAccumDuration[i];
						dIntensity[i] += dAccumIntensity[i];

						dAccumDuration[i] = 0.0;
						dAccumIntensity[i] = 0.0;
					}
				}
			}
		}

		// Loop through all times
		for (int t = 0; t < vecTimes.size(); t++) {

			// Load the search variable data
			Variable & varTag = varreg.Get(varixTag);
			vecFiles.SetTime(vecTimes[t]);
			varTag.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataTag = varTag.GetData();

			Variable & varTemp = varreg.Get(varixTemp);
			vecFiles.SetTime(vecTimes[t]);
			varTemp.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataTemp = varTemp.GetData();

			// Perform actual metrics calculation
			for (size_t i = 0; i < dCount.GetRows(); i++) {
				if (dataTag[i] == 0) {
					if (dAccumDuration[i] != 0) {
						dCount[i] += 1.0;

						dDuration[i] += dAccumDuration[i];
						dIntensity[i] += dAccumIntensity[i];

						dAccumDuration[i] = 0.0;
						dAccumIntensity[i] = 0.0;
					}
				} else {
					dAccumDuration[i] += 1.0;
					dAccumIntensity[i] += dataTemp[i];
				}
			}
		}
	}

	// Final step
	for (size_t i = 0; i < dCount.GetRows(); i++) {
		if (dAccumDuration[i] != 0) {
			dCount[i] += 1.0;

			dDuration[i] += dAccumDuration[i];
			dIntensity[i] += dAccumIntensity[i];
		}

		dDuration[i] /= dCount[i];
		dIntensity[i] /= dCount[i];
	}

	AnnounceEndBlock("Done");

	// Write results
	{
		AnnounceStartBlock("Writing results");

		NcDim * dim0 = NULL;
		NcDim * dim1 = NULL;

		// Unstructured grids
		if (grid.m_nGridDim.size() == 1) {
			dim0 = ncOutput.get_dim("ncol");
			if (dim0 == NULL) {
				dim0 = ncOutput.add_dim("ncol", grid.GetSize());
			}
			if (dim0 == NULL) {
				_EXCEPTION1("Error creating dim \"ncol\" in file \"%s\"",
					strOutputFile.c_str());
			}
			if (dim0->size() != grid.GetSize()) {
				_EXCEPTION4("Dimension \"%s\" in file \"%s\" has inconsistent length (%i vs %i)",
					dim0->name(),
					strOutputFile.c_str(),
					dim0->size(),
					grid.GetSize());
			}

			NcVar * varCount = ncOutput.add_var("count", ncDouble, dim0);
			if (varCount == NULL) {
				_EXCEPTION1("Error creating variable \"count\" in file \"%s\"",
					strOutputFile.c_str());
			}
			varCount->put(&(dCount[0]), dim0->size());

			NcVar * varAvgDuration = ncOutput.add_var("avg_duration", ncDouble, dim0);
			if (varAvgDuration == NULL) {
				_EXCEPTION1("Error creating variable \"avg_duration\" in file \"%s\"",
					strOutputFile.c_str());
			}
			varAvgDuration->put(&(dDuration[0]), dim0->size());

			NcVar * varAvgIntensity = ncOutput.add_var("avg_intensity", ncDouble, dim0);
			if (varAvgIntensity == NULL) {
				_EXCEPTION1("Error creating variable \"avg_intensity\" in file \"%s\"",
					strOutputFile.c_str());
			}
			varAvgIntensity->put(&(dIntensity[0]), dim0->size());

		// Structured grids
		} else {

			// Get output dimensions
			NcVar * varLat = ncOutput.get_var(param.strLatitudeName.c_str());
			if ((varLat != NULL) && (varLat->num_dims() == 2)) {
				dim0 = varLat->get_dim(0);
				dim1 = varLat->get_dim(1);

			} else {
				dim0 = ncOutput.get_dim(param.strLatitudeName.c_str());
				if (dim0 == NULL) {
					_EXCEPTION1("Unable to copy variable \"%s\" to output file",
						param.strLatitudeName.c_str());
				}
				dim1 = ncOutput.get_dim(param.strLongitudeName.c_str());
				if (dim1 == NULL) {
					_EXCEPTION1("Unable to copy variable \"%s\" to output file",
						param.strLongitudeName.c_str());
				}
			}

			NcVar * varCount = ncOutput.add_var("count", ncDouble, dim0, dim1);
			if (varCount == NULL) {
				_EXCEPTION1("Error creating variable \"count\" in file \"%s\"",
					strOutputFile.c_str());
			}
			varCount->put(&(dCount[0]), dim0->size(), dim1->size());

			NcVar * varAvgDuration = ncOutput.add_var("avg_duration", ncDouble, dim0, dim1);
			if (varAvgDuration == NULL) {
				_EXCEPTION1("Error creating variable \"avg_duration\" in file \"%s\"",
					strOutputFile.c_str());
			}
			varAvgDuration->put(&(dDuration[0]), dim0->size(), dim1->size());

			NcVar * varAvgIntensity = ncOutput.add_var("avg_intensity", ncDouble, dim0, dim1);
			if (varAvgIntensity == NULL) {
				_EXCEPTION1("Error creating variable \"avg_intensity\" in file \"%s\"",
					strOutputFile.c_str());
			}
			varAvgIntensity->put(&(dIntensity[0]), dim0->size(), dim1->size());

		}

		AnnounceEndBlock("Done");
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	NcError error(NcError::silent_nonfatal);

try {
	// Parameters for PersistentBlobs
	HeatWaveMetricsParam hwparam;

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

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineBool(hwparam.fRegional, "regional");
		CommandLineBool(hwparam.fSplitYears, "split_years");
		CommandLineString(hwparam.strTagVar, "tagvar", "binary_tag");
		CommandLineString(hwparam.strTempVar, "t2var", "t2");
		CommandLineString(hwparam.strLongitudeName, "lonname", "lon");
		CommandLineString(hwparam.strLatitudeName, "latname", "lat");
		CommandLineInt(hwparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(hwparam.iVerbosityLevel);

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

	// Perform HeatWaveMetrics
	HeatWaveMetrics(
		vecInputFiles,
		strOutputFile,
		strConnectivity,
		varreg,
		hwparam);

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

}

///////////////////////////////////////////////////////////////////////////////

