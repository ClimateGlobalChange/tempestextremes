///////////////////////////////////////////////////////////////////////////////
///
///	\file    Climatology.cpp
///	\author  Paul Ullrich
///	\version June 3, 2020
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
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "FilenameList.h"
#include "Variable.h"
#include "TimeObj.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"

#include <complex>
#include <vector>
#include <set>
#include <map>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Apply a Fourier filter to a given data sequence.
///	</summary>

// TODO: Test this function
void fourier_filter(
	double * const data,
	size_t sCount,
	size_t sStride,
	size_t sModes,
	DataArray1D<double> & an,
	DataArray1D<double> & bn
) {
	_ASSERT(an.GetRows() >= sModes);
	_ASSERT(bn.GetRows() >= sModes);

	an.Zero();
	bn.Zero();

	{
		for (int n = 0; n < sCount; n++) {
			an[0] += data[n*sStride];
		}

		double dProd0 = 2.0 * M_PI / static_cast<double>(sCount);
		for (int k = 1; k < sModes; k++) {
			double dProd1 = dProd0 * static_cast<double>(k);
			for (int n = 0; n < sCount; n++) {
				double dProd2 = dProd1 * static_cast<double>(n);
				an[k] += data[n*sStride] * cos(dProd2);
				bn[k] -= data[n*sStride] * sin(dProd2);
			}

			//an[sCount-k] = an[k];
			//bn[sCount-k] = -bn[k];
		}
	}
	{
		for (int n = 0; n < sCount; n++) {
			data[n*sStride] = 0.0;
		}

		double dProd0 = 2.0 * M_PI / static_cast<double>(sCount);
		for (int n = 0; n < sCount; n++) {
			data[n*sStride] += an[0];

			double dProd1 = dProd0 * static_cast<double>(n);
			for (int k = 1; k < sModes; k++) {
				double dProd2 = dProd1 * static_cast<double>(k);
				double dProd3 = dProd1 * static_cast<double>(sCount-k);
				data[n*sStride] += an[k] * cos(dProd2) - bn[k] * sin(dProd2);
				data[n*sStride] += an[k] * cos(dProd3) + bn[k] * sin(dProd3);

				//printf("%1.15e %1.15e : ", cos(dProd2), cos(dProd3));
				//printf("%1.15e %1.15e\n", sin(dProd2), sin(dProd3));
			}
			//for (int k = sCount-sModes+1; k < sCount; k++) {
			//	double dProd2 = dProd1 * static_cast<double>(k);
			//	data[n*sStride] += an[k] * cos(dProd2) - bn[k] * sin(dProd2);
			//}
		}

		for (int n = 0; n < sCount; n++) {
			data[n*sStride] /= static_cast<double>(sCount);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a number of bytes into a user-readable string.
///	</summary>
std::string MemoryBytesToString(
	size_t sTotalOutputStorage
) {
	std::string strStorageUnit;
	double dTotalOutputStorage;
	if (sTotalOutputStorage > 1024 * 1024 * 1024) {
		strStorageUnit = "G";
		dTotalOutputStorage =
			static_cast<double>(sTotalOutputStorage)
			/ static_cast<double>(1024 * 1024 * 1024);
	} else if (sTotalOutputStorage > 1024 * 1024) {
		strStorageUnit = "M";
		dTotalOutputStorage =
			static_cast<double>(sTotalOutputStorage)
			/ static_cast<double>(1024 * 1024);
	} else {
		strStorageUnit = "K";
		dTotalOutputStorage =
			static_cast<double>(sTotalOutputStorage)
			/ static_cast<double>(1024);
	}

	char szBuffer[20];
	sprintf(szBuffer, "%1.2f", dTotalOutputStorage);


	return (std::string(szBuffer) + strStorageUnit);
}

///////////////////////////////////////////////////////////////////////////////

void CalculatePartialDataPosSize(
	const DataArray1D<long> & lDimSize,
	long lSkipDims,
	long long llNcArrayPos,
	long long llNcArraySize,
	DataArray1D<long> & lPos,
	DataArray1D<long> & lSize
) {
	_ASSERT(lDimSize.GetRows() == lPos.GetRows());
	_ASSERT(lDimSize.GetRows() == lSize.GetRows());

	// llAccumulatedSize[d+1] stores the size of one read of an element
	// in that dimension
	DataArray1D<long long> llAccumulatedSize(lDimSize.GetRows()+1);
	llAccumulatedSize[lDimSize.GetRows()] = 1;
	for (long d = lDimSize.GetRows()-1; d >= lSkipDims; d--) {
		llAccumulatedSize[d] =
			llAccumulatedSize[d+1]
			* static_cast<long long>(lDimSize[d]);
	}

	// Given a 1D llNcArrayPos find the multi-dimensional index corresponding
	// to this offset in the larger array.
	long long llNcArrayPosTemp = llNcArrayPos;
	for (long d = lDimSize.GetRows()-1; d >= lSkipDims; d--) {
		lPos[d] = llNcArrayPosTemp % lDimSize[d];
		llNcArrayPosTemp /= lDimSize[d];
	}
	if (llNcArrayPosTemp > 0) {
		_EXCEPTIONT("NcVar 1D index out of range");
	}

	for (long d = lSkipDims; d < lDimSize.GetRows(); d++) {
		if ((llNcArraySize <= llAccumulatedSize[d]) &&
			(llNcArraySize > llAccumulatedSize[d+1])
		) {
			lSize[d] = llNcArraySize / llAccumulatedSize[d+1];
			for (d++; d < lDimSize.GetRows(); d++) {
				lSize[d] = lDimSize[d];
			}
			break;

		} else {
			lSize[d] = 1;
		}
	}

	for (long d = lDimSize.GetRows()-1; d >= lSkipDims; d--) {
		if (lSize[d] + lPos[d] > lDimSize[d]) {
			lSize[d] = lDimSize[d] - lPos[d];
			for (d--; d >= lSkipDims; d--) {
				lSize[d] = 1;
			}
			break;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

enum NetCDF_GetPut {
	NetCDF_Get,
	NetCDF_Put
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
long long NcGetPutSpecifiedDataSize(
	NcVar * ncvar,
	DataArray1D<long> lPos0,
	long lSizePerIteration,
	long lIteration,
	DataArray1D<T> & data,
	NetCDF_GetPut eGetPut
) {
	_ASSERT(lSizePerIteration <= data.GetRows());
	_ASSERT(ncvar->num_dims() > lPos0.GetRows());

	long long llNcArrayPos =
		static_cast<long long>(lSizePerIteration)
		* static_cast<long long>(lIteration);

	DataArray1D<long> lDimSize(ncvar->num_dims());
	for (long d = 0; d < lDimSize.GetRows(); d++) {
		lDimSize[d] = ncvar->get_dim(d)->size();
	}

	long long llTotalReadWriteSize = 0;

	// Position and size of read/write
	DataArray1D<long> lPos(lDimSize.GetRows());
	DataArray1D<long> lSize(lDimSize.GetRows());

	for (int d = 0; d < lPos0.GetRows(); d++) {
		lPos[d] = lPos0[d];
		lSize[d] = 1;
	}

	// Iterate until nothing left to be read/write
	while(lSizePerIteration > 0) {
		CalculatePartialDataPosSize(
			lDimSize,
			lPos0.GetRows(),
			llNcArrayPos,
			lSizePerIteration,
			lPos,
			lSize);

		// Number of elements to read
		long lReadWriteSize = 1;
		for (long d = lPos0.GetRows(); d < lDimSize.GetRows(); d++) {
			lReadWriteSize *= lSize[d];
		}
		_ASSERT(lReadWriteSize <= lSizePerIteration);

		// Update data remaining
		lSizePerIteration -= lReadWriteSize;
		llNcArrayPos += lReadWriteSize;

		// Read data
/*
		if (eGetPut == NetCDF_Get) {
			printf("get ");
		} else {
			printf("put ");
		}
		for (long d = 0; d < lPos.GetRows(); d++) {
			printf("%lu ", lPos[d]);
		}
		printf(": ");
		for (long d = 0; d < lSize.GetRows(); d++) {
			printf("%lu ", lSize[d]);
		}
		printf(" :: %lu read/written %lu remaining", lReadWriteSize, lSizePerIteration);
		printf("\n");
*/
		ncvar->set_cur(&(lPos[0]));
		if (eGetPut == NetCDF_Get) {
			ncvar->get(&(data[llTotalReadWriteSize]), &(lSize[0]));
		} else if (eGetPut == NetCDF_Put) {
			ncvar->put(&(data[llTotalReadWriteSize]), &(lSize[0]));
		} else {
			_EXCEPTIONT("Logic error");
		}

		llTotalReadWriteSize += lReadWriteSize;
	}

	return llTotalReadWriteSize;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

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

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input data file
	std::string strInputFile;

	// List of input data files
	std::string strInputFileList;

	// Output data file
	std::string strOutputFile;

	// Variable to use
	std::string strVariable;

	// Maximum memory allocation per thread
	std::string strMemoryMax;

	// Type of mean climatology
	std::string strMeanType;

	// Include leap days
	bool fIncludeLeapDays;

	// Number of Fourier modes to retain
	int nFourierModes;

	// Dataset has missing data
	bool fMissingData;

	// Verbose
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strOutputFile, "out_data", "");
		CommandLineString(strVariable, "var", "");
		CommandLineStringD(strMemoryMax, "memmax", "2G", "[#K,#M,#G]");	
		CommandLineStringD(strMeanType, "mean", "daily", "[daily|monthly|seasonal|annual]");
		CommandLineBool(fIncludeLeapDays, "include_leap_days");
		CommandLineInt(nFourierModes, "time_modes", 0);
		CommandLineBool(fMissingData, "missingdata");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Validate arguments
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output file (--out) specified");
	}
	if (strVariable.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}
	if (fIncludeLeapDays) {
		_EXCEPTIONT("--include_leap_days not implemented");
	}

	// Parse input file list
	FilenameList vecInputFileList;
	if (strInputFile.length() != 0) {
		vecInputFileList.push_back(strInputFile);
		if (strInputFile.find(';') != std::string::npos) {
			_EXCEPTIONT("Only one filename allowed in --in_data");
		}
	} else {
		vecInputFileList.FromFile(strInputFileList, false);
	}
	_ASSERT(vecInputFileList.size() != 0);

	// Parse variable list
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVariable, vecVariableStrings);

	// Memory maximum
	size_t sMemoryMax = (-1);
	if (strMemoryMax != "") {
		char cType = tolower(strMemoryMax[strMemoryMax.length()-1]);
		std::string strMemoryMaxValue = strMemoryMax.substr(0,strMemoryMax.length()-1);
		if (!STLStringHelper::IsInteger(strMemoryMaxValue)) {
			_EXCEPTIONT("Invalid --memmax value, expected integer [#K,#M,#G]");
		}
		std::stringstream sstream(strMemoryMaxValue);
		sstream >> sMemoryMax;
		if (cType == 'k') {
			sMemoryMax *= 1024;
		} else if (cType == 'm') {
			sMemoryMax *= 1024 * 1024;
		} else if (cType == 'g') {
			sMemoryMax *= 1024 * 1024 * 1024;
		} else {
			_EXCEPTIONT("Invalid --memmax unit, expected [#K,#M,#G]");
		}
	}
	if (sMemoryMax == 0) {
		_EXCEPTIONT("--memmax must be non-zero");
	}

	AnnounceBanner();

	// Time units
	std::string strTimeUnits = "days since 0001-01-01";

	// Open output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputFile.c_str());
	}

	// Determine dimensionality of each variable

	std::vector<DimInfoVector> vecVarDimInfo;
	vecVarDimInfo.resize(vecVariableStrings.size());

	std::vector<NcVar*> vecNcVarOut;
	vecNcVarOut.resize(vecVariableStrings.size());

	std::string strTimeCalendar;
	size_t sOutputTimes = 0;

	std::vector<size_t> vecOutputAuxSize;
	vecOutputAuxSize.resize(vecVariableStrings.size(), 0);

	std::vector<size_t> vecOutputAuxCount;
	vecOutputAuxCount.resize(vecVariableStrings.size(), 0);

	{
		NcFile ncinfile(vecInputFileList[0].c_str());
		if (!ncinfile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"",
				vecInputFileList[0].c_str());
		}

		// Get information from the time variable
		AnnounceStartBlock("Loading time information from input datafile");
		NcVar * varInTime = ncinfile.get_var("time");
		if (varInTime == NULL) {
			_EXCEPTION1("File \"%s\" does not contain variable \"time\"",
				vecInputFileList[0].c_str());
		}
		if (varInTime->num_dims() != 1) {
			_EXCEPTION1("Variable \"time\" in NetCDF file \"%s\" has more than 1 dimension",
				vecInputFileList[0].c_str());
		}

		{

			NcAtt * attCalendar = varInTime->get_att("calendar");
			if (attCalendar == NULL) {
				Announce("Dataset \"time\" does not specify \"calendar\" attribute.  Assuming \"standard\".");
				strTimeCalendar = "standard";
			} else {
				strTimeCalendar = attCalendar->as_string(0);
			}

			Time::CalendarType caltype =
				Time::CalendarTypeFromString(strTimeCalendar);
			if (caltype == Time::CalendarUnknown) {
				_EXCEPTION1("Unknown calendar name \"%s\"; cannot determine number of days per year", strTimeCalendar.c_str());

			} else if (caltype == Time::Calendar360Day) {
				Announce("Calendar \"%s\" contains 360 days per year (30 days per month)",
					strTimeCalendar.c_str());
				sOutputTimes = 360;

			} else {
				Announce("Calendar contains 365 days per year",
					strTimeCalendar.c_str());
				sOutputTimes = 365;
			}
		}
		AnnounceEndBlock("Done");

		// Get the auxiliary dimensions of each variable
		AnnounceStartBlock("Loading variable information from input datafile");
		std::set<DimInfo> setDimInfo;
		for (int v = 0; v < vecVariableStrings.size(); v++) {
			NcVar * var = ncinfile.get_var(vecVariableStrings[v].c_str());
			if (var == NULL) {
				_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
					vecInputFileList[0].c_str(),
					vecVariableStrings[v].c_str());
			}

			if (var->num_dims() < 2) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have at least 2 dimensions",
					vecInputFileList[0].c_str(),
					vecVariableStrings[v].c_str());
			}

			NcDim * dim0 = var->get_dim(0);
			_ASSERT(dim0 != NULL);

			std::string strDimTimeName = dim0->name();
			if (strDimTimeName != "time") {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have leftmost dimension \"time\"",
					vecInputFileList[0].c_str(),
					vecVariableStrings[v].c_str());
			}

			for (int d = 1; d < var->num_dims(); d++) {
				NcDim * dimVar = var->get_dim(d);
				DimInfo diminfo(dimVar->name(), dimVar->size());
				vecVarDimInfo[v].push_back(diminfo);
				setDimInfo.insert(diminfo);
			}

			// Output storage requirements
			size_t sTotalAuxDims = 1;
			for (int d = 0; d < vecVarDimInfo[v].size(); d++) {
				sTotalAuxDims *= vecVarDimInfo[v][d].size;
			}

			// Number of auxiliary values to store at each iteration
			size_t sClimatologyStorageSize = sOutputTimes * sTotalAuxDims * sizeof(double);
			vecOutputAuxCount[v] = ((sClimatologyStorageSize + 1) / sMemoryMax) + 1;
			if (vecOutputAuxCount[v] > sTotalAuxDims) {
				_EXCEPTIONT("Insufficient memory for operation, increase --memmax");
			}
			vecOutputAuxSize[v] = sTotalAuxDims / vecOutputAuxCount[v];

			// Output information to user
			std::string strPerTimesliceStorage =
				MemoryBytesToString(sTotalAuxDims * sizeof(double));

			std::string strTotalOutputStorage =
				MemoryBytesToString(sOutputTimes * sTotalAuxDims * sizeof(double));

			AnnounceStartBlock("Variable \"%s\" has %lu distinct time series",
				vecVariableStrings[v].c_str(),
				sTotalAuxDims);
			Announce("Memory requirement: per-timeslice %s / total %s",
				strPerTimesliceStorage.c_str(),
				strTotalOutputStorage.c_str());
			Announce("%lu time series will be loaded per iteration over %lu iteration(s)",
				vecOutputAuxSize[v],
				vecOutputAuxCount[v]);
			AnnounceEndBlock(NULL);
		}
		AnnounceEndBlock("Done");

		// Write time dimension to output file
		AnnounceStartBlock("Initializing output file");

		Announce("Initialize time dimension and variable");
		NcDim * dimOutTime = ncoutfile.add_dim("time", (long)sOutputTimes);
		if (dimOutTime == NULL) {
			_EXCEPTIONT("Unable to create output dimension \"time\" in output file");
		}

		NcVar * varOutTime = ncoutfile.add_var("time", varInTime->type(), dimOutTime);
		if (varOutTime == NULL) {
			_EXCEPTIONT("Unable to create output dimension \"time\" in output file");
		}

		varOutTime->add_att("calendar", strTimeCalendar.c_str());
		varOutTime->add_att("units", strTimeUnits.c_str());

		if (varOutTime->type() == ncInt) {
			DataArray1D<int> nTimes(sOutputTimes);
			for (size_t s = 0; s < sOutputTimes; s++) {
				nTimes[s] = static_cast<int>(s);
			}
			varOutTime->put(&(nTimes[0]), sOutputTimes);

		} else if (varOutTime->type() == ncDouble) {
			DataArray1D<double> dTimes(sOutputTimes);
			for (size_t s = 0; s < sOutputTimes; s++) {
				dTimes[s] = static_cast<int>(s);
			}
			varOutTime->put(&(dTimes[0]), sOutputTimes);

		// TODO: Add option for int64 time
		} else {
			_EXCEPTIONT("Invalid type of \"time\" variable, expected \"int\" or \"double\"");
		}

		// Copy dimension information to output file
		Announce("Copy auxiliary dimension information from input file");
		std::map<std::string, NcDim *> mapOutputNcDim;
		for (auto iter = setDimInfo.begin(); iter != setDimInfo.end(); iter++) {
			NcDim * dimOut = ncoutfile.add_dim(iter->name.c_str(), iter->size);
			if (dimOut == NULL) {
				_EXCEPTION1("Unable to create output dimension \"%s\" in output file",
					iter->name.c_str());
			}

			CopyNcVarIfExists(ncinfile, ncoutfile, iter->name);

			mapOutputNcDim.insert(
				std::map<std::string, NcDim *>::value_type(
					iter->name, dimOut));
		}

		// Create output variables
		Announce("Creating output variables");
		for (int v = 0; v < vecVariableStrings.size(); v++) {
			std::vector<NcDim *> vecVarNcDims;
			vecVarNcDims.push_back(dimOutTime);
			for (int d = 0; d < vecVarDimInfo[v].size(); d++) {
				auto itNcDim = mapOutputNcDim.find(vecVarDimInfo[v][d].name);
				_ASSERT(itNcDim != mapOutputNcDim.end());
				vecVarNcDims.push_back(itNcDim->second);
			}
			NcVar * varOut =
				ncoutfile.add_var(
					vecVariableStrings[v].c_str(),
					ncFloat,
					vecVarNcDims.size(),
					const_cast<const NcDim**>(&(vecVarNcDims[0])));

			if (varOut == NULL) {
				_EXCEPTION1("Unable to create output variable \"%s\" in output file",
					vecVariableStrings[v].c_str());
			}

			NcVar * varIn = ncinfile.get_var(vecVariableStrings[v].c_str());
			_ASSERT(varIn != NULL);
			CopyNcVarAttributes(varIn, varOut);

			varOut->add_att("climatology", "daily mean");

			vecNcVarOut[v] = varOut;
		}

		AnnounceEndBlock("Done");
	}

	// Storage array for the time dimension
	DataArray1D<long> lPos0(1);

	// Loop through each variable
	for (int v = 0; v < vecVariableStrings.size(); v++) {
		AnnounceStartBlock("Calculating climatology of variable \"%s\"",
			vecVariableStrings[v].c_str());

		DataArray1D<float> dDataIn(vecOutputAuxSize[v]);

		// Loop through each iteration of this variable
		for (int c = 0; c < vecOutputAuxCount[v]; c++) {

			if (vecOutputAuxCount[v] > 1) {
				AnnounceStartBlock("Iteration %i/%i", c, vecOutputAuxCount[v]);
			}

			// Number of time slices
			DataArray1D<int> nTimeSlices(sOutputTimes);

			// Accumulated data array for mean
			DataArray2D<double> dAccumulatedData(sOutputTimes, vecOutputAuxSize[v]);

			// If missing data is present need to count number of data points
			// at each location.
			DataArray2D<int> nTimeSlicesGrid;
			if (fMissingData) {
				nTimeSlicesGrid.Allocate(sOutputTimes, vecOutputAuxSize[v]);
			}

			// Loop through each input file
			for (int f = 0; f < vecInputFileList.size(); f++) {

				// Open input file
				NcFile ncinfile(vecInputFileList[f].c_str());
				if (!ncinfile.is_valid()) {
					_EXCEPTION1("Unable to open NetCDF file \"%s\"",
						vecInputFileList[f].c_str());
				}

				// Get time information
				NcVar * varInTime = ncinfile.get_var("time");
				if (varInTime == NULL) {
					_EXCEPTION1("Cannot load variable \"time\" from NetCDF file \"%s\"",
						vecInputFileList[f].c_str());
				}
				if (varInTime->num_dims() != 1) {
					_EXCEPTION1("Variable \"time\" in NetCDF file \"%s\" has "
						"more than 1 dimension",
						vecInputFileList[f].c_str());
				}

				long lThisFileTimeCount = varInTime->get_dim(0)->size();

				// Verify calendar
				std::string strThisFileCalendar;
				NcAtt * attCalendar = varInTime->get_att("calendar");
				if (attCalendar == NULL) {
					strThisFileCalendar = "standard";
				} else {
					strThisFileCalendar = attCalendar->as_string(0);
				}
				if (strThisFileCalendar != strTimeCalendar) {
					_EXCEPTION3("Inconsistency in calendar in file \"%s\", "
						"expected \"%s\" found \"%s\"",
						vecInputFileList[f].c_str(),
						strTimeCalendar.c_str(),
						strThisFileCalendar.c_str());
				}

				// Load Times from file
				std::vector<Time> vecTimes;
				ReadCFTimeDataFromNcFile(
					&ncinfile, 
					vecInputFileList[f],
					vecTimes,
					false);

				// Get the variable
				NcVar * varIn = ncinfile.get_var(vecVariableStrings[v].c_str());
				if (varIn == NULL) {
					_EXCEPTION2("Cannot load variable \"%s\" from NetCDF file \"%s\"",
						vecVariableStrings[v].c_str(),
						vecInputFileList[f].c_str());
				}
				
				// Get the fill value
				float dFillValue;
				if (fMissingData) {
					NcAtt * attFillValue = varIn->get_att("_FillValue");
					if (attFillValue == NULL) {
						_EXCEPTION2("Variable \"%s\" missing _FillValue attribute, "
							"needed for --missingdata in NetCDF file \"%s\"",
							vecVariableStrings[v].c_str(),
							vecInputFileList[f].c_str());
					}

					dFillValue = attFillValue->as_float(0);
				}

				// TODO: Verify variable dimensionality

				// Daily mean
				for (int t = 0; t < vecTimes.size(); t++) {
					if (vecTimes[t].IsLeapDay()) {
						if (fVerbose) {
							Announce("Time %s (leap day; skipping)",
								vecTimes[t].ToString().c_str());
						}
						continue;
					} else {
						if (fVerbose) {
							Announce("Time %s",
								vecTimes[t].ToString().c_str());
						}
					}

					lPos0[0] = t;
					size_t sGetRows =
						NcGetPutSpecifiedDataSize<float>(
							varIn,
							lPos0,
							vecOutputAuxSize[v],
							c,
							dDataIn,
							NetCDF_Get);

					//std::cout << sGetRows << " read" << std::endl;

					Time timeJan01 = vecTimes[t];
					timeJan01.SetMonth(1);
					timeJan01.SetDay(1);

					int iDay = vecTimes[t].DayNumber() - timeJan01.DayNumber();
					if ((vecTimes[t].IsLeapYear()) && (iDay >= 60)) {
						iDay--;
					}

					_ASSERT((iDay >= 0) && (iDay < 365));

					if (fMissingData) {
						for (size_t i = 0; i < sGetRows; i++) {
							if (dDataIn[i] != dFillValue) {
								nTimeSlicesGrid[iDay][i]++;
								dAccumulatedData[iDay][i] += static_cast<double>(dDataIn[i]);
							}
						}
					} else {
						nTimeSlices[iDay]++;
						for (size_t i = 0; i < sGetRows; i++) {
							dAccumulatedData[iDay][i] += static_cast<double>(dDataIn[i]);
						}
					}
				}
			}

			// Calculate the mean over all data
			for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {

				// ..with missing data
				if (fMissingData) {
					for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
						if (nTimeSlicesGrid[t][i] != 0) {
							dAccumulatedData[t][i] /=
								static_cast<double>(nTimeSlicesGrid[t][i]);
						} else {
							dAccumulatedData[t][i] = 0.0;
						}
					}

				// ..without missing data
				} else {
					if (nTimeSlices[t] == 0) {
						continue;
					}
					for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
						dAccumulatedData[t][i] /= static_cast<double>(nTimeSlices[t]);
					}
				}
			}

			// Apply Fourier smoothing in time
			if ((nFourierModes != 0) &&
			    (nFourierModes < dAccumulatedData.GetColumns() / 2)
			) {
				AnnounceStartBlock("Performing Fourier filtering");
				DataArray1D<double> an(dAccumulatedData.GetRows());
				DataArray1D<double> bn(dAccumulatedData.GetRows());

				for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
					fourier_filter(
						&(dAccumulatedData[0][i]),
						dAccumulatedData.GetRows(),
						dAccumulatedData.GetColumns(),
						nFourierModes,
						an,
						bn);
				}
				AnnounceEndBlock("Done");
			}

			// Write data
			AnnounceStartBlock("Writing data to file");
			for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {
				// Convert to float and write to disk
				DataArray1D<float> dMeanData(dAccumulatedData.GetColumns());
				for (size_t i = 0; i < dMeanData.GetRows(); i++) {
					dMeanData[i] = static_cast<float>(dAccumulatedData[t][i]);
				}

				lPos0[0] = t;
				NcGetPutSpecifiedDataSize<float>(
					vecNcVarOut[v],
					lPos0,
					vecOutputAuxSize[v],
					c,
					dMeanData,
					NetCDF_Put);
			}
			AnnounceEndBlock("Done");

			if (vecOutputAuxCount[v] > 1) {
				AnnounceEndBlock("Done");
			}
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

