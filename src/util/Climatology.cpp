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
#include "FourierTransforms.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

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

///	<summary>
///		Calculate the number of iterations that need to be performed to
///		stay within the memory limit and the read size per iteration.
///	</summary>
template<typename Type>
bool CalculateIterationCountSize(
	size_t sMemoryMax,
	size_t sDataInstances,
	size_t sTotalValuesPerInstance,
	size_t & sAuxCount,
	size_t & sAuxSize
) {
	// Memory max too low; can't even store one value per data instance
	if (sMemoryMax < sDataInstances * sizeof(Type)) {
		return false;
	}

	// Total amount of space needed for storage
	size_t sClimatologyStorageSize = sDataInstances * sTotalValuesPerInstance * sizeof(Type);

	// Adjust memory maximum to be a multiple of sDataInstances * sizeof(Type)
	sMemoryMax -= sMemoryMax % (sDataInstances * sizeof(Type));

	if (sClimatologyStorageSize <= sMemoryMax) {
		sAuxCount = 1;
		sAuxSize = sTotalValuesPerInstance;

	} else if (sClimatologyStorageSize % sMemoryMax == 0) {
		sAuxCount = sClimatologyStorageSize / sMemoryMax;
		sAuxSize = sMemoryMax / (sDataInstances * sizeof(Type));

	} else {
		sAuxCount = sClimatologyStorageSize / sMemoryMax + 1;
		sAuxSize = sMemoryMax / (sDataInstances * sizeof(Type));
	}

	//printf("%lu %lu\n", sDataInstances, sTotalValuesPerInstance);
	//printf("%lu %lu\n", sClimatologyStorageSize, sMemoryMax);
	//printf("%lu %lu (%lu)\n", sAuxCount, sAuxSize, sAuxSize);
	//_EXCEPTION();
	return true;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the next nD NetCDF file position (lPos) and size (lSize)
///		given a 1D position (llNcArrayPos) and size (llNcArraySize).
///	</summary>
bool CalculatePartialDataPosSize(
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

	if (llNcArrayPos == llAccumulatedSize[lSkipDims]) {
		return false;
	}

	// Given a 1D llNcArrayPos find the multi-dimensional index corresponding
	// to this offset in the larger array.
	long long llNcArrayPosTemp = llNcArrayPos;
	for (long d = lDimSize.GetRows()-1; d >= lSkipDims; d--) {
		lPos[d] = llNcArrayPosTemp % lDimSize[d];
		llNcArrayPosTemp /= lDimSize[d];
	}
	if (llNcArrayPosTemp > 0) {
		_EXCEPTION2("NcVar 1D index (%lu) out of range (%lu)", llNcArrayPos, llNcArrayPosTemp);
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

	return true;
}

///////////////////////////////////////////////////////////////////////////////

enum NetCDF_GetPut {
	NetCDF_Get,
	NetCDF_Put
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load in exactly lSizePerIteration contiguous elements from a NcVar,
///		starting at the specified position.
///	</summary>
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
		bool fDataRemaining =
			CalculatePartialDataPosSize(
				lDimSize,
				lPos0.GetRows(),
				llNcArrayPos,
				lSizePerIteration,
				lPos,
				lSize);

		if (!fDataRemaining) {
			break;
		}

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

		NcError err;
		if (err.get_err() != NC_NOERR) {
			_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
		}

		if (eGetPut == NetCDF_Get) {
			ncvar->get(&(data[llTotalReadWriteSize]), &(lSize[0]));
		} else if (eGetPut == NetCDF_Put) {
			ncvar->put(&(data[llTotalReadWriteSize]), &(lSize[0]));
		} else {
			_EXCEPTIONT("Logic error");
		}
		if (err.get_err() != NC_NOERR) {
			_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
		}

		llTotalReadWriteSize += lReadWriteSize;
	}

	return llTotalReadWriteSize;
}

///////////////////////////////////////////////////////////////////////////////

void SplitVariableStrings(
	const std::vector<std::string> & vecVariableStrings,
	std::vector<std::string> & vecVariableNames,
	std::vector< std::vector<long> > & vecVariableSpecifiedDims
) {
	vecVariableNames.clear();
	vecVariableSpecifiedDims.clear();

	vecVariableNames.resize(vecVariableStrings.size());
	vecVariableSpecifiedDims.resize(vecVariableStrings.size());

	// Loop through all variables
	for (int v = 0; v < vecVariableStrings.size(); v++) {

		// Split variable string into variable name and any specified dimensions
		int iBeginParentheses = (-1);
		int iEndParentheses = (-1);
		for (int i = 0; i < vecVariableStrings[v].length(); i++) {
			if (vecVariableStrings[v][i] == '(') {
				if (iBeginParentheses != (-1)) {
					_EXCEPTIONT("Multiple open parentheses in --var");
				} else {
					iBeginParentheses = i;
				}
			}
			if (vecVariableStrings[v][i] == ')') {
				if (iEndParentheses != (-1)) {
					_EXCEPTIONT("Multiple closed parentheses in --var");
				} else {
					iEndParentheses = i;
				}
			}
		}
	
		// Extract variable name and specified dimensions
		if (iBeginParentheses != (-1)) {
			if ((iBeginParentheses != (-1)) && (iEndParentheses == (-1))) {
				_EXCEPTIONT("Unbalanced open parentheses in --var");
			}
			if ((iBeginParentheses == (-1)) && (iEndParentheses != (-1))) {
				_EXCEPTIONT("Unbalanced closed parentheses in --var");
			}
			if (iBeginParentheses >= iEndParentheses) {
				_EXCEPTIONT("Unbalanced closed parentheses in --var");
			}
	
			vecVariableNames[v] =
				vecVariableStrings[v].substr(0, iBeginParentheses);
	
			int iLast = iBeginParentheses+1;
			for (int i = iBeginParentheses+1; i <= iEndParentheses; i++) {
				if ((i == iEndParentheses) ||
				    (vecVariableStrings[v][i] == ',')
				) {
					std::string strDimValue =
						vecVariableStrings[v].substr(iLast, i-iLast);
					//if (strDimValue == "*") {
					//	vecVariableSpecifiedDims[v].push_back(-1);
					//	iLast = i+1;
					//}
					if (!STLStringHelper::IsInteger(strDimValue)) {
						_EXCEPTION1("Invalid dimension index \"%s\" in --var; expected positive integer.",
							strDimValue.c_str());
					}
					long lDimValue = std::stol(strDimValue);
					if (lDimValue < 0) {
						_EXCEPTION1("Invalid dimension index \"%s\" in --var; expected positive integer.",
							strDimValue.c_str());
					}
					vecVariableSpecifiedDims[v].push_back(lDimValue);
					iLast = i+1;
				}
			}
		} else {
			vecVariableNames[v] = vecVariableStrings[v];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Possible climatology types that can be calculated with Climatology().
///	</summary>
enum ClimatologyType {
	ClimatologyType_Mean,
	ClimatologyType_Sq
};

///	<summary>
///		Produce a climatology over the given input files.
///	</summary>
void Climatology(
	const std::vector<std::string> & vecInputFileList,
	const std::string & strOutputFile,
	const std::vector<std::string> & vecVariableNames,
	const std::vector< std::vector<long> > & vecVariableSpecifiedDims,
	size_t sMemoryMax,
	bool fIncludeLeapDays,
	ClimatologyType eClimoType,
	int nFourierModes,
	bool fMissingData,
	bool fOutputSliceCounts,
	bool fVerbose
) {
	_ASSERT(vecInputFileList.size() > 0);
	_ASSERT(vecVariableNames.size() > 0);
	_ASSERT(vecVariableNames.size() == vecVariableSpecifiedDims.size());

	// Climatology type
	if ((eClimoType != ClimatologyType_Mean) &&
	    (eClimoType != ClimatologyType_Sq)
	) {
		_EXCEPTIONT("Invalid eClimoType");
	}

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
	vecVarDimInfo.resize(vecVariableNames.size());

	std::vector<NcVar*> vecNcVarOut;
	vecNcVarOut.resize(vecVariableNames.size());

	std::vector<NcVar*> vecNcVarCount;
	vecNcVarCount.resize(vecVariableNames.size());

	std::string strTimeCalendar;
	size_t sOutputTimes = 0;

	std::vector<size_t> vecOutputAuxSize;
	vecOutputAuxSize.resize(vecVariableNames.size(), 0);

	std::vector<size_t> vecOutputAuxCount;
	vecOutputAuxCount.resize(vecVariableNames.size(), 0);

	/////////////////////////////////////////////
	// Initialization
	{
		// Copy over details from the input file
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
		for (int v = 0; v < vecVariableNames.size(); v++) {

			// Get the variable
			NcVar * var = ncinfile.get_var(vecVariableNames[v].c_str());
			if (var == NULL) {
				_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
					vecInputFileList[0].c_str(),
					vecVariableNames[v].c_str());
			}

			if (var->num_dims() - vecVariableSpecifiedDims[v].size() < 2) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have at least 2 dimensions (in addition to specified dimensions)",
					vecInputFileList[0].c_str(),
					vecVariableNames[v].c_str());
			}

			NcDim * dim0 = var->get_dim(0);
			_ASSERT(dim0 != NULL);

			std::string strDimTimeName = dim0->name();
			if (strDimTimeName != "time") {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have leftmost dimension \"time\"",
					vecInputFileList[0].c_str(),
					vecVariableNames[v].c_str());
			}

			for (int d = 0; d < vecVariableSpecifiedDims[v].size(); d++) {
				if ((vecVariableSpecifiedDims[v][d] < (-1)) ||
				    (vecVariableSpecifiedDims[v][d] >= var->get_dim(d+1)->size())
				) {
					_EXCEPTION4("File \"%s\" variable \"%s\" specified dimension index out of range (%lu/%lu).",
					vecInputFileList[0].c_str(),
					vecVariableNames[v].c_str(),
					vecVariableSpecifiedDims[v][d],
					var->get_dim(d+1)->size());
				}
			}

			for (int d = 1; d < var->num_dims(); d++) {
				NcDim * dimVar = var->get_dim(d);
				DimInfo diminfo(dimVar->name(), dimVar->size());
				vecVarDimInfo[v].push_back(diminfo);
				setDimInfo.insert(diminfo);
			}

			// Output storage requirements
			size_t sTotalAuxDims = 1;
			for (int d = vecVariableSpecifiedDims[v].size(); d < vecVarDimInfo[v].size(); d++) {
				sTotalAuxDims *= vecVarDimInfo[v][d].size;
			}

			// Number of auxiliary values to store at each iteration
			// We require storage for sOutputTimes persistent arrays and 1 temporary array
			bool fSufficientMemory =
				CalculateIterationCountSize<double>(
					sMemoryMax,
					sOutputTimes + 1,
					sTotalAuxDims,
					vecOutputAuxCount[v],
					vecOutputAuxSize[v]);
/*
			// Only allow one read per iteration
			size_t sAuxCount = 1;
			size_t sAuxSize = 1;
			for (int d = var->num_dims()-1; d >= 1; d--) {
				size_t sDimSize = var->get_dim(d)->size();
				if (sAuxSize * sDimSize * sOutputTimes * sizeof(double) > sMemoryMax) {
					for (int di = d-1; d >= 1; d--) {
						sAuxCount *= sDimSize;
					}
					break;
				}
				sAuxSize *= sDimSize;
			}

			Announce("Multi-read: %lu (%lu) / Single-read: %lu (%lu)", vecOutputAuxSize[v], vecOutputAuxCount[v], sAuxSize, sAuxCount);
			vecOutputAuxSize[v] = sAuxSize;
			vecOutputAuxCount[v] = sAuxCount;
*/
			// Output information to user
			std::string strPerTimesliceStorage =
				MemoryBytesToString(sTotalAuxDims * sizeof(double));

			std::string strTotalOutputStorage =
				MemoryBytesToString(sOutputTimes * sTotalAuxDims * sizeof(double));

			AnnounceStartBlock("Variable \"%s\" has %lu distinct time series",
				vecVariableNames[v].c_str(),
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

		if (fIncludeLeapDays) {
			_EXCEPTIONT("--include_leap_days not implemented");
		} else {
			varOutTime->add_att("calendar", "noleap");
		}
		varOutTime->add_att("units", strTimeUnits.c_str());
		varOutTime->add_att("type", "daily mean climatology");

		if (varOutTime->type() == ncInt) {
			DataArray1D<int> nTimes(sOutputTimes);
			for (size_t s = 0; s < sOutputTimes; s++) {
				nTimes[s] = static_cast<int>(s);
			}
			varOutTime->put(&(nTimes[0]), sOutputTimes);

		} else if (varOutTime->type() == ncDouble) {
			DataArray1D<double> dTimes(sOutputTimes);
			for (size_t s = 0; s < sOutputTimes; s++) {
				dTimes[s] = static_cast<double>(s);
			}
			varOutTime->put(&(dTimes[0]), sOutputTimes);

		} else if (varOutTime->type() == ncInt64) {
			DataArray1D<ncint64> dTimes(sOutputTimes);
			for (size_t s = 0; s < sOutputTimes; s++) {
				dTimes[s] = static_cast<ncint64>(s);
			}
			varOutTime->put(&(dTimes[0]), sOutputTimes);

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
		for (int v = 0; v < vecVariableNames.size(); v++) {

			// Create the vector NcDims for the output variable
			std::vector<NcDim *> vecVarNcDims;
			vecVarNcDims.push_back(dimOutTime);
			for (int d = vecVariableSpecifiedDims[v].size(); d < vecVarDimInfo[v].size(); d++) {
				auto itNcDim = mapOutputNcDim.find(vecVarDimInfo[v][d].name);
				_ASSERT(itNcDim != mapOutputNcDim.end());
				vecVarNcDims.push_back(itNcDim->second);
			}
			std::string strOutVar = std::string("dailymean_") + vecVariableNames[v];
			NcVar * varOut =
				ncoutfile.add_var(
					strOutVar.c_str(),
					ncFloat,
					vecVarNcDims.size(),
					const_cast<const NcDim**>(&(vecVarNcDims[0])));

			if (varOut == NULL) {
				_EXCEPTION1("Unable to create output variable \"%s\" in output file",
					strOutVar.c_str());
			}

			// Create output variable for slice count
			if (fOutputSliceCounts) {
				if (fMissingData) {
					_EXCEPTIONT("Output slice counts with missing data not implemented");
				}
				std::string strOutCountVar = strOutVar + "_count";
				NcVar * varOutCount =
					ncoutfile.add_var(
						strOutCountVar.c_str(),
						ncInt,
						dimOutTime);

				if (varOutCount == NULL) {
					_EXCEPTION1("Unable to create output variable \"%s\" in output file",
						strOutCountVar.c_str());
				}

				vecNcVarCount[v] = varOutCount;
			}

			// Add attributes to variable describing specified hyperslab indices
			for (int d = 0; d < vecVariableSpecifiedDims[v].size(); d++) {
				NcVar * varDimIn = ncinfile.get_var(vecVarDimInfo[v][d].name.c_str());
				if (varDimIn == NULL) {
					std::string strAttName = vecVarDimInfo[v][d].name + "_index";
					varOut->add_att(
						strAttName.c_str(),
						vecVariableSpecifiedDims[v][d]);
				} else {
					_ASSERT(varDimIn->num_dims() == 1);
					_ASSERT(varDimIn->get_dim(0)->size() > vecVariableSpecifiedDims[v][d]);
					varDimIn->set_cur(vecVariableSpecifiedDims[v][d]);
					double dValue;
					varDimIn->get(&dValue, (long)1);
					varOut->add_att(
						vecVarDimInfo[v][d].name.c_str(),
						dValue);
				}
			}

			// Get the input NcVar
			NcVar * varIn = ncinfile.get_var(vecVariableNames[v].c_str());
			_ASSERT(varIn != NULL);
			CopyNcVarAttributes(varIn, varOut);

			vecNcVarOut[v] = varOut;
		}

		AnnounceEndBlock("Done");
	}

	/////////////////////////////////////////////
	// Actual climatology calculation

	// Loop through each variable
	for (int v = 0; v < vecVariableNames.size(); v++) {
		AnnounceStartBlock("Calculating climatology of variable \"%s\"",
			vecVariableNames[v].c_str());

		DataArray1D<float> dDataIn(vecOutputAuxSize[v]);

		// Storage array for auxiliary dimensions
		DataArray1D<long> lPos0get(1 + vecVariableSpecifiedDims[v].size());
		for (int d = 0; d < vecVariableSpecifiedDims[v].size(); d++) {
			lPos0get[d+1] = vecVariableSpecifiedDims[v][d];
		}
		DataArray1D<long> lPos0put(1);

		// Loop through each iteration of this variable
		for (int c = 0; c < vecOutputAuxCount[v]; c++) {

			if (vecOutputAuxCount[v] > 1) {
				AnnounceStartBlock("Iteration %i/%i", c+1, vecOutputAuxCount[v]);
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
				NcVar * varIn = ncinfile.get_var(vecVariableNames[v].c_str());
				if (varIn == NULL) {
					_EXCEPTION2("Cannot load variable \"%s\" from NetCDF file \"%s\"",
						vecVariableNames[v].c_str(),
						vecInputFileList[f].c_str());
				}

				// Get the fill value
				float dFillValue;
				if (fMissingData) {
					NcAtt * attFillValue = varIn->get_att("_FillValue");
					if (attFillValue == NULL) {
						_EXCEPTION2("Variable \"%s\" missing _FillValue attribute, "
							"needed for --missingdata in NetCDF file \"%s\"",
							vecVariableNames[v].c_str(),
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

					lPos0get[0] = t;
					size_t sGetRows =
						NcGetPutSpecifiedDataSize<float>(
							varIn,
							lPos0get,
							vecOutputAuxSize[v],
							c,
							dDataIn,
							NetCDF_Get);

					Time timeJan01 = vecTimes[t];
					timeJan01.SetMonth(1);
					timeJan01.SetDay(1);

					int iDay = vecTimes[t].DayNumber() - timeJan01.DayNumber();
					if ((vecTimes[t].IsLeapYear()) && (iDay >= 60)) {
						iDay--;
					}

					_ASSERT((iDay >= 0) && (iDay < 365));

					if (fMissingData) {
						if (eClimoType == ClimatologyType_Mean) {
							for (size_t i = 0; i < sGetRows; i++) {
								if (dDataIn[i] != dFillValue) {
									nTimeSlicesGrid[iDay][i]++;
									dAccumulatedData[iDay][i] += static_cast<double>(dDataIn[i]);
								}
							}

						} else if (eClimoType == ClimatologyType_Sq) {
							for (size_t i = 0; i < sGetRows; i++) {
								if (dDataIn[i] != dFillValue) {
									nTimeSlicesGrid[iDay][i]++;
									dAccumulatedData[iDay][i] +=
										static_cast<double>(dDataIn[i])
										* static_cast<double>(dDataIn[i]);
								}
							}
						}

					} else {
						nTimeSlices[iDay]++;
						if (eClimoType == ClimatologyType_Mean) {
							for (size_t i = 0; i < sGetRows; i++) {
								dAccumulatedData[iDay][i] += static_cast<double>(dDataIn[i]);
							}

						} else if (eClimoType == ClimatologyType_Sq) {
							for (size_t i = 0; i < sGetRows; i++) {
								dAccumulatedData[iDay][i] +=
									static_cast<double>(dDataIn[i])
									* static_cast<double>(dDataIn[i]);
							}
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
/*
			// Apply Fourier smoothing in time
			if ((nFourierModes != 0) &&
			    (nFourierModes < dAccumulatedData.GetColumns() / 2)
			) {
				AnnounceStartBlock("Performing Fourier filtering");
				DataArray1D<double> an(dAccumulatedData.GetRows());
				DataArray1D<double> bn(dAccumulatedData.GetRows());

				for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
					fourier_filter<double>(
						&(dAccumulatedData[0][i]),
						dAccumulatedData.GetRows(),
						dAccumulatedData.GetColumns(),
						nFourierModes,
						an,
						bn);
				}
				AnnounceEndBlock("Done");
			}
*/
			// Write data
			AnnounceStartBlock("Writing data to file");
			for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {

				// Convert to float and write to disk
				DataArray1D<float> dMeanData(dAccumulatedData.GetColumns());
				for (size_t i = 0; i < dMeanData.GetRows(); i++) {
					dMeanData[i] = static_cast<float>(dAccumulatedData[t][i]);
				}

				lPos0put[0] = t;
				NcGetPutSpecifiedDataSize<float>(
					vecNcVarOut[v],
					lPos0put,
					vecOutputAuxSize[v],
					c,
					dMeanData,
					NetCDF_Put);
			}

			// Write slice counts
			if (fOutputSliceCounts) {
				AnnounceStartBlock("Writing slice counts");
				if (fMissingData) {
					_EXCEPTION();

				// If no missing data only need to write this once
				} else {
					if (c == 0) {
						vecNcVarCount[v]->set_cur((long)0);
						vecNcVarCount[v]->put(&(nTimeSlices[0]), dAccumulatedData.GetRows());
					}
				}
				AnnounceEndBlock("Done");
			}

			AnnounceEndBlock("Done");

			if (vecOutputAuxCount[v] > 1) {
				AnnounceEndBlock("Done");
			}
		}
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

	// Period of the climatology
	std::string strPeriod;

	// Type of climatology
	std::string strClimoType;

	// Include leap days
	bool fIncludeLeapDays;

	// Number of Fourier modes to retain
	int nFourierModes;

	// Dataset has missing data
	bool fMissingData;

	// Path to write temporary cliamtology files
	std::string strTempFilePath;

	// Do not delete temp files after completing climatology
	bool fKeepTempFiles;

	// Verbose
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strOutputFile, "out_data", "");
		CommandLineString(strVariable, "var", "");
		CommandLineStringD(strMemoryMax, "memmax", "2G", "[#K,#M,#G]");
		CommandLineStringD(strPeriod, "period", "daily", "[daily|monthly|seasonal|annual]");
		CommandLineStringD(strClimoType, "type", "mean", "[mean|sq]");
		CommandLineBool(fIncludeLeapDays, "include_leap_days");
		//CommandLineInt(nFourierModes, "time_modes", 0);
		CommandLineBool(fMissingData, "missingdata");
		CommandLineString(strTempFilePath, "temp_file_path", ".");
		CommandLineBool(fKeepTempFiles, "keep_temp_files");
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

	// Climatology type
	ClimatologyType eClimoType;
	STLStringHelper::ToLower(strClimoType);
	if (strClimoType == "mean") {
		eClimoType = ClimatologyType_Mean;
	} else if (strClimoType == "sq") {
		eClimoType = ClimatologyType_Sq;
	} else {
		_EXCEPTIONT("--type invalid; expected \"mean\" or \"sq\"");
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

	std::vector<std::string> vecVariableNames;
	std::vector< std::vector<long> > vecVariableSpecifiedDims;

	SplitVariableStrings(
		vecVariableStrings,
		vecVariableNames,
		vecVariableSpecifiedDims
	);

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

#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
#else
	int nMPIRank = 0;
	int nMPISize = 1;
#endif

	// Distribute workload over processors
	std::vector< std::string > vecMyInputFileList;
	std::string strMyOutputFile;

	bool fOutputSliceCounts = false;

	if (nMPISize != 1) {
		fOutputSliceCounts = true;

		size_t sFilesPerProc = vecInputFileList.size() / nMPISize;

		// Split the file list among processors
		for (size_t f = 0; f < sFilesPerProc; f++) {
			vecMyInputFileList.push_back(
				vecInputFileList[sFilesPerProc * nMPIRank + f]);
		}
		if (nMPIRank == nMPISize-1) {
			for (size_t f = sFilesPerProc * nMPISize; f < vecInputFileList.size(); f++) {
				vecMyInputFileList.push_back(vecInputFileList[f]);
			}
		}

		char szBuffer[16];
		sprintf(szBuffer, "%06i", nMPIRank);

		strMyOutputFile =
			strTempFilePath + "/tempClimatology" + szBuffer + ".nc";

		std::string strLogFile =
			strTempFilePath + "/tempLog" + szBuffer + ".txt";

		FILE * fp = fopen(strLogFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Unable to open log file \"%s\"", strLogFile.c_str());
		}

		Announce("Logs will be written to \"%s/tempLog######.txt\"", strTempFilePath.c_str());
		Announce("Temporary climatologies will be written to \"%s/tempClimatology######.nc\"", strTempFilePath.c_str());
		AnnounceSetOutputBuffer(fp);
		AnnounceOutputOnAllRanks();

	} else {
		vecMyInputFileList = vecInputFileList;
		strMyOutputFile = strOutputFile;
	}

	// Calculate the climatology
	Climatology(
		vecMyInputFileList,
		strMyOutputFile,
		vecVariableNames,
		vecVariableSpecifiedDims,
		sMemoryMax,
		fIncludeLeapDays,
		eClimoType,
		nFourierModes,
		fMissingData,
		fOutputSliceCounts,
		fVerbose
	);

#if defined(TEMPEST_MPIOMP)
	if (nMPISize != 1) {
		MPI_Barrier(MPI_COMM_WORLD);
		AnnounceSetOutputBuffer(stdout);
		AnnounceOnlyOutputOnRankZero();
	}
#endif

	if ((nMPISize != 1) && (nMPIRank == 0)) {
		AnnounceBanner();

		if (nMPIRank == 0) {
			Announce("TODO: Combine files");
			DataArray1D<float> dData;
		}
	}

	AnnounceBanner();

} catch(Exception & e) {
	AnnounceOutputOnAllRanks();
	AnnounceSetOutputBuffer(stdout);
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

