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
#include "Units.h"

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
///	<returns>
///		true if sufficient memory is available after splitting, false
///		otherwise.
///	</returns>
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
	//printf("%lu %lu\n", sAuxCount, sAuxSize);
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
	std::vector< std::vector<std::string> > & vecVariableSpecifiedDims
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

					std::string strValue;
					std::string strUnits;
					SplitIntoValueAndUnits(strDimValue, strValue, strUnits);

					if (!STLStringHelper::IsFloat(strValue)) {
						_EXCEPTION1("Invalid dimension index \"%s\" in --var; expected positive integer or value index.",
							strDimValue.c_str());
					}

					//if (strDimValue == "*") {
					//	vecVariableSpecifiedDims[v].push_back(-1);
					//	iLast = i+1;
					//}
					//long lDimValue = std::stol(strDimValue);
					//if (lDimValue < 0) {
					//	_EXCEPTION1("Invalid dimension index \"%s\" in --var; expected positive integer.",
					//		strDimValue.c_str());
					//}
					vecVariableSpecifiedDims[v].push_back(strDimValue);
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
///		Possible climatology periods that can be calculated with Climatology().
///	</summary>
enum ClimatologyPeriod {
	ClimatologyPeriod_Daily,
	ClimatologyPeriod_Monthly,
	ClimatologyPeriod_Seasonal,
	ClimatologyPeriod_Annual
};

///	<summary>
///		Possible climatology types that can be calculated with Climatology().
///	</summary>
enum ClimatologyType {
	ClimatologyType_Mean,
	ClimatologyType_MeanSq
};

///	<summary>
///		Produce a climatology over the given input files.
///	</summary>
void Climatology(
	const std::vector<std::string> & vecInputFileList,
	const std::string & strOutputFile,
	const std::vector<std::string> & vecVariableNames,
	const std::vector< std::vector<std::string> > & vecVariableSpecifiedDims,
	size_t sMemoryMax,
	bool fIncludeLeapDays,
	ClimatologyPeriod eClimoPeriod,
	ClimatologyType eClimoType,
	int nFourierModes,
	bool fMissingData,
	bool fOutputSliceCounts,
	bool fVerbose
) {
	_ASSERT(vecInputFileList.size() > 0);
	_ASSERT(vecVariableNames.size() > 0);
	_ASSERT(vecVariableNames.size() == vecVariableSpecifiedDims.size());

	// Climatology period
	if ((eClimoPeriod != ClimatologyPeriod_Daily) &&
		(eClimoPeriod != ClimatologyPeriod_Monthly) &&
		(eClimoPeriod != ClimatologyPeriod_Seasonal) &&
		(eClimoPeriod != ClimatologyPeriod_Annual)
	) {
		_EXCEPTIONT("Invalid eClimoPeriod");
	}

	// Climatology type
	if ((eClimoType != ClimatologyType_Mean) &&
	    (eClimoType != ClimatologyType_MeanSq)
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
	std::vector< std::vector<long> > vecVariableSpecifiedDimIxs;
	vecVariableSpecifiedDimIxs.resize(vecVariableSpecifiedDims.size());

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

			// Other climatologies
			if (eClimoPeriod == ClimatologyPeriod_Monthly) {
				sOutputTimes = 12;
			} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
				sOutputTimes = 4;
			} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
				sOutputTimes = 1;
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

			// Get the dimension indices
			vecVariableSpecifiedDimIxs[v].resize(vecVariableSpecifiedDims[v].size());
			for (int d = 0; d < vecVariableSpecifiedDims[v].size(); d++) {

				vecVariableSpecifiedDimIxs[v][d] =
					GetIntegerIndexFromValueBasedIndex(
						&ncinfile,
						vecInputFileList[0],
						vecVariableNames[v],
						d+1,
						vecVariableSpecifiedDims[v][d]
					);

				if ((vecVariableSpecifiedDimIxs[v][d] < (-1)) ||
				    (vecVariableSpecifiedDimIxs[v][d] >= var->get_dim(d+1)->size())
				) {
					_EXCEPTION4("File \"%s\" variable \"%s\" specified dimension index out of range (%lu/%lu).",
						vecInputFileList[0].c_str(),
						vecVariableNames[v].c_str(),
						vecVariableSpecifiedDimIxs[v][d],
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

			std::string strFullVariableName = vecVariableNames[v];
			if (vecVariableSpecifiedDims[v].size() > 0) {
				strFullVariableName += "(";
				for (int d = 0; d < vecVariableSpecifiedDims[v].size(); d++) {
					strFullVariableName += vecVariableSpecifiedDims[v][d];
					if (d != vecVariableSpecifiedDims.size()-1) {
						strFullVariableName += ",";
					}
				}
				strFullVariableName += ")";
			}

			AnnounceStartBlock("Variable \"%s\" has %lu distinct time series",
				strFullVariableName.c_str(),
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
		//NcDim * dimOutTime = ncoutfile.add_dim("time", (long)sOutputTimes);
		//if (dimOutTime == NULL) {
		//	_EXCEPTIONT("Unable to create output dimension \"time\" in output file");
		//}

		//NcVar * varOutTime = ncoutfile.add_var("time", varInTime->type(), dimOutTime);
		//if (varOutTime == NULL) {
		//	_EXCEPTIONT("Unable to create output dimension \"time\" in output file");
		//}

		//if (fIncludeLeapDays) {
		//	_EXCEPTIONT("--include_leap_days not implemented");
		//} else {
		//	varOutTime->add_att("calendar", "noleap");
		//}
		//varOutTime->add_att("units", strTimeUnits.c_str());

		NcTimeDimension vecTimes;
		vecTimes.m_nctype = varInTime->type();
		vecTimes.m_units = strTimeUnits;

		if (eClimoPeriod == ClimatologyPeriod_Daily) {
			Time timeOut("0001-01-01-00000-000000", Time::CalendarNoLeap);
			for (size_t s = 0; s < sOutputTimes; s++) {
				vecTimes.push_back(timeOut);
				timeOut.AddDays(1);
			}

		} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
			Time timeOut("0001-01-01-00000-000000", Time::CalendarNoLeap);
			for (size_t s = 0; s < sOutputTimes; s++) {
				vecTimes.push_back(timeOut);
				timeOut.AddMonths(1);
			}

		} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
			Time timeOut("0001-01-01-00000-000000", Time::CalendarNoLeap);
			for (size_t s = 0; s < sOutputTimes; s++) {
				vecTimes.push_back(timeOut);
				timeOut.AddMonths(3);
			}

		} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
			Time timeOut("0001-01-01-00000-000000", Time::CalendarNoLeap);
			vecTimes.push_back(timeOut);
		}

		WriteCFTimeDataToNcFile(
			&ncoutfile,
			strOutputFile,
			vecTimes,
			true);

		NcDim * dimOutTime = ncoutfile.get_dim("time");
		_ASSERT(dimOutTime != NULL);

		NcVar * varOutTime = ncoutfile.get_var("time");
		_ASSERT(varOutTime != NULL);
		if (eClimoPeriod == ClimatologyPeriod_Daily) {
			varOutTime->add_att("type", "daily mean climatology");
		} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
			varOutTime->add_att("type", "monthly mean climatology");
		} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
			varOutTime->add_att("type", "seasonal mean climatology");
		} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
			varOutTime->add_att("type", "annual mean climatology");
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
			std::string strOutVar;
			if (eClimoType == ClimatologyType_Mean) {
				strOutVar = std::string("dailymean_") + vecVariableNames[v];
			} else if (eClimoType == ClimatologyType_MeanSq) {
				strOutVar = std::string("dailymeansq_") + vecVariableNames[v];
			} else {
				_EXCEPTIONT("Invalid eClimoType");
			}
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
				std::string strOutCountVar = strOutVar + "_count";

				if (fMissingData) {
					vecNcVarCount[v] =
						ncoutfile.add_var(
							strOutCountVar.c_str(),
							ncInt,
							vecVarNcDims.size(),
							const_cast<const NcDim**>(&(vecVarNcDims[0])));

				} else {
					vecNcVarCount[v] =
						ncoutfile.add_var(
							strOutCountVar.c_str(),
							ncInt,
							dimOutTime);
				}

				if (vecNcVarCount[v] == NULL) {
					_EXCEPTION1("Unable to create output variable \"%s\" in output file",
						strOutCountVar.c_str());
				}
			}

			// Add attributes to variable describing specified hyperslab indices
			for (int d = 0; d < vecVariableSpecifiedDims[v].size(); d++) {
				NcVar * varDimIn = ncinfile.get_var(vecVarDimInfo[v][d].name.c_str());
				if (varDimIn == NULL) {
					std::string strAttName = vecVarDimInfo[v][d].name + "_index";
					varOut->add_att(
						strAttName.c_str(),
						vecVariableSpecifiedDims[v][d].c_str());
				} else {
					_ASSERT(varDimIn->num_dims() == 1);
					_ASSERT(varDimIn->get_dim(0)->size() > vecVariableSpecifiedDimIxs[v][d]);
					varDimIn->set_cur(vecVariableSpecifiedDimIxs[v][d]);
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
		DataArray1D<long> lPos0get(1 + vecVariableSpecifiedDimIxs[v].size());
		for (int d = 0; d < vecVariableSpecifiedDimIxs[v].size(); d++) {
			lPos0get[d+1] = vecVariableSpecifiedDimIxs[v][d];
		}
		DataArray1D<long> lPos0put(1);

		float dFillValue = 0.0;

		// Loop through each iteration of this variable
		for (int c = 0; c < vecOutputAuxCount[v]; c++) {

			if (vecOutputAuxCount[v] > 1) {
				AnnounceStartBlock("Iteration %i/%i", c+1, vecOutputAuxCount[v]);
			}

			// Number of time slices
			DataArray1D<int> nTimeSlices;
			if (!fMissingData) {
				nTimeSlices.Allocate(sOutputTimes);
			}

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

				if (fVerbose) {
					Announce("%s", vecInputFileList[f].c_str());
				}

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
				NcTimeDimension vecTimes;
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

					int iTimeIndex;
					if (eClimoPeriod == ClimatologyPeriod_Daily) {
						iTimeIndex = vecTimes[t].DayNumber() - timeJan01.DayNumber();
						if ((vecTimes[t].IsLeapYear()) && (iTimeIndex >= 60)) {
							iTimeIndex--;
						}

						_ASSERT((iTimeIndex >= 0) && (iTimeIndex < 365));

					} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
						iTimeIndex = vecTimes[t].GetZeroIndexedMonth();
						_ASSERT((iTimeIndex >= 0) && (iTimeIndex < 12));

					} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
						int iMonth = vecTimes[t].GetZeroIndexedMonth();
						_ASSERT((iMonth >= 0) && (iMonth < 12));
						if ((iMonth == 11) || (iMonth == 0) || (iMonth == 1)) {
							iTimeIndex = 0;
						} else if ((iMonth == 2) || (iMonth == 3) || (iMonth == 4)) {
							iTimeIndex = 1;
						} else if ((iMonth == 5) || (iMonth == 6) || (iMonth == 7)) {
							iTimeIndex = 2;
						} else if ((iMonth == 8) || (iMonth == 9) || (iMonth == 10)) {
							iTimeIndex = 3;
						}

					} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
						iTimeIndex = 0;
					}

					// Count number of time slices at each point and accumulate data
					if (fMissingData) {
						if (eClimoType == ClimatologyType_Mean) {
							for (size_t i = 0; i < sGetRows; i++) {
								if (dDataIn[i] != dFillValue) {
									nTimeSlicesGrid[iTimeIndex][i]++;
									dAccumulatedData[iTimeIndex][i] += static_cast<double>(dDataIn[i]);
								}
							}

						} else if (eClimoType == ClimatologyType_MeanSq) {
							for (size_t i = 0; i < sGetRows; i++) {
								if (dDataIn[i] != dFillValue) {
									nTimeSlicesGrid[iTimeIndex][i]++;
									dAccumulatedData[iTimeIndex][i] +=
										static_cast<double>(dDataIn[i])
										* static_cast<double>(dDataIn[i]);
								}
							}
						}

					// Count number of time slices at each time and accumulate data
					} else {
						nTimeSlices[iTimeIndex]++;
						if (eClimoType == ClimatologyType_Mean) {
							for (size_t i = 0; i < sGetRows; i++) {
								dAccumulatedData[iTimeIndex][i] += static_cast<double>(dDataIn[i]);
							}

						} else if (eClimoType == ClimatologyType_MeanSq) {
							for (size_t i = 0; i < sGetRows; i++) {
								dAccumulatedData[iTimeIndex][i] +=
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
							dAccumulatedData[t][i] = dFillValue;
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
			if (fMissingData) {
				AnnounceStartBlock("Writing data and slice counts to file");
			} else {
				AnnounceStartBlock("Writing data to file");
			}
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

				// If missing data is present write slice counts
				if (fOutputSliceCounts && fMissingData) {

					DataArray1D<int> nTimeSlicesGrid1D(nTimeSlicesGrid.GetColumns(), false);
					nTimeSlicesGrid1D.AttachToData(&(nTimeSlicesGrid[t][0]));

					NcGetPutSpecifiedDataSize<int>(
						vecNcVarCount[v],
						lPos0put,
						vecOutputAuxSize[v],
						c,
						nTimeSlicesGrid1D,
						NetCDF_Put);

					AnnounceEndBlock("Done");
				}
			}

			// Write slice counts
			if (fOutputSliceCounts && (!fMissingData) && (c == 0)) {
				AnnounceStartBlock("Writing slice counts");
				vecNcVarCount[v]->set_cur((long)0);
				vecNcVarCount[v]->put(&(nTimeSlices[0]), dAccumulatedData.GetRows());
				AnnounceEndBlock("Done");
			}

			AnnounceEndBlock("Done");

			if (vecOutputAuxCount[v] > 1) {
				AnnounceEndBlock("Done");
			}
		}
	}
	ncoutfile.close();

	AnnounceEndBlock("Done");
}

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Average over NetCDF files.
///	</summary>
void AverageOverNcFiles(
	const std::vector<std::string> & vecInputFilenames,
	const std::string & strOutputFile,
	const std::vector<std::string> & vecVariableNames,
	size_t sMemoryMax
) {
	_ASSERT(vecInputFilenames.size() > 0);
	_ASSERT(vecVariableNames.size() > 0);

	AnnounceStartBlock("Averaging over files");

	// Open output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputFile.c_str());
	}

	// Open input files
	std::vector<NcFile *> vecInputNcFiles;
	for (int f = 0; f < vecInputFilenames.size(); f++) {
		NcFile * pncinfile = new NcFile(vecInputFilenames[f].c_str(), NcFile::ReadOnly);
		if (!pncinfile->is_valid()) {
			_EXCEPTION1("Unable to open input datafile \"%s\"",
				vecInputFilenames[f].c_str());
		}

		vecInputNcFiles.push_back(pncinfile);
	}
	_ASSERT(vecInputNcFiles.size() == vecInputFilenames.size());

	// Loop through all variables
	for (int v = 0; v < vecVariableNames.size(); v++) {

		AnnounceStartBlock("Variable \"%s\"", vecVariableNames[v].c_str());

		// Variable dimension information
		std::vector<std::string> vecVariableDimNames;
		std::vector<long> vecVariableDimSizes;

		float dFillValueFloat = 0.0;
		double dFillValueDouble = 0.0;

		DataArray1D<float> dDataFloat;
		DataArray1D<double> dDataDouble;

		DataArray1D<float> dAccumDataFloat;
		DataArray1D<double> dAccumDataDouble;

		long lCountDims = (-1);
		size_t sAuxDimSize;
		DataArray1D<int> nDataCount;
		DataArray1D<int> nAccumDataCount;

		NcType nctype;

		// Number of iterations needed to average all data
		size_t sAuxCount = (-1);
		size_t sAuxSize = (-1);

		NcVar * varOut = NULL;

		// Loop through input files, initialize allocation size from file 0
		for (int f = 0; f < vecInputNcFiles.size(); f++) {

			NcVar * var = vecInputNcFiles[f]->get_var(vecVariableNames[v].c_str());
			if (var == NULL) {
				_EXCEPTION2("Input datafile \"%s\" does not contain variable \"%s\"",
					vecInputFilenames[f].c_str(),
					vecVariableNames[v].c_str());
			}

			if (var->num_dims() < 1) {
				_EXCEPTION2("Variable \"%s\" in file \"%s\" must have at least one dimension",
					vecVariableNames[v].c_str(),
					vecInputFilenames[f].c_str());
			}

			// Check for counts
			std::string strVariableCount = vecVariableNames[v] + "_count";
			NcVar * varCount = vecInputNcFiles[f]->get_var(strVariableCount.c_str());

			// Store variable dimension names and size, initialize variable in output file
			if (f == 0) {
				nctype = var->type();

				NcAtt * attFill = var->get_att("_FillValue");
				if (attFill != NULL) {
					if (attFill->type() != nctype) {
						_EXCEPTION1("Variable \"%s\" attribute \"_FillValue\" has incompatible type",
							vecVariableNames[v].c_str());
					}
					if (nctype == ncFloat) {
						dFillValueFloat = attFill->as_float(0);
					}
					if (nctype == ncDouble) {
						dFillValueDouble = attFill->as_double(0);
					}
				}

				CopyNcVar(*(vecInputNcFiles[f]), ncoutfile, vecVariableNames[v], true, false);

				varOut = ncoutfile.get_var(vecVariableNames[v].c_str());
				if (varOut == NULL) {
					_EXCEPTION3("Error copying variable \"%s\" from datafile \"%s\" to datafile \"%s\"",
						vecVariableNames[v].c_str(),
						vecInputFilenames[f].c_str(),
						strOutputFile.c_str());
				}

				size_t sTotalValues = 1;
				vecVariableDimNames.resize(var->num_dims());
				vecVariableDimSizes.resize(var->num_dims());
				for (int d = 0; d < var->num_dims(); d++) {
					vecVariableDimNames[d] = var->get_dim(d)->name();
					vecVariableDimSizes[d] = var->get_dim(d)->size();

					sTotalValues *= static_cast<size_t>(vecVariableDimSizes[d]);
				}
				sAuxDimSize = 1;
				for (int d = 1; d < var->num_dims(); d++) {
					sAuxDimSize *= static_cast<size_t>(vecVariableDimSizes[d]);
				}

				// Calculate number of iterations needed to stay within memory limit
				bool fSufficientMemory = false;
				if (nctype == ncFloat) {
					fSufficientMemory =
						CalculateIterationCountSize<float>(
							sMemoryMax,
							2,
							sTotalValues,
							sAuxCount,
							sAuxSize);

					dDataFloat.Allocate(sAuxSize);
					dAccumDataFloat.Allocate(sAuxSize);

				} else if (nctype == ncDouble) {
					fSufficientMemory =
						CalculateIterationCountSize<double>(
							sMemoryMax,
							2,
							sTotalValues,
							sAuxCount,
							sAuxSize);

					dDataDouble.Allocate(sAuxSize);
					dAccumDataDouble.Allocate(sAuxSize);

				} else {
					_EXCEPTION2("Variable \"%s\" has unsupported nctype (%i)",
						var->name(),
						var->type());
				}

				if (!fSufficientMemory) {
					_EXCEPTIONT("Insufficient memory to perform averaging");
				}

				// Allocate count array
				if (varCount == NULL) {
					nAccumDataCount.Allocate(1);

				} else {
					lCountDims = varCount->num_dims();

					if (varCount->get_dim(0)->size() != var->get_dim(0)->size()) {
						_EXCEPTION4("Variable \"%s\" not compatible with variable \"%s\": "
							"Lead dimension must have same size (%i found vs. %i expected)",
							vecVariableNames[v].c_str(),
							strVariableCount.c_str(),
							var->get_dim(0)->size(),
							varCount->get_dim(0)->size());
					}

					if (varCount->num_dims() == 1) {
						nDataCount.Allocate(varCount->get_dim(0)->size());
						nAccumDataCount.Allocate(varCount->get_dim(0)->size());

					} else if (varCount->num_dims() == var->num_dims()) {
						nDataCount.Allocate(sAuxSize);
						nAccumDataCount.Allocate(sAuxSize);

					} else {
						_EXCEPTION3("Variable \"%s\" not compatible with variable \"%s\": "
							"Either 1 dimension or %i dimensions expected",
							vecVariableNames[v].c_str(),
							strVariableCount.c_str(),
							var->num_dims());
					}
				}

				Announce("%lu values will be loaded per iteration over %lu iteration(s)",
					sAuxSize, sAuxCount);

			// Verify variable dimension names and size
			} else {
				if (var->num_dims() != vecVariableDimSizes.size()) {
					_EXCEPTION3("Input datafile \"%s\" dimension count mismatch. Found (%i) expected (%i)",
						vecInputFilenames[f].c_str(),
						var->num_dims(),
						vecVariableDimSizes.size());
				}
				for (int d = 0; d < var->num_dims(); d++) {
					if (var->get_dim(d)->name() != vecVariableDimNames[d]) {
						_EXCEPTION4("Input datafile \"%s\" dimension name mismatch in dimension %i. Found (%s) expected (%s)",
							vecInputFilenames[f].c_str(),
							d,
							var->get_dim(d)->name(),
							vecVariableDimNames[d].c_str());
					}
					if (var->get_dim(d)->size() != vecVariableDimSizes[d]) {
						_EXCEPTION4("Input datafile \"%s\" dimension size mismatch in dimension \"%s\". Found (%i) expected (%i)",
							vecInputFilenames[f].c_str(),
							var->get_dim(d)->name(),
							var->get_dim(d)->size(),
							vecVariableDimSizes[d]);
					}
				}

				// Verify count array
				if ((varCount == NULL) && (lCountDims >= 0)) {
					_EXCEPTION2("Variable \"%s\" missing from file \"%s\"",
						strVariableCount.c_str(),
						vecInputFilenames[f].c_str());
				}
				if ((varCount != NULL) && (lCountDims == (-1))) {
					_EXCEPTION2("Variable \"%s\" present in file \"%s\" but missing from earlier files",
						strVariableCount.c_str(),
						vecInputFilenames[f].c_str());
				}

				if (varCount != NULL) {
					if (varCount->get_dim(0)->size() != var->get_dim(0)->size()) {
						_EXCEPTION4("Variable \"%s\" not compatible with variable \"%s\": "
							"Lead dimension must have same size (%i found vs. %i expected)",
							vecVariableNames[v].c_str(),
							strVariableCount.c_str(),
							var->get_dim(0)->size(),
							varCount->get_dim(0)->size());
					}

					if (varCount->num_dims() != lCountDims) {
						_EXCEPTION2("Variable \"%s\" inconsistent across files: "
							"%i dimensions expected",
							strVariableCount.c_str(),
							var->num_dims());
					}
				}
			}
		}

		_ASSERT(sAuxCount != (-1));
		_ASSERT(sAuxSize != (-1));

		// Perform averaging
		AnnounceStartBlock("Performing averaging");

		size_t sCurrentArrayPos = 0;

		DataArray1D<long> lPos0;
		for (size_t c = 0; c < sAuxCount; c++) {
			if (sAuxCount > 1) {
				Announce("Iteration %lu/%lu", c+1, sAuxCount);
			}

			long long llCountReadSize = 0;
			long long llCurrentReadSize = 0;
			long long llCurrentWriteSize = 0;

			// Accumulate data
			if (nAccumDataCount.IsAttached()) {
				nAccumDataCount.Zero();
			}
			if (dAccumDataFloat.IsAttached()) {
				dAccumDataFloat.Zero();
			}
			if (dAccumDataDouble.IsAttached()) {
				dAccumDataDouble.Zero();
			}

			for (int f = 0; f < vecInputFilenames.size(); f++) {

				NcVar * var = vecInputNcFiles[f]->get_var(vecVariableNames[v].c_str());
				if (var == NULL) {
					_EXCEPTION2("Input datafile \"%s\" does not contain variable \"%s\"",
						vecInputFilenames[f].c_str(),
						vecVariableNames[v].c_str());
				}

				// Check and read counts
				std::string strVariableCount = vecVariableNames[v] + "_count";
				NcVar * varCount = vecInputNcFiles[f]->get_var(strVariableCount.c_str());

				if (varCount != NULL) {
					if (lCountDims == 1) {
						varCount->set_cur((long)0);
						varCount->get(&(nDataCount[0]), nDataCount.GetRows());

						for (size_t s = 0; s < nDataCount.GetRows(); s++) {
							nAccumDataCount[s] += nDataCount[s];
						}

					} else {
						llCountReadSize =
							NcGetPutSpecifiedDataSize<int>(
								varCount,
								lPos0,
								sAuxSize,
								c,
								nDataCount,
								NetCDF_Get);

						for (size_t s = 0; s < llCountReadSize; s++) {
							nAccumDataCount[s] += nDataCount[s];
						}
					}
				}

				// Read data
				if (nctype == ncFloat) {
					llCurrentReadSize =
						NcGetPutSpecifiedDataSize<float>(
							var,
							lPos0,
							sAuxSize,
							c,
							dDataFloat,
							NetCDF_Get);

					if (lCountDims == (-1)) {
						nAccumDataCount[0]++;
						for (size_t s = 0; s < llCurrentReadSize; s++) {
							dAccumDataFloat[s] += dDataFloat[s];
						}

					} else if (lCountDims == 1) {
						for (size_t s = 0; s < llCurrentReadSize; s++) {
							size_t sDim0Index = (sCurrentArrayPos + s) / sAuxDimSize;
							dAccumDataFloat[s] += dDataFloat[s] * static_cast<float>(nDataCount[sDim0Index]);
						}

					} else {
						_ASSERT(llCountReadSize == llCurrentReadSize);
						_ASSERT(dDataFloat.GetRows() == nAccumDataCount.GetRows());

						for (size_t s = 0; s < llCurrentReadSize; s++) {
							dAccumDataFloat[s] += dDataFloat[s] * static_cast<float>(nDataCount[s]);
						}
					}

				} else if (nctype == ncDouble) {
					llCurrentReadSize =
						NcGetPutSpecifiedDataSize<double>(
							var,
							lPos0,
							sAuxSize,
							c,
							dDataDouble,
							NetCDF_Get);

					if (lCountDims == (-1)) {
						nAccumDataCount[0]++;
						for (size_t s = 0; s < llCurrentReadSize; s++) {
							dAccumDataDouble[s] += dDataDouble[s];
						}

					} else if (lCountDims == 1) {
						for (size_t s = 0; s < llCurrentReadSize; s++) {
							size_t sDim0Index = (sCurrentArrayPos + s) / sAuxDimSize;
							dAccumDataDouble[s] += dDataDouble[s] * static_cast<float>(nAccumDataCount[sDim0Index]);
						}

					} else {
						_ASSERT(llCountReadSize == llCurrentReadSize);
						_ASSERT(dDataDouble.GetRows() == nAccumDataCount.GetRows());
						for (size_t s = 0; s < llCurrentReadSize; s++) {
							dAccumDataDouble[s] += dDataDouble[s] * static_cast<float>(nAccumDataCount[s]);
						}
					}
				}
			}

			// Write averaged data (float)
			if (nctype == ncFloat) {
				if (lCountDims == (-1)) {
					if (nAccumDataCount[0] == 0) {
						for (size_t s = 0; s < dAccumDataFloat.GetRows(); s++) {
							dAccumDataFloat[s] = dFillValueFloat;
						}
					} else {
						for (size_t s = 0; s < dAccumDataFloat.GetRows(); s++) {
							dAccumDataFloat[s] /= static_cast<float>(nAccumDataCount[0]);
						}
					}

				} else if (lCountDims == 1) {
					for (size_t s = 0; s < llCurrentReadSize; s++) {
						size_t sDim0Index = (sCurrentArrayPos + s) / sAuxDimSize;
						if (nAccumDataCount[sDim0Index] == 0) {
							dAccumDataFloat[s] = dFillValueFloat;
						} else {
							dAccumDataFloat[s] /= static_cast<float>(nAccumDataCount[sDim0Index]);
						}
					}

				} else {
					for (size_t s = 0; s < dAccumDataFloat.GetRows(); s++) {
						if (nAccumDataCount[s] == 0) {
							dAccumDataFloat[s] = dFillValueFloat;
						} else {
							dAccumDataFloat[s] /= static_cast<float>(nAccumDataCount[s]);
						}
					}
				}

				llCurrentWriteSize =
					NcGetPutSpecifiedDataSize<float>(
						varOut,
						lPos0,
						sAuxSize,
						c,
						dAccumDataFloat,
						NetCDF_Put);

			// Write averaged data (double)
			} else if (nctype == ncDouble) {
				if (lCountDims == (-1)) {
					if (nAccumDataCount[0] == 0) {
						for (size_t s = 0; s < dAccumDataDouble.GetRows(); s++) {
							dAccumDataDouble[s] = dFillValueDouble;
						}
					} else {
						for (size_t s = 0; s < dAccumDataDouble.GetRows(); s++) {
							dAccumDataDouble[s] /= static_cast<float>(nAccumDataCount[0]);
						}
					}

				} else if (lCountDims == 1) {
					for (size_t s = 0; s < llCurrentReadSize; s++) {
						size_t sDim0Index = (sCurrentArrayPos + s) / sAuxDimSize;
						if (nAccumDataCount[sDim0Index] == 0) {
							dAccumDataDouble[s] = dFillValueDouble;
						} else {
							dAccumDataDouble[s] /= static_cast<float>(nAccumDataCount[sDim0Index]);
						}
					}

				} else {
					for (size_t s = 0; s < dAccumDataDouble.GetRows(); s++) {
						if (nAccumDataCount[s] == 0) {
							dAccumDataDouble[s] = dFillValueDouble;
						} else {
							dAccumDataDouble[s] /= static_cast<float>(nAccumDataCount[s]);
						}
					}
				}

				llCurrentWriteSize =
					NcGetPutSpecifiedDataSize<double>(
						varOut,
						lPos0,
						sAuxSize,
						c,
						dAccumDataDouble,
						NetCDF_Put);
			}

			if (llCurrentReadSize != llCurrentWriteSize) {
				_EXCEPTION2("Read/write size mismatch: %ld / %ld",
					llCurrentReadSize,
					llCurrentWriteSize);
			}

			// Update 1D array position
			sCurrentArrayPos += static_cast<size_t>(llCurrentReadSize);
		}

		AnnounceEndBlock("Done");
		AnnounceEndBlock(NULL);
	}

	for (int f = 0; f < vecInputNcFiles.size(); f++) {
		vecInputNcFiles[f]->close();
		delete vecInputNcFiles[f];
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
	std::string strClimoPeriod;

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
		CommandLineStringD(strClimoPeriod, "period", "daily", "[daily|monthly|seasonal|annual]");
		CommandLineStringD(strClimoType, "type", "mean", "[mean|meansq]");
		CommandLineBool(fIncludeLeapDays, "include_leap_days");
		//CommandLineInt(nFourierModes, "time_modes", 0);
		CommandLineBool(fMissingData, "missingdata");
		CommandLineString(strTempFilePath, "temp_file_path", "/tmp");
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

	// Climatology period
	ClimatologyPeriod eClimoPeriod;
	STLStringHelper::ToLower(strClimoPeriod);
	if (strClimoPeriod == "daily") {
		eClimoPeriod = ClimatologyPeriod_Daily;
	} else if (strClimoPeriod == "monthly") {
		eClimoPeriod = ClimatologyPeriod_Monthly;
	} else if (strClimoPeriod == "seasonal") {
		eClimoPeriod = ClimatologyPeriod_Seasonal;
	} else if (strClimoPeriod == "annual") {
		eClimoPeriod = ClimatologyPeriod_Annual;
	} else {
		_EXCEPTIONT("--period invalid; expected \"daily\", \"monthly\", \"seasonal\" or \"annual\"");
	}

	// Climatology type
	ClimatologyType eClimoType;
	STLStringHelper::ToLower(strClimoType);
	if (strClimoType == "mean") {
		eClimoType = ClimatologyType_Mean;
	} else if (strClimoType == "meansq") {
		eClimoType = ClimatologyType_MeanSq;
	} else {
		_EXCEPTIONT("--type invalid; expected \"mean\" or \"meansq\"");
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
	std::vector< std::vector<std::string> > vecVariableSpecifiedDims;

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
		eClimoPeriod,
		eClimoType,
		nFourierModes,
		fMissingData,
		fOutputSliceCounts,
		fVerbose
	);

#if defined(TEMPEST_MPIOMP)
	if (nMPISize != 1) {

		// Barrier
		MPI_Barrier(MPI_COMM_WORLD);
		AnnounceSetOutputBuffer(stdout);
		AnnounceOnlyOutputOnRankZero();

		// Perform averaging
		if (nMPIRank == 0) {
			AnnounceBanner();

			if (nMPIRank == 0) {
				std::vector< std::string > vecAllOutputFileList;
				for (int f = 0; f < nMPISize; f++) {
					char szBuffer[16];
					sprintf(szBuffer, "%06i", f);

					std::string strTempOutputFile =
						strTempFilePath + "/tempClimatology" + szBuffer + ".nc";

					vecAllOutputFileList.push_back(strTempOutputFile);
				}

				std::vector<std::string> vecTempVariableNames;

				std::string strVariablePrefix;
				if (eClimoPeriod == ClimatologyPeriod_Daily) {
					strVariablePrefix = "daily";
				} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
					strVariablePrefix = "monthly";
				} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
					strVariablePrefix = "seasonal";
				} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
					strVariablePrefix = "annual";
				}
				if (eClimoType == ClimatologyType_Mean) {
					strVariablePrefix += "mean_";
				} else if (eClimoType == ClimatologyType_MeanSq) {
					strVariablePrefix += "meansq_";
				}

				for (int v = 0; v < vecVariableNames.size(); v++) {
					vecTempVariableNames.push_back(strVariablePrefix + vecVariableNames[v]);
				}

				AverageOverNcFiles(
					vecAllOutputFileList,
					strOutputFile,
					vecTempVariableNames,
					sMemoryMax);
			}
		}
	}
#endif

	AnnounceBanner();

} catch(Exception & e) {
	AnnounceOutputOnAllRanks();
	AnnounceSetOutputBuffer(stdout);
	Announce(e.ToString().c_str());

#if defined(TEMPEST_MPIOMP)
	MPI_Abort(MPI_COMM_WORLD, -1);
#endif
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

