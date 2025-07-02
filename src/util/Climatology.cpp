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
#include "TimeMatch.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "FourierTransforms.h"
#include "Units.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <ctime>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Possible climatology periods that can be calculated with Climatology().
///	</summary>
enum ClimatologyPeriod {
	ClimatologyPeriod_Daily,
	ClimatologyPeriod_Monthly,
	ClimatologyPeriod_Seasonal,
	ClimatologyPeriod_Annual,
	ClimatologyPeriod_All
};

///	<summary>
///		Possible climatology types that can be calculated with Climatology().
///	</summary>
enum ClimatologyType {
	ClimatologyType_Mean,
	ClimatologyType_MeanSq,
	ClimatologyType_StdDev,
	ClimatologyType_AutoCor,
	ClimatologyType_Min,
	ClimatologyType_Max,
	ClimatologyType_AvgMin,
	ClimatologyType_AvgMax,
	ClimatologyType_AvgCount,
	ClimatologyType_ThreshSum,
	ClimatologyType_TimeUntil,
	ClimatologyType_MaxConsec
};

///////////////////////////////////////////////////////////////////////////////

int GetClimatologyTimeIndex(
	const Time & time,
	ClimatologyPeriod eClimoPeriod
) {
	int iTimeIndex;
	if (eClimoPeriod == ClimatologyPeriod_Daily) {
		Time timeJan01 = time;
		timeJan01.SetMonth(1);
		timeJan01.SetDay(1);

		iTimeIndex = time.DayNumber() - timeJan01.DayNumber();
		if ((time.IsLeapYear()) && (iTimeIndex >= 60)) {
			iTimeIndex--;
		}

		_ASSERT((iTimeIndex >= 0) && (iTimeIndex < 365));

	} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
		iTimeIndex = time.GetZeroIndexedMonth();
		_ASSERT((iTimeIndex >= 0) && (iTimeIndex < 12));

	} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
		iTimeIndex = (int)(time.GetSeasonIndex());
		_ASSERT((iTimeIndex >= 0) && (iTimeIndex < 4));

	} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
		iTimeIndex = 0;

	} else if (eClimoPeriod == ClimatologyPeriod_All) {
		iTimeIndex = 0;

	} else {
		_EXCEPTIONT("Invalid eClimoPeriod");
	}

	return iTimeIndex;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Climatology threshold types.
///	</summary>
enum ThresholdType {
	ThresholdType_GreaterThan,
	ThresholdType_GreaterThanOrEqualTo,
	ThresholdType_LessThan,
	ThresholdType_LessThanOrEqualTo,
	ThresholdType_EqualTo,
	ThresholdType_NotEqualTo
};

///	<summary>
///		Parse a threshold string.
///	</summary>
void ParseClimatologyThreshold(
	std::string strThreshold,
	ThresholdType & eThreshType,
	float & dThreshValue,
	std::string & strThreshUnits
) {
	STLStringHelper::RemoveWhitespaceInPlace(strThreshold);
	if (strThreshold.length() < 2) {
		_EXCEPTION1("Invalid threshold string \"%s\": expected form \"<op><value>\"",
			strThreshold.c_str());
	}
	if (strThreshold[0] == '>') {
		if (strThreshold[1] == '=') {
			eThreshType = ThresholdType_GreaterThanOrEqualTo;
			strThreshold = strThreshold.substr(2);
		} else {
			eThreshType = ThresholdType_GreaterThan;
			strThreshold = strThreshold.substr(1);
		}

	} else if (strThreshold[0] == '<') {
		if (strThreshold[1] == '=') {
			eThreshType = ThresholdType_LessThanOrEqualTo;
			strThreshold = strThreshold.substr(2);
		} else {
			eThreshType = ThresholdType_LessThan;
			strThreshold = strThreshold.substr(1);
		}

	} else if (strThreshold[0] == '=') {
		eThreshType = ThresholdType_EqualTo;
		strThreshold = strThreshold.substr(1);

	} else if ((strThreshold[0] == '!') && (strThreshold[1] == '=')) {
		eThreshType = ThresholdType_NotEqualTo;
		strThreshold = strThreshold.substr(1);

	} else {
		_EXCEPTION1("Invalid threshold operator in string \"%s\": operator must be one of >,>=,<,<=,=,!=",
			strThreshold.c_str());
	}

	std::string strThreshValue;
	SplitIntoValueAndUnits(strThreshold, strThreshValue, strThreshUnits);
	dThreshValue = std::stof(strThreshValue);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A structure holding requested climatology info.
///	</summary>
struct ClimatologyInfo {
public:
	ClimatologyInfo() :
		type(ClimatologyType_Mean),
		lag(0),
		threshtype(ThresholdType_GreaterThan),
		threshvalue(0.0),
		threshunits(""),
		threshorigunits(""),
		threshstring("")
	{ }

public:
	ClimatologyType type;
	int lag;
	ThresholdType threshtype;
	float threshvalue;
	std::string threshunits;
	std::string threshorigunits;
	std::string threshstring;
};

///	<summary>
///		Parse a list of climatology types.
///	</summary>
void ParseClimatologyTypes(
	const std::string & strClimatologyTypes,
	std::vector<ClimatologyInfo> & vecClimatologyInfo
) {
	vecClimatologyInfo.clear();

	// Break up into individual types
	std::vector<std::string> vecParsedClimatologyList;
	STLStringHelper::ParseVariableList(strClimatologyTypes, vecParsedClimatologyList);

	for (int i = 0; i < vecParsedClimatologyList.size(); i++) {
		STLStringHelper::RemoveWhitespaceInPlace(vecParsedClimatologyList[i]);

		for (int j = 0; j < vecParsedClimatologyList[i].length(); j++) {
			if (std::isalpha(vecParsedClimatologyList[i][j])) {
				vecParsedClimatologyList[i][j] = std::tolower(vecParsedClimatologyList[i][j]);
			} else {
				break;
			}
		}

		const std::string & strClimoType = vecParsedClimatologyList[i];

		ClimatologyInfo climoinfo;
		if (vecParsedClimatologyList[i] == "mean") {
			climoinfo.type = ClimatologyType_Mean;
		} else if (strClimoType == "meansq") {
			climoinfo.type = ClimatologyType_MeanSq;
		} else if (strClimoType == "stddev") {
			climoinfo.type = ClimatologyType_StdDev;
		} else if (strClimoType == "min") {
			climoinfo.type = ClimatologyType_Min;
		} else if (strClimoType == "max") {
			climoinfo.type = ClimatologyType_Max;
		} else if (strClimoType == "avgmin") {
			climoinfo.type = ClimatologyType_AvgMin;
		} else if (strClimoType == "avgmax") {
			climoinfo.type = ClimatologyType_AvgMax;
		} else if (strClimoType == "avgmin") {
			climoinfo.type = ClimatologyType_AvgMin;
		} else if ((strClimoType.length() >= 9) && (strClimoType.substr(0,7) == "autocor")) {
			climoinfo.type = ClimatologyType_AutoCor;
			if (strClimoType[7] != ':') {
				_EXCEPTIONT("--type invalid; expected \"autocor:<lag>\"");
			}
			std::string strAutoCorLag = strClimoType.substr(8);
			if (!STLStringHelper::IsInteger(strAutoCorLag)) {
				_EXCEPTIONT("--type invalid; expected \"autocor:<lag>\"");
			}
			climoinfo.lag = std::stoi(strAutoCorLag);
			if (climoinfo.lag < 1) {
				_EXCEPTIONT("--type \"autocor:<lag>\" must provide a non-negative lag time");
			}

		} else if ((strClimoType.length() >= 10) && (strClimoType.substr(0,8) == "avgcount")) {
			climoinfo.type = ClimatologyType_AvgCount;
			climoinfo.threshstring = strClimoType.substr(8);
			ParseClimatologyThreshold(
				climoinfo.threshstring,
				climoinfo.threshtype,
				climoinfo.threshvalue,
				climoinfo.threshunits);

		} else if ((strClimoType.length() >= 11) && (strClimoType.substr(0,9) == "threshsum")) {
			climoinfo.type = ClimatologyType_ThreshSum;
			climoinfo.threshstring = strClimoType.substr(9);
			ParseClimatologyThreshold(
				climoinfo.threshstring,
				climoinfo.threshtype,
				climoinfo.threshvalue,
				climoinfo.threshunits);

			if ((climoinfo.threshtype != ThresholdType_GreaterThan) &&
			    (climoinfo.threshtype != ThresholdType_LessThan)
			) {
				_EXCEPTIONT("--type \"threshsum\" can only be combined with > or < operators");
			}

		} else if ((strClimoType.length() >= 11) && (strClimoType.substr(0,9) == "timeuntil")) {
			climoinfo.type = ClimatologyType_TimeUntil;
			climoinfo.threshstring = strClimoType.substr(9);
			ParseClimatologyThreshold(
				climoinfo.threshstring,
				climoinfo.threshtype,
				climoinfo.threshvalue,
				climoinfo.threshunits);

		} else if ((strClimoType.length() >= 11) && (strClimoType.substr(0,9) == "maxconsec")) {
			climoinfo.type = ClimatologyType_MaxConsec;
			climoinfo.threshstring = strClimoType.substr(9);
			ParseClimatologyThreshold(
				climoinfo.threshstring,
				climoinfo.threshtype,
				climoinfo.threshvalue,
				climoinfo.threshunits);

		} else {
			_EXCEPTIONT("--type invalid; expected \"mean\", \"meansq\", \"stddev\", \"autocor:<lag>\", \"min\", \"max\", \"avgmin\", \"avgmax\", \"avgcount<op><value>\", \"threshsum<op><value>\", \"timeuntil<op><value>\", \"maxconsec<op><value>\"");
		}

		climoinfo.threshorigunits = climoinfo.threshunits;
		vecClimatologyInfo.push_back(climoinfo);
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string ClimoVariablePrefix(
	ClimatologyPeriod eClimoPeriod,
	const ClimatologyInfo & aClimoInfo
) {
	ClimatologyType eClimoType = aClimoInfo.type;

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
	} else if (eClimoType == ClimatologyType_StdDev) {
		strVariablePrefix += "stddev_";
	} else if (eClimoType == ClimatologyType_AutoCor) {
		
		if (eClimoPeriod == ClimatologyPeriod_All) {
			strVariablePrefix += "autocor" + std::to_string((long)aClimoInfo.lag) + "_";
		} else {
			strVariablePrefix += "avgautocor" + std::to_string((long)aClimoInfo.lag) + "_";
		}
	} else if (eClimoType == ClimatologyType_Min) {
		strVariablePrefix += "min_";
	} else if (eClimoType == ClimatologyType_Max) {
		strVariablePrefix += "max_";
	} else if (eClimoType == ClimatologyType_AvgMin) {
		strVariablePrefix += "avgmin_";
	} else if (eClimoType == ClimatologyType_AvgMax) {
		strVariablePrefix += "avgmax_";
	} else if ((eClimoType == ClimatologyType_AvgCount) ||
	           (eClimoType == ClimatologyType_ThreshSum) ||
	           (eClimoType == ClimatologyType_TimeUntil) ||
	           (eClimoType == ClimatologyType_MaxConsec)
	) {
		if (eClimoType == ClimatologyType_AvgCount) {
			strVariablePrefix += "avgcount_";
		} else if (eClimoType == ClimatologyType_ThreshSum) {
			strVariablePrefix += "threshsum_";
		} else if (eClimoType == ClimatologyType_TimeUntil) {
			strVariablePrefix += "timeuntil_";
		} else {
			strVariablePrefix += "maxconsec_";
		}

		if (aClimoInfo.threshtype == ThresholdType_GreaterThan) {
			strVariablePrefix += "gt_";
		} else if (aClimoInfo.threshtype == ThresholdType_GreaterThanOrEqualTo) {
			strVariablePrefix += "gte_";
		} else if (aClimoInfo.threshtype == ThresholdType_LessThan) {
			strVariablePrefix += "lt_";
		} else if (aClimoInfo.threshtype == ThresholdType_LessThanOrEqualTo) {
			strVariablePrefix += "lte_";
		} else if (aClimoInfo.threshtype == ThresholdType_EqualTo) {
			strVariablePrefix += "eq_";
		} else if (aClimoInfo.threshtype == ThresholdType_NotEqualTo) {
			strVariablePrefix += "neq_";
		}

		char szThresholdValueBuffer[15];
		snprintf(szThresholdValueBuffer, 15, "%g", aClimoInfo.threshvalue);
		for (int i = 0; i < 15; i++) {
			if (szThresholdValueBuffer[i] == '.') {
				szThresholdValueBuffer[i] = 'p';
			}
		}

		strVariablePrefix += std::string(szThresholdValueBuffer) + "_";
	}
	return strVariablePrefix;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Produce a climatology over the given input files.
///	</summary>
void Climatology(
	const std::vector<std::string> & vecInputFileList,
	const std::string & strConnectivityFile,
	const std::string & strOutputFile,
	const std::string & strVariables,
	const std::string & strOutputVariables,
	bool fIncludeLeapDays,
	const std::string & strStartTime,
	const std::string & strEndTime,
	ClimatologyPeriod eClimoPeriod,
	std::vector<ClimatologyInfo> vecClimoInfo,
	bool fMissingData,
	const std::string & strFillValueOverride,
	const std::string & strTimeFilter,
	std::string strLatitudeName,
	std::string strLongitudeName,
	bool fDiagonalConnectivity,
	bool fRegional,
	bool fOutputSliceCounts,
	bool fVerbose
) {
	_ASSERT(vecInputFileList.size() > 0);
	_ASSERT(vecClimoInfo.size() > 0);

	// Climatology period
	if ((eClimoPeriod != ClimatologyPeriod_Daily) &&
		(eClimoPeriod != ClimatologyPeriod_Monthly) &&
		(eClimoPeriod != ClimatologyPeriod_Seasonal) &&
		(eClimoPeriod != ClimatologyPeriod_Annual) &&
		(eClimoPeriod != ClimatologyPeriod_All)
	) {
		_EXCEPTIONT("Invalid eClimoPeriod");
	}

	// Check vecClimoInfo input
	bool fPeriodAccumulationNeeded = false;
	for (int ct = 0; ct < vecClimoInfo.size(); ct++) {
		// Climatology type
		if ((vecClimoInfo[ct].type != ClimatologyType_Mean) &&
		    (vecClimoInfo[ct].type != ClimatologyType_MeanSq) &&
		    (vecClimoInfo[ct].type != ClimatologyType_StdDev) &&
			(vecClimoInfo[ct].type != ClimatologyType_Min) &&
			(vecClimoInfo[ct].type != ClimatologyType_Max) &&
			(vecClimoInfo[ct].type != ClimatologyType_AvgMin) &&
			(vecClimoInfo[ct].type != ClimatologyType_AvgMax) &&
	   		(vecClimoInfo[ct].type != ClimatologyType_AutoCor) &&
			(vecClimoInfo[ct].type != ClimatologyType_AvgCount) &&
			(vecClimoInfo[ct].type != ClimatologyType_ThreshSum) &&
			(vecClimoInfo[ct].type != ClimatologyType_TimeUntil) &&
			(vecClimoInfo[ct].type != ClimatologyType_MaxConsec)
		) {
			_EXCEPTIONT("Invalid eClimoType");
		}

		// Autocorrelation only available without missing data
		if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
			if (fMissingData) {
				_EXCEPTIONT("At present, --type \"autocor\" cannot be used with --missingdata");
			}
			if (vecClimoInfo[ct].lag < 1) {
				_EXCEPTIONT("Nonnegative lag required with climatology type \"autocor\"");
			}
		}

		// AvgMin and AvgMax with ClimatologyPeriod_All is the same as Min and Max
		if (eClimoPeriod == ClimatologyPeriod_All) {
			if (vecClimoInfo[ct].type == ClimatologyType_AvgMin) {
				vecClimoInfo[ct].type = ClimatologyType_Min;
			}
			if (vecClimoInfo[ct].type == ClimatologyType_AvgMax) {
				vecClimoInfo[ct].type = ClimatologyType_Max;
			}
		}

		// Min and Max with ClimatologyPeriod_Annual is invalid; override to AvgMin and AvgMax
		if (eClimoPeriod == ClimatologyPeriod_Annual) {
			if (vecClimoInfo[ct].type == ClimatologyType_Min) {
				vecClimoInfo[ct].type = ClimatologyType_AvgMin;
				std::cout << "WARNING: For --period \"annual\" --type \"min\" is invalid; setting --type \"avgmin\"" << std::endl;
			}
			if (vecClimoInfo[ct].type == ClimatologyType_Max) {
				vecClimoInfo[ct].type = ClimatologyType_AvgMax;
				std::cout << "WARNING: For --period \"annual\" --type \"max\" is invalid; setting --type \"avgmax\"" << std::endl;
			}
		}

		// Check if this climatology type requires data to be accumulated after each period
		if ((vecClimoInfo[ct].type == ClimatologyType_AvgMin) ||
		    (vecClimoInfo[ct].type == ClimatologyType_AvgMax) ||
		    (vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
		    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) ||
		    (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) ||
		    (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) ||
		    ((vecClimoInfo[ct].type == ClimatologyType_AutoCor) && (eClimoPeriod != ClimatologyPeriod_All))
		) {
			fPeriodAccumulationNeeded = true;
		}
	}

	// Time units
	std::string strTimeUnits = "days since 0001-01-01";

	// Create VariableRegistry
	VariableRegistry varreg;

	// Parse --var argument
	std::vector<std::string> vecVariableNames;
	std::vector<VariableIndex> vecVarIxIn;
	if (strVariables != "") {
		std::string strVariablesTemp = strVariables;
		for (;;) {
			VariableIndex varix;
			int iLast = varreg.FindOrRegisterSubStr(strVariablesTemp, &varix) + 1;

			vecVarIxIn.push_back(varix);
			vecVariableNames.push_back( strVariablesTemp.substr(0,iLast-1) );
			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}
	}
 
	// Parse --varout argument
	std::vector<std::string> vecVariableOutNames;
	if (strOutputVariables != "") {
		STLStringHelper::ParseVariableList(strOutputVariables, vecVariableOutNames);
		if (vecVariableOutNames.size() != vecVariableNames.size()) {
			_EXCEPTION2("Inconsistent number of variables in --var (%lu) and --varout (%lu)",
				vecVariableNames.size(), vecVariableOutNames.size());
		}
	}

	// Open output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputFile.c_str());
	}

	// Flag indicating output file has been initialized
	bool fOutputInitialized = false;

	// Parse timefilter
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

		std::string strTimeFilterTemp = strTimeFilter;
		if (strTimeFilterTemp == "3hr") {
			strTimeFilterTemp = "....-..-.. (00|03|06|09|12|15|18|21):00:00";
		}
		if (strTimeFilterTemp == "6hr") {
			strTimeFilterTemp = "....-..-.. (00|06|12|18):00:00";
		}
		if (strTimeFilterTemp == "daily") {
			strTimeFilterTemp = "....-..-.. 00:00:00";
		}

		try {
			reTimeSubset.assign(strTimeFilterTemp);
		} catch(std::regex_error & reerr) {
			_EXCEPTION2("Parse error in --timefilter regular expression \"%s\" (code %i)",
				strTimeFilterTemp.c_str(), reerr.code());
		}
	}
#endif
	// Vector of pointers to output NcVars for climatology
	std::vector<NcVar*> vecNcVarOut;
	vecNcVarOut.resize(vecClimoInfo.size() * vecVariableNames.size());

	// Vector of pointers to output NcVars for slice counts
	std::vector<NcVar*> vecNcVarCount;
	vecNcVarCount.resize(vecVariableNames.size());

	// Time information
	Time::CalendarType caltype;
	std::string strTimeCalendar;
	size_t sOutputTimes = 0;

	Time timeStartTime;
	Time timeEndTime;

	// Vector of output NcDim * for each variable
	std::vector< std::vector<NcDim *> > vecOutputNcDim;
	vecOutputNcDim.resize(vecVariableNames.size());

	// SimpleGrid
	SimpleGrid grid;

	// Grid dimensions in output file
	NcDim * dimGrid0 = NULL;
	NcDim * dimGrid1 = NULL;

	/////////////////////////////////////////////
	// Initialization
	{
		// Initialize the NcFileVector
		NcFileVector vecncDataFiles;
		vecncDataFiles.ParseFromString(vecInputFileList[0].c_str());
		_ASSERT(vecncDataFiles.size() > 0);

		NcFile & ncinfile = *(vecncDataFiles[0]);

		// Initialize the SimpleGrid

		// Check for connectivity file
		if (strConnectivityFile != "") {
			AnnounceStartBlock("Generating grid information from connectivity file");
			grid.FromFile(strConnectivityFile);
			AnnounceEndBlock("Done");

		// No connectivity file; check for latitude/longitude dimension
		} else {
			AnnounceStartBlock("No connectivity file specified");
			Announce("Attempting to generate latitude-longitude grid from data file");

			if (vecInputFileList.size() < 1) {
				_EXCEPTIONT("No data files specified; unable to generate grid");
			}

			grid.GenerateLatitudeLongitude(
				vecncDataFiles[0],
				strLatitudeName,
				strLongitudeName,
				fRegional,
				fDiagonalConnectivity);

			if (grid.m_nGridDim.size() != 2) {
				_EXCEPTIONT("Logic error when generating connectivity");
			}
			AnnounceEndBlock("Done");
		}

		// Write grid information to output file
		if (grid.m_nGridDim.size() == 1) {
			dimGrid0 = ncoutfile.add_dim("ncol", grid.m_nGridDim[0]);
		} else if (grid.m_nGridDim.size() == 2) {
			dimGrid0 = ncoutfile.add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]);
			dimGrid1 = ncoutfile.add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]);

			CopyNcVarIfExists(
				ncinfile,
				ncoutfile,
				strLatitudeName.c_str(),
				true,
				true);

			CopyNcVarIfExists(
				ncinfile,
				ncoutfile,
				strLongitudeName.c_str(),
				true,
				true);

		} else {
			_EXCEPTIONT("Only 1D or 2D spatial data supported");
		}

		// Get information from the time variable
		AnnounceStartBlock("Loading time information from input datafile");
		NcVar * varInTime = NcGetTimeVariable(ncinfile);
		if (varInTime == NULL) {
			_EXCEPTION1("File \"%s\" does not contain variable \"time\"",
				vecncDataFiles.GetFilename(0).c_str());
		}
		if (varInTime->num_dims() != 1) {
			_EXCEPTION1("Variable \"time\" in NetCDF file \"%s\" has more than 1 dimension",
				vecncDataFiles.GetFilename(0).c_str());
		}

		{
			NcAtt * attCalendar = varInTime->get_att("calendar");
			if (attCalendar == NULL) {
				Announce("Dataset \"time\" does not specify \"calendar\" attribute.  Assuming \"standard\".");
				strTimeCalendar = "standard";
			} else {
				strTimeCalendar = attCalendar->as_string(0);
			}

			caltype = Time::CalendarTypeFromString(strTimeCalendar);
			if (caltype == Time::CalendarUnknown) {
				_EXCEPTION1("Unknown calendar name \"%s\"; cannot determine number of days per year", strTimeCalendar.c_str());

			} else if (caltype == Time::Calendar360Day) {
				Announce("Calendar \"%s\" contains 360 days per year (30 days per month)",
					strTimeCalendar.c_str());
				sOutputTimes = 360;

			} else {
				Announce("Calendar \"%s\" contains 365 days per year",
					strTimeCalendar.c_str());
				sOutputTimes = 365;
			}

			// Start and end time
			if (strStartTime != "") {
				timeStartTime = Time(caltype);
				timeStartTime.FromFormattedString(strStartTime);
			}
			if (strEndTime != "") {
				timeEndTime = Time(caltype);
				timeEndTime.FromFormattedString(strEndTime);
			}

			// Other climatologies
			if (eClimoPeriod == ClimatologyPeriod_Monthly) {
				sOutputTimes = 12;
			} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
				sOutputTimes = 4;
			} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
				sOutputTimes = 1;
			} else if (eClimoPeriod == ClimatologyPeriod_All) {
				sOutputTimes = 1;
			}
		}
		AnnounceEndBlock("Done");

		// Write time dimension to output file
		AnnounceStartBlock("Initializing output file");

		Announce("Initialize time dimension and variable");

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

		} else if (eClimoPeriod == ClimatologyPeriod_All) {
			Time timeOut("0001-01-01-00000-000000", Time::CalendarNoLeap);
			vecTimes.push_back(timeOut);
		}

		WriteCFTimeDataToNcFile(
			&ncoutfile,
			strOutputFile,
			vecTimes,
			true);

		NcDim * dimOutTime = NcGetTimeDimension(ncoutfile);
		_ASSERT(dimOutTime != NULL);

		NcVar * varOutTime = NcGetTimeVariable(ncoutfile);
		_ASSERT(varOutTime != NULL);

		// Generate output attribute for time dimension
		std::string strTimeTypeAtt;
		if (eClimoPeriod == ClimatologyPeriod_Daily) {
			strTimeTypeAtt = "daily ";
		} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
			strTimeTypeAtt = "monthly ";
		} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
			strTimeTypeAtt = "seasonal ";
		} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
			strTimeTypeAtt = "annual ";
		}

		strTimeTypeAtt += "climatology";

		varOutTime->add_att("type", strTimeTypeAtt.c_str());

		// Create output variables
		Announce("Creating output variables");
		for (int v = 0; v < vecVariableNames.size(); v++) {

			// Get auxiliary dimension info and verify consistency
			DimInfoVector vecAuxDimInfo;

			varreg.AppendVariableToProcessingQueue(
				vecncDataFiles,
				grid,
				vecVarIxIn[v],
				&vecAuxDimInfo);

			// Copy auxiliary dimension variables from input to output
			vecOutputNcDim[v].resize(vecAuxDimInfo.size() + 1);
			vecOutputNcDim[v][0] = dimOutTime;
			for (int d = 0; d < vecAuxDimInfo.size(); d++) {
				vecOutputNcDim[v][d] =
					ncoutfile.get_dim(vecAuxDimInfo[d].name.c_str());
				if (vecOutputNcDim[v][d] == NULL) {
					vecOutputNcDim[v][d] = ncoutfile.add_dim(
						vecAuxDimInfo[d].name.c_str(),
						vecAuxDimInfo[d].size);

					if (vecOutputNcDim[v][d] == NULL) {
						_EXCEPTION2("Error adding dimension \"%s\" (%li) to output file",
							vecAuxDimInfo[d].name.c_str(),
							vecAuxDimInfo[d].size);
					}

				} else {
					if (vecOutputNcDim[v][d]->size() != vecAuxDimInfo[d].size) {
						std::string strVarName =
							varreg.GetVariableString(vecVarIxIn[v]);

						_EXCEPTION4("Dimension size mismatch when initializing variable \"%s\": Expected dimension \"%s\" to have size \"%li\" (found \"%li\")",
							strVarName.c_str(),
							vecAuxDimInfo[d].name.c_str(),
							vecAuxDimInfo[d].size,
							vecOutputNcDim[v][d]->size());
					}
				}
			}
			if (grid.m_nGridDim.size() == 1) {
				vecOutputNcDim[v].push_back(dimGrid0);
			} else {
				vecOutputNcDim[v].push_back(dimGrid0);
				vecOutputNcDim[v].push_back(dimGrid1);
			}

			// One output variable for each climatology type
			for (int ct = 0; ct < vecClimoInfo.size(); ct++) {

				// Create output variable
				std::string strOutVar;
				if (vecVariableOutNames.size() == 0) {
		   			strOutVar =
						ClimoVariablePrefix(eClimoPeriod, vecClimoInfo[ct])
						+ vecVariableNames[v];
				} else {
					_ASSERT(vecClimoInfo.size() == 1);
					_ASSERT(v < vecVariableOutNames.size());
					strOutVar =
						ClimoVariablePrefix(eClimoPeriod, vecClimoInfo[ct])
						+ vecVariableOutNames[v];
				}

				NcVar * varOut =
					ncoutfile.add_var(
						strOutVar.c_str(),
						ncFloat,
						vecOutputNcDim[v].size(),
						const_cast<const NcDim**>(&(vecOutputNcDim[v][0])));

				if (varOut == NULL) {
					_EXCEPTION1("Unable to create output variable \"%s\" in output file",
						strOutVar.c_str());
				}
/*
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
*/

				// Copy attributes from input variable to output variable
				NcVar * varIn = ncinfile.get_var(vecVariableNames[v].c_str());
				if (varIn != NULL) {
					std::vector<std::string> vecDoNotCopyNames;
					if ((vecClimoInfo[ct].type == ClimatologyType_MeanSq) ||
					    (vecClimoInfo[ct].type == ClimatologyType_AutoCor) ||
					    (vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
					    (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) ||
					    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) ||
					    (vecClimoInfo[ct].type == ClimatologyType_MaxConsec)
					) {
						vecDoNotCopyNames.push_back("units");
						vecDoNotCopyNames.push_back("standard_name");
						vecDoNotCopyNames.push_back("long_name");
						vecDoNotCopyNames.push_back("original_units");

						std::string strUnits;
						NcAtt * attUnits = varIn->get_att("units");
						if (attUnits != NULL) {
							strUnits = attUnits->as_string(0);
						} else {
							strUnits = "unknown";
						}
						std::string strStandardName;
						NcAtt * attStandardName = varIn->get_att("standard_name");
						if (attStandardName != NULL) {
							strStandardName = attStandardName->as_string(0);
						} else {
							strStandardName = "variable";
						}
						std::string strLongName;
						NcAtt * attLongName = varIn->get_att("long_name");
						if (attLongName != NULL) {
							strLongName = attLongName->as_string(0);
						} else {
							strLongName = "variable";
						}

						varOut->add_att("original_units", strUnits.c_str());

						if (vecClimoInfo[ct].type == ClimatologyType_MeanSq) {
							strUnits = std::string("(") + strUnits + std::string(")2");
							strStandardName = strStandardName + std::string(" squared");
							strLongName = strLongName + std::string(" squared");

						} else if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
							strUnits = "unitless";
							strStandardName = std::string("autocorrelation of ") + strStandardName;
							strLongName = std::string("autocorrelation of ") + strLongName;

						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgCount) {
							strUnits = "times";
							strStandardName = "average count of " + strStandardName + " satisfying threshold";
							strLongName = "average count of " + strLongName + " satisfying threshold";

						} else if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
							strUnits = "times";
							strStandardName = "average time until " + strStandardName + " satisfies threshold";
							strLongName = "average time until " + strLongName + " satisfies threshold";

						} else if (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) {
							strUnits = vecClimoInfo[ct].threshorigunits;
							strStandardName = "threshold sum of " + strStandardName;
							strLongName = "threshold sum of " + strLongName;

						} else if (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) {
							strUnits = "times";
							strStandardName = "average max consecutive times " + strStandardName + " satisfies threshold";
							strLongName = "average max consecutive times " + strLongName + " satisfies threshold";
						}

						varOut->add_att("units", strUnits.c_str());
						varOut->add_att("standard_name", strStandardName.c_str());
						varOut->add_att("long_name", strLongName.c_str());
					}

					CopyNcVarAttributes(varIn, varOut, vecDoNotCopyNames);
				}

				// Add threshold attribute
				if ((vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
				    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) ||
				    (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) ||
					(vecClimoInfo[ct].type == ClimatologyType_MaxConsec)
				) {
					varOut->add_att("threshold", vecClimoInfo[ct].threshstring.c_str());
				}

				// Store output variable
				vecNcVarOut[v * vecClimoInfo.size() + ct] = varOut;
			}

			// Create output variable for slice count
			if (fOutputSliceCounts) {
				std::string strOutCountVar;
				if (vecVariableOutNames.size() != 0) {
					strOutCountVar = vecVariableOutNames[v];
				} else {
					strOutCountVar = vecVariableNames[v];
				}
				strOutCountVar += "_count";

				if (fMissingData) {
					vecNcVarCount[v] =
						ncoutfile.add_var(
							strOutCountVar.c_str(),
							ncInt,
							vecOutputNcDim[v].size(),
							const_cast<const NcDim**>(&(vecOutputNcDim[v][0])));

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
		}

		AnnounceEndBlock("Done");
	}

	/////////////////////////////////////////////
	// Actual climatology calculation

	// Loop through all Variables in ProcessingQueue
	varreg.ResetProcessingQueue();
	while(varreg.AdvanceProcessingQueue()) {

		// ProcessingQueue Variable and position
		size_t v = varreg.GetProcessingQueueVarPos();
		_ASSERT(v < vecVarIxIn.size());

		Variable & var = varreg.GetProcessingQueueVariable();

		AnnounceStartBlock("Calculating climatology of variable \"%s\"",
			vecVariableNames[v].c_str());

		// FillValue for this variable
		float dFillValue = std::numeric_limits<float>::max();
		if (strFillValueOverride != "") {
			if (!STLStringHelper::IsFloat(strFillValueOverride)) {
				_EXCEPTION1("FillValue override must be of type float: given \"%s\"",
					strFillValueOverride.c_str());
			}
			dFillValue = std::stof(strFillValueOverride);
		}

		// Prepare output auxiliary index and size
		VariableAuxIndex auxixOut = varreg.GetProcessingQueueAuxIx();
		auxixOut.insert(auxixOut.begin(), 0);
		VariableAuxSize auxsizeOut(auxixOut.size(),1);

		for (size_t d = 0; d < grid.DimCount(); d++) {
			auxixOut.push_back(0);
			auxsizeOut.push_back(grid.m_nGridDim[d]);
		}

		// Number of time slices within the current period.  If there is missing
		// data this is the number of time slices at each grid point.
		DataArray1D<int> nTimeSlices;
		DataArray2D<int> nTimeSlicesGrid;
		if (fMissingData) {
			nTimeSlicesGrid.Allocate(sOutputTimes, grid.GetSize());
		} else {
			nTimeSlices.Allocate(sOutputTimes);
		}

		// Number of time periods.  If there is missing data this is the number
		// of periods at each grid point.
		DataArray1D<int> nTimePeriods;
		DataArray2D<int> nTimePeriodsGrid;
		if (fMissingData) {
			nTimePeriodsGrid.Allocate(sOutputTimes, grid.GetSize());
		} else {
			nTimePeriods.Allocate(sOutputTimes);
		}

		// Accumulated data arrays
		std::vector< DataArray2D<double> > vecAccumulatedData(vecClimoInfo.size());
		for (int ct = 0; ct < vecClimoInfo.size(); ct++) {
			vecAccumulatedData[ct].Allocate(sOutputTimes, grid.GetSize());

			// Initialize the Min and Max arrays when first file is loaded
			DataArray2D<double> & dAccumulatedData = vecAccumulatedData[ct];
			if (vecClimoInfo[ct].type == ClimatologyType_Min) {
				dAccumulatedData.Set(std::numeric_limits<double>::max());
			}
			if (vecClimoInfo[ct].type == ClimatologyType_Max) {
				dAccumulatedData.Set(-std::numeric_limits<double>::max());
			}
		}

		// Scratch data
		int iScratchIndex = 0;
		std::vector< DataArray3D<double> > vecScratchData(vecClimoInfo.size());
		std::vector< DataArray3D<float> > vecScratchDataFloat(vecClimoInfo.size());

		for (int ct = 0; ct < vecClimoInfo.size(); ct++) {
			DataArray3D<double> & dScratchData = vecScratchData[ct];
			DataArray3D<float> & dScratchDataFloat = vecScratchDataFloat[ct];

			if (vecClimoInfo[ct].type == ClimatologyType_StdDev) {
				dScratchData.Allocate(1, sOutputTimes, grid.GetSize());

			} else if (vecClimoInfo[ct].type == ClimatologyType_AvgMin) {
				dScratchData.Allocate(1, 1, grid.GetSize());
				for (size_t i = 0; i < dScratchData.GetSubColumns(); i++) {
					dScratchData(0,0,i) = std::numeric_limits<double>::max();
				}

			} else if (vecClimoInfo[ct].type == ClimatologyType_AvgMax) {
				dScratchData.Allocate(1, 1, grid.GetSize());
				for (size_t i = 0; i < dScratchData.GetSubColumns(); i++) {
					dScratchData(0,0,i) = -std::numeric_limits<double>::max();
				}

			} else if (vecClimoInfo[ct].type == ClimatologyType_AvgCount) {
				dScratchData.Allocate(1, 1, grid.GetSize());

			} else if (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) {
				dScratchData.Allocate(1, 1, grid.GetSize());

			} else if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
				dScratchData.Allocate(1, 1, grid.GetSize());
				for (size_t i = 0; i < dScratchData.GetSubColumns(); i++) {
					dScratchData(0,0,i) = -1.0;
				}

			} else if (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) {
				dScratchData.Allocate(2, 1, grid.GetSize());

			} else if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
				dScratchData.Allocate(5, 1, grid.GetSize());
				dScratchDataFloat.Allocate(vecClimoInfo[ct].lag, 1, grid.GetSize());
			}
		}

		// Time when processing started
		Time timeProcessingStart(Time::CalendarUnknown);

		// Last time processed
		Time timeProcessingLast(Time::CalendarUnknown);

		// Number of times processed
		size_t nProcessingTimes = 0;

		// Time when filtering started
		Time timeFilterStart(Time::CalendarUnknown);

		// Last time filtered
		Time timeFilterLast(Time::CalendarUnknown);

		// Reason for filtering times
		std::string strFilterReason;

		// Number of times filtered
		size_t nFilterTimes = 0;

		// Current time (used for avgmin and avgmax calculations)
		Time timeCurrent(Time::CalendarUnknown);

		// Loop through each input file
		for (int f = 0; f < vecInputFileList.size(); f++) {

			// Initialize the NcFileVector
			NcFileVector vecncDataFiles;
			vecncDataFiles.ParseFromString(vecInputFileList[f].c_str());
			_ASSERT(vecncDataFiles.size() > 0);

			if (fVerbose) {
				Announce("%s", vecInputFileList[f].c_str());
			}
/*
			// Open input file
			NcFile ncinfile(vecInputFileList[f].c_str());
			if (!ncinfile.is_valid()) {
				_EXCEPTION1("Unable to open NetCDF file \"%s\"",
					vecInputFileList[f].c_str());
			}
*/
			// Get time information
			NcVar * varInTime = NcGetTimeVariable(*(vecncDataFiles[0]));
			if (varInTime == NULL) {
				_EXCEPTION1("Cannot load variable \"time\" from NetCDF file \"%s\"",
					vecncDataFiles.GetFilename(0).c_str());
			}
			if (varInTime->num_dims() != 1) {
				_EXCEPTION1("Variable \"time\" in NetCDF file \"%s\" has "
					"more than 1 dimension",
					vecncDataFiles.GetFilename(0).c_str());
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
			NcTimeDimension vecTimesIn;
			ReadCFTimeDataFromNcFile(
				vecncDataFiles[0], 
				vecncDataFiles.GetFilename(0),
				vecTimesIn,
				false);

			// Loop through all times
			for (int t = 0; t < vecTimesIn.size(); t++) {

				bool fFilterTime = false;

#ifndef TEMPEST_NOREGEX
				if (strTimeFilter != "") {
					std::string strTime = vecTimesIn[t].ToString();
					std::smatch match;
					if (!std::regex_search(strTime, match, reTimeSubset)) {
						fFilterTime = true;
						if (strFilterReason.length() == 0) {
							strFilterReason = "timefilter";
						} else if (strFilterReason.find("timefilter") == std::string::npos) {
							strFilterReason += "; timefilter";
						}
					}
				}
#endif
				// Filter by --time_start
				if ((!fFilterTime) && (strStartTime != "")) {
					double dDeltaSeconds = timeStartTime - vecTimesIn[t];
					if (dDeltaSeconds > 0.0) {
						fFilterTime = true;
						if (strFilterReason.length() == 0) {
							strFilterReason = "time_start";
						} else if (strFilterReason.find("time_start") == std::string::npos) {
							strFilterReason += "; time_start";
						}
					}
				}

				// Filter by --time_end
				if ((!fFilterTime) && (strEndTime != "")) {
					double dDeltaSeconds = vecTimesIn[t] - timeEndTime;
					if (dDeltaSeconds > 0.0) {
						fFilterTime = true;
						if (strFilterReason.length() == 0) {
							strFilterReason = "time_end";
						} else if (strFilterReason.find("time_end") == std::string::npos) {
							strFilterReason += "; time_end";
						}
					}
				}

				// Filter leap days
				if ((!fFilterTime) && (vecTimesIn[t].IsLeapDay())) {
					fFilterTime = true;
					if (strFilterReason.length() == 0) {
						strFilterReason = "leap day";
					} else if (strFilterReason.find("leap day") == std::string::npos) {
						strFilterReason += "; leap day";
					}
				}

				// Notify user of this time slice being filtered
				if (fFilterTime) {
					if (fVerbose) {
						Announce("Time %s (%s; skipping)",
							vecTimesIn[t].ToString().c_str(),
							strFilterReason.c_str());
					} else {
						if (timeFilterStart.GetCalendarType() == Time::CalendarUnknown) {
							timeFilterStart = vecTimesIn[t];
						}
						nFilterTimes++;
						timeFilterLast = vecTimesIn[t];
					}
					if (timeProcessingStart.GetCalendarType() != Time::CalendarUnknown) {
						Announce("Processed %lu times %s to %s",
							nProcessingTimes,
							timeProcessingStart.ToString().c_str(),
							timeProcessingLast.ToString().c_str());

						timeProcessingStart = Time(Time::CalendarUnknown);
						timeProcessingLast = Time(Time::CalendarUnknown);
						nProcessingTimes = 0;
					}
					continue;
				}

				// Notify user of this time slice being processed
				if (fVerbose) {
					Announce("Time %s", vecTimesIn[t].ToString().c_str());
				} else {
					if (timeFilterStart.GetCalendarType() != Time::CalendarUnknown) {
						Announce("Filtered %lu times %s to %s (%s)",
							nFilterTimes,
							timeFilterStart.ToString().c_str(),
							timeFilterLast.ToString().c_str(),
							strFilterReason.c_str());
					
						timeFilterStart = Time(Time::CalendarUnknown);
						timeFilterLast = Time(Time::CalendarUnknown);
						nFilterTimes = 0;
						strFilterReason = "";
					}
					if (timeProcessingStart.GetCalendarType() == Time::CalendarUnknown) {
						timeProcessingStart = vecTimesIn[t];
					}
					nProcessingTimes++;
					timeProcessingLast = vecTimesIn[t];
				}

				// Load the data for the search variable
				vecncDataFiles.SetTime(vecTimesIn[t]);
				var.LoadGridData(varreg, vecncDataFiles, grid);
				const DataArray1D<float> & dDataIn = var.GetData();
				_ASSERT(dDataIn.GetRows() == grid.GetSize());

				// Get FillValue
				if (strFillValueOverride != "") {
					if (dDataIn.HasFillValue()) {
						if (std::isnan(dFillValue)) {
							if (!std::isnan(dDataIn.GetFillValue())) {
								_EXCEPTION1("_FillValue changes value across files from NaN to %g", dFillValue);
							}
						} else {
							if (dFillValue == std::numeric_limits<float>::max()) {
								dFillValue = dDataIn.GetFillValue();
							} else if (std::isnan(dDataIn.GetFillValue())) {
								_EXCEPTION1("_FillValue changes value across files from %g to NaN", dFillValue);
							} else if (dFillValue != dDataIn.GetFillValue()) {
								_EXCEPTION2("_FillValue changes value across files from %g to %g",
									dFillValue, dDataIn.GetFillValue());
							}
						}
					} else if (dFillValue != std::numeric_limits<float>::max()) {
						_EXCEPTION1("_FillValue found in earlier file (%g) missing in subsequent file", dFillValue);
					}
				}

				// Check if we are in a new period and accumulate data
				if (fPeriodAccumulationNeeded) {

					if (timeCurrent.GetCalendarType() == Time::CalendarUnknown) {
						timeCurrent = vecTimesIn[t];
					}

					if (vecTimesIn[t] - timeCurrent < 0.0) {
						_EXCEPTIONT("For --type \"avgmin\", \"avgmax\" and \"autocor\" time across files must be monotone increasing");
					}

					// Check if we are in a new period
					bool fNewPeriod = false;
					if (eClimoPeriod == ClimatologyPeriod_Daily) {
						if ((vecTimesIn[t].GetDay() != timeCurrent.GetDay()) ||
						    (vecTimesIn[t].GetMonth() != timeCurrent.GetMonth()) ||
						    (vecTimesIn[t].GetYear() != timeCurrent.GetYear())
						) {
							fNewPeriod = true;
						}

					} else if (eClimoPeriod == ClimatologyPeriod_Monthly) {
						if ((vecTimesIn[t].GetMonth() != timeCurrent.GetMonth()) ||
						    (vecTimesIn[t].GetYear() != timeCurrent.GetYear())
						) {
							fNewPeriod = true;
						}

					} else if (eClimoPeriod == ClimatologyPeriod_Seasonal) {
						if (vecTimesIn[t].GetSeasonIndex() != timeCurrent.GetSeasonIndex()) {
							fNewPeriod = true;
						}

					} else if (eClimoPeriod == ClimatologyPeriod_Annual) {
						if (vecTimesIn[t].GetYear() != timeCurrent.GetYear()) {
							fNewPeriod = true;
						}
					}

					// If we're in a new period accumulate data
					int iCurrentTimeIndex = GetClimatologyTimeIndex(timeCurrent, eClimoPeriod);
					if (fNewPeriod) {
						Announce("Accumulating index %i; Next time %s",
							iCurrentTimeIndex, vecTimesIn[t].ToString().c_str());

						// Increment the number of time periods at each grid point
						if (fMissingData) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if (nTimeSlicesGrid(iCurrentTimeIndex,i) != 0) {
									nTimePeriodsGrid(iCurrentTimeIndex,i)++;
								}
							}
						}

						// Loop through all climatology types
						for (int ct = 0; ct < vecClimoInfo.size(); ct++) {

							DataArray3D<double> & dScratchData = vecScratchData[ct];
							DataArray2D<double> & dAccumulatedData = vecAccumulatedData[ct];

							// Accumulate avgmin, avgmax, avgcount and threshsum
							if ((vecClimoInfo[ct].type == ClimatologyType_AvgMin) ||
							    (vecClimoInfo[ct].type == ClimatologyType_AvgMax) ||
							    (vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
							    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum)
							) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (fabs(dScratchData(0,0,i)) == std::numeric_limits<double>::max()) {
										continue;
									}
									if (fMissingData) {
										if (dScratchData(0,0,i) == dFillValue) {
											continue;
										}
										if (dScratchData(0,0,i) != dScratchData(0,0,i)) {
											continue;
										}
									}
									dAccumulatedData(iCurrentTimeIndex,i) += dScratchData(0,0,i);
								}
							}

							// Accumulate autocorrelation
							if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
								double dModCount = 1.0 / static_cast<double>(nTimeSlices[0] - vecClimoInfo[ct].lag);

								for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
									double dEX2 = dScratchData(0,0,i) * dModCount;
									double dEY2 = dScratchData(1,0,i) * dModCount;

									double dEX = dScratchData(2,0,i) * dModCount;
									double dEY = dScratchData(3,0,i) * dModCount;

									double dEXY = dScratchData(4,0,i) * dModCount;

									double dVarX = dEX2 - dEX * dEX;
									double dVarY = dEY2 - dEY * dEY;

									if ((fabs(dVarX) <= 1.0e-12 * dEX2) || (fabs(dVarY) <= 1.0e-12 * dEY2)) {
										dAccumulatedData(iCurrentTimeIndex,i) += 1.0;
									} else {
										dAccumulatedData(iCurrentTimeIndex,i) +=
											(dEXY - dEX * dEY) / sqrt(dVarX) / sqrt(dVarY);
									}
								}
							}

							// Accumulate timeuntil
							if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
									}
									dAccumulatedData(iCurrentTimeIndex,i) += dScratchData(0,0,i);
								}
							}

							// Accumulate maxconsec (check both length of current period
							// and longest period to date)
							if (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
										dAccumulatedData(iCurrentTimeIndex,i) += dScratchData(1,0,i);
									} else {
										dAccumulatedData(iCurrentTimeIndex,i) += dScratchData(0,0,i);
									}
								}
							}

							// Reset scratch data
							if (vecClimoInfo[ct].type == ClimatologyType_AvgMin) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									dScratchData(0,0,i) = std::numeric_limits<double>::max();
								}
							}
							if (vecClimoInfo[ct].type == ClimatologyType_AvgMax) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									dScratchData(0,0,i) = -std::numeric_limits<double>::max();
								}
							}
							if ((vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
							    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) ||
								(vecClimoInfo[ct].type == ClimatologyType_MaxConsec)
							) {
								dScratchData.Zero();
							}
							if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
								dScratchData.Zero();
								iScratchIndex = 0;
								nTimeSlices[0] = 0;
							}
							if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
								for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
									dScratchData(0,0,i) = -1.0;
								}
							}
						}

						// Increment the number of periods (if no missing data)
						if (!fMissingData) {
							nTimePeriods(iCurrentTimeIndex)++;
						}

						// Update the beginning of the accumulation time
						timeCurrent = vecTimesIn[t];
					}
				}

				// Current time index
				int iTimeIndex = GetClimatologyTimeIndex(vecTimesIn[t], eClimoPeriod);

				// Time slice count
				if (fMissingData) {
					for (size_t i = 0; i < dDataIn.GetRows(); i++) {
						if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
							nTimeSlicesGrid(iTimeIndex,i)++;
						}
					}
				} else {
					nTimeSlices[iTimeIndex]++;
				}

				// Loop through all climatology types
				for (int ct = 0; ct < vecClimoInfo.size(); ct++) {

					DataArray3D<double> & dScratchData = vecScratchData[ct];
					DataArray3D<float> & dScratchDataFloat = vecScratchDataFloat[ct];
					DataArray2D<double> & dAccumulatedData = vecAccumulatedData[ct];

					// Get the threshold
					float dClimoThreshold = vecClimoInfo[ct].threshvalue;

					// Convert threshold units, if appropriate
					if (t == 0) {
						if (vecClimoInfo[ct].threshunits == "") {
							continue;
						}
						if (vecClimoInfo[ct].threshunits != dDataIn.GetUnits()) {
							float dThreshValueBak = vecClimoInfo[ct].threshvalue;

							bool fSuccess =
								ConvertUnits(
									vecClimoInfo[ct].threshvalue,
									vecClimoInfo[ct].threshunits,
									dDataIn.GetUnits());

							if (!fSuccess) {
								_EXCEPTION2("Unable to convert threshold from units %s to %s: No known conversion",
									vecClimoInfo[ct].threshunits.c_str(), dDataIn.GetUnits().c_str());
							}

							Announce("Converted threshold %g %s to %g %s",
								dThreshValueBak, vecClimoInfo[ct].threshunits.c_str(),
								vecClimoInfo[ct].threshvalue, dDataIn.GetUnits().c_str());

							vecClimoInfo[ct].threshunits = dDataIn.GetUnits();
						}
					}

					// Count number of time slices at each point and accumulate data
					if (fMissingData) {

						// Mean
						if (vecClimoInfo[ct].type == ClimatologyType_Mean) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									dAccumulatedData(iTimeIndex,i) += static_cast<double>(dDataIn[i]);
								}
							}

						// Expectation of the square
						} else if (vecClimoInfo[ct].type == ClimatologyType_MeanSq) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									dAccumulatedData(iTimeIndex,i) +=
										static_cast<double>(dDataIn[i])
										* static_cast<double>(dDataIn[i]);
								}
							}

						// Standard deviation
						// + dAccumulatedData stores the mean
						// + dScratchData(0) stores the expectation of the square
						} else if (vecClimoInfo[ct].type == ClimatologyType_StdDev) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									dAccumulatedData(iTimeIndex,i) +=
										static_cast<double>(dDataIn[i]);

									dScratchData(0,iTimeIndex,i) +=
										static_cast<double>(dDataIn[i])
										* static_cast<double>(dDataIn[i]);
								}
							}

						// Autocorrelation
						} else if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
							_EXCEPTION();

						// Minimum
						} else if (vecClimoInfo[ct].type == ClimatologyType_Min) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									if (dDataIn[i] < dAccumulatedData(iTimeIndex,i)) {
										dAccumulatedData(iTimeIndex,i) = dDataIn[i];
									}
								}
							}

						// Maximum
						} else if (vecClimoInfo[ct].type == ClimatologyType_Max) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									if (dDataIn[i] > dAccumulatedData(iTimeIndex,i)) {
										dAccumulatedData(iTimeIndex,i) = dDataIn[i];
									}
								}
							}

						// Average of the minimum
						// + dScratchData(0) stores the minimum over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgMin) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									if (dDataIn[i] < dScratchData(0,0,i)) {
										dScratchData(0,0,i) = dDataIn[i];
									}
								}
							}

						// Average of the maximum
						// + dScratchData(0) stores the maximum over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgMax) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
									if (dDataIn[i] > dScratchData(0,0,i)) {
										dScratchData(0,0,i) = dDataIn[i];
									}
								}
							}

						// Average of the count
						// + dScratchData(0) stores the count over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgCount) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] > dClimoThreshold) {
											dScratchData(0,0,i)++;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] >= dClimoThreshold) {
											dScratchData(0,0,i)++;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] < dClimoThreshold) {
											dScratchData(0,0,i)++;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] <= dClimoThreshold) {
											dScratchData(0,0,i)++;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_EqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] == dClimoThreshold) {
											dScratchData(0,0,i)++;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_NotEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] != dClimoThreshold) {
											dScratchData(0,0,i)++;
										}
									}
								}
							}

						// Sum over/under threshold
						// + dScratchData(0) stores the sum over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] > dClimoThreshold) {
											dScratchData(0,0,i) += (dDataIn[i] - dClimoThreshold);
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] < dClimoThreshold) {
											dScratchData(0,0,i) += (dClimoThreshold - dDataIn[i]);
										}
									}
								}
							}

						// Time until a condition is met
						// + dScratchData(0) stores -1.0 minus the count over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										
											if (dDataIn[i] > dClimoThreshold) {
												dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
											} else {
												dScratchData(0,0,i) -= 1.0;
											}
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
											if (dDataIn[i] >= dClimoThreshold) {
												dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
											} else {
												dScratchData(0,0,i) -= 1.0;
											}
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
											if (dDataIn[i] < dClimoThreshold) {
												dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
											} else {
												dScratchData(0,0,i) -= 1.0;
											}
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
											if (dDataIn[i] <= dClimoThreshold) {
												dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
											} else {
												dScratchData(0,0,i) -= 1.0;
											}
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_EqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
											if (dDataIn[i] == dClimoThreshold) {
												dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
											} else {
												dScratchData(0,0,i) -= 1.0;
											}
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_NotEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
											if (dDataIn[i] != dClimoThreshold) {
												dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
											} else {
												dScratchData(0,0,i) -= 1.0;
											}
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}
							}

						// Maximum consecutive days meeting a criteria
						// + dScratchData(0) stores the current maximum in this period
						// + dScratchData(1) stores the number of time slices in the current accumulation
						} else if (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									bool fSatisfiesThreshold = false;
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] > dClimoThreshold) {
											fSatisfiesThreshold = true;
										}
									}
									if (fSatisfiesThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									bool fSatisfiesThreshold = false;
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] >= dClimoThreshold) {
											fSatisfiesThreshold = true;
										}
									}
									if (fSatisfiesThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									bool fSatisfiesThreshold = false;
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] < dClimoThreshold) {
											fSatisfiesThreshold = true;
										}
									}
									if (fSatisfiesThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									bool fSatisfiesThreshold = false;
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] <= dClimoThreshold) {
											fSatisfiesThreshold = true;
										}
									}
									if (fSatisfiesThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_EqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									bool fSatisfiesThreshold = false;
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] == dClimoThreshold) {
											fSatisfiesThreshold = true;
										}
									}
									if (fSatisfiesThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_NotEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									bool fSatisfiesThreshold = false;
									if ((dDataIn[i] != dFillValue) && (!std::isnan(dDataIn[i]))) {
										if (dDataIn[i] != dClimoThreshold) {
											fSatisfiesThreshold = true;
										}
									}
									if (fSatisfiesThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}
							}
						}

					// Count number of time slices at each time and accumulate data
					} else {

						// Mean
						if (vecClimoInfo[ct].type == ClimatologyType_Mean) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								dAccumulatedData(iTimeIndex,i) += static_cast<double>(dDataIn[i]);
							}

						// Expectation of the square
						} else if (vecClimoInfo[ct].type == ClimatologyType_MeanSq) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								dAccumulatedData(iTimeIndex,i) +=
									static_cast<double>(dDataIn[i])
									* static_cast<double>(dDataIn[i]);
							}

						// Standard deviation
						} else if (vecClimoInfo[ct].type == ClimatologyType_StdDev) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								dAccumulatedData(iTimeIndex,i) +=
									static_cast<double>(dDataIn[i]);

								dScratchData(0,iTimeIndex,i) +=
									static_cast<double>(dDataIn[i])
									* static_cast<double>(dDataIn[i]);
							}

						// Autocorrelation
						// + dScratchData(0) stores the sum of the square of i=[0,n-k)
						// + dScratchData(1) stores the sum of the square of i=[k,n)
						// + dScratchData(2) stores the sum of i=[0,n-k)
						// + dScratchData(3) stores the sum of i=[k,n)
						// + dScratchData(4) stores sum_{i=k}^{n-1}(y_i y_{i-k})
						} else if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {

							if (nTimeSlices[0] >= vecClimoInfo[ct].lag) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									dScratchData(0,0,i) +=
										static_cast<double>(dScratchDataFloat(iScratchIndex,0,i))
										* static_cast<double>(dScratchDataFloat(iScratchIndex,0,i));

									dScratchData(1,0,i) +=
										static_cast<double>(dDataIn[i])
										* static_cast<double>(dDataIn[i]);

									dScratchData(2,0,i) +=
										static_cast<double>(dScratchDataFloat(iScratchIndex,0,i));

									dScratchData(3,0,i) +=
										static_cast<double>(dDataIn[i]);

									dScratchData(4,0,i) +=
										static_cast<double>(dDataIn[i])
										* static_cast<double>(dScratchDataFloat(iScratchIndex,0,i));
								}
							}
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								dScratchDataFloat(iScratchIndex,0,i) = dDataIn[i];
							}
							iScratchIndex++;
							if (iScratchIndex >= vecClimoInfo[ct].lag) {
								iScratchIndex = 0;
							}
							nTimeSlices[0]++;

						// Minimum
						} else if (vecClimoInfo[ct].type == ClimatologyType_Min) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if (dDataIn[i] < dAccumulatedData(iTimeIndex,i)) {
									dAccumulatedData(iTimeIndex,i) = dDataIn[i];
								}
							}

						// Maximum
						} else if (vecClimoInfo[ct].type == ClimatologyType_Max) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if (dDataIn[i] > dAccumulatedData(iTimeIndex,i)) {
									dAccumulatedData(iTimeIndex,i) = dDataIn[i];
								}
							}

						// Average of the minimum
						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgMin) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if (dDataIn[i] < dScratchData(0,0,i)) {
									dScratchData(0,0,i) = dDataIn[i];
								}
							}

						// Average of the maximum
						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgMax) {
							for (size_t i = 0; i < dDataIn.GetRows(); i++) {
								if (dDataIn[i] > dScratchData(0,0,i)) {
									dScratchData(0,0,i) = dDataIn[i];
								}
							}

						// Average of the count
						// + dScratchData(0) stores the count over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_AvgCount) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] > dClimoThreshold) {
										dScratchData(0,0,i)++;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] >= dClimoThreshold) {
										dScratchData(0,0,i)++;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] < dClimoThreshold) {
										dScratchData(0,0,i)++;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] <= dClimoThreshold) {
										dScratchData(0,0,i)++;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_EqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] == dClimoThreshold) {
										dScratchData(0,0,i)++;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_NotEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] != dClimoThreshold) {
										dScratchData(0,0,i)++;
									}
								}
							}

						// Sum over/under threshold
						// + dScratchData(0) stores the sum over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] > dClimoThreshold) {
										dScratchData(0,0,i) += (dDataIn[i] - dClimoThreshold);
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] < dClimoThreshold) {
										dScratchData(0,0,i) += (dClimoThreshold - dDataIn[i]);
									}
								}
							}

						// Time until a condition is met
						// + dScratchData(0) stores -1.0 minus the count over the current period
						} else if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if (dDataIn[i] > dClimoThreshold) {
											dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if (dDataIn[i] >= dClimoThreshold) {
											dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if (dDataIn[i] < dClimoThreshold) {
											dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if (dDataIn[i] <= dClimoThreshold) {
											dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_EqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if (dDataIn[i] == dClimoThreshold) {
											dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_NotEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dScratchData(0,0,i) < 0.0) {
										if (dDataIn[i] != dClimoThreshold) {
											dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
										} else {
											dScratchData(0,0,i) -= 1.0;
										}
									}
								}
							}

						// Maximum consecutive days meeting a criteria
						// + dScratchData(0) stores the current maximum in this period
						// + dScratchData(1) stores the number of time slices in the current accumulation
						} else if (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) {
							if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] > dClimoThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_GreaterThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] >= dClimoThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThan) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] < dClimoThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_LessThanOrEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] <= dClimoThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_EqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] == dClimoThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
										dScratchData(0,0,i) = dScratchData(1,0,i);
										dScratchData(1,0,i) = 0.0;
									}
								}

							} else if (vecClimoInfo[ct].threshtype == ThresholdType_NotEqualTo) {
								for (size_t i = 0; i < dDataIn.GetRows(); i++) {
									if (dDataIn[i] != dClimoThreshold) {
										dScratchData(1,0,i) += 1.0;
									} else {
										if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
											dScratchData(0,0,i) = dScratchData(1,0,i);
										}
										dScratchData(1,0,i) = 0.0;
									}
								}
							}
						}
					}
				}
			}
		}

		// Flag indicating number of time periods has been updated
		bool fUpdatedNumberOfTimePeriods = false;

		// Loop through all climatology types and perform post-processing on climatology
		for (int ct = 0; ct < vecClimoInfo.size(); ct++) {

			DataArray3D<double> & dScratchData = vecScratchData[ct];
			DataArray2D<double> & dAccumulatedData = vecAccumulatedData[ct];

			// Standard deviation: Summarize
			if (vecClimoInfo[ct].type == ClimatologyType_StdDev) {
				for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {
					if (fMissingData) {
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							double dTimeSliceCount = static_cast<double>(nTimeSlicesGrid(t,i));
							double dBesselCorrection = 1.0 / (dTimeSliceCount - 1.0);

							if (nTimeSlicesGrid(t,i) != 0) {
								double dSumX = dAccumulatedData(t,i);
								double dSumX2 = dScratchData(0,t,i);
								dAccumulatedData(t,i) =
									sqrt(dBesselCorrection * (dSumX2 - dSumX * dSumX / dTimeSliceCount));
							} else {
								dAccumulatedData(t,i) = dFillValue;
							}
						}

					} else {
						if (nTimeSlices[t] == 0) {
							continue;
						}
						double dTimeSliceCount = static_cast<double>(nTimeSlices(t));
						double dBesselCorrection = 1.0 / (dTimeSliceCount - 1.0);

						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							double dSumX = dAccumulatedData(t,i);
							double dSumX2 = dScratchData(0,t,i);
							dAccumulatedData(t,i) =
								sqrt(dBesselCorrection * (dSumX2 - dSumX * dSumX / dTimeSliceCount));
						}
					}
				}
			}

			// AvgMin, AvgMax, AvgCount, ThreshSum, TimeUntil, MaxConsec:
			// Perform one final accumulation of scratch data into accumulated data
			if ((vecClimoInfo[ct].type == ClimatologyType_AvgMin) ||
			    (vecClimoInfo[ct].type == ClimatologyType_AvgMax) ||
			    (vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
			    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) ||
				(vecClimoInfo[ct].type == ClimatologyType_TimeUntil) ||
				(vecClimoInfo[ct].type == ClimatologyType_MaxConsec)
			) {
				if (timeCurrent.GetCalendarType() == Time::CalendarUnknown) {
					_EXCEPTIONT("Data contains no valid time steps given command-line constraints");
				}

				if (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) {
					for (size_t i = 0; i < dScratchData.GetSubColumns(); i++) {
						if (dScratchData(0,0,i) < 0.0) {
							dScratchData(0,0,i) = -1.0 - dScratchData(0,0,i);
						}
					}
				}
				if (vecClimoInfo[ct].type == ClimatologyType_MaxConsec) {
					for (size_t i = 0; i < dScratchData.GetSubColumns(); i++) {
						if (dScratchData(1,0,i) > dScratchData(0,0,i)) {
							dScratchData(0,0,i) = dScratchData(1,0,i);
						}
					}
				}
				int iCurrentTimeIndex = GetClimatologyTimeIndex(timeCurrent, eClimoPeriod);

				if (!fUpdatedNumberOfTimePeriods) {
					AnnounceStartBlock("Accumulating index %i; Final", iCurrentTimeIndex);
				}
				if (fMissingData) {
					for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
						if (nTimeSlicesGrid(iCurrentTimeIndex,i) == 0) {
							continue;
						}
						if (!fUpdatedNumberOfTimePeriods) {
							nTimePeriodsGrid(iCurrentTimeIndex,i)++;
						}
						dAccumulatedData(iCurrentTimeIndex,i) += dScratchData(0,0,i);
					}

				} else {
					if (nTimeSlices(iCurrentTimeIndex) != 0) {
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							dAccumulatedData(iCurrentTimeIndex,i) += dScratchData(0,0,i);
						}
						if (!fUpdatedNumberOfTimePeriods) {
							nTimePeriods(iCurrentTimeIndex)++;
						}
					}
				}
				if (!fUpdatedNumberOfTimePeriods) {
					fUpdatedNumberOfTimePeriods = true;
					AnnounceEndBlock(NULL);
				}
			}

			// AutoCor: Perform one final accumulation of scratch data into accumulated data
			if (vecClimoInfo[ct].type == ClimatologyType_AutoCor) {
				int iCurrentTimeIndex = GetClimatologyTimeIndex(timeCurrent, eClimoPeriod);

				if (!fUpdatedNumberOfTimePeriods) {
					AnnounceStartBlock("Accumulating index %i; Final", iCurrentTimeIndex);
				}
				if (nTimeSlices(iCurrentTimeIndex) != 0) {
					double dModCount = 1.0 / static_cast<double>(nTimeSlices[0] - vecClimoInfo[ct].lag);

					for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
						double dEX2 = dScratchData(0,0,i) * dModCount;
						double dEY2 = dScratchData(1,0,i) * dModCount;

						double dEX = dScratchData(2,0,i) * dModCount;
						double dEY = dScratchData(3,0,i) * dModCount;

						double dEXY = dScratchData(4,0,i) * dModCount;

						double dVarX = dEX2 - dEX * dEX;
						double dVarY = dEY2 - dEY * dEY;

						if ((fabs(dVarX) <= 1.0e-12 * dEX2) || (fabs(dVarY) <= 1.0e-12 * dEY2)) {
							dAccumulatedData(iCurrentTimeIndex,i) += 1.0;
						} else {
							dAccumulatedData(iCurrentTimeIndex,i) +=
								(dEXY - dEX * dEY) / sqrt(dVarX) / sqrt(dVarY);
						}
					}

					nTimePeriods(iCurrentTimeIndex)++;
				}
				if (!fUpdatedNumberOfTimePeriods) {
					fUpdatedNumberOfTimePeriods = true;
					AnnounceEndBlock(NULL);
				}
			}

			// Calculate mean over time slices
			if ((vecClimoInfo[ct].type == ClimatologyType_Mean) ||
			    (vecClimoInfo[ct].type == ClimatologyType_MeanSq)
			) {
				for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {

					// ..with missing data
					if (fMissingData) {
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							if (nTimeSlicesGrid(t,i) != 0) {
								dAccumulatedData(t,i) /=
									static_cast<double>(nTimeSlicesGrid(t,i));
							} else {
								dAccumulatedData(t,i) = dFillValue;
							}
						}

					// ..without missing data
					} else {
						if (nTimeSlices[t] == 0) {
							continue;
						}
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							dAccumulatedData(t,i) /= static_cast<double>(nTimeSlices[t]);
						}
					}
				}
			}

			// Calculate mean over time periods
			if ((vecClimoInfo[ct].type == ClimatologyType_AutoCor) ||
			    (vecClimoInfo[ct].type == ClimatologyType_AvgMin) ||
			    (vecClimoInfo[ct].type == ClimatologyType_AvgMax) ||
			    (vecClimoInfo[ct].type == ClimatologyType_AvgCount) ||
			    (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) || 
			    (vecClimoInfo[ct].type == ClimatologyType_TimeUntil) ||
			    (vecClimoInfo[ct].type == ClimatologyType_MaxConsec)
			) {
				for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {

					// ..with missing data
					if (fMissingData) {
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							if (nTimePeriodsGrid(t,i) != 0) {
								dAccumulatedData(t,i) /=
									static_cast<double>(nTimePeriodsGrid(t,i));
							} else {
								dAccumulatedData(t,i) = dFillValue;
							}
						}

					// ..without missing data
					} else {
						if (nTimePeriods[t] == 0) {
							continue;
						}
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							dAccumulatedData(t,i) /= static_cast<double>(nTimePeriods[t]);
						}
					}
				}
			}

			// Convert max and min to FillValue
			if ((vecClimoInfo[ct].type == ClimatologyType_Min) ||
			    (vecClimoInfo[ct].type == ClimatologyType_Max)
			) {
				if (fMissingData) {
					for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {
					for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
						if (fabs(dAccumulatedData(t,i)) == std::numeric_limits<double>::max()) {
							dAccumulatedData(t,i) = dFillValue;
						}
					}
					}
				}
			}

			// Notify user of filtered or processed times
			if (timeFilterStart.GetCalendarType() != Time::CalendarUnknown) {
				Announce("Filtered %lu times %s to %s (%s)",
					nFilterTimes,
					timeFilterStart.ToString().c_str(),
					timeFilterLast.ToString().c_str(),
					strFilterReason.c_str());

				timeFilterStart = Time(Time::CalendarUnknown);
				timeFilterLast = Time(Time::CalendarUnknown);
				strFilterReason = "";
				nFilterTimes = 0;
			}
			if (timeProcessingStart.GetCalendarType() != Time::CalendarUnknown) {
				Announce("Processed %lu times %s to %s",
					nProcessingTimes,
					timeProcessingStart.ToString().c_str(),
					timeProcessingLast.ToString().c_str());

				timeProcessingStart = Time(Time::CalendarUnknown);
				timeProcessingLast = Time(Time::CalendarUnknown);
				nProcessingTimes = 0;
			}

			// Convert threshold sum back to original units
			if (vecClimoInfo[ct].type == ClimatologyType_ThreshSum) {
				if (vecClimoInfo[ct].threshorigunits != vecClimoInfo[ct].threshunits) {
					Announce("Converting data units from \"%s\" to \"%s\"",
						vecClimoInfo[ct].threshunits.c_str(), 
						vecClimoInfo[ct].threshorigunits.c_str());

					if (fMissingData) {
						for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							if (!std::isnan(dAccumulatedData(t,i)) &&
							    (dAccumulatedData(t,i) != dFillValue)
							) {
								ConvertUnits(
									dAccumulatedData(t,i),
									vecClimoInfo[ct].threshunits,
									vecClimoInfo[ct].threshorigunits,
									true);
							}
						}
						}
					} else {
						for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {
						for (size_t i = 0; i < dAccumulatedData.GetColumns(); i++) {
							ConvertUnits(
								dAccumulatedData(t,i),
								vecClimoInfo[ct].threshunits,
								vecClimoInfo[ct].threshorigunits,
								true);
						}
						}
					}
				}
			}

			// Write data
			AnnounceStartBlock("Writing \"%s\" to file", vecNcVarOut[v * vecClimoInfo.size() + ct]->name());
			for (size_t t = 0; t < dAccumulatedData.GetRows(); t++) {

				// Convert to float and write to disk
				DataArray1D<float> dOutputData(dAccumulatedData.GetColumns());
				for (size_t i = 0; i < dOutputData.GetRows(); i++) {
					dOutputData[i] = static_cast<float>(dAccumulatedData(t,i));
				}

				auxixOut[0] = t;
				vecNcVarOut[v * vecClimoInfo.size() + ct]->set_cur(&(auxixOut[0]));
				vecNcVarOut[v * vecClimoInfo.size() + ct]->put(&(dOutputData[0]), &(auxsizeOut[0]));

				// Write slice counts
				if ((fOutputSliceCounts) && (ct == 0)) {
					if (fMissingData) {
						vecNcVarCount[v]->set_cur(&(auxixOut[0]));
						vecNcVarCount[v]->put(&(nTimeSlicesGrid(t,0)), &(auxsizeOut[0]));
					} else {
						vecNcVarCount[v]->set_cur(t);
						vecNcVarCount[v]->put(&(nTimeSlices(t)), 1);
					}
				}
			}
			AnnounceEndBlock("Done");
		}
		AnnounceEndBlock("Done");
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
	size_t sMemoryMax,
	const std::string & strFillValueOverride
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

				if (strFillValueOverride != "") {
					if (strFillValueOverride == "nan") {
						if (var->type() == ncFloat) {
							dFillValueFloat = NAN;
						} else if (var->type() == ncDouble) {
							dFillValueDouble = NAN;
						} else {
							_EXCEPTION();
						}
					} else if (STLStringHelper::IsFloat(strFillValueOverride)) {
						if (var->type() == ncFloat) {
							dFillValueFloat = std::stof(strFillValueOverride);
						} else if (var->type() == ncDouble) {
							dFillValueDouble = std::stod(strFillValueOverride);
						} else {
							_EXCEPTION();
						}
					} else {
						_EXCEPTIONT("Invalid --fillvalue argument");
					}

				} else {
					NcAtt * attFill = var->get_att("_FillValue");
					if (attFill == NULL) {
						attFill = var->get_att("missing_value");
					}
					if (attFill != NULL) {
						if (attFill->type() != nctype) {
							_EXCEPTION2("Variable \"%s\" attribute \"%s\" has incompatible type",
								vecVariableNames[v].c_str(), attFill->name());
						}
						if (nctype == ncFloat) {
							dFillValueFloat = attFill->as_float(0);
						}
						if (nctype == ncDouble) {
							dFillValueDouble = attFill->as_double(0);
						}
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
/*
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
*/
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
/*
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
*/
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
/*
						llCountReadSize =
							NcGetPutSpecifiedDataSize<int>(
								varCount,
								lPos0,
								sAuxSize,
								c,
								nDataCount,
								NetCDF_Get);
*/
						for (size_t s = 0; s < llCountReadSize; s++) {
							nAccumDataCount[s] += nDataCount[s];
						}
					}
				}

				// Read data
				if (nctype == ncFloat) {
/*
					llCurrentReadSize =
						NcGetPutSpecifiedDataSize<float>(
							var,
							lPos0,
							sAuxSize,
							c,
							dDataFloat,
							NetCDF_Get);
*/
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
/*
					llCurrentReadSize =
						NcGetPutSpecifiedDataSize<double>(
							var,
							lPos0,
							sAuxSize,
							c,
							dDataDouble,
							NetCDF_Get);
*/
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
/*
				llCurrentWriteSize =
					NcGetPutSpecifiedDataSize<float>(
						varOut,
						lPos0,
						sAuxSize,
						c,
						dAccumDataFloat,
						NetCDF_Put);
*/
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
/*
				llCurrentWriteSize =
					NcGetPutSpecifiedDataSize<double>(
						varOut,
						lPos0,
						sAuxSize,
						c,
						dAccumDataDouble,
						NetCDF_Put);
*/
			}
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

	// Connectivity file
	std::string strConnectivity;

	// Output data file
	std::string strOutputFile;

	// Variable to use
	std::string strVariables;

	// Output variable names
	std::string strOutputVariables;

	// Period of the climatology
	std::string strClimoPeriod;

	// Type of climatology
	std::string strClimoType;

	// Thresholds
	std::string strThresholds;

	// Include leap days
	bool fIncludeLeapDays;

	// Start time
	std::string strStartTime;

	// End time
	std::string strEndTime;

	// Dataset has missing data
	bool fMissingData;

	// Output slice counts
	bool fOutputSliceCounts;

	// Path to write temporary cliamtology files
	std::string strTempFilePath;

	// Override the fillvalue
	std::string strFillValueOverride;

	// Time filter
	std::string strTimeFilter;

	// Latitude dimension name
	std::string strLatitudeName;

	// Longitude dimension name
	std::string strLongitudeName;

	// Diagonal connectivity
	bool fDiagonalConnectivity;

	// Regional grid
	bool fRegional;

	// Do not delete temp files after completing climatology
	bool fKeepTempFiles;

	// Verbose
	bool fVerbose;

	// Command line
	std::string strCommandLine = GetCommandLineAsString(argc, argv);

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutputFile, "out_data", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputVariables, "varout", "");
		CommandLineStringD(strClimoPeriod, "period", "daily", "[daily|monthly|seasonal|annual|all]");
		CommandLineStringD(strClimoType, "type", "mean", "[mean|meansq|stddev|autocor|min|max|avgmin|avgmax|avgcount|threshsum|timeuntil|maxconsec]");
		CommandLineBool(fIncludeLeapDays, "include_leap_days");
		CommandLineString(strStartTime, "time_start", "");
		CommandLineString(strEndTime, "time_end", "");
		CommandLineBool(fMissingData, "missingdata");
		CommandLineBool(fOutputSliceCounts, "write_counts");
		CommandLineString(strTempFilePath, "temp_file_path", "/tmp");
		CommandLineString(strFillValueOverride, "fillvalue", "");
		CommandLineString(strTimeFilter, "timefilter", "");
		CommandLineString(strLongitudeName, "lonname", "[auto]");
		CommandLineString(strLatitudeName, "latname", "[auto]");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");
		CommandLineBool(fKeepTempFiles, "keep_temp_files");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

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
	if (strVariables.length() == 0) {
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
	} else if (strClimoPeriod == "all") {
		eClimoPeriod = ClimatologyPeriod_All;
	} else {
		_EXCEPTIONT("--period invalid; expected \"daily\", \"monthly\", \"seasonal\", \"annual\" or \"all\"");
	}

	// Parse climatology types
	std::vector<ClimatologyInfo> vecClimoInfo;
	ParseClimatologyTypes(strClimoType, vecClimoInfo);

	// Parse input file list
	FilenameList vecInputFileList;
	if (strInputFile.length() != 0) {
		std::vector<std::string> vecInputFilesParsed;
		STLStringHelper::ParseVariableList(strInputFile, vecInputFilesParsed, ";");
		for (auto & strInputFileParsed : vecInputFilesParsed) {
			vecInputFileList.push_back(strInputFileParsed);
		}
	} else {
		vecInputFileList.FromFile(strInputFileList, false);
	}
	_ASSERT(vecInputFileList.size() != 0);
/*
	// Parse variable list
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVariable, vecVariableStrings);

	std::vector<std::string> vecVariableNames;
	std::vector< std::vector<std::string> > vecVariableSpecifiedDims;

	STLStringHelper::SplitVariableStrings(
		vecVariableStrings,
		vecVariableNames,
		vecVariableSpecifiedDims
	);

	// Parse output variable list
	std::vector<std::string> vecVariableOutNames;
	if (strVariableOut != "") {
		STLStringHelper::ParseVariableList(strVariableOut, vecVariableOutNames);
		if (vecVariableOutNames.size() != vecVariableStrings.size()) {
			_EXCEPTIONT("If specified, --varout and --var must have the same length");
		}
	}
*/
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

	if (nMPISize != 1) {
		_EXCEPTIONT("Sorry! MPI Support in Climatology is currently disabled. It may be added in later.");
/*
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
		snprintf(szBuffer, 16, "%06i", nMPIRank);

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
*/
	} else {
		vecMyInputFileList = vecInputFileList;
		strMyOutputFile = strOutputFile;
	}

	// Calculate the climatology
	Climatology(
		vecMyInputFileList,
		strConnectivity,
		strMyOutputFile,
		strVariables,
		strOutputVariables,
		fIncludeLeapDays,
		strStartTime,
		strEndTime,
		eClimoPeriod,
		vecClimoInfo,
		fMissingData,
		strFillValueOverride,
		strTimeFilter,
		strLatitudeName,
		strLongitudeName,
		fDiagonalConnectivity,
		fRegional,
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
					snprintf(szBuffer, 16, "%06i", f);

					std::string strTempOutputFile =
						strTempFilePath + "/tempClimatology" + szBuffer + ".nc";

					vecAllOutputFileList.push_back(strTempOutputFile);
				}

				std::vector<std::string> vecTempVariableNames;

				std::string strVariablePrefix =
					ClimoVariablePrefix(eClimoPeriod, vecClimoInfo[0]);

				for (int v = 0; v < vecVariableNames.size(); v++) {
					vecTempVariableNames.push_back(strVariablePrefix + vecVariableNames[v]);
				}

				AverageOverNcFiles(
					vecAllOutputFileList,
					strOutputFile,
					vecTempVariableNames,
					sMemoryMax,
					strFillValueOverride);
			}
		}
	}
#endif

	// Add provenance information
	if (nMPIRank == 0) {
		const std::time_t timetNow = std::time(nullptr);
		std::string strProvenance = std::asctime(std::localtime(&timetNow));
		strProvenance += ": " + strCommandLine;

		NcFile ncoutfile(strOutputFile.c_str(), NcFile::Write);
		if (!ncoutfile.is_valid()) {
			_EXCEPTION1("Unable to open output datafile \"%s\"",
				strOutputFile.c_str());
		}
		ncoutfile.add_att("history", strProvenance.c_str());
	}

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

