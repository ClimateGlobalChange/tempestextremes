///////////////////////////////////////////////////////////////////////////////
///
///	\file    AccumulateERA5Forecast.cpp
///	\author  Paul Ullrich
///	\version June 20, 2021
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

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "NetCDFUtilities.h"
#include "Variable.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "Variable.h"
#include "FilenameList.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>
#include <ctime>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

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

#if defined(TEMPEST_MPIOMP)
	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
	if (nMPISize > 1) {
		_EXCEPTIONT("At present QuantileCalculator only supports serial execution.");
	}
#endif

	// Input directory
	std::string strInputDir;

	// Input data
	std::string strInputData;

	// Input data list
	std::string strInputDataList;

	// Input year
	int nInputYear;

	// Input month
	int nInputMonth;

	// Output data
	std::string strOutputData;

	// Accumulation frequency
	std::string strAccumFrequency;

	// File string
	std::string strFileString;

	// Variable to use for quantile calculation
	std::string strVariableName;

	// Output variable
	std::string strVariableOutName;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Command line
	std::string strCommandLine = GetCommandLineAsString(argc, argv);

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputDir, "in_dir", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineInt(nInputYear, "year", -1);
		CommandLineInt(nInputMonth, "month", -1);
		CommandLineString(strOutputData, "out_data", "");
		CommandLineStringD(strAccumFrequency, "accumfreq", "6h", "[1h|3h|6h|daily]");
		CommandLineString(strFileString, "file", "*235_055*");
		CommandLineString(strVariableName, "var", "MTPR");
		CommandLineString(strVariableOutName, "varout", "tp");

		CommandLineString(strLongitudeName, "lonname", "longitude");
		CommandLineString(strLatitudeName, "latname", "latitude");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	int nInputCount =
		  ((strInputDir.length() == 0)?(0):(1))
		+ ((strInputData.length() == 0)?(0):(1))
		+ ((strInputDataList.length() == 0)?(0):(1));

	if (nInputCount == 0) {
		_EXCEPTIONT("No input directory (--in_dir), file (--in_data) or (--in_data_list) specified");
	}
	if (nInputCount > 1) {
		_EXCEPTIONT("Only one of (--in_dir), (--in_data) or (--in_data_list) may be specified");
	}
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}

	if (strVariableName == "") {
		_EXCEPTIONT("No variable (--var) specified");
	}
	if (strVariableOutName == "") {
		strVariableOutName = strVariableName;
	}

	// Accumulation frequency
	long lAccumFrequency = 0;
	if (strAccumFrequency == "1h") {
		lAccumFrequency = 1;
	} else if (strAccumFrequency == "3h") {
		lAccumFrequency = 3;
	} else if (strAccumFrequency == "6h") {
		lAccumFrequency = 6;
	} else if ((strAccumFrequency == "daily") || (strAccumFrequency == "24h")) {
		lAccumFrequency = 24;
	} else {
		_EXCEPTIONT("--accumfreq must be \"1h\", \"3h\", \"6h\" or \"daily\"");
	}

	// Input files
	FilenameList vecInputFiles;

	if (strInputData != "") {
		vecInputFiles.push_back(strInputData);
	}
	if (strInputDataList != "") {
		vecInputFiles.FromFile(strInputDataList);
	}
	if (strInputDir != "") {
		int nInputMonthPrev = nInputMonth - 1;
		int nInputYearPrev = nInputYear;
		if (nInputMonthPrev == 0) {
			nInputYearPrev--;
			nInputMonthPrev = 12;
		}

		char szYearMonth[32];
		snprintf(szYearMonth, 32, "%04i%02i", nInputYearPrev, nInputMonthPrev);

		std::string strTempFile = std::string("/tmp/accum") + std::string(szYearMonth);

		std::string strSystemString = std::string("rm -rf ") + strTempFile;

		Announce(strSystemString.c_str());
		system(strSystemString.c_str());

		std::string strYearMonthDir = strInputDir + "/" + std::string(szYearMonth);

		strSystemString = std::string("ls ") + strYearMonthDir + std::string("/") + strFileString + std::string(" >> ") + strTempFile;

		Announce(strSystemString.c_str());
		system(strSystemString.c_str());

		snprintf(szYearMonth, 32, "%04i%02i", nInputYear, nInputMonth);
		strYearMonthDir = strInputDir + "/" + std::string(szYearMonth);

		strSystemString = std::string("ls ") + strYearMonthDir + std::string("/") + strFileString + std::string(" >> ") + strTempFile;

		Announce(strSystemString.c_str());
		system(strSystemString.c_str());

		vecInputFiles.FromFile(strTempFile);
	}
	if (vecInputFiles.size() == 0) {
		_EXCEPTIONT("No input data files found");
	}

	// Output file
	NcFile ncfileout(strOutputData.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
	if (!ncfileout.is_valid()) {
		_EXCEPTION1("Cannot open file \"%s\" for writing", strOutputData.c_str());
	}
	NcDim * dimTimeOut = ncfileout.add_dim("time", 0);
	if (dimTimeOut == NULL) {
		_EXCEPTION1("Cannot create dimension \"time\" in \"%s\"", strOutputData.c_str());
	}

	// Build time array across files
	AnnounceStartBlock("Building time array across files");
	std::vector<int> vecTimeInt;

	std::string strTimeCalendar;
	Time::CalendarType eCalType;
	std::string strTimeUnits;

	for (int f = 0; f < vecInputFiles.size(); f++) {

		NcFile ncfilein(vecInputFiles[f].c_str(), NcFile::ReadOnly);
		if (!ncfilein.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"", vecInputFiles[f].c_str());
		}

		NcDim * dimTimeInitial = ncfilein.get_dim("forecast_initial_time");
		if (dimTimeInitial == NULL) {
			_EXCEPTIONT("Unable to open dimension \"forecast_initial_time\"");
		}

		NcDim * dimTimeHour = ncfilein.get_dim("forecast_hour");
		if (dimTimeHour == NULL) {
			_EXCEPTIONT("Unable to open dimension \"forecast_hour\"");
		}

		NcVar * varTimeInitial = ncfilein.get_var("forecast_initial_time");
		if (varTimeInitial == NULL) {
			_EXCEPTION1("Cannot find variable \"forecast_initial_time\" in \"%s\"", vecInputFiles[f].c_str());
		}

		NcVar * varTimeHour = ncfilein.get_var("forecast_hour");
		if (varTimeHour == NULL) {
			_EXCEPTION1("Cannot find variable \"forecast_hour\" in \"%s\"", vecInputFiles[f].c_str());
		}

		NcAtt * attTimeCalendar = varTimeInitial->get_att("calendar");
		if (attTimeCalendar == NULL) {
			_EXCEPTION1("Variable \"forecast_initial_time\" missing \"calendar\" attribute in \"%s\"", vecInputFiles[f].c_str());
		}
		strTimeCalendar = attTimeCalendar->as_string(0);
		eCalType = Time::CalendarTypeFromString(strTimeCalendar);

		NcAtt * attTimeUnits = varTimeInitial->get_att("units");
		if (attTimeUnits == NULL) {
			_EXCEPTION1("Variable \"forecast_initial_time\" missing \"units\" attribute in \"%s\"", vecInputFiles[f].c_str());
		}
		strTimeUnits = attTimeUnits->as_string(0);

		if (lAccumFrequency <= dimTimeHour->size()) {
			for (long lFI = 0; lFI < dimTimeInitial->size(); lFI++) {
				Time timeInitial(eCalType);

				int nTimeInitial;
				varTimeInitial->set_cur(lFI);
				varTimeInitial->get(&nTimeInitial, 1);

				timeInitial.FromCFCompliantUnitsOffsetInt(strTimeUnits, nTimeInitial);

				for (long lFH = 0; lFH < dimTimeHour->size(); lFH += lAccumFrequency) {
					Time timeForecast = timeInitial;
					timeForecast.AddHours(lFH + lAccumFrequency);

					if ((nInputYear == (-1)) || (timeForecast.GetYear() == nInputYear)) {
						if ((nInputMonth == (-1)) || (timeForecast.GetMonth() == nInputMonth)) {
							vecTimeInt.push_back(nTimeInitial + lFH + lAccumFrequency);
						}
					}
				}
			}

		} else {
			_EXCEPTIONT("Not implemented");
		}
	}

	AnnounceEndBlock("Done");

	// Write relevant variables to output file
	AnnounceStartBlock("Initializing output file");
	NcVar * varOut;

	long lLatitudeCount;
	long lLongitudeCount;
	{
		NcFile ncfilein(vecInputFiles[0].c_str(), NcFile::ReadOnly);
		if (!ncfilein.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"", vecInputFiles[0].c_str());
		}

		NcVar * varIn = ncfilein.get_var(strVariableName.c_str());
		if (varIn == NULL) {
			_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
				vecInputFiles[0].c_str(), strVariableName.c_str());
		}
	
		NcVar * varTimeOut = ncfileout.add_var("time", ncInt, dimTimeOut);
		if (varTimeOut == NULL) {
			_EXCEPTION1("Cannot create variable \"time\" in \"%s\"", strOutputData.c_str());
		}

		varTimeOut->put(&(vecTimeInt[0]), vecTimeInt.size());

		varTimeOut->add_att("long_name", "time");
		varTimeOut->add_att("short_name", "time");
		varTimeOut->add_att("units", strTimeUnits.c_str());
		varTimeOut->add_att("calendar", strTimeCalendar.c_str());
		varTimeOut->add_att("bounds", "time_bnds");

		NcDim * dimTimeBnds = ncfileout.add_dim("bounds", 2);
		if (dimTimeBnds == NULL) {
			_EXCEPTIONT("Unable to create dimension \"bounds\" in output file");
		}

		NcVar * varTimeBnds = ncfileout.add_var("time_bnds", ncInt, dimTimeOut, dimTimeBnds);
		if (dimTimeBnds == NULL) {
			_EXCEPTIONT("Unable to create dimension \"time_bnds\" in output file");
		}
	
		{
			DataArray2D<int> nTimeBnds(vecTimeInt.size(), 2);
			for (int t = 0; t < vecTimeInt.size(); t++) {
				nTimeBnds(t,0) = vecTimeInt[t] - static_cast<int>(lAccumFrequency);
				nTimeBnds(t,1) = vecTimeInt[t];
			}
			varTimeBnds->set_cur(0,0);
			varTimeBnds->put(&(nTimeBnds(0,0)), vecTimeInt.size(), 2);

			varTimeBnds->add_att("long_name", "time bounds");
			varTimeBnds->add_att("short_name", "time bounds");
			varTimeBnds->add_att("units", strTimeUnits.c_str());
			varTimeBnds->add_att("calendar", strTimeCalendar.c_str());
		}

		CopyNcVar(
			ncfilein,
			ncfileout,
			strLongitudeName);

		CopyNcVar(
			ncfilein,
			ncfileout,
			strLatitudeName);

		NcDim * dimLongitude = ncfileout.get_dim(strLongitudeName.c_str());
		if (dimLongitude == NULL) {
			_EXCEPTION1("Unable to get dimension \"%s\" from output file", strLongitudeName.c_str());
		}
		
		NcDim * dimLatitude = ncfileout.get_dim(strLatitudeName.c_str());
		if (dimLatitude == NULL) {
			_EXCEPTION1("Unable to get dimension \"%s\" from output file", strLatitudeName.c_str());
		}

		lLatitudeCount = dimLatitude->size();
		lLongitudeCount = dimLongitude->size();

		varOut =
			ncfileout.add_var(
				strVariableOutName.c_str(),
				ncFloat,
				dimTimeOut,
				dimLatitude,
				dimLongitude);

		if (varOut == NULL) {
			_EXCEPTION1("Unable to add variable \"%s\" to output file", strVariableOutName.c_str());
		}

		NcAtt * attLongName = varIn->get_att("long_name");
		if (attLongName != NULL) {
			varOut->add_att("long_name", attLongName->as_string(0));
		}
		NcAtt * attShortName = varIn->get_att("short_name");
		if (attShortName != NULL) {
			varOut->add_att("short_name", attShortName->as_string(0));
		}
		NcAtt * attUnits = varIn->get_att("units");
		if (attUnits != NULL) {
			varOut->add_att("units", attUnits->as_string(0));
		}
		
		varOut->add_att("grid_specification", "0.25 degree x 0.25 degree from 90N to 90S and 0E to 359.75E (721 x 1440 Latitude/Longitude)");
		varOut->add_att("rda_dataset", "ds633.0");
		varOut->add_att("rda_dataset_url", "https:/rda.ucar.edu/datasets/ds633.0/");
		varOut->add_att("rda_dataset_doi", "DOI: 10.5065/BH6N-5N20");
		varOut->add_att("rda_dataset_group", "ERA5 atmospheric surface forecast (mean rates or fluxes) [netCDF4]");
		varOut->add_att("number_of_significant_digits", 7);
		varOut->add_att("time_frequency", strAccumFrequency.c_str());
		varOut->add_att("computed_using", strVariableName.c_str());

		// Add provenance information
		const std::time_t timetNow = std::time(nullptr);
		std::string strProvenance = std::asctime(std::localtime(&timetNow));
		strProvenance += ": " + strCommandLine;

		ncfileout.add_att("history", strProvenance.c_str());
	}
	AnnounceEndBlock("Done");

	// Accumulate data
	AnnounceStartBlock("Accumulating data");
	long lTimeIx = 0;

	for (int f = 0; f < vecInputFiles.size(); f++) {
		NcFile ncfilein(vecInputFiles[f].c_str(), NcFile::ReadOnly);
		if (!ncfilein.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"", vecInputFiles[f].c_str());
		}

		NcVar * varTimeInitial = ncfilein.get_var("forecast_initial_time");
		if (varTimeInitial == NULL) {
			_EXCEPTION1("Cannot find variable \"forecast_initial_time\" in \"%s\"", vecInputFiles[f].c_str());
		}

		NcVar * varTimeHour = ncfilein.get_var("forecast_hour");
		if (varTimeHour == NULL) {
			_EXCEPTION1("Cannot find variable \"forecast_hour\" in \"%s\"", vecInputFiles[f].c_str());
		}

		NcVar * varIn = ncfilein.get_var(strVariableName.c_str());
		if (varIn == NULL) {
			_EXCEPTION2("Input file \"%s\" missing variable \"%s\"",
				strInputData.c_str(), strVariableName.c_str());
		}
		if (varIn->num_dims() != 4) {
			_EXCEPTION2("Expected variable \"%s\" to have 4 dimensions (%li found)",
				varIn->name(), varIn->num_dims());
		}

		DataArray1D<float> dData(lLongitudeCount * lLatitudeCount);
		DataArray1D<float> dAccum(lLongitudeCount * lLatitudeCount);

		long lForecastInitSize = varTimeInitial->get_dim(0)->size();
		long lForecastHourSize = varTimeHour->get_dim(0)->size();

		// Accumulate data
		for (long lFI = 0; lFI < lForecastInitSize; lFI++) {
			Time timeInitial(eCalType);

			int nTimeInitial;
			varTimeInitial->set_cur(lFI);
			varTimeInitial->get(&nTimeInitial, 1);

			timeInitial.FromCFCompliantUnitsOffsetInt(strTimeUnits, nTimeInitial);

			for (long lFH = 0; lFH < lForecastHourSize; lFH += lAccumFrequency) {

				// Determine if the end of this forecast is in the right year/month
				Time timeForecast = timeInitial;
				timeForecast.AddHours(lFH + lAccumFrequency);

				if ((nInputYear != (-1)) && (timeForecast.GetYear() != nInputYear)) {
					continue;
				}
				if ((nInputMonth != (-1)) && (timeForecast.GetMonth() != nInputMonth)) {
					continue;
				}

				// Accumulate data
				dAccum.Zero();
				AnnounceStartBlock("%s", timeForecast.ToString().c_str());
				for (long lFHsub = 0; lFHsub < lAccumFrequency; lFHsub++) {
					Announce("File %i FI %li FH %li", f, lFI, lFH + lFHsub);

					varIn->set_cur(lFI, lFH + lFHsub, 0, 0);
					varIn->get(&(dData[0]), 1, 1, lLatitudeCount, lLongitudeCount);

					for (int i = 0; i < dData.GetRows(); i++) {
						dAccum[i] += dData[i];
					}
				}
				AnnounceEndBlock(NULL);

				// Write data
				if (lTimeIx >= varOut->get_dim(0)->size()) {
					_EXCEPTIONT("Time index exceeds output");
				}
				varOut->set_cur(lTimeIx, 0, 0);
				varOut->put(&(dAccum[0]), 1, lLatitudeCount, lLongitudeCount);
				lTimeIx++;
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

