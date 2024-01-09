///////////////////////////////////////////////////////////////////////////////
///
///	\file    AutoCuratorTool.cpp
///	\author  Paul Ullrich
///	\version January 9, 2024
///
///	<remarks>
///		Copyright 2024 Paul Ullrich
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
#include "AutoCurator.h"

#include <fstream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input data
	std::string strInputData;

	// Input data list
	std::string strInputDataList;

	// Add a calendar attribute
	std::string strCalendarName;

	// Output file
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strCalendarName, "add_calendar", "");
		CommandLineString(strOutputFile, "out_index", "");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check command line arguments
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output index file (--out_index) specified");
	}

	// Create autocurator
	AutoCurator autocurator;

	// Set calendar manually
	if (strCalendarName.length() != 0) {
		autocurator.SetCalendar(Time::CalendarTypeFromString(strCalendarName));
	}

	// Curate input data
	if (strInputData.length() != 0) {
		AnnounceStartBlock("Autocurating in_data");
		autocurator.IndexFiles(strInputData);

	} else if (strInputDataList.length() != 0) {
		AnnounceStartBlock("Autocurating in_data_list");
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
			Announce(strFileLine.c_str());
			autocurator.IndexFiles(strFileLine);
		}
	}
	AnnounceEndBlock("Done");

	// Write to file
	autocurator.ToYAMLFile(strOutputFile);
/*
	// Read from file
	AutoCurator autocurator2;
	autocurator2.FromYAMLFile(strOutputFile);
	autocurator2.ToYAMLFile("test.txt");
*/
	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

