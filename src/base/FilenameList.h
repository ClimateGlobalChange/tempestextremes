///////////////////////////////////////////////////////////////////////////////
///
///	\file    FilenameList.h
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

#ifndef _FILENAMELIST_H_
#define _FILENAMELIST_H_

#include "Exception.h"

#include <string>
#include <fstream>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

class FilenameList : public std::vector<std::string> {

public:
	///	<summary>
	///		Parse the filename list from a file containing a list of filenames.
	///	</summary>
	void FromFile(
		const std::string & strFileListFile,
		bool fAllowMultipleFilesPerLine = true
	) {
		std::ifstream ifFileList(strFileListFile.c_str());
		if (!ifFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strFileListFile.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			if (!fAllowMultipleFilesPerLine) {
				if (strFileLine.find(';') != std::string::npos) {
					_EXCEPTION1("Only one filename allowed per line in \"%s\"",
						strFileListFile.c_str());
				}
			}
			push_back(strFileLine);
		}
		if (size() == 0) {
			_EXCEPTION1("No filenames found in \"%s\"",
				strFileListFile.c_str());

		}
	}
};

///////////////////////////////////////////////////////////////////////////////

#endif

