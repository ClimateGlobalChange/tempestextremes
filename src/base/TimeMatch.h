///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimeMatch.h
///	\author  Paul Ullrich
///	\version August 25, 2020
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

#ifndef _TIMEMATCH_H_
#define _TIMEMATCH_H_

#include "Exception.h"

#include <string>

#ifndef TEMPEST_NOREGEX
#include <regex>
#endif

///////////////////////////////////////////////////////////////////////////////

#ifndef TEMPEST_NOREGEX
void TestRegex() {
	std::string strExp("(00|06):..");
	std::string strString1("00:55");
	std::string strString2("01:55");

	try {
		std::regex reTest;
		reTest.assign(strExp);

		std::smatch match1;
		if (!std::regex_search(strString1, match1, reTest)) {
			_EXCEPTIONT("std::regex failure; std::regex does not appear to work on this system");
		}

		std::smatch match2;
		if (std::regex_search(strString2, match2, reTest)) {
			_EXCEPTIONT("std::regex failure; std::regex does not appear to work on this system");
		}

	} catch(std::regex_error & reerr) {
		_EXCEPTION1("std::regex failure; std::regex does not appear to work on this system (%s)",
			reerr.what());
	}
}
#endif

///////////////////////////////////////////////////////////////////////////////

#endif
