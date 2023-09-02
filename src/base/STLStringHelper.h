///////////////////////////////////////////////////////////////////////////////
///
///	\file    STLStringHelper.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _STLSTRINGHELPER_H_
#define _STLSTRINGHELPER_H_

#include "Units.h"

#include <string>
#include <vector>

#include <cstring>

///	<summary>
///		This class exposes additional functionality which can be used to
///		supplement the STL string class.
///	</summary>
class STLStringHelper {

///////////////////////////////////////////////////////////////////////////////

private:
STLStringHelper() { }

public:

///////////////////////////////////////////////////////////////////////////////

inline static void ToLower(std::string &str) {
	for(size_t i = 0; i < str.length(); i++) {
		str[i] = tolower(str[i]);
	}
}

///////////////////////////////////////////////////////////////////////////////

inline static void ToUpper(std::string &str) {
	for(size_t i = 0; i < str.length(); i++) {
		str[i] = toupper(str[i]);
	}
}

///////////////////////////////////////////////////////////////////////////////

inline static bool IsIntegerIndex(const std::string &str) {
	if (str.length() == 0) {
		return false;
	}
	for(size_t i = 0; i < str.length(); i++) {
		if ((str[i] < '0') || (str[i] > '9')) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

inline static bool IsInteger(const std::string &str) {
	if (str.length() == 0) {
		return false;
	}
	for(size_t i = 0; i < str.length(); i++) {
		if ((i == 0) && ((str[i] == '-') || (str[i] == '+'))) {
			if (str.length() == 1) {
				return false;
			}
			continue;
		}
		if ((str[i] < '0') || (str[i] > '9')) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

inline static bool IsFloat(const std::string &str) {
	bool fIsFloat = false;
	bool fHasExponent = false;
	bool fHasDecimal = false;
	for(size_t i = 0; i < str.length(); i++) {
		if ((str[i] < '0') || (str[i] > '9')) {
			if (str[i] == '.') {
				if (fHasDecimal) {
					return false;
				}
				if (fHasExponent) {
					return false;
				}
				fHasDecimal = true;
				continue;
			}
			if (str[i] == 'e') {
				if (fHasExponent) {
					return false;
				}
				fHasExponent = true;
				continue;
			}
			if ((str[i] == '-') || (str[i] == '+')) {
				if (i == 0) {
					continue;
				} else if (str[i-1] == 'e') {
					continue;
				} else {
					return false;
				}
			}
			if (str[i] == 'f') {
				if (i != str.length()-1) {
					return false;
				}
			}
			return false;

		} else {
			fIsFloat = true;
		}
	}
	return fIsFloat;
}

///////////////////////////////////////////////////////////////////////////////

static void RemoveWhitespaceInPlace(
	std::string & strString
) {
	size_t sBegin = strString.length();
	for (size_t s = 0; s < strString.length(); s++) {
		if ((strString[s] != ' ') && (strString[s] != '\t')) {
			sBegin = s;
			break;
		}
	}
	if (sBegin == strString.length()) {
		strString = "";
		return;
	}

	size_t sEnd = strString.length();
	for (size_t s = sEnd-1; s > sBegin; s--) {
		if ((strString[s] != ' ') && (strString[s] != '\t')) {
			sEnd = s+1;
			break;
		}
	}

	strString = strString.substr(sBegin, sEnd - sBegin);
}

///////////////////////////////////////////////////////////////////////////////

static void ParseVariableList(
	const std::string & strVariables,
	std::vector< std::string > & vecVariableStrings,
	const std::string & strDelimiters = std::string(" ,;")
) {
	int iVarBegin = 0;
	int iVarCurrent = 0;

	int nParenthesesLevel = 0;

	if (strVariables.length() == 0) {
		return;
	}

	for (;;) {
		if ((iVarCurrent >= strVariables.length()) ||
		    (strDelimiters.find(strVariables[iVarCurrent]) != std::string::npos)
		) {
			if (nParenthesesLevel == 0) {
				_ASSERT(iVarCurrent > iVarBegin);
				if (iVarCurrent == iVarBegin) {
					_EXCEPTION1("Zero length variable string in \"%s\"",
						strVariables.c_str());
				}

				vecVariableStrings.push_back(
					strVariables.substr(iVarBegin, iVarCurrent - iVarBegin));

				iVarBegin = iVarCurrent + 1;
			}
			if (iVarCurrent >= strVariables.length()) {
				break;
			}
		}
		if (strVariables[iVarCurrent] == '(') {
			nParenthesesLevel++;
		}
		if (strVariables[iVarCurrent] == ')') {
			nParenthesesLevel--;
			if (nParenthesesLevel < 0) {
				_EXCEPTIONT("Unmatched closed parenthesis in variable list");
			}
		}

		iVarCurrent++;
	}

	if (nParenthesesLevel != 0) {
		_EXCEPTIONT("Unmatched open parenthesis in variable list");
	}
}


///////////////////////////////////////////////////////////////////////////////

static void SplitVariableStrings(
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

static std::string ConcatenateStringVector(
	const std::vector< std::string > & vecStrings,
	std::string strDelimiter
) {
	std::string strConcat;
	for (int v = 0; v < vecStrings.size(); v++) {
		strConcat += vecStrings[v];
		if (v != vecStrings.size() - 1) {
			strConcat += strDelimiter;
		}
	}
	return strConcat;
}

///////////////////////////////////////////////////////////////////////////////

};

#endif

