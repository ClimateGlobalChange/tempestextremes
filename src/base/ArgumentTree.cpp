///////////////////////////////////////////////////////////////////////////////
///
///	\file    ArgumentTree.cpp
///	\author  Paul Ullrich
///	\version September 7, 2018
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "ArgumentTree.h"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////

ArgumentTree::~ArgumentTree() {
	for (int i = 0; i < m_vecSubArguments.size(); i++) {
		if (m_vecSubArguments[i] != NULL) {
			delete m_vecSubArguments[i];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ArgumentTree::SubParse(
	const std::string & str,
	char szListType
) {
	if (str.length() == 0) {
		return;
	}

	m_strArgumentList.push_back(str);

	ArgumentTree * pSubTree = new ArgumentTree(false);
	if (pSubTree == NULL) {
		_EXCEPTION();
	}

	pSubTree->m_szListType = szListType;

	bool fHasSubArguments =
		pSubTree->Parse(
			m_strArgumentList[m_strArgumentList.size()-1]);

	if (fHasSubArguments) {
		m_vecSubArguments.push_back(pSubTree);
	} else {
		m_vecSubArguments.push_back(NULL);
		delete pSubTree;
	}
}

///////////////////////////////////////////////////////////////////////////////

bool ArgumentTree::Parse(
	const std::string & str
) {
	//std::cout << str << std::endl;
	if (m_strArgumentList.size() != 0) {
		_EXCEPTIONT("Attempting to initialize already initialized ArgumentTree");
	}

	// Check length
	if (str.length() == 0) {
		m_strArgumentList.push_back("");
		m_vecSubArguments.push_back(NULL);
		return false;
	}

	// Last character where something interesting happened
	int iLast = 0;

	// Search for semi-colon separated strings at top level
	if (m_fSemiColonDelimited) {
		for (size_t i = 0; i <= str.length(); i++) {

			// Semi-colon separator
			if ((i == str.length()) || (str[i] == ';')) {
				SubParse(str.substr(iLast,i-iLast), ';');
				iLast = i+1;
			}
		}

		return true;
	}

	// Parse delimeters
	bool fHasSubArguments = false;

	for (size_t i = 0; i <= str.length(); i++) {

		// End-of-string
		if (i == str.length()) {
			if (iLast != 0) {
				SubParse(str.substr(iLast,i-iLast), '\0');
			}
			break;
		}

		// Comma separator
		if (str[i] == ',') {
			SubParse(str.substr(iLast,i-iLast), '\0');
			iLast = i+1;
			continue;
		}

		// Equal sign
		if (str[i] == '=') {
			if (m_szListType != '(') {
				SubParse(str.substr(iLast,i-iLast), '\0');
				m_strArgumentList.push_back("=");
				m_vecSubArguments.push_back(NULL);
				iLast = i+1;
				continue;
			}
		}

		// Open parentheses
		if (str[i] == '(') {
			SubParse(str.substr(iLast,i-iLast), '\0');
			iLast = i+1;

			// Find matching parentheses
			int iLevel = 1;
			for (i++; i < str.length(); i++) {
				if (str[i] == '(') {
					iLevel++;
				}
				if (str[i] == ')') {
					iLevel--;
					if (iLevel == 0) {
						SubParse(str.substr(iLast,i-iLast), '(');
						iLast = i+1;
						break;
					}
				}
			}
			if (iLevel != 0) {
				_EXCEPTION1("Input string \"%s\" missing closed parentheses \')\'",
					str.c_str());
			}
			continue;
		}

		// Close parentheses (unbalanced)
		if (str[i] == ')') {
			_EXCEPTION1("Input string \"%s\" contains unbalanced parentheses \')\'",
				str.c_str());
		}

		// Semi-colon (not allowed)
		if (str[i] == ';') {
			_EXCEPTION1("Input string \"%s\" contains semi-colons",
				str.c_str());
		}
	}

	return (m_strArgumentList.size() > 1);
}

///////////////////////////////////////////////////////////////////////////////

