///////////////////////////////////////////////////////////////////////////////
///
///	\file    CalculationList.cpp
///	\author  Paul Ullrich
///	\version May 25, 2025

#include "CalculationList.h"
#include "STLStringHelper.h"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////

void CalculationCommand::Parse(
	const std::string & strCalculateCommand
) {
	if (strCalculateCommand.length() < 4) {
		_EXCEPTION1("Malformed calculate command \"%s\"", strCalculateCommand.c_str());
	}

	// Find the equal sign and populate LHS
	size_t sPos = 0;
	if (!std::isalpha(strCalculateCommand[0])) {
		_EXCEPTION1("LHS target variable must start with an alphabetic "
			"character in calculate command \"%s\"", strCalculateCommand.c_str());
	}
	bool fWhitespace = false;
	for (; sPos < strCalculateCommand.length(); sPos++) {
		if (strCalculateCommand[sPos] == ' ') {
			if (!fWhitespace) {
				fWhitespace = true;
			}
			continue;
		}
		if (strCalculateCommand[sPos] == '=') {
			break;
		}
		if (fWhitespace) {
			_EXCEPTION1("Invalid LHS target variable name in calculate "
				   "command \"%s\"", strCalculateCommand.c_str());
		}
		if (!std::isalnum(strCalculateCommand[sPos])) {
			_EXCEPTION1("LHS target variable must only contain alphanumeric "
				"characters in calculate command \"%s\"", strCalculateCommand.c_str());
		}
		lhs += strCalculateCommand[sPos];
	}
	if (sPos == strCalculateCommand.length()) {
		_EXCEPTION1("Missing = in calculate command \"%s\"", strCalculateCommand.c_str());
	}
	if (sPos == strCalculateCommand.length()-1) {
		_EXCEPTION1("Missing RHS in calculate command \"%s\"", strCalculateCommand.c_str());
	}

	// Get command
	sPos++;
	size_t sCommandStartPos = sPos;
	for (; sPos < strCalculateCommand.length(); sPos++) {
		if ((strCalculateCommand[sPos] == '\"') ||
		    (strCalculateCommand[sPos] == '\'') ||
			(strCalculateCommand[sPos] == ')')
		) {
			_EXCEPTION2("RHS command name may not contain %c character \"%s\"",
				strCalculateCommand[sPos], strCalculateCommand.c_str());
		}
		if (strCalculateCommand[sPos] == '(') {
			if (sPos == sCommandStartPos) {
				_EXCEPTION1("Missing RHS command name in calculate command \"%s\"", strCalculateCommand.c_str());
			}
			rhs = strCalculateCommand.substr(sCommandStartPos, sPos - sCommandStartPos);
			STLStringHelper::RemoveWhitespaceInPlace(rhs);
			if (rhs.length() == 0) {
				_EXCEPTION1("Missing RHS command name in calculate command \"%s\"", strCalculateCommand.c_str());
			}
			break;
		}
	}

	// Straight assignment operation; no arguments
	if (strCalculateCommand[sPos] != '(') {
		rhs = strCalculateCommand.substr(sCommandStartPos, sPos - sCommandStartPos);
		STLStringHelper::RemoveWhitespaceInPlace(rhs);

	// Parse arguments
	} else {
		int nParenthesisLevel = 1;
		sPos++;
		size_t sArgBegin = sPos;
		bool fInQuote = false;
		bool fFinalParenthesis = false;
		for (; sPos < strCalculateCommand.length(); sPos++) {
			if (fFinalParenthesis && (strCalculateCommand[sPos] != ' ')) {
				_EXCEPTION1("Additional characters after arguments in calculate command \"%s\"", strCalculateCommand.c_str());
			}
			if (strCalculateCommand[sPos] == '\"') {
				fInQuote = !fInQuote;
				continue;
			}
			if (fInQuote) {
				continue;
			}
			if (strCalculateCommand[sPos] == '(') {
				nParenthesisLevel++;
				continue;
			}
			if (strCalculateCommand[sPos] == ')') {
				if (nParenthesisLevel == 0) {
					_EXCEPTION1("Unbalanced parentheses in calculate command \"%s\"", strCalculateCommand.c_str());
				}
				nParenthesisLevel--;
				if (nParenthesisLevel == 0) {
					arg.push_back(strCalculateCommand.substr(sArgBegin, sPos - sArgBegin));
					STLStringHelper::RemoveWhitespaceInPlace(arg[arg.size()-1]);
					fFinalParenthesis = true;
				}
				continue;
			}
			if ((nParenthesisLevel == 1) && (strCalculateCommand[sPos] == ',')) {
				arg.push_back(strCalculateCommand.substr(sArgBegin, sPos - sArgBegin));
				STLStringHelper::RemoveWhitespaceInPlace(arg[arg.size()-1]);
				sArgBegin = sPos + 1;
				continue;
			}
		}
	}
/*
	std::cout << "LHS: " << lhs << std::endl;
	std::cout << "RHS: " << rhs << std::endl;
	if (arg.size() != 0) {
		std::cout << "ARG: [";
		for (size_t a = 0; a < arg.size(); a++) {
			std::cout << arg[a];
			if (a != arg.size()-1) {
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

std::string CalculationCommand::ToString() const {
	std::string strCommand = lhs + "=" + rhs;
	if (arg.size() != 0) {
		strCommand += "(";
		for (size_t a = 0; a < arg.size(); a++) {
			strCommand += arg[a];
			if (a != arg.size()-1) {
				strCommand += ",";
			}
		}
		strCommand += ")";
	}
	return strCommand;
}

///////////////////////////////////////////////////////////////////////////////

void CalculationList::Parse(
	const std::string & strCalculateCommand
) {

	// Vector storing semi-colon positions
	std::vector<size_t> vecCommandPos;
	vecCommandPos.push_back(0);

	// Identify occurrence of semi-colons in string, avoiding quoted semi-colons
	size_t sLast = 0;
	bool fInQuote = false;

	for (size_t s = 0; s < strCalculateCommand.length(); s++) {
		if (s == strCalculateCommand.length()) {
			vecCommandPos.push_back(s+1);
			break;
		}
		if (strCalculateCommand[s] == '\"') {
			fInQuote = !fInQuote;
			continue;
		}
		if (!fInQuote && (strCalculateCommand[s] == ';')) {
			vecCommandPos.push_back(s+1);
			sLast = s+1;
		}
	}
	vecCommandPos.push_back(strCalculateCommand.length()+1);

	// Parse each CalculateCommand
	m_vecCommands.resize(vecCommandPos.size()-1);
	for (size_t s = 0; s < m_vecCommands.size(); s++) {
		m_vecCommands[s].Parse(
			strCalculateCommand.substr(
				vecCommandPos[s],
				vecCommandPos[s+1] - vecCommandPos[s] - 1));
	}
}

///////////////////////////////////////////////////////////////////////////////

