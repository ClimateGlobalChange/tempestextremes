///////////////////////////////////////////////////////////////////////////////
///
///	\file    CalculationList.h
///	\author  Paul Ullrich
///	\version May 25, 2025

#ifndef _CALCULATIONLIST_H_
#define _CALCULATIONLIST_H_

#include "Exception.h"

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A single calculation command, including lhs, rhs and arguments.
///	</summary>
class CalculationCommand {

public:
	///	<summary>
	///		Populate the CalculationCommand from a string.
	///	</summary>
	void Parse(
		const std::string & strCalculateCommand
	);

	///	<summary>
	///		Get a string representation of this command.
	///	</summary>
	std::string ToString() const;

public:
	///	<summary>
	///		Variable receiving the output.
	///	</summary>
	std::string lhs;

	///	<summary>
	///		Command to be executed.
	///	</summary>
	std::string rhs;

	///	<summary>
	///		Command arguments.
	///	</summary>
	std::vector<std::string> arg;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A list of CalculationCommands.
///	</summary>
class CalculationList {

public:
	///	<summary>
	///		Parse a list of calculate commands.
	///	</summary>
	void Parse(
		const std::string & strCalculateCommand
	);

	///	<summary>
	///		Number of CalculationCommands.
	///	</summary>
	size_t size() const {
		return m_vecCommands.size();
	}

	///	<summary>
	///		Get the specified CalculationCommand.
	///	</summary>
	const CalculationCommand & operator[](size_t s) const {
		if (s >= m_vecCommands.size()) {
			_EXCEPTION2("Index out of range %lu >= %lu)", s, m_vecCommands.size());
		}
		return m_vecCommands[s];
	}

private:
	///	<summary>
	///		Vector of CalculationCommands.
	///	</summary>
	std::vector<CalculationCommand> m_vecCommands;
};

///////////////////////////////////////////////////////////////////////////////

#endif // _CALCULATIONLIST_H_

