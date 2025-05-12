///////////////////////////////////////////////////////////////////////////////
///
///	\file    CommandLine.h
///	\author  Paul Ullrich
///	\version February 24, 2013
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

#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_

#include "TimeObj.h"
#include "Announce.h"
#include "Exception.h"

#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////
//
// If set then invalid command line arguments will only throw warnings.
//
//#define WARN_ONLY_ON_INVALID_ARGUMENT

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Types of parameters.
///	</summary>
enum ArgumentType {
	ArgumentTypeNone,
	ArgumentTypeBool,
	ArgumentTypeString,
	ArgumentTypeInt,
	ArgumentTypeDouble,
	ArgumentTypeTime
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line parameter.
///	</summary>
class CommandLineArgument {
public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	CommandLineArgument(
		std::string strName,
		std::string strDescription
	) :
		m_strName(std::string("--") + strName),
		m_strDescription(strDescription),
		m_fHidden(false)
	{
		if ((strName.length() > 0) && (strName[0] == '*')) {
			m_strName = std::string("--") + strName.substr(1);
			m_fHidden = true;
		}
	}

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~CommandLineArgument() {
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ArgumentType GetArgumentType() {
		return ArgumentTypeNone;
	}

	///	<summary>
	///		Number of values required.
	///	</summary>
	virtual int GetValueCount() const = 0;

	///	<summary>
	///		Print the usage information of this parameter.
	///	</summary>
	virtual void PrintUsage() const = 0;

	///	<summary>
	///		Activate this parameter.
	///	</summary>
	virtual void Activate() {
	}

	///	<summary>
	///		Set the value from a string.
	///	</summary>
	virtual void SetValue(
		int ix,
		std::string strValue
	) {
		_EXCEPTIONT("Invalid value index.");
	}

public:
	///	<summary>
	///		Name of this parameter.
	///	</summary>
	std::string m_strName;

	///	<summary>
	///		Description of this parameter.
	///	</summary>
	std::string m_strDescription;

	///	<summary>
	///		Flag indicating this argument is hidden.
	///	</summary>
	bool m_fHidden;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line boolean.
///	</summary>
class CommandLineArgumentBool : public CommandLineArgument {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineArgumentBool(
		bool & ref,
		std::string strName,
		std::string strDescription
	) :
		CommandLineArgument(strName, strDescription),
		m_fValue(ref)
	{
		m_fValue = false;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ArgumentType GetArgumentType() {
		return ArgumentTypeBool;
	}

	///	<summary>
	///		Number of values required.
	///	</summary>
	virtual int GetValueCount() const {
		return (0);
	}

	///	<summary>
	///		Print the usage information of this parameter.
	///	</summary>
	virtual void PrintUsage() const {
		if (!m_fHidden) {
			Announce("  %s <bool> [%s] %s",
				m_strName.c_str(),
				(m_fValue)?("true"):("false"),
				m_strDescription.c_str());
		}
	}

	///	<summary>
	///		Activate this parameter.
	///	</summary>
	virtual void Activate() {
		m_fValue = true;
	}

public:
	///	<summary>
	///		Argument value.
	///	</summary>
	bool & m_fValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line string.
///	</summary>
class CommandLineArgumentString : public CommandLineArgument {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineArgumentString(
		std::string & ref,
		std::string strName,
	 	std::string strDefaultValue,
		std::string strDescription
	) :
		CommandLineArgument(strName, strDescription),
		m_strValue(ref)
	{
		m_strValue = strDefaultValue;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ArgumentType GetArgumentType() {
		return ArgumentTypeString;
	}

	///	<summary>
	///		Number of values required.
	///	</summary>
	virtual int GetValueCount() const {
		return (1);
	}

	///	<summary>
	///		Print the usage information of this parameter.
	///	</summary>
	virtual void PrintUsage() const {
		if (!m_fHidden) {
			Announce("  %s <string> [\"%s\"] %s",
				m_strName.c_str(),
				m_strValue.c_str(),
				m_strDescription.c_str());
		}
	}

	///	<summary>
	///		Set the value from a string.
	///	</summary>
	virtual void SetValue(
		int ix,
		std::string strValue
	) {
		if (ix != 0) {
			_EXCEPTIONT("Invalid value index.");
		}
		m_strValue = strValue.c_str();
	}

public:
	///	<summary>
	///		Argument value.
	///	</summary>
	std::string & m_strValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line integer.
///	</summary>
class CommandLineArgumentInt : public CommandLineArgument {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineArgumentInt(
		int & ref,
		std::string strName,
		int dDefaultValue,
		std::string strDescription
	) :
		CommandLineArgument(strName, strDescription),
		m_dValue(ref)
	{
		m_dValue = dDefaultValue;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ArgumentType GetArgumentType() const {
		return ArgumentTypeInt;
	}

	///	<summary>
	///		Number of values required.
	///	</summary>
	virtual int GetValueCount() const {
		return (1);
	}

	///	<summary>
	///		Print the usage information of this parameter.
	///	</summary>
	virtual void PrintUsage() const {
		if (!m_fHidden) {
			Announce("  %s <integer> [%i] %s",
				m_strName.c_str(),
				m_dValue,
				m_strDescription.c_str());
		}
	}

	///	<summary>
	///		Set the value from a string.
	///	</summary>
	virtual void SetValue(
		int ix,
		std::string strValue
	) {
		if (ix != 0) {
			_EXCEPTIONT("Invalid value index.");
		}
		m_dValue = atoi(strValue.c_str());
	}

public:
	///	<summary>
	///		Argument value.
	///	</summary>
	int & m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line double.
///	</summary>
class CommandLineArgumentDouble : public CommandLineArgument {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineArgumentDouble(
		double & ref,
		std::string strName,
		double dDefaultValue,
		std::string strDescription
	) :
		CommandLineArgument(strName, strDescription),
		m_dValue(ref)
	{
		m_dValue = dDefaultValue;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ArgumentType GetArgumentType() const {
		return ArgumentTypeDouble;
	}

	///	<summary>
	///		Number of values required.
	///	</summary>
	virtual int GetValueCount() const {
		return (1);
	}

	///	<summary>
	///		Print the usage information of this parameter.
	///	</summary>
	virtual void PrintUsage() const {
		if (!m_fHidden) {
			if (fabs(m_dValue) < 1.0e6) {
				Announce("  %s <double> [%f] %s",
					m_strName.c_str(),
					m_dValue,
					m_strDescription.c_str());
			} else {
				Announce("  %s <double> [%e] %s",
					m_strName.c_str(),
					m_dValue,
					m_strDescription.c_str());
			}
		}
	}

	///	<summary>
	///		Set the value from a string.
	///	</summary>
	virtual void SetValue(
		int ix,
		std::string strValue
	) {
		if (ix != 0) {
			_EXCEPTIONT("Invalid value index.");
		}
		m_dValue = atof(strValue.c_str());
	}

public:
	///	<summary>
	///		Argument value.
	///	</summary>
	double & m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line Time
///	</summary>
class CommandLineArgumentTime : public CommandLineArgument {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineArgumentTime(
		Time & ref,
		std::string strName,
	 	std::string strDefaultValue,
		std::string strDescription,
		Time::TimeType eTimeType = Time::TypeFixed
	) :
		CommandLineArgument(strName, strDescription),
		m_timeValue(ref),
		m_eTimeType(eTimeType)
	{
		m_timeValue.FromFormattedString(strDefaultValue);
		m_timeValue.SetTimeType(m_eTimeType);
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ArgumentType GetArgumentType() {
		return ArgumentTypeTime;
	}

	///	<summary>
	///		Number of values required.
	///	</summary>
	virtual int GetValueCount() const {
		return (1);
	}

	///	<summary>
	///		Print the usage information of this parameter.
	///	</summary>
	virtual void PrintUsage() const {
		if (!m_fHidden) {
			if (m_eTimeType == Time::TypeFixed) {
				Announce("  %s <time> [%s] %s",
					m_strName.c_str(),
					m_timeValue.ToFreeString().c_str(),
					m_strDescription.c_str());

			} else if (m_eTimeType == Time::TypeDelta) {
				Announce("  %s <dtime> [%s] %s",
					m_strName.c_str(),
					m_timeValue.ToFreeString().c_str(),
					m_strDescription.c_str());

			} else {
				_EXCEPTIONT("Invalid TimeType");
			}
		}
	}

	///	<summary>
	///		Set the value from a string.
	///	</summary>
	virtual void SetValue(
		int ix,
		std::string strValue
	) {
		if (ix != 0) {
			_EXCEPTIONT("Invalid value index.");
		}
		m_timeValue.FromFormattedString(strValue);
		m_timeValue.SetTimeType(m_eTimeType);
	}

public:
	///	<summary>
	///		Argument value.
	///	</summary>
	Time & m_timeValue;

	///	<summary>
	///		Time type.
	///	</summary>
	Time::TimeType m_eTimeType;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Begin the definition of command line parameters.
///	</summary>
#define BeginCommandLine() \
	{ bool _errorCommandLine = false; \
	  bool _invalidArgument = false; \
	  std::vector<CommandLineArgument*> _vecArguments;

///	<summary>
///		Define a new command line boolean parameter.
///	</summary>
#define CommandLineBool(ref, name) \
	_vecArguments.push_back( \
		new CommandLineArgumentBool(ref, name, ""));

///	<summary>
///		Define a new command line boolean parameter with description.
///	</summary>
#define CommandLineBoolD(ref, name, desc) \
	_vecArguments.push_back( \
		new CommandLineArgumentBool(ref, name, desc));

///	<summary>
///		Define a new command line string parameter.
///	</summary>
#define CommandLineString(ref, name, value) \
	_vecArguments.push_back( \
		new CommandLineArgumentString(ref, name, value, ""));

///	<summary>
///		Define a new command line string parameter with description.
///	</summary>
#define CommandLineStringD(ref, name, value, desc) \
	_vecArguments.push_back( \
		new CommandLineArgumentString(ref, name, value, desc));

///	<summary>
///		Define a new command line integer parameter.
///	</summary>
#define CommandLineInt(ref, name, value) \
	_vecArguments.push_back( \
		new CommandLineArgumentInt(ref, name, value, ""));

///	<summary>
///		Define a new command line integer parameter with description.
///	</summary>
#define CommandLineIntD(ref, name, value, desc) \
	_vecArguments.push_back( \
		new CommandLineArgumentInt(ref, name, value, desc));

///	<summary>
///		Define a new command line double parameter.
///	</summary>
#define CommandLineDouble(ref, name, value) \
	_vecArguments.push_back( \
		new CommandLineArgumentDouble(ref, name, value, ""));

///	<summary>
///		Define a new command line double parameter with description.
///	</summary>
#define CommandLineDoubleD(ref, name, value, desc) \
	_vecArguments.push_back( \
		new CommandLineArgumentDouble(ref, name, value, desc));

///	<summary>
///		Define a new command line Time parameter.
///	</summary>
#define CommandLineFixedTime(ref, name, value) \
	_vecArguments.push_back( \
		new CommandLineArgumentTime(ref, name, value, "", Time::TypeFixed));

///	<summary>
///		Define a new command line Time parameter with description.
///	</summary>
#define CommandLineFixedTimeD(ref, name, value, desc) \
	_vecArguments.push_back( \
		new CommandLineArgumentTime(ref, name, value, desc, Time::TypeFixed));

///	<summary>
///		Define a new command line Time parameter.
///	</summary>
#define CommandLineDeltaTime(ref, name, value) \
	_vecArguments.push_back( \
		new CommandLineArgumentTime(ref, name, value, "", Time::TypeDelta));

///	<summary>
///		Define a new command line Time parameter with description.
///	</summary>
#define CommandLineDeltaTimeD(ref, name, value, desc) \
	_vecArguments.push_back( \
		new CommandLineArgumentTime(ref, name, value, desc, Time::TypeDelta));

///	<summary>
///		Begin the loop for command line parameters.
///	</summary>
#ifdef WARN_ONLY_ON_INVALID_ARGUMENT
#define INVALID_ARGUMENT_STRING "WARNING: Invalid argument \"%s\""
#else
#define INVALID_ARGUMENT_STRING "ERROR: Invalid argument \"%s\""
#endif

#define ParseCommandLine(argc, argv) \
    for(int _command = 1; _command < argc; _command++) { \
		bool _found = false; \
		for(int _p = 0; _p < _vecArguments.size(); _p++) { \
			if (_vecArguments[_p]->m_strName == argv[_command]) { \
				_found = true; \
				_vecArguments[_p]->Activate(); \
				int _z; \
				if (_vecArguments[_p]->GetValueCount() >= argc - _command) { \
					Announce("Error: Insufficient values for option %s", \
						argv[_command]); \
					_errorCommandLine = true; \
					_command = argc; \
					break; \
				} \
				for (_z = 0; _z < _vecArguments[_p]->GetValueCount(); _z++) { \
					if ((_command + _z + 1 < argc) && \
						(strlen(argv[_command + _z + 1]) > 2) && \
						(argv[_command + _z + 1][0] == '-') && \
						(argv[_command + _z + 1][1] == '-') \
					) break; \
					_command++; \
					_vecArguments[_p]->SetValue(_z, argv[_command]); \
				} \
				if ((_vecArguments[_p]->GetValueCount() >= 0) && \
					(_z != _vecArguments[_p]->GetValueCount()) \
				) { \
					Announce("Error: Insufficient values for option %s", \
						argv[_command]); \
					_errorCommandLine = true; \
					_command = argc; \
					break; \
				} \
			} \
		} \
		if (!_found) { \
			_invalidArgument = true; \
			Announce(INVALID_ARGUMENT_STRING, argv[_command]); \
		} \
	}

///	<summary>
///		Print usage information.
///	</summary>
#ifdef WARN_ONLY_ON_INVALID_ARGUMENT
#define PrintCommandLineUsage(argv) \
	if (_errorCommandLine) \
		Announce("\nUsage: %s <Argument List>", argv[0]); \
	Announce("Arguments:"); \
	for (int _p = 0; _p < _vecArguments.size(); _p++) \
		_vecArguments[_p]->PrintUsage(); \
	if (_errorCommandLine) \
		exit(-1);
#else
#define PrintCommandLineUsage(argv) \
	if ((_errorCommandLine) || (_invalidArgument)) \
		Announce("\nUsage: %s <Argument List>", argv[0]); \
	Announce("Arguments:"); \
	for (int _p = 0; _p < _vecArguments.size(); _p++) \
		_vecArguments[_p]->PrintUsage(); \
	if ((_errorCommandLine) || (_invalidArgument)) \
		exit(-1);
#endif

///	<summary>
///		End the definition of command line parameters.
///	</summary>
#define EndCommandLine(argv) \
		PrintCommandLineUsage(argv); \
		for (int _p = 0; _p < _vecArguments.size(); _p++) { \
			delete _vecArguments[_p]; \
		} \
	}

///	<summary>
///		Concatenate the command line into a string.
///	</summary>
inline std::string GetCommandLineAsString(int argc, char ** argv) {
	std::string strCommandLine;
	for (int i = 0; i < argc; i++) {
		strCommandLine += argv[i];
		if (i != argc-1) {
			strCommandLine += " ";
		}
	}
	return strCommandLine;
}

///////////////////////////////////////////////////////////////////////////////

#endif
