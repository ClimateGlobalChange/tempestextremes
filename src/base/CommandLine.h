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

///	<summary>
///		Types of parameters.
///	</summary>
enum ParameterType {
	ParameterTypeNone,
	ParameterTypeBool,
	ParameterTypeString,
	ParameterTypeInt,
	ParameterTypeDouble,
	ParameterTypeTime
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line parameter.
///	</summary>
class CommandLineParameter {
public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	CommandLineParameter(
		std::string strName,
		std::string strDescription
	) :
		m_strName(std::string("--") + strName),
		m_strDescription(strDescription)
	{ }

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~CommandLineParameter() {
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ParameterType GetParameterType() {
		return ParameterTypeNone;
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
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line boolean.
///	</summary>
class CommandLineParameterBool : public CommandLineParameter {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineParameterBool(
		bool & ref,
		std::string strName,
		std::string strDescription
	) :
		CommandLineParameter(strName, strDescription),
		m_fValue(ref)
	{
		m_fValue = false;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ParameterType GetParameterType() {
		return ParameterTypeBool;
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
		Announce("  %s <bool> [%s] %s",
			m_strName.c_str(),
			(m_fValue)?("true"):("false"),
			m_strDescription.c_str());
	}

	///	<summary>
	///		Activate this parameter.
	///	</summary>
	virtual void Activate() {
		m_fValue = true;
	}

public:
	///	<summary>
	///		Parameter value.
	///	</summary>
	bool & m_fValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line string.
///	</summary>
class CommandLineParameterString : public CommandLineParameter {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineParameterString(
		std::string & ref,
		std::string strName,
	 	std::string strDefaultValue,
		std::string strDescription
	) :
		CommandLineParameter(strName, strDescription),
		m_strValue(ref)
	{
		m_strValue = strDefaultValue;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ParameterType GetParameterType() {
		return ParameterTypeString;
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
		Announce("  %s <string> [\"%s\"] %s",
			m_strName.c_str(),
			m_strValue.c_str(),
			m_strDescription.c_str());
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
	///		Parameter value.
	///	</summary>
	std::string & m_strValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line integer.
///	</summary>
class CommandLineParameterInt : public CommandLineParameter {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineParameterInt(
		int & ref,
		std::string strName,
		int dDefaultValue,
		std::string strDescription
	) :
		CommandLineParameter(strName, strDescription),
		m_dValue(ref)
	{
		m_dValue = dDefaultValue;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ParameterType GetParameterType() const {
		return ParameterTypeInt;
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
		Announce("  %s <integer> [%i] %s",
			m_strName.c_str(),
			m_dValue,
			m_strDescription.c_str());
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
	///		Parameter value.
	///	</summary>
	int & m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line double.
///	</summary>
class CommandLineParameterDouble : public CommandLineParameter {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineParameterDouble(
		double & ref,
		std::string strName,
		double dDefaultValue,
		std::string strDescription
	) :
		CommandLineParameter(strName, strDescription),
		m_dValue(ref)
	{
		m_dValue = dDefaultValue;
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ParameterType GetParameterType() const {
		return ParameterTypeDouble;
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
	///		Parameter value.
	///	</summary>
	double & m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A command line Time
///	</summary>
class CommandLineParameterTime : public CommandLineParameter {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CommandLineParameterTime(
		Time & ref,
		std::string strName,
	 	std::string strDefaultValue,
		std::string strDescription,
		Time::TimeType eTimeType = Time::TypeFixed
	) :
		CommandLineParameter(strName, strDescription),
		m_timeValue(ref),
		m_eTimeType(eTimeType)
	{
		m_timeValue.FromFormattedString(strDefaultValue);
		m_timeValue.SetTimeType(m_eTimeType);
	}

	///	<summary>
	///		Identify the type of parameter.
	///	</summary>
	virtual ParameterType GetParameterType() {
		return ParameterTypeTime;
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
	///		Parameter value.
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
	  std::vector<CommandLineParameter*> _vecParameters;

///	<summary>
///		Define a new command line boolean parameter.
///	</summary>
#define CommandLineBool(ref, name) \
	_vecParameters.push_back( \
		new CommandLineParameterBool(ref, name, ""));

///	<summary>
///		Define a new command line boolean parameter with description.
///	</summary>
#define CommandLineBoolD(ref, name, desc) \
	_vecParameters.push_back( \
		new CommandLineParameterBool(ref, name, desc));

///	<summary>
///		Define a new command line string parameter.
///	</summary>
#define CommandLineString(ref, name, value) \
	_vecParameters.push_back( \
		new CommandLineParameterString(ref, name, value, ""));

///	<summary>
///		Define a new command line string parameter with description.
///	</summary>
#define CommandLineStringD(ref, name, value, desc) \
	_vecParameters.push_back( \
		new CommandLineParameterString(ref, name, value, desc));

///	<summary>
///		Define a new command line integer parameter.
///	</summary>
#define CommandLineInt(ref, name, value) \
	_vecParameters.push_back( \
		new CommandLineParameterInt(ref, name, value, ""));

///	<summary>
///		Define a new command line integer parameter with description.
///	</summary>
#define CommandLineIntD(ref, name, value, desc) \
	_vecParameters.push_back( \
		new CommandLineParameterInt(ref, name, value, desc));

///	<summary>
///		Define a new command line double parameter.
///	</summary>
#define CommandLineDouble(ref, name, value) \
	_vecParameters.push_back( \
		new CommandLineParameterDouble(ref, name, value, ""));

///	<summary>
///		Define a new command line double parameter with description.
///	</summary>
#define CommandLineDoubleD(ref, name, value, desc) \
	_vecParameters.push_back( \
		new CommandLineParameterDouble(ref, name, value, desc));

///	<summary>
///		Define a new command line Time parameter.
///	</summary>
#define CommandLineFixedTime(ref, name, value) \
	_vecParameters.push_back( \
		new CommandLineParameterTime(ref, name, value, "", Time::TypeFixed));

///	<summary>
///		Define a new command line Time parameter with description.
///	</summary>
#define CommandLineFixedTimeD(ref, name, value, desc) \
	_vecParameters.push_back( \
		new CommandLineParameterTime(ref, name, value, desc, Time::TypeFixed));

///	<summary>
///		Define a new command line Time parameter.
///	</summary>
#define CommandLineDeltaTime(ref, name, value) \
	_vecParameters.push_back( \
		new CommandLineParameterTime(ref, name, value, "", Time::TypeDelta));

///	<summary>
///		Define a new command line Time parameter with description.
///	</summary>
#define CommandLineDeltaTimeD(ref, name, value, desc) \
	_vecParameters.push_back( \
		new CommandLineParameterTime(ref, name, value, desc, Time::TypeDelta));

///	<summary>
///		Begin the loop for command line parameters.
///	</summary>
#define ParseCommandLine(argc, argv) \
    for(int _command = 1; _command < argc; _command++) { \
		bool _found = false; \
		for(int _p = 0; _p < _vecParameters.size(); _p++) { \
			if (_vecParameters[_p]->m_strName == argv[_command]) { \
				_found = true; \
				_vecParameters[_p]->Activate(); \
				int _z; \
				if (_vecParameters[_p]->GetValueCount() >= argc - _command) { \
					Announce("Error: Insufficient values for option %s", \
						argv[_command]); \
					_errorCommandLine = true; \
					_command = argc; \
					break; \
				} \
				for (_z = 0; _z < _vecParameters[_p]->GetValueCount(); _z++) { \
					if ((_command + _z + 1 < argc) && \
						(strlen(argv[_command + _z + 1]) > 2) && \
						(argv[_command + _z + 1][0] == '-') && \
						(argv[_command + _z + 1][1] == '-') \
					) break; \
					_command++; \
					_vecParameters[_p]->SetValue(_z, argv[_command]); \
				} \
				if ((_vecParameters[_p]->GetValueCount() >= 0) && \
					(_z != _vecParameters[_p]->GetValueCount()) \
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
			Announce("Warning: Invalid parameter \"%s\"", argv[_command]); \
		} \
	}

///	<summary>
///		Print usage information.
///	</summary>
#define PrintCommandLineUsage(argv) \
	if (_errorCommandLine) \
		Announce("\nUsage: %s <Parameter List>", argv[0]); \
	Announce("Parameters:"); \
	for (int _p = 0; _p < _vecParameters.size(); _p++) \
		_vecParameters[_p]->PrintUsage(); \
	if (_errorCommandLine) \
		exit(-1);

///	<summary>
///		End the definition of command line parameters.
///	</summary>
#define EndCommandLine(argv) \
	PrintCommandLineUsage(argv); \
	for (int _p = 0; _p < _vecParameters.size(); _p++) \
		delete _vecParameters[_p]; \
	}

///////////////////////////////////////////////////////////////////////////////

#endif
