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

#include "Exception.h"

#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <fstream>

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
	virtual ~CommandLineParameter()
	{ }

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
		printf("  %s <bool> [%s] %s\n",
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
		printf("  %s <string> [\"%s\"] %s\n",
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
		printf("  %s <integer> [%i] %s\n",
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
		printf("  %s <double> [%f] %s\n",
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
///		Parse command line parameters.
///	</summary>
void _ParseCommandLine(
	int argc,
	char ** argv,
	std::vector<CommandLineParameter*> & vecParameters,
	bool & errorCommandLine
) {
	for (int command = 1; command < argc; command++) {
		bool found = false;

		// Include command line arguments from file
		if ((command == 1) &&
			(strlen(argv[command]) > 0) &&
			(argv[command][0] != '-')
		) {
			std::vector<std::string> strFileCommands;
			
			std::ifstream filein(argv[command]);
			if (filein.fail()) {
				_EXCEPTION1("No command line file \"%s\" found",
					argv[command]);
			}

			for (;;) {
				std::string strInput;
				filein >> strInput;
				if (filein.eof()) {
					break;
				}
				strFileCommands.push_back(strInput);
			}

			// Build command
			if (strFileCommands.size() != 0) {
				if (strFileCommands[0][0] != '-') {
					_EXCEPTIONT("Command line file must only contain commands");
				}
				int argcin = (int)strFileCommands.size()+1;
				std::vector<char*> argvin;
				argvin.reserve(argcin);
				argvin.push_back(argv[0]);
				for (int d = 0; d < argcin; d++) {
					argvin.push_back(&(strFileCommands[d][0]));
				}

				// Recursively call this routine
				_ParseCommandLine(
					argcin, &(argvin[0]), vecParameters, errorCommandLine);

				continue;
			}
		}

		// Parse parameters
		for (int p = 0; p < vecParameters.size(); p++) {
			if (vecParameters[p]->m_strName == argv[command]) {
				found = true;
				vecParameters[p]->Activate();
				int z;
				if (vecParameters[p]->GetValueCount() >= argc - command) {
					printf("Error: Insufficient values for option %s\n",
						argv[command]);
					errorCommandLine = true;
					command = argc;
					break;
				}
				for (z = 0; z < vecParameters[p]->GetValueCount(); z++) {
					if ((command + z + 1 < argc) &&
						(strlen(argv[command + z + 1]) > 2) &&
						(argv[command + z + 1][0] == '-') &&
						(argv[command + z + 1][1] == '-')
					) break;
					command++;
					vecParameters[p]->SetValue(z, argv[command]);
				}
				if ((vecParameters[p]->GetValueCount() >= 0) &&
					(z != vecParameters[p]->GetValueCount())
				) {
					printf("Error: Insufficient values for option %s\n",
						argv[command]);
					errorCommandLine = true;
					command = argc;
					break;
				}
			}
		}
		if (!found) {
			printf("Error: Invalid parameter \"%s\"\n", argv[command]);
			errorCommandLine = true;
			break;
		}
	}
}

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
///		Begin the loop for command line parameters.
///	</summary>
#define ParseCommandLine(argc, argv) \
	_ParseCommandLine(argc, argv, _vecParameters, _errorCommandLine);
/*
    for(int _command = 1; _command < argc; _command++) { \
		bool _found = false; \
		for(int _p = 0; _p < _vecParameters.size(); _p++) { \
			if (_vecParameters[_p]->m_strName == argv[_command]) { \
				_found = true; \
				_vecParameters[_p]->Activate(); \
				int _z; \
				if (_vecParameters[_p]->GetValueCount() >= argc - _command) { \
					printf("Error: Insufficient values for option %s\n", \
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
					printf("Error: Insufficient values for option %s\n", \
						argv[_command]); \
					_errorCommandLine = true; \
					_command = argc; \
					break; \
				} \
			} \
		} \
		if (!_found) { \
			printf("Error: Invalid parameter \"%s\"\n", argv[_command]); \
			_errorCommandLine = true; \
			break; \
		} \
	}
*/
///	<summary>
///		Print usage information.
///	</summary>
#define PrintCommandLineUsage(argv) \
	if (_errorCommandLine) \
		printf("\nUsage: %s <Parameter List>\n", argv[0]); \
	printf("Parameters:\n"); \
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
