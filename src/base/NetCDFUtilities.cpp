///////////////////////////////////////////////////////////////////////////////
///
///	\file    NetCDFUtilities.cpp
///	\author  Paul Ullrich
///	\version August 14, 2014
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "NetCDFUtilities.h"
#include "Exception.h"
#include "Announce.h"
#include "DataArray1D.h"
#include "netcdfcpp.h"
#include "Units.h"

#include <cstring>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

void CopyNcFileAttributes(
	NcFile * fileIn,
	NcFile * fileOut
) {
	for (int a = 0; a < fileIn->num_atts(); a++) {
		NcAtt * att = fileIn->get_att(a);
		long num_vals = att->num_vals();

		NcValues * pValues = att->values();
		if (pValues == NULL) {
			_EXCEPTION1("Invalid attribute type in file \"%s\"",
				att->name());
		}

		if (att->type() == ncByte) {
			fileOut->add_att(att->name(), num_vals,
				(const ncbyte*)(pValues->base()));

		} else if (att->type() == ncChar) {
			fileOut->add_att(att->name(), num_vals,
				(const char*)(pValues->base()));

		} else if (att->type() == ncShort) {
			fileOut->add_att(att->name(), num_vals,
				(const short*)(pValues->base()));

		} else if (att->type() == ncInt) {
			fileOut->add_att(att->name(), num_vals,
				(const int*)(pValues->base()));

		} else if (att->type() == ncFloat) {
			fileOut->add_att(att->name(), num_vals,
				(const float*)(pValues->base()));

		} else if (att->type() == ncDouble) {
			fileOut->add_att(att->name(), num_vals,
				(const double*)(pValues->base()));

		} else if (att->type() == ncString) {
			fileOut->add_att(att->name(), strlen((const char *)pValues->base()),
				(const char*)(pValues->base()));

		} else {
			_EXCEPTIONT("Invalid attribute type");
		}

		delete pValues;
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVarAttributes(
	NcVar * varIn,
	NcVar * varOut
) {
	for (int a = 0; a < varIn->num_atts(); a++) {
		NcAtt * att = varIn->get_att(a);
		long num_vals = att->num_vals();

		NcValues * pValues = att->values();
		if (pValues == NULL) {
			_EXCEPTION2("Invalid attribute type \"%s::%s\"",
				varIn->name(), att->name());
		}

		if (att->type() == ncByte) {
			varOut->add_att(att->name(), num_vals,
				(const ncbyte*)(pValues->base()));

		} else if (att->type() == ncChar) {
			varOut->add_att(att->name(), num_vals,
				(const char*)(pValues->base()));

		} else if (att->type() == ncShort) {
			varOut->add_att(att->name(), num_vals,
				(const short*)(pValues->base()));

		} else if (att->type() == ncInt) {
			varOut->add_att(att->name(), num_vals,
				(const int*)(pValues->base()));

		} else if (att->type() == ncFloat) {
			varOut->add_att(att->name(), num_vals,
				(const float*)(pValues->base()));

		} else if (att->type() == ncDouble) {
			varOut->add_att(att->name(), num_vals,
				(const double*)(pValues->base()));

		} else if (att->type() == ncString) {
			varOut->add_att(att->name(), strlen((const char *)pValues->base()),
				(const char*)(pValues->base()));

		} else {
			_EXCEPTIONT("Invalid attribute type");
		}

		delete pValues;
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVar(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes,
	bool fCopyData
) {
	if (!ncIn.is_valid()) {
		_EXCEPTIONT("Invalid input file specified");
	}
	if (!ncOut.is_valid()) {
		_EXCEPTIONT("Invalid output file specified");
	}
	NcVar * var = ncIn.get_var(strVarName.c_str());
	if (var == NULL) {
		_EXCEPTION1("NetCDF file does not contain variable \"%s\"",
			strVarName.c_str());
	}

	NcVar * varOut;

	std::vector<NcDim *> dimOut;
	dimOut.resize(var->num_dims());

	std::vector<long> counts;
	counts.resize(var->num_dims());

	long nDataSize = 1;

	for (int d = 0; d < var->num_dims(); d++) {
		NcDim * dimA = var->get_dim(d);

		dimOut[d] = ncOut.get_dim(dimA->name());

		if (dimOut[d] == NULL) {
			if (dimA->is_unlimited()) {
				dimOut[d] = ncOut.add_dim(dimA->name());
			} else {
				dimOut[d] = ncOut.add_dim(dimA->name(), dimA->size());
			}

			if (dimOut[d] == NULL) {
				_EXCEPTION2("Failed to add dimension \"%s\" (%i) to file",
					dimA->name(), dimA->size());
			}

			if (strVarName != dimA->name()) {
				CopyNcVarIfExists(
					ncIn,
					ncOut,
					dimA->name(),
					true,
					true);
			}
		}
		if (dimOut[d]->size() != dimA->size()) {
			if (dimA->is_unlimited() && !dimOut[d]->is_unlimited()) {
				_EXCEPTION2("Mismatch between input file dimension \"%s\" and "
					"output file dimension (UNLIMITED / %i)",
					dimA->name(), dimOut[d]->size());
			} else if (!dimA->is_unlimited() && dimOut[d]->is_unlimited()) {
				_EXCEPTION2("Mismatch between input file dimension \"%s\" and "
					"output file dimension (%i / UNLIMITED)",
					dimA->name(), dimA->size());
			} else if (!dimA->is_unlimited() && !dimOut[d]->is_unlimited()) {
				_EXCEPTION3("Mismatch between input file dimension \"%s\" and "
					"output file dimension (%i / %i)",
					dimA->name(), dimA->size(), dimOut[d]->size());
			}
		}

		counts[d] = dimA->size();
		nDataSize *= counts[d];
	}

	// Check for existence of variable
	if (ncOut.get_var(var->name()) != NULL) {
		_EXCEPTION1("Variable \"%s\" already exists in output file", var->name());
	}

	// ncByte / ncChar type
	if ((var->type() == ncByte) || (var->type() == ncChar)) {
		DataArray1D<char> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// ncShort type
	if (var->type() == ncShort) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			DataArray1D<short> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncInt type
	if (var->type() == ncInt) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			DataArray1D<int> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncFloat type
	if (var->type() == ncFloat) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			DataArray1D<float> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncDouble type
	if (var->type() == ncDouble) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			DataArray1D<double> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncInt64 type
	if (var->type() == ncInt64) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			DataArray1D<ncint64> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// Check output variable exists
	if (varOut == NULL) {
		_EXCEPTION1("Unable to create output variable \"%s\"",
			var->name());
	}

	// Copy attributes
	if (fCopyAttributes) {
		CopyNcVarAttributes(var, varOut);
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVarIfExists(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes,
	bool fCopyData
) {
	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	NcVar * var = ncIn.get_var(strVarName.c_str());
	if (var != NULL) {
		CopyNcVar(ncIn, ncOut, strVarName, fCopyAttributes, fCopyData);
	}
}

////////////////////////////////////////////////////////////////////////////////

void ReadCFTimeDataFromNcFile(
	NcFile * ncfile,
	const std::string & strFilename,
	std::vector<Time> & vecTimes,
	bool fWarnOnMissingCalendar
) {
	// Empty existing Time vector
	vecTimes.clear();

	// Get time dimension
	NcDim * dimTime = ncfile->get_dim("time");
	if (dimTime == NULL) {
		_EXCEPTION1("Dimension \"time\" not found in file \"%s\"",
			strFilename.c_str());
	}

	// Get time variable
	NcVar * varTime = ncfile->get_var("time");
	if (varTime == NULL) {
		_EXCEPTION1("Variable \"time\" not found in file \"%s\"",
			strFilename.c_str());
	}
	if (varTime->num_dims() != 1) {
		_EXCEPTION1("Variable \"time\" has more than one dimension in file \"%s\"",
			strFilename.c_str());
	}
	if (strcmp(varTime->get_dim(0)->name(), "time") != 0) {
		_EXCEPTION1("Variable \"time\" does not have dimension \"time\" in file \"%s\"",
			strFilename.c_str());
	}

	// Calendar attribute
	NcAtt * attTimeCal = varTime->get_att("calendar");
	std::string strCalendar;
	if (attTimeCal == NULL) {
		if (fWarnOnMissingCalendar) {
			Announce("WARNING: Variable \"time\" is missing \"calendar\" attribute; assuming \"standard\"");
		}
		strCalendar = "standard";
	} else {
		strCalendar = attTimeCal->as_string(0);
	}
	Time::CalendarType eCalendarType =
		Time::CalendarTypeFromString(strCalendar);

	// Units attribute
	NcAtt * attTimeUnits = varTime->get_att("units");
	if (attTimeUnits == NULL) {
		_EXCEPTION1("Variable \"time\" is missing \"units\" attribute in file \"%s\"",
			strFilename.c_str());
	}
	std::string strTimeUnits = attTimeUnits->as_string(0);

	// Load in time data
	DataArray1D<int> vecTimeInt;
	DataArray1D<float> vecTimeFloat;
	DataArray1D<double> vecTimeDouble;
	DataArray1D<ncint64> vecTimeInt64;

	if (varTime->type() == ncInt) {
		vecTimeInt.Allocate(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeInt[0]), dimTime->size());

	} else if (varTime->type() == ncFloat) {
		vecTimeFloat.Allocate(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeFloat[0]), dimTime->size());

	} else if (varTime->type() == ncDouble) {
		vecTimeDouble.Allocate(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeDouble[0]), dimTime->size());

	} else if (varTime->type() == ncInt64) {
		vecTimeInt64.Allocate(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeInt64[0]), dimTime->size());

	} else {
		_EXCEPTION1("Variable \"time\" has invalid type "
			"(expected \"int\", \"int64\", \"float\" or \"double\")"
			" in file \"%s\"", strFilename.c_str());
	}

	for (int t = 0; t < dimTime->size(); t++) {
		Time time(eCalendarType);
		if (varTime->type() == ncInt) {
			time.FromCFCompliantUnitsOffsetInt(
				strTimeUnits,
				vecTimeInt[t]);

		} else if (varTime->type() == ncFloat) {
			time.FromCFCompliantUnitsOffsetDouble(
				strTimeUnits,
				static_cast<double>(vecTimeFloat[t]));

		} else if (varTime->type() == ncDouble) {
			time.FromCFCompliantUnitsOffsetDouble(
				strTimeUnits,
				vecTimeDouble[t]);

		} else if (varTime->type() == ncInt64) {
			time.FromCFCompliantUnitsOffsetInt(
				strTimeUnits,
				(int)(vecTimeInt64[t]));

		}

		vecTimes.push_back(time);
	}
}

////////////////////////////////////////////////////////////////////////////////

long GetIntegerIndexFromValueBasedIndex(
	NcFile * ncfile,
	const std::string & strFilename,
	const std::string & strVariableName,
	long lDim,
	const std::string & strValueIndex
) {
	_ASSERT(ncfile != NULL);

	long lDimIndex = (-1);

	// Get the value and units from the argument string
	std::string strValue;
	std::string strUnits;
	SplitIntoValueAndUnits(strValueIndex, strValue, strUnits);

	// This is already a value-based index
	if (STLStringHelper::IsIntegerIndex(strValue) && (strUnits == "")) {
		return std::stol(strValue);
	}

	// Get the variable
	NcVar * var = ncfile->get_var(strVariableName.c_str());
	if (var == NULL) {
		_EXCEPTION2("Unable to load variable \"%s\" in file \"%s\"",
			strVariableName.c_str(), strFilename.c_str());
	}
	if ((lDim < 0) || (lDim >= var->num_dims())) {
		_EXCEPTION3("Specified dimension index (%li) exceeds number of dimensions in variable \"%s\" (%li)",
			lDim, strVariableName.c_str(), var->num_dims());
	}


	// Get dimension variable
	NcDim * dimArg = var->get_dim(lDim);
	if (dimArg == NULL) {
		_EXCEPTION3("Failure to read dimension %li of variable \"%s\" in file \"%s\"",
			lDim,
			strVariableName.c_str(),
			strFilename.c_str());
	}

	std::string strDimName = dimArg->name();

	NcVar * varDim = ncfile->get_var(strDimName.c_str());
	if (varDim == NULL) {
		_EXCEPTION3("Dimension \"%s\" in file \"%s\" has no corresponding variable (unable to use value notation for variable \"%s\")",
			dimArg->name(),
			strFilename.c_str(),
			strVariableName.c_str());
	}
	if (varDim->num_dims() != 1) {
		_EXCEPTION2("Dimension variable \"%s\" in file \"%s\" has more than one dimension",
			dimArg->name(),
			strFilename.c_str());
	}
	if (strDimName != varDim->get_dim(0)->name()) {
		_EXCEPTION3("Dimension variable \"%s\" in file \"%s\" instead refers to dimension \"%s\"",
			dimArg->name(),
			strFilename.c_str(),
			varDim->get_dim(0)->name());
	}

	NcAtt * attUnits = varDim->get_att("units");
	if (attUnits == NULL) {
		_EXCEPTION3("Dimension variable \"%s\" in file \"%s\" has no units attribute (unable to use value notation for variable \"%s\")",
			dimArg->name(),
			strFilename.c_str(),
			strVariableName.c_str());
	}

	std::string strNcDimUnits = attUnits->as_string(0);

	// Integer type dimension variable
	if (varDim->type() == ncInt) {
		if (!STLStringHelper::IsInteger(strValue)) {
			_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is not of type integer, as expected for dimension \"%s\"",
				strValueIndex.c_str(),
				strVariableName.c_str(),
				dimArg->name());
		}
		long lValue = std::stol(strValue);

		bool fConvert =
			ConvertUnits<long>(
				lValue,
				strUnits,
				strNcDimUnits);

		if (!fConvert) {
			_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is incompatible with dimension units \"%s\"",
				strValueIndex.c_str(),
				strVariableName.c_str(),
				strNcDimUnits.c_str());
		}

		DataArray1D<int> nDimValues(dimArg->size());
		varDim->get(&(nDimValues[0]), dimArg->size());

		NcError err;
		if (err.get_err() != NC_NOERR) {
			_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
		}

		long vd = 0;
		for (; vd < dimArg->size(); vd++) {
			if (nDimValues[vd] == lValue) {
				lDimIndex = vd;
				break;
			}
		}
		if (vd == dimArg->size()) {
			_EXCEPTION2("Dimension \"%s\" does not contain coordinate value \"%s\"",
				dimArg->name(),
				strValueIndex.c_str());
		}


	// double type dimension
	} else if ((varDim->type() == ncFloat) || (varDim->type() == ncDouble)) {
		if (!STLStringHelper::IsFloat(strValue)) {
			_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is not of type float/double, as expected for dimension \"%s\"",
				strValueIndex.c_str(),
				strVariableName.c_str(),
				dimArg->name());
		}
		double dValue = std::stof(strValue);

		bool fConvert =
			ConvertUnits<double>(
				dValue,
				strUnits,
				strNcDimUnits);

		if (!fConvert) {
			_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is incompatible with dimension units \"%s\"",
				strValueIndex.c_str(),
				strVariableName.c_str(),
				strNcDimUnits.c_str());
		}

		DataArray1D<double> dDimValues(dimArg->size());
		if (varDim->type() == ncFloat) {
			DataArray1D<float> flDimValues(dimArg->size());
			varDim->get(&(flDimValues[0]), dimArg->size());
			for (long vd = 0; vd < dimArg->size(); vd++) {
				dDimValues[vd] = static_cast<double>(flDimValues[vd]);
			}
		} else {
			varDim->get(&(dDimValues[0]), dimArg->size());
		}

		NcError err;
		if (err.get_err() != NC_NOERR) {
			_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
		}

		long vd = 0;
		for (; vd < dimArg->size(); vd++) {
			if (fabs(dDimValues[vd] - dValue) < 1.0e-12 * fabs(dValue)) {
				lDimIndex = vd;
				break;
			}
		}
		if (vd == dimArg->size()) {
			_EXCEPTION2("Dimension \"%s\" does not contain coordinate value \"%s\"",
				dimArg->name(),
				strValueIndex.c_str());
		}

		//printf("\"%s\" corresponds to dimension index %li", m_strArg[a].c_str(), vd);

	// Unknown dimension type
	} else {
		_EXCEPTION2("Dimension variable \"%s\" in file \"%s\" must be of type int, float, or double to use value indexing",
			varDim->name(),
			strFilename.c_str());
	}

	_ASSERT(lDimIndex != (-1));

	return lDimIndex;
}

////////////////////////////////////////////////////////////////////////////////

