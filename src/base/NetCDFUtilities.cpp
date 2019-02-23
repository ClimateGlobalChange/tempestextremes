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
#include "DataVector.h"
#include "netcdfcpp.h"

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

	// ncByte / ncChar type
	if ((var->type() == ncByte) || (var->type() == ncChar)) {
		DataVector<char> data;
		data.Initialize(nDataSize);

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
		DataVector<short> data;
		data.Initialize(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncInt type
	if (var->type() == ncInt) {
		DataVector<int> data;
		data.Initialize(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncFloat type
	if (var->type() == ncFloat) {
		DataVector<float> data;
		data.Initialize(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncDouble type
	if (var->type() == ncDouble) {
		DataVector<double> data;
		data.Initialize(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncInt64 type
	if (var->type() == ncInt64) {
		DataVector<ncint64> data;
		data.Initialize(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyData) {
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
	DataVector<int> vecTimeInt;
	DataVector<float> vecTimeFloat;
	DataVector<double> vecTimeDouble;
	DataVector<ncint64> vecTimeInt64;

	if (varTime->type() == ncInt) {
		vecTimeInt.Initialize(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeInt[0]), dimTime->size());

	} else if (varTime->type() == ncFloat) {
		vecTimeFloat.Initialize(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeFloat[0]), dimTime->size());

	} else if (varTime->type() == ncDouble) {
		vecTimeDouble.Initialize(dimTime->size());
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeDouble[0]), dimTime->size());

	} else if (varTime->type() == ncInt64) {
		vecTimeInt64.Initialize(dimTime->size());
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

