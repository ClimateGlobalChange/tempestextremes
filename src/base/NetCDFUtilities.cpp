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

#include "Defines.h"
#include "NetCDFUtilities.h"
#include "Exception.h"
#include "Announce.h"
#include "DataArray1D.h"
#include "netcdfcpp.h"
#include "Units.h"

#include <cstring>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

size_t NcGetVarFromList(
	NcFile & ncFile,
	const std::vector<std::string> & vecVarNames,
	NcVar ** pvar,
	NcDim ** pdim
) {
	for (size_t v = 0; v < vecVarNames.size(); v++) {
		NcVar * var = ncFile.get_var(vecVarNames[v].c_str());
		if (var != NULL) {
			if (pvar != NULL) {
				(*pvar) = var;
			}

			if (pdim != NULL) {
				NcDim * dim = ncFile.get_dim(vecVarNames[v].c_str());
				if (dim != NULL) {
					(*pdim) = dim;
				} else {
					(*pdim) = NULL;
				}
			}
			return v;
		}
	}
	return vecVarNames.size();
}

////////////////////////////////////////////////////////////////////////////////

void NcGetLatitudeLongitudeName(
	NcFile & ncFile,
	std::string & strLatitudeName,
	std::string & strLongitudeName
) {
	// Check possible latitude names
	std::vector<std::string> vecLatitudeNames;
	if (strLatitudeName == std::string("[auto]")) {
		vecLatitudeNames.push_back("lat");
		vecLatitudeNames.push_back("latitude");
		vecLatitudeNames.push_back("LAT");
		vecLatitudeNames.push_back("latitude0");
		vecLatitudeNames.push_back("Latitude");
		vecLatitudeNames.push_back("XLAT");
	} else {
		vecLatitudeNames.push_back(strLatitudeName);
	}

	// Check possible longitude names
	std::vector<std::string> vecLongitudeNames;
	if (strLatitudeName == std::string("[auto]")) {
		vecLongitudeNames.push_back("lon");
		vecLongitudeNames.push_back("longitude");
		vecLongitudeNames.push_back("LON");
		vecLongitudeNames.push_back("longitude0");
		vecLongitudeNames.push_back("Longitude");
		vecLongitudeNames.push_back("XLONG");
	} else {
		vecLongitudeNames.push_back(strLongitudeName);
	}
}

////////////////////////////////////////////////////////////////////////////////

bool NcIsTimeDimension(
	NcDim * dim
) {
	if (strcmp(dim->name(), "time") == 0) {
		return true;
	}
	if (strcmp(dim->name(), "Time") == 0) {
		return true;
	}
	if (strcmp(dim->name(), "xtime") == 0) {
		return true;
	}
	if (strcmp(dim->name(), "initial_time0_hours") == 0) {
		return true;
	}
	if (strcmp(dim->name(), "valid_time") == 0) {
		return true;
	}
	if (strcmp(dim->name(), "day") == 0) {
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////

NcDim * NcGetTimeDimension(
	NcFile & ncFile
) {
	NcDim * dim = ncFile.get_dim("time");
	if (dim != NULL) {
		return dim;
	}

	dim = ncFile.get_dim("Time");
	if (dim != NULL) {
		return dim;
	}

	dim = ncFile.get_dim("xtime");
	if (dim != NULL) {
		return dim;
	}

	dim = ncFile.get_dim("initial_time0_hours");
	if (dim != NULL) {
		return dim;
	}

	dim = ncFile.get_dim("valid_time");
	if (dim != NULL) {
		return dim;
	}

	dim = ncFile.get_dim("day");
	if (dim != NULL) {
		return dim;
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////////////

NcVar * NcGetTimeVariable(
	NcFile & ncFile
) {
	NcVar * var = ncFile.get_var("time");
	if (var != NULL) {
		return var;
	}

	var = ncFile.get_var("Time");
	if (var != NULL) {
		return var;
	}

	var = ncFile.get_var("xtime");
	if (var != NULL) {
		return var;
	}

	var = ncFile.get_var("initial_time0_hours");
	if (var != NULL) {
		return var;
	}

	var = ncFile.get_var("valid_time");
	if (var != NULL) {
		return var;
	}

	var = ncFile.get_var("day");
	if (var != NULL) {
		return var;
	}

	return NULL;
}

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

		} else if (att->type() == ncUShort) {
			fileOut->add_att(att->name(), num_vals,
				(const ushort*)(pValues->base()));

		} else if (att->type() == ncUInt) {
			fileOut->add_att(att->name(), num_vals,
				(const uint*)(pValues->base()));

		} else if (att->type() == ncInt64) {
			fileOut->add_att(att->name(), num_vals,
				(const ncint64*)(pValues->base()));

		} else if (att->type() == ncUInt64) {
			fileOut->add_att(att->name(), num_vals,
				(const ncuint64*)(pValues->base()));

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
	std::vector<std::string> vecDoNotCopyNames;
	CopyNcVarAttributes(varIn, varOut, vecDoNotCopyNames);
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVarAttributes(
	NcVar * varIn,
	NcVar * varOut,
	const std::vector<std::string> & vecDoNotCopyNames
) {
	for (int a = 0; a < varIn->num_atts(); a++) {
		NcAtt * att = varIn->get_att(a);
		long num_vals = att->num_vals();
		std::string strAttName = att->name();

		// Check names against list of attributes not to copy
		bool fSkipAttribute = false;
		for (int b = 0; b < vecDoNotCopyNames.size(); b++) {
			if (strAttName == vecDoNotCopyNames[b]) {
				fSkipAttribute = true;
			}
		}
		if (fSkipAttribute) {
			continue;
		}

		// Do not copy over add_offset or scale_factor
		if (strAttName == "add_offset") {
			continue;
		}
		if (strAttName == "scale_factor") {
			continue;
		}

		// Sometimes we can have strings of length zero.  It seems that this
		// check isn't actually necessary to prevent segfaults, so it's
		// commented out for now.
		//_ASSERT(num_vals > 0);

		NcValues * pValues = att->values();
		if (pValues == NULL) {
			_EXCEPTION2("Invalid attribute type \"%s::%s\"",
				varIn->name(), att->name());
		}

		// Do not copy over _FillValue if input/output variable types are different
		if (strAttName == "_FillValue") {
			if (varIn->type() != varOut->type()) {
				continue;
			}
		}

		// Otherwise copy over attributes
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

		} else if (att->type() == ncUShort) {
			varOut->add_att(att->name(), num_vals,
				(const ushort*)(pValues->base()));

		} else if (att->type() == ncUInt) {
			varOut->add_att(att->name(), num_vals,
				(const uint*)(pValues->base()));

		} else if (att->type() == ncInt64) {
			varOut->add_att(att->name(), num_vals,
				(const ncint64*)(pValues->base()));

		} else if (att->type() == ncUInt64) {
			varOut->add_att(att->name(), num_vals,
				(const ncuint64*)(pValues->base()));

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

NcDim * AddNcDimOrUseExisting(
	NcFile & ncFile,
	const std::string & strDimName,
	long lDimSize
) {
	NcDim * dim = ncFile.get_dim(strDimName.c_str());
	if (dim == NULL) {
		dim = ncFile.add_dim(strDimName.c_str(), lDimSize);
		if (dim == NULL) {
			_EXCEPTION2("Error adding dimension \"%s\" (%li) to file",
				strDimName.c_str(), lDimSize);
		}
	} else if (dim->size() != lDimSize) {
		_EXCEPTION3("Attempting to redefine dimension \"%s\" from size %li to %li",
			strDimName.c_str(), dim->size(), lDimSize);
	}
	return dim;
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
				if (!NcIsTimeDimension(dimA)) {
					_EXCEPTION2("Mismatch between input file dimension \"%s\" and "
						"output file dimension (%i / UNLIMITED)",
						dimA->name(), dimA->size());
				}
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

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
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

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
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

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
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

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
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

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		if (fCopyData) {
			DataArray1D<double> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncUShort type
	if (var->type() == ncUShort) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		if (fCopyData) {
			DataArray1D<ushort> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncInt type
	if (var->type() == ncUInt) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		if (fCopyData) {
			DataArray1D<uint> data(nDataSize);
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

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		if (fCopyData) {
			DataArray1D<ncint64> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// ncUInt64 type
	if (var->type() == ncUInt64) {
		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		if (fCopyData) {
			DataArray1D<ncuint64> data(nDataSize);
			var->get(&(data[0]), &(counts[0]));
			varOut->put(&(data[0]), &(counts[0]));
		}
	}

	// Check output variable exists
	if (varOut == NULL) {
		_EXCEPTION1("Unable to create output variable \"%s\"",
			var->name());
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

void CopyNcVarTimeSubset(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	const std::vector<Time> & vecOutputTimes
) {
	if (!ncIn.is_valid()) {
		_EXCEPTIONT("Invalid input file specified");
	}
	if (!ncOut.is_valid()) {
		_EXCEPTIONT("Invalid output file specified");
	}

	NcDim * dimIn = ncIn.get_dim(strVarName.c_str());
	if (dimIn == NULL) {
		_EXCEPTION1("NetCDF file does not contain dimension \"%s\"",
			strVarName.c_str());
	}

	// Check if this variable is the record dimension
	bool fRecDim = false;
	NcDim * dimRec = ncIn.rec_dim();
	if ((dimRec != NULL) && (strVarName == std::string(dimRec->name())))  {
		fRecDim = true;
	}

	NcVar * varIn = ncIn.get_var(strVarName.c_str());
	if (varIn == NULL) {
		_EXCEPTION1("NetCDF file does not contain variable \"%s\"",
			strVarName.c_str());
	}

	NcDim * dimOut;
	if (fRecDim) {
		dimOut = ncOut.add_dim(strVarName.c_str(), 0);
		if (dimOut == NULL) {
			_EXCEPTION1("Error creating dimension \"%s\" (0) in output file",
				strVarName.c_str());
		}

	} else {
		dimOut = ncOut.add_dim(strVarName.c_str(), vecOutputTimes.size());
		if (dimOut == NULL) {
			_EXCEPTION2("Error creating dimension \"%s\" (%lu) in output file",
				strVarName.c_str(), vecOutputTimes.size());
		}
	}

	NcVar * varOut = ncOut.add_var(strVarName.c_str(), varIn->type(), dimOut);
	if (varOut == NULL) {
		_EXCEPTION1("Error creating variable \"%s\" in output file",
			strVarName.c_str());
	}

	CopyNcVarAttributes(varIn, varOut);

	NcAtt * attTimeUnits = varOut->get_att("units");
	if (attTimeUnits == NULL) {
		_EXCEPTIONT("Variable \"time\" is missing \"units\" attribute");
	}
	std::string strFormattedTime = attTimeUnits->as_string(0);

	// ncInt type
	if (varOut->type() == ncInt) {
		DataArray1D<int> data(vecOutputTimes.size());
		for (int t = 0; t < vecOutputTimes.size(); t++) {
			data[t] = static_cast<int>(vecOutputTimes[t].GetCFCompliantUnitsOffsetDouble(strFormattedTime));
		}
		varOut->put(&(data[0]), data.GetRows());

	// ncFloat type
	} else if (varOut->type() == ncFloat) {
		DataArray1D<float> data(vecOutputTimes.size());
		for (int t = 0; t < vecOutputTimes.size(); t++) {
			data[t] = static_cast<float>(vecOutputTimes[t].GetCFCompliantUnitsOffsetDouble(strFormattedTime));
		}
		varOut->put(&(data[0]), data.GetRows());

	// ncDouble type
	} else if (varOut->type() == ncDouble) {
		DataArray1D<double> data(vecOutputTimes.size());
		for (int t = 0; t < vecOutputTimes.size(); t++) {
			data[t] = vecOutputTimes[t].GetCFCompliantUnitsOffsetDouble(strFormattedTime);
		}
		varOut->put(&(data[0]), data.GetRows());

	// ncInt64 type
	} else if (varOut->type() == ncInt64) {
		DataArray1D<ncint64> data(vecOutputTimes.size());
		for (int t = 0; t < vecOutputTimes.size(); t++) {
			data[t] = static_cast<ncint64>(vecOutputTimes[t].GetCFCompliantUnitsOffsetDouble(strFormattedTime));
		}
		varOut->put(&(data[0]), data.GetRows());

	// Invalid time type
	} else {
		_EXCEPTION1("Invalid time type (%i)", varOut->type());
	}

}

////////////////////////////////////////////////////////////////////////////////

void ReadCFTimeDataFromNcFile(
	NcFile * ncfile,
	const std::string & strFilename,
	NcTimeDimension & vecTimes,
	bool fWarnOnMissingCalendar
) {
	_ASSERT(ncfile != NULL);

	// Empty existing Time vector
	vecTimes.clear();

	// Get time dimension and  variable
	long lTimeCount = 0;
	NcDim * dimTime = NcGetTimeDimension(*ncfile);
	NcVar * varTime = NcGetTimeVariable(*ncfile);

	if (varTime == NULL) {
		_EXCEPTION1("Variable \"time\" not found in file \"%s\"",
			strFilename.c_str());
	}
	if (dimTime == NULL) {
		if (varTime->num_dims() != 0) {
			_EXCEPTION1("Dimension \"time\" not found in file \"%s\"",
				strFilename.c_str());
		}
		lTimeCount = 1;

	} else {
		if (varTime->num_dims() != 1) {
			_EXCEPTION1("Variable \"time\" has more than one dimension in file \"%s\"",
				strFilename.c_str());
		}
		if (!NcIsTimeDimension(varTime->get_dim(0))) {
			_EXCEPTION1("Variable \"time\" does not have time dimension in file \"%s\"",
				strFilename.c_str());
		}
		lTimeCount = dimTime->size();
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

	// Type attribute
	NcAtt * attTimeType = varTime->get_att("type");
	if (attTimeType != NULL) {
		std::string strTimeType = attTimeType->as_string(0);
		if (strTimeType.find("climatology") != std::string::npos) {
			if (strTimeType.find("daily") != std::string::npos) {
				vecTimes.m_dimtype = NcTimeDimension::TimeDimType_DailyMean;
			} else if (strTimeType.find("monthly") != std::string::npos) {
				vecTimes.m_dimtype = NcTimeDimension::TimeDimType_MonthlyMean;
			} else if (strTimeType.find("seasonal") != std::string::npos) {
				vecTimes.m_dimtype = NcTimeDimension::TimeDimType_SeasonalMean;
			} else if (strTimeType.find("annual") != std::string::npos) {
				vecTimes.m_dimtype = NcTimeDimension::TimeDimType_AnnualMean;
			}
		}
	}

	// Store details of time dimension
	vecTimes.m_nctype = varTime->type();
	vecTimes.m_units = attTimeUnits->as_string(0);

	// Load in time data
	DataArray1D<int> vecTimeInt;
	DataArray1D<float> vecTimeFloat;
	DataArray1D<double> vecTimeDouble;
	DataArray1D<ncint64> vecTimeInt64;

	if (varTime->type() == ncInt) {
		vecTimeInt.Allocate(lTimeCount);
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeInt[0]), lTimeCount);

	} else if (varTime->type() == ncFloat) {
		vecTimeFloat.Allocate(lTimeCount);
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeFloat[0]), lTimeCount);

	} else if (varTime->type() == ncDouble) {
		vecTimeDouble.Allocate(lTimeCount);
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeDouble[0]), lTimeCount);

	} else if (varTime->type() == ncInt64) {
		vecTimeInt64.Allocate(lTimeCount);
		varTime->set_cur((long)0);
		varTime->get(&(vecTimeInt64[0]), lTimeCount);

	} else {
		_EXCEPTION1("Variable \"time\" has invalid type "
			"(expected \"int\", \"int64\", \"float\" or \"double\")"
			" in file \"%s\"", strFilename.c_str());
	}

	for (int t = 0; t < lTimeCount; t++) {
		Time time(eCalendarType);
		if (varTime->type() == ncInt) {
			time.FromCFCompliantUnitsOffsetInt(
				vecTimes.m_units,
				vecTimeInt[t]);

		} else if (varTime->type() == ncFloat) {
			time.FromCFCompliantUnitsOffsetDouble(
				vecTimes.m_units,
				static_cast<double>(vecTimeFloat[t]));

		} else if (varTime->type() == ncDouble) {
			time.FromCFCompliantUnitsOffsetDouble(
				vecTimes.m_units,
				vecTimeDouble[t]);

		} else if (varTime->type() == ncInt64) {
			time.FromCFCompliantUnitsOffsetInt(
				vecTimes.m_units,
				(int)(vecTimeInt64[t]));

		}

#if defined(ROUND_TIMES_TO_NEAREST_MINUTE)
		time.RoundToNearestMinute();
#endif

		vecTimes.push_back(time);
	}
}

////////////////////////////////////////////////////////////////////////////////

void WriteCFTimeDataToNcFile(
	NcFile * ncfile,
	const std::string & strFilename,
	NcTimeDimension & vecTimes,
	bool fRecordDim
) {
	_ASSERT(ncfile != NULL);

	if (vecTimes.size() == 0) {
		_EXCEPTIONT("NcTimeDimension has zero size");
	}

	NcDim * dimTime = NcGetTimeDimension(*ncfile);
	if (dimTime == NULL) {
		if (fRecordDim) {
			dimTime = ncfile->add_dim("time", 0);
		} else {
			dimTime = ncfile->add_dim("time", vecTimes.size());
		}
		if (dimTime == NULL) {
			_EXCEPTION1("Unable to create dimension \"time\" in file \"%s\"",
				strFilename.c_str());
		}
	} else {
		if (dimTime->size() != vecTimes.size()) {
			_EXCEPTION3("File \"%s\" already contains dimension \"time\" with unexpected size (%li != %li)",
				strFilename.c_str(), dimTime->size(), vecTimes.size());
		}
	}

	NcVar * varTime = ncfile->add_var("time", vecTimes.nctype(), dimTime);
	if (varTime == NULL) {
		NcError ncerr;
		_EXCEPTION3("Unable to create variable \"time\" in file \"%s\" (%i: %s)",
			strFilename.c_str(), ncerr.get_err(), ncerr.get_errmsg());
	}

	if (vecTimes.nctype() == ncInt) {
		DataArray1D<int> nTimes(vecTimes.size());
		for (int t = 0; t < nTimes.GetRows(); t++) {
			nTimes[t] =
				static_cast<int>(
					vecTimes[t].GetCFCompliantUnitsOffsetDouble(vecTimes.units()));
		}
		varTime->put(&(nTimes[0]), nTimes.GetRows());

	} else if (vecTimes.nctype() == ncFloat) {
		DataArray1D<float> dTimes(vecTimes.size());
		for (int t = 0; t < dTimes.GetRows(); t++) {
			dTimes[t] =
				static_cast<float>(
					vecTimes[t].GetCFCompliantUnitsOffsetDouble(vecTimes.units()));
		}
		varTime->put(&(dTimes[0]), dTimes.GetRows());

	} else if (vecTimes.nctype() == ncDouble) {
		DataArray1D<double> dTimes(vecTimes.size());
		for (int t = 0; t < dTimes.GetRows(); t++) {
			dTimes[t] =
				vecTimes[t].GetCFCompliantUnitsOffsetDouble(vecTimes.units());
		}
		varTime->put(&(dTimes[0]), dTimes.GetRows());

	} else if (vecTimes.nctype() == ncInt64) {
		DataArray1D<ncint64> nTimes(vecTimes.size());
		for (int t = 0; t < nTimes.GetRows(); t++) {
			nTimes[t] =
				static_cast<ncint64>(
					vecTimes[t].GetCFCompliantUnitsOffsetDouble(vecTimes.units()));
		}
		varTime->put(&(nTimes[0]), nTimes.GetRows());

	} else {
		_EXCEPTION1("Invalid \"time\" type (%i)", varTime->type());
	}

	varTime->add_att("long_name", "time");
	varTime->add_att("calendar", vecTimes[0].GetCalendarName().c_str());
	varTime->add_att("units", vecTimes.units().c_str());

	if (vecTimes.dimtype() != NcTimeDimension::TimeDimType_Standard) {
		if (vecTimes.dimtype() == NcTimeDimension::TimeDimType_DailyMean) {
			varTime->add_att("type", "daily climatology");
		} else if (vecTimes.dimtype() == NcTimeDimension::TimeDimType_MonthlyMean) {
			varTime->add_att("type", "monthly climatology");
		} else if (vecTimes.dimtype() == NcTimeDimension::TimeDimType_SeasonalMean) {
			varTime->add_att("type", "seasonal climatology");
		} else if (vecTimes.dimtype() == NcTimeDimension::TimeDimType_AnnualMean) {
			varTime->add_att("type", "annual climatology");
		}
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

