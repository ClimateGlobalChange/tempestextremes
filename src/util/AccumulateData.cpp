///////////////////////////////////////////////////////////////////////////////
///
///	\file    AccumulateData.cpp
///	\author  Paul Ullrich
///	\version March 9, 2023
///
///	<remarks>
///		Copyright 2023 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "NetCDFUtilities.h"
#include "Variable.h"
#include "DataArray1D.h"
#include "Variable.h"
#include "FilenameList.h"
#include "NetCDFUtilities.h"
#include "STLStringHelper.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>
#include <limits>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

enum AccumulationFrequency {
	AccumulationFrequency_1h,
	AccumulationFrequency_3h,
	AccumulationFrequency_6h,
	AccumulationFrequency_24h,
	AccumulationFrequency_Monthly,
	AccumulationFrequency_Annual,
	AccumulationFrequency_WaterYear
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Operators
	enum AccumulateDataOp {
		AccumulateDataOp_Sum,
		AccumulateDataOp_Avg,
		AccumulateDataOp_Min,
		AccumulateDataOp_Max
	};

	// Input data
	std::string strInputData;

	// Input data list
	std::string strInputDataList;

	// Output data
	std::string strOutputData;

	// Output data
	std::string strOutputDataList;

	// Operation
	std::string strOp;

	// Accumulation frequency
	std::string strAccumFrequency;

	// Accumulate from previous times
	bool fAccumBackward;

	// Check for missing data
	bool fMissingData;

	// Variable to use for quantile calculation
	std::string strVariable;

	// Output variable
	std::string strVariableOut;

	// Variables to preserve
	std::string strPreserveVar;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputDataList, "out_data_list", "");
		CommandLineStringD(strOp, "op", "sum", "[sum|avg|min|max]");
		CommandLineStringD(strAccumFrequency, "accumfreq", "6h", "[1h|3h|6h|daily|monthly|annual|wateryear]");
		CommandLineBool(fAccumBackward, "accumbackward");
		CommandLineBool(fMissingData, "missingdata");
		CommandLineString(strVariable, "var", "");
		CommandLineString(strVariableOut, "varout", "");
		CommandLineString(strPreserveVar, "preserve", "");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	int nInputCount =
		  ((strInputData.length() == 0)?(0):(1))
		+ ((strInputDataList.length() == 0)?(0):(1));

	if (nInputCount == 0) {
		_EXCEPTIONT("No input file (--in_data) or (--in_data_list) specified");
	}
	if (nInputCount > 1) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list) may be specified");
	}

	if ((strInputData.length() != 0) && (strOutputData.length() == 0)) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}
	if ((strInputDataList.length() != 0) && (strOutputDataList.length() == 0) && (strOutputData.length() == 0)) {
		_EXCEPTIONT("No output file (--out_data) or output data list (--out_data_list) specified");
	}

	// Parse preserve list
	std::vector<std::string> vecPreserveVariableStrings;
	STLStringHelper::ParseVariableList(strPreserveVar, vecPreserveVariableStrings);

	// Get filename lists
	FilenameList vecInputFiles;
	if (strInputData != "") {
		vecInputFiles.push_back(strInputData);
	}
	if (strInputDataList != "") {
		vecInputFiles.FromFile(strInputDataList);
	}

	FilenameList vecOutputFiles;
	if (strOutputData != "") {
		vecOutputFiles.push_back(strOutputData);
	}
	if (strOutputDataList != "") {
		vecOutputFiles.FromFile(strOutputDataList);
	}

	if ((strOutputDataList.length() != 0) && (vecInputFiles.size() != vecOutputFiles.size())) {
		_EXCEPTION2("Input data list must be same length as output data list (%lu vs %lu)",
			vecInputFiles.size(), vecOutputFiles.size());
	}

	// Operator (--op)
	AccumulateDataOp m_eAccumOp;
	STLStringHelper::ToLower(strOp);
	if (strOp == "sum") {
		m_eAccumOp = AccumulateDataOp_Sum;
	} else if (strOp == "avg") {
		m_eAccumOp = AccumulateDataOp_Avg;
	} else if (strOp == "min") {
		m_eAccumOp = AccumulateDataOp_Min;
	} else if (strOp == "max") {
		m_eAccumOp = AccumulateDataOp_Max;
	} else {
		_EXCEPTIONT("Invalid value for --op.  Expected \"sum\", \"avg\", \"min\", or \"max\"");
	}

	// Parse variable list (--var)
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVariable, vecVariableStrings);

	std::vector<std::string> vecVariableNames;
	std::vector< std::vector<std::string> > vecVariableSpecifiedDims;

	STLStringHelper::SplitVariableStrings(
		vecVariableStrings,
		vecVariableNames,
		vecVariableSpecifiedDims
	);

	if (vecVariableNames.size() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}

	// Parse output variable list (--varout)
	std::vector<std::string> vecVariableOutNames;
	if (strVariableOut == "") {
		vecVariableOutNames = vecVariableNames;
	} else {
		STLStringHelper::ParseVariableList(strVariableOut, vecVariableOutNames);
		if (vecVariableOutNames.size() != vecVariableStrings.size()) {
			_EXCEPTIONT("If specified, --varout and --var must have the same length");
		}
	}

	// Accumulation frequency (--accumfreq)
	AccumulationFrequency eAccumulationFreq;
	if (strAccumFrequency == "1h") {
		eAccumulationFreq = AccumulationFrequency_1h;
	} else if (strAccumFrequency == "3h") {
		eAccumulationFreq = AccumulationFrequency_3h;
	} else if (strAccumFrequency == "6h") {
		eAccumulationFreq = AccumulationFrequency_6h;
	} else if ((strAccumFrequency == "daily") || (strAccumFrequency == "24h")) {
		eAccumulationFreq = AccumulationFrequency_24h;
	} else if (strAccumFrequency == "monthly") {
		eAccumulationFreq = AccumulationFrequency_Monthly;
	} else if (strAccumFrequency == "annual") {
		eAccumulationFreq = AccumulationFrequency_Annual;
	} else if (strAccumFrequency == "wateryear") {
		eAccumulationFreq = AccumulationFrequency_WaterYear;
	} else {
		_EXCEPTIONT("--accumfreq must be \"1h\", \"3h\", \"6h\", \"daily\", \"monthly\", \"annual\", or \"wateryear\"");
	}

	// Data
	size_t sTotalSize = 1;
	std::vector<long> vecPos;
	std::vector<long> vecSize;

	DataArray1D<float> dataIn;
	DataArray1D<float> dataAccum;
	DataArray1D<int> dataAccumCount;

	// Number of accumulated timesteps
	int nAccumulatedTimesteps = 0;

	// Start accumulation
	AnnounceStartBlock("Start accumulation");

	for (size_t v = 0; v < vecVariableStrings.size(); v++) {

		AnnounceStartBlock("Variable \"%s\"", vecVariableStrings[v].c_str());

		sTotalSize = 1;

		vecPos.clear();
		vecSize.clear();

		if (dataAccum.IsAttached()) {
			if ((m_eAccumOp == AccumulateDataOp_Sum) ||
			    (m_eAccumOp == AccumulateDataOp_Avg)
			) {
				dataAccum.Zero();
	
			} else if (m_eAccumOp == AccumulateDataOp_Min) {
				for (size_t i = 0; i < sTotalSize; i++) {
					dataAccum[i] = std::numeric_limits<float>::max();
				}
	
			} else if (m_eAccumOp == AccumulateDataOp_Max) {
				for (size_t i = 0; i < sTotalSize; i++) {
					dataAccum[i] = -std::numeric_limits<float>::max();
				}
			}
		}

		if (dataAccumCount.IsAttached()) {
			dataAccumCount.Zero();
		}

		NcTimeDimension vecTimesOut;
		Time timeLastFileTime(Time::CalendarUnknown);
		NcFile * pncfileout = NULL;
		NcVar * varout = NULL;
		int iCurrentAccumIndex = -1;

		Time timeLastProcessed(Time::CalendarUnknown);

		for (size_t f = 0; f < vecInputFiles.size(); f++) {
	
			AnnounceStartBlock("%s", vecInputFiles[f].c_str());

			// Open input files
			NcFile ncfilein(vecInputFiles[f].c_str(), NcFile::ReadOnly);
			if (!ncfilein.is_valid()) {
				_EXCEPTION1("Unable to open input file \"%s\"", vecInputFiles[f].c_str());
			}

			// Get variable
			NcVar * var = ncfilein.get_var(vecVariableStrings[v].c_str());
			if (var == NULL) {
				_EXCEPTION2("Unable to open variable \"%s\" in file \"%s\"",
					vecVariableStrings[v].c_str(), vecInputFiles[f].c_str());
			}
			if ((var->num_dims() == 0) || !NcIsTimeDimension(var->get_dim(0))) {
				_EXCEPTION2("First dimension of variable \"%s\" in file \"%s\" must be \"time\"",
					vecVariableStrings[v].c_str(), vecInputFiles[f].c_str());
			}
	
			// Get FillValue
			float dFillValue = std::numeric_limits<float>::max();
			if (fMissingData) {
				NcAtt * attFillValue = var->get_att("_FillValue");
				if (attFillValue == NULL) {
					attFillValue = var->get_att("missing_value");
				}
				if (attFillValue != NULL) {
					dFillValue = attFillValue->as_float(0);
				}
			}

			// Get scale_factor and add_offset
			float dScaleFactor = 1.0f;
			float dAddOffset = 0.0f;
			NcAtt * attScaleFactor = var->get_att("scale_factor");
			if (attScaleFactor != NULL) {
				dScaleFactor = attScaleFactor->as_float(0);
			}
			NcAtt * attAddOffset = var->get_att("add_offset");
			if (attAddOffset != NULL) {
				dAddOffset = attAddOffset->as_float(0);
			}
	
			// Allocate space for accumulation
			if (f == 0) {
				vecSize.resize(var->num_dims());
				vecSize[0] = 1;
				for (long d = 1; d < var->num_dims(); d++) {
					vecSize[d] = var->get_dim(d)->size();
					sTotalSize *= vecSize[d];
				}
				vecPos.resize(var->num_dims(), 0);
				dataIn.Allocate(sTotalSize);
				dataAccum.Allocate(sTotalSize);
	
				if (fMissingData) {
					dataAccumCount.Allocate(sTotalSize);
				}
			}
	
			// Get times in this input file
			NcTimeDimension vecTimesIn;
			ReadCFTimeDataFromNcFile(&ncfilein, vecInputFiles[f], vecTimesIn, false);
	
			if (vecTimesIn.size() < 1) {
				_EXCEPTION1("Time variable in file \"%s\" has zero length",
					vecInputFiles[f].c_str());
			}

			// Get output times
			std::vector<int> iAccumIndex(vecTimesIn.size());
			NcTimeDimension vecTimesOutNew;
			vecTimesOutNew.m_nctype = vecTimesIn.m_nctype;
			vecTimesOutNew.m_units = vecTimesIn.m_units;
			vecTimesOutNew.m_dimtype = vecTimesIn.m_dimtype;

			if (vecOutputFiles.size() > 1) {
				vecTimesOut.clear();
				vecTimesOut.m_nctype = vecTimesIn.m_nctype;
				vecTimesOut.m_units = vecTimesIn.m_units;
				vecTimesOut.m_dimtype = vecTimesIn.m_dimtype;
			}

			for (size_t t = 0; t < vecTimesIn.size(); t++) {
				if ((t == 0) && (timeLastProcessed.GetCalendarType() != Time::CalendarUnknown)) {
					_ASSERT(f != 0);
					if (vecTimesIn[0] <= timeLastProcessed) {
						_EXCEPTION2("Time values across files \"%s\" and \"%s\" are not monotone increasing",
							vecInputFiles[f-1].c_str(), vecInputFiles[f].c_str());
					}
				}
				if ((t != 0) && (vecTimesIn[t] <= vecTimesIn[t-1])) {
					_EXCEPTION1("Time values in file \"%s\" are not monotone increasing",
						vecInputFiles[f].c_str());
				}

				Time timeRounded = vecTimesIn[t];
				if (eAccumulationFreq == AccumulationFrequency_1h) {
					timeRounded.SetSecond((timeRounded.GetSecond() / 3600) * 3600);
				} else if (eAccumulationFreq == AccumulationFrequency_3h) {
					timeRounded.SetSecond((timeRounded.GetSecond() / 10800) * 10800);
				} else if (eAccumulationFreq == AccumulationFrequency_6h) {
					timeRounded.SetSecond((timeRounded.GetSecond() / 21600) * 21600);
				} else if (eAccumulationFreq == AccumulationFrequency_24h) {
					timeRounded.SetSecond(0);
				} else if (eAccumulationFreq == AccumulationFrequency_Monthly) {
					timeRounded.SetSecond(0);
					timeRounded.SetDay(1);
				} else if (eAccumulationFreq == AccumulationFrequency_Annual) {
					timeRounded.SetSecond(0);
					timeRounded.SetDay(1);
					timeRounded.SetMonth(1);
				} else if (eAccumulationFreq == AccumulationFrequency_WaterYear) {
					if (timeRounded.GetMonth() >= 10) {
						timeRounded.SetSecond(0);
						timeRounded.SetDay(1);
						timeRounded.SetMonth(10);
					} else {
						timeRounded.SetSecond(0);
						timeRounded.SetDay(1);
						timeRounded.SetMonth(10);
						timeRounded.SetYear(timeRounded.GetYear()-1);
					}
				}
				if (vecTimesOut.size() == 0) {
					vecTimesOut.push_back(timeRounded);
					vecTimesOutNew.push_back(timeRounded);
				} else if (timeRounded != vecTimesOut[vecTimesOut.size()-1]) {
					vecTimesOut.push_back(timeRounded);
					vecTimesOutNew.push_back(timeRounded);
				}
				iAccumIndex[t] = vecTimesOut.size()-1;
			}

			// End of output file indicator
			if ((f == vecInputFiles.size()-1) || (vecOutputFiles.size() > 1)) {
				iAccumIndex.push_back(-1);
			}

			// Last processed time from previous file, when reducing multiple input files to one output file
			if ((vecInputFiles.size() > 1) && (vecOutputFiles.size() == 1)) {
				timeLastFileTime = vecTimesIn[vecTimesIn.size()-1];
			}

			// Open output file
			if ((f == 0) || (vecOutputFiles.size() > 1)) {

				_ASSERT(f < vecOutputFiles.size());
				NcFile::FileMode ncfileoutmode = NcFile::Replace;
				if (v != 0) {
					ncfileoutmode = NcFile::Write;
				}

				pncfileout = new NcFile(vecOutputFiles[f].c_str(), ncfileoutmode, NULL, 0, NcFile::Netcdf4);
				if (!pncfileout->is_valid()) {
					_EXCEPTION1("Unable to open output file \"%s\"", vecOutputFiles[f].c_str());
				}
			}
			_ASSERT(pncfileout != NULL);

			// Write new output times to output file
			if (v == 0) {
				if (vecOutputFiles.size() == 1) {
					WriteCFTimeDataToNcFile(pncfileout, vecOutputFiles[0], vecTimesOutNew, true, true);
				} else {
					WriteCFTimeDataToNcFile(pncfileout, vecOutputFiles[f], vecTimesOutNew);
				}
			}

			// Create output variables
			if ((f == 0) || (vecOutputFiles.size() > 1)) {

				// Copy dimensions to output file
				std::vector<NcDim *> vecVarDims;
				vecVarDims.resize(var->num_dims());
				vecVarDims[0] = NcGetTimeDimension(*pncfileout);
				if (vecVarDims[0] == NULL) {
					_EXCEPTIONT("Error writing dimension \"time\" to output file");
				}
				for (long d = 1; d < var->num_dims(); d++) {
					if (var->get_dim(d)->size() != vecSize[d]) {
						_EXCEPTION4("Variable \"%s\" in file \"%s\" has mismatched dimension %li size: expected %li",
							var->name(), vecInputFiles[f].c_str(), d, vecSize[d]);
					}
					NcDim * dimOut = pncfileout->get_dim(var->get_dim(d)->name());
					if (dimOut == NULL) {
						vecVarDims[d] = pncfileout->add_dim(var->get_dim(d)->name(), var->get_dim(d)->size());
						if (vecVarDims[d] == NULL) {
							_EXCEPTION1("Error writing dimension \"%s\" to output file", var->get_dim(d)->name());
						}
						CopyNcVarIfExists(ncfilein, *pncfileout, var->get_dim(d)->name());

					} else {
						vecVarDims[d] = dimOut;
					}
				}
	
				// Create output variable
				varout =
					pncfileout->add_var(
						vecVariableOutNames[v].c_str(),
						ncFloat,
						vecVarDims.size(),
						const_cast<const NcDim**>(&(vecVarDims[0])));

				if (varout == NULL) {
					_EXCEPTION1("Error adding variable \"%s\" to output file", var->name());
				}
	
				if (m_eAccumOp == AccumulateDataOp_Min) {
					varout->add_att("_FillValue", std::numeric_limits<float>::max());
				} else if (m_eAccumOp == AccumulateDataOp_Max) {
					varout->add_att("_FillValue", -std::numeric_limits<float>::max());
				}

				CopyNcVarAttributes(var, varout);

				// Copy preserve variables
				if ((v == 0) && (strPreserveVar != "")) {
					AnnounceStartBlock("Preserving variables");
	
					for (int v = 0; v < vecPreserveVariableStrings.size(); v++) {
						Announce("Variable %s", vecPreserveVariableStrings[v].c_str());
						CopyNcVar(
							ncfilein,
							*pncfileout,
							vecPreserveVariableStrings[v]);
					}
	
					AnnounceEndBlock("Done");
				}

				// Reset accumulation index
				iCurrentAccumIndex = iAccumIndex[0];
			}
	
			// Loop through all input times
			for (size_t t = 0; t < iAccumIndex.size(); t++) {
	
				// Write to file
				if (iCurrentAccumIndex != iAccumIndex[t]) {
					Announce("%s (writing)", vecTimesOut[iCurrentAccumIndex].ToString().c_str());
	
					if (m_eAccumOp == AccumulateDataOp_Avg) {
						if (!fMissingData) {
							for (size_t i = 0; i < sTotalSize; i++) {
								dataAccum[i] /= static_cast<float>(nAccumulatedTimesteps);
							}
						} else {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (dataAccumCount[i] == 0) {
									dataAccum[i] = dFillValue;
								} else {
									dataAccum[i] /= static_cast<float>(dataAccumCount[i]);
								}
							}
						}

					} else {
						if (fMissingData) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (dataAccumCount[i] == 0) {
									dataAccum[i] = dFillValue;
								}
							}
						}
					}
	
					vecPos[0] = iCurrentAccumIndex;
					varout->set_cur(&(vecPos[0]));
					varout->put(&(dataAccum[0]), &(vecSize[0]));
	
					// Reset accumulation arrays
					if ((m_eAccumOp == AccumulateDataOp_Sum) ||
					    (m_eAccumOp == AccumulateDataOp_Avg)
					) {
						dataAccum.Zero();
	
					} else if (m_eAccumOp == AccumulateDataOp_Min) {
						for (size_t i = 0; i < sTotalSize; i++) {
							dataAccum[i] = std::numeric_limits<float>::max();
						}
	
					} else if (m_eAccumOp == AccumulateDataOp_Max) {
						for (size_t i = 0; i < sTotalSize; i++) {
							dataAccum[i] = -std::numeric_limits<float>::max();
						}
					}

					if (fMissingData) {
						dataAccumCount.Zero();
					}

					iCurrentAccumIndex = iAccumIndex[t];
					nAccumulatedTimesteps = 0;
				}

				// Accumulate data
				if (iCurrentAccumIndex != (-1)) {
					Announce("%s => %s (accumulating)",
						vecTimesIn[t].ToString().c_str(),
						vecTimesOut[iCurrentAccumIndex].ToString().c_str());
	
					nAccumulatedTimesteps++;
	
					vecPos[0] = t;
					var->set_cur(&(vecPos[0]));
					var->get(&(dataIn[0]), &(vecSize[0]));

					if (dScaleFactor != 1.0f) {
						for (size_t i = 0; i < sTotalSize; i++) {
							dataIn[i] *= dScaleFactor;
						}
					}
					if (dAddOffset != 0.0f) {
						for (size_t i = 0; i < sTotalSize; i++) {
							dataIn[i] += dAddOffset;
						}
					}
	
					if (!fMissingData) {
						if ((m_eAccumOp == AccumulateDataOp_Sum) ||
						    (m_eAccumOp == AccumulateDataOp_Avg)
						) {
							for (size_t i = 0; i < sTotalSize; i++) {
								dataAccum[i] += dataIn[i];
							}
	
						} else if (m_eAccumOp == AccumulateDataOp_Min) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (dataIn[i] < dataAccum[i]) {
									dataAccum[i] = dataIn[i];
								}
							}
	
						} else if (m_eAccumOp == AccumulateDataOp_Max) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (dataIn[i] > dataAccum[i]) {
									dataAccum[i] = dataIn[i];
								}
							}
						}
	
					} else {
						if (m_eAccumOp == AccumulateDataOp_Sum) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (std::isnan(dataIn[i]) || (dataIn[i] == dFillValue)) {
									continue;
								}
								dataAccum[i] += dataIn[i];
								dataAccumCount[i]++;
							}
	
						} else if (m_eAccumOp == AccumulateDataOp_Avg) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (std::isnan(dataIn[i]) || (dataIn[i] == dFillValue)) {
									continue;
								}
								dataAccum[i] += dataIn[i];
								dataAccumCount[i]++;
							}
	
						} else if (m_eAccumOp == AccumulateDataOp_Min) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (std::isnan(dataIn[i]) || (dataIn[i] == dFillValue)) {
									continue;
								}
								if (dataIn[i] < dataAccum[i]) {
									dataAccum[i] = dataIn[i];
								}
								dataAccumCount[i]++;
							}
	
						} else if (m_eAccumOp == AccumulateDataOp_Max) {
							for (size_t i = 0; i < sTotalSize; i++) {
								if (std::isnan(dataIn[i]) || (dataIn[i] == dFillValue)) {
									continue;
								}
								if (dataIn[i] > dataAccum[i]) {
									dataAccum[i] = dataIn[i];
								}
								dataAccumCount[i]++;
							}
						}
					}
				}
			}
	
			// Close the file
			if (vecOutputFiles.size() > 1) {
				delete pncfileout;
			}

			AnnounceEndBlock("Done");
		}

		// Close the file
		if ((vecInputFiles.size() > 1) && (vecOutputFiles.size() == 1)) {
			delete pncfileout;
		}

		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

