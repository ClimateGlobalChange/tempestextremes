///////////////////////////////////////////////////////////////////////////////
///
///	\file    QuantileCalculator.cpp
///	\author  Paul Ullrich
///	\version June 19, 2021
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
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
#include "DataArray2D.h"
#include "Variable.h"
#include "FilenameList.h"
#include "TimeMatch.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>
#include <limits>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {

#if defined(TEMPEST_MPIOMP)
	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
	if (nMPISize > 1) {
		_EXCEPTIONT("At present QuantileCalculator only supports serial execution.");
	}
#endif

	// Input data
	std::string strInputData;

	// Input data list
	std::string strInputDataList;

	// Output data
	std::string strOutputData;

	// Connectivity
	std::string strConnectivity;

	// Data is regional
	bool fRegional;

	// Input includes missing data
	bool fMissingData;

	// Start time
	std::string strStartTime;
	
	// End time
	std::string strEndTime;

	// Time filter
	std::string strTimeFilter;

	// Dimension to use for quantile calculation (currently can only be "time")
	std::string strDimensionName = "time";

	// Variable to use for quantile calculation
	std::string strVariableName;

	// Output variable
	std::string strVariableOutName;

	// Minimum value of the field incorporated in quantile calculation
	std::string strMinThreshold;

	// Maximum value of the field incorporated in quantile calculation
	std::string strMaxThreshold;

	// Quantile to calculate
	double dQuantile;

	// Number of quantile bins
	int nQuantileBins;

	// Number of iterations
	int nIterations;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");
		CommandLineBool(fMissingData, "missingdata");
		CommandLineString(strStartTime, "time_start", "");
		CommandLineString(strEndTime, "time_end", "");
		CommandLineString(strTimeFilter, "timefilter", "");
		//CommandLineString(strDimensionName, "dimname", "time");
		CommandLineString(strVariableName, "var", "");
		CommandLineString(strVariableOutName, "varout", "");
		CommandLineString(strMinThreshold, "minthreshold", "");
		CommandLineString(strMaxThreshold, "maxthreshold", "");
		CommandLineDouble(dQuantile, "quantile", 0.95);
		CommandLineInt(nQuantileBins, "bins", 16);
		CommandLineInt(nIterations, "iter", 8);

		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input file (--in_data) or (--in_data_list) specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list) may be specified");
	}
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}
	if (nQuantileBins < 2) {
		_EXCEPTIONT("--bins must be at least 2");
	}
	if ((dQuantile < 0.0) || (dQuantile > 1.0)) {
		_EXCEPTIONT("--quantile must be between 0 and 1, inclusive");
	}

	if (nIterations < 1) {
		_EXCEPTIONT("--iter must be at least 1");
	}

	if (strVariableOutName == "") {
		strVariableOutName = strVariableName;
	}

	// Thresholds (automatically turn on --missingdata)
	double dMinThreshold = -std::numeric_limits<float>::max();
	double dMaxThreshold = std::numeric_limits<float>::max();
	if (strMinThreshold != "") {
		if (!STLStringHelper::IsFloat(strMinThreshold)) {
			_EXCEPTIONT("--minthreshold must be a floating point number");
		}
		dMinThreshold = std::stof(strMinThreshold);
		fMissingData = true;
	}
	if (strMaxThreshold != "") {
		if (!STLStringHelper::IsFloat(strMaxThreshold)) {
			_EXCEPTIONT("--maxthreshold must be a floating point number");
		}
		dMaxThreshold = std::stof(strMaxThreshold);
		fMissingData = true;
	}
	if (dMaxThreshold < dMinThreshold) {
		_EXCEPTIONT("--maxthreshold must be greater than or equal to --minthreshold");
	}

	// Start and end time (initialized on first pass through files)
	Time timeStartTime;
	Time timeEndTime;

	// Time filter
	bool fTimeFilterSpecified = (strTimeFilter != "");
#ifdef TEMPEST_NOREGEX
	if (strTimeFilter != "") {
		_EXCEPTIONT("Cannot use --timefilter with -DTEMPEST_NOREGEX compiler flag");
	}
#endif
#ifndef TEMPEST_NOREGEX
	std::regex reTimeSubset;
	if (strTimeFilter != "") {
		
		// Test regex support
		TestRegex();

		if (strTimeFilter == "3hr") {
			strTimeFilter = "....-..-.. (00|03|06|09|12|15|18|21):00:00";
		}
		if (strTimeFilter == "6hr") {
			strTimeFilter = "....-..-.. (00|06|12|18):00:00";
		}
		if (strTimeFilter == "daily") {
			strTimeFilter = "....-..-.. 00:00:00";
		}

		try {
			reTimeSubset.assign(strTimeFilter);
		} catch(std::regex_error & reerr) {
			_EXCEPTION2("Parse error in --timefilter regular expression \"%s\" (code %i)",
				strTimeFilter.c_str(), reerr.code());
		}
	}
#endif

	// Create Variable registry
	VariableRegistry varreg;
	VariableIndex varix = varreg.FindOrRegister(strVariableName);
	_ASSERT(varix != InvalidVariableIndex);

	// Input file list
	FilenameList vecInputFiles;

	if (strInputData != "") {
		vecInputFiles.push_back(strInputData);
	}
	if (strInputDataList != "") {
		vecInputFiles.FromFile(strInputDataList);
	}
	if (vecInputFiles.size() == 0) {
		_EXCEPTIONT("No input data files found");
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		// Load in file vector
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[0]);
		_ASSERT(vecNcFiles.size() > 0);

		grid.GenerateLatitudeLongitude(
			vecNcFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			true);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	size_t sGridSize = grid.GetSize();

	// Output grid information
	NcFile ncfileout(strOutputData.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
	std::vector<NcDim *> vecOutputGridDim;
	NcVar * varOutputQuantile = NULL;

	// Fill value
	float dFillValue = std::numeric_limits<float>::max();

	// Allocate space for counts
	DataArray1D<int> nBinCountPreHist(grid.GetSize());
	DataArray2D<int> nBinCounts(grid.GetSize(), nQuantileBins);
	DataArray2D<float> dBinEdges(grid.GetSize(), nQuantileBins+1);

	DataArray1D<int> nLocalPointCounts;
	if (fMissingData) {
		nLocalPointCounts.Allocate(grid.GetSize());
	}

	// Initialize first pass bin sizes from first file
	{
		AnnounceStartBlock("Using first file to set array bounds");
		NcFileVector ncfilevec;
		ncfilevec.ParseFromString(vecInputFiles[0]);
		if (ncfilevec.size() == 0) {
			_EXCEPTION1("No NetCDF files found in \"%s\"", vecInputFiles[0].c_str());
		}

		NcDim * dimTime = ncfilevec[0]->get_dim(strDimensionName.c_str());
		if (dimTime == NULL) {
			_EXCEPTION2("Unable to load dimension \"%s\" from \"%s\"",
				strDimensionName.c_str(),
				ncfilevec.GetFilename(0).c_str());
		}

		// Get grid variables from first input file
		{
			if (grid.m_nGridDim.size() == 1) {
				NcDim * dim0 = ncfileout.add_dim("ncol", grid.m_nGridDim[0]);
				vecOutputGridDim.push_back(dim0);

			} else if (grid.m_nGridDim.size() == 2) {
				NcDim * dim0 = ncfileout.add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]);
				NcDim * dim1 = ncfileout.add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]);

				vecOutputGridDim.push_back(dim0);
				vecOutputGridDim.push_back(dim1);

				CopyNcVarIfExists(
					*(ncfilevec[0]),
					ncfileout,
					strLatitudeName.c_str(),
					true,
					true);

				CopyNcVarIfExists(
					*(ncfilevec[0]),
					ncfileout,
					strLongitudeName.c_str(),
					true,
					true);

			} else {
				_EXCEPTIONT("Only 1D or 2D spatial data supported");
			}
		}

		// Variable for computing quantile
		Variable & varQuantile = varreg.Get(varix);

		// Create output variable
		{
			// Open output file
			varOutputQuantile =
				ncfileout.add_var(
					strVariableOutName.c_str(),
					ncFloat,
					vecOutputGridDim.size(),
					const_cast<const NcDim**>(&(vecOutputGridDim[0])));

			if (varOutputQuantile == NULL) {
				_EXCEPTION2("Unable to create variable \"%s\" in output file \"%s\"",
					strVariableOutName.c_str(), strOutputData.c_str());
			}

			// Copy attributes from input file for output variable
			if (!varQuantile.IsOp()) {
				NcVar * ncvarIn = ncfilevec[0]->get_var(strVariableName.c_str());
				if (ncvarIn != NULL) {
					CopyNcVarAttributes(ncvarIn, varOutputQuantile);
				}
			}

			// Add attributes describing quantile
			varOutputQuantile->add_att("quantile", dQuantile);
			varOutputQuantile->add_att("quantile_info", "via TempestExtremes QuantileCalculator");
		}

		// Get minimum and maximum at each point
		int nRetainedTimeSlices = 0;

		for (long t = 0; t < dimTime->size(); t++) {
			//Announce("Time %li/%li", t+1, dimTime->size());
			ncfilevec.SetConstantTimeIx(t);
			varQuantile.LoadGridData(varreg, ncfilevec, grid);

#ifndef TEMPEST_NOREGEX
			// Time filter
			if (fTimeFilterSpecified) {
				const NcTimeDimension & vecTimes = ncfilevec.GetNcTimeDimension(0);
				_ASSERT(vecTimes.size() == dimTime->size());
				std::string strTime = vecTimes[t].ToString();
				std::smatch match;
				
				if (!std::regex_search(strTime, match, reTimeSubset)) {
					continue;
				}
			}
#endif
			// Apply time bounds
			if ((strStartTime != "") || (strEndTime != "")) {
				const NcTimeDimension & vecTimes = ncfilevec.GetNcTimeDimension(0);
				_ASSERT(vecTimes.size() == dimTime->size());
				if (strStartTime != "") {
					if (t == 0) {
						timeStartTime = Time(vecTimes[0].GetCalendarType());
						timeStartTime.FromFormattedString(strStartTime);
					}
					double dDeltaSeconds = timeStartTime - vecTimes[t];
					if (dDeltaSeconds > 0.0) {
						continue;
					}
				}
				if (strEndTime != "") {
					if (t == 0) {
						timeEndTime = Time(vecTimes[0].GetCalendarType());
						timeEndTime.FromFormattedString(strEndTime);
					}
					double dDeltaSeconds = vecTimes[t] - timeEndTime;
					if (dDeltaSeconds > 0.0) {
						continue;
					}
				}
			}

			// This time slice is retained
			nRetainedTimeSlices++;

			// Get _FillValue if it exists
			if (dFillValue == std::numeric_limits<float>::max()) {
				dFillValue = varQuantile.GetFillValueFloat();
			}

			// Load data
			const DataArray1D<float> & data = varQuantile.GetData();

			// Estimate first bin edges with missing data
			if (fMissingData) {
				if (t == 0) {
					for (size_t i = 0; i < sGridSize; i++) {
						if ((data[i] == dFillValue) || (data[i] != data[i])) {
							dBinEdges(i,0) = dFillValue;
						} else if ((data[i] < dMinThreshold) || (data[i] > dMaxThreshold)) {
							dBinEdges(i,0) = dFillValue;
						} else {
							dBinEdges(i,0) = data[i];
						}
						dBinEdges(i,nQuantileBins) = dFillValue;
					}

				} else {
					for (size_t i = 0; i < sGridSize; i++) {
						if ((data[i] == dFillValue) || (data[i] != data[i])) {
							continue;
						}
						if ((data[i] < dMinThreshold) || (data[i] > dMaxThreshold)) {
							continue;
						}
						if ((dBinEdges(i,0) == dFillValue) || (dBinEdges(i,0) != dBinEdges(i,0))) {
							dBinEdges(i,0) = data[i];
							continue;
						}
						if ((dBinEdges(i,nQuantileBins) == dFillValue) ||
						    (dBinEdges(i,nQuantileBins) != dBinEdges(i,nQuantileBins))
						) {
							if (data[i] < dBinEdges(i,0)) {
								dBinEdges(i,nQuantileBins) = dBinEdges(i,0);
								dBinEdges(i,0) = data[i];
							} else {
								dBinEdges(i,nQuantileBins) = data[i];
							}
							continue;
						}
						if (dBinEdges(i,0) > data[i]) {
							dBinEdges(i,0) = data[i];
						}
						if (dBinEdges(i,nQuantileBins) < data[i]) {
							dBinEdges(i,nQuantileBins) = data[i];
						}
					}
				}

			// Estimate first bin edges without missing data
			} else {
				if (t == 0) {
					for (size_t i = 0; i < sGridSize; i++) {
						dBinEdges(i,0) = data[i];
						dBinEdges(i,nQuantileBins) = data[i];
					}

				} else {
					for (size_t i = 0; i < sGridSize; i++) {
						if (dBinEdges(i,0) > data[i]) {
							dBinEdges(i,0) = data[i];
						}
						if (dBinEdges(i,nQuantileBins) < data[i]) {
							dBinEdges(i,nQuantileBins) = data[i];
						}
					}
				}
			}
		}

		if ((!fMissingData) && (nRetainedTimeSlices < 2)) {
			_EXCEPTIONT("At present, at least 2 time slices must be available in the first file provided.  Reorganize file list or rerun with --missingdata.");
		}

		// Use a linear fit for bin edges
		for (int i = 0; i < grid.GetSize(); i++) {
			for (int b = 1; b < nQuantileBins; b++) {
				double dAlpha = static_cast<float>(b) / static_cast<float>(nQuantileBins);
				dBinEdges(i,b) = (1.0 - dAlpha) * dBinEdges(i,0) + dAlpha * dBinEdges(i,nQuantileBins);
			}
		}

		AnnounceEndBlock("Done");
	}

	// Run the quantile calculation
	for (int it = 0; it < nIterations; it++) {
		AnnounceStartBlock("Running quantiles: Iteration %i / %i", it+1, nIterations);

		// Get the Variable
		Variable & varQuantile = varreg.Get(varix);

		// Total count of points
		int nTotalPointCount = 0;

		// Populate bins from all files
		nBinCountPreHist.Zero();
		nBinCounts.Zero();

		for (size_t f = 0; f < vecInputFiles.size(); f++) {
			AnnounceStartBlock("File \"%s\"", vecInputFiles[f].c_str());

			NcFileVector ncfilevec;
			ncfilevec.ParseFromString(vecInputFiles[f]);
			if (ncfilevec.size() == 0) {
				_EXCEPTION1("No NetCDF files found in \"%s\"", vecInputFiles[f].c_str());
			}

			NcDim * dimTime = ncfilevec[0]->get_dim(strDimensionName.c_str());
			if (dimTime == NULL) {
				_EXCEPTION2("Unable to load dimension \"%s\" from \"%s\"",
					strDimensionName.c_str(),
					ncfilevec.GetFilename(0).c_str());
			}

			int nRetainedTimeSlices = 0;

			for (long t = 0; t < dimTime->size(); t++) {

#ifndef TEMPEST_NOREGEX
				// Time filter
				if (fTimeFilterSpecified) {
					const NcTimeDimension & vecTimes = ncfilevec.GetNcTimeDimension(0);
					_ASSERT(vecTimes.size() == dimTime->size());
					std::string strTime = vecTimes[t].ToString();
					std::smatch match;
				
					if (!std::regex_search(strTime, match, reTimeSubset)) {
						continue;
					}
				}
#endif
				// Apply time bounds
				if ((strStartTime != "") || (strEndTime != "")) {
					const NcTimeDimension & vecTimes = ncfilevec.GetNcTimeDimension(0);
					_ASSERT(vecTimes.size() == dimTime->size());
					if (strStartTime != "") {
						double dDeltaSeconds = timeStartTime - vecTimes[t];
						if (dDeltaSeconds > 0.0) {
							continue;
						}
					}
					if (strEndTime != "") {
						double dDeltaSeconds = vecTimes[t] - timeEndTime;
						if (dDeltaSeconds > 0.0) {
							continue;
						}
					}
				}

				// This time slice is retained
				nRetainedTimeSlices++;

				// Load data
				ncfilevec.SetConstantTimeIx(t);
				varQuantile.LoadGridData(varreg, ncfilevec, grid);

				const DataArray1D<float> & data = varQuantile.GetData();

				// Get _FillValue if it exists
				float dLocalFillValue = varQuantile.GetFillValueFloat();
				if (dLocalFillValue != dFillValue) {
					_EXCEPTIONT("_FillValue attribute appears to change across files");
				}

				// Build histograms with missing data
				if (fMissingData) {

					// On first iteration expand bins if needed
					if (it == 0) {
						for (size_t i = 0; i < sGridSize; i++) {
							if ((data[i] == dFillValue) || (data[i] != data[i])) {
								continue;
							}
							if ((data[i] < dMinThreshold) || (data[i] > dMaxThreshold)) {
								continue;
							}

							nLocalPointCounts[i]++;

							if ((dBinEdges(i,0) == dFillValue) || (dBinEdges(i,0) != dBinEdges(i,0))) {
								dBinEdges(i,0) = data[i];
								nBinCounts(i,0)++;
								continue;
							}
							if ((dBinEdges(i,nQuantileBins) == dFillValue) ||
							    (dBinEdges(i,nQuantileBins) != dBinEdges(i,nQuantileBins))
							) {
								if (data[i] == dBinEdges(i,0)) {
									nBinCounts(i,0)++;
									continue;

								} else if (data[i] < dBinEdges(i,0)) {
									dBinEdges(i,nQuantileBins) = dBinEdges(i,0);
									dBinEdges(i,0) = data[i];
									nBinCounts(i,nQuantileBins-1) = nBinCounts(i,0);
									nBinCounts(i,0) = 0;

								} else {
									dBinEdges(i,nQuantileBins) = data[i];
								}

								for (int b = 1; b < nQuantileBins; b++) {
									double dAlpha = static_cast<float>(b) / static_cast<float>(nQuantileBins);
									dBinEdges(i,b) = (1.0 - dAlpha) * dBinEdges(i,0) + dAlpha * dBinEdges(i,nQuantileBins);
								}
							}

							if (dBinEdges(i,0) >= data[i]) {
								dBinEdges(i,0) = data[i];
								nBinCounts(i,0)++;
								continue;
							}

							int b = 1;
							for (; b < nQuantileBins; b++) {
								if (dBinEdges(i,b) >= data[i]) {
									nBinCounts(i,b-1)++;
									break;
								}
							}

							if (b == nQuantileBins) {
								if (data[i] > dBinEdges(i,b)) {
									dBinEdges(i,b) = data[i];
								}
								nBinCounts(i,b-1)++;
							}
						}

					// ..otherwise discard values outside of histogram range
					} else {
						for (size_t i = 0; i < sGridSize; i++) {
							if ((data[i] == dFillValue) || (data[i] != data[i])) {
								continue;
							}
							if ((data[i] < dMinThreshold) || (data[i] > dMaxThreshold)) {
								continue;
							}
							if ((dBinEdges(i,nQuantileBins) == dFillValue) ||
							    (dBinEdges(i,nQuantileBins) != dBinEdges(i,nQuantileBins))
							) {
								continue;
							}

							if (dBinEdges(i,0) > data[i]) {
								nBinCountPreHist[i]++;
								continue;
							}

							for (int b = 1; b <= nQuantileBins; b++) {
								if (dBinEdges(i,b) >= data[i]) {
									nBinCounts(i,b-1)++;
									break;
								}
							}
						}
					}

				// Build histograms without missing data
				} else {
					if (it == 0) {
						for (size_t i = 0; i < grid.GetSize(); i++) {
							if (dBinEdges(i,0) >= data[i]) {
								dBinEdges(i,0) = data[i];
								nBinCounts(i,0)++;
								continue;
							}

							int b = 1;
							for (; b < nQuantileBins; b++) {
								if (dBinEdges(i,b) >= data[i]) {
									nBinCounts(i,b-1)++;
									break;
								}
							}

							if (b == nQuantileBins) {
								if (data[i] > dBinEdges(i,b)) {
									dBinEdges(i,b) = data[i];
								}
								nBinCounts(i,b-1)++;
							}
						}

					// ..otherwise discard values outside of histogram range
					} else {
						for (size_t i = 0; i < grid.GetSize(); i++) {
							if (dBinEdges(i,0) > data[i]) {
								nBinCountPreHist[i]++;
								continue;
							}

							for (int b = 1; b <= nQuantileBins; b++) {
								if (dBinEdges(i,b) >= data[i]) {
									nBinCounts(i,b-1)++;
									break;
								}
							}
						}
					}
				}
			}

			nTotalPointCount += nRetainedTimeSlices;

			if (fTimeFilterSpecified || (strStartTime != "") || (strEndTime != "")) {
				Announce("%i/%li time slices retained after time filtering", nRetainedTimeSlices, dimTime->size());
			}

			AnnounceEndBlock(NULL);
		}

		// Determine which bin contains the correct quantile
		Announce("Identifying correct bin for associated quantile");

		double dAverageBinWidth = 0.0;
		double dMaximumBinWidth = 0.0;

		int nQuantileIndex = static_cast<int>(static_cast<double>(nTotalPointCount-1) * dQuantile + 0.5);
		_ASSERT((nQuantileIndex >= 0) && (nQuantileIndex < nTotalPointCount));

		size_t sActiveGridPoints = sGridSize;

		for (size_t i = 0; i < sGridSize; i++) {

			// If there is missing data the number of samples may vary by grid point
			if (fMissingData) {
				if ((dBinEdges(i,nQuantileBins) == dFillValue) ||
				    (dBinEdges(i,nQuantileBins) != dBinEdges(i,nQuantileBins))
				) {
					sActiveGridPoints--;
					continue;
				}

				nQuantileIndex = static_cast<int>(static_cast<double>(nLocalPointCounts[i]-1) * dQuantile + 0.5);
				_ASSERT((nQuantileIndex >= 0) && (nQuantileIndex < nLocalPointCounts[i]));
			}

			// Find the bin with the right sample number
			int nAccumulatedCount = nBinCountPreHist[i];

			int b = 0;
			for (; b < nQuantileBins; b++) {
				nAccumulatedCount += nBinCounts(i,b);

				if (nAccumulatedCount > nQuantileIndex) {
					dBinEdges(i,0) = dBinEdges(i,b);
					dBinEdges(i,nQuantileBins) = dBinEdges(i,b+1);
					break;
				}
			}
			if (b == nQuantileBins) {
				_EXCEPTION1("Index %i", i);
			}

			// Calculate the average bin width
			double dBinWidth = static_cast<double>(dBinEdges(i,nQuantileBins)) - static_cast<double>(dBinEdges(i,0));
			dAverageBinWidth += dBinWidth;
			if (dBinWidth > dMaximumBinWidth) {
				dMaximumBinWidth = dBinWidth;
			}

			// Set new bin edges by subdividing this bin
			for (b = 1; b < nQuantileBins; b++) {
				double dAlpha = static_cast<float>(b) / static_cast<float>(nQuantileBins);
				dBinEdges(i,b) = (1.0 - dAlpha) * dBinEdges(i,0) + dAlpha * dBinEdges(i,nQuantileBins);
			}
		}

		Announce("Average bin width: %1.7e", dAverageBinWidth / static_cast<double>(sActiveGridPoints));
		Announce("Maximum bin width: %1.7e", dMaximumBinWidth);
		AnnounceEndBlock("Done");
	}

	// Cleanup prior to transposition
	nBinCountPreHist.Deallocate();
	nBinCounts.Deallocate();

	// Write quantiles
	{
		AnnounceStartBlock("Write quantiles");

		// Get quantiles
		DataArray1D<float> dQuantile(grid.GetSize());

		if (fMissingData) {
			for (int i = 0; i < grid.GetSize(); i++) {
				if ((dBinEdges(i,nQuantileBins) == dFillValue) ||
				    (dBinEdges(i,nQuantileBins) != dBinEdges(i,nQuantileBins))
				) {
					dQuantile[i] = dFillValue;
				} else {
					dQuantile[i] = 0.5 * (dBinEdges(i,0) + dBinEdges(i,nQuantileBins));
				}
			}

		} else {
			for (int i = 0; i < grid.GetSize(); i++) {
				dQuantile[i] = 0.5 * (dBinEdges(i,0) + dBinEdges(i,nQuantileBins));
			}
		}

		// Write to file
		if (vecOutputGridDim.size() == 1) {
			varOutputQuantile->set_cur((long)0);
			varOutputQuantile->put(&(dQuantile[0]), vecOutputGridDim[0]->size());

		} else if (vecOutputGridDim.size() == 2) {
			varOutputQuantile->set_cur(0, 0);
			varOutputQuantile->put(&(dQuantile[0]), vecOutputGridDim[0]->size(), vecOutputGridDim[1]->size());

		} else {
			_EXCEPTION();
		}

		AnnounceEndBlock("Done");
	}

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

