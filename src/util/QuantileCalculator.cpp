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

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////
/*
static unsigned int mylog2 (unsigned int val) {
	if (val == 0) return UINT_MAX;
	if (val == 1) return 0;
	unsigned int ret = 0;
	while (val > 1) {
		val >>= 1;
		ret++;
	}
	return ret;
}
*/
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

	// Dimension to use for quantile calculation (currently can only be "time")
	std::string strDimensionName = "time";

	// Variable to use for quantile calculation
	std::string strVariableName;

	// Output variable
	std::string strVariableOutName;

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
		//CommandLineString(strDimensionName, "dimname", "time");
		CommandLineString(strVariableName, "var", "");
		CommandLineString(strVariableOutName, "varout", "");
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
/*
	unsigned int nLog2QuantileBins = mylog2((unsigned int) nQuantileBins);
	if ((1 << nLog2QuantileBins) != nQuantileBins) {
		_EXCEPTIONT("--bins must be a power of 2");
	}
*/
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

	// Output grid information
	NcFile ncfileout(strOutputData.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
	std::vector<NcDim *> vecOutputGridDim;

	// Allocate space for counts
	DataArray1D<int> nBinCountPreHist(grid.GetSize());
	DataArray2D<int> nBinCounts(grid.GetSize(), nQuantileBins);
	DataArray2D<float> dBinEdges(grid.GetSize(), nQuantileBins+1);
	//DataArray2D<float> dBinEdges(grid.GetSize(), 2);

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

		Variable & varQuantile = varreg.Get(varix);

		// Get minimum and maximum at each point
		//printf("[");
		for (long t = 0; t < dimTime->size(); t++) {
			//Announce("Time %li/%li", t+1, dimTime->size());
			ncfilevec.SetConstantTimeIx(t);
			varQuantile.LoadGridData(varreg, ncfilevec, grid);

			const DataArray1D<float> & data = varQuantile.GetData();

			//printf("%1.14f", data[0]);
			//if (t != dimTime->size()-1) {
			//	printf(",");
			//}

			if (t == 0) {
				for (size_t i = 0; i < grid.GetSize(); i++) {
					dBinEdges(i,0) = data[i];
					dBinEdges(i,nQuantileBins) = data[i];
				}
			} else {
				for (size_t i = 0; i < grid.GetSize(); i++) {
					if (dBinEdges(i,0) > data[i]) {
						dBinEdges(i,0) = data[i];
					}
					if (dBinEdges(i,nQuantileBins) < data[i]) {
						dBinEdges(i,nQuantileBins) = data[i];
					}
				}
			}
		}
		//printf("]\n");

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
			Announce("File \"%s\"", vecInputFiles[f].c_str());

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

			nTotalPointCount += dimTime->size();

			for (long t = 0; t < dimTime->size(); t++) {

				//Announce("File %lu/%lu Time %li/%li", f+1, vecInputFiles.size(), t+1, dimTime->size());
				ncfilevec.SetConstantTimeIx(t);
				varQuantile.LoadGridData(varreg, ncfilevec, grid);

				const DataArray1D<float> & data = varQuantile.GetData();

				// On first iteration expand bins if needed
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
/*
		{
			printf("\n");
			for (int b = 0; b < nQuantileBins; b++) {
				printf("(%i %i) ", b, nBinCounts(27,b));
			}
			printf("\n");
		}
*/
		// Determine which bin contains the correct quantile
		Announce("Identifying correct bin for associated quantile");

		int nQuantileIndex = static_cast<int>(static_cast<double>(nTotalPointCount-1) * dQuantile + 0.5);
		_ASSERT((nQuantileIndex >= 0) && (nQuantileIndex < nTotalPointCount));

		double dAverageBinWidth = 0.0;
		double dMaximumBinWidth = 0.0;

		for (size_t i = 0; i < grid.GetSize(); i++) {
			int nAccumulatedCount = nBinCountPreHist[i];

			int b = 0;
			for (; b < nQuantileBins; b++) {
				nAccumulatedCount += nBinCounts(i,b);
/*
				if (i == 0) {
					printf("%i %i / %i\n", b, nAccumulatedCount, nQuantileIndex);
				}
*/
				if (nAccumulatedCount > nQuantileIndex) {
/*
					if (i == 0) {
						printf("Prev [%1.15f %1.15f] : New [%1.15f %1.15f] : Count %i\n",
							dBinEdges(i,0), dBinEdges(i,nQuantileBins),
							dBinEdges(i,b), dBinEdges(i,b+1),
							nBinCounts(i,b));
					}
*/
					dBinEdges(i,0) = dBinEdges(i,b);
					dBinEdges(i,nQuantileBins) = dBinEdges(i,b+1);
					break;
				}
			}
			if (b == nQuantileBins) {
				_EXCEPTION1("Index %i", i);
			}

			double dBinWidth = static_cast<double>(dBinEdges(i,nQuantileBins)) - static_cast<double>(dBinEdges(i,0));
			dAverageBinWidth += dBinWidth;
			if (dBinWidth > dMaximumBinWidth) {
				dMaximumBinWidth = dBinWidth;
			}

			for (b = 1; b < nQuantileBins; b++) {
				double dAlpha = static_cast<float>(b) / static_cast<float>(nQuantileBins);
				dBinEdges(i,b) = (1.0 - dAlpha) * dBinEdges(i,0) + dAlpha * dBinEdges(i,nQuantileBins);
			}
		}

		Announce("Average bin width: %1.7e", dAverageBinWidth / static_cast<double>(grid.GetSize()));
		Announce("Maximum bin width: %1.7e", dMaximumBinWidth);
		AnnounceEndBlock("Done");
	}

	// Cleanup prior to transposition
	nBinCountPreHist.Deallocate();
	nBinCounts.Deallocate();

	// Write quantiles
	{
		AnnounceStartBlock("Write quantiles");

		// Open output file
		NcVar * varQuantile =
			ncfileout.add_var(
				strVariableOutName.c_str(),
				ncFloat,
				vecOutputGridDim.size(),
				const_cast<const NcDim**>(&(vecOutputGridDim[0])));

		// Transpose quantiles
		DataArray1D<float> dQuantile(grid.GetSize());
		for (int i = 0; i < grid.GetSize(); i++) {
			dQuantile[i] = 0.5 * (dBinEdges(i,0) + dBinEdges(i,nQuantileBins));
		}

		// Write to file
		if (vecOutputGridDim.size() == 1) {
			varQuantile->set_cur((long)0);
			varQuantile->put(&(dQuantile[0]), vecOutputGridDim[0]->size());

		} else if (vecOutputGridDim.size() == 2) {
			varQuantile->set_cur(0, 0);
			varQuantile->put(&(dQuantile[0]), vecOutputGridDim[0]->size(), vecOutputGridDim[1]->size());

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

