///////////////////////////////////////////////////////////////////////////////
///
///	\file    CalculateIDF.cpp
///	\author  Paul Ullrich
///	\version August 6, 2026
///
///	<remarks>
///		Copyright 2026 Paul Ullrich
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
#include "SimpleGrid.h"
#include "STLStringHelper.h"
#include "FilenameList.h"
#include "NetCDFUtilities.h"
#include "NcFileVector.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "NcVarSpatialIterator.h"

#include "netcdfcpp.h"

#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
/*
///////////////////////////////////////////////////////////////////////////////

double GEVNegativeLogLikelihood(
	const std::vector<double> & vecData,
	const GEVParameters & p)
{
	if (p.sigma <= 0.0) {
		return std::numeric_limits<double>::infinity();
	}

	const double dEps = 1.0e-10;
	double dL = 0.0;

	if (std::abs(p.xi) < 1e-6) {
		// Gumbel limit
		for (double dValue : vecData) {
			double dZ = (dValue - p.mu) / p.sigma;
			dL += std::log(p.sigma) + dZ + std::exp(-dZ);
		}

		return dL;
	}

	for (double dValue : vecData) {
		double dT = 1.0 + p.xi * (dValue - p.mu) / p.sigma;

		if (dT <= dEps) {
			return std::numeric_limits<double>::infinity();
		}

		dL += std::log(p.sigma)
		   + (1.0 + 1.0/p.xi) * std::log(dT)
		   + std::pow(dT, -1.0/p.xi);
	}

	return dL;
}
*/
///////////////////////////////////////////////////////////////////////////////

void CalculateIDF(
	DataArray2D<float> & dData,
	DataArray2D<float> & dIDFParameters,
	size_t sPointsRef,
	float dFillValue
) {
	size_t sTimeCount = dData.GetColumns();
	for (size_t i = 0; i < sPointsRef; i++) {

		// Don't operate if no data in this column
		size_t sNonFillValueCount = 0;
		for (size_t j = 0; j < sTimeCount; j++) {
			if (std::isnan(dData(i,j)) || (dData(i,j) == dFillValue)) {
				continue;
			}
			sNonFillValueCount++;
		}
		if (sNonFillValueCount == 0) {
			dIDFParameters(0,i) = -std::numeric_limits<float>::max();
			dIDFParameters(1,i) = -std::numeric_limits<float>::max();
			dIDFParameters(2,i) = -std::numeric_limits<float>::max();
			continue;
		}

		// Sort data
		std::sort(&(dData(i,0)), &(dData(i,0)) + sTimeCount);

		// Compute Probability Weighted Moments
		double dB0 = 0.0;
		double dB1 = 0.0;
		double dB2 = 0.0;

		double dTimeCount = static_cast<double>(sTimeCount);
		for (size_t j = 0; j < sTimeCount; j++) {
			dB0 += dData(i,j);
			if (sTimeCount > 1) {
				dB1 += static_cast<double>(j) / (dTimeCount - 1.0) * dData(i,j);
			}
			if (sTimeCount > 2) {
				dB2 += static_cast<double>(j * (j-1)) / ((dTimeCount - 1.0) * (dTimeCount - 2.0)) * dData(i,j);
			}
		}
		dB0 /= dTimeCount;
		dB1 /= dTimeCount;
		dB2 /= dTimeCount;

		// Convert to L-moments
		double dL1 = dB0;
		double dL2 = 2.0*dB1 - dB0;
		double dL3 = 6.0*dB2 - 6.0*dB1 + dB0;
		double dTau3 = dL3/dL2;

		// Calculate GEV parameters
		double dC = 2.0/(3.0+dTau3) - std::log(2.0)/std::log(3.0);
		double dXi = 7.8590*dC + 2.9554*dC*dC;
		double dGamma = std::tgamma(1.0 + dXi);
		double dSigma = dL2 * dXi / (dGamma * (1.0 - std::pow(2.0,-dXi)));
		double dMu = dL1 - dSigma/dXi * (1.0 - dGamma);

		// Store parameters
		dIDFParameters(0,i) = dMu;
		dIDFParameters(1,i) = dSigma;
		dIDFParameters(2,i) = dXi;
/*
		for (size_t j = 0; j < sTimeCount; j++) {
			printf("%1.5f\n", dData(i,j));
		}
		printf("%1.5f %1.5f %1.5f : mean %1.5f\n", dMu, dSigma, dXi, dMu + (std::tgamma(1.0 - dXi) - 1.0) / dXi);

		_EXCEPTION();
*/
	}
}

///////////////////////////////////////////////////////////////////////////////

double IDFReturnLevel(
	double dMu,
	double dSigma,
	double dXi,
	double dT
) {
    double dPnon = 1.0 - 1.0 / dT;

    if (std::abs(dXi) < 1.0e-6) {
        return dMu - dSigma * std::log(-std::log(dPnon));
    }

    return dMu
         + dSigma / dXi *
           (std::pow(-std::log(dPnon), -dXi) - 1.0);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input data file
	std::string strInputData;

	// Input list of data files
	std::string strInputDataList;

	// Connectivity file
	std::string strConnectivity;

	// Grid file
	std::string strGridFile;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Data is regional
	bool fRegional;

	// Output data file
	std::string strOutputData;

	// Variable to process
	std::string strVariable;

	// Output the parameters of the GEV fit
	bool fOutputParameters;

	// List of frequencies of the GEV fit
	std::string strFrequencies;

	// Units for the frequency variable
	std::string strFrequencyUnit;

	// Number of spatial points to load in each block
	int nBlockSize;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strGridFile, "in_grid", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");
		CommandLineString(strOutputData, "out_data", "");

		CommandLineString(strVariable, "var", "");
		CommandLineBool(fOutputParameters, "output_param");
		CommandLineString(strFrequencies, "freq", "");
		CommandLineString(strFrequencyUnit, "freq_unit", "year");

		CommandLineInt(nBlockSize, "block_size", 0);
		CommandLineString(strLongitudeName, "lonname", "[auto]");
		CommandLineString(strLatitudeName, "latname", "[auto]");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file --out_data specified");
	}
	if (strVariable.length() == 0) {
		_EXCEPTIONT("No variable --var specified");
	}
	if ((strFrequencies == "") && !fOutputParameters) {
		_EXCEPTIONT("One of --freq or --output_param must be set");
	}

	// Store input data
	FilenameList vecInputFileList;

	if (strInputData.length() != 0) {
		vecInputFileList.push_back(strInputData);

	} else {
		AnnounceStartBlock("Building input data list");
		vecInputFileList.FromFile(strInputDataList);
		AnnounceEndBlock("Done");
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

		if (vecInputFileList.size() < 1) {
			_EXCEPTIONT("No data files specified; unable to generate grid");
		}

		// Get the list of data files
		NcFileVector vecFiles;
		if (strGridFile == "") {
			vecFiles.ParseFromString(vecInputFileList[0]);
		} else {
			vecFiles.ParseFromString(strGridFile);
		}
		_ASSERT(vecFiles.size() > 0);

		grid.GenerateLatitudeLongitude(
			vecFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Parse freqs
	std::vector<float> vecFrequencies;
	if (strFrequencies != "") {
		std::vector<std::string> vecFrequenciesStr;
		STLStringHelper::ParseVariableList(strFrequencies, vecFrequenciesStr);
		for (size_t i = 0; i < vecFrequenciesStr.size(); i++) {
			if (!STLStringHelper::IsFloat(vecFrequenciesStr[i])) {
				_EXCEPTIONT("--freq must be a list of floating point numbers");
			}
			vecFrequencies.push_back(std::stof(vecFrequenciesStr[i]));
		}
	}

	/////////////////////////////////////////////
	// Number of blocks
	if (nBlockSize < 0) {
		_EXCEPTIONT("--block_size must be nonnegative");
	}

	size_t sBlockSize = static_cast<size_t>(nBlockSize);

	/////////////////////////////////////////////
	// Check chunk information

	// Open output file
	NcFile ncout(strOutputData.c_str(), NcFile::Replace);
	if (!ncout.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputData.c_str());
	}

	// Chunk size of variable to be processed
	std::vector<NcDim *> vecDimOut;
	std::vector<size_t> vecChunkSize;
	size_t sTimeCount = 0;
	float dFillValue = -std::numeric_limits<float>::max();

	// Variables for storying IDF parameters
	NcVar * varMu = NULL;
	NcVar * varSigma = NULL;
	NcVar * varXi = NULL;

	// Freuqency variables
	NcDim * dimFreq = NULL;
	NcVar * varFreq = NULL;

	{
		AnnounceStartBlock("Initializing output file");

		// Get the list of data files
		NcFileVector vecFiles;
		vecFiles.ParseFromString(vecInputFileList[0]);

		_ASSERT(vecFiles.size() > 0);

		// Get the variable
		NcVar * varIn = NULL;
		vecFiles.FindContainingVariable(strVariable, &varIn);
		if (varIn == NULL) {
			_EXCEPTION2("Cannot find variable \"%s\" in file \"%s\"", 
				strVariable.c_str(), vecFiles.GetFilename(0).c_str());
		}

		if ((varIn->num_dims() != 2) && (varIn->num_dims() != 3)) {
			_EXCEPTION1("Variable \"%s\" must have 2 or 3 dimensions",
				strVariable.c_str());
		}

		const NcTimeDimension & vecTimes =
			vecFiles.GetNcTimeDimension(0);

		sTimeCount += vecTimes.size();

		// Get FillValue
		NcAtt * attFillValue = varIn->get_att("_FillValue");
		if (attFillValue != NULL) {
			dFillValue = attFillValue->as_float(0);
		}

		// Get count of times across all files
		for (int f = 1; f < vecInputFileList.size(); f++) {
			NcFileVector vecFilesSeq;
			vecFilesSeq.ParseFromString(vecInputFileList[f]);

			_ASSERT(vecFilesSeq.size() > 0);

			NcVar * varInSeq = NULL;
			vecFiles.FindContainingVariable(strVariable, &varInSeq);
			if (varInSeq == NULL) {
				_EXCEPTION2("Cannot find variable \"%s\" in file \"%s\"", 
					strVariable.c_str(), vecFiles.GetFilename(0).c_str());
			}
			if (varInSeq->num_dims() != varIn->num_dims()) {
				_EXCEPTION4("Variable \"%s\" in file \"%s\" does not match structure of variable \"%s\" in file \"%s\"",
					strVariable.c_str(), vecFiles.GetFilename(f).c_str(),
					strVariable.c_str(), vecFiles.GetFilename(0).c_str());
			}

			const NcTimeDimension & vecTimesSeq =
				vecFilesSeq.GetNcTimeDimension(0);

			sTimeCount += vecTimesSeq.size();
		}

		// Initialize vector of output dimensions
		if (grid.m_nGridDim.size() == 1) {
			vecDimOut.push_back(ncout.add_dim("ncol", grid.m_nGridDim[0]));
		} else if (grid.m_nGridDim.size() == 2) {
			vecDimOut.push_back(ncout.add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]));
			vecDimOut.push_back(ncout.add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]));
		} else {
			_EXCEPTIONT("Only 1D or 2D spatial data supported");
		}

		// Copy latitude and longitude variables to output file
		CopyNcVarIfExists(
			*(vecFiles[0]),
			ncout,
			strLatitudeName,
			true,
			true);

		CopyNcVarIfExists(
			*(vecFiles[0]),
			ncout,
			strLongitudeName,
			true,
			true);

		// Get chunking information
		if (nBlockSize != 0) {
			vecChunkSize = varIn->get_chunk_size();

			if (vecChunkSize.size() == 0) {
				Announce("Variable %s is not chunked", strVariable.c_str());
			} else if (vecChunkSize.size() == 2) {
				Announce("Variable %s has chunking pattern [%lu,%lu]",
					strVariable.c_str(), vecChunkSize[0], vecChunkSize[1]);
			} else if (vecChunkSize.size() == 3) {
				Announce("Variable %s has chunking pattern [%lu,%lu,%lu]",
					strVariable.c_str(), vecChunkSize[0], vecChunkSize[1], vecChunkSize[2]);
			}
		}

		// Create output variables for IDF parameters
		if (fOutputParameters) {
			Announce("Initializing parameter variables");

			std::string strVariableMu = strVariable + "_mu";
			std::string strVariableSigma = strVariable + "_sigma";
			std::string strVariableXi = strVariable + "_xi";

			varMu =
				ncout.add_var(
					strVariableMu.c_str(),
					ncFloat,
					vecDimOut.size(),
					const_cast<const NcDim**>(&(vecDimOut[0])));

			if (varMu == NULL) {
				_EXCEPTION2("Cannot create variable \"%s\" in output file \"%s\"",
					strVariableMu.c_str(), strOutputData.c_str());
			}

			varSigma =
				ncout.add_var(
					strVariableSigma.c_str(),
					ncFloat,
					vecDimOut.size(),
					const_cast<const NcDim**>(&(vecDimOut[0])));

			if (varSigma == NULL) {
				_EXCEPTION2("Cannot create variable \"%s\" in output file \"%s\"",
					strVariableSigma.c_str(), strOutputData.c_str());
			}

			varXi =
				ncout.add_var(
					strVariableXi.c_str(),
					ncFloat,
					vecDimOut.size(),
					const_cast<const NcDim**>(&(vecDimOut[0])));

			if (varXi == NULL) {
				_EXCEPTION2("Cannot create variable \"%s\" in output file \"%s\"",
					strVariableXi.c_str(), strOutputData.c_str());
			}

			varMu->add_att("_FillValue", -std::numeric_limits<float>::max());
			varSigma->add_att("_FillValue", -std::numeric_limits<float>::max());
			varXi->add_att("_FillValue", -std::numeric_limits<float>::max());

			if (vecChunkSize.size() != 0) {
				std::vector<size_t> vecChunkSizeOut(vecChunkSize.begin()+1, vecChunkSize.end());

				varMu->set_chunk_size(vecChunkSizeOut);
				varSigma->set_chunk_size(vecChunkSizeOut);
				varXi->set_chunk_size(vecChunkSizeOut);
			}
		}

		// Create frequency variable
		if (vecFrequencies.size() != 0) {
			Announce("Initializing frequency variable");
		}

		std::string strVariableFreq = strVariable + "_freq";
		if (vecFrequencies.size() == 1) {
			varFreq =
				ncout.add_var(
					strVariableFreq.c_str(),
					ncFloat,
					vecDimOut.size(),
					const_cast<const NcDim**>(&(vecDimOut[0])));

			if (varFreq == NULL) {
				_EXCEPTION2("Cannot create variable \"%s\" in output file \"%s\"",
					strVariableFreq.c_str(), strOutputData.c_str());
			}

			varFreq->add_att("_FillValue", -std::numeric_limits<float>::max());

			std::vector<size_t> vecChunkSizeOut(vecChunkSize.begin()+1, vecChunkSize.end());
			varFreq->set_chunk_size(vecChunkSizeOut);

		} else if (vecFrequencies.size() > 1) {
			dimFreq = ncout.add_dim("Frequency", static_cast<long>(vecFrequencies.size()));

			std::vector<NcDim*> vecDimFreqOut = vecDimOut;
			vecDimFreqOut.insert(vecDimFreqOut.begin(), dimFreq);

			varFreq =
				ncout.add_var(
					strVariableFreq.c_str(),
					ncFloat,
					vecDimFreqOut.size(),
					const_cast<const NcDim**>(&(vecDimFreqOut[0])));

			if (varFreq == NULL) {
				_EXCEPTION2("Cannot create variable \"%s\" in output file \"%s\"",
					strVariableFreq.c_str(), strOutputData.c_str());
			}

			varFreq->add_att("_FillValue", -std::numeric_limits<float>::max());

			std::vector<size_t> vecChunkSizeOut(vecChunkSize.begin()+1, vecChunkSize.end());
			vecChunkSizeOut.insert(vecChunkSizeOut.begin(), 1);
			varFreq->set_chunk_size(vecChunkSizeOut);

			NcVar * varFreqCoord =
				ncout.add_var(
					"Frequency",
					ncFloat,
					dimFreq);

			varFreqCoord->put(
				&(vecFrequencies[0]),
				static_cast<long>(vecFrequencies.size()));

			varFreqCoord->add_att("units", strFrequencyUnit.c_str());
		}

		if (vecFrequencies.size() != 0) {
			CopyNcVarAttributes(varIn, varFreq);
		}

		AnnounceEndBlock(NULL);
	}

	/////////////////////////////////////////////
	// Actual calculation
	AnnounceStartBlock("Performing calculation");

	// Data
	DataArray1D<float> dDataTemp;
	DataArray2D<float> dData;
	DataArray2D<float> dIDFParameters;
	DataArray2D<float> dFrequencies;
	if (sBlockSize == 0) {
		dDataTemp.Allocate(grid.GetSize());
		dData.Allocate(grid.GetSize(), sTimeCount);
		dIDFParameters.Allocate(3, grid.GetSize());
	} else {
		dDataTemp.Allocate(sBlockSize);
		dData.Allocate(sBlockSize, sTimeCount);
		dIDFParameters.Allocate(3, sBlockSize);
	}
	if (vecFrequencies.size() != 0) {
		dFrequencies.Allocate(vecFrequencies.size(), sBlockSize);
	}

	// Loop over number of blocks
	NcVarSpatialIterator iter(grid.m_nGridDim.size(), sBlockSize, true);
	while(true) {
		AnnounceStartBlock("Block %lu", iter.GetCurrentBlockIx());

		size_t sPointsReadRef = 0;

		// Loop over all files to be processed
		size_t tGlobal = 0;
		for (size_t f = 0; f < vecInputFileList.size(); f++) {

			// Get the list of data files
			NcFileVector vecFiles;
			vecFiles.ParseFromString(vecInputFileList[f]);

			_ASSERT(vecFiles.size() > 0);

			const NcTimeDimension & vecTimes =
				vecFiles.GetNcTimeDimension(0);

			// Get the variable
			NcVar * varIn = NULL;
			vecFiles.FindContainingVariable(strVariable, &varIn);
			if (varIn == NULL) {
				_EXCEPTION2("Cannot find variable \"%s\" in file list (line %lu)",
					strVariable.c_str(), f);
			}

			// Compare chunk information
			if (sBlockSize != 0) {
				std::vector<size_t> vecThisChunkSize = varIn->get_chunk_size();
				if (vecThisChunkSize != vecChunkSize) {
					_EXCEPTION2("Chunk sizes change in file \"%s\" (line %lu);"
						" --block_size must be set to 0",
						vecFiles.GetFilename(f).c_str(), f);
				}
			}

			// Get data
			for (size_t t = 0; t < vecTimes.size(); t++,tGlobal++) {
				size_t sPointsRead = iter.GetBlockFloat(varIn, t, dDataTemp);
				if ((f == 0) && (t == 0)) {
					sPointsReadRef = sPointsRead;
				} else if (sPointsRead != sPointsReadRef) {
					_EXCEPTION4("Block size mismatch in file \"%s\" time %lu (%lu != %lu);"
						" chunking structure changes between files",
						vecFiles.GetFilename(f).c_str(), t, sPointsRead, sPointsReadRef);
				}

				//printf("Loaded %lu points\n", sPoints);

				// Store in transposed array
				for (size_t i = 0; i < sPointsRead; i++) {
					dData(i,tGlobal) = dDataTemp[i];
				}
			}
		}

		// Calculate IDF
		CalculateIDF(dData, dIDFParameters, sPointsReadRef, dFillValue);

		// Write data
		if (fOutputParameters) {
			size_t sPointsWritten;
			sPointsWritten = iter.PutBlockFloat(varMu, -1, &(dIDFParameters(0,0)));
			if (sPointsWritten != sPointsReadRef) {
				_EXCEPTIONT("Logic error: mismatch in points read/written");
			}
			sPointsWritten = iter.PutBlockFloat(varSigma, -1, &(dIDFParameters(1,0)));
			if (sPointsWritten != sPointsReadRef) {
				_EXCEPTIONT("Logic error: mismatch in points read/written");
			}
			sPointsWritten = iter.PutBlockFloat(varXi, -1, &(dIDFParameters(2,0)));
			if (sPointsWritten != sPointsReadRef) {
				_EXCEPTIONT("Logic error: mismatch in points read/written");
			}
		}

		// Compute return levels and write
		for (size_t fr = 0; fr < vecFrequencies.size(); fr++) {
			for (size_t i = 0; i < sPointsReadRef; i++) {
				dFrequencies(fr,i) =
					IDFReturnLevel(
						dIDFParameters(0,i),
						dIDFParameters(1,i),
						dIDFParameters(2,i),
						vecFrequencies[fr]);
			}

			size_t sPointsWritten;
			if (vecFrequencies.size() == 1) {
				sPointsWritten = iter.PutBlockFloat(varFreq, -1, &(dFrequencies(0,0)));
				if (sPointsWritten != sPointsReadRef) {
					_EXCEPTIONT("Logic error: mismatch in points read/written");
				}
			} else {
				sPointsWritten = iter.PutBlockFloat(varFreq, fr, &(dFrequencies(fr,0)));
				if (sPointsWritten != sPointsReadRef) {
					_EXCEPTIONT("Logic error: mismatch in points read/written");
				}
			}
		}

		// Next block
		if(!iter.NextBlock()) {
			break;
		}

		AnnounceEndBlock(NULL);
	}

	AnnounceEndBlock("Done");
	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

}

///////////////////////////////////////////////////////////////////////////////

