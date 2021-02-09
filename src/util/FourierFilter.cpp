///////////////////////////////////////////////////////////////////////////////
///
///	\file    FourierFilter.cpp
///	\author  Paul Ullrich
///	\version June 3, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
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
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "DataArray1D.h"
#include "FourierTransforms.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>

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
		_EXCEPTIONT("At present FourierFilter only supports serial execution.");
	}
#endif

	// Input data file
	std::string strInputFile;

	// Output data file
	std::string strOutputFile;

	// Variable to apply Fourier transform to
	std::string strVarName;

	// Variables to preserve
	std::string strPreserveVarName;

	// Dimension along which to perform Fourier transform
	std::string strDimName;

	// Number of Fourier modes
	int nFourierModes;

	// Output phase and magnitude
	bool fOutputCoeffs;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strOutputFile, "out_data", "");
		CommandLineString(strVarName, "var", "");
		CommandLineString(strPreserveVarName, "preserve", "");
		CommandLineString(strDimName, "dim", "");
		CommandLineInt(nFourierModes, "modes", 4);
		CommandLineBool(fOutputCoeffs, "output_coeffs");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Validate arguments
	if (strInputFile.length() == 0) {
		_EXCEPTIONT("No input data file (--in_data) specified");
	}
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output data file (--out_data) specified");
	}
	if (strVarName.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}
	if (strDimName.length() == 0) {
		_EXCEPTIONT("No dimension name (--dim) specified");
	}
	if (nFourierModes < 1) {
		_EXCEPTIONT("--modes must be a positive integer");
	}

	// Parse variable list
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVarName, vecVariableStrings);

	// Prase preserve list
	std::vector<std::string> vecPreserveVariableStrings;
	STLStringHelper::ParseVariableList(strPreserveVarName, vecPreserveVariableStrings);

	// Begin processing
	AnnounceStartBlock("Initializing output file");

	// Open the input file
	NcFile ncinfile(strInputFile.c_str());
	if (!ncinfile.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for reading",
			strInputFile.c_str());
	}

	// Open the output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\"",
			strOutputFile.c_str());
	}

	// Copy attributes
	CopyNcFileAttributes(
		&ncinfile,
		&ncoutfile);

	// Loop through all dimensions in the input file
	NcDim * dimPrimary = ncinfile.get_dim(strDimName.c_str());
	if (dimPrimary == NULL) {
		_EXCEPTION2("Input NetCDF file \"%s\" missing dimension \"%s\"",
			strInputFile.c_str(),
			strDimName.c_str());
	}

	AnnounceEndBlock("Done");

	// Copy preserve variables
	if (strPreserveVarName != "") {
		AnnounceStartBlock("Preserving variables");

		for (int v = 0; v < vecPreserveVariableStrings.size(); v++) {
			Announce("Variable %s", vecPreserveVariableStrings[v].c_str());
			CopyNcVar(
				ncinfile,
				ncoutfile,
				vecPreserveVariableStrings[v]);
		}

		AnnounceEndBlock("Done");
	}

	// Copy variables from input file into output file
	AnnounceStartBlock("Processing");
	for (int v = 0; v < ncinfile.num_vars(); v++) {
		NcVar * varIn = ncinfile.get_var(v);

		// Check that the variable is in the variable list
		bool fFound = false;
		for (int s = 0; s < vecVariableStrings.size(); s++) {
			if (vecVariableStrings[s] == varIn->name()) {
				fFound = true;
			}
		}
		if (!fFound) {
			continue;
		}

		AnnounceStartBlock("Variable \"%s\"", varIn->name());

		// Check if the variable has the Fourier dimension
		int iFourierDimension = (-1);
		for (int d = 0; d < varIn->num_dims(); d++) {
			NcDim * dim = varIn->get_dim(d);
			if (strDimName == dim->name()) {
				iFourierDimension = d;
				break;
			}
		}

		CopyNcVar(
			ncinfile,
			ncoutfile,
			varIn->name(),
			true,
			iFourierDimension != (-1));

		if (iFourierDimension == (-1)) {
			Announce("Variable has no dimension \"%s\"; not processed", strDimName.c_str());
			AnnounceEndBlock(NULL);
			continue;
		}

		// Vector of dimensions for phase and magnitude outputs
		std::vector<long> vecdimCoeffOutSize;
		std::vector<NcDim*> vecdimCoeffOut;
		if (fOutputCoeffs) {
			NcDim * dimMode = ncoutfile.add_dim("fourier_mode", nFourierModes);
			if (dimMode == NULL) {
				_EXCEPTIONT("Unable to create dimension \"fourier_mode\" in output file");
			}
			vecdimCoeffOut.push_back(dimMode);
			vecdimCoeffOutSize.push_back(nFourierModes);
			for (int d = 0; d < varIn->num_dims(); d++) {
				NcDim * dimIn = varIn->get_dim(d);
				if (strDimName == dimIn->name()) {
					continue;
				}
				NcDim * dimOut = ncoutfile.get_dim(dimIn->name());
				if (dimOut == NULL) {
					_EXCEPTION1("Missing dimension \"%s\" in output file",
						dimIn->name());
				}
				vecdimCoeffOut.push_back(dimOut);
				vecdimCoeffOutSize.push_back(dimOut->size());
			}
		}

		// Apply Fourier filter
		Announce("Applying Fourier filter to dimension %i", iFourierDimension);

		// Load in data and Fourier filter
		if ((varIn->type() != ncFloat) && (varIn->type() != ncDouble)) {
			_EXCEPTIONT("Only variables of type \"float\" or \"double\" currently supported");
		}

		NcVar * varOut = ncoutfile.get_var(varIn->name());
		if (varOut == NULL) {
			_EXCEPTION1("Unable to create variable \"%s\" in output file.",
				varIn->name());
		}

		size_t sTotalSize = 1;
		size_t sFourierStride = 1;
		size_t sFourierCount;
		size_t sFourierOuterLoops = 1;
		std::vector<long> lSize;
		for (int d = 0; d < varOut->num_dims(); d++) {
			long lDimSize = varOut->get_dim(d)->size();
			lSize.push_back(lDimSize);
			sTotalSize *= static_cast<size_t>(lDimSize);

			if (d > iFourierDimension) {
				sFourierStride *= static_cast<size_t>(lDimSize);
			} else if (d == iFourierDimension) {
				sFourierCount = static_cast<size_t>(lDimSize);
			} else {
				sFourierOuterLoops *= static_cast<size_t>(lDimSize);
			}
		}

		// Perform filtering on variables of type float
		if (varIn->type() == ncFloat) {

			// Output phase and magnitude
			DataArray1D<float> dMagnitude;
			DataArray1D<float> dPhase;

			NcVar * varMagnitude;
			NcVar * varPhase;

			if (fOutputCoeffs) {
				dMagnitude.Allocate(nFourierModes * sTotalSize);
				dPhase.Allocate(nFourierModes * sTotalSize);

				std::string strMagnitudeName = std::string(varIn->name()) + "_mag";
				varMagnitude =
					ncoutfile.add_var(
						strMagnitudeName.c_str(),
						ncFloat,
						vecdimCoeffOut.size(),
						const_cast<const NcDim**>(&(vecdimCoeffOut[0])));

				std::string strPhaseName = std::string(varIn->name()) + "_phase";
				varPhase =
					ncoutfile.add_var(
							strPhaseName.c_str(),
							ncFloat,
							vecdimCoeffOut.size(),
							const_cast<const NcDim**>(&(vecdimCoeffOut[0])));
			}

			DataArray1D<float> data(sTotalSize);
			varIn->get(&(data[0]), &(lSize[0]));

			DataArray1D<float> an(sFourierCount);
			DataArray1D<float> bn(sFourierCount);

			//std::cout << sFourierCount << std::endl;
			//std::cout << sFourierStride << std::endl;
			for (size_t s = 0; s < sFourierOuterLoops; s++) {
			for (size_t i = 0; i < sFourierStride; i++) {
				size_t sDataIx = s * sFourierStride * sFourierCount + i;

				fourier_filter<float>(
					&(data[sDataIx]),
					sFourierCount,
					sFourierStride,
					static_cast<size_t>(nFourierModes),
					an, bn);

				if (fOutputCoeffs) {
					for (int m = 0; m < nFourierModes; m++) {
						size_t sModeSize = sTotalSize / sFourierCount;
						dMagnitude[m * sModeSize + sDataIx] =
							sqrt(an[m] * an[m] + bn[m] * bn[m]) / static_cast<float>(sFourierCount);
						dPhase[m * sModeSize + sDataIx] =
							atan2(bn[m], an[m]);
					}
				}
			}
			}

			varOut->put(&(data[0]), &(lSize[0]));

			if (fOutputCoeffs) {
				varMagnitude->put(&(dMagnitude[0]), &(vecdimCoeffOutSize[0]));
				varPhase->put(&(dPhase[0]), &(vecdimCoeffOutSize[0]));
			}

		// Perform filtering on variables of type double
		} else {
			DataArray1D<double> data(sTotalSize);
			varIn->get(&(data[0]), &(lSize[0]));

			DataArray1D<double> an(sFourierCount);
			DataArray1D<double> bn(sFourierCount);

			//std::cout << sFourierCount << std::endl;
			//std::cout << sFourierStride << std::endl;
			for (size_t s = 0; s < sFourierOuterLoops; s++) {
			for (size_t i = 0; i < sFourierStride; i++) {
				fourier_filter<double>(
					&(data[s * sFourierStride * sFourierCount + i]),
					sFourierCount,
					sFourierStride,
					static_cast<size_t>(nFourierModes),
					an, bn);
			}
			}

			varOut->put(&(data[0]), &(lSize[0]));
		}
		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

