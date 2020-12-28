///////////////////////////////////////////////////////////////////////////////
///
///	\file    IntegrateDimension.cpp
///	\author  Paul Ullrich
///	\version December 26th, 2020
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
#include "NcFileVector.h"

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
		_EXCEPTIONT("At present IntegrateDimension only supports serial execution.");
	}
#endif

	// Input data file
	std::string strInputFiles;

	// Output data file
	std::string strOutputFile;

	// Variable to apply Fourier transform to
	std::string strVarName;

	// Variables to preserve
	std::string strPreserveVarName;

	// Dimension along which to perform Fourier transform
	std::string strDimName;

	// Variable name for surface pressure
	std::string strPSVariableName;

	// Reference pressure
	double dReferencePressure;

	// Scale factor
	double dScaleFactor;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFiles, "in_data", "");
		CommandLineString(strOutputFile, "out_data", "");
		CommandLineString(strVarName, "var", "");
		CommandLineString(strPreserveVarName, "preserve", "");
		CommandLineString(strDimName, "dim", "lev");
		CommandLineString(strPSVariableName, "psvarname", "PS");
		CommandLineDouble(dReferencePressure, "psref", 1.0e5);
		CommandLineDouble(dScaleFactor, "scale", 1.0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Validate arguments
	if (strInputFiles.length() == 0) {
		_EXCEPTIONT("No input data file(s) (--in_data) specified");
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

	// Parse input file list (--in_data)
	NcFileVector vecInputFileList;
	vecInputFileList.ParseFromString(strInputFiles);
	_ASSERT(vecInputFileList.size() != 0);

	// Parse variable list (--var)
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVarName, vecVariableStrings);

	// Parse preserve list (--preserve)
	std::vector<std::string> vecPreserveVariableStrings;

	vecPreserveVariableStrings.push_back("time");
	vecPreserveVariableStrings.push_back("lon");
	vecPreserveVariableStrings.push_back("lat");

	STLStringHelper::ParseVariableList(strPreserveVarName, vecPreserveVariableStrings);

	// Begin processing
	AnnounceBanner();
	AnnounceStartBlock("Initializing output file");

	// Open the output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\"",
			strOutputFile.c_str());
	}

	AnnounceEndBlock("Done");

	// Copy preserve variables
	if (vecPreserveVariableStrings.size() != 0) {
		AnnounceStartBlock("Preserving variables");

		for (int v = 0; v < vecPreserveVariableStrings.size(); v++) {
			Announce("Variable %s", vecPreserveVariableStrings[v].c_str());

			NcVar * var;
			size_t sFileIx = vecInputFileList.FindContainingVariable(vecPreserveVariableStrings[v], &var);
			if (var == NULL) {
				_EXCEPTION1("Unable to find variable \"%s\" in specified files",
					vecPreserveVariableStrings[v].c_str());
			}

			NcFile * ncinfile = vecInputFileList[sFileIx];
			_ASSERT(ncinfile != NULL);

			CopyNcVarIfExists(
				(*ncinfile),
				ncoutfile,
				vecPreserveVariableStrings[v]);
		}

		AnnounceEndBlock("Done");
	}

	// Load hybrid coefficients (if present)
	DataArray1D<double> dHYAI;
	DataArray1D<double> dHYBI;

	{
		// Load in hybrid coefficients
		NcVar * varHYAI;
		NcVar * varHYBI;

		size_t sFileIx = vecInputFileList.FindContainingVariable("hyai", &varHYAI);
		if (varHYAI != NULL) {
			NcFile * ncinfile = vecInputFileList[sFileIx];
			_ASSERT(ncinfile != NULL);

			varHYBI = ncinfile->get_var("hybi");
			if (varHYBI == NULL) {
				Announce("WARNING: File \"%s\" contains variable \"hyai\" but not \"hybi\"; unable to load hybrid coefficients",
					vecInputFileList.GetFilename(sFileIx).c_str());
			} else if (varHYAI->num_dims() != 1) {
				Announce("WARNING: File \"%s\" variable \"hyai\" has more than one dimension; unable to load hybrid coefficients",
					vecInputFileList.GetFilename(sFileIx).c_str());
			} else if (varHYBI->num_dims() != 1) {
				Announce("WARNING: File \"%s\" variable \"hybi\" has more than one dimension; unable to load hybrid coefficients",
					vecInputFileList.GetFilename(sFileIx).c_str());
			} else if (varHYAI->get_dim(0)->size() < 2) {
				Announce("WARNING: File \"%s\" variable \"hyai\" must have size greater than 1; unable to load hybrid coefficients",
					vecInputFileList.GetFilename(sFileIx).c_str());
			} else if (varHYAI->get_dim(0)->size() != varHYBI->get_dim(0)->size()) {
				Announce("WARNING: File \"%s\" variable \"hyai\" has size %li but \"hybi\" has size %li; unable to load hybrid coefficients",
					vecInputFileList.GetFilename(sFileIx).c_str(),
					varHYAI->get_dim(0)->size(),
					varHYBI->get_dim(0)->size());
			} else {
				dHYAI.Allocate(varHYAI->get_dim(0)->size());
				dHYBI.Allocate(varHYBI->get_dim(0)->size());

				varHYAI->set_cur((long)0);
				varHYAI->get(&(dHYAI[0]), (long)dHYAI.GetRows());
				varHYBI->set_cur((long)0);
				varHYBI->get(&(dHYBI[0]), (long)dHYBI.GetRows());
			}

			// Verify monotonicity of hybrid pressure coefficients
			// dSign indicates direction of coordinate (-1.0 = low to high)
			double dFirstSign = ((dHYAI[1] + dHYBI[1]) > (dHYAI[0] + dHYBI[0]))?(-1.0):(+1.0);
			for (size_t i = 0; i < dHYAI.GetRows()-1; i++) {
				double dThisSign = ((dHYAI[i+1] + dHYBI[i+1]) > (dHYAI[i] + dHYBI[i]))?(-1.0):(+1.0);
				if (dFirstSign != dThisSign) {
					Announce("WARNING: File \"%s\" variable \"hyai\" and \"hybi\" are non-monotone; unable to load hybrid coefficients",
						vecInputFileList.GetFilename(sFileIx).c_str());
				}
			}

			// Change the scale factor to reflect if data is top -> bottom so
			// integral is always from the surface to top of atmosphere.
			dScaleFactor *= dFirstSign;
		}
	}

	_ASSERT(dHYAI.GetRows() == dHYBI.GetRows());
	if (dHYAI.GetRows() == 0) {
		_EXCEPTIONT("NOT IMPLEMENTED: Currently this executable is only for computing integrals over hybrid coordinates");
	}

	// Load PS variable
	NcVar * varPS;
	{
		size_t sFileIx = vecInputFileList.FindContainingVariable(strPSVariableName, &varPS);
		if (varPS == NULL) {
			_EXCEPTION1("Variable \"%s\" not found among data files",
				strPSVariableName.c_str());
		}
	}

	// Begin processing
	AnnounceStartBlock("Processing");
	for (int v = 0; v < vecVariableStrings.size(); v++) {

		AnnounceStartBlock("Variable \"%s\"", vecVariableStrings[v].c_str());

		// Find the file containing this variable
		NcVar * varIn;
		size_t sFileIx = vecInputFileList.FindContainingVariable(vecVariableStrings[v], &varIn);
		if (varIn == NULL) {
			_EXCEPTION1("Unable to find variable \"%s\" among input files",
				vecVariableStrings[v].c_str());
		}

		// Find the integral dimension
		long lIntegralDimSize = 0;
		long lIntegralDimIx = (-1);
		for (long d = 0; d < varIn->num_dims(); d++) {
			if (strDimName == std::string(varIn->get_dim(d)->name())) {
				lIntegralDimIx = d;
				break;
			}
		}
		if (lIntegralDimIx == (-1)) {
			_EXCEPTION3("Variable \"%s\" in file \"%s\" does not contain dimension \"%s\"",
				vecVariableStrings[v].c_str(),
				vecInputFileList.GetFilename(sFileIx).c_str(),
				strDimName.c_str());
		}
		lIntegralDimSize = varIn->get_dim(lIntegralDimIx)->size();

		if (dHYAI.GetRows()-1 != lIntegralDimSize) {
			_EXCEPTION2("Integral dimension size (%li) must be one less than interface hybrid coefficient array size (%li)",
				lIntegralDimSize,
				dHYAI.GetRows());
		}

		// Get the number of auxiliary dimensions (dimensions preceding the integral dimension)
		long lAuxSize = 1;
		std::vector<long> vecAuxDimSize;
		if (lIntegralDimIx != 0) {
			vecAuxDimSize.resize(lIntegralDimIx);
			for (long d = 0; d < lIntegralDimIx; d++) {
				vecAuxDimSize[d] = varIn->get_dim(d)->size();
				lAuxSize *= vecAuxDimSize[d];
			}
		}

		// Get the number of grid dimensions (dimensions after the integral dimension)
		long lGridSize = 1;
		std::vector<long> vecGridDimSize;
		if (lIntegralDimIx != varIn->num_dims()-1) {
			vecGridDimSize.resize(varIn->num_dims() - lIntegralDimIx - 1);
			for (long d = lIntegralDimIx+1; d < varIn->num_dims(); d++) {
				vecGridDimSize[d-lIntegralDimIx-1] = varIn->get_dim(d)->size();
				lGridSize *= vecGridDimSize[d-lIntegralDimIx-1];
			}
		}

		// Verify PS auxiliary dimension count matches
		_ASSERT(varPS != NULL);
		if (varPS->num_dims() != varIn->num_dims()-1) {
			_EXCEPTION4("Dimension mismatch: Variable \"%s\" has %li dimensions, but \"%s\" has %li dimensions (should be 1 less)",
				vecVariableStrings[v].c_str(),
				varIn->num_dims(),
				strPSVariableName.c_str(),
				varPS->num_dims());
		}
		for (long d = 0; d < vecAuxDimSize.size(); d++) {
			if (varPS->get_dim(d)->size() != varIn->get_dim(d)->size()) {
				_EXCEPTION6("Dimension mismatch: Variable \"%s\" dimension %li has size %li, but \"%s\" dimension %li has size %li",
					vecVariableStrings[v].c_str(),
					d,
					varIn->get_dim(d)->size(),
					strPSVariableName.c_str(),
					d,
					varPS->get_dim(d)->size());
			}
		}
		for (long d = 0; d < vecGridDimSize.size(); d++) {
			long dPS = d + vecAuxDimSize.size();
			long dIn = d + vecAuxDimSize.size() + 1;
			if (varPS->get_dim(dPS)->size() != varIn->get_dim(dIn)->size()) {
				_EXCEPTION6("Dimension mismatch: Variable \"%s\" dimension %li has size %li, but \"%s\" dimension %li has size %li",
					vecVariableStrings[v].c_str(),
					dIn,
					varIn->get_dim(dIn)->size(),
					strPSVariableName.c_str(),
					dPS,
					varPS->get_dim(dPS)->size());
			}
		}

		// Copy dimensions
		std::vector<NcDim *> vecDimOut;
		for (long d = 0; d < varIn->num_dims(); d++) {
			if (d == lIntegralDimIx) {
				continue;
			}
			NcDim * dimIn = varIn->get_dim(d);
			NcDim * dimOut = ncoutfile.get_dim(dimIn->name());
			if (dimOut != NULL) {
				if (dimOut->size() != dimIn->size()) {
					_EXCEPTION3("Dimension \"%s\" has incompatible size in input (%li) and output (%li)",
						dimIn->name(), dimIn->size(), dimOut->size());
				}
			} else {
				dimOut = ncoutfile.add_dim(dimIn->name(), dimIn->size());
				if (dimOut == NULL) {
					_EXCEPTION1("Unable to create dimension \"%s\" in output file",
						dimIn->name());
				}
			}
			vecDimOut.push_back(dimOut);
		}

		// Create output variable
		NcVar * varOut =
			ncoutfile.add_var(
				varIn->name(),
				ncFloat,
				vecDimOut.size(),
				const_cast<const NcDim**>(&(vecDimOut[0])));
		if (varOut == NULL) {
			_EXCEPTION1("Unable to create variable \"%s\" in output file",
				varIn->name());
		}

		// Allocate data
		DataArray1D<double> dDataPS(lGridSize);
		DataArray1D<double> dDataVar(lGridSize);
		DataArray1D<double> dDataOut(lGridSize);

		// Loop through auxiliary dimensions
		for (long lAux = 0; lAux < lAuxSize; lAux++) {

			dDataOut.Zero();

			// Get the data index
			std::vector<long> lPos(varIn->num_dims(), 0);
			std::vector<long> lSize(varIn->num_dims(), 1);
			long lAuxTemp = lAux;

			for (long d = vecAuxDimSize.size()-1; d >= 0; d--) {
				lPos[d] = lAuxTemp % vecAuxDimSize[d];
				lAuxTemp /= vecAuxDimSize[d];
			}
			for (long d = 0; d < vecGridDimSize.size(); d++) {
				lSize[lIntegralDimIx+d+1] = vecGridDimSize[d];
			}

			if (lAuxSize != 1) {
				char szPos[10];
				std::string strPos;
				for (long d = 0; d < vecAuxDimSize.size(); d++) {
					sprintf(szPos, "%li", lPos[d]);
					strPos += szPos;
					if (d != vecAuxDimSize.size()-1) {
						strPos += ",";
					}
				}
				Announce("Processing %s(%s)",
					vecVariableStrings[v].c_str(),
					strPos.c_str());
			}

			// Load PS data
			std::vector<long> lPosPS = lPos;
			std::vector<long> lSizePS = lSize;
			lPosPS.erase(lPosPS.begin() + lIntegralDimIx);
			lSizePS.erase(lSizePS.begin() + lIntegralDimIx);
			_ASSERT(lPosPS.size() == varPS->num_dims());
			_ASSERT(lSizePS.size() == varPS->num_dims());

			varPS->set_cur(&(lPosPS[0]));
			varPS->get(&(dDataPS[0]), &(lSizePS[0]));

			// Loop through integral dimension
			_ASSERT(lPos.size() > lIntegralDimIx);
			_ASSERT(varIn->get_dim(lIntegralDimIx)->size() == lIntegralDimSize);
			for (long lLev = 0; lLev < lIntegralDimSize; lLev++) {
				lPos[lIntegralDimIx] = lLev;
				varIn->set_cur(&(lPos[0]));
				varIn->get(&(dDataVar[0]), &(lSize[0]));

				for (long lGrid = 0; lGrid < lGridSize; lGrid++) {
					double dPl = dHYAI[lLev] * dReferencePressure + dHYBI[lLev] * dDataPS[lGrid];
					double dPu = dHYAI[lLev+1] * dReferencePressure + dHYBI[lLev+1] * dDataPS[lGrid];

					dDataOut[lGrid] += (dPu - dPl) * dDataVar[lGrid];
				}
			}

			if (dScaleFactor != 1.0) {
				for (long lGrid = 0; lGrid < lGridSize; lGrid++) {
					dDataOut[lGrid] *= dScaleFactor;
				}
			}

			// Write output data
			varOut->set_cur(&(lPosPS[0]));
			varOut->put(&(dDataOut[0]), &(lSizePS[0]));
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

