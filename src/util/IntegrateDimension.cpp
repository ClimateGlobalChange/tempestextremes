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
#include "Constants.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "NcFileVector.h"
#include "MathExpression.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

template <typename T>
class FieldUnion {

public:
	///	<summary>
	///		Type stored in this FieldUnion.
	///	</summary>
	enum class Type {
		Unknown,
		Scalar,
		BoundsVector,
		Field
	};

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	FieldUnion() :
		type(Type::Unknown),
		fieldvar(NULL)
	{ }

public:
	///	<summary>
	///		Type being stored.
	///	</summary>
	Type type;

	///	<summary>
	///		Scalar value.
	///	</summary>
	T scalar;

	///	<summary>
	///		Interfacial bounds.
	///	</summary>
	DataArray2D<T> bounds;

	///	<summary>
	///		Field.
	///	</summary>
	NcVar * fieldvar;

	///	<summary>
	///		Field data.
	///	</summary>
	DataArray1D<T> fielddata;
};

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

	// Variable to integrate
	std::string strVarName;

	// Names for output variables after integration
	std::string strVarOutName;

	// Variables to preserve
	std::string strPreserveVarName;

	// Dimension along which to perform integration
	std::string strDimName;

	// Variable name for surface pressure
	std::string strPSVariableName;

	// Variable name for reference pressure
	std::string strHybridExpr;

	// Type of hybrid coordinate
	std::string strHybridCoordType;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFiles, "in_data", "");
		CommandLineString(strOutputFile, "out_data", "");
		CommandLineString(strVarName, "var", "");
		CommandLineString(strVarOutName, "varout", "");
		CommandLineString(strPreserveVarName, "preserve", "");
		CommandLineString(strDimName, "dim", "lev");
		CommandLineString(strHybridExpr, "hybridexpr", "a*p0+b*ps");
		CommandLineStringD(strHybridCoordType, "hybridtype", "p", "[p|z]");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

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
	if ((strHybridCoordType != "p") && (strHybridCoordType != "z")) {
		_EXCEPTIONT("Hybrid coordinate type (--hytype) must be either \"p\" or \"z\"");
	}

	// Parse input file list (--in_data)
	NcFileVector vecInputFileList;
	vecInputFileList.ParseFromString(strInputFiles);
	_ASSERT(vecInputFileList.size() != 0);

	// Parse variable list (--var)
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVarName, vecVariableStrings);

	// Parse output variable list (--outvar)
	std::vector<std::string> vecOutputVariableStrings;
	if (strVarOutName.length() == 0) {
		vecOutputVariableStrings = vecVariableStrings;
	} else {
		STLStringHelper::ParseVariableList(strVarOutName, vecOutputVariableStrings);
		if (vecVariableStrings.size() != vecOutputVariableStrings.size()) {
			_EXCEPTION2("Inconsistent number of variables in --var (%lu) and --varout (%lu)",
				vecVariableStrings.size(), vecOutputVariableStrings.size());
		}
	}

	// Parse preserve list (--preserve)
	std::vector<std::string> vecPreserveVariableStrings;

	vecPreserveVariableStrings.push_back("time");
	vecPreserveVariableStrings.push_back("lon");
	vecPreserveVariableStrings.push_back("lat");

	STLStringHelper::ParseVariableList(strPreserveVarName, vecPreserveVariableStrings);

	// Parse the hybrid expression
	enum class HybridExprType {
		ValuePlusDoubleProduct,
		DoubleProductPlusDoubleProduct,
	};

	MathExpression exprHybridExpr(strHybridExpr);

	// Begin processing
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

	// Build another array that stores the values of the hybrid coordinate array
	AnnounceStartBlock("Loading interfacial variables");
	std::vector< FieldUnion<double> > vecExprContents(exprHybridExpr.size());
	for(size_t t = 0; t < exprHybridExpr.size(); t++) {
		const MathExpression::Token & token = exprHybridExpr[t];
		if (token.type == MathExpression::Token::Type::Variable) {
			NcVar * var;

			std::string strBnds = token.str + "_bnds";
			std::string strHy = "hy" + token.str + "i";

			size_t sFileIx = vecInputFileList.FindContainingVariable(token.str, &var);

			// Interfacial variables need to be of the form $_bnds or hy$i
			if (sFileIx != NcFileVector::InvalidIndex) {
				if (var->num_dims() == 1) {
					sFileIx = NcFileVector::InvalidIndex;
				}
				if ((var->num_dims() == 2) && (var->get_dim(1)->size() == 2)) {
					sFileIx = NcFileVector::InvalidIndex;
				}
			}

			if (sFileIx == NcFileVector::InvalidIndex) {
				sFileIx = vecInputFileList.FindContainingVariable(strBnds, &var);
			}

			if (sFileIx == NcFileVector::InvalidIndex) {
				sFileIx = vecInputFileList.FindContainingVariable(strHy, &var);
			}

			if (sFileIx == NcFileVector::InvalidIndex) {
				_EXCEPTION3("Cannot file variables \"%s\" or \"%s\" in input files (or missing corresponding interfacial variables)",
					token.str.c_str(), strBnds.c_str(), strHy.c_str());
			}

			_ASSERT(var != NULL);

			// Constants
			if (var->num_dims() == 0) {
				Announce("%s (constant)", var->name());
				double dValue;
				var->get(&dValue, 1);

				vecExprContents[t].type = FieldUnion<double>::Type::Scalar;
				vecExprContents[t].scalar = dValue;

			// Interfacial values (single array)
			} else if (var->num_dims() == 1) {
				Announce("%s (1D bounds)", var->name());

				long lDimSize = var->get_dim(0)->size();
				if (lDimSize < 2) {
					_EXCEPTION2("Interfacial variable \"%s\" dimension \"%s\" must have size >= 2",
						var->name(), var->get_dim(0)->name());
				}

				vecExprContents[t].type = FieldUnion<double>::Type::BoundsVector;
				vecExprContents[t].bounds.Allocate(lDimSize-1, 2);

				DataArray1D<double> data(lDimSize);
				var->get(&(data[0]), lDimSize);
				for (long l = 0; l < lDimSize-1; l++) {
					vecExprContents[t].bounds(l,0) = data(l);
					vecExprContents[t].bounds(l,1) = data(l+1);
				}

			// Interfacial values (bounds)
			} else if ((var->num_dims() == 2) && (var->get_dim(1)->size() == 2)) {
				Announce("%s (2D bounds)", var->name());

				long lDimSize = var->get_dim(0)->size();
				if (lDimSize < 1) {
					_EXCEPTION2("Interfacial variable \"%s\" dimension \"%s\" must have size >= 1",
						var->name(), var->get_dim(0)->name());
				}

				vecExprContents[t].type = FieldUnion<double>::Type::BoundsVector;
				vecExprContents[t].bounds.Allocate(lDimSize, 2);

				var->get(&(vecExprContents[t].bounds(0,0)), lDimSize, 2);

				//for (int i = 0; i < (*parray).GetRows(); i++) {
				//for (int j = 0; j < (*parray).GetColumns(); j++) {
				//	printf("(%i %i) %1.15f\n", i, j, (*parray)(i,j));
				//}
				//}

			// Interfacial values (bounds) in reverse order (why?)
			} else if ((var->num_dims() == 2) && (var->get_dim(0)->size() == 2)) {
				Announce("%s (2D bounds)", var->name());
				Announce("WARNING: \"bnds\" dimension appears first. This may indicate something is incorrect with your vertical level array ordering");

				long lDimSize = var->get_dim(1)->size();
				if (lDimSize < 1) {
					_EXCEPTION2("Interfacial variable \"%s\" dimension \"%s\" must have size >= 1",
						var->name(), var->get_dim(0)->name());
				}

				vecExprContents[t].type = FieldUnion<double>::Type::BoundsVector;
				vecExprContents[t].bounds.Allocate(lDimSize, 2);

				DataArray2D<double> dBounds(lDimSize, 2);

				var->get(&(dBounds(0,0)), 2, lDimSize);

				for (int i = 0; i < lDimSize; i++) {
					vecExprContents[t].bounds(i,0) = dBounds(i,0);
					vecExprContents[t].bounds(i,1) = dBounds(i,1);
				}

			// Fields
			} else {
				Announce("%s (field)", token.str.c_str());

				vecExprContents[t].type = FieldUnion<double>::Type::Field;
				vecExprContents[t].fieldvar = var;
			}
		}
	}
	AnnounceEndBlock("Done");
/*
	_EXCEPTION();

	double dScaleFactor = 1.0;
	double dReferencePressure = 1.0;

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
*/
	// Begin processing
	AnnounceStartBlock("Processing");
	_ASSERT(vecVariableStrings.size() == vecOutputVariableStrings.size());
	for (int v = 0; v < vecVariableStrings.size(); v++) {

		if (vecVariableStrings[v] == vecOutputVariableStrings[v]) {
			AnnounceStartBlock("Variable \"%s\"",
				vecVariableStrings[v].c_str());
		} else {
			AnnounceStartBlock("Variable \"%s\" -> \"%s\"",
				vecVariableStrings[v].c_str(),
				vecOutputVariableStrings[v].c_str());
		}

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
/*
		if (dHYAI.GetRows()-1 != lIntegralDimSize) {
			_EXCEPTION2("Integral dimension size (%li) must be one less than interface hybrid coefficient array size (%li)",
				lIntegralDimSize,
				dHYAI.GetRows());
		}
*/
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

		// Verify auxiliary dimension count matches field variables
		_ASSERT(exprHybridExpr.size() == vecExprContents.size());
		for (int t = 0; t < vecExprContents.size(); t++) {

			if (vecExprContents[t].type == FieldUnion<double>::Type::BoundsVector) {
				if (vecExprContents[t].bounds.GetRows() != lIntegralDimSize) {
					_EXCEPTION3("Mismatch in vertical bounds for variable \"%s\" (%lu != %lu)",
						exprHybridExpr[t].str.c_str(),
						vecExprContents[t].bounds.GetRows(),
						lIntegralDimSize);
				}
			}

			if (vecExprContents[t].type != FieldUnion<double>::Type::Field) {
				continue;
			}

			std::string strF = exprHybridExpr[t].str;
			NcVar * varF = vecExprContents[t].fieldvar;

			if (varF->num_dims() != varIn->num_dims()-1) {
				_EXCEPTION4("Dimension mismatch: Variable \"%s\" has %li dimensions, but \"%s\" has %li dimensions (should be 1 less)",
					vecVariableStrings[v].c_str(),
					varIn->num_dims(),
					strF.c_str(),
					varF->num_dims());
			}
			for (long d = 0; d < vecAuxDimSize.size(); d++) {
				if (varF->get_dim(d)->size() != varIn->get_dim(d)->size()) {
					_EXCEPTION6("Dimension mismatch: Variable \"%s\" dimension %li has size %li, but \"%s\" dimension %li has size %li",
						vecVariableStrings[v].c_str(),
						d,
						varIn->get_dim(d)->size(),
						strF.c_str(),
						d,
						varF->get_dim(d)->size());
				}
			}
			for (long d = 0; d < vecGridDimSize.size(); d++) {
				long dPS = d + vecAuxDimSize.size();
				long dIn = d + vecAuxDimSize.size() + 1;
				if (varF->get_dim(dPS)->size() != varIn->get_dim(dIn)->size()) {
					_EXCEPTION6("Dimension mismatch: Variable \"%s\" dimension %li has size %li, but \"%s\" dimension %li has size %li",
						vecVariableStrings[v].c_str(),
						dIn,
						varIn->get_dim(dIn)->size(),
						strF.c_str(),
						dPS,
						varF->get_dim(dPS)->size());
				}
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
				vecOutputVariableStrings[v].c_str(),
				ncFloat,
				vecDimOut.size(),
				const_cast<const NcDim**>(&(vecDimOut[0])));
		if (varOut == NULL) {
			_EXCEPTION1("Unable to create variable \"%s\" in output file",
				varIn->name());
		}

		// Allocate data
		DataArray1D<double> dDataVar(lGridSize);
		DataArray1D<double> dDataOut(lGridSize);

		for (int t = 0; t < vecExprContents.size(); t++) {
			if (vecExprContents[t].type != FieldUnion<double>::Type::Field) {
				continue;
			}
			vecExprContents[t].fielddata.Allocate(lGridSize);
		}

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
					sprintf(szPos, "%s=%li", varIn->get_dim(d)->name(), lPos[d]);
					strPos += szPos;
					if (d != vecAuxDimSize.size()-1) {
						strPos += ",";
					}
				}
				Announce("Processing %s(%s)",
					vecVariableStrings[v].c_str(),
					strPos.c_str());
			}

			// Load field data
			std::vector<long> lPosF = lPos;
			std::vector<long> lSizeF = lSize;
			lPosF.erase(lPosF.begin() + lIntegralDimIx);
			lSizeF.erase(lSizeF.begin() + lIntegralDimIx);

			for (int t = 0; t < vecExprContents.size(); t++) {
				if (vecExprContents[t].type != FieldUnion<double>::Type::Field) {
					continue;
				}

				_ASSERT(lPosF.size() == vecExprContents[t].fieldvar->num_dims());
				_ASSERT(lSizeF.size() == vecExprContents[t].fieldvar->num_dims());

				vecExprContents[t].fieldvar->set_cur(&(lPosF[0]));
				vecExprContents[t].fieldvar->get(&(vecExprContents[t].fielddata[0]), &(lSizeF[0]));
			}

			// Loop through integral dimension
			_ASSERT(lPos.size() > lIntegralDimIx);
			_ASSERT(varIn->get_dim(lIntegralDimIx)->size() == lIntegralDimSize);

			// "a*p0+b*ps" style hybrid expression
			if ((exprHybridExpr.size() == 7) &&
				(vecExprContents[0].type == FieldUnion<double>::Type::BoundsVector) &&
				(exprHybridExpr[1].str == "*") &&
				(vecExprContents[2].type == FieldUnion<double>::Type::Scalar) &&
				(exprHybridExpr[3].str == "+") &&
				(vecExprContents[4].type == FieldUnion<double>::Type::BoundsVector) &&
				(exprHybridExpr[5].str == "*") &&
				(vecExprContents[6].type == FieldUnion<double>::Type::Field)
			) {
				double dRefValue = vecExprContents[2].scalar;

				for (long lLev = 0; lLev < lIntegralDimSize; lLev++) {

					lPos[lIntegralDimIx] = lLev;
					varIn->set_cur(&(lPos[0]));
					varIn->get(&(dDataVar[0]), &(lSize[0]));

					for (long lGrid = 0; lGrid < lGridSize; lGrid++) {
						double dPl =
							vecExprContents[0].bounds(lLev,0) * dRefValue
							+ vecExprContents[4].bounds(lLev,0) * vecExprContents[6].fielddata(lGrid);

						double dPu =
							vecExprContents[0].bounds(lLev,1) * dRefValue
							+ vecExprContents[4].bounds(lLev,1) * vecExprContents[6].fielddata(lGrid);

						dDataOut[lGrid] += fabs(dPu - dPl) * dDataVar[lGrid];
					}
				}

			// "ap+b*ps" style hybrid expression
			} else if ((exprHybridExpr.size() == 5) &&
				(vecExprContents[0].type == FieldUnion<double>::Type::BoundsVector) &&
				(exprHybridExpr[1].str == "+") &&
				(vecExprContents[2].type == FieldUnion<double>::Type::BoundsVector) &&
				(exprHybridExpr[3].str == "*") &&
				(vecExprContents[4].type == FieldUnion<double>::Type::Field)
			) {

				for (long lLev = 0; lLev < lIntegralDimSize; lLev++) {

					lPos[lIntegralDimIx] = lLev;
					varIn->set_cur(&(lPos[0]));
					varIn->get(&(dDataVar[0]), &(lSize[0]));

					for (long lGrid = 0; lGrid < lGridSize; lGrid++) {
						double dPl =
							vecExprContents[0].bounds(lLev,0)
							+ vecExprContents[2].bounds(lLev,0) * vecExprContents[4].fielddata(lGrid);

						double dPu =
							vecExprContents[0].bounds(lLev,1)
							+ vecExprContents[2].bounds(lLev,1) * vecExprContents[4].fielddata(lGrid);

						//if (lGrid == 0) {
						//	printf("%1.8f %1.8f : %1.8f\n", dPu, dPl, vecExprContents[4].fielddata(lGrid));
						//}

						dDataOut[lGrid] += fabs(dPu - dPl) * dDataVar[lGrid];
					}
				}

			// Unknown expression
			} else {
				_EXCEPTION1("Unimplemented --hybridexpr \"%s\"", strHybridExpr.c_str());
			}

			// Rescale by -1/g
			if (strHybridCoordType == "p") {
				for (long lGrid = 0; lGrid < lGridSize; lGrid++) {
					dDataOut[lGrid] *= 1.0 / EarthGravity;
				}
			}
/*
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
*/
			// Write output data
			varOut->set_cur(&(lPosF[0]));
			varOut->put(&(dDataOut[0]), &(lSizeF[0]));
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

