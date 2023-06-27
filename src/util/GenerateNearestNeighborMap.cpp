///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateNearestNeighborMap.cpp
///	\author  Paul Ullrich
///	\version March 30, 2023
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
#include "CoordTransforms.h"

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

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Source data
	std::string strSourceData;

	// Source connectivity file
	std::string strSourceConnect;

	// Source variable (used for masking)
	std::string strSourceVar;

	// Target data
	std::string strTargetData;

	// Target connectivity file
	std::string strTargetConnect;

	// Target variable (used for masking)
	std::string strTargetVar;

	// Regional
	bool fRegional = true;

	// Output map
	std::string strOutputMap;

	// Name of latitude dimension on source grid
	std::string strSourceLatitudeName;

	// Name of longitude dimension on source grid
	std::string strSourceLongitudeName;

	// Name of latitude dimension on target grid
	std::string strTargetLatitudeName;

	// Name of longitude dimension on target grid
	std::string strTargetLongitudeName;

	// Maximum allowed distance from nearest neighbor
	double dMaximumDistDeg;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strSourceData, "src_data", "");
		CommandLineString(strSourceConnect, "src_connect", "");
		CommandLineString(strSourceVar, "src_var", "");
		CommandLineString(strTargetData, "tgt_data", "");
		CommandLineString(strTargetConnect, "tgt_connect", "");
		CommandLineString(strTargetVar, "tgt_var", "");
		CommandLineString(strOutputMap, "out_map", "");
		CommandLineString(strSourceLongitudeName, "src_lonname", "lon");
		CommandLineString(strSourceLatitudeName, "src_latname", "lat");
		CommandLineString(strTargetLongitudeName, "tgt_lonname", "lon");
		CommandLineString(strTargetLatitudeName, "tgt_latname", "lat");
		CommandLineDoubleD(dMaximumDistDeg, "maxdist", 180.0, "(degrees GCD)");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if ((strSourceData.length() == 0) && (strSourceConnect.length() == 0)) {
		_EXCEPTIONT("No source data file (--src_data) or source connectivity file (--src_connect) specified");
	}
	if ((strSourceData.length() != 0) && (strSourceConnect.length() != 0)) {
		_EXCEPTIONT("Only one source data file (--src_data) or source connectivity file (--src_connect) may be specified");
	}
	if ((strSourceData.length() == 0) && (strSourceVar.length() != 0)) {
		_EXCEPTIONT("Argument (--src_var) must be combined with (--src_data)");
	}
	if ((strTargetData.length() == 0) && (strTargetConnect.length() == 0)) {
		_EXCEPTIONT("No target data file (--tgt_data) or target connectivity file (--tgt_connect) specified");
	}
	if ((strTargetData.length() != 0) && (strTargetConnect.length() != 0)) {
		_EXCEPTIONT("Only one target data file (--tgt_data) or target connectivity file (--tgt_connect) may be specified");
	}
	if ((strTargetData.length() == 0) && (strTargetVar.length() != 0)) {
		_EXCEPTIONT("Argument (--tgt_var) must be combined with (--tgt_data)");
	}
	if (strOutputMap.length() == 0) {
		_EXCEPTIONT("No output map (--out_map) specified");
	}

	// Define the source SimpleGrid
	SimpleGrid gridSource;

	// Check for connectivity file
	if (strSourceConnect != "") {
		AnnounceStartBlock("Generating source grid information from connectivity file");
		gridSource.FromFile(strSourceConnect);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		Announce("Generating latitude-longitude grid from source data file");

		// Load in file vector
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(strSourceData);
		_ASSERT(vecNcFiles.size() > 0);

		gridSource.GenerateLatitudeLongitude(
			vecNcFiles[0],
			strSourceLatitudeName,
			strSourceLongitudeName,
			fRegional,
			true);

		if (gridSource.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Define the target SimpleGrid
	SimpleGrid gridTarget;

	// Check for connectivity file
	if (strTargetConnect != "") {
		AnnounceStartBlock("Generating source grid information from connectivity file");
		gridTarget.FromFile(strTargetConnect);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		Announce("Generating latitude-longitude grid from target data file");

		// Load in file vector
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(strTargetData);
		_ASSERT(vecNcFiles.size() > 0);

		gridTarget.GenerateLatitudeLongitude(
			vecNcFiles[0],
			strTargetLatitudeName,
			strTargetLongitudeName,
			fRegional,
			true);

		if (gridTarget.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Build the KD tree on the source grid
	if (strSourceVar == "") {
		AnnounceStartBlock("Building kd-tree on source grid");
		gridSource.BuildKDTree();
		AnnounceEndBlock("Done");

	} else {
		AnnounceStartBlock("Building kd-tree on source grid (masked)");
		size_t sGridSizeSrc = gridSource.m_dLon.GetRows();

		// Load the data for the src grid
		NcFileVector vecNcFilesSrc;
		vecNcFilesSrc.ParseFromString(strSourceData);
		vecNcFilesSrc.SetConstantTimeIx(0);

		VariableRegistry varregSrc;
		VariableIndex ixSrcData = varregSrc.FindOrRegister(strSourceVar);
		Variable & varSrcData = varregSrc.Get(ixSrcData);
		varSrcData.LoadGridData(varregSrc, vecNcFilesSrc, gridSource);

		const DataArray1D<float> & dataSrc = varSrcData.GetData();

		float dFillValueSrc = varSrcData.GetFillValueFloat();

		DataArray1D<bool> fMask(sGridSizeSrc);
		for (size_t i = 0; i < sGridSizeSrc; i++) {
			if ((dataSrc[i] == dFillValueSrc) || (dataSrc[i] != dataSrc[i])) {
				fMask[i] = false;
			} else {
				fMask[i] = true;
			}
		}
		gridSource.BuildMaskedKDTree(fMask);
		AnnounceEndBlock("Done");
	}

	// Load target data
	VariableRegistry varregTgt;
	float dFillValueTgt;
	const DataArray1D<float> * pTargetData = NULL;

	if (strTargetVar != "") {
		size_t sGridSizeTgt = gridTarget.m_dLon.GetRows();

		// Load the data for the tgt grid
		NcFileVector vecNcFilesTgt;
		vecNcFilesTgt.ParseFromString(strTargetData);
		vecNcFilesTgt.SetConstantTimeIx(0);

		VariableIndex ixTgtData = varregTgt.FindOrRegister(strTargetVar);
		Variable & varTgtData = varregTgt.Get(ixTgtData);
		varTgtData.LoadGridData(varregTgt, vecNcFilesTgt, gridTarget);

		pTargetData = &(varTgtData.GetData());
		dFillValueTgt = varTgtData.GetFillValueFloat();

		_ASSERT(pTargetData->GetRows() == gridTarget.GetSize());
	}

	// Build the map
	AnnounceStartBlock("Building map");
	long nA = gridSource.GetSize();
	long nB = gridTarget.GetSize();

	std::vector<int> vecCol(nB, 1);
	std::vector<int> vecRow(nB, 1);
	std::vector<double> vecS(nB, 0.0);

	std::vector<double> vecSrcLonDeg(nA);
	std::vector<double> vecSrcLatDeg(nA);

	std::vector<double> vecTgtLonDeg(nB);
	std::vector<double> vecTgtLatDeg(nB);

	for (int i = 0; i < nA; i++) {
		vecSrcLonDeg[i] = RadToDeg(gridSource.m_dLon[i]);
		vecSrcLatDeg[i] = RadToDeg(gridSource.m_dLat[i]);
	}
	for (int j = 0; j < nB; j++) {
		if ((j % (nB / 10) == 0) && (j != 0)) {
			int nPctComplete = j / (nB / 10);
			Announce("%i%% complete", nPctComplete * 10);
		}

		vecTgtLonDeg[j] = RadToDeg(gridTarget.m_dLon[j]);
		vecTgtLatDeg[j] = RadToDeg(gridTarget.m_dLat[j]);

		if (pTargetData != NULL) {
			if (((*pTargetData)[j] == dFillValueTgt) || ((*pTargetData)[j] != (*pTargetData)[j])) {
				continue;
			}
		}

		size_t sNearestNode = 
			gridSource.NearestNode(
				gridTarget.m_dLon[j],
				gridTarget.m_dLat[j]);

		if (dMaximumDistDeg < 180.0) {
			double dR =
				GreatCircleDistance_Deg(
					DegToRad(vecSrcLonDeg[sNearestNode]),
					DegToRad(vecSrcLatDeg[sNearestNode]),
					DegToRad(vecTgtLonDeg[j]),
					DegToRad(vecTgtLatDeg[j]));

			if (dR > dMaximumDistDeg) {
				vecRow[j] = 1;
				vecCol[j] = 1;
				vecS[j] = 0.0;
				continue;
			}
		}

		_ASSERT(sNearestNode < static_cast<size_t>(nA));

		vecRow[j] = j + 1;
		vecCol[j] = static_cast<int>(sNearestNode) + 1;
		vecS[j] = 1.0;
	}
	AnnounceEndBlock(NULL);

	// Write the map
	AnnounceStartBlock("Writing map");
	NcFile ncmap(strOutputMap.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
	if (!ncmap.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\"", strOutputMap.c_str());
	}

	ncmap.add_att("Title", "TempestExtremes GenerateNearestNeighborMap");

	// Write output dimensions entries
	int nSrcGridDims = (int)(gridSource.m_nGridDim.size());
	int nDstGridDims = (int)(gridTarget.m_nGridDim.size());

	NcDim * dimSrcGridRank = ncmap.add_dim("src_grid_rank", gridSource.m_nGridDim.size());
	NcDim * dimDstGridRank = ncmap.add_dim("dst_grid_rank", gridTarget.m_nGridDim.size());

	NcVar * varSrcGridDims = ncmap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
	NcVar * varDstGridDims = ncmap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

	if ((nSrcGridDims == 1) && (gridSource.m_nGridDim[0] != nA)) {
		varSrcGridDims->put(&nA, 1);
		varSrcGridDims->add_att("name0", "ncol");

	} else {
		for (int i = 0; i < gridSource.m_nGridDim.size(); i++) {
			int iGridDim = static_cast<int>(gridSource.m_nGridDim[i]);
			varSrcGridDims->set_cur(nSrcGridDims - i - 1);
			varSrcGridDims->put(&iGridDim, 1);
		}

		for (int i = 0; i < gridSource.m_nGridDim.size(); i++) {
			varSrcGridDims->add_att("name0", strSourceLongitudeName.c_str());
			varSrcGridDims->add_att("name1", strSourceLatitudeName.c_str());
		}
	}

	if ((nDstGridDims == 1) && (gridTarget.m_nGridDim[0] != nA)) {
		varDstGridDims->put(&nA, 1);
		varDstGridDims->add_att("name0", "ncol");

	} else {
		for (int i = 0; i < gridTarget.m_nGridDim.size(); i++) {
			int iGridDim = static_cast<int>(gridTarget.m_nGridDim[i]);
			varDstGridDims->set_cur(nDstGridDims - i - 1);
			varDstGridDims->put(&iGridDim, 1);
		}

		for (int i = 0; i < gridTarget.m_nGridDim.size(); i++) {
			varDstGridDims->add_att("name0", strTargetLongitudeName.c_str());
			varDstGridDims->add_att("name1", strTargetLatitudeName.c_str());
		}
	}

	// Source and Target mesh resolutions
	NcDim * dimNA = ncmap.add_dim("n_a", nA);
	NcDim * dimNB = ncmap.add_dim("n_b", nB);

	NcDim * dimNVA = ncmap.add_dim("nv_a", 1);
	NcDim * dimNVB = ncmap.add_dim("nv_b", 1);

	// Coordinates
	NcVar * varYCA = ncmap.add_var("yc_a", ncDouble, dimNA);
	NcVar * varYCB = ncmap.add_var("yc_b", ncDouble, dimNB);

	NcVar * varXCA = ncmap.add_var("xc_a", ncDouble, dimNA);
	NcVar * varXCB = ncmap.add_var("xc_b", ncDouble, dimNB);

	NcVar * varYVA = ncmap.add_var("yv_a", ncDouble, dimNA, dimNVA);
	NcVar * varYVB = ncmap.add_var("yv_b", ncDouble, dimNB, dimNVB);

	NcVar * varXVA = ncmap.add_var("xv_a", ncDouble, dimNA, dimNVA);
	NcVar * varXVB = ncmap.add_var("xv_b", ncDouble, dimNB, dimNVB);

	std::vector<double> dAreaA(nA, 1.0);
	std::vector<double> dAreaB(nA, 1.0);

	NcVar * varAreaA = ncmap.add_var("area_a", ncDouble, dimNA);
	NcVar * varAreaB = ncmap.add_var("area_b", ncDouble, dimNB);

	varAreaA->put(&(dAreaA[0]), nA);
	varAreaB->put(&(dAreaB[0]), nB);

	varAreaA->add_att("units", "unitless");
	varAreaB->add_att("units", "unitless");

	std::vector<double> dFracA(nA, 1.0);
	std::vector<double> dFracB(nA, 1.0);

	NcVar * varFracA = ncmap.add_var("frac_a", ncDouble, dimNA);
	NcVar * varFracB = ncmap.add_var("frac_b", ncDouble, dimNB);

	varFracA->add_att("name", "fraction of target coverage of source dof");
	varFracB->add_att("name", "fraction of target coverage of target dof");

	varFracA->add_att("units", "unitless");
	varFracB->add_att("units", "unitless");

	varFracA->put(&(dFracA[0]), nA);
	varFracB->put(&(dFracB[0]), nB);

	varYCA->add_att("units", "degrees");
	varYCB->add_att("units", "degrees");

	varXCA->add_att("units", "degrees");
	varXCB->add_att("units", "degrees");

	varYVA->add_att("units", "degrees");
	varYVB->add_att("units", "degrees");

	varXVA->add_att("units", "degrees");
	varXVB->add_att("units", "degrees");

	varYCA->put(&(vecSrcLatDeg[0]), vecSrcLatDeg.size());
	varYCB->put(&(vecTgtLatDeg[0]), vecTgtLatDeg.size());

	varXCA->put(&(vecSrcLonDeg[0]), vecSrcLonDeg.size());
	varXCB->put(&(vecTgtLonDeg[0]), vecTgtLonDeg.size());

	varYVA->put(&(vecSrcLatDeg[0]), vecSrcLatDeg.size(), 1);
	varYVB->put(&(vecTgtLatDeg[0]), vecTgtLatDeg.size(), 1);

	varXVA->put(&(vecSrcLonDeg[0]), vecSrcLonDeg.size(), 1);
	varXVB->put(&(vecTgtLonDeg[0]), vecTgtLonDeg.size(), 1);

	// 1D longitude and latitude
	if ((gridTarget.m_nGridDim.size() != 1) && (gridTarget.m_nGridDim[0] != gridTarget.GetSize())) {
		NcDim * dimBnds = ncmap.add_dim("bnds", 2);

		NcDim * dimLatB = ncmap.add_dim("lat_b", gridTarget.m_nGridDim[0]);
		NcDim * dimLonB = ncmap.add_dim("lon_b", gridTarget.m_nGridDim[1]);

		NcVar * varLatBnds = ncmap.add_var("lat_bnds", ncDouble, dimLatB, dimBnds);
		NcVar * varLonBnds = ncmap.add_var("lon_bnds", ncDouble, dimLonB, dimBnds);

		NcVar * varLatB = ncmap.add_var("latc_b", ncDouble, dimLatB);
		NcVar * varLonB = ncmap.add_var("lonc_b", ncDouble, dimLonB);

		std::vector<double> dLat1D(gridTarget.m_nGridDim[0]);
		for (int i = 0; i < gridTarget.m_nGridDim[0]; i++) {
			dLat1D[i] = RadToDeg(gridTarget.m_dLat[i * gridTarget.m_nGridDim[1]]);
		}
		std::vector<double> dLon1D(gridTarget.m_nGridDim[1]);
		for (int i = 0; i < gridTarget.m_nGridDim[1]; i++) {
			dLon1D[i] = RadToDeg(gridTarget.m_dLon[i]);
		}

		varLatB->put(&(dLat1D[0]), dLat1D.size());
		varLonB->put(&(dLon1D[0]), dLon1D.size());
	}

	// Write out data
	int nS = vecS.size();
	NcDim * dimNS = ncmap.add_dim("n_s", nS);

	NcVar * varRow = ncmap.add_var("row", ncInt, dimNS);
	varRow->add_att("name", "sparse matrix target dof index");
	varRow->add_att("first_index", "1");

	NcVar * varCol = ncmap.add_var("col", ncInt, dimNS);
	varCol->add_att("name", "sparse matrix source dof index");
	varCol->add_att("first_index", "1");

	NcVar * varS = ncmap.add_var("S", ncDouble, dimNS);
	varS->add_att("name", "sparse matrix coefficient");

	varRow->set_cur((long)0);
	varRow->put(&(vecRow[0]), nS);

	varCol->set_cur((long)0);
	varCol->put(&(vecCol[0]), nS);

	varS->set_cur((long)0);
	varS->put(&(vecS[0]), nS);

	AnnounceStartBlock("Done");
	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

}
