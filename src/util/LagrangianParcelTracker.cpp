///////////////////////////////////////////////////////////////////////////////
///
///	\file    LagrangianParcelTracker.cpp
///	\author  Paul Ullrich
///	\version April 11, 2021
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

#include "Constants.h"
#include "CoordTransforms.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "NetCDFUtilities.h"
#include "DataArray1D.h"
#include "Variable.h"
#include "NodeFileUtilities.h"
#include "FilenameList.h"
#include "GridElements.h"

#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

struct LLPPosition {
	int ix;
	double lon;
	double lat;
	double pres;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate remapping weights for a particular location using inverse
///		distance to nearby grid points.
///	</summary>
void CalculateRemapWeightsByWeightedDistance(
	const SimpleGrid & grid,
	double dLonRad,
	double dLatRad,
	std::map<size_t, double> & mapDist
) {

	// Find nearest neighbor and calculate distance
	size_t sNearestIx = grid.NearestNode(dLonRad, dLatRad);

	double dGCD0Rad =
		GreatCircleDistance_Rad(
			dLonRad,
			dLatRad,
			grid.m_dLon[sNearestIx],
			grid.m_dLat[sNearestIx]);

	if (dGCD0Rad < 1.0e-12) {
		mapDist.insert(std::pair<size_t, double>(sNearestIx, 1.0));
		return;
	}

	mapDist.insert(std::pair<size_t, double>(sNearestIx, 1.0/dGCD0Rad));

	// Find distance to grid cells connected to this point
	const std::vector<int> & vecConnectivity = grid.m_vecConnectivity[sNearestIx];
	for (size_t i = 0; i < vecConnectivity.size(); i++) {

		size_t sNeighborIx = vecConnectivity[i];

		double dGCD1Rad =
			GreatCircleDistance_Rad(
				grid.m_dLon[sNearestIx],
				grid.m_dLat[sNearestIx],
				grid.m_dLon[sNeighborIx],
				grid.m_dLat[sNeighborIx]);

		double dGCD2Rad =
			GreatCircleDistance_Rad(
				dLonRad,
				dLatRad,
				grid.m_dLon[sNeighborIx],
				grid.m_dLat[sNeighborIx]);

		if (dGCD2Rad >= dGCD1Rad) {
			continue;
		}

		mapDist.insert(std::pair<size_t, double>(sNeighborIx, 1.0/dGCD2Rad));
	}

	// Reweight by total distance to enforce consistency
	double dTotalDist = 0.0;
	for (auto itDist = mapDist.begin(); itDist != mapDist.end(); itDist++) {
		dTotalDist += itDist->second;
	}
	for (auto itDist = mapDist.begin(); itDist != mapDist.end(); itDist++) {
		itDist->second /= dTotalDist;
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		At nodeX0 generate an orthonormal basis tangent to the sphere and
///		store in vectors nodeX1 and nodeX2.
///	</summary>
void GenerateOrthonormalBasis(
	const Node & nodeX0,
	Node & nodeX1,
	Node & nodeX2
) {
	if (nodeX0.z > 0.9) {
		nodeX1.Set(1.0, 0.0, 0.0);
		nodeX1 -= (nodeX0 * DotProduct(nodeX0, nodeX1));
		nodeX1 /= nodeX1.Magnitude();

		nodeX2 = CrossProduct(nodeX0, nodeX1);

	} else {
		nodeX1.Set(0.0, 0.0, 1.0);
		nodeX1 -= (nodeX0 * DotProduct(nodeX0, nodeX1));
		nodeX1 /= nodeX1.Magnitude();

		nodeX2 = CrossProduct(nodeX0, nodeX1);
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Project the velocity (dUlon,dUlat) at point (dLonRad, dLatRad) onto
///		the orthogonal basis defined at nodeX0 consisting of vectors
///		(nodeX1, nodeX2) and store in (dU1,dU2).
///	</summary>
void ProjectVelocityOntoOrthonormalFields(
	const Node & nodeX0,
	const Node & nodeX1,
	const Node & nodeX2,
	const double dLonRad,
	const double dLatRad,
	const double dUlon,
	const double dUlat,
	double & dU1,
	double & dU2
) {
	// Cartesian velocity at point
	Node nodeU;
	VecTransRLL2DtoXYZ_Rad(
		dLonRad, dLatRad,
		dUlon, dUlat,
		nodeU.x, nodeU.y, nodeU.z);

	// nodeY1 and nodeY2 are the projections of this point onto the (X0,X1) and (X0,X2) planes
	Node nodeY;
	RLLtoXYZ_Rad(dLonRad, dLatRad, nodeY.x, nodeY.y, nodeY.z);

	Node nodeY1 = nodeY;
	Node nodeY2 = nodeY;

	_ASSERT(fabs(DotProduct(nodeU, nodeY1)) < 1.0e-12);

	nodeY1 -= nodeX2 * DotProduct(nodeX2, nodeY1);
	nodeY2 -= nodeX1 * DotProduct(nodeX1, nodeY1);

	// Calculate angle along the great circle in X1 direction
	double dA1 = DotProduct(nodeY1, nodeX0);
	double dB1 = DotProduct(nodeY1, nodeX1);
	_ASSERT((fabs(dA1) > 1.0e-12) || (fabs(dB1) > 1.0e-12));

	double dSinTheta1 = dB1 / sqrt(dB1 * dB1 + dA1 * dA1);
	double dCosTheta1 = dA1 / sqrt(dB1 * dB1 + dA1 * dA1);

	//printf("%1.15e %1.15e : %1.15e %1.15e\n", dSinTheta1, sin(dTheta1), dCosTheta1, cos(dTheta1));

	// Calculate local tangent vector in X1 direction
	Node nodeT1;
	nodeT1.x = - dSinTheta1 * nodeX0.x + dCosTheta1 * nodeX1.x;
	nodeT1.y = - dSinTheta1 * nodeX0.y + dCosTheta1 * nodeX1.y;
	nodeT1.z = - dSinTheta1 * nodeX0.z + dCosTheta1 * nodeX1.z;

	Node nodeC = CrossProduct(nodeY, nodeT1);

	_ASSERT(fabs(DotProduct(nodeY, nodeT1)) < 1.0e-12);
	_ASSERT(fabs(DotProduct(nodeT1, nodeT1) - 1.0) < 1.0e-12);
	_ASSERT(fabs(DotProduct(nodeC, nodeC) - 1.0) < 1.0e-12);

	// Calculate X1 contribution to velocity
	dU1 = DotProduct(nodeU, nodeT1);

	// Repeat calculation in X2 direction
	double dA2 = DotProduct(nodeY2, nodeX0);
	double dB2 = DotProduct(nodeY2, nodeX2);
	_ASSERT((fabs(dA2) > 1.0e-12) || (fabs(dB2) > 1.0e-12));

	double dSinTheta2 = dB2 / sqrt(dB2 * dB2 + dA2 * dA2);
	double dCosTheta2 = dA2 / sqrt(dB2 * dB2 + dA2 * dA2);

	Node nodeT2;
	nodeT2.x = - dSinTheta2 * nodeX0.x + dCosTheta2 * nodeX2.x;
	nodeT2.y = - dSinTheta2 * nodeX0.y + dCosTheta2 * nodeX2.y;
	nodeT2.z = - dSinTheta2 * nodeX0.z + dCosTheta2 * nodeX2.z;

	_ASSERT(fabs(DotProduct(nodeY, nodeT2)) < 1.0e-12);
	_ASSERT(fabs(DotProduct(nodeT2, nodeT2) - 1.0) < 1.0e-12);

	dU2 = DotProduct(nodeU, nodeT2);

	//nodeC.Print("C ");
	//nodeT2.Print("T2");
}

///////////////////////////////////////////////////////////////////////////////

void SwitchFiles(
	NcFileVector & vecNcFiles,
	const std::string & strInputFiles,
	const std::string & strUVariable,
	const std::string & strVVariable,
	const std::string & strWVariable,
	DataArray1D<double> & dLevel,
	NcVar ** pvarU,
	NcVar ** pvarV,
	NcVar ** pvarW
) {
	vecNcFiles.clear();
	vecNcFiles.ParseFromString(strInputFiles);
	_ASSERT(vecNcFiles.size() > 0);

	// Load the level variable
	NcVar * varLevel = vecNcFiles[0]->get_var("level");
	if (varLevel == NULL) {
		_EXCEPTIONT("Unable to load variable \"level\"");
	}
	if (varLevel->num_dims() != 1) {
		_EXCEPTIONT("Variable \"level\" has more than one dimension");
	}
	dLevel.Allocate(varLevel->get_dim(0)->size());
	varLevel->get(&(dLevel[0]), varLevel->get_dim(0)->size());

	// TODO: Verify level strict monotonicity

	// Convert units
	NcAtt * attLevUnits = varLevel->get_att("units");
	if (attLevUnits != NULL) {
		std::string strLevelUnits = attLevUnits->as_string(0);
		if ((strLevelUnits == "hPa") || (strLevelUnits == "mb")) {
			for (size_t k = 0; k < dLevel.GetRows(); k++) {
				dLevel[k] *= 100.0;
			}
		}
	}

	// Find the file with the right variables
	vecNcFiles.FindContainingVariable(strUVariable, pvarU);
	if ((*pvarU) == NULL) {
		_EXCEPTION2("Unable to find variable \"%s\" in \"%s\"",
			strUVariable.c_str(), strInputFiles.c_str());
	}
	vecNcFiles.FindContainingVariable(strVVariable, pvarV);
	if ((*pvarV) == NULL) {
		_EXCEPTION2("Unable to find variable \"%s\" in \"%s\"",
			strVVariable.c_str(), strInputFiles.c_str());
	}
	vecNcFiles.FindContainingVariable(strWVariable, pvarW);
	if ((*pvarW) == NULL) {
		_EXCEPTION2("Unable to find variable \"%s\" in \"%s\"",
			strWVariable.c_str(), strInputFiles.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load all grid data needed for Lagrangian parcel tracking.
///	</summary>
void LoadGridData(
	const SimpleGrid & grid,
	const long f,
	const long t,
	const DataArray1D<double> & dLevel,
	const std::set<int> & setActivePathIxs,
	const std::vector<LLPPosition> & vecX,
	std::map<int, std::pair<int,int> > & mapPathIxToBoundingLevelIx,
	std::map<int,int> & mapLevelIxToLoadIx,
	NcVar * varU,
	NcVar * varV,
	NcVar * varW,
	std::vector<DataArray1D<float> *> & vecdataU,
	std::vector<DataArray1D<float> *> & vecdataV,
	std::vector<DataArray1D<float> *> & vecdataW,
	bool fKeepLoadedData
) {
	if (!fKeepLoadedData) {
		mapPathIxToBoundingLevelIx.clear();
		mapLevelIxToLoadIx.clear();
	}

	// Determine monotonic order of levels
	double dLevelSign;
	if (dLevel.GetRows() <= 1) {
		dLevelSign = 1.0;
	} else if (dLevel[1] > dLevel[0]) {
		dLevelSign = 1.0;
	} else {
		dLevelSign = -1.0;
	}

	// Determine which levels are needed
	std::set<int> setNeededLevels;
	for (auto itPath = setActivePathIxs.begin(); itPath != setActivePathIxs.end(); itPath++) {
		if (vecX[*itPath].pres * dLevelSign <= dLevel[0] * dLevelSign) {
			setNeededLevels.insert(0);
			mapPathIxToBoundingLevelIx.insert(
				std::map<int, std::pair<int,int> >::value_type(
					*itPath, std::pair<int,int>(-1,0)));
			continue;
		}
		if (vecX[*itPath].pres * dLevelSign >= dLevel[dLevel.GetRows()-1] * dLevelSign) {
			setNeededLevels.insert(dLevel.GetRows()-1);
			mapPathIxToBoundingLevelIx.insert(
				std::map<int, std::pair<int,int> >::value_type(
					*itPath, std::pair<int,int>(dLevel.GetRows()-1,-1)));
			continue;
		}
		for (size_t lev = 0; lev < dLevel.GetRows()-1; lev++) {
			if ((vecX[*itPath].pres * dLevelSign >= dLevel[lev] * dLevelSign) &&
			    (vecX[*itPath].pres * dLevelSign <= dLevel[lev+1] * dLevelSign)
			) {
				setNeededLevels.insert(lev);
				setNeededLevels.insert(lev+1);
				mapPathIxToBoundingLevelIx.insert(
					std::map<int, std::pair<int,int> >::value_type(
						*itPath, std::pair<int,int>(lev,lev+1)));
			}
		}
	}

	// Load velocity data at this time slice
	if (setActivePathIxs.size() != 0) {
		std::vector<long> nSize(2+grid.m_nGridDim.size(), 1);
		std::vector<long> nPos(2+grid.m_nGridDim.size(), 0);
		nPos[0] = t;
		for (int d = 0; d < grid.m_nGridDim.size(); d++) {
			nSize[d+2] = static_cast<long>(grid.m_nGridDim[d]);
		}

		_ASSERT(nSize.size() == nPos.size());
		_ASSERT(nSize.size() == varU->num_dims());
		_ASSERT(nSize.size() == varV->num_dims());
		_ASSERT(nSize.size() == varW->num_dims());

		int iLoadIx = mapLevelIxToLoadIx.size();
		for (auto itLev = setNeededLevels.begin(); itLev != setNeededLevels.end(); itLev++) {

			// Don't load data if it's already loaded
			if (mapLevelIxToLoadIx.find(*itLev) != mapLevelIxToLoadIx.end()) {
				continue;
			}

			Announce("Loading velocity (file %i, time %i, level %i) into data index %i", f, t, *itLev, iLoadIx);

			if (iLoadIx >= vecdataU.size()) {
				DataArray1D<float> * pdataU = new DataArray1D<float>(grid.GetSize());
				DataArray1D<float> * pdataV = new DataArray1D<float>(grid.GetSize());
				DataArray1D<float> * pdataW = new DataArray1D<float>(grid.GetSize());

				if ((pdataU == NULL) || (pdataV == NULL) || (pdataW == NULL)) {
					_EXCEPTIONT("Out of memory");
				}

				vecdataU.push_back(pdataU);
				vecdataV.push_back(pdataV);
				vecdataW.push_back(pdataW);
			}

			nPos[1] = *itLev;

			DataArray1D<float> & dataU = *(vecdataU[iLoadIx]);
			DataArray1D<float> & dataV = *(vecdataV[iLoadIx]);
			DataArray1D<float> & dataW = *(vecdataW[iLoadIx]);

			varU->set_cur(&(nPos[0]));
			varU->get(&(dataU[0]), &(nSize[0]));

			varV->set_cur(&(nPos[0]));
			varV->get(&(dataV[0]), &(nSize[0]));

			varW->set_cur(&(nPos[0]));
			varW->get(&(dataW[0]), &(nSize[0]));
/*
			// Solid body rotation with 1 day cycle
			for (int i = 0; i < grid.GetSize(); i++) {
				dataU[i] = 2.0 * M_PI * EarthRadius * cos(grid.m_dLat[i]) / 86400.0;
				dataV[i] = 0.0;
				dataW[i] = 0.0;
			}
*/
			mapLevelIxToLoadIx.insert(std::pair<int,int>(*itLev, iLoadIx));

			iLoadIx++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Update the path position using equation [Xs = Xs - dDeltaTime * u(Xn)]
///	</summary>
void UpdatePaths(
	const SimpleGrid & grid,
	const DataArray1D<double> & dLevel,
	const std::set<int> & setActivePathIxs,
	const std::map<int, std::pair<int,int> > & mapPathIxToBoundingLevelIx,
	const std::map<int,int> & mapLevelIxToLoadIx,
	const std::vector<DataArray1D<float> *> & vecdataU,
	const std::vector<DataArray1D<float> *> & vecdataV,
	const std::vector<DataArray1D<float> *> & vecdataW,
	double dDeltaTime,
	const std::vector<LLPPosition> & vecXn,
	std::vector<LLPPosition> & vecXs
) {
	for (auto itPath = setActivePathIxs.begin(); itPath != setActivePathIxs.end(); itPath++) {

		// Find nearest neighbors
		std::map<size_t, double> mapDist;
		CalculateRemapWeightsByWeightedDistance(
			grid,
			DegToRad(vecXn[*itPath].lon),
			DegToRad(vecXn[*itPath].lat),
			mapDist);

		// Build local orthonormal basis
		// nodeX0 is the local node
		// nodeX1 and nodeX2 are orthonormal tangent vectors to the sphere at this point
		Node nodeX0, nodeX1, nodeX2;
		RLLtoXYZ_Deg(vecXn[*itPath].lon, vecXn[*itPath].lat, nodeX0.x, nodeX0.y, nodeX0.z);

		GenerateOrthonormalBasis(nodeX0, nodeX1, nodeX2);

		_ASSERT(fabs(DotProduct(nodeX0, nodeX1)) < 1.0e-12);
		_ASSERT(fabs(DotProduct(nodeX0, nodeX2)) < 1.0e-12);
		_ASSERT(fabs(DotProduct(nodeX1, nodeX2)) < 1.0e-12);
		_ASSERT(fabs(DotProduct(nodeX0, nodeX0) - 1.0) < 1.0e-12);
		_ASSERT(fabs(DotProduct(nodeX1, nodeX1) - 1.0) < 1.0e-12);
		_ASSERT(fabs(DotProduct(nodeX2, nodeX2) - 1.0) < 1.0e-12);

		// Calculate velocity components
		double dTotalU1 = 0.0;
		double dTotalU2 = 0.0;
		double dTotalUp = 0.0;

		double dTotalUe = 0.0;
		double dTotalUn = 0.0;

		// Vertical weighting
		auto itBoundingLevels = mapPathIxToBoundingLevelIx.find(*itPath);
		_ASSERT(itBoundingLevels != mapPathIxToBoundingLevelIx.end());

		int iUpperLev = itBoundingLevels->second.second;
		int iLowerLev = itBoundingLevels->second.first;
		double dUpperWt;
		double dLowerWt;
		if (iUpperLev == (-1)) {
			dLowerWt = 1.0;
			dUpperWt = 0.0;

		} else if (iLowerLev == (-1)) {
			dLowerWt = 0.0;
			dUpperWt = 1.0;

		} else {
			double dPres = vecXn[*itPath].pres;
			dLowerWt = (dLevel[iUpperLev] - dPres) / (dLevel[iUpperLev] - dLevel[iLowerLev]);
			dUpperWt = (dPres - dLevel[iLowerLev]) / (dLevel[iUpperLev] - dLevel[iLowerLev]);
		}

		// Contributions to velocity from lower levels
		if (iLowerLev != (-1)) {
			auto itLevelIxToLoadIx = mapLevelIxToLoadIx.find(iLowerLev);
			_ASSERT(itLevelIxToLoadIx != mapLevelIxToLoadIx.end());
			int iLoadIx = itLevelIxToLoadIx->second;

			DataArray1D<float> * dataU = vecdataU[iLoadIx];
			DataArray1D<float> * dataV = vecdataV[iLoadIx];
			DataArray1D<float> * dataW = vecdataW[iLoadIx];

			_ASSERT(dataU != NULL);
			_ASSERT(dataV != NULL);
			_ASSERT(dataW != NULL);

			for (auto itDist = mapDist.begin(); itDist != mapDist.end(); itDist++) {

				double dU1;
				double dU2;

				ProjectVelocityOntoOrthonormalFields(
					nodeX0, nodeX1, nodeX2,
					grid.m_dLon(itDist->first),
					grid.m_dLat(itDist->first),
					(*dataU)[itDist->first],
					(*dataV)[itDist->first],
					dU1, dU2);

				//printf("U1 U2: %1.5f %1.5f\n", dU1, dU2);

				double dCombinedWeight = itDist->second * dLowerWt;
				dTotalU1 += dCombinedWeight * dU1;
				dTotalU2 += dCombinedWeight * dU2;
				dTotalUp += dCombinedWeight * (*dataW)[itDist->first];

				dTotalUe += dCombinedWeight * (*dataU)[itDist->first];
				dTotalUn += dCombinedWeight * (*dataV)[itDist->first];
			}
		}

		// Contributions to velocity from upper levels
		if (iUpperLev != (-1)) {
			auto itLevelIxToLoadIx = mapLevelIxToLoadIx.find(iUpperLev);
			_ASSERT(itLevelIxToLoadIx != mapLevelIxToLoadIx.end());
			int iLoadIx = itLevelIxToLoadIx->second;

			DataArray1D<float> * dataU = vecdataU[iLoadIx];
			DataArray1D<float> * dataV = vecdataV[iLoadIx];
			DataArray1D<float> * dataW = vecdataW[iLoadIx];

			_ASSERT(dataU != NULL);
			_ASSERT(dataV != NULL);
			_ASSERT(dataW != NULL);

			for (auto itDist = mapDist.begin(); itDist != mapDist.end(); itDist++) {

				double dU1;
				double dU2;

				ProjectVelocityOntoOrthonormalFields(
					nodeX0, nodeX1, nodeX2,
					grid.m_dLon(itDist->first),
					grid.m_dLat(itDist->first),
					(*dataU)[itDist->first],
					(*dataV)[itDist->first],
					dU1, dU2);

				double dCombinedWeight = itDist->second * dUpperWt;
				dTotalU1 += dCombinedWeight * dU1;
				dTotalU2 += dCombinedWeight * dU2;
				dTotalUp += dCombinedWeight * (*dataW)[itDist->first];

				dTotalUe += dCombinedWeight * (*dataU)[itDist->first];
				dTotalUn += dCombinedWeight * (*dataV)[itDist->first];
			}
		}

		//double dLonS0 = vecXs[*itPath].lon;
		//double dLatS0 = vecXs[*itPath].lat;

		//double dLonS = vecXs[*itPath].lon - RadToDeg(dDeltaTime * dTotalUe / (cos(DegToRad(vecXs[*itPath].lat)) * EarthRadius));
		//double dLatS = vecXs[*itPath].lat - RadToDeg(dDeltaTime * dTotalUn / EarthRadius);

		// Build local orthonormal basis
		// nodeX0 is the local node
		// nodeX1 and nodeX2 are orthonormal tangent vectors to the sphere at this point
		Node nodeU = (nodeX1 * dTotalU1) + (nodeX2 * dTotalU2);

		double dTangentDist = 0.0;
		double dUmag = nodeU.Magnitude();
		if (dUmag > 1.0e-12) {
			dTangentDist = 1.0 / dUmag * tan(dUmag * dDeltaTime / EarthRadius);
		}

		//nodeU.Print("U");
		//printf("Tangent Dist: %1.5f\n", dTangentDist * dUmag);

		Node nodeXs;
		RLLtoXYZ_Deg(vecXs[*itPath].lon, vecXs[*itPath].lat, nodeXs.x, nodeXs.y, nodeXs.z);

		nodeXs -= nodeU * dTangentDist;
		nodeXs /= nodeXs.Magnitude();

		XYZtoRLL_Deg(nodeXs.x, nodeXs.y, nodeXs.z, vecXs[*itPath].lon, vecXs[*itPath].lat);

		vecXs[*itPath].pres -= dDeltaTime * dTotalUp;
/*
		printf("lonlat: %1.5f %1.5f : %1.5f %1.5f : %1.5f %1.5f\n",
			dLonS0, dLatS0, dLonS, dLatS, vecXs[*itPath].lon, vecXs[*itPath].lat);

		printf("speed * time: %1.5f : %1.5f : %1.5f\n",
			dDeltaTime * sqrt(dTotalUe * dTotalUe + dTotalUn * dTotalUn),
			dDeltaTime * sqrt(dTotalU1 * dTotalU1 + dTotalU2 * dTotalU2),
			dDeltaTime * nodeU.Magnitude());

		printf("distance: %1.5f : %1.5f\n",
			EarthRadius * GreatCircleDistance_Rad(DegToRad(dLonS0), DegToRad(dLatS0), DegToRad(dLonS), DegToRad(dLatS)),
			EarthRadius * GreatCircleDistance_Rad(DegToRad(dLonS0), DegToRad(dLatS0), DegToRad(vecXs[*itPath].lon), DegToRad(vecXs[*itPath].lat)));
*/
		//vecXs[*itPath].lon = dLonS;
		//vecXs[*itPath].lat = dLatS;


	}
}

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
		_EXCEPTIONT("At present LagrangianParcelTracker only supports serial execution.");
	}
#endif

	// Input nodefile
	std::string strInputNodeFile;

	// Input nodefile column format
	std::string strInputFormat;

	// Input data
	std::string strInputData;

	// Input data list
	std::string strInputDataList;

	// Connectivity
	std::string strConnectivity;

	// Data is regional
	bool fRegional;

	// Output nodefile
	std::string strOutputNodeFile;

	// Output nodefile column format
	std::string strOutputFormat;

	// Output nodefile format
	std::string strOutputFileFormat;

	// U variable
	std::string strUVariable;

	// V variable
	std::string strVVariable;

	// Omega variable
	std::string strWVariable;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_nodefile", "");
		CommandLineString(strInputFormat, "in_fmt", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");
		CommandLineString(strOutputNodeFile, "out_nodefile", "");
		CommandLineStringD(strOutputFileFormat, "out_nodefile_format", "gfdl", "[gfdl|csv|csvnohead]");
		CommandLineString(strUVariable, "uvar", "U");
		CommandLineString(strVVariable, "vvar", "V");
		CommandLineString(strWVariable, "omegavar", "OMEGA");

		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if (strInputNodeFile.length() == 0) {
		_EXCEPTIONT("No input nodefile (--in_nodefile) specified");
	}
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input file (--in_data) or (--in_data_list) specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list) may be specified");
	}
	if (strOutputNodeFile.length() == 0) {
		_EXCEPTIONT("No output nodefile (--out_nodefile) specified");
	}
	if ((strOutputFileFormat != "gfdl") &&
		(strOutputFileFormat != "csv") &&
		(strOutputFileFormat != "csvnohead")
	) {
		_EXCEPTIONT("Output format must be either \"gfdl\", \"csv\", or \"csvnohead\"");
	}
/*
	// Create Variable registry and register velocity variables
	VariableRegistry varreg;

	VariableIndex varixU = varreg.FindOrRegister(strUVariable);
	if (varixU == InvalidVariableIndex) {
		_EXCEPTION1("Unable to register variable \"%s\"", strUVariable.c_str());
	}

	VariableIndex varixV = varreg.FindOrRegister(strVVariable);
	if (varixV == InvalidVariableIndex) {
		_EXCEPTION1("Unable to register variable \"%s\"", strVVariable.c_str());
	}

	VariableIndex varixOmega = varreg.FindOrRegister(strWVariable);
	if (varixOmega == InvalidVariableIndex) {
		_EXCEPTION1("Unable to register variable \"%s\"", strWVariable.c_str());
	}
*/
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

	// Input file type
	NodeFile::PathType iftype = NodeFile::PathTypeSN;

	// Parse --in_fmt string
	ColumnDataHeader cdhInput;
	cdhInput.Parse(strInputFormat);

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

	// Build the KD tree on this grid
	AnnounceStartBlock("Building KD tree");
	grid.BuildKDTree();
	AnnounceEndBlock("Done");

	// Get time dimension over all files
	// strOutTimeUnits is either predetermined or set at the command line
	AnnounceStartBlock("Loading time array from files");
	std::map<Time, std::pair<int, int> > mapGlobalTimeIxToFileTimeIx;

	for (int f = 0; f < vecInputFiles.size(); f++){

		// Load in file vector
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[f]);
		_ASSERT(vecNcFiles.size() > 0);

		// Load in CF-compliant time data
		const NcTimeDimension & vecTimes = vecNcFiles.GetNcTimeDimension(0);
		if (vecTimes.size() == 0) {
			_EXCEPTION1("WARNING: File group does not contain any time data (%s)",
				vecInputFiles[f].c_str());
		}

		// Store times as a Time -> <fileix,timeix> map
		for (int t = 0; t < vecTimes.size(); t++) {
			mapGlobalTimeIxToFileTimeIx.insert(
				std::map<Time, std::pair<int, int> >::value_type(
					vecTimes[t], std::pair<int, int>(f, t)));
		}
	}
	AnnounceEndBlock("Done");

	// Get the calendar type
	Time::CalendarType caltype;
	{
		if (mapGlobalTimeIxToFileTimeIx.size() == 0) {
			_EXCEPTIONT("No times found in input data files");
		}
		auto it = mapGlobalTimeIxToFileTimeIx.begin();
		caltype = it->first.GetCalendarType();
	}

	// Read contents of NodeFile
	NodeFile nodefilein;
	std::multimap<Time, LLPPosition> mapStartPositions;

	AnnounceStartBlock("Reading nodefile");
	nodefilein.Read(
		strInputNodeFile,
		NodeFile::PathTypeSN,
		cdhInput,
		grid,
		caltype);

	if (nodefilein.m_pathvec.size() != 1) {
		_EXCEPTIONT("--in_nodefile must contain exactly one path");
	}

	{
		int ixLonCDH = cdhInput.GetIndexFromString("lon");
		if (ixLonCDH == ColumnDataHeader::InvalidIndex) {
			_EXCEPTIONT("Unable to find header \"lon\" in --in_fmt");
		}
		int ixLatCDH = cdhInput.GetIndexFromString("lat");
		if (ixLatCDH == ColumnDataHeader::InvalidIndex) {
			_EXCEPTIONT("Unable to find header \"lat\" in --in_fmt");
		}
		int ixPresCDH = cdhInput.GetIndexFromString("pres");
		if (ixPresCDH == ColumnDataHeader::InvalidIndex) {
			_EXCEPTIONT("Unable to find header \"pres\" in --in_fmt");
		}

		for (size_t n = 0; n < nodefilein.m_pathvec[0].size(); n++) {
			LLPPosition llp;
			llp.ix = n;
			llp.lon = nodefilein.m_pathvec[0][n].GetColumnDataAsDouble(ixLonCDH);
			llp.lat = nodefilein.m_pathvec[0][n].GetColumnDataAsDouble(ixLatCDH);
			llp.pres = nodefilein.m_pathvec[0][n].GetColumnDataAsDouble(ixPresCDH);
			mapStartPositions.insert(
				std::multimap<Time, LLPPosition>::value_type(
					nodefilein.m_pathvec[0][n].m_time, llp));
		}
	}

	AnnounceEndBlock("Done");

	// Output NodeFile
	NodeFile nodefileout;
	nodefileout.m_pathvec.resize(mapStartPositions.size());
	nodefileout.m_cdh.push_back("lon");
	nodefileout.m_cdh.push_back("lat");
	nodefileout.m_cdh.push_back("pres");

	// Pressure array
	DataArray1D<double> dLevel;

	// Begin tracking
	AnnounceStartBlock("Beginning tracking");

	// Current position in the start position map
	auto itStartPos = mapStartPositions.rbegin();

	// Set of active path indices and coordinates of these parcels
	std::set<int> setActivePathIxs;

	// TODO: Convert LLPPosition to (X,Y,Z,pres) to reduce trigonometric operations
	std::vector<LLPPosition> vecX(mapStartPositions.size());
	std::vector<LLPPosition> vecXs(mapStartPositions.size());

	// Map from path index to the bounding levels of the parcel
	std::map<int, std::pair<int,int> > mapPathIxToBoundingLevelIx;

	// Map from level index to the storage index in the storage arrays
	std::map<int,int> mapLevelIxToLoadIx;

	std::vector<DataArray1D<float> *> vecdataU;
	std::vector<DataArray1D<float> *> vecdataV;
	std::vector<DataArray1D<float> *> vecdataW;

	NcFileVector vecNcFiles;
	int fcurix = (-1);

	// NetCDF variables storing velocity fields
	NcVar * varU;
	NcVar * varV;
	NcVar * varW;

	// Loop through all time intervals
	size_t sTimeStep = 0;
	for (auto itTime = mapGlobalTimeIxToFileTimeIx.rbegin(); itTime != mapGlobalTimeIxToFileTimeIx.rend(); itTime++, sTimeStep++) {

		// Next time
		auto itTimeNext = itTime;
		itTimeNext++;
		if (itTimeNext == mapGlobalTimeIxToFileTimeIx.rend()) {
			break;
		}
		double dDeltaTime = itTimeNext->first.DeltaSeconds(itTime->first);

		AnnounceStartBlock("Time %s => %s (%.0f seconds)",
			itTime->first.ToString().c_str(),
			itTimeNext->first.ToString().c_str(),
			dDeltaTime);

		// Generate new paths
		while ((itStartPos != mapStartPositions.rend()) && (itStartPos->first >= itTime->first)) {

			int ixNode = itStartPos->second.ix;

			Announce("Starting trajectory %i at lon %3.5f lat %3.5f pres %3.5f",
				ixNode,
				itStartPos->second.lon,
				itStartPos->second.lat,
				itStartPos->second.pres);

			if (itStartPos->first != itTime->first) {
				_EXCEPTIONT("Not Implemented: Currently we require starting point of trajectory to have a corresponding timeslice in data");
			}

			setActivePathIxs.insert(ixNode);

			vecX[ixNode] = itStartPos->second;

			PathNode pn;
			pn.m_gridix = 0;
			pn.m_fileix = 0;
			pn.m_time = itStartPos->first;
			pn.m_vecColumnData.push_back(new ColumnDataDouble(itStartPos->second.lon));
			pn.m_vecColumnData.push_back(new ColumnDataDouble(itStartPos->second.lat));
			pn.m_vecColumnData.push_back(new ColumnDataDouble(itStartPos->second.pres));
			nodefileout.m_pathvec[ixNode].push_back(pn);

			itStartPos++;
		}

		// No active paths to track
		if (setActivePathIxs.size() == 0) {
			AnnounceEndBlock("No parcels (continuing)");
			continue;
		}

		// Load data
		int f = itTime->second.first;
		int t = itTime->second.second;

		// Switch files
		if (fcurix != f) {
			fcurix = f;

			SwitchFiles(
				vecNcFiles,
				vecInputFiles[f],
				strUVariable,
				strVVariable,
				strWVariable,
				dLevel,
				&varU,
				&varV,
				&varW);
		}

		// Load grid data needed for calculating velocity vectors
		LoadGridData(
			grid,
			f, t,
			dLevel,
			setActivePathIxs,
			vecX,
			mapPathIxToBoundingLevelIx,
			mapLevelIxToLoadIx,
			varU, varV, varW,
			vecdataU, vecdataV, vecdataW,
			true); // Keep any existing data levels from last iteration

		// Calculate the predicted position
		{
			Announce("Calculating predicted position");
			vecXs = vecX;

			UpdatePaths(
				grid,
				dLevel,
				setActivePathIxs,
				mapPathIxToBoundingLevelIx,
				mapLevelIxToLoadIx,
				vecdataU, vecdataV, vecdataW,
				dDeltaTime,
				vecX,
				vecXs);

			//vecX = vecXs;

			// Calculate half position and store in vecXs
			for (auto itPath = setActivePathIxs.begin(); itPath != setActivePathIxs.end(); itPath++) {
				Node nodeX, nodeXs;
				RLLtoXYZ_Deg(vecX[*itPath].lon, vecX[*itPath].lat, nodeX.x, nodeX.y, nodeX.z);
				RLLtoXYZ_Deg(vecXs[*itPath].lon, vecXs[*itPath].lat, nodeXs.x, nodeXs.y, nodeXs.z);

				nodeX.x = 0.5 * (nodeX.x + nodeXs.x);
				nodeX.y = 0.5 * (nodeX.y + nodeXs.y);
				nodeX.z = 0.5 * (nodeX.z + nodeXs.z);

				nodeX /= nodeX.Magnitude();

				XYZtoRLL_Deg(nodeX.x, nodeX.y, nodeX.z, vecX[*itPath].lon, vecX[*itPath].lat);
			}
		}

		// Update file position
		f = itTimeNext->second.first;
		t = itTimeNext->second.second;

		// Switch files
		if (fcurix != f) {
			fcurix = f;

			SwitchFiles(
				vecNcFiles,
				vecInputFiles[f],
				strUVariable,
				strVVariable,
				strWVariable,
				dLevel,
				&varU,
				&varV,
				&varW);
		}

		// Load grid data needed for calculating velocity vectors
		LoadGridData(
			grid,
			f, t,
			dLevel,
			setActivePathIxs,
			vecX,
			mapPathIxToBoundingLevelIx,
			mapLevelIxToLoadIx,
			varU, varV, varW,
			vecdataU, vecdataV, vecdataW,
			false); // Always reload data

		// Calculate the corrected position
		{

			Announce("Calculating corrected position");
			UpdatePaths(
				grid,
				dLevel,
				setActivePathIxs,
				mapPathIxToBoundingLevelIx,
				mapLevelIxToLoadIx,
				vecdataU, vecdataV, vecdataW,
				0.5 * dDeltaTime,
				vecXs,
				vecX);

			// Insert into path array
			for (auto itPath = setActivePathIxs.begin(); itPath != setActivePathIxs.end(); itPath++) {
				PathNode pn;
				pn.m_gridix = 0;
				pn.m_fileix = 0;
				pn.m_time = itTimeNext->first;
				pn.m_vecColumnData.push_back(new ColumnDataDouble(vecX[*itPath].lon));
				pn.m_vecColumnData.push_back(new ColumnDataDouble(vecX[*itPath].lat));
				pn.m_vecColumnData.push_back(new ColumnDataDouble(vecX[*itPath].pres));
				nodefileout.m_pathvec[*itPath].push_back(pn);
			}
			//if (sTimeStep == 4) {
			//	break;
			//}
		}

		AnnounceEndBlock("Done");
	}
	AnnounceEndBlock("Done");

	// Clear data
	for (int i = 0; i < vecdataU.size(); i++) {
		delete vecdataU[i];
		delete vecdataV[i];
		delete vecdataW[i];
	}

	// Write nodefile
	AnnounceStartBlock("Write nodefile");
	nodefileout.Write(strOutputNodeFile);
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

