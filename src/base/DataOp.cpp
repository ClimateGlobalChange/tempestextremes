///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataOp.cpp
///	\author  Paul Ullrich
///	\version July 22, 2018
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "Exception.h"
#include "Announce.h"
#include "DataOp.h"
#include "Variable.h"
#include "SimpleGrid.h"
#include "kdtree.h"

#include <cstdlib>
#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////
// DataOpManager
///////////////////////////////////////////////////////////////////////////////

DataOpManager::~DataOpManager() {
	for (iterator iter = begin(); iter != end(); iter++) {
		delete iter->second;
	}
}

///////////////////////////////////////////////////////////////////////////////

DataOp * DataOpManager::Add(DataOp * pdo) {
	if (pdo == NULL) {
		_EXCEPTIONT("Invalid pointer");
	}

	iterator iter = find(pdo->GetName());
	if (iter != end()) {
		_EXCEPTION1("DataOp with given name \"%s\" already exists",
			pdo->GetName().c_str());
	}

	insert(DataOpMapPair(pdo->GetName(), pdo));

	return pdo;
}

///////////////////////////////////////////////////////////////////////////////

DataOp * DataOpManager::Add(
	const std::string & strName
) {
	if (strName == "_VECMAG") {
		return Add(new DataOp_VECMAG);

	} else if (strName == "_ABS") {
		return Add(new DataOp_ABS);

	} else if (strName == "_SIGN") {
		return Add(new DataOp_SIGN);

	} else if (strName == "_AVG") {
		return Add(new DataOp_AVG);

	} else if (strName == "_DIFF") {
		return Add(new DataOp_DIFF);
		
	} else if (strName == "_MULT") {
		return Add(new DataOp_MULT);

	} else if (strName == "_DIV") {
		return Add(new DataOp_DIV);
	
	} else if (strName == "_LAT") {
		return Add(new DataOp_LAT);

	} else if (strName == "_F") {
		return Add(new DataOp_F);

	} else if (strName.substr(0,10) == "_LAPLACIAN") {
		int nPoints = 0;
		double dDist = 0.0;

		int iMode = 0;
		int iLast = 0;
		for (int i = 0; i < strName.length(); i++) {
			if (iMode == 0) {
				if (strName[i] == '{') {
					iLast = i;
					iMode = 1;
				}
				continue;

			} else if (iMode == 1) {
				if (strName[i] == ',') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _LAPLACIAN{npts,dist} name");
					}
					nPoints = atoi(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 2;
				}

			} else if (iMode == 2) {
				if (strName[i] == '}') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _LAPLACIAN{npts,dist} name");
					}
					dDist = atof(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 3;
				}

			} else if (iMode == 3) {
				_EXCEPTIONT("Malformed _LAPLACIAN{npts,dist} name");
			}
		}
		if (iMode != 3) {
			_EXCEPTIONT("Malformed _LAPLACIAN{npts,dist} name");
		}

		return Add(new DataOp_LAPLACIAN(strName, nPoints, dDist));
	}

	return NULL;
}

///////////////////////////////////////////////////////////////////////////////

DataOp * DataOpManager::Find(
	const std::string & strName
) {
	iterator iter = find(strName);
	if (iter == end()) {
		return NULL;
	}
	return iter->second;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp
///////////////////////////////////////////////////////////////////////////////

bool DataOp::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & out
) {
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_VECMAG
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_VECMAG::name = "_VECMAG";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_VECMAG::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataVector<float> & dataLeft  = *(vecArgData[0]);
	const DataVector<float> & dataRight = *(vecArgData[1]);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] =
			sqrt(dataLeft[i] * dataLeft[i]
				+ dataRight[i] * dataRight[i]);
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_ABS
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_ABS::name = "_ABS";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_ABS::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (vecArgData[0] == NULL) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataVector<float> & data = *(vecArgData[0]);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = fabs(data[i]);
	}
	
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_SIGN
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_SIGN::name = "_SIGN";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_SIGN::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (vecArgData[0] == NULL) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataVector<float> & data = *(vecArgData[0]);

	for (int i = 0; i < dataout.GetRows(); i++) {
		if (data[i] > 0.0) {
			dataout[i] = 1.0;
		} else if (data[i] < 0.0) {
			dataout[i] = -1.0;
		} else {
			dataout[i] = 0.0;
		}
	}
	
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_AVG
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_AVG::name = "_AVG";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_AVG::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() <= 1) {
		_EXCEPTION2("%s expects at least two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			_EXCEPTION1("Arguments to %s must be data variables",
				m_strName.c_str());
		}
	}

	dataout.Zero();
	for (int v = 0; v < vecArgData.size(); v++) {
		const DataVector<float> & data  = *(vecArgData[v]);

		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] += data[i];
		}
	}

	const double dScale = 1.0 / static_cast<double>(strArg.size());
	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] *= dScale;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_DIFF
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_DIFF::name = "_DIFF";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_DIFF::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataVector<float> & dataLeft  = *(vecArgData[0]);
	const DataVector<float> & dataRight = *(vecArgData[1]);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = dataLeft[i] - dataRight[i];
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_MULT
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_MULT::name = "_MULT";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_MULT::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataVector<float> & dataLeft  = *(vecArgData[0]);
	const DataVector<float> & dataRight = *(vecArgData[1]);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = dataLeft[i] * dataRight[i];
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_DIV
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_DIV::name = "_DIV";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_DIV::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataVector<float> & dataLeft  = *(vecArgData[0]);
	const DataVector<float> & dataRight = *(vecArgData[1]);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = dataLeft[i] / dataRight[i];
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_LAT
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_LAT::name = "_LAT";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_LAT::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 0) {
		_EXCEPTION2("%s expects zero arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = grid.m_dLat[i] * 180.0 / M_PI;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_F
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_F::name = "_F";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_F::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	static const double Omega = 7.2921e-5;

	if (strArg.size() != 0) {
		_EXCEPTION2("%s expects zero arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = 2.0 * Omega * sin(grid.m_dLat[i]);
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_LAPLACIAN
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate an array of points that are of equal distance along the
///		surface of the sphere from a given central point.
///	</summary>
///	<param name="dDist">
///		Great circle distance in degrees.
///	</param>
void GenerateEqualDistanceSpherePoints(
	double dX0,
	double dY0,
	double dZ0,
	int nPoints,
	double dDist,
	std::vector<double> & dXout,
	std::vector<double> & dYout,
	std::vector<double> & dZout
) { 
	dXout.resize(nPoints);
	dYout.resize(nPoints);
	dZout.resize(nPoints);

	// Pick a quasi-arbitrary reference direction
	double dX1;
	double dY1;
	double dZ1;
	if ((fabs(dX0) >= fabs(dY0)) && (fabs(dX0) >= fabs(dZ0))) {
		dX1 = dX0;
		dY1 = dY0 + 1.0;
		dZ1 = dZ0;
	} else if ((fabs(dY0) >= fabs(dX0)) && (fabs(dY0) >= fabs(dZ0))) {
		dX1 = dX0;
		dY1 = dY0;
		dZ1 = dZ0 + 1.0;
	} else {
		dX1 = dX0 + 1.0;
		dY1 = dY0;
		dZ1 = dZ0;
	}

	// Project perpendicular to detection location
	double dDot = dX1 * dX0 + dY1 * dY0 + dZ1 * dZ0;

	dX1 -= dDot * dX0;
	dY1 -= dDot * dY0;
	dZ1 -= dDot * dZ0;

	// Normalize
	double dMag1 = sqrt(dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1);

	if (dMag1 < 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	double dScale1 = tan(dDist * M_PI / 180.0) / dMag1;

	dX1 *= dScale1;
	dY1 *= dScale1;
	dZ1 *= dScale1;

	// Verify dot product is zero
	dDot = dX0 * dX1 + dY0 * dY1 + dZ0 * dZ1;
	if (fabs(dDot) > 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	// Cross product (magnitude automatically 
	double dCrossX = dY0 * dZ1 - dZ0 * dY1;
	double dCrossY = dZ0 * dX1 - dX0 * dZ1;
	double dCrossZ = dX0 * dY1 - dY0 * dX1;

	// Generate all points
	for (int j = 0; j < nPoints; j++) {

		// Angle of rotation
		double dAngle = 2.0 * M_PI
			* static_cast<double>(j)
			/ static_cast<double>(nPoints);

		// Calculate new rotated vector
		double dX2 = dX0 + dX1 * cos(dAngle) + dCrossX * sin(dAngle);
		double dY2 = dY0 + dY1 * cos(dAngle) + dCrossY * sin(dAngle);
		double dZ2 = dZ0 + dZ1 * cos(dAngle) + dCrossZ * sin(dAngle);

		double dMag2 = sqrt(dX2 * dX2 + dY2 * dY2 + dZ2 * dZ2);

		dX2 /= dMag2;
		dY2 /= dMag2;
		dZ2 /= dMag2;

		dXout[j] = dX2;
		dYout[j] = dY2;
		dZout[j] = dZ2;
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build a sparse Laplacian operator on an unstructured SimpleGrid.
///	</summary>
///	<param name="dLaplacianDist">
///		Great circle radius of the Laplacian operator, in degrees.
///	</param>
void BuildLaplacianOperator(
	const SimpleGrid & grid,
	int nLaplacianPoints,
	double dLaplacianDist,
	SparseMatrix<float> & opLaplacian
) {
	opLaplacian.Clear();

	int iRef = 0;

	// Scaling factor used in Laplacian calculation
	const double dScale = 4.0 / static_cast<double>(nLaplacianPoints);

	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataVector<double> dXi(grid.GetSize());
	DataVector<double> dYi(grid.GetSize());
	DataVector<double> dZi(grid.GetSize());
	for (int i = 0; i < grid.GetSize(); i++) {
		double dLat = grid.m_dLat[i];
		double dLon = grid.m_dLon[i];

		dXi[i] = cos(dLon) * cos(dLat);
		dYi[i] = sin(dLon) * cos(dLat);
		dZi[i] = sin(dLat);

		kd_insert3(kdGrid, dXi[i], dYi[i], dZi[i], (void*)((&iRef)+i));
	}

	// Construct the Laplacian operator using SPH
	for (int i = 0; i < grid.GetSize(); i++) {

		// Generate points for the Laplacian
		std::vector<double> dXout;
		std::vector<double> dYout;
		std::vector<double> dZout;

		GenerateEqualDistanceSpherePoints(
			dXi[i], dYi[i], dZi[i],
			nLaplacianPoints,
			dLaplacianDist,
			dXout, dYout, dZout);

		dXout.push_back(dXi[i]);
		dYout.push_back(dYi[i]);
		dZout.push_back(dZi[i]);
/*
		kdres * kdr = kd_nearest_range3(kdGrid, dXi[i], dYi[i], dZi[i], dMaxDist);
		if (kdr == NULL) {
			_EXCEPTIONT("Error in kd_nearest_range3");
		}
		int nNodes = kd_res_size(kdr);

		std::cout << nNodes << " found at distance " << dMaxDist << std::endl;
*/
		//std::vector<int> vecCols;
		//std::vector<double> vecWeight;
		std::set<int> setPoints;
		setPoints.insert(i);

		double dAccumulatedDiff = 0.0;

		for (int j = 0; j < dXout.size(); j++) {
			// Find the nearest grid point to the output point
			kdres * kdr = kd_nearest3(kdGrid, dXout[j], dYout[j], dZout[j]);
			if (kdr == NULL) {
				_EXCEPTIONT("NULL return value in call to kd_nearest3");
			}

			void* pData = kd_res_item_data(kdr);
			if (pData == NULL) {
				_EXCEPTIONT("NULL data index");
			}
			int k = ((int*)(pData)) - (&iRef);

			kd_res_free(kdr);

			if (k == i) {
				continue;
			}

			if (k > dXi.GetRows()) {
				_EXCEPTIONT("Invalid point index");
			}

			// Ensure points are not duplicated
			if (setPoints.find(k) != setPoints.end()) {
				continue;
			} else {
				setPoints.insert(k);
			}
	/*
			// ============== BEGIN DEBUGGING =============================
			double dLon0 = atan2(dYi[i], dXi[i]) * 180.0 / M_PI;
			double dLat0 = asin(dZi[i]) * 180.0 / M_PI;

			double dLon1 = atan2(dYout[j], dXout[j]) * 180.0 / M_PI;
			double dLat1 = asin(dZout[j]) * 180.0 / M_PI;

			double dDist0 =
				sqrt(
					(dXout[j] - dXi[i]) * (dXout[j] - dXi[i])
					+ (dYout[j] - dYi[i]) * (dYout[j] - dYi[i])
					+ (dZout[j] - dZi[i]) * (dZout[j] - dZi[i]));

			std::cout << 2.0 * sin(dDist0 / 2.0) * 180.0 / M_PI << std::endl;

			double dDist1 =
				sqrt(
					(dXi[k] - dXi[i]) * (dXi[k] - dXi[i])
					+ (dYi[k] - dYi[i]) * (dYi[k] - dYi[i])
					+ (dZi[k] - dZi[i]) * (dZi[k] - dZi[i]));

			std::cout << 2.0 * sin(dDist1 / 2.0) * 180.0 / M_PI << std::endl;

			double dLon2 = atan2(dYi[k], dXi[k]) * 180.0 / M_PI;
			double dLat2 = asin(dZi[k]) * 180.0 / M_PI;

			printf("XY: %1.3f %1.3f :: %1.3f %1.3f :: %1.3f %1.3f\n", dXi[i], dYi[i], dXout[j], dYout[j], dXi[k], dYi[k]);
			printf("LL: %1.2f %1.2f :: %1.2f %1.2f\n", dLon0, dLat0, dLon1, dLat1);
			// ============== END DEBUGGING =============================
*/

			double dX1 = dXi[k] - dXi[i];
			double dY1 = dYi[k] - dYi[i];
			double dZ1 = dZi[k] - dZi[i];

			double dChordDist2 = dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1;

			double dSurfDist2 = 2.0 * asin(0.5 * sqrt(dChordDist2));
			dSurfDist2 *= dSurfDist2;

			//double dTanDist2 = (2.0 - dChordDist2);
			//dTanDist2 = dChordDist2 * (4.0 - dChordDist2) / (dTanDist2 * dTanDist2);

			opLaplacian(i,k) = static_cast<float>(dScale / dSurfDist2);

			//printf("(%1.2e %1.2e %1.2e) (%1.2e %1.2e %1.2e) %1.5e\n", dXout[j], dYout[j], dZout[j], dXi[k], dYi[k], dZi[k], dSurfDist2 * 180.0 / M_PI);
			//printf("%1.5e %i %i %1.5e\n", sqrt(dSurfDist2) * 180.0 / M_PI, i, k, opLaplacian(i,k));

			dAccumulatedDiff += dScale / dSurfDist2;
		}

		opLaplacian(i,i) = static_cast<float>(- dAccumulatedDiff);

		//printf("%i %i %1.5e\n", i, i, opLaplacian(i,i));

		if (setPoints.size() < 5) {
			Announce("WARNING: Only %i points used for Laplacian in cell %i"
				" -- accuracy may be affected", setPoints.size(), i);
		}
	}
/*
	for (
		SparseMatrix<double>::SparseMapIterator iter = opLaplacian.begin();
		iter != opLaplacian.end(); iter++
	) {
		std::cout << iter->first.first << " " << iter->first.second << " " << iter->second << std::endl;
	}
*/
	kd_free(kdGrid);
}

///////////////////////////////////////////////////////////////////////////////

DataOp_LAPLACIAN::DataOp_LAPLACIAN(
	const std::string & strName,
	int nLaplacianPoints,
	double dLaplacianDist
) :
	DataOp(strName),
	m_nLaplacianPoints(nLaplacianPoints),
	m_dLaplacianDist(dLaplacianDist),
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

bool DataOp_LAPLACIAN::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataVector<float> const *> & vecArgData,
	DataVector<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (vecArgData[0] == NULL) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	if (!m_fInitialized) {
		Announce("Building Laplacian operator %s (%i, %1.2f)",
			m_strName.c_str(), m_nLaplacianPoints, m_dLaplacianDist);

		BuildLaplacianOperator(
			grid,
			m_nLaplacianPoints,
			m_dLaplacianDist,
			m_opLaplacian);

		m_fInitialized = true;
	}

	m_opLaplacian.Apply(*(vecArgData[0]), dataout);

	return true;
}

///////////////////////////////////////////////////////////////////////////////

