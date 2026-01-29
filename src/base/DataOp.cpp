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
#include "STLStringHelper.h"
#include "kdtree.h"
#include "Constants.h"
#include "CoordTransforms.h"
#include "Units.h"

#include <cstdlib>
#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

double SatVapPres_FromCC_degC(
	double dTempDegC
) {
	return 6.1094 * exp(17.625 * dTempDegC / (dTempDegC + 243.04));
}

///////////////////////////////////////////////////////////////////////////////

double SatVapPres_FromCC_K(
	double dTempK
) {
	return 6.1094 * exp(17.625 * (dTempK - 273.15) / (dTempK - 30.11));
}

///////////////////////////////////////////////////////////////////////////////

double SatVapPres_FromCC_degF(
	double dTempDegF
) {
	return SatVapPres_FromCC_degC((5.0/9.0)*(dTempDegF - 32.0));
}

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

	} else if (strName == "_ALLPOS") {
		return Add(new DataOp_ALLPOS);

	} else if (strName == "_SUM") {
		return Add(new DataOp_SUM);

	} else if (strName == "_AVG") {
		return Add(new DataOp_AVG);

	} else if (strName == "_DIFF") {
		return Add(new DataOp_DIFF);
		
	} else if (strName == "_PROD") {
		return Add(new DataOp_PROD);

	} else if (strName == "_DIV") {
		return Add(new DataOp_DIV);
		
	} else if (strName == "_MIN") {
		return Add(new DataOp_MIN);

	} else if (strName == "_MAX") {
		return Add(new DataOp_MAX);

	} else if (strName == "_COND") {
		return Add(new DataOp_COND);

	} else if (strName == "_EQUALS") {
		return Add(new DataOp_EQUALS);

	} else if (strName == "_SQRT") {
		return Add(new DataOp_SQRT);
	
	} else if (strName == "_POW") {
		return Add(new DataOp_POW);

	} else if (strName == "_LAT") {
		return Add(new DataOp_LAT);

	} else if (strName == "_LON") {
		return Add(new DataOp_LON);

	} else if (strName == "_AREA") {
		return Add(new DataOp_AREA);

	} else if (strName == "_F") {
		return Add(new DataOp_F);

	// Laplacian operator
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

	// Curl operator
	} else if (strName.substr(0,5) == "_CURL") {
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
						_EXCEPTIONT("Malformed _CURL{npts,dist} name");
					}
					nPoints = atoi(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 2;
				}

			} else if (iMode == 2) {
				if (strName[i] == '}') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _CURL{npts,dist} name");
					}
					dDist = atof(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 3;
				}

			} else if (iMode == 3) {
				_EXCEPTIONT("Malformed _CURL{npts,dist} name");
			}
		}
		if (iMode != 3) {
			_EXCEPTIONT("Malformed _CURL{npts,dist} name");
		}

		return Add(new DataOp_CURL(strName, nPoints, dDist));

	// Divergence operator
	} else if (strName.substr(0,11) == "_DIVERGENCE") {
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
						_EXCEPTIONT("Malformed _DIVERGENCE{npts,dist} name");
					}
					nPoints = atoi(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 2;
				}

			} else if (iMode == 2) {
				if (strName[i] == '}') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _DIVERGENCE{npts,dist} name");
					}
					dDist = atof(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 3;
				}

			} else if (iMode == 3) {
				_EXCEPTIONT("Malformed _DIVERGENCE{npts,dist} name");
			}
		}
		if (iMode != 3) {
			_EXCEPTIONT("Malformed _DIVERGENCE{npts,dist} name");
		}

		return Add(new DataOp_DIVERGENCE(strName, nPoints, dDist));

	// Gradient magnitude operator
	} else if (strName.substr(0,8) == "_GRADMAG") {
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
						_EXCEPTIONT("Malformed _GRADMAG{npts,dist} name");
					}
					nPoints = atoi(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 2;
				}

			} else if (iMode == 2) {
				if (strName[i] == '}') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _GRADMAG{npts,dist} name");
					}
					dDist = atof(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 3;
				}

			} else if (iMode == 3) {
				_EXCEPTIONT("Malformed _GRADMAG{npts,dist} name");
			}
		}
		if (iMode != 3) {
			_EXCEPTIONT("Malformed _GRADMAG{npts,dist} name");
		}

		return Add(new DataOp_GRADMAG(strName, nPoints, dDist));

	// Vector dot gradient operator
	} else if (strName.substr(0,11) == "_VECDOTGRAD") {
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
						_EXCEPTIONT("Malformed _VECDOTGRAD{npts,dist} name");
					}
					nPoints = atoi(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 2;
				}

			} else if (iMode == 2) {
				if (strName[i] == '}') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _VECDOTGRAD{npts,dist} name");
					}
					dDist = atof(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 3;
				}

			} else if (iMode == 3) {
				_EXCEPTIONT("Malformed _VECDOTGRAD{npts,dist} name");
			}
		}
		if (iMode != 3) {
			_EXCEPTIONT("Malformed _VECDOTGRAD{npts,dist} name");
		}

		return Add(new DataOp_VECDOTGRAD(strName, nPoints, dDist));

	// Mean operator
	} else if (strName.substr(0,5) == "_MEAN") {
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
				if (strName[i] == '}') {
					if ((i-iLast-1) < 1) {
						_EXCEPTIONT("Malformed _MEAN{dist} name");
					}
					dDist = atof(strName.substr(iLast+1, i-iLast-1).c_str());
					iLast = i;
					iMode = 2;
				}

			} else if (iMode == 2) {
				_EXCEPTIONT("Malformed _MEAN{dist} name");
			}
		}
		if (iMode != 2) {
			_EXCEPTIONT("Malformed _MEAN{dist} name");
		}

		return Add(new DataOp_MEAN(strName, dDist));

	// Chill hours calculator using tasmin and tasmax
	} else if (strName == "_DAILYCHILLHOURS") {
		return Add(new DataOp_DAILYCHILLHOURS);

	// Relative humidity using T2d and T2m
	} else if (strName == "_RELHUMFROMTDTA") {
		return Add(new DataOp_RELHUMFROMTDTA);

	// VPD from Ta and RH
	} else if (strName == "_VPDFROMTAHUR") {
		return Add(new DataOp_VPDFROMTAHUR);

	// VPD from Ta and SH
	} else if (strName == "_VPDFROMTAHUSPRES") {
		return Add(new DataOp_VPDFROMTAHUSPRES);

	} else {
		_EXCEPTION1("Invalid DataOp \"%s\"", strName.c_str());
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & out
) {
	return true;
}

///////////////////////////////////////////////////////////////////////////////

std::string DataOp::GetUnits_Common(
	const std::vector<std::string> & vecUnits
) {
	if (vecUnits.size() == 0) {
		return std::string("");
	}
	for (size_t i = 1; i < vecUnits.size(); i++) {
		if (vecUnits[i] != vecUnits[0]) {
			return std::string("");
		}
	}
	return vecUnits[0];
}

///////////////////////////////////////////////////////////////////////////////

std::string DataOp::GetUnits(
	const std::vector<std::string> & vecUnits
) {
	return "";
}

///////////////////////////////////////////////////////////////////////////////

bool DataOp::HasFillValue(
	const std::vector<DataArray1D<float> const *> vecArgData
) {
	for (int i = 0; i < vecArgData.size(); i++) {
		if (vecArgData[i] != NULL) {
			if (vecArgData[i]->HasFillValue()) {
				return true;
			}
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////

float DataOp::GetFillValue_Common(
	const std::vector<DataArray1D<float> const *> vecArgData
) {
	float dCommonFillValue = DefaultFillValue;
	bool fCommonFillValueDefined = false;
	for (int i = 0; i < vecArgData.size(); i++) {
		if (vecArgData[i] != NULL) {
			if (vecArgData[i]->HasFillValue()) {
				if (!fCommonFillValueDefined) {
					dCommonFillValue = vecArgData[i]->GetFillValue();
					fCommonFillValueDefined = true;

				} else if (dCommonFillValue != vecArgData[i]->GetFillValue()) {
					dCommonFillValue = DefaultFillValue;
				}
			}
		}
	}
	if (!fCommonFillValueDefined) {
		return DefaultFillValue;
	}
	return dCommonFillValue;
}

///////////////////////////////////////////////////////////////////////////////

float DataOp::GetFillValue(
	const std::vector<DataArray1D<float> const *> vecArgData
) {
	return DefaultFillValue;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_VECMAG
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_VECMAG::name = "_VECMAG";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_VECMAG::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataArray1D<float> & dataLeft  = *(vecArgData[0]);
	const DataArray1D<float> & dataRight = *(vecArgData[1]);

	dataout.SetFillValue(DefaultFillValue);

	if (dataLeft.HasFillValue() || dataRight.HasFillValue()) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataLeft[i] == dataLeft.GetFillValue()) ||
			    (dataRight[i] == dataRight.GetFillValue())
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] =
					sqrt(dataLeft[i] * dataLeft[i]
						+ dataRight[i] * dataRight[i]);
			}
		}

	} else {
		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] =
				sqrt(dataLeft[i] * dataLeft[i]
					+ dataRight[i] * dataRight[i]);
		}
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (vecArgData[0] == NULL) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataArray1D<float> & data = *(vecArgData[0]);

	dataout.SetFillValue(DefaultFillValue);

	if (data.HasFillValue()) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if (data[i] == data.GetFillValue()) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = fabs(data[i]);
			}
		}

	} else {
		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] = fabs(data[i]);
		}
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (vecArgData[0] == NULL) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	const DataArray1D<float> & data = *(vecArgData[0]);

	dataout.SetFillValue(DefaultFillValue);

	if (data.HasFillValue()) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((std::isnan(data[i])) || (data[i] == data.GetFillValue())) {
				dataout[i] = dataout.GetFillValue();
			} else if (data[i] > 0.0) {
				dataout[i] = 1.0;
			} else if (data[i] < 0.0) {
				dataout[i] = -1.0;
			} else {
				dataout[i] = 0.0;
			}
		}

	} else {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if (std::isnan(data[i])) {
				dataout[i] = dataout.GetFillValue();
			} else if (data[i] > 0.0) {
				dataout[i] = 1.0;
			} else if (data[i] < 0.0) {
				dataout[i] = -1.0;
			} else {
				dataout[i] = 0.0;
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_ALLPOS
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_ALLPOS::name = "_ALLPOS";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_ALLPOS::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
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

	dataout.SetFillValue(DefaultFillValue);
	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = 1.0;
	}

	for (int v = 0; v < vecArgData.size(); v++) {
		const DataArray1D<float> & data  = *(vecArgData[v]);

		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else if ((data[i] <= 0.0) && (dataout[i] != dataout.GetFillValue())) {
					dataout[i] = 0.0;
				}
			}

		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else if ((data[i] <= 0.0) && (dataout[i] != dataout.GetFillValue())) {
					dataout[i] = 0.0;
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_SUM
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_SUM::name = "_SUM";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_SUM::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() <= 1) {
		_EXCEPTION2("%s expects at least two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			if (!STLStringHelper::IsFloat(strArg[v])) {
				_EXCEPTION1("Arguments to %s must be data variables or floats",
					m_strName.c_str());
			}
		}
	}

	dataout.Zero();
	dataout.SetFillValue(DefaultFillValue);

	for (int v = 0; v < vecArgData.size(); v++) {

		if (vecArgData[v] == NULL) {
			float dValue = atof(strArg[v].c_str());
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (dataout[i] != dataout.GetFillValue()) {
					dataout[i] += dValue;
				}
			}

		} else {
			const DataArray1D<float> & data  = *(vecArgData[v]);

			if (data.HasFillValue()) {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if ((std::isnan(data[i])) ||
					    (data[i] == data.GetFillValue()) ||
					    (dataout[i] == dataout.GetFillValue())
					) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] += data[i];
					}
				}

			} else {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if ((std::isnan(data[i])) || (dataout[i] == dataout.GetFillValue())) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] += data[i];
					}
				}
			}
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
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
	dataout.SetFillValue(DefaultFillValue);

	for (int v = 0; v < vecArgData.size(); v++) {

		if (vecArgData[v] == NULL) {
			float dValue = atof(strArg[v].c_str());
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (dataout[i] != dataout.GetFillValue()) {
					dataout[i] += dValue;
				}
			}

		} else {
			const DataArray1D<float> & data  = *(vecArgData[v]);

			if (data.HasFillValue()) {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if ((std::isnan(data[i])) ||
					    (data[i] == data.GetFillValue()) ||
					    (dataout[i] == dataout.GetFillValue())
					) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] += data[i];
					}
				}

			} else {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if ((std::isnan(data[i])) || (dataout[i] == dataout.GetFillValue())) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] += data[i];
					}
				}
			}
		}
	}

	const double dScale = 1.0 / static_cast<double>(strArg.size());
	for (int i = 0; i < dataout.GetRows(); i++) {
		if (dataout[i] != dataout.GetFillValue()) {
			dataout[i] *= dScale;
		}
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) && (vecArgData[1] == NULL)) {
		_EXCEPTION1("At least one arguments to %s must be a data variable",
			m_strName.c_str());
	}
	if (vecArgData[0] == NULL) {
		if (!STLStringHelper::IsFloat(strArg[0])) {
			_EXCEPTION1("Arguments to %s must be data variables or floats",
				m_strName.c_str());
		}
	}
	if (vecArgData[1] == NULL) {
		if (!STLStringHelper::IsFloat(strArg[1])) {
			_EXCEPTION1("Arguments to %s must be data variables or floats",
				m_strName.c_str());
		}
	}

	dataout.Zero();
	dataout.SetFillValue(DefaultFillValue);

	if (vecArgData[0] == NULL) {
		float dValue = atof(strArg[0].c_str());
		const DataArray1D<float> & data = *(vecArgData[1]);
		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dValue - data[i];
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dValue - data[i];
				}
			}
		}

	} else if (vecArgData[1] == NULL) {
		const DataArray1D<float> & data = *(vecArgData[0]);
		float dValue = atof(strArg[1].c_str());
		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i] - dValue;
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i] - dValue;
				}
			}
		}

	} else {
		const DataArray1D<float> & dataLeft  = *(vecArgData[0]);
		const DataArray1D<float> & dataRight = *(vecArgData[1]);

		if (dataLeft.HasFillValue() && dataRight.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) ||
					(dataLeft[i] == dataLeft.GetFillValue()) ||
				    std::isnan(dataRight[i]) ||
					(dataRight[i] == dataRight.GetFillValue())
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] - dataRight[i];
				}
			}

		} else if (dataLeft.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) ||
					(dataLeft[i] == dataLeft.GetFillValue()) ||
				    std::isnan(dataRight[i])
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] - dataRight[i];
				}
			}

		} else if (dataRight.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) ||
				    std::isnan(dataRight[i]) ||
					(dataRight[i] == dataRight.GetFillValue())
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] - dataRight[i];
				}
			}

		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) || std::isnan(dataRight[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] - dataRight[i];
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_PROD
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_PROD::name = "_PROD";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_PROD::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() <= 1) {
		_EXCEPTION2("%s expects at least two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			if (!STLStringHelper::IsFloat(strArg[v])) {
				_EXCEPTION1("Arguments to %s must be data variables or floats",
					m_strName.c_str());
			}
		}
	}

	dataout.SetFillValue(DefaultFillValue);
	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = 1.0;
	}

	for (int v = 0; v < vecArgData.size(); v++) {

		if (vecArgData[v] == NULL) {
			float dValue = atof(strArg[v].c_str());
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (dataout[i] != dataout.GetFillValue()) {
					dataout[i] *= dValue;
				}
			}

		} else {
			const DataArray1D<float> & data  = *(vecArgData[v]);
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) ||
				    (data.HasFillValue() && (data[i] == data.GetFillValue())) ||
				    (dataout[i] == dataout.GetFillValue())
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] *= data[i];
				}
			}
		}
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) && (vecArgData[1] == NULL)) {
		_EXCEPTION1("At least one arguments to %s must be a data variable",
			m_strName.c_str());
	}
	if (vecArgData[0] == NULL) {
		if (!STLStringHelper::IsFloat(strArg[0])) {
			_EXCEPTION1("Arguments to %s must be data variables or floats",
				m_strName.c_str());
		}
	}
	if (vecArgData[1] == NULL) {
		if (!STLStringHelper::IsFloat(strArg[1])) {
			_EXCEPTION1("Arguments to %s must be data variables or floats",
				m_strName.c_str());
		}
	}

	dataout.Zero();
	dataout.SetFillValue(DefaultFillValue);

	if (vecArgData[0] == NULL) {
		float dValue = atof(strArg[0].c_str());
		const DataArray1D<float> & data = *(vecArgData[1]);
		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dValue / data[i];
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dValue / data[i];
				}
			}
		}

	} else if (vecArgData[1] == NULL) {
		const DataArray1D<float> & data = *(vecArgData[0]);
		float dValue = atof(strArg[1].c_str());
		if (dValue == 0.0) {
			_EXCEPTION1("Division by zero in %s", m_strName.c_str());
		}
		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i] / dValue;
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i] / dValue;
				}
			}
		}

	} else {
		const DataArray1D<float> & dataLeft  = *(vecArgData[0]);
		const DataArray1D<float> & dataRight = *(vecArgData[1]);

		if (dataLeft.HasFillValue() && dataRight.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) ||
					(dataLeft[i] == dataLeft.GetFillValue()) ||
				    std::isnan(dataRight[i]) ||
					(dataRight[i] == dataRight.GetFillValue())
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] / dataRight[i];
				}
			}

		} else if (dataLeft.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) ||
					(dataLeft[i] == dataLeft.GetFillValue()) ||
				    std::isnan(dataRight[i])
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] / dataRight[i];
				}
			}

		} else if (dataRight.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) ||
				    std::isnan(dataRight[i]) ||
					(dataRight[i] == dataRight.GetFillValue())
				) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] / dataRight[i];
				}
			}

		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(dataLeft[i]) || std::isnan(dataRight[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = dataLeft[i] / dataRight[i];
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_MIN
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_MIN::name = "_MIN";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_MIN::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() <= 1) {
		_EXCEPTION2("%s expects at least two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			if (!STLStringHelper::IsFloat(strArg[v])) {
				_EXCEPTION1("Arguments to %s must be data variables or floats",
					m_strName.c_str());
			}
		}
	}

	dataout.SetFillValue(DefaultFillValue);

	if (vecArgData[0] == NULL) {
		float dValue = atof(strArg[0].c_str());
		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] = dValue;
		}

	} else {
		const DataArray1D<float> & data = *(vecArgData[0]);

		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i];
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i];
				}
			}
		}
	}

	for (int v = 1; v < vecArgData.size(); v++) {

		if (vecArgData[v] == NULL) {
			float dValue = atof(strArg[v].c_str());
			for (int i = 0; i < dataout.GetRows(); i++) {
				if ((dValue < dataout[i]) && (dataout[i] != dataout.GetFillValue())) {
					dataout[i] = dValue;
				}
			}

		} else {
			const DataArray1D<float> & data = *(vecArgData[v]);
			if (data.HasFillValue()) {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if (std::isnan(data[i]) ||
					    (data[i] == data.GetFillValue()) ||
					    (dataout[i] == dataout.GetFillValue())
					) {
						dataout[i] = dataout.GetFillValue();
					} else if (data[i] < dataout[i]) {
						dataout[i] = data[i];
					}
				}

			} else {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if (std::isnan(data[i]) ||
					    (dataout[i] == dataout.GetFillValue())
					) {
						dataout[i] = dataout.GetFillValue();
					} else if (data[i] < dataout[i]) {
						dataout[i] = data[i];
					}
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_MAX
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_MAX::name = "_MAX";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_MAX::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() <= 1) {
		_EXCEPTION2("%s expects at least two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			if (!STLStringHelper::IsFloat(strArg[v])) {
				_EXCEPTION1("Arguments to %s must be data variables or floats",
					m_strName.c_str());
			}
		}
	}

	dataout.SetFillValue(DefaultFillValue);

	if (vecArgData[0] == NULL) {
		float dValue = atof(strArg[0].c_str());
		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] = dValue;
		}

	} else {
		const DataArray1D<float> & data = *(vecArgData[0]);

		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i]) || (data[i] == data.GetFillValue())) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i];
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data[i])) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = data[i];
				}
			}
		}
	}

	for (int v = 1; v < vecArgData.size(); v++) {

		if (vecArgData[v] == NULL) {
			float dValue = atof(strArg[v].c_str());
			for (int i = 0; i < dataout.GetRows(); i++) {
				if ((dValue > dataout[i]) && (dataout[i] != dataout.GetFillValue())) {
					dataout[i] = dValue;
				}
			}

		} else {
			const DataArray1D<float> & data = *(vecArgData[v]);
			if (data.HasFillValue()) {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if (std::isnan(data[i]) ||
					    (data[i] == data.GetFillValue()) ||
					    (dataout[i] == dataout.GetFillValue())
					) {
						dataout[i] = dataout.GetFillValue();
					} else if (data[i] > dataout[i]) {
						dataout[i] = data[i];
					}
				}

			} else {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if (std::isnan(data[i]) ||
					    (dataout[i] == dataout.GetFillValue())
					) {
						dataout[i] = dataout.GetFillValue();
					} else if (data[i] > dataout[i]) {
						dataout[i] = data[i];
					}
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_COND
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_COND::name = "_COND";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_COND::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 3) {
		_EXCEPTION2("%s expects three arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			if (!STLStringHelper::IsFloat(strArg[v])) {
				_EXCEPTION1("Arguments to %s must be data variables or floats",
					m_strName.c_str());
			}
		}
	}

	dataout.SetFillValue(DefaultFillValue);

	// First conditional is a float
	if (vecArgData[0] == NULL) {
		float dValue0 = atof(strArg[0].c_str());

		int ix = 2;
		if (dValue0 > 0.0) {
			ix = 1;
		}
		if (vecArgData[ix] == NULL) {
			float dValueRHS = atof(strArg[ix].c_str());
			for (int i = 0; i < dataout.GetRows(); i++) {
				dataout[i] = dValueRHS;
			}
		} else {
			const DataArray1D<float> & data = *(vecArgData[ix]);
			if (data.HasFillValue()) {
				for (int i = 0; i < dataout.GetRows(); i++) {
					if (data[i] == data.GetFillValue()) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] = data[i];
					}
				}
			} else {
				for (int i = 0; i < dataout.GetRows(); i++) {
					dataout[i] = data[i];
				}
			}
		}

	// First conditional is a field
	} else {
		const DataArray1D<float> & datacond = *(vecArgData[0]);

		// Outputs are floats
		if ((vecArgData[1] == NULL) && (vecArgData[2] == NULL)) {
			float dValue1 = atof(strArg[1].c_str());
			float dValue2 = atof(strArg[2].c_str());
			for (int i = 0; i < datacond.GetRows(); i++) {
				if (std::isnan(datacond[i]) || (datacond.HasFillValue() && (datacond[i] == datacond.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (datacond[i] > 0.0) {
					dataout[i] = dValue1;
				} else {
					dataout[i] = dValue2;
				}
			}

		// First output is a float
		} else if (vecArgData[1] == NULL) {
			float dValue1 = atof(strArg[1].c_str());
			const DataArray1D<float> & data2 = *(vecArgData[2]);
			for (int i = 0; i < datacond.GetRows(); i++) {
				if (std::isnan(datacond[i]) || (datacond.HasFillValue() && (datacond[i] == datacond.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (datacond[i] > 0.0) {
					dataout[i] = dValue1;
				} else {
					if (std::isnan(data2[i]) || (data2.HasFillValue() && (data2[i] == data2.GetFillValue()))) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] = data2[i];
					}
				}
			}

		// Second output is a float
		} else if (vecArgData[2] == NULL) {
			const DataArray1D<float> & data1 = *(vecArgData[1]);
			float dValue2 = atof(strArg[2].c_str());
			for (int i = 0; i < datacond.GetRows(); i++) {
				if (std::isnan(datacond[i]) || (datacond.HasFillValue() && (datacond[i] == datacond.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (datacond[i] > 0.0) {
					if (std::isnan(data1[i]) || (data1.HasFillValue() && (data1[i] == data1.GetFillValue()))) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] = data1[i];
					}
				} else {
					dataout[i] = dValue2;
				}
			}

		// Both outputs are fields
		} else {
			const DataArray1D<float> & data1 = *(vecArgData[1]);
			const DataArray1D<float> & data2 = *(vecArgData[2]);
			for (int i = 0; i < datacond.GetRows(); i++) {
				if (std::isnan(datacond[i]) || (datacond.HasFillValue() && (datacond[i] == datacond.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (datacond[i] > 0.0) {
					if (std::isnan(data1[i]) || (data1.HasFillValue() && (data1[i] == data1.GetFillValue()))) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] = data1[i];
					}
				} else {
					if (std::isnan(data2[i]) || (data2.HasFillValue() && (data2[i] == data2.GetFillValue()))) {
						dataout[i] = dataout.GetFillValue();
					} else {
						dataout[i] = data2[i];
					}
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_EQUALS
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_EQUALS::name = "_EQUALS";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_EQUALS::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	for (int v = 0; v < vecArgData.size(); v++) {
		if (vecArgData[v] == NULL) {
			if (!STLStringHelper::IsFloat(strArg[v])) {
				_EXCEPTION1("Arguments to %s must be data variables or floats",
					m_strName.c_str());
			}
		}
	}

	dataout.SetFillValue(DefaultFillValue);

	// First argument is a float
	if (vecArgData[0] == NULL) {
		float dValue0 = atof(strArg[0].c_str());

		// Second argument is a float
		if (vecArgData[1] == NULL) {
			float dValue1 = atof(strArg[1].c_str());

			if (dValue0 == dValue1) {
				for (int i = 0; i < dataout.GetRows(); i++) {
					dataout[i] = 1.0;
				}
			} else {
				for (int i = 0; i < dataout.GetRows(); i++) {
					dataout[i] = 0.0;
				}
			}

		// Second argument is a field
		} else {
			const DataArray1D<float> & data1 = *(vecArgData[1]);
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data1[i]) || (data1.HasFillValue() && (data1[i] == data1.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (dValue0 == data1[i]) {
					dataout[i] = 1.0;
				} else {
					dataout[i] = 0.0;
				}
			}
		}

	// First argument is a field
	} else {
		const DataArray1D<float> & data0 = *(vecArgData[0]);

		// Second argument is a float
		if (vecArgData[1] == NULL) {
			float dValue1 = atof(strArg[1].c_str());

			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data0[i]) || (data0.HasFillValue() && (data0[i] == data0.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (data0[i] == dValue1) {
					dataout[i] = 1.0;
				} else {
					dataout[i] = 0.0;
				}
			}

		// Second argument is a field
		} else {
			const DataArray1D<float> & data1 = *(vecArgData[1]);
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (std::isnan(data0[i]) || (data0.HasFillValue() && (data0[i] == data0.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (std::isnan(data1[i]) || (data1.HasFillValue() && (data1[i] == data1.GetFillValue()))) {
					dataout[i] = dataout.GetFillValue();
				} else if (data0[i] == data1[i]) {
					dataout[i] = 1.0;
				} else {
					dataout[i] = 0.0;
				}
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_EQUALS
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_SQRT::name = "_SQRT";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_SQRT::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}

	dataout.SetFillValue(DefaultFillValue);

	if (vecArgData[0] == NULL) {
		float dValue = atof(strArg[0].c_str());
		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] = sqrt(dValue);
		}

	} else {
		const DataArray1D<float> & data = *(vecArgData[0]);

		if (data.HasFillValue()) {
			for (int i = 0; i < dataout.GetRows(); i++) {
				if (data[i] == data.GetFillValue()) {
					dataout[i] = dataout.GetFillValue();
				} else {
					dataout[i] = sqrt(data[i]);
				}
			}
		} else {
			for (int i = 0; i < dataout.GetRows(); i++) {
				dataout[i] = sqrt(data[i]);
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_POW
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_POW::name = "_POW";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_POW::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	dataout.SetFillValue(DefaultFillValue);

	if ((vecArgData[0] == NULL) && (vecArgData[1] == NULL)) {
		float dValue = atof(strArg[0].c_str());
		float dExponent = atof(strArg[1].c_str());
		for (int i = 0; i < dataout.GetRows(); i++) {
			dataout[i] = pow(dValue, dExponent);
		}

	} else if ((vecArgData[0] == NULL) && (vecArgData[1] != NULL)) {
		float dValue = atof(strArg[0].c_str());
		const DataArray1D<float> & data = *(vecArgData[1]);

		for (int i = 0; i < dataout.GetRows(); i++) {
			if (data.HasFillValue() && (data[i] == data.GetFillValue())) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = pow(dValue, data[i]);
			}
		}

	} else if ((vecArgData[0] != NULL) && (vecArgData[1] == NULL)) {
		const DataArray1D<float> & data = *(vecArgData[0]);
		float dExponent = atof(strArg[1].c_str());

		for (int i = 0; i < dataout.GetRows(); i++) {
			if (data.HasFillValue() && (data[i] == data.GetFillValue())) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = pow(data[i],dExponent);
			}
		}

	} else if ((vecArgData[0] != NULL) && (vecArgData[1] == NULL)) {
		const DataArray1D<float> & dataValue = *(vecArgData[0]);
		const DataArray1D<float> & dataExponent = *(vecArgData[1]);

		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataValue.HasFillValue() && (dataValue[i] == dataValue.GetFillValue())) ||
			    (dataExponent.HasFillValue() && (dataExponent[i] == dataExponent.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = pow(dataValue[i],dataExponent[i]);
			}
		}
	
	} else {
		_EXCEPTION();
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 0) {
		_EXCEPTION2("%s expects zero arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	dataout.SetFillValue(DefaultFillValue);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = grid.m_dLat[i] * 180.0 / M_PI;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_LAT
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_LON::name = "_LON";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_LON::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 0) {
		_EXCEPTION2("%s expects zero arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	dataout.SetFillValue(DefaultFillValue);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = grid.m_dLon[i] * 180.0 / M_PI;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_AREA
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_AREA::name = "_AREA";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_AREA::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() > 1) {
		_EXCEPTION2("%s expects at most one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (!grid.HasAreas()) {
		_EXCEPTION1("Grid area not available in %s operator",
			m_strName.c_str());
	}

	dataout.SetFillValue(DefaultFillValue);

	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = grid.m_dArea[i] * EarthRadius * EarthRadius;
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	static const double Omega = 7.2921e-5;

	if (strArg.size() != 0) {
		_EXCEPTION2("%s expects zero arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	dataout.SetFillValue(DefaultFillValue);

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
	double dDistDeg,
	std::vector<double> & dXout,
	std::vector<double> & dYout,
	std::vector<double> & dZout
) { 
	dXout.resize(nPoints);
	dYout.resize(nPoints);
	dZout.resize(nPoints);

	double dX1;
	double dY1;
	double dZ1;
/*
	// Pick a quasi-arbitrary reference direction
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
*/
	// Pick eastward reference direction
	if (fabs(fabs(dZ0) - 1.0) > 1.0e-10) {
		dX1 = dY0;
		dY1 = - dX0;
		dZ1 = 0.0;

	// At poles point to grand meridian
	} else {
		dX1 = 1.0;
		dY1 = 0.0;
		dZ1 = 0.0;
	}
/*
	// Project perpendicular to detection location
	double dDot = dX1 * dX0 + dY1 * dY0 + dZ1 * dZ0;

	dX1 -= dDot * dX0;
	dY1 -= dDot * dY0;
	dZ1 -= dDot * dZ0;
*/
	// Normalize
	double dMag1 = sqrt(dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1);

	if (dMag1 < 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	double dScale1 = tan(DegToRad(dDistDeg)) / dMag1;

	dX1 *= dScale1;
	dY1 *= dScale1;
	dZ1 *= dScale1;

	// Verify dot product is zero
	double dDot = dX0 * dX1 + dY0 * dY1 + dZ0 * dZ1;
	if (fabs(dDot) > 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	// Cross product to obtain northward vector (magnitude automatically set to 1)
	double dCrossX = dZ0 * dY1 - dY0 * dZ1;
	double dCrossY = dX0 * dZ1 - dZ0 * dX1;
	double dCrossZ = dY0 * dX1 - dX0 * dY1;

	if (fabs(dZ0) < 1.0 - 1.0e-12) {
		_ASSERT(dCrossZ >= 0.0);
	}

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
	double dLaplacianDistDeg,
	SparseMatrix<float> & opLaplacian
) {
	opLaplacian.Clear();

	int iRef = 0;

	// Scaling factor used in Laplacian calculation
	const double dScale = 4.0 / static_cast<double>(nLaplacianPoints);

	// TODO: Use SimpleGrid's built-in functionality
	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataArray1D<double> dXi(grid.GetSize());
	DataArray1D<double> dYi(grid.GetSize());
	DataArray1D<double> dZi(grid.GetSize());
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
			dLaplacianDistDeg,
			dXout, dYout, dZout);
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
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
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
// DataOp_CURL
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build a sparse Curl operator on an unstructured SimpleGrid.
///	</summary>
///	<param name="dCurlDist">
///		Great circle radius of the Curl operator, in degrees.
///	</param>
void BuildCurlOperator(
	const SimpleGrid & grid,
	int nCurlPoints,
	double dCurlDistDeg,
	SparseMatrix<float> & opCurlE,
	SparseMatrix<float> & opCurlN
) {
	opCurlE.Clear();
	opCurlN.Clear();

	int iRef = 0;

	// Scaling factor used in Curl calculation
	const double dScale = 2.0 / (EarthRadius * static_cast<double>(nCurlPoints) * DegToRad(dCurlDistDeg));

	// TODO: Use SimpleGrid's built-in functionality
	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataArray1D<double> dXi(grid.GetSize());
	DataArray1D<double> dYi(grid.GetSize());
	DataArray1D<double> dZi(grid.GetSize());
	for (int i = 0; i < grid.GetSize(); i++) {
		double dLat = grid.m_dLat[i];
		double dLon = grid.m_dLon[i];

		dXi[i] = cos(dLon) * cos(dLat);
		dYi[i] = sin(dLon) * cos(dLat);
		dZi[i] = sin(dLat);

		kd_insert3(kdGrid, dXi[i], dYi[i], dZi[i], (void*)((&iRef)+i));
	}

	// Construct the Curl operator using SPH
	for (int i = 0; i < grid.GetSize(); i++) {

		// Generate points for the Curl
		std::vector<double> dXout;
		std::vector<double> dYout;
		std::vector<double> dZout;

		GenerateEqualDistanceSpherePoints(
			dXi[i], dYi[i], dZi[i],
			nCurlPoints,
			dCurlDistDeg,
			dXout, dYout, dZout);
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

		// Center point in spherical coordinates
		double dLonRadI, dLatRadI;
		XYZtoRLL_Rad(dXi[i], dYi[i], dZi[i], dLonRadI, dLatRadI);

		//printf("%1.5e %1.5e\n", RadToDeg(grid.m_dLon[i]), RadToDeg(grid.m_dLat[i]));
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

			// Local point in spherical coordinates
			double dLonRadK, dLatRadK;
			XYZtoRLL_Rad(dXi[k], dYi[k], dZi[k], dLonRadK, dLatRadK);

			// Calculate local radial vector from central lat/lon
			// i.e. project \vec{X}_k - \vec{X}_i to the surface of the
			//      sphere and normalize to unit length.
			double dRx = dXi[k] - dXi[i];
			double dRy = dYi[k] - dYi[i];
			double dRz = dZi[k] - dZi[i];

			double dChordLength = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

			if (dChordLength < 1.0e-12) {
				_EXCEPTIONT("Curl sample point not distinct from center point: increase radius of operator.");
			}

			// Calculate azimuthal velocity vector at point I
			double dDot = dRx * dXi[i] + dRy * dYi[i] + dRz * dZi[i];

			dRx -= dDot * dXi[i];
			dRy -= dDot * dYi[i];
			dRz -= dDot * dZi[i];

			double dMag = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

			_ASSERT(dMag >= 1.0e-12);

			dRx /= dMag;
			dRy /= dMag;
			dRz /= dMag;

			double dAxI = dYi[i] * dRz - dZi[i] * dRy;
			double dAyI = dZi[i] * dRx - dXi[i] * dRz;
			double dAzI = dXi[i] * dRy - dYi[i] * dRx;

			// Calculate azimuthal velocity vector at point K
			dDot = dRx * dXi[k] + dRy * dYi[k] + dRz * dZi[k];

			dRx -= dDot * dXi[k];
			dRy -= dDot * dYi[k];
			dRz -= dDot * dZi[k];

			dMag = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

			_ASSERT(dMag >= 1.0e-12);

			dRx /= dMag;
			dRy /= dMag;
			dRz /= dMag;

			double dAxK = dYi[k] * dRz - dZi[k] * dRy;
			double dAyK = dZi[k] * dRx - dXi[k] * dRz;
			double dAzK = dXi[k] * dRy - dYi[k] * dRx;

			// Great circle distance between points I and K
			double dRdeg =
				RadToDeg(GreatCircleDistanceFromChordLength_Rad(dChordLength));

			// Coefficients for linear interpolation / extrapolation to distance dCurlDistDeg
			double dCoeff0 = (dRdeg - dCurlDistDeg) / dRdeg;
			double dCoeff1 = dCurlDistDeg / dRdeg;

			// Multiplicative coefficients for converting spherical velocities to azimuthal velocities
			double dAzLonCoeffI = - sin(dLonRadI) * dAxI + cos(dLonRadI) * dAyI;
			double dAzLatCoeffI = - sin(dLatRadI) * (cos(dLonRadI) * dAxI + sin(dLonRadI) * dAyI) + cos(dLatRadI) * dAzI;

			double dAzLonCoeffK = - sin(dLonRadK) * dAxK + cos(dLonRadK) * dAyK;
			double dAzLatCoeffK = - sin(dLatRadK) * (cos(dLonRadK) * dAxK + sin(dLonRadK) * dAyK) + cos(dLatRadK) * dAzK;

			// Add contributions to sparse matrix operator
			opCurlE(i,i) += dScale * dCoeff0 * dAzLonCoeffI;
			opCurlN(i,i) += dScale * dCoeff0 * dAzLatCoeffI;

			opCurlE(i,k) += dScale * dCoeff1 * dAzLonCoeffK;
			opCurlN(i,k) += dScale * dCoeff1 * dAzLatCoeffK;
		}

		if (setPoints.size() < 4) {
			Announce("WARNING: Only %i points used for Curl in cell %i"
				" -- accuracy may be affected", setPoints.size(), i);
		}
	}

	kd_free(kdGrid);
}

///////////////////////////////////////////////////////////////////////////////

DataOp_CURL::DataOp_CURL(
	const std::string & strName,
	int nCurlPoints,
	double dCurlDist
) :
	DataOp(strName),
	m_nCurlPoints(nCurlPoints),
	m_dCurlDist(dCurlDist),
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

bool DataOp_CURL::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	if (!m_fInitialized) {
		Announce("Building Curl operator %s (%i, %1.2f)",
			m_strName.c_str(), m_nCurlPoints, m_dCurlDist);

		BuildCurlOperator(
			grid,
			m_nCurlPoints,
			m_dCurlDist,
			m_opCurlE,
			m_opCurlN);

		m_fInitialized = true;
	}

	m_opCurlE.Apply(*(vecArgData[0]), dataout, true);
	m_opCurlN.Apply(*(vecArgData[1]), dataout, false); 

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_DIVERGENCE
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build a sparse Div operator on an unstructured SimpleGrid.
///	</summary>
///	<param name="dDivDist">
///		Great circle radius of the Div operator, in degrees.
///	</param>
void BuildDivergenceOperator(
	const SimpleGrid & grid,
	int nDivPoints,
	double dDivDistDeg,
	SparseMatrix<float> & opDivE,
	SparseMatrix<float> & opDivN
) {
	opDivE.Clear();
	opDivN.Clear();

	int iRef = 0;

	// Scaling factor used in Div calculation
	const double dScale = 2.0 / (EarthRadius * static_cast<double>(nDivPoints) * DegToRad(dDivDistDeg));

	// TODO: Use SimpleGrid's built-in functionality
	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataArray1D<double> dXi(grid.GetSize());
	DataArray1D<double> dYi(grid.GetSize());
	DataArray1D<double> dZi(grid.GetSize());
	for (int i = 0; i < grid.GetSize(); i++) {
		double dLat = grid.m_dLat[i];
		double dLon = grid.m_dLon[i];

		dXi[i] = cos(dLon) * cos(dLat);
		dYi[i] = sin(dLon) * cos(dLat);
		dZi[i] = sin(dLat);

		kd_insert3(kdGrid, dXi[i], dYi[i], dZi[i], (void*)((&iRef)+i));
	}

	// Construct the Div operator using SPH
	for (int i = 0; i < grid.GetSize(); i++) {

		// Generate points for the Div
		std::vector<double> dXout;
		std::vector<double> dYout;
		std::vector<double> dZout;

		GenerateEqualDistanceSpherePoints(
			dXi[i], dYi[i], dZi[i],
			nDivPoints,
			dDivDistDeg,
			dXout, dYout, dZout);
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

		// Center point in spherical coordinates
		double dLonRadI, dLatRadI;
		XYZtoRLL_Rad(dXi[i], dYi[i], dZi[i], dLonRadI, dLatRadI);

		//printf("%1.5e %1.5e\n", RadToDeg(grid.m_dLon[i]), RadToDeg(grid.m_dLat[i]));
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

			// Local point in spherical coordinates
			double dLonRadK, dLatRadK;
			XYZtoRLL_Rad(dXi[k], dYi[k], dZi[k], dLonRadK, dLatRadK);

			// Calculate local radial vector from central lat/lon
			// i.e. project \vec{X}_k - \vec{X}_i to the surface of the
			//      sphere and normalize to unit length.
			double dRx = dXi[k] - dXi[i];
			double dRy = dYi[k] - dYi[i];
			double dRz = dZi[k] - dZi[i];

			double dChordLength = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

			if (dChordLength < 1.0e-12) {
				_EXCEPTIONT("Div sample point not distinct from center point: increase radius of operator.");
			}

			// Calculate radial unit vector at point I
			double dDot = dRx * dXi[i] + dRy * dYi[i] + dRz * dZi[i];

			double dRxI = dRx - dDot * dXi[i];
			double dRyI = dRy - dDot * dYi[i];
			double dRzI = dRz - dDot * dZi[i];

			double dMagI = sqrt(dRxI * dRxI + dRyI * dRyI + dRzI * dRzI);

			_ASSERT(dMagI >= 1.0e-12);

			dRxI /= dMagI;
			dRyI /= dMagI;
			dRzI /= dMagI;

			// Calculate radial unit vector at point K
			dDot = dRx * dXi[k] + dRy * dYi[k] + dRz * dZi[k];

			double dRxK = dRx - dDot * dXi[k];
			double dRyK = dRy - dDot * dYi[k];
			double dRzK = dRz - dDot * dZi[k];

			double dMagK = sqrt(dRxK * dRxK + dRyK * dRyK + dRzK * dRzK);

			_ASSERT(dMagK >= 1.0e-12);

			dRxK /= dMagK;
			dRyK /= dMagK;
			dRzK /= dMagK;

			// Great circle distance between points I and K
			double dRdeg =
				RadToDeg(GreatCircleDistanceFromChordLength_Rad(dChordLength));

			// Coefficients for linear interpolation / extrapolation to distance dDivDistDeg
			double dCoeff0 = (dRdeg - dDivDistDeg) / dRdeg;
			double dCoeff1 = dDivDistDeg / dRdeg;

			// Multiplicative coefficients for converting spherical velocities to azimuthal velocities
			double dRaLonCoeffI = - sin(dLonRadI) * dRxI + cos(dLonRadI) * dRyI;
			double dRaLatCoeffI = - sin(dLatRadI) * (cos(dLonRadI) * dRxI + sin(dLonRadI) * dRyI) + cos(dLatRadI) * dRzI;

			double dRaLonCoeffK = - sin(dLonRadK) * dRxK + cos(dLonRadK) * dRyK;
			double dRaLatCoeffK = - sin(dLatRadK) * (cos(dLonRadK) * dRxK + sin(dLonRadK) * dRyK) + cos(dLatRadK) * dRzK;

			// Add contributions to sparse matrix operator
			opDivE(i,i) += dScale * dCoeff0 * dRaLonCoeffI;
			opDivN(i,i) += dScale * dCoeff0 * dRaLatCoeffI;

			opDivE(i,k) += dScale * dCoeff1 * dRaLonCoeffK;
			opDivN(i,k) += dScale * dCoeff1 * dRaLatCoeffK;
		}

		if (setPoints.size() < 4) {
			Announce("WARNING: Only %i points used for Div in cell %i"
				" -- accuracy may be affected", setPoints.size(), i);
		}
	}
/*
	for (
		SparseMatrix<double>::SparseMapIterator iter = opDiv.begin();
		iter != opDiv.end(); iter++
	) {
		std::cout << iter->first.first << " " << iter->first.second << " " << iter->second << std::endl;
	}
*/
	kd_free(kdGrid);
}

///////////////////////////////////////////////////////////////////////////////

DataOp_DIVERGENCE::DataOp_DIVERGENCE(
	const std::string & strName,
	int nDivPoints,
	double dDivDist
) :
	DataOp(strName),
	m_nDivPoints(nDivPoints),
	m_dDivDist(dDivDist),
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

bool DataOp_DIVERGENCE::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 2) {
		_EXCEPTION2("%s expects two arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	if (!m_fInitialized) {
		Announce("Building Div operator %s (%i, %1.2f)",
			m_strName.c_str(), m_nDivPoints, m_dDivDist);

		BuildDivergenceOperator(
			grid,
			m_nDivPoints,
			m_dDivDist,
			m_opDivE,
			m_opDivN);

		m_fInitialized = true;
	}

	m_opDivE.Apply(*(vecArgData[0]), dataout, true);
	m_opDivN.Apply(*(vecArgData[1]), dataout, false); 

	return true;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build a sparse gradient magnitude operator on an unstructured SimpleGrid.
///	</summary>
///	<param name="dGradMagDistDeg">
///		Great circle radius of the gradient magnitude operator, in degrees.
///	</param>
void BuildGradOperator(
	const SimpleGrid & grid,
	int nGradPoints,
	double dGradDistDeg,
	SparseMatrix< DataOp_GRADMAG::pair_with_plus_minus<float> > & opGrad
) {
	opGrad.Clear();

	int iRef = 0;

	// Scaling factor used in Laplacian calculation
	const double dScale = 2.0 / (EarthRadius * static_cast<double>(nGradPoints) * DegToRad(dGradDistDeg));

	// TODO: Use SimpleGrid's built-in functionality
	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataArray1D<double> dXi(grid.GetSize());
	DataArray1D<double> dYi(grid.GetSize());
	DataArray1D<double> dZi(grid.GetSize());
	for (int i = 0; i < grid.GetSize(); i++) {
		double dLat = grid.m_dLat[i];
		double dLon = grid.m_dLon[i];

		dXi[i] = cos(dLon) * cos(dLat);
		dYi[i] = sin(dLon) * cos(dLat);
		dZi[i] = sin(dLat);

		kd_insert3(kdGrid, dXi[i], dYi[i], dZi[i], (void*)((&iRef)+i));
	}

	// Construct the gradient magnitude operator
	for (int i = 0; i < grid.GetSize(); i++) {

		// Generate points for the Laplacian
		std::vector<double> dXout;
		std::vector<double> dYout;
		std::vector<double> dZout;

		GenerateEqualDistanceSpherePoints(
			dXi[i], dYi[i], dZi[i],
			nGradPoints,
			dGradDistDeg,
			dXout, dYout, dZout);
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

			// Calculate local radial vector from central lat/lon
			// i.e. project \vec{X}_k - \vec{X}_i to the surface of the
			//      sphere and normalize to unit length.
			double dRx = dXi[k] - dXi[i];
			double dRy = dYi[k] - dYi[i];
			double dRz = dZi[k] - dZi[i];

			double dChordLength = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

			if (dChordLength < 1.0e-12) {
				_EXCEPTIONT("Grad sample point not distinct from center point: increase radius of operator.");
			}

			// Generate perpendicular vector fields in Cartesian coords
			double dEx, dEy, dEz;
			double dNx, dNy, dNz;

			// Pick eastward and northward reference directions
			if (fabs(fabs(dZi[i]) - 1.0) > 1.0e-10) {
				double dEmag = sqrt(dXi[i] * dXi[i] + dYi[i] * dYi[i]);

				dEx = - dYi[i] / dEmag;
				dEy = dXi[i] / dEmag;
				dEz = 0.0;

				dNx = - dEy * dZi[i];
				dNy = dEx * dZi[i];
				dNz = dEmag;

			// At poles point
			} else {
				dEx = 1.0;
				dEy = 0.0;
				dEz = 0.0;

				dNx = 0.0;
				dNy = 1.0;
				dNz = 0.0;
			}

			// Calculate radial unit vector at point K
			double dDot = dRx * dXi[k] + dRy * dYi[k] + dRz * dZi[k];

			double dRxK = dRx - dDot * dXi[k];
			double dRyK = dRy - dDot * dYi[k];
			double dRzK = dRz - dDot * dZi[k];

			double dMagK = sqrt(dRxK * dRxK + dRyK * dRyK + dRzK * dRzK);

			_ASSERT(dMagK >= 1.0e-12);

			dRxK /= dMagK;
			dRyK /= dMagK;
			dRzK /= dMagK;

			// Coefficients of gradient operators
			double dEDotRK = dRxK * dEx + dRyK * dEy + dRzK * dEz;
			double dNDotRK = dRxK * dNx + dRyK * dNy + dRzK * dNz;

			opGrad(i,k) += DataOp_GRADMAG::pair_with_plus_minus<float>(dScale * dEDotRK, dScale * dNDotRK);
			opGrad(i,i) -= DataOp_GRADMAG::pair_with_plus_minus<float>(dScale * dEDotRK, dScale * dNDotRK);
		}

		if (setPoints.size() < 4) {
			Announce("WARNING: Only %i points used for gradient in cell %i"
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

DataOp_GRADMAG::DataOp_GRADMAG(
	const std::string & strName,
	int nGradPoints,
	double dGradDist
) :
	DataOp(strName),
	m_nGradPoints(nGradPoints),
	m_dGradDist(dGradDist),
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

bool DataOp_GRADMAG::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
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
		Announce("Building gradient operator %s (%i, %1.2f)",
			m_strName.c_str(), m_nGradPoints, m_dGradDist);

		BuildGradOperator(
			grid,
			m_nGradPoints,
			m_dGradDist,
			m_opGrad);

		m_fInitialized = true;
	}

	DataArray1D<float> const & data = *(vecArgData[0]);
	DataArray1D<float> datatemp(dataout.GetRows());
	dataout.Zero();
	for (auto it = m_opGrad.begin(); it != m_opGrad.end(); it++) {
		int row = it->first.first;
		int col = it->first.second;
		dataout[row] += it->second.first * data[col];
		datatemp[row] += it->second.second * data[col];
	}
	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = sqrt(dataout[i] * dataout[i] + datatemp[i] * datatemp[i]);
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

DataOp_VECDOTGRAD::DataOp_VECDOTGRAD(
	const std::string & strName,
	int nGradPoints,
	double dGradDist
) :
	DataOp(strName),
	m_nGradPoints(nGradPoints),
	m_dGradDist(dGradDist),
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

bool DataOp_VECDOTGRAD::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 3) {
		_EXCEPTION2("%s expects three arguments: %i given",
			m_strName.c_str(), strArg.size());
	}
	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL) || (vecArgData[2] == NULL)) {
		_EXCEPTION1("Arguments to %s must be data variables",
			m_strName.c_str());
	}

	if (!m_fInitialized) {
		Announce("Building gradient operator %s (%i, %1.2f)",
			m_strName.c_str(), m_nGradPoints, m_dGradDist);

		BuildGradOperator(
			grid,
			m_nGradPoints,
			m_dGradDist,
			m_opGrad);

		m_fInitialized = true;
	}

	DataArray1D<float> const & dataU = *(vecArgData[0]);
	DataArray1D<float> const & dataV = *(vecArgData[1]);

	DataArray1D<float> const & data = *(vecArgData[2]);
	DataArray1D<float> datatemp(dataout.GetRows());
	dataout.Zero();
	for (auto it = m_opGrad.begin(); it != m_opGrad.end(); it++) {
		int row = it->first.first;
		int col = it->first.second;
		dataout[row] += it->second.first * data[col];
		datatemp[row] += it->second.second * data[col];
	}
	for (int i = 0; i < dataout.GetRows(); i++) {
		dataout[i] = dataU[i] * dataout[i] + dataV[i] * datatemp[i];
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_DIVERGENCE
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build a sparse Mean operator on an unstructured SimpleGrid.
///	</summary>
///	<param name="dMeanDist">
///		Great circle radius of the Mean operator, in degrees.
///	</param>
void BuildMeanOperator(
	const SimpleGrid & grid,
	double dMeanDistDeg,
	SparseMatrix<float> & opMean
) {
	opMean.Clear();

	int iRef = 0;

	if (dMeanDistDeg <= 0.0) {
		_EXCEPTION1("dist (%1.2f) in _MEAN{dist} must be positive", dMeanDistDeg);
	}
	if (grid.m_dArea.GetRows() != grid.GetSize()) {
		_EXCEPTIONT("SimpleGrid must have precomputed areas to call _MEAN{}");
	}

	// Convert great circle dist to chord dist
	double dMeanDistXYZ = ChordLengthFromGreatCircleDistance_Deg(dMeanDistDeg);

	// TODO: Use SimpleGrid's built-in functionality
	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataArray1D<double> dXi(grid.GetSize());
	DataArray1D<double> dYi(grid.GetSize());
	DataArray1D<double> dZi(grid.GetSize());
	for (int i = 0; i < grid.GetSize(); i++) {
		double dLat = grid.m_dLat[i];
		double dLon = grid.m_dLon[i];

		dXi[i] = cos(dLon) * cos(dLat);
		dYi[i] = sin(dLon) * cos(dLat);
		dZi[i] = sin(dLat);

		kd_insert3(kdGrid, dXi[i], dYi[i], dZi[i], (void*)((&iRef)+i));
	}

	// Area of accumulated nodes
	std::vector<int> vecNodeIx;
	double dAccumulatedArea;
	std::vector<double> vecNodeArea;

	// Construct the Mean operator
	for (int i = 0; i < grid.GetSize(); i++) {

		vecNodeIx.clear();
		vecNodeArea.clear();
		dAccumulatedArea = 0.0;

		// Query kd-tree
		kdres * kdr = kd_nearest_range3(kdGrid, dXi[i], dYi[i], dZi[i], dMeanDistXYZ);
		if (kdr == NULL) {
			_EXCEPTIONT("Error in kd_nearest_range3");
		}
		int nNodes = kd_res_size(kdr);

		//if (i % 1000 == 0) {
		//	std::cout << dMeanDistXYZ << " : " << nNodes << std::endl;
		//}

		// Loop through kd-tree
		while (!kd_res_end(kdr)) {

			// Get data index
			void* pData = kd_res_item_data(kdr);
			if (pData == NULL) {
				_EXCEPTIONT("NULL data index");
			}
			int k = ((int*)(pData)) - (&iRef);
			_ASSERT((k >= 0) && (k < grid.GetSize()));

			// Insert this node
			vecNodeIx.push_back(k);
			vecNodeArea.push_back(grid.m_dArea[k]);
			dAccumulatedArea += grid.m_dArea[k];

			kd_res_next(kdr);
		}

		_ASSERT(dAccumulatedArea > 0.0);

		// Delete the kdres
		kd_res_free(kdr);

		// Insert new row into sparse matrix
		for (int k = 0; k < vecNodeIx.size(); k++) {
			opMean(i,vecNodeIx[k]) = static_cast<float>(vecNodeArea[k] / dAccumulatedArea);
		}
	}

	// Cleanup kd-tree
	kd_free(kdGrid);
}

///////////////////////////////////////////////////////////////////////////////

DataOp_MEAN::DataOp_MEAN(
	const std::string & strName,
	double dMeanDist
) :
	DataOp(strName),
	m_dMeanDist(dMeanDist),
	m_fInitialized(false)
{ }

///////////////////////////////////////////////////////////////////////////////

bool DataOp_MEAN::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 1) {
		_EXCEPTION2("%s expects one argument: %i given",
			m_strName.c_str(), strArg.size());
	}
	if (vecArgData[0] == NULL) {
		_EXCEPTION1("Argument to %s must be data variable",
			m_strName.c_str());
	}

	if (!m_fInitialized) {
		Announce("Building MEAN operator %s (%1.2f)",
			m_strName.c_str(), m_dMeanDist);

		BuildMeanOperator(
			grid,
			m_dMeanDist,
			m_opMean);

		m_fInitialized = true;
	}

	m_opMean.Apply(*(vecArgData[0]), dataout, true);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_DAILYCHILLHOURS
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_DAILYCHILLHOURS::name = "_DAILYCHILLHOURS";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_DAILYCHILLHOURS::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 3) {
		_EXCEPTION2("%s expects three arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("%s expects first argument (tasmin) and second argument (tasmax) to be data variables",
			m_strName.c_str());
	}
	const DataArray1D<float> & dataTasmin = *(vecArgData[0]);
	const DataArray1D<float> & dataTasmax = *(vecArgData[1]);

	// Freezing and chilling temperature are defined in degrees F
	float dFreezingTemp = 32.0;
	float dChillingTemp = 45.0;

	// Convert to local unit
	bool fSuccess;
	fSuccess = ConvertUnits<float>(dFreezingTemp, "degF", strArg[2], false);
	if (!fSuccess) {
		_EXCEPTION2("%s cannot convert freezing temperature to provided units \"%s\"",
			m_strName.c_str(), strArg[2].c_str());
	}

	fSuccess = ConvertUnits<float>(dChillingTemp, "degF", strArg[2], false);
	if (!fSuccess) {
		_EXCEPTION2("%s cannot convert chilling temperature to provided units \"%s\"",
			m_strName.c_str(), strArg[2].c_str());
	}

	// Calculate daily chill hours using similar triangles
	bool fWarning = false;
	for (int i = 0; i < dataout.GetRows(); i++) {
		if (dataTasmin[i] > dataTasmax[i]) {
			dataout[i] = 0.0;
			fWarning = true;
		} else if (dataTasmax[i] < dFreezingTemp) {
			dataout[i] = 0.0;
		} else if (dataTasmax[i] < dChillingTemp) {
			if (dataTasmin[i] < dFreezingTemp) {
				dataout[i] = 24.0 * (dataTasmax[i] - dFreezingTemp) / (dataTasmax[i] - dataTasmin[i]);
			} else {
				dataout[i] = 24.0;
			}
		} else {
			if (dataTasmin[i] < dFreezingTemp) {
				dataout[i] = 24.0 * (dChillingTemp - dFreezingTemp) / (dataTasmax[i] - dataTasmin[i]);
			} else if (dataTasmin[i] < dChillingTemp) {
				dataout[i] = 24.0 * (dChillingTemp - dataTasmin[i]) / (dataTasmax[i] - dataTasmin[i]);
			} else {
				dataout[i] = 0.0;
			}
		}
	}

	if (fWarning) {
		Announce("WARNING: Some grid points have tasmax < tasmin");
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_RELHUMFROMTDTA
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_RELHUMFROMTDTA::name = "_RELHUMFROMTDTA";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_RELHUMFROMTDTA::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 3) {
		_EXCEPTION2("%s expects three arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("%s expects first argument (td) and second argument (ta) to be data variables",
			m_strName.c_str());
	}
	const DataArray1D<float> & dataTd = *(vecArgData[0]);
	const DataArray1D<float> & dataTa = *(vecArgData[1]);

	dataout.SetFillValue(DefaultFillValue);
	dataout.SetUnits("percent");

	// Calculate relative humidity from Td and Ta

	// Calculation with variables provided in degrees Celsius
	if ((strArg[2] == "degC") || (strArg[2] == "C")) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataTd.HasFillValue() && (dataTd[i] == dataTd.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = 100.0 * exp((17.1 * dataTd[i]) / (235.0 + dataTd[i]) - (17.1 * dataTa[i]) / (235.0 + dataTa[i]));
			}
		}

	// Calculation with variables provided in Kelvin
	} else if (strArg[2] == "K") {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataTd.HasFillValue() && (dataTd[i] == dataTd.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = 100.0 * exp((17.1 * (dataTd[i] - 273.15)) / (dataTd[i] - 38.15) - (17.1 * (dataTa[i] - 273.15)) / (dataTa[i] - 38.15));
			}
		}

	// Calculation with variables provided in degrees Fahrenheit
	} else if ((strArg[2] == "degF") || (strArg[2] == "F")) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataTd.HasFillValue() && (dataTd[i] == dataTd.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				double dTdDegC = (5.0/9.0) * (dataTd[i] - 32.0);
				double dTaDegC = (5.0/9.0) * (dataTa[i] - 32.0);

				dataout[i] = 100.0 * exp((17.1 * dTdDegC) / (235.0 + dTdDegC) - (17.1 * dTaDegC) / (235.0 + dTaDegC));
			}
		}
	
	// Invalid unit
	} else {
		_EXCEPTION1("Invalid third argument units in _RELHUMFROMTDT2M (%s): Expected \"degC\", \"degF\", \"K\"",
			strArg[2].c_str());
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_VPDFROMTAHUR
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_VPDFROMTAHUR::name = "_VPDFROMTAHUR";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_VPDFROMTAHUR::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 3) {
		_EXCEPTION2("%s expects three arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("%s expects first argument (ta) and second argument (hur) to be data variables",
			m_strName.c_str());
	}
	const DataArray1D<float> & dataTa = *(vecArgData[0]);
	const DataArray1D<float> & dataHur = *(vecArgData[1]);

	dataout.SetFillValue(DefaultFillValue);
	dataout.SetUnits("hPa");

	// Calculation with variables provided in degrees Celsius
	if ((strArg[2] == "degC") || (strArg[2] == "C")) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataHur.HasFillValue() && (dataHur[i] == dataHur.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = SatVapPres_FromCC_degC(dataTa[i]) * (1.0 - dataHur[i] / 100.0);
			}
		}

	// Calculation with variables provided in Kelvin
	} else if (strArg[2] == "K") {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataHur.HasFillValue() && (dataHur[i] == dataHur.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = SatVapPres_FromCC_K(dataTa[i]) * (1.0 - dataHur[i] / 100.0);
			}
		}

	// Calculation with variables provided in degrees Fahrenheit
	} else if ((strArg[2] == "degF") || (strArg[2] == "F")) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataHur.HasFillValue() && (dataHur[i] == dataHur.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				dataout[i] = SatVapPres_FromCC_degF(dataTa[i]) * (1.0 - dataHur[i] / 100.0);
			}
		}

	// Invalid unit
	} else {
		_EXCEPTION1("Invalid third argument units in _VPDFROMTARH (%s): Expected \"degC\", \"degF\", \"K\"",
			strArg[2].c_str());
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// DataOp_VPDFROMTAHUSPRES
///////////////////////////////////////////////////////////////////////////////

const char * DataOp_VPDFROMTAHUSPRES::name = "_VPDFROMTAHUSPRES";

///////////////////////////////////////////////////////////////////////////////

bool DataOp_VPDFROMTAHUSPRES::Apply(
	const SimpleGrid & grid,
	const std::vector<std::string> & strArg,
	const std::vector<DataArray1D<float> const *> & vecArgData,
	DataArray1D<float> & dataout
) {
	if (strArg.size() != 4) {
		_EXCEPTION2("%s expects four arguments: %i given",
			m_strName.c_str(), strArg.size());
	}

	if ((vecArgData[0] == NULL) || (vecArgData[1] == NULL)) {
		_EXCEPTION1("%s expects first argument (ta), second argument (hus) and third argument (pr) to be data variables",
			m_strName.c_str());
	}
	const DataArray1D<float> & dataTa = *(vecArgData[0]);
	const DataArray1D<float> & dataHus = *(vecArgData[1]);
	const DataArray1D<float> & dataPr = *(vecArgData[2]);

	dataout.SetFillValue(DefaultFillValue);
	dataout.SetUnits("hPa");

	// Calculation with variables provided in degrees Celsius
	if ((strArg[3] == "degC") || (strArg[3] == "C")) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataHus.HasFillValue() && (dataHus[i] == dataHus.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue())) ||
			    (dataPr.HasFillValue() && (dataPr[i] == dataPr.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				// Actual vapor pressure = q * p / ((Mw/Md) + (1-Mw/Md) * q)
				// Assumes pressure is given in Pa
				double dAVP = dataHus[i] * (dataPr[i] / 100.0) / (0.6221 + 0.3779 * dataHus[i]);
				dataout[i] = SatVapPres_FromCC_degC(dataTa[i]) - dAVP;
			}
		}

	// Calculation with variables provided in Kelvin
	} else if (strArg[3] == "K") {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataHus.HasFillValue() && (dataHus[i] == dataHus.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue())) ||
			    (dataPr.HasFillValue() && (dataPr[i] == dataPr.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				double dAVP = dataHus[i] * (dataPr[i] / 100.0) / (0.6221 + 0.3779 * dataHus[i]);
				dataout[i] = SatVapPres_FromCC_K(dataTa[i]) - dAVP;
			}
		}

	// Calculation with variables provided in degrees Fahrenheit
	} else if ((strArg[3] == "degF") || (strArg[3] == "F")) {
		for (int i = 0; i < dataout.GetRows(); i++) {
			if ((dataHus.HasFillValue() && (dataHus[i] == dataHus.GetFillValue())) ||
			    (dataTa.HasFillValue() && (dataTa[i] == dataTa.GetFillValue())) ||
			    (dataPr.HasFillValue() && (dataPr[i] == dataPr.GetFillValue()))
			) {
				dataout[i] = dataout.GetFillValue();
			} else {
				double dAVP = dataHus[i] * (dataPr[i] / 100.0) / (0.6221 + 0.3779 * dataHus[i]) / (1000.0);
				dataout[i] = SatVapPres_FromCC_degF(dataTa[i]) - dAVP;
			}
		}

	// Invalid unit
	} else {
		_EXCEPTION1("Invalid third argument units in _VPDFROMTAHUSPRES (%s): Expected \"degC\", \"degF\", \"K\"",
			strArg[2].c_str());
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

