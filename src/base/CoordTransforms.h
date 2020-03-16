///////////////////////////////////////////////////////////////////////////////
///
///	\file    CoordElements.h
///	\author  Paul Ullrich
///	\version March 14, 2020
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _COORDTRANSFORMS_H_
#define _COORDTRANSFORMS_H_

///////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "Defines.h"
#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate latitude and longitude from normalized 3D Cartesian
///		coordinates, in degrees.
///	</summary>
inline void XYZtoRLL_Deg(
	const double & dX,
	const double & dY,
	const double & dZ,
	double & dLonDeg,
	double & dLatDeg
) {
	_ASSERT(fabs(dX * dX + dY * dY + dZ * dZ - 1.0) < HighTolerance);

	if (fabs(dZ) < 1.0 - ReferenceTolerance) {
		dLonDeg = atan2(dY, dX);
		dLatDeg = asin(dZ);

		if (dLonDeg < 0.0) {
			dLonDeg += 2.0 * M_PI;
		}

		dLonDeg = dLonDeg / M_PI * 180.0;
		dLatDeg = dLatDeg / M_PI * 180.0;

	} else if (dZ > 0.0) {
		dLonDeg = 0.0;
		dLatDeg = 90.0;

	} else {
		dLonDeg = 0.0;
		dLatDeg = -90.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate an average longitude from two given longitudes (in radians).
///	</summary>
inline double AverageLongitude_Rad(
	double dLonRad1,
	double dLonRad2
) {
	double dDeltaLonRad;
	if (dLonRad2 > dLonRad1) {
		dDeltaLonRad = fmod(dLonRad2 - dLonRad1, 2.0 * M_PI);
		if (dDeltaLonRad > M_PI) {
			dDeltaLonRad = dDeltaLonRad - 2.0 * M_PI;
		}
	} else {
		dDeltaLonRad = - fmod(dLonRad1 - dLonRad2, 2.0 * M_PI);
		if (dDeltaLonRad < -M_PI) {
			dDeltaLonRad = dDeltaLonRad + 2.0 * M_PI;
		}
	}

	double dLonRadAvg = dLonRad1 + 0.5 * dDeltaLonRad;

	if ((dLonRadAvg < 0.0) && (dLonRad1 >= 0.0) && (dLonRad2 >= 0.0)) {
		dLonRadAvg += 2.0 * M_PI;
	}
	if ((dLonRadAvg > 2.0 * M_PI) && (dLonRad1 <= 2.0 * M_PI) && (dLonRad2 <= 2.0 * M_PI)) {
		dLonRadAvg -= 2.0 * M_PI;
	}

	return dLonRadAvg;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the great circle distance (in radians) between two RLL points.
///	</summary>
inline double GreatCircleDistance_Rad(
	double dLonRad1,
	double dLatRad1,
	double dLonRad2,
	double dLatRad2
) {
	double dR =
		sin(dLatRad1) * sin(dLatRad2)
		+ cos(dLatRad1) * cos(dLatRad2) * cos(dLonRad2 - dLonRad1);

	if (dR >= 1.0) {
		dR = 0.0;
	} else if (dR <= -1.0) {
		dR = M_PI;
	} else {
		dR = acos(dR);
	}
	if (dR != dR) {
		_EXCEPTIONT("NaN value detected");
	}

	return dR;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the great circle distance (in degrees) between two RLL points.
///	</summary>
inline double GreatCircleDistance_Deg(
	double dLonRad1,
	double dLatRad1,
	double dLonRad2,
	double dLatRad2
) {
	return (180.0 / M_PI) *
		GreatCircleDistance_Rad(dLonRad1, dLatRad1, dLonRad2, dLatRad2);
}

///////////////////////////////////////////////////////////////////////////////

inline void StereographicProjection(
	double dLonRad0,
	double dLatRad0,
	double dLonRad,
	double dLatRad,
	double & dXs,
	double & dYs
) {
	// Forward projection using equations (1)-(3)
	// http://mathworld.wolfram.com/StereographicProjection.html
	double dK = 2.0 / (1.0 + sin(dLatRad0) * sin(dLatRad) + cos(dLatRad0) * cos(dLatRad) * cos(dLonRad - dLonRad0));
	dXs = dK * cos(dLatRad) * sin(dLonRad - dLonRad0);
	dYs = dK * (cos(dLatRad0) * sin(dLatRad) - sin(dLatRad0) * cos(dLatRad) * cos(dLonRad - dLonRad0));
}

///////////////////////////////////////////////////////////////////////////////

inline void StereographicProjectionInv(
	double dLonRad0,
	double dLatRad0,
	double dXs,
	double dYs,
	double & dLonRad,
	double & dLatRad
) {
	// Forward projection using equations (3)-(5)
	// http://mathworld.wolfram.com/StereographicProjection.html
	double dRho = sqrt(dXs * dXs + dYs * dYs);
	double dC = 2.0 * atan(0.5 * dRho);

	if (dRho < 1.0e-14) {
		dLatRad = dLatRad0;
		dLonRad = dLonRad0;
		return;
	}

	dLatRad = asin(cos(dC) * sin(dLatRad0) + dYs * sin(dC) * cos(dLatRad0) / dRho);
	dLonRad = dLonRad0 + atan2(
			dXs * sin(dC),
			dRho * cos(dLatRad0) * cos(dC) - dYs * sin(dLatRad0) * sin(dC));
}

///////////////////////////////////////////////////////////////////////////////

#endif // _COORDTRANSFORMS_H_

