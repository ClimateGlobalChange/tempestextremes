///////////////////////////////////////////////////////////////////////////////
///
///	\file    CoordTransforms.h
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
///		Convert radians to degrees.
///	</summary>
inline double RadToDeg(
	double dRad
) {
	return (dRad * 180.0 / M_PI);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert degrees to radians.
///	</summary>
inline double DegToRad(
	double dDeg
) {
	return (dDeg * M_PI / 180.0);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Translate a longitude value to the range [0,360)
///	</summary>
inline double LonDegToStandardRange(
	double dLonDeg
) {
	dLonDeg = (dLonDeg - 360.0 * floor(dLonDeg / 360.0));
	if ((dLonDeg < 0.0) || (dLonDeg >= 360.0)) {
		return 0.0;
	}
	return dLonDeg;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Translate a longitude value to the range [0,2pi)
///	</summary>
inline double LonRadToStandardRange(
	double dLonRad
) {
	dLonRad = (dLonRad - (2.0 * M_PI) * floor(dLonRad / (2.0 * M_PI)));
	if ((dLonRad < 0.0) || (dLonRad >= (2.0 * M_PI))) {
		return 0.0;
	}
	return dLonRad;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the chord length from the great circle distance (in degrees).
///	</summary>
inline double ChordLengthFromGreatCircleDistance_Deg(
	double dGCDDeg
) {
	return 2.0 * sin(0.5 * DegToRad(dGCDDeg));
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the great circle distance (in radians) from the chord length.
///	</summary>
inline double GreatCircleDistanceFromChordLength_Rad(
	double dChordLength
) {
	double dHalfDist = 0.5 * dChordLength;
	if (dHalfDist > 1.0) {
		_ASSERT(dHalfDist < 1.0 + 1.0e-8);
		dHalfDist = 1.0;
	}

	return 2.0 * asin(dHalfDist);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate 3D Cartesian coordinates from latitude and longitude,
///		in radians.
///	</summary>
inline void RLLtoXYZ_Rad(
	double dLonRad,
	double dLatRad,
	double & dX,
	double & dY,
	double & dZ
) {
	if (fabs(dLatRad) > 0.5 * M_PI + HighTolerance) {
		_EXCEPTION1("Latitude out of range (%2.14f)", dLatRad);
	}

	dX = cos(dLonRad) * cos(dLatRad);
	dY = sin(dLonRad) * cos(dLatRad);
	dZ = sin(dLatRad);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate 3D Cartesian coordinates from latitude and longitude,
///		in degrees.
///	</summary>
inline void RLLtoXYZ_Deg(
	double dLonDeg,
	double dLatDeg,
	double & dX,
	double & dY,
	double & dZ
) {
	return RLLtoXYZ_Rad(
		DegToRad(dLonDeg),
		DegToRad(dLatDeg),
		dX, dY, dZ);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate latitude and longitude from normalized 3D Cartesian
///		coordinates, in degrees.
///	</summary>
inline void XYZtoRLL_Deg(
	double dX,
	double dY,
	double dZ,
	double & dLonDeg,
	double & dLatDeg
) {
	double dMag2 = dX * dX + dY * dY + dZ * dZ;

	if (fabs(dMag2 - 1.0) >= 0.01) {
		_EXCEPTION4("Grid point has non-unit magnitude: "
			"(%1.15e, %1.15e, %1.15e) (magnitude %1.15e)",
			dX, dY, dZ, fabs(dX * dX + dY * dY + dZ * dZ));

	}

	double dMag = sqrt(dMag2);

	dX /= dMag;
	dY /= dMag;
	dZ /= dMag;

	if (fabs(dZ) < 1.0 - ReferenceTolerance) {
		dLonDeg = RadToDeg(atan2(dY, dX));
		dLatDeg = RadToDeg(asin(dZ));

		if (dLonDeg < 0.0) {
			dLonDeg += 360.0;
		}

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
///		Calculate latitude and longitude from normalized 3D Cartesian
///		coordinates, in radians.
///	</summary>
inline void XYZtoRLL_Rad(
	double dX,
	double dY,
	double dZ,
	double & dLonRad,
	double & dLatRad
) {
	double dMag2 = dX * dX + dY * dY + dZ * dZ;

	if (fabs(dMag2 - 1.0) >= 0.01) {
		_EXCEPTION4("Grid point has non-unit magnitude: "
			"(%1.15e, %1.15e, %1.15e) (magnitude %1.15e)",
			dX, dY, dZ, fabs(dX * dX + dY * dY + dZ * dZ));

	}

	double dMag = sqrt(dMag2);

	dX /= dMag;
	dY /= dMag;
	dZ /= dMag;

	if (fabs(dZ) < 1.0 - ReferenceTolerance) {
		dLonRad = atan2(dY, dX);
		dLatRad = asin(dZ);

		if (dLonRad < 0.0) {
			dLonRad += 2.0 * M_PI;
		}

	} else if (dZ > 0.0) {
		dLonRad = 0.0;
		dLatRad = 0.5 * M_PI;

	} else {
		dLonRad = 0.0;
		dLatRad = -0.5 * M_PI;
	}

}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a vector in RLL coordinates to XYZ.
///	</summary>
inline void VecTransRLL2DtoXYZ_Rad(
	double dLonRad,
	double dLatRad,
	double dUlon,
	double dUlat,
	double & dUx,
	double & dUy,
	double & dUz
) {
	double dSinLon = sin(dLonRad);
	double dCosLon = cos(dLonRad);
	double dSinLat = sin(dLatRad);
	double dCosLat = cos(dLatRad);

	dUx = - dSinLon * dUlon - dCosLon * dSinLat * dUlat;
	dUy =   dCosLon * dUlon - dSinLon * dSinLat * dUlat;
	dUz =                               dCosLat * dUlat;
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
///		Calculate the great circle distance between two points on the sphere,
///		with return value in radians.
///	</summary>
inline double GreatCircleDistanceXYZ_Rad(
	double dX0,
	double dY0,
	double dZ0,
	double dX1,
	double dY1,
	double dZ1
) {
	double dDX = dX1 - dX0;
	double dDY = dY1 - dY0;
	double dDZ = dZ1 - dZ0;

	double dChordLength = sqrt(dDX * dDX + dDY * dDY + dDZ * dDZ);

	return GreatCircleDistanceFromChordLength_Rad(dChordLength);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the great circle distance between two points on the sphere,
///		with return value in degrees.
///	</summary>
inline double GreatCircleDistanceXYZ_Deg(
	double dX0,
	double dY0,
	double dZ0,
	double dX1,
	double dY1,
	double dZ1
) {
	return RadToDeg(GreatCircleDistanceXYZ_Rad(dX0, dY0, dZ0, dX1, dY1, dZ1));
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
	return
		RadToDeg(
			GreatCircleDistance_Rad(
				dLonRad1,
				dLatRad1,
				dLonRad2,
				dLatRad2));
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the coordinates of a point (dLonRad, dLatRad) on the
///		stereographic projection centered at (dLonRad0, dLatRad0).  The Xs
///		coordinate is to the east and the Ys coordinate is to the north.
///		At the poles the Y coordinate is a continuation of the constant
///		longitude line dLonRad0 and the X coordinate points to its right.
///	</summary>
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

///	<summary>
///		Calculate the inverse of StereographicProjection().
///	</summary>
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

///	<summary>
///		Calculate the spherical direction between two RLL points, specified
///		in radians.  The result is returned in radians.
///	</summary>
inline void GreatCircleDirection_Rad(
	double dLon1Rad,
	double dLat1Rad,
	double dLon2Rad,
	double dLat2Rad,
	double & dZonalDirRad,
	double & dMeridDirRad
) {
	if ((fabs(dLon1Rad - dLon2Rad) < 1.0e-14) && (fabs(dLat1Rad - dLat2Rad) < 1.0e-14)) {
		dZonalDirRad = 0.0;
		dMeridDirRad = 0.0;
		return;
	}

	// Use consistent stereographic projection at poles
	if ((dLat1Rad > 0.5 * M_PI - 1.0e-7) || (dLat1Rad < -0.5 * M_PI + 1.0e-7))  {
		dLon1Rad = 0.0;
	}
	if ((dLat2Rad > 0.5 * M_PI - 1.0e-7) || (dLat2Rad < -0.5 * M_PI + 1.0e-7))  {
		dLon2Rad = 0.0;
	}

	double dXs;
	double dYs;
	StereographicProjection(
		dLon1Rad,
		dLat1Rad,
		dLon2Rad,
		dLat2Rad,
		dXs,
		dYs);

	double dDeltaX = sqrt(dXs * dXs + dYs * dYs);

	double dGreatCircleDistRad = atan(dDeltaX);

	if (fabs(dDeltaX) < 1.0e-12) {
		dZonalDirRad = 0.0;
		dMeridDirRad = 0.0;
		return;
	}

	dZonalDirRad = dXs * dGreatCircleDistRad / dDeltaX;
	dMeridDirRad = dYs * dGreatCircleDistRad / dDeltaX;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the spherical direction between two RLL points, specified
///		in radians. The result is returned in degrees.
///	</summary>
inline void GreatCircleDirection_Deg(
	double dLon1Rad,
	double dLat1Rad,
	double dLon2Rad,
	double dLat2Rad,
	double & dZonalDirDeg,
	double & dMeridDirDeg
) {
	double dZonalDirRad;
	double dMeridDirRad;

	GreatCircleDirection_Rad(
		dLon1Rad,
		dLat1Rad,
		dLon2Rad,
		dLat2Rad,
		dZonalDirRad,
		dMeridDirRad);

	dZonalDirDeg = RadToDeg(dZonalDirRad);
	dMeridDirDeg = RadToDeg(dMeridDirRad);
}

///////////////////////////////////////////////////////////////////////////////

#endif // _COORDTRANSFORMS_H_

