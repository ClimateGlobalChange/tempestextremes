///////////////////////////////////////////////////////////////////////////////
///
///	\file    Units.h
///	\author  Paul Ullrich
///	\version August 14, 2020
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

#ifndef _UNITS_H_
#define _UNITS_H_

#include <string>
#include <cmath>

#include "Constants.h"

///////////////////////////////////////////////////////////////////////////////

template <typename T>
bool ConvertUnits(
	T & dValue,
	const std::string & strUnit,
	const std::string & strTargetUnit,
	bool fIsDelta = false
) {
	// Unit is equal to TargetUnit
	if (strUnit == strTargetUnit) {

	// Perform unit conversion from great circle distance (degrees)
	} else if ((strUnit == "deg") || (strUnit == "degrees_north") || (strUnit == "degrees_east")) {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {

		} else if (strTargetUnit == "rad") {
			dValue *= M_PI / 180.0;

		} else if (strTargetUnit == "m") {
			dValue *= EarthRadius * M_PI / 180.0;

		} else if (strTargetUnit == "km") {
			dValue *= EarthRadius * M_PI / 180.0;

		} else {
			return false;
		}

	// Perform unit conversion from great circle distance (radians)
	} else if (strUnit == "rad") {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			dValue *= 180.0 / M_PI;

		} else if (strTargetUnit == "m") {
			dValue *= EarthRadius;

		} else if (strTargetUnit == "km") {
			dValue *= EarthRadius;

		} else {
			return false;
		}

	// Perform unit conversion from great circle distance (meters)
	// or altitude (meters)
	} else if (strUnit == "m") {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			dValue *= 180.0 / M_PI / EarthRadius;

		} else if (strTargetUnit == "rad") {
			dValue /= EarthRadius;

		} else if (strTargetUnit == "km") {
			dValue /= 1000.0;

		} else if (strTargetUnit == "m2/s2") {
			dValue *= EarthGravity;

		} else {
			return false;
		}

	// Perform unit conversion from great circle distance (kilometers)
	} else if (strUnit == "km") {
		if ((strTargetUnit == "deg") ||
		    (strTargetUnit == "degrees_north") ||
		    (strTargetUnit == "degrees_east")
		) {
			dValue *= 180.0 / M_PI / (EarthRadius / 1000.0);

		} else if (strTargetUnit == "rad") {
			dValue /= (EarthRadius / 1000.0);

		} else if (strTargetUnit == "m") {
			dValue *= 1000.0;

		} else if (strTargetUnit == "m2/s2") {
			dValue *= 1000.0 * EarthGravity;

		} else {
			return false;
		}

	// Perform unit conversion from temperature (K)
	} else if (strUnit == "K") {
		if ((strTargetUnit == "degC") || (strTargetUnit == "C")) {
			if (!fIsDelta) {
				dValue -= 273.15;
			}

		} else if ((strTargetUnit == "degF") || (strTargetUnit == "F")) {
			if (!fIsDelta) {
				dValue = 32.0 + 9.0 / 5.0 * (dValue - 273.15);
			} else {
				dValue *= 9.0 / 5.0;
			}

		} else {
			return false;
		}

	// Perform unit conversion from temperature (degC)
	} else if ((strUnit == "degC") || (strUnit == "C")) {
		if (strTargetUnit == "K") {
			if (!fIsDelta) {
				dValue += 273.15;
			}

		} else if ((strTargetUnit == "degF") || (strTargetUnit == "F")) {
			if (!fIsDelta) {
				dValue = 32.0 + 9.0 / 5.0 * dValue;
			} else {
				dValue *= 9.0 / 5.0;
			}

		} else {
			return false;
		}

	// Perform unit conversion from temperature (degF)
	} else if ((strUnit == "degF") || (strUnit == "F")) {
		if (strTargetUnit == "K") {
			if (!fIsDelta) {
				dValue = 273.15 + (dValue - 32.0) * 5.0 / 9.0;
			} else {
				dValue *= 5.0 / 9.0;
			}

		} else if ((strTargetUnit == "degC") || (strTargetUnit == "C")) {
			if (!fIsDelta) {
				dValue = (dValue - 32.0) * 5.0 / 9.0;
			} else {
				dValue *= 5.0 / 9.0;
			}

		} else {
			return false;
		}

	// Perform unit conversion from pressure (Pa)
	} else if (strUnit == "Pa") {
		if ((strTargetUnit == "hPa") ||
		    (strTargetUnit == "mb") ||
		    (strTargetUnit == "mbar")
		) {
			dValue /= 100.0;

		} else if (strTargetUnit == "atm") {
			dValue /= EarthAtmosphericPressure;

		} else {
			return false;
		}

	// Perform unit conversion from pressure (hPa,mb,mbar)
	} else if ((strUnit == "hPa") || (strUnit == "mb") || (strUnit == "mbar")) {
		if (strTargetUnit == "Pa") {
			dValue *= 100.0;

		} else if (
		    (strTargetUnit == "hPa") ||
		    (strTargetUnit == "mb") ||
		    (strTargetUnit == "mbar")
		) {

		} else if (strTargetUnit == "atm") {
			dValue /= (EarthAtmosphericPressure / 100.0);

		} else {
			return false;
		}

	// Perform unit conversion from pressure (atm)
	} else if (strUnit == "atm") {
		if (strTargetUnit == "Pa") {
			dValue *= EarthAtmosphericPressure;

		} else if (
		    (strTargetUnit == "hPa") ||
		    (strTargetUnit == "mb") ||
		    (strTargetUnit == "mbar")
		) {
			dValue *= (EarthAtmosphericPressure / 100.0);

		} else {
			return false;
		}

	// Perform unit conversion from geopotential (m2/s2)
	} else if (strUnit == "m2/s2") {
		if (strTargetUnit == "m") {
			dValue /= EarthGravity;

		} else if (strUnit == "km") {
			dValue /= (EarthGravity * 1000.0);

		} else {
			return false;
		}

	// Perform unit conversion from steradians (sr)
	} else if (strUnit == "sr") {
		if (strTargetUnit == "m2") {
			dValue *= EarthRadius * EarthRadius;

		} else if (strTargetUnit == "km2") {
			dValue *= EarthRadius * EarthRadius / (1000.0 * 1000.0);

		} else {
			return false;
		}

	// Perform unit conversion from square meters (m2)
	} else if (strUnit == "m2") {
		if (strTargetUnit == "sr") {
			dValue /= (EarthRadius * EarthRadius);

		} else if (strTargetUnit == "km2") {
			dValue /= (1000.0 * 1000.0);

		} else {
			return false;
		}

	// Perform unit conversion from square kilometers (km2)
	} else if (strUnit == "km2") {
		if (strTargetUnit == "sr") {
			dValue *= 1000.0 * 1000.0 / (EarthRadius * EarthRadius);

		} else if (strTargetUnit == "m2") {
			dValue *= 1000.0 * 1000.0;

		} else {
			return false;
		}

	// Perform unit conversion from kg m-2 s-1
	} else if (strUnit == "kg m-2 s-1"){
		if ((strTargetUnit == "mm/day") || (strTargetUnit == "mm/d") || (strTargetUnit == "mm d-1")) {
			dValue *= 86400.0;

		} else if ((strTargetUnit == "in/day") || (strTargetUnit == "in/d") || (strTargetUnit == "in d-1")) {
			dValue *= 86400.0 / 25.4;

		} else if ((strTargetUnit == "m s-1") || (strTargetUnit == "m/s")) {
			dValue /= 1000.0;

		} else {
			return false;
		}

	// Perform unit conversion from mm/day
	} else if ((strUnit == "mm/day") || (strUnit == "mm/d") || (strUnit == "mm d-1")) {
		if (strTargetUnit == "kg m-2 s-1") {
			dValue /= 86400.0;

		} else if ((strTargetUnit == "in/day") || (strTargetUnit == "in/d") || (strTargetUnit == "in d-1")) {
			dValue /= 25.4;

		} else if ((strTargetUnit == "m s-1") || (strTargetUnit == "m/s")) {
			dValue /= (1000.0 * 86400.0);

		} else {
			return false;
		}

	// Perform unit conversion from in/day
	} else if ((strUnit == "in/day") || (strUnit == "in/d") || (strUnit == "in d-1")) {
		if (strTargetUnit == "kg m-2 s-1") {
			dValue *= 25.4 / 86400.0;

		} else if ((strTargetUnit == "mm/day") || (strTargetUnit == "mm/d") || (strTargetUnit == "mm d-1")) {
			dValue *= 25.4;

		} else if ((strTargetUnit == "m s-1") || (strTargetUnit == "m/s")) {
			dValue *= 25.4 / (1000.0 * 86400.0);

		} else {
			return false;
		}

	// Perform unit conversion from m/s
	} else if ((strUnit == "m s-1") || (strUnit == "m/s")) {
		if (strTargetUnit == "kg m-2 s-1") {
			dValue *= 1000.0;

		} else if ((strTargetUnit == "mm/day") || (strTargetUnit == "mm/d") || (strTargetUnit == "mm d-1")) {
			dValue *= 1000.0 * 86400.0;

		} else if ((strUnit == "in/day") || (strUnit == "in/d") || (strUnit == "in d-1")) {
			dValue *= 1000.0 * 86400.0 / 25.4;

		} else {
			return false;
		}

	} else {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

inline static void SplitIntoValueAndUnits(
	const std::string & str,
	std::string & strValue,
	std::string & strUnits
) {
	size_t iFirstNonValue = std::string::npos;

	bool fIsFloat = false;
	bool fHasExponent = false;
	bool fHasDecimal = false;
	for(size_t i = 0; i < str.length(); i++) {
		if ((str[i] < '0') || (str[i] > '9')) {
			if (str[i] == '.') {
				if (fHasDecimal) {
					iFirstNonValue = i;
					break;
				}
				if (fHasExponent) {
					iFirstNonValue = i;
					break;
				}
				fHasDecimal = true;
				continue;
			}
			if (str[i] == 'e') {
				if (fHasExponent) {
					if ((i != 0) && (str[i-1] == 'e')) {
						iFirstNonValue = i-1;
					} else {
						iFirstNonValue = i;
					}
					break;
				}
				fHasExponent = true;
				continue;
			}
			if ((str[i] == '-') || (str[i] == '+')) {
				if (i == 0) {
					continue;
				} else if (str[i-1] == 'e') {
					continue;
				} else {
					iFirstNonValue = i;
					break;
				}
			}
			iFirstNonValue = i;
			break;
		}
	}

	strValue = str.substr(0, iFirstNonValue);
	if (iFirstNonValue == std::string::npos) {
		strUnits = "";
	} else {
		strUnits = str.substr(iFirstNonValue, std::string::npos);
	}
}

///////////////////////////////////////////////////////////////////////////////

#endif // _UNITS_H_
