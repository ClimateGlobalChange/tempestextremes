///////////////////////////////////////////////////////////////////////////////
///
///	\file    Constants.h
///	\author  Paul Ullrich
///	\version September 17, 2019
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

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

///////////////////////////////////////////////////////////////////////////////
// Standard Earth Radius in meters
static const double EarthRadius = 6.37122e6;

///////////////////////////////////////////////////////////////////////////////
// Rotation rate of the planet, in inverse seconds
static const double EarthOmega = 7.2921e-5;

///////////////////////////////////////////////////////////////////////////////
// Gravitational acceleration at the surface, in meters per second squared
static const double EarthGravity = 9.80616;

///////////////////////////////////////////////////////////////////////////////
// Earth standard atmospheric pressure, in Pascals
static const double EarthAtmosphericPressure = 101325.0;

///////////////////////////////////////////////////////////////////////////////
// Knots per meter/second
static const double KnotsPerMetersPerSecond = 1.94384;

#endif // _CONSTANTS_H_

