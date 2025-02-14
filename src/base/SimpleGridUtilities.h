///////////////////////////////////////////////////////////////////////////////
///
///	\file    SimpleGridUtilities.h
///	\author  Paul Ullrich
///	\version August 14, 2018
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

#ifndef _SIMPLEGRIDUTILITIES_H_
#define _SIMPLEGRIDUTILITIES_H_

#include "SimpleGrid.h"
#include "DataArray1D.h"

#include <set>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the minimum/maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDistDeg">
///		Maximum distance from the initial point in degrees.
///	</param>
///	<param name="ixExtremum">
///		Output node index at which the extremum occurs.
///	</param>
///	<param name="dExtremumValue">
///		Output value of the field taken at the extremum point.
///	</param>
///	<param name="dRMax">
///		Output distance from the centerpoint at which the extremum occurs
///		in great circle distance (degrees).
///	</param>
template <typename real>
void FindLocalMinMax(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataArray1D<real> & data,
	int ix0,
	double dMaxDistDeg,
	int & ixExtremum,
	real & dExtremumValue,
	float & dRMaxDeg
);

///	<summary>
///		Find the locations of all minima in the given DataArray1D.
///	</summary>
template <typename real>
void FindAllLocalMinima(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	std::set<int> & setMinima
);

///	<summary>
///		Find the locations of all maxima in the given DataArray1D.
///	</summary>
template <typename real>
void FindAllLocalMaxima(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	std::set<int> & setMaxima
);

///	<summary>
///		Find the locations of all minima in the given DataArray1D
///		with a prescribed threshold.
///	</summary>
template <typename real>
void FindAllLocalMinMaxWithThreshold(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	bool fMinima,
	const std::string & strThreshold,
	std::set<int> & setMinima
);

///	<summary>
///		Find the locations of all local min/max in the given DataArray1D
///		for a given search distance.
///	</summary>
template <typename real>
void FindAllLocalMinMaxWithGraphDistance(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	bool fMinima,
	int nMaxGraphDistance,
	std::set<int> & setMinMax
);

///	<summary>
///		Find the local average of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
template <typename real>
void FindLocalAverage(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	int ix0,
	double dMaxDistDeg,
	real & dAverage
);

///	<summary>
///		Find the largest positive or negative value of the closed contour
///		delta in the given field at the given grid point.
///	</summary>
template <typename real>
void FindMaxClosedContourDelta(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	int ix0,
	double dDistDeg,
	double dMinMaxDistDeg,
	bool fMaxClosedContourDeltaSign,
	real & dMaxClosedContourDelta
);

///	<summary>
///		Find the weighted area of positive values in a given radius
///		minus the weighted area of negative values in that radius.
///	</summary>
template <typename real>
void PositiveMinusNegativeWeightedArea(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	int ix0,
	double dDistDeg,
	real & dValue);

///	<summary>
///		Find the maximum value of the data field that is poleward of
///		the given candidate within a given longitude swath.
///	</summary>
template <typename real>
void MaxPolewardValue(
	const SimpleGrid & grid,
	const DataArray1D<real> & data,
	int ix0,
	double dDistDeg,
	real & dValue);

///////////////////////////////////////////////////////////////////////////////

#endif // _SIMPLEGRIDUTILITIES_H_

