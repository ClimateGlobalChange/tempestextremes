///////////////////////////////////////////////////////////////////////////////
///
///	\file    LatLonBox.h
///	\author  Paul Ullrich
///	\version August 29, 2023
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

#ifndef _LATLONBOX_H_
#define _LATLONBOX_H_

#include <sstream>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A structure for storing a bounding box in latitude / longitude space.
///		The box is treated as singularly periodic in the longitudinal direction
///		if lon_periodic is set to true.
///	</summary>
template <typename Type>
class LatLonBox {

public:
	///	<summary>
	///		Flag indicating this is a null box.
	///	</summary>
	bool is_null;

	///	<summary>
	///		Flag indicating this is a regional.
	///	</summary>
	bool lon_periodic;

	///	<summary>
	///		Width of the longitude variable.
	///	</summary>
	Type lon_width;

	///	<summary>
	///		Bounding longitudes (endpoints are included).
	///	</summary>
	Type lon[2];

	///	<summary>
	///		Bounding latitudes (endpoints are included).
	///	</summar>
	Type lat[2];

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	LatLonBox(
		bool a_lon_periodic = true,
		Type a_lon_width = static_cast<Type>(360)
	) :
		is_null(true),
		lon_periodic(a_lon_periodic),
		lon_width(a_lon_width)
	{
		lon[0] = static_cast<Type>(0);
		lon[1] = static_cast<Type>(0);
		lat[0] = static_cast<Type>(0);
		lat[1] = static_cast<Type>(0);
	}

	///	<summary>
	///		Constructor with longitude-latitude coordinates.
	///	</summary>
	LatLonBox(
		Type a_lon0,
		Type a_lon1,
		Type a_lat0,
		Type a_lat1,
		bool a_lon_periodic = true,
		Type a_lon_width = static_cast<Type>(360)
	) :
		is_null(false),
		lon_periodic(a_lon_periodic),
		lon_width(a_lon_width)
	{
		set(a_lon0, a_lon1, a_lat0, a_lat1);
	}

	///	<summary>
	///		Set the latlonbox bounds to the specified values.
	///	</summary>
	void set(
		Type a_lon0,
		Type a_lon1,
		Type a_lat0,
		Type a_lat1
	) {
		lon[0] = a_lon0;
		lon[1] = a_lon1;
		lat[0] = a_lat0;
		lat[1] = a_lat1;

		if (lon[0] < static_cast<Type>(0)) {
			std::stringstream strError;
			strError << "lon_pt0 out of range (" << lon[0] << " < 0)" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}
		if (lon[0] > lon_width) {
			std::stringstream strError;
			strError << "lon_pt0 out of range (" << lon[0] << " > ";
			strError << lon_width << ")" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}
		if (lon[1] < static_cast<Type>(0)) {
			std::stringstream strError;
			strError << "lon_pt1 out of range (" << lon[1] << " < 0)" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}
		if (lon[1] > lon_width) {
			std::stringstream strError;
			strError << "lon_pt1 out of range (" << lon[1] << " > ";
			strError << lon_width << ")" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}
	}

	///	<summary>
	///		Width of this LatLonBox.
	///	</summary>
	Type width() const {
		if (is_null) {
			return 0;
		}

		if (!lon_periodic) {
			return (lon[1] - lon[0]);
		}

		if (lon[0] == lon[1]) {
			return static_cast<Type>(0);
		} else if (lon[0] <= lon[1]) {
			return (lon[1] - lon[0]);
		} else {
			return (lon[1] - lon[0] + lon_width);
		}
	}

	///	<summary>
	///		Insert a point into this LatLonBox.
	///	</summary>
	void insert(
		const Type & lat_pt,
		const Type & lon_pt
	) {
		if (is_null) {
			is_null = false;
			lat[0] = lat_pt;
			lat[1] = lat_pt;
			lon[0] = lon_pt;
			lon[1] = lon_pt;
			return;
		}
		if (lon_pt < static_cast<Type>(0)) {
			std::stringstream strError;
			strError << "lon_pt out of range (" << lon_pt << " < 0)" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}
		if (lon_pt > lon_width) {
			std::stringstream strError;
			strError << "lon_pt out of range (" << lon_pt << " > ";
			strError << lon_width << ")" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}

		// Expand latitudes
		if (lat_pt > lat[1]) {
			lat[1] = lat_pt;
		}
		if (lat_pt < lat[0]) {
			lat[0] = lat_pt;
		}

		// Expand longitude, if non-periodic
		if (!lon_periodic) {
			if (lon_pt > lon[1]) {
				lon[1] = lon_pt;
			}
			if (lon_pt < lon[0]) {
				lon[0] = lon_pt;
			}
			return;
		}

		// New longitude lies within existing range
		if (lon[0] <= lon[1]) {
			if ((lon_pt >= lon[0]) && (lon_pt <= lon[1])) {
				return;
			}
		} else {
			if ((lon_pt >= lon[0]) || (lon_pt <= lon[1])) {
				return;
			}
		}

		// New longitude lies outside of existing range
		LatLonBox boxA(*this);
		boxA.lon[0] = lon_pt;

		LatLonBox boxB(*this);
		boxB.lon[1] = lon_pt;

		// The updated box is the box of minimum width
		Type dWidthNow = width();
		Type dWidthA = boxA.width();
		Type dWidthB = boxB.width();

		if ((dWidthA - dWidthNow < -1.0e-14) || (dWidthB - dWidthNow < -1.0e-14)) {
			_EXCEPTIONT("Logic error");
		}

		if (dWidthA < dWidthB) {
			(*this) = boxA;
		} else {
			(*this) = boxB;
		}
	}

	///	<summary>
	///		Determine if this LatLonBox contains the given point.
	///	</summary>
	bool contains(
		const Type & lat_pt,
		const Type & lon_pt
	) {
		// Sanity check
		if (!lon_periodic && (lon[0] > lon[1])) {
			_EXCEPTION2("Maximum longitude (%1.5f) is smaller than minimum longitude (%1.5f) in non-periodic LatLonBox",
				lon[0], lon[1]);
		}

		// Check latitudes
		if (lat[0] > lat_pt) {
			return false;
		}
		if (lat[1] < lat_pt) {
			return false;
		}

		// Check longitudes, if non-periodic
		if (!lon_periodic) {
			if (lon[0] > lon_pt) {
				return false;
			}
			if (lon[1] < lon_pt) {
				return false;
			}
			return true;
		}

		// This box crosses lon 360
		if (lon[0] > lon[1]) {
			if (lon_pt >= lon[0]) {
				return true;
			}
			if (lon_pt <= lon[1]) {
				return true;
			}
			return false;
		}

		// This box does not cross lon 360
		if ((lon_pt >= lon[0]) && (lon_pt <= lon[1])) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Determine if this LatLonBox overlaps with box.
	///	</summary>
	bool overlaps(
		const LatLonBox<Type> & box
	) const {
		// Sanity check
		if (!lon_periodic && (lon[0] > lon[1])) {
			_EXCEPTION2("Maximum longitude (%1.5f) is smaller than minimum longitude (%1.5f) in non-periodic LatLonBox",
				lon[0], lon[1]);
		}

		// Check flags
		if ((is_null) || (box.is_null)) {
			return false;
		}
		if (lon_periodic != box.lon_periodic) {
			_EXCEPTION2("Inconsistent periodicity flag (%s/%s)",
				(lon_periodic)?("true"):("false"),
				(box.lon_periodic)?("true"):("false"));
		}
		if (lon_width != box.lon_width) {
			_EXCEPTIONT("Inconsistent box lon_width in comparison");
		}

		// Check latitudes
		if (lat[0] > box.lat[1]) {
			return false;
		}
		if (lat[1] < box.lat[0]) {
			return false;
		}

		// Cases when longitude is periodic
		if (lon_periodic) {

			// Both boxes cross lon 360
			if ((lon[0] > lon[1]) && (box.lon[0] > box.lon[1])) {
				return true;
			}

			// This box crosses lon 360
			if (lon[0] > lon[1]) {
				if (box.lon[1] >= lon[0]) {
					return true;
				}
				if (box.lon[0] <= lon[1]) {
					return true;
				}
				return false;
			}

			// That box crosses lon 360
			if (box.lon[0] > box.lon[1]) {
				if (lon[1] >= box.lon[0]) {
					return true;
				}
				if (lon[0] <= box.lon[1]) {
					return true;
				}
				return false;
			}
		}

		// No boxes cross lon 360
		if (box.lon[1] < lon[0]) {
			return false;
		}
		if (box.lon[0] > lon[1]) {
			return false;
		}
		return true;
	}
};

///////////////////////////////////////////////////////////////////////////////

#endif // _LATLONBOX_H_

