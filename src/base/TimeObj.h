///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimeObj.h
///	\author  Paul Ullrich
///	\version April 18, 2018
///
///	<remarks>
///		Copyright 2000- Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TIMEOBJ_H_
#define _TIMEOBJ_H_

#include "Exception.h"
#include "STLStringHelper.h"

#include <string>
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for storing a time as year, month, day, seconds and useconds.
///	</summary>
class Time {

public:
	///	<summary>
	///		Type of calendar.
	///	</summary>
	enum CalendarType {
		CalendarUnknown,
		CalendarNone,
		CalendarNoLeap,
		CalendarStandard,
		CalendarGregorian,
		Calendar360Day,
	};

	///	<summary>
	///		Type of time.
	///	</summary>
	enum TimeType {
		TypeFixed,
		TypeDelta
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Time(
		CalendarType eCalendarType = CalendarNoLeap,
		TimeType eTimeType = TypeFixed
	) :
		m_iYear(0),
		m_iMonth(0),
		m_iDay(0),
		m_iSecond(0),
		m_iMicroSecond(0),
		m_eCalendarType(eCalendarType),
		m_eTimeType(eTimeType)
	{
		if (m_eCalendarType == CalendarUnknown) {
			_EXCEPTIONT("Invalid CalendarType");
		}
	}

	///	<summary>
	///		Constructor.
	///	</summary>
	Time(
		int iYear,
		int iMonth,
		int iDay,
		int iSecond,
		int iMicroSecond,
		CalendarType eCalendarType = CalendarNoLeap,
		TimeType eTimeType = TypeFixed
	) :
		m_iYear(iYear),
		m_iMonth(iMonth),
		m_iDay(iDay),
		m_iSecond(iSecond),
		m_iMicroSecond(iMicroSecond),
		m_eCalendarType(eCalendarType),
		m_eTimeType(eTimeType)
	{
		if (m_eCalendarType == CalendarUnknown) {
			_EXCEPTIONT("Invalid CalendarType");
		}
		NormalizeTime();
	}

	///	<summary>
	///		Constructor from std::string.
	///	</summary>
	Time(
		const std::string & strTime,
		CalendarType eCalendarType = CalendarNoLeap,
		TimeType eTimeType = TypeFixed
	) :
		m_iYear(0),
		m_iMonth(0),
		m_iDay(0),
		m_iSecond(0),
		m_iMicroSecond(0),
		m_eCalendarType(eCalendarType),
		m_eTimeType(eTimeType)
	{
		if (m_eCalendarType == CalendarUnknown) {
			_EXCEPTIONT("Invalid CalendarType");
		}
		FromLongString(strTime);
	}

public:
	///	<summary>
	///		Returns the CalendarType associated with the given string.
	///	</summary>
	static CalendarType CalendarTypeFromString(
		const std::string & strCalendar
	) {
		std::string strCalendarTemp = strCalendar;
		STLStringHelper::ToLower(strCalendarTemp);

		if (strCalendarTemp == "none") {
			return CalendarNone;
		} else if (strCalendarTemp == "noleap") {
			return CalendarNoLeap;
		} else if (strCalendarTemp == "365_day"){
			return CalendarNoLeap;
		} else if (strCalendarTemp == "standard") {
			return CalendarStandard;
		} else if (strCalendarTemp == "gregorian") {
			return CalendarGregorian;
		} else if (strCalendarTemp == "proleptic_gregorian") {
			return CalendarGregorian;
		} else if (strCalendarTemp == "360_day") {
			return Calendar360Day;
		} else {
			return CalendarUnknown;
		}
	}

public:
	///	<summary>
	///		Determine if two Times occur on the same date.
	///	</summary>
	inline bool IsSameDate(const Time & time) {
		if (m_eTimeType != time.m_eTimeType) {
			return false;
		}
		if (m_eCalendarType != time.m_eCalendarType) {
			return false;
		}

		if ((m_iYear  == time.m_iYear) &&
		    (m_iMonth == time.m_iMonth) &&
		    (m_iDay   == time.m_iDay)
		) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Equality between Times.
	///	</summary>
	bool operator==(const Time & time) const;

	///	<summary>
	///		Inequality between Times.
	///	</summary>
	bool operator!=(const Time & time) const {
		return !((*this) == time);
	}

	///	<summary>
	///		Less-than between Times.
	///	</summary>
	bool operator<(const Time & time) const;

	///	<summary>
	///		Greater-than comparator between Times.
	///	</summary>
	bool operator>(const Time & time) const;

	///	<summary>
	///		Less-than-or-equal comparator between Times.
	///	</summary>
	bool operator<=(const Time & time) const {
		return !((*this) > time);
	}

	///	<summary>
	///		Greater-than-or-equal comparator between Times.
	///	</summary>
	bool operator>=(const Time & time) const {
		return !((*this) < time);
	}

protected:
	///	<summary>
	///		Verify that the Time is in accordance with the calendar.
	///	</summary>
	void VerifyTime();

	///	<summary>
	///		Normalize the Time, in accordance with the Calendar.
	///	</summary>
	void NormalizeTime();

public:
	///	<summary>
	///		Add Time to this object.
	///	</summary>
	void AddTime(const Time & timeDelta);

	///	<summary>
	///		Add Time to this object.
	///	</summary>
	inline Time & operator+=(const Time & timeDelta) {
		AddTime(timeDelta);
		return (*this);
	}

	///	<summary>
	///		Add a number of seconds to the Time.
	///	</summary>
	inline void AddSeconds(int nSeconds) {
		m_iSecond += nSeconds;

		NormalizeTime();
	}

	///	<summary>
	///		Add a number of seconds to the Time.
	///	</summary>
	inline Time & operator+=(int nSeconds) {
		AddSeconds(nSeconds);
		return (*this);
	}

	///	<summary>
	///		Add a number of minutes to the Time.
	///	</summary>
	inline void AddMinutes(int nMinutes) {
		AddHours(nMinutes / 60);
		AddSeconds((nMinutes % 60) * 60);
	}

	///	<summary>
	///		Add a number of hours to the Time.
	///	</summary>
	inline void AddHours(int nHours) {
		AddDays(nHours / 24);
		AddSeconds((nHours % 24) * 3600);
	}

	///	<summary>
	///		Add a number of days to the Time.
	///	</summary>
	inline void AddDays(int nDays) {
		m_iDay += nDays;

		NormalizeTime();
	}

	///	<summary>
	///		Add a number of months to the Time.
	///	</summary>
	inline void AddMonths(int nMonths) {
		m_iMonth += nMonths;

		NormalizeTime();
	}

	///	<summary>
	///		Add a number of years to the Time.
	///	</summary>
	inline void AddYears(int nYears) {
		m_iYear += nYears;

		NormalizeTime();
	}

public:
/*
	///	<summary>
	///		Add a number of seconds to the Time to produce a new Time value.
	///	</summary>
	Time operator+(double dSeconds) const;
*/
	///	<summary>
	///		Calculate the day number for this Time.
	///	</summary>
	int DayNumber() const;

	///	<summary>
	///		Determine the number of seconds between two Times.
	///	</summary>
	double operator-(const Time & time) const;

	///	<summary>
	///		Determine the number of seconds from this time to the given Time.
	///	</summary>
	double DeltaSeconds(const Time & time) const;

	///	<summary>
	///		Determine the number of minutes from this time to the given Time.
	///	</summary>
	double DeltaMinutes(const Time & time) const;

	///	<summary>
	///		Determine the number of hours from this time to the given Time.
	///	</summary>
	double DeltaHours(const Time & time) const;

	///	<summary>
	///		Determine the number of days from this time to the given Time.
	///	</summary>
	double DeltaDays(const Time & time) const;

public:
	///	<summary>
	///		Check if this object is zero.
	///	</summary>
	inline bool IsZero() const {
		if ((m_iYear == 0) &&
			(m_iMonth == 0) &&
			(m_iDay == 0) &&
			(m_iSecond == 0) &&
			(m_iMicroSecond == 0)
		) {
			return true;
		}

		return false;
	}

	///	<summary>
	///		Get the year.
	///	</summary>
	inline int GetYear() const {
		return m_iYear;
	}

	///	<summary>
	///		Get the month.
	///	</summary>
	inline int GetMonth() const {
		if ((m_eTimeType == TypeFixed) && (m_eCalendarType != CalendarNone)) {
			return m_iMonth + 1;
		} else {
			return m_iMonth;
		}
	}

	///	<summary>
	///		Get the day.
	///	</summary>
	inline int GetDay() const {
		if ((m_eTimeType == TypeFixed) && (m_eCalendarType != CalendarNone)) {
			return m_iDay + 1;
		} else {
			return m_iDay;
		}
	}

	///	<summary>
	///		Get the number of seconds.
	///	</summary>
	inline int GetSecond() const {
		return m_iSecond;
	}

	///	<summary>
	///		Get the number of microseconds.
	///	</summary>
	inline int GetMicroSecond() const {
		return m_iMicroSecond;
	}

	///	<summary>
	///		Get the CalendarType.
	///	</summary>
	inline CalendarType GetCalendarType() const {
		return m_eCalendarType;
	}

	///	<summary>
	///		Get the TimeType.
	///	</summary>
	inline TimeType GetTimeType() const {
		return m_eTimeType;
	}

	///	<summary>
	///		Set the year.
	///	</summary>
	inline void SetYear(int iYear) {
		m_iYear = iYear;
	}

	///	<summary>
	///		Set the month.
	///	</summary>
	inline void SetMonth(int iMonth) {
		if ((m_eTimeType == TypeFixed) && (m_eCalendarType != CalendarNone)) {
			m_iMonth = iMonth - 1;
		} else {
			m_iMonth = iMonth;
		}
	}

	///	<summary>
	///		Set the day.
	///	</summary>
	inline void SetDay(int iDay) {
		if ((m_eTimeType == TypeFixed) && (m_eCalendarType != CalendarNone)) {
			m_iDay = iDay - 1;
		} else {
			m_iDay = iDay;
		}
	}

	///	<summary>
	///		Set the number of seconds.
	///	</summary>
	inline void SetSecond(int iSecond) {
		m_iSecond = iSecond;
	}

	///	<summary>
	///		Set the number of microseconds.
	///	</summary>
	inline void SetMicroSecond(int iMicroSecond) {
		m_iMicroSecond = iMicroSecond;
	}

	///	<summary>
	///		Set the time type.
	///	</summary>
	inline void SetTimeType(TimeType eTimeType) {
		m_eTimeType = eTimeType;
	}

public:
	///	<summary>
	///		Get the Time as a string, only including the date: "YYYY-MM-DD"
	///	</summary>
	std::string ToDateString() const;

	///	<summary>
	///		Get the Time as a short string: "YYYY-MM-DD-SSSSS"
	///	</summary>
	std::string ToShortString() const;

	///	<summary>
	///		Get the Time as a long string: "YYYY-MM-DD-SSSSS-UUUUUU"
	///	</summary>
	std::string ToLongString() const;

	///	<summary>
	///		Get the Time as a full-length string showing hours, days
	///		and seconds: "YYYY-MM-DD hh:mm:ss[.uuuuuu]"
	///	</summary>
	std::string ToString() const;

	///	<summary>
	///		Get the Time as a free string.
	///	</summary>
	std::string ToFreeString() const;

	///	<summary>
	///		Set the Time using a long string: "YYYY-MM-DD-SSSSS-UUUUUU"
	///	</summary>
	void FromLongString(const std::string & strLongTime);

	///	<summary>
	///		Set the Time using a formatted string:
	///		- 15h (15 hours)
	///		- 15h30s (15 hours 30 seconds)
	///		- 7y (7 years)
	///		- 1500s (1500 seconds)
	///		- 2y1500u (2 years 1500 microseconds)
	///		- yyyy
	///		- yyyy-MM
	///		- yyyy-MM-dd
	///		- yyyy-MM-dd-sssss
	///		- yyyy-MM-dd-hh:mm:ss.uuuuuu
	///		- HH:MM:SS.UUUU
	///	</summary>
	void FromFormattedString(const std::string & strFormattedTime);

	///	<summary>
	///		Set the Time using a CF-compliant time unit string.
	///		- "hours since ..."
	///		- "seconds since ..."
	///	</summary>
	void FromCFCompliantUnitsOffsetInt(
		const std::string & strFormattedTime,
		int nOffset
	);

	///	<summary>
	///		Set the Time using a CF-compliant time unit string.
	///		- "hours since ..."
	///		- "seconds since ..."
	///	</summary>
	void FromCFCompliantUnitsOffsetDouble(
		const std::string & strFormattedTime,
		double dOffset
	);

	///	<summary>
	///		Get the Time using a CF-compliant time unit string.
	///		- "hours since ..."
	///		- "seconds since ..."
	///	</summary>
	double GetCFCompliantUnitsOffsetDouble(
		const std::string & strFormattedTime
	);

	///	<summary>
	///		Get the name of the calendar.
	///	</summary>
	std::string GetCalendarName() const;

private:
	///	<summary>
	///		The year.
	///	</summary>
	int m_iYear;

	///	<summary>
	///		The month.
	///	</summary>
	int m_iMonth;

	///	<summary>
	///		The day.
	///	</summary>
	int m_iDay;

	///	<summary>
	///		The second count.
	///	</summary>
	int m_iSecond;

	///	<summary>
	///		The microsecond count.
	///	</summary>
	int m_iMicroSecond;

	///	<summary>
	///		Calendar type.
	///	</summary>
	CalendarType m_eCalendarType;

	///	<summary>
	///		Time type.
	///	</summary>
	TimeType m_eTimeType;
};

///////////////////////////////////////////////////////////////////////////////

#endif

