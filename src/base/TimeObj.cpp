///////////////////////////////////////////////////////////////////////////////
///
///	\file    TimeObj.cpp
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

#include "TimeObj.h"
#include "Exception.h"

#include <iostream>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////

bool Time::operator==(const Time & time) const {
	if ((m_eCalendarType != time.m_eCalendarType) ||
	    (m_eTimeType     != time.m_eTimeType)     ||
	    (m_iYear        != time.m_iYear)          ||
	    (m_iMonth       != time.m_iMonth)         ||
	    (m_iDay         != time.m_iDay)           ||
		(m_iSecond      != time.m_iSecond)        ||
		(m_iMicroSecond != time.m_iMicroSecond)
	) {
		return false;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

bool Time::operator<(const Time & time) const {
	if (m_eTimeType != time.m_eTimeType) {
		_EXCEPTIONT("Cannot compare Time objects with different types");
	}
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTION2("Cannot compare Time objects with different calendars (%s/%s)",
			GetCalendarName().c_str(),
			time.GetCalendarName().c_str());
	}

	if (m_iYear < time.m_iYear) {
		return true;
	} else if (m_iYear > time.m_iYear) {
		return false;
	}

	if (m_iMonth < time.m_iMonth) {
		return true;
	} else if (m_iMonth > time.m_iMonth) {
		return false;
	}

	if (m_iDay < time.m_iDay) {
		return true;
	} else if (m_iDay > time.m_iDay) {
		return false;
	}

	if (m_iSecond < time.m_iSecond) {
		return true;
	} else if (m_iSecond > time.m_iSecond) {
		return false;
	}

	if (m_iMicroSecond < time.m_iMicroSecond) {
		return true;
	} else if (m_iMicroSecond > time.m_iMicroSecond) {
		return false;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

bool Time::operator>(const Time & time) const {
	if (m_eTimeType != time.m_eTimeType) {
		_EXCEPTIONT("Cannot compare Time objects with different types");
	}
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT(
			"Cannot compare Time objects with different calendars");
	}

	if (m_iYear > time.m_iYear) {
		return true;
	} else if (m_iYear < time.m_iYear) {
		return false;
	}

	if (m_iMonth > time.m_iMonth) {
		return true;
	} else if (m_iMonth < time.m_iMonth) {
		return false;
	}

	if (m_iDay > time.m_iDay) {
		return true;
	} else if (m_iDay < time.m_iDay) {
		return false;
	}

	if (m_iSecond > time.m_iSecond) {
		return true;
	} else if (m_iSecond < time.m_iSecond) {
		return false;
	}

	if (m_iMicroSecond > time.m_iMicroSecond) {
		return true;
	} else if (m_iMicroSecond < time.m_iMicroSecond) {
		return false;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

void Time::VerifyTime() {

	// Verification only for known CalendarTypes
	if ((m_eCalendarType == CalendarNoLeap) || 
		(m_eCalendarType == CalendarStandard) ||
		(m_eCalendarType == CalendarGregorian) ||
		(m_eCalendarType == Calendar360Day) ||
		(m_eCalendarType == Calendar365Day)
	) {
		int nDaysPerMonth[]
			= {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

		if ((m_eCalendarType == CalendarStandard) ||
			(m_eCalendarType == CalendarGregorian)
		) {
			if ((m_iYear % 4) == 0) {
				nDaysPerMonth[1] = 29;
				if (((m_iYear % 100) == 0) && ((m_iYear % 400) != 0)) {
					nDaysPerMonth[1] = 28;
				}
			}
		}
		if (m_eCalendarType == Calendar360Day) {
			for (int i = 0; i < 12; i++) {
				nDaysPerMonth[i] = 30;
			}
		}

		if ((m_iYear < 0) || (m_iYear > 9999)) {
			_EXCEPTIONT("Year out of range");
		}
		if ((m_iMonth < 0) || (m_iMonth > 11)) {
			_EXCEPTIONT("Month out of range");
		}
		if ((m_iDay < 0) || (m_iDay >= nDaysPerMonth[m_iMonth])) {
			_EXCEPTIONT("Day out of range");
		}
		if ((m_iSecond < 0) || (m_iSecond >= 86400)) {
			_EXCEPTIONT("Seconds out of range");
		}
		if ((m_iMicroSecond < 0) || (m_iMicroSecond >= 1000000)) {
			_EXCEPTIONT("MicroSecond out of range");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Time::NormalizeTime() {

	// Normalization only for fixed time
	if (m_eTimeType == TypeDelta) {
		return;
	}

	// Normalization only for known CalendarTypes
	if ((m_eCalendarType == CalendarNoLeap) || 
		(m_eCalendarType == CalendarStandard) ||
		(m_eCalendarType == CalendarGregorian) ||
		(m_eCalendarType == Calendar360Day) ||
		(m_eCalendarType == Calendar365Day)
	) {
		int nDaysPerMonth[]
			= {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

		if ((m_eCalendarType == CalendarStandard) ||
			(m_eCalendarType == CalendarGregorian)
		) {
			if ((m_iYear % 4) == 0) {
				nDaysPerMonth[1] = 29;
				if (((m_iYear % 100) == 0) && ((m_iYear % 400) != 0)) {
					nDaysPerMonth[1] = 28;
				}
			}
		}
		if (m_eCalendarType == Calendar360Day) {
			for (int i = 0; i < 12; i++) {
				nDaysPerMonth[i] = 30;
			}
		}

		// Add seconds
		int nAddedSeconds = 0;
		if (m_iMicroSecond >= 1000000) {
			nAddedSeconds = m_iMicroSecond / 1000000;

		} else if (m_iMicroSecond < 0) {
			nAddedSeconds = - (1000000 - m_iSecond) / 1000000;
		}
		m_iMicroSecond -= nAddedSeconds * 1000000;
		m_iSecond += nAddedSeconds;

		// Add days
		int nAddedDays = 0;
		if (m_iSecond >= 86400) {
			nAddedDays = m_iSecond / 86400;

		} else if (m_iSecond < 0) {
			nAddedDays = - (86400 - m_iSecond) / 86400;
		}

		m_iSecond -= nAddedDays * 86400;
		m_iDay += nAddedDays;

		// Add years
		int nAddedYears = 0;
		if (m_iMonth >= 12) {
			nAddedYears = m_iMonth / 12;
		} else if (m_iMonth < 0) {
			nAddedYears = - (12 - m_iMonth) / 12;
		}
		m_iMonth -= nAddedYears * 12;
		m_iYear += nAddedYears;

		if ((m_iMonth < 0) || (m_iMonth >= 12)) {
			_EXCEPTIONT("Logic error");
		}

		// Subtract months
		while (m_iDay < 0) {
			m_iMonth--;
			if (m_iMonth < 0) {
				m_iMonth = 11;
				m_iYear--;

				// Adjust number of days per month
				if ((m_eCalendarType == CalendarStandard) ||
					(m_eCalendarType == CalendarGregorian)
				) {
					nDaysPerMonth[1] = 28;
					if ((m_iYear % 4) == 0) {
						nDaysPerMonth[1] = 29;
						if (((m_iYear % 100) == 0) && ((m_iYear % 400) != 0)) {
							nDaysPerMonth[1] = 28;
						}
					}
				}
			}
			m_iDay += nDaysPerMonth[m_iMonth];
		}

		// Add months
		while (m_iDay >= nDaysPerMonth[m_iMonth]) {
			m_iDay -= nDaysPerMonth[m_iMonth];
			m_iMonth++;
			if (m_iMonth > 11) {
				m_iMonth = 0;
				m_iYear++;

				// Adjust number of days per month
				if ((m_eCalendarType == CalendarStandard) ||
					(m_eCalendarType == CalendarGregorian)
				) {
					nDaysPerMonth[1] = 28;
					if ((m_iYear % 4) == 0) {
						nDaysPerMonth[1] = 29;
						if (((m_iYear % 100) == 0) && ((m_iYear % 400) != 0)) {
							nDaysPerMonth[1] = 28;
						}
					}
				}
			}
		}

		// Check that the result is ok
		if ((m_iMonth < 0) || (m_iMonth >= 12) ||
		    (m_iDay < 0) || (m_iDay >= nDaysPerMonth[m_iMonth]) ||
		    (m_iSecond < 0) || (m_iSecond >= 86400) ||
			(m_iMicroSecond < 0) || (m_iMicroSecond >= 1e6)
		) {
			_EXCEPTION5("Logic error: %i %i %i %i %i",
				m_iYear, m_iMonth, m_iDay, m_iSecond, m_iMicroSecond);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Time::AddTime(const Time & timeDelta) {

	if (timeDelta.m_eTimeType != TypeDelta) {
		_EXCEPTIONT("Argument to AddTime() is not a TimeDelta");
	}

	m_iYear        += timeDelta.m_iYear;
	m_iMonth       += timeDelta.m_iMonth;
	m_iDay         += timeDelta.m_iDay;
	m_iSecond      += timeDelta.m_iSecond;
	m_iMicroSecond += timeDelta.m_iMicroSecond;

	NormalizeTime();
}

///////////////////////////////////////////////////////////////////////////////

void Time::SubtractTime(const Time & timeDelta) {

	if (timeDelta.m_eTimeType != TypeDelta) {
		_EXCEPTIONT("Argument to AddTime() is not a TimeDelta");
	}

	m_iYear        -= timeDelta.m_iYear;
	m_iMonth       -= timeDelta.m_iMonth;
	m_iDay         -= timeDelta.m_iDay;
	m_iSecond      -= timeDelta.m_iSecond;
	m_iMicroSecond -= timeDelta.m_iMicroSecond;

	NormalizeTime();
}

///////////////////////////////////////////////////////////////////////////////

void Time::RoundToNearestSecond() {

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("RoundToNearestMinute() only implemented for fixed times");
	}

	if (m_iMicroSecond < 500000) {
		m_iMicroSecond = 0;
	} else {
		m_iSecond++;
		m_iMicroSecond = 0;
		NormalizeTime();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Time::RoundToNearestMinute() {

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("RoundToNearestMinute() only implemented for fixed times");
	}

	int iMinute = m_iSecond / 60;
	int iSecond0to59 = m_iSecond % 60;

	if (iSecond0to59 < 30) {
		m_iSecond = iMinute * 60;
		m_iMicroSecond = 0;

	} else if (m_iSecond >= 30) {
		m_iSecond = (iMinute + 1) * 60;
		m_iMicroSecond = 0;
		NormalizeTime();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Time::RoundToNearestHour() {

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("RoundToNearestMinute() only implemented for fixed times");
	}

	int iHour = m_iSecond / 3600;
	int iSecond0to3599 = m_iSecond % 3600;

	if (iSecond0to3599 < 1800) {
		m_iSecond = iHour * 3600;
		m_iMicroSecond = 0;

	} else if (m_iSecond >= 1800) {
		m_iSecond = (iHour + 1) * 3600;
		m_iMicroSecond = 0;
		NormalizeTime();
	}
}

///////////////////////////////////////////////////////////////////////////////

int Time::DayNumber() const {

	// Based on https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
	// but modified since m_iMonth and m_iDay are zero-indexed
	if ((m_eCalendarType == CalendarNoLeap) ||
	    (m_eCalendarType == Calendar365Day)
	) {
		int nM = (m_iMonth + 10) % 12;
		int nY = m_iYear - nM/10;
		int nDay = 365 * nY + (nM * 306 + 5) / 10 + m_iDay;

		return nDay;

	} else if (
		(m_eCalendarType == CalendarStandard) ||
		(m_eCalendarType == CalendarGregorian)
	) {
		int nM = (m_iMonth + 10) % 12;
		int nY = m_iYear - nM/10;
		int nDay = 365 * nY + nY / 4 - nY / 100 + nY / 400
			+ (nM * 306 + 5) / 10 + m_iDay;

		return nDay;

	} else if (m_eCalendarType == Calendar360Day) {
		return (360 * m_iYear + 30 * m_iMonth + m_iDay);

	} else {
		_EXCEPTIONT("Not implemented");
	}
}

///////////////////////////////////////////////////////////////////////////////

bool Time::IsLeapDay() const {
	if ((m_eCalendarType == CalendarUnknown) ||
	    (m_eCalendarType == CalendarNone) ||
	    (m_eCalendarType == CalendarNoLeap) ||
	    (m_eCalendarType == Calendar360Day) ||
	    (m_eCalendarType == Calendar365Day)
	) {
		return false;
	}

	if ((m_iMonth == 1) && (m_iDay == 28)) {
		return true;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

bool Time::IsLeapYear() const {
	if ((m_eCalendarType == CalendarUnknown) ||
	    (m_eCalendarType == CalendarNone) ||
	    (m_eCalendarType == CalendarNoLeap) ||
	    (m_eCalendarType == Calendar360Day) ||
	    (m_eCalendarType == Calendar365Day)
	) {
		return false;
	}
	if ((m_iYear % 4) == 0) {
		if (((m_iYear % 100) == 0) && ((m_iYear % 400) != 0)) {
			return false;
		}
		return true;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

double Time::AsSeconds() const {
	if (m_eTimeType != TypeDelta) {
		_EXCEPTIONT("AsSeconds() only valid for operands of Time::TypeDelta");
	}
	if ((m_iYear != 0) || (m_iMonth != 0)) {
		_EXCEPTIONT("AsSeconds() not supported for non-zero year / month as result is dependent on calendar");
	}

	return (
		86400.0 * static_cast<double>(m_iDay)
		+ static_cast<double>(m_iSecond)
		+ 1.0e-6 * static_cast<double>(m_iMicroSecond));
}

///////////////////////////////////////////////////////////////////////////////

double Time::operator-(const Time & time) const {
	return -DeltaSeconds(time);
}

///////////////////////////////////////////////////////////////////////////////

double Time::DeltaSeconds(const Time & time) const {

	if ((m_eTimeType != TypeFixed) || (time.m_eTimeType != TypeFixed)) {
		_EXCEPTIONT("DeltaSeconds() only valid for operands of Time::TypeFixed");
	}
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT("Incompatible calendar types");
	}

	int nDayNumber1 = DayNumber();
	int nDayNumber2 = time.DayNumber();

	return static_cast<double>(nDayNumber2 - nDayNumber1) * 86400.0
		+ static_cast<double>(time.m_iSecond - m_iSecond)
		+ static_cast<double>(time.m_iMicroSecond - m_iMicroSecond) / 1.0e6;
}

///////////////////////////////////////////////////////////////////////////////

double Time::DeltaMinutes(const Time & time) const {

	if ((m_eTimeType != TypeFixed) || (time.m_eTimeType != TypeFixed)) {
		_EXCEPTIONT("DeltaMinutes() only valid for operands of Time::TypeFixed");
	}
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT("Incompatible calendar types");
	}

	int nDayNumber1 = DayNumber();
	int nDayNumber2 = time.DayNumber();

	return static_cast<double>(nDayNumber2 - nDayNumber1) * 1440.0
		+ static_cast<double>(time.m_iSecond - m_iSecond) / 60.0
		+ static_cast<double>(time.m_iMicroSecond - m_iMicroSecond) / 6.0e7;
}

///////////////////////////////////////////////////////////////////////////////

double Time::DeltaHours(const Time & time) const {

	if ((m_eTimeType != TypeFixed) || (time.m_eTimeType != TypeFixed)) {
		_EXCEPTIONT("DeltaHours() only valid for operands of Time::TypeFixed");
	}
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT("Incompatible calendar types");
	}

	int nDayNumber1 = DayNumber();
	int nDayNumber2 = time.DayNumber();

	return static_cast<double>(nDayNumber2 - nDayNumber1) * 24.0
		+ static_cast<double>(time.m_iSecond - m_iSecond) / 3600.0
		+ static_cast<double>(time.m_iMicroSecond - m_iMicroSecond) / 3.6e9;
}

///////////////////////////////////////////////////////////////////////////////

double Time::DeltaDays(const Time & time) const {

	if ((m_eTimeType != TypeFixed) || (time.m_eTimeType != TypeFixed)) {
		_EXCEPTIONT("DeltaDays() only valid for operands of Time::TypeFixed");
	}
	if (m_eCalendarType != time.m_eCalendarType) {
		_EXCEPTIONT("Incompatible calendar types");
	}

	int nDayNumber1 = DayNumber();
	int nDayNumber2 = time.DayNumber();

	return static_cast<double>(nDayNumber2 - nDayNumber1)
		+ static_cast<double>(time.m_iSecond - m_iSecond) / 86400.0
		+ static_cast<double>(time.m_iMicroSecond - m_iMicroSecond) / 8.64e10;
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToDateString() const {
	char szBuffer[100];

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("ToDateString() only valid for Time::TypeFixed");
	}

	snprintf(szBuffer, 100, "%04i-%02i-%02i",
		m_iYear, m_iMonth+1, m_iDay+1);

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToShortString() const {
	char szBuffer[100];

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("ToShortString() only valid for Time::TypeFixed");
	}

	snprintf(szBuffer, 100, "%04i-%02i-%02i-%05i",
		m_iYear, m_iMonth+1, m_iDay+1, m_iSecond);

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToLongString() const {
	char szBuffer[100];

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("ToLongString() only valid for Time::TypeFixed");
	}

	int iReindexedMonth = m_iMonth;
	int iReindexedDay = m_iDay;

	if ((m_eCalendarType != CalendarNone) && (m_eCalendarType != CalendarUnknown)) {
		iReindexedMonth++;
		iReindexedDay++;
	}

	snprintf(szBuffer, 100, "%04i-%02i-%02i-%05i-%06i",
		m_iYear, iReindexedMonth, iReindexedDay, m_iSecond, m_iMicroSecond);

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToString() const {
	char szBuffer[100];

	if (m_eTimeType != TypeFixed) {
		_EXCEPTIONT("ToString() only valid for Time::TypeFixed");
	}

	int iReindexedMonth = m_iMonth;
	int iReindexedDay = m_iDay;

	if ((m_eCalendarType != CalendarNone) && (m_eCalendarType != CalendarUnknown)) {
		iReindexedMonth++;
		iReindexedDay++;
	}

	if (m_iMicroSecond == 0) {
		snprintf(szBuffer, 100, "%04i-%02i-%02i %02i:%02i:%02i",
			m_iYear,
			iReindexedMonth,
			iReindexedDay,
			m_iSecond / 3600,
			(m_iSecond % 3600) / 60,
			(m_iSecond % 60));

	} else {
		snprintf(szBuffer, 100, "%04i-%02i-%02i %02i:%02i:%02i.%06i",
			m_iYear,
			iReindexedMonth,
			iReindexedDay,
			m_iSecond / 3600,
			(m_iSecond % 3600) / 60,
			(m_iSecond % 60),
			m_iMicroSecond);
	}

	return std::string(szBuffer);
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::ToFreeString() const {
	std::string strFreeString;

	char szBuffer[100];

	if (m_iYear != 0) {
		snprintf(szBuffer, 100, "%iy", m_iYear);
		strFreeString += szBuffer;
	}
	if (m_iMonth != 0) {
		snprintf(szBuffer, 100, "%iM", m_iMonth);
		strFreeString += szBuffer;
	}
	if (m_iDay != 0) {
		snprintf(szBuffer, 100, "%id", m_iDay);
		strFreeString += szBuffer;
	}
	if (m_iSecond != 0) {
		snprintf(szBuffer, 100, "%is", m_iSecond);
		strFreeString += szBuffer;
	}
	if (m_iMicroSecond != 0) {
		snprintf(szBuffer, 100, "%iu", m_iMicroSecond);
		strFreeString += szBuffer;
	}

	return strFreeString;
}

///////////////////////////////////////////////////////////////////////////////

void Time::FromLongString(const std::string & strLongTime) {
	if ((m_eCalendarType == CalendarNone) && (m_eCalendarType == CalendarUnknown)) {
		_EXCEPTIONT("Not implemented");
	}
	if (strLongTime.length() != 23) {
		_EXCEPTIONT("Invalid Time LongString");
	}
	m_iYear = atoi(strLongTime.c_str());
	m_iMonth = atoi(strLongTime.c_str()+5) - 1;
	m_iDay = atoi(strLongTime.c_str()+8) - 1;
	m_iSecond = atoi(strLongTime.c_str()+11);
	m_iMicroSecond = atoi(strLongTime.c_str()+17);

	VerifyTime();
}

///////////////////////////////////////////////////////////////////////////////

void Time::FromFormattedString(
	const std::string & strFormattedTime
) {
	const char * szYear = NULL;
	const char * szMonth = NULL;
	const char * szDay = NULL;
	const char * szHour = NULL;
	const char * szMinute = NULL;
	const char * szSecond = NULL;
	const char * szMicroSecond = NULL;

	// Reset values
	m_iYear = 0;
	m_iMonth = 0;
	m_iDay = 0;
	m_iSecond = 0;
	m_iMicroSecond = 0;

	// Check for length zero string
	if (strFormattedTime.length() == 0) {
		return;
	}

	// Current state
	enum FormatState {
		FormatState_Date,
		FormatState_Time,
		FormatState_Free
	};

	FormatState state = FormatState_Date;

	// Loop over all characters in string
	int i = 0;
	int j = 0;
	for (; i < strFormattedTime.length(); i++) {

		// Digit (ignore)
		if ((strFormattedTime[i] >= '0') &&
			(strFormattedTime[i] <= '9')
		) {
			continue;
		}

		// Record in Date format (yyyy-MM-dd-sssss)
		if ((strFormattedTime[i] == '-') ||
		    (strFormattedTime[i] == ' ') ||
			(strFormattedTime[i] == 'T') ||
			(strFormattedTime[i] == '_')
		) {
			if (state != FormatState_Date) {
				break;
/*
				_EXCEPTION1(
					"Malformed Time string (%s): "
						"Cannot return to Date format",
						strFormattedTime.c_str());
*/
			}

			if (szYear == NULL) {
				szYear = &(strFormattedTime[j]);
				j = i + 1;

			} else if (szMonth == NULL) {
				szMonth = &(strFormattedTime[j]);
				j = i + 1;

			} else if (szDay == NULL) {
				szDay = &(strFormattedTime[j]);
				j = i + 1;

			} else {
				_EXCEPTION1("Malformed Time string (%s): "
					"Unexpected \'-\'",
					strFormattedTime.c_str());
			}

		// Record in Time format (hh:mm:ss.uuuuuu)
		} else if (
			(strFormattedTime[i] == ':') ||
			(strFormattedTime[i] == '.')
		) {
			if (state == FormatState_Free) {
				_EXCEPTION1(
					"Malformed Time string (%s): "
						"Cannot mix formatting types",
						strFormattedTime.c_str());
			}

			state = FormatState_Time;

			if (szHour == NULL) {
				szHour = &(strFormattedTime[j]);
				j = i + 1;

			} else if (szMinute == NULL) {
				szMinute = &(strFormattedTime[j]);
				j = i + 1;

			} else if (szSecond == NULL) {
				szSecond = &(strFormattedTime[j]);
				j = i + 1;

			} else {
				_EXCEPTION1(
					"Malformed Time string (%s): "
						"Unexpected \':\' or \'.\'",
						strFormattedTime.c_str());
			}

		// Record in Free format (##y##M##d##m##s##u)
		} else {
			if (state == FormatState_Time) {
				_EXCEPTION1(
					"Malformed Time string (%s): "
						"Cannot mix Free and Time format",
						strFormattedTime.c_str());
			}
			if ((state == FormatState_Date) && (j != 0)) {
				_EXCEPTION1(
					"Malformed Time string (%s): "
						"Cannot mix Free and Date format",
						strFormattedTime.c_str());
			}

			state = FormatState_Free;

			// Record year
			if (strFormattedTime[i] == 'y') {
				if (szYear == NULL) {
					szYear = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'y\'",
						strFormattedTime.c_str());
				}

			// Record month
			} else if (strFormattedTime[i] == 'M') {
				if (szMonth == NULL) {
					szMonth = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'M\'",
						strFormattedTime.c_str());
				}

			// Record day
			} else if (strFormattedTime[i] == 'd') {
				if (szDay == NULL) {
					szDay = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'d\'",
						strFormattedTime.c_str());
				}

			// Record hour
			} else if (strFormattedTime[i] == 'h') {
				if (szHour == NULL) {
					szHour = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'h\'",
						strFormattedTime.c_str());
				}

			// Record minute
			} else if (strFormattedTime[i] == 'm') {
				if (szMinute == NULL) {
					szMinute = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'m\'",
						strFormattedTime.c_str());
				}

			// Record second 
			} else if (strFormattedTime[i] == 's') {
				if (szSecond == NULL) {
					szSecond = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'s\'",
						strFormattedTime.c_str());
				}

			// Record microsecond 
			} else if (strFormattedTime[i] == 'u') {
				if (szMicroSecond == NULL) {
					szMicroSecond = &(strFormattedTime[j]);
				} else {
					_EXCEPTION1("Malformed Time string (%s): "
						"Two instances of \'s\'",
						strFormattedTime.c_str());
				}

			// Unexpected character
			} else {
				_EXCEPTION2("Malformed Time string (%s): "
					"Unexpected character \'%c\'",
					strFormattedTime.c_str(),
					strFormattedTime[i]);
			}

			// Increment location pointer
			j = i + 1;
		}
	}

	// Finalize formatting
	if (i != j) {

		// Date format
		if (state == FormatState_Date) {
			if (szYear == NULL) {
				_EXCEPTION1("Malformed time string (%s): "
					"No formatting characters detected",
					strFormattedTime.c_str());

			} else if (szMonth == NULL) {
				szMonth = &(strFormattedTime[j]);

			} else if (szDay == NULL) {
				szDay = &(strFormattedTime[j]);

			} else if (szSecond == NULL) {
				szSecond = &(strFormattedTime[j]);

			} else {
				_EXCEPTION1("Malformed time string (%s): "
					"Extraneous value at end",
					strFormattedTime.c_str());
			}

		// Time format
		} else if (state == FormatState_Time) {
			if (szMinute == NULL) {
				szMinute = &(strFormattedTime[j]);

			} else if (szSecond == NULL) {
				szSecond = &(strFormattedTime[j]);

			} else if (szMicroSecond == NULL) {
				szMicroSecond = &(strFormattedTime[j]);

			} else {
				_EXCEPTION1("Malformed time string (%s): "
					"Extraneous value at end",
					strFormattedTime.c_str());
			}

		// Free format
		} else {
			_EXCEPTION1("Malformed time string (%s): "
				"Dangling values not allowed in Free formatting",
				strFormattedTime.c_str());
		} 
	}

	// Type
	if ((state == FormatState_Date) || (state == FormatState_Time)) {
		m_eTimeType = TypeFixed;
	} else {
		m_eTimeType = TypeDelta;
	}

	// Convert
	m_iSecond = 0;

	if (szYear != NULL) {
		m_iYear = atoi(szYear);
	}
	if (szMonth != NULL) {
		if (state == FormatState_Free) {
			m_iMonth = atoi(szMonth);
		} else {
			m_iMonth = atoi(szMonth) - 1;
		}
	}
	if (szDay != NULL) {
		if (state == FormatState_Free) {
			m_iDay = atoi(szDay);
		} else {
			m_iDay = atoi(szDay) - 1;
		}
	}
	if (szHour != NULL) {
		m_iSecond += 3600 * atoi(szHour);
	}
	if (szMinute != NULL) {
		m_iSecond += 60 * atoi(szMinute);
	}
	if (szSecond != NULL) {
		m_iSecond += atoi(szSecond);
	}
	if (szMicroSecond != NULL) {
		m_iMicroSecond = atoi(szMicroSecond);
	}

	// Normalize
	NormalizeTime();
}

///////////////////////////////////////////////////////////////////////////////

void Time::FromCFCompliantUnitsOffsetInt(
	const std::string & strFormattedTime,
	int nOffset
) {
	// Time format is "days since ..."
	if ((strFormattedTime.length() >= 11) &&
	    (strncmp(strFormattedTime.c_str(), "days since ", 11) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(11);
		FromFormattedString(strSubStr);

		AddDays(nOffset);

	// Time format is "hours since ..."
	} else if (
	    (strFormattedTime.length() >= 12) &&
	    (strncmp(strFormattedTime.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(12);
		FromFormattedString(strSubStr);

		AddHours(nOffset);

	// Time format is "minutes since ..."
	} else if (
	    (strFormattedTime.length() >= 14) &&
	    (strncmp(strFormattedTime.c_str(), "minutes since ", 14) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(14);
		FromFormattedString(strSubStr);

		AddMinutes(nOffset);

	// Time format is "seconds since ..."
	} else if (
	    (strFormattedTime.length() >= 14) &&
	    (strncmp(strFormattedTime.c_str(), "seconds since ", 14) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(14);
		FromFormattedString(strSubStr);

		AddSeconds(nOffset);

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
}

///////////////////////////////////////////////////////////////////////////////

void Time::FromCFCompliantUnitsOffsetDouble(
	const std::string & strFormattedTime,
	double dOffset
) {
	// Time format is "days since ..."
	if ((strFormattedTime.length() >= 11) &&
	    (strncmp(strFormattedTime.c_str(), "days since ", 11) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(11);
		FromFormattedString(strSubStr);

		int nDays = static_cast<int>(dOffset);
		AddDays(nDays);

		int nSeconds = static_cast<int>(fmod(dOffset, 1.0) * 86400.0);
		AddSeconds(nSeconds);

	// Time format is "hours since ..."
	} else if (
	    (strFormattedTime.length() >= 12) &&
	    (strncmp(strFormattedTime.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(12);
		FromFormattedString(strSubStr);

		int nHours = static_cast<int>(dOffset);
		AddHours(nHours);

		int nSeconds = static_cast<int>(fmod(dOffset, 1.0) * 3600.0);
		AddSeconds(nSeconds);

	// Time format is "minutes since ..."
	} else if (
	    (strFormattedTime.length() >= 14) &&
	    (strncmp(strFormattedTime.c_str(), "minutes since ", 14) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(14);
		FromFormattedString(strSubStr);

		int nMinutes = static_cast<int>(dOffset);
		AddMinutes(nMinutes);

		int nSeconds = static_cast<int>(fmod(dOffset, 1.0) * 60.0);
		AddSeconds(nSeconds);

	// Time format is "seconds since ..."
	} else if (
	    (strFormattedTime.length() >= 14) &&
	    (strncmp(strFormattedTime.c_str(), "seconds since ", 14) == 0)
	) {
		std::string strSubStr = strFormattedTime.substr(14);
		FromFormattedString(strSubStr);

		int nDays = static_cast<int>(dOffset / 86400.0);
		AddDays(nDays);

		int nSeconds = static_cast<int>(fmod(dOffset, 86400.0));
		AddSeconds(nSeconds);

	// Time format is "day as %Y%m%d.%f"
	} else if (
	    (strFormattedTime.length() >= 16) &&
	    (strncmp(strFormattedTime.c_str(), "day as %%Y%%m%%d.%%f", 16) == 0)
	) {
		m_iSecond = static_cast<int>(fmod(dOffset, 1.0) * 86400.0);

		int nRem = static_cast<int>(dOffset);

		m_iDay = (nRem % 100) - 1;

		m_iMonth = ((nRem / 100) % 100) - 1;

		m_iYear = (nRem / 10000);

		NormalizeTime();

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
}

///////////////////////////////////////////////////////////////////////////////

double Time::GetCFCompliantUnitsOffsetDouble(
	const std::string & strFormattedTime
) const {
	// Time format is "days since ..."
	if ((strFormattedTime.length() >= 11) &&
	    (strncmp(strFormattedTime.c_str(), "days since ", 11) == 0)
	) {
		Time timeBuf(m_eCalendarType);
		std::string strSubStr = strFormattedTime.substr(11);
		timeBuf.FromFormattedString(strSubStr);

		return timeBuf.DeltaDays(*this);

	// Time format is "hours since ..."
	} else if (
	    (strFormattedTime.length() >= 12) &&
	    (strncmp(strFormattedTime.c_str(), "hours since ", 12) == 0)
	) {
		Time timeBuf(m_eCalendarType);
		std::string strSubStr = strFormattedTime.substr(12);
		timeBuf.FromFormattedString(strSubStr);

		return timeBuf.DeltaHours(*this);

	// Time format is "minutes since ..."
	} else if (
	    (strFormattedTime.length() >= 14) &&
	    (strncmp(strFormattedTime.c_str(), "minutes since ", 14) == 0)
	) {
		Time timeBuf(m_eCalendarType);
		std::string strSubStr = strFormattedTime.substr(14);
		timeBuf.FromFormattedString(strSubStr);

		return timeBuf.DeltaMinutes(*this);

	// Time format is "seconds since ..."
	} else if (
	    (strFormattedTime.length() >= 14) &&
	    (strncmp(strFormattedTime.c_str(), "seconds since ", 14) == 0)
	) {
		Time timeBuf(m_eCalendarType);
		std::string strSubStr = strFormattedTime.substr(14);
		timeBuf.FromFormattedString(strSubStr);

		return timeBuf.DeltaSeconds(*this);

	// Time format is "day as %Y%m%d.%f"
	} else if (
	    (strFormattedTime.length() >= 16) &&
	    (strncmp(strFormattedTime.c_str(), "day as %%Y%%m%%d.%%f", 16) == 0)
	) {
		double dOffset =
			static_cast<double>(m_iYear * 10000)
			+ static_cast<double>((m_iMonth+1) * 100)
			+ static_cast<double>(m_iDay+1)
			+ (static_cast<double>(m_iSecond) / 86400.0);

		return dOffset;

	} else {
		_EXCEPTION1("Unknown \"time::units\" format \"%s\"", strFormattedTime.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

std::string Time::GetCalendarName() const {
	if (m_eCalendarType == CalendarNone) {
		return std::string("none");
	} else if (m_eCalendarType == CalendarNoLeap) {
		return std::string("noleap");
	} else if (m_eCalendarType == CalendarStandard) {
		return std::string("standard");
	} else if (m_eCalendarType == CalendarGregorian) {
		return std::string("gregorian");
	} else if (m_eCalendarType == Calendar360Day) {
		return std::string("360_day");
	} else if (m_eCalendarType == Calendar365Day) {
		return std::string("365_day");
	} else {
		_EXCEPTIONT("Invalid CalendarType");
	}
}

///////////////////////////////////////////////////////////////////////////////

