///////////////////////////////////////////////////////////////////////////////
///
///	\file    Exception.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<summary>
///		This file provides functionality for formatted Exceptions.
///	</summary>
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

///////////////////////////////////////////////////////////////////////////////

#include <cstdarg>

///////////////////////////////////////////////////////////////////////////////

#define _EXCEPTION() \
throw Exception(__FILE__, __LINE__)

#define _EXCEPTIONT(text) \
throw Exception(__FILE__, __LINE__, text)

#define _EXCEPTION1(text, var1) \
throw Exception(__FILE__, __LINE__, text, var1)

#define _EXCEPTION2(text, var1, var2) \
throw Exception(__FILE__, __LINE__, text, var1, var2)

#define _EXCEPTION3(text, var1, var2, var3) \
throw Exception(__FILE__, __LINE__, text, var1, var2, var3)

#define _EXCEPTION4(text, var1, var2, var3, var4) \
throw Exception(__FILE__, __LINE__, text, var1, var2, var3, var4)

#define _EXCEPTION5(text, var1, var2, var3, var4, var5) \
throw Exception(__FILE__, __LINE__, text, var1, var2, var3, var4, var5)

#define _EXCEPTION6(text, var1, var2, var3, var4, var5, var6) \
throw Exception(__FILE__, __LINE__, text, var1, var2, var3, var4, var5, var6)

#define _EXCEPTION7(text, var1, var2, var3, var4, var5, var6, var7) \
throw Exception(__FILE__, __LINE__, text, var1, var2, var3, var4, var5, var6, var7)

#define _EXCEPTION8(text, var1, var2, var3, var4, var5, var6, var7, var8) \
throw Exception(__FILE__, __LINE__, text, var1, var2, var3, var4, var5, var6, var7, var8)

///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <cstdio>
#include <stdarg.h>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An Exception is a formatted error message that is generated from a
///		throw directive.  This class is automatically generated when using
///		the _EXCEPTION macros.
///	</summary>
class Exception {
	
	public:
		///	<summary>
		///		Maximum buffer size for exception strings.
		///	</summary>
		static const int ExceptionBufferSize = 1024;

	public:
		///	<summary>
		///		Generic constructor.
		///	</summary>
		Exception(
			const char * szFile,
			unsigned int uiLine
		) :
			m_strText("General exception"),
			m_strFile(szFile),
			m_uiLine(uiLine)
		{ }
	
		///	<summary>
		///		Constructor with text and variables.
		///	</summary>
		Exception(
			const char * szFile,
			unsigned int uiLine,
			const char * szText,
			...
		) :
			m_strFile(szFile),
			m_uiLine(uiLine)
		{
			char szBuffer[ExceptionBufferSize];

			va_list arguments;

			// Initialize the argument list
			va_start(arguments, szText);

			// Write to string
			vsprintf(szBuffer, szText, arguments);

			m_strText = szBuffer;

			// Cleans up the argument list
			va_end(arguments);
		}
	
	public:
		///	<summary>
		///		Get a string representation of this exception.
		///	</summary>
		std::string ToString() const {
			std::string strReturn;

			char szBuffer[128];

			// Preamble
			sprintf(szBuffer, "EXCEPTION (");
			strReturn.append(szBuffer);

			// File name
			strReturn.append(m_strFile);

			// Line number
			sprintf(szBuffer, ", Line %u) ", m_uiLine);
			strReturn.append(szBuffer);

			// Text
			strReturn.append(m_strText);

			return strReturn;
		}

	private:
		///	<summary>
		///		A string denoting the error in question.
		///	</summary>
		std::string m_strText;
	
		///	<summary>
		///		A string containing the filename where the exception occurred.
		///	</summary>
		std::string m_strFile;
	
		///	<summary>
		///		A constant containing the line number where the exception
		///		occurred.
		///	</summary>
		unsigned int m_uiLine;
};

///////////////////////////////////////////////////////////////////////////////

#endif

