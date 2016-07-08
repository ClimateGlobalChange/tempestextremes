///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataVector.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _DATAVECTOR_H_
#define _DATAVECTOR_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data vector is a datatype that stores data in a 1D structure.
///		Arithmetic operations are supported for this datatype.
///	</summary>
///	<warning>
///		Memory is allocated using malloc and deallocated using free, so no
///		calls will be made to the constructor or destructor of DataType.  This
///		class is primarily designed to efficiently handle primitive types.
///	</warning>

template <typename DataType>
class DataVector {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		DataVector() :
			m_sRows(0),
			m_data(NULL)
		{ }

		///	<summary>
		///		Constructor.
		///	</summary>
		DataVector(
			unsigned int sRows
		) :
			m_sRows(0),
			m_data(NULL)
		{
			Initialize(sRows);
		}

		///	<summary>
		///		Copy constructor.
		///	</summary>
		DataVector(
			const DataVector<DataType> & dv
		) :
			m_sRows(0),
			m_data(NULL)
		{
			Assign(dv);
		}

		///	<summary>
		///		Destructor.
		///	</summary>
		virtual ~DataVector() {
			if (m_data != NULL) {
				free(reinterpret_cast<void*>(m_data));
			}
		}

	public:
		///	<summary>
		///		Determine if this DataVector is initialized.
		///	</summary>
		bool IsInitialized() const {
			if (m_data == NULL) {
				return false;
			} else {
				return true;
			}
		}

	public:
		///	<summary>
		///		Deallocate data for this object.
		///	</summary>
		void Deinitialize() {
			if (m_data != NULL) {
				free(reinterpret_cast<void*>(m_data));
			}

			m_data = NULL;
			m_sRows = 0;
		}

		///	<summary>
		///		Allocate data for this object.
		///	</summary>
		void Initialize(
			unsigned int sRows,
			bool fAutoZero = true
		) {

			// Check for zero size
			if (sRows == 0) {
				Deinitialize();
				return;
			}

			// No need to reallocate memory if this vector already has
			// the correct dimensions.
			if (m_sRows == sRows) {
				// Auto zero
				if (fAutoZero) {
					Zero();
				}

				return;
			}

			// Deinitialize existing content
			Deinitialize();

			// Allocate memory
			m_data = reinterpret_cast<DataType*>(
				malloc(sRows * sizeof(DataType)));

			if (m_data == NULL) {
				_EXCEPTIONT("Out of memory.");
			}

			// Assign dimensions
			m_sRows = sRows;

			// Auto zero
			if (fAutoZero) {
				Zero();
			}
		}

	public:
		///	<summary>
		///		Assignment operator.
		///	</summary>
		void Assign(const DataVector<DataType> & dv) {

			// Check initialization status
			if (!dv.IsInitialized()) {
				Deinitialize();
				return;
			}

			// Allocate memory
			Initialize(dv.m_sRows, false);

			// Copy data
			memcpy(m_data, dv.m_data, m_sRows * sizeof(DataType));
		}

		///	<summary>
		///		Assignment operator.
		///	</summary>
		DataVector & operator= (const DataVector<DataType> & dv) {
			Assign(dv);
			return (*this);
		}

		///	<summary>
		///		Zero the data content of this object.
		///	</summary>
		void Zero() {

			// Check initialization status
			if (!IsInitialized()) {
				_EXCEPTIONT(
					"Attempted operation on uninitialized DataVector.");
			}

			// Set content to zero
			memset(m_data, 0, m_sRows * sizeof(DataType));
		}

	public:
		///	<summary>
		///		Get the number of rows in this matrix.
		///	</summary>
		inline unsigned int GetRows() const {
			return m_sRows;
		}

	public:
		///	<summary>
		///		Cast to an array.
		///	</summary>
		inline operator DataType*() {
			return m_data;
		}

		///	<summary>
		///		Cast to an array.
		///	</summary>
		inline operator const DataType*() const {
			return m_data;
		}

	private:
		///	<summary>
		///		The number of rows in this matrix.
		///	</summary>
		unsigned int m_sRows;

		///	<summary>
		///		A pointer to the data associated with this matrix.
		///	</summary>
		DataType* m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

