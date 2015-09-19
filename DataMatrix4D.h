///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataMatrix4D.h
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

#ifndef _DATAMATRIX4D_H_
#define _DATAMATRIX4D_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

#include <iostream>
#include <cstdlib>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A 4D data matrix is a datatype that stores data in a 4D structure.
///		Arithmatic operations are not supported for this datatype.
///	</summary>

template <typename DataType>
class DataMatrix4D {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		DataMatrix4D() :
			m_data(NULL)
		{
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
			m_sSize[3] = 0;
		}

		///	<summary>
		///		Constructor.
		///	</summary>
		DataMatrix4D(
			unsigned int sSize0,
			unsigned int sSize1,
			unsigned int sSize2,
			unsigned int sSize3
		) :
			m_data(NULL)
		{
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
			m_sSize[3] = 0;

			Initialize(sSize0, sSize1, sSize2, sSize3);
		}

		///	<summary>
		///		Copy constructor.
		///	</summary>
		DataMatrix4D(const DataMatrix4D<DataType> & dm) 
			: m_data(NULL)
		{
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
			m_sSize[3] = 0;

			Assign(dm);
		}

		///	<summary>
		///		Destructor.
		///	</summary>
		virtual ~DataMatrix4D() {
			if (m_data != NULL) {
				free(reinterpret_cast<void*>(m_data));
			}
		}

	public:
		///	<summary>
		///		Determine if this DataMatrix4D is initialized.
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
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
			m_sSize[3] = 0;
		}

		///	<summary>
		///		Allocate data for this object.
		///	</summary>
		void Initialize(
			unsigned int sSize0,
			unsigned int sSize1,
			unsigned int sSize2,
			unsigned int sSize3,
			bool fAutoZero = true
		) {
			unsigned int sI;
			unsigned int sJ;
			unsigned int sK;

			// Check for zero size
			if ((sSize0 == 0) || (sSize1 == 0) ||
				(sSize2 == 0) || (sSize3 == 0)
			) {
				Deinitialize();
				return;
			}

			// No need to reallocate memory if this matrix already has
			// the correct dimensions.
			if ((m_sSize[0] == sSize0) &&
				(m_sSize[1] == sSize1) &&
				(m_sSize[2] == sSize2) &&
				(m_sSize[3] == sSize3)
			) {
				// Auto zero
				if (fAutoZero) {
					Zero();
				}

				return;
			}

			// Deinitialize existing content
			Deinitialize();

			// Calculate the footprint in each direction
			unsigned int sDim0PtrFootprint = sSize0 * sizeof(DataType ***);
			unsigned int sDim1PtrFootprint = sSize1 * sizeof(DataType **);
			unsigned int sDim2PtrFootprint = sSize2 * sizeof(DataType *);

			unsigned int sDim3DataFootprint = sSize3 * sizeof(DataType);

			unsigned int sTotalSize =
				sSize0 * sizeof(DataType ***) +
				sSize0 * sSize1 * sizeof(DataType **) +
				sSize0 * sSize1 * sSize2 * sizeof(DataType *) +
				sSize0 * sSize1 * sSize2 * sSize3 * sizeof(DataType);

			// Allocate memory
			char *rawdata = reinterpret_cast<char*>(malloc(sTotalSize));

			if (rawdata == NULL) {
				_EXCEPTIONT("Out of memory.");
			}

			// Assign memory pointers
			char *pDim0Ptrs = rawdata;

			char *pDim1Ptrs = pDim0Ptrs + sDim0PtrFootprint;

			char *pDim2Ptrs = pDim1Ptrs + sSize0 * sDim1PtrFootprint;

			char *pDataStart = pDim2Ptrs + sSize0 * sSize1 * sDim2PtrFootprint;

			m_data = reinterpret_cast<DataType****>(rawdata);

			for (sI = 0; sI < sSize0; sI++) {
				m_data[sI] = reinterpret_cast<DataType***>(
					pDim1Ptrs + sI * sDim1PtrFootprint);

			for (sJ = 0; sJ < sSize1; sJ++) {
				m_data[sI][sJ] = reinterpret_cast<DataType**>(
					pDim2Ptrs + (sI * sSize1 + sJ) * sDim2PtrFootprint);

			for (sK = 0; sK < sSize2; sK++) {
				m_data[sI][sJ][sK] = reinterpret_cast<DataType*>(
					pDataStart +
					((sI * sSize1 + sJ) * sSize2 + sK)
						* sDim3DataFootprint);
			}
			}
			}

			// Assign dimensions
			m_sSize[0] = sSize0;
			m_sSize[1] = sSize1;
			m_sSize[2] = sSize2;
			m_sSize[3] = sSize3;

			// Auto zero
			if (fAutoZero) {
				Zero();
			}
		}

	public:
		///	<summary>
		///		Assignment operator.
		///	</summary>
		void Assign(const DataMatrix4D<DataType> & dm) {

			// Check initialization status
			if (!dm.IsInitialized()) {
				Deinitialize();
				return;
			}

			// Allocate memory
			Initialize(
				dm.m_sSize[0],
				dm.m_sSize[1],
				dm.m_sSize[2],
				dm.m_sSize[3],
				false);

			// Copy data
			unsigned int nOffset =
				m_sSize[0] * sizeof(DataType ***)
				+ m_sSize[0] * m_sSize[1] * sizeof(DataType **)
				+ m_sSize[0] * m_sSize[1] * m_sSize[2] * sizeof(DataType *);

			unsigned int nDataSize = sizeof(DataType) *
				m_sSize[0] * m_sSize[1] * m_sSize[2] * m_sSize[3];

			memcpy(
				reinterpret_cast<char*>(m_data) + nOffset,
				reinterpret_cast<char*>(dm.m_data) + nOffset,
				nDataSize
			);
		}

		///	<summary>
		///		Assignment operator.
		///	</summary>
		DataMatrix4D & operator= (const DataMatrix4D<DataType> & dm) {
			Assign(dm);
			return (*this);
		}

		///	<summary>
		///		Zero the data content of this object.
		///	</summary>
		void Zero() {

			// Check initialization status
			if (!IsInitialized()) {
				_EXCEPTIONT(
					"Attempted operation on uninitialized DataMatrix3D.");
			}

			// Set content to zero
			unsigned int nOffset =
				m_sSize[0] * sizeof(DataType ***)
				+ m_sSize[0] * m_sSize[1] * sizeof(DataType **)
				+ m_sSize[0] * m_sSize[1] * m_sSize[2] * sizeof(DataType *);

			unsigned int nDataSize = sizeof(DataType) *
				m_sSize[0] * m_sSize[1] * m_sSize[2] * m_sSize[3];

			memset(
				reinterpret_cast<char*>(m_data) + nOffset,
				0,
				nDataSize
			);
		}

	public:
		///	<summary>
		///		Get the number of elements along the specified dimension.
		///	</summary>
		inline unsigned int GetSize(int dim) const {
			return m_sSize[dim];
		}

		///	<summary>
		///		Get the total number of elements in this matrix.
		///	</summary>
		inline unsigned int GetTotalElements() const {
			return m_sSize[0] * m_sSize[1] * m_sSize[2] * m_sSize[3];
		}

	public:
		///	<summary>
		///		Cast to an array.
		///	</summary>
		inline operator DataType****() const {
			return m_data;
		}

	private:
		///	<summary>
		///		The number of elements in each dimension of this matrix.
		///	</summary>
		unsigned int m_sSize[4];

		///	<summary>
		///		A pointer to the data associated with this matrix.
		///	</summary>
		DataType**** m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

