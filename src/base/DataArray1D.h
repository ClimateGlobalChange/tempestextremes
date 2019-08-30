///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataArray1D.h
///	\author  Paul Ullrich
///	\version June 26, 2015
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

#ifndef _DATAARRAY1D_H_
#define _DATAARRAY1D_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

#include <cstdlib>
#include <cstring>

template <typename T>
class DataArray1D {

public:
	typedef T ValueType;

	///	<summary>
	///		Constructor.
	///	</summary>
	DataArray1D() :
		m_fOwnsData(true),
		m_sSize(0),
		m_data(NULL)
	{ }

	///	<summary>
	///		Constructor allowing specification of size.
	///	</summary>
	DataArray1D(
		size_t sSize,
		bool fAllocate = true
	) :
		m_fOwnsData(true),
		m_sSize(sSize),
		m_data(NULL)
	{
		if (fAllocate) {
			Allocate(sSize);
		}
	}

	///	<summary>
	///		Copy constructor.
	///	</summary>
	DataArray1D(const DataArray1D<T> & da) :
		m_fOwnsData(true),
		m_sSize(0),
		m_data(NULL)
	{
		Assign(da);
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~DataArray1D() {
		Detach();
	}

	///	<summary>
	///		Get the size of the data, in bytes.
	///	</summary>
	virtual size_t GetByteSize() const {

		// Verify data aligns on word boundaries
		if ((m_sSize * sizeof(T)) % sizeof(size_t) == 0) {
			return (m_sSize * sizeof(T));
		} else {
			return ((m_sSize * sizeof(T)) / sizeof(size_t) + 1)
				* sizeof(size_t);
		}
	}

	///	<summary>
	///		Allocate data in this DataArray1D.
	///	</summary>
	void Allocate(
		size_t sSize
	) {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting to Allocate() on attached DataArray1D");
		}

		Detach();

		if (sSize == 0) {
			m_sSize = 0;

			return;
		}
		if ((m_data == NULL) || (m_sSize != sSize)) {
			m_sSize = sSize;

			m_data = reinterpret_cast<T *>(malloc(GetByteSize()));

			if (m_data == NULL) {
				_EXCEPTION1("Failed malloc call (%lu bytes)", GetByteSize());
			}
		}

		Zero();
	}

	///	<summary>
	///		Set the number of rows in this DataArray1D.
	///	</summary>
	inline void SetSize(
		size_t sSize
	) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting to SetSize() on attached DataArray1D");
		}

		m_sSize = sSize;
	}

public:
	///	<summary>
	///		Determine if this DataChunk is attached to a data array.
	///	</summary>
	virtual bool IsAttached() const {
		return (m_data != NULL);
	}

	///	<summary>
	///		Attach this DataChunk to an array of pre-allocated data.
	///	</summary>
	virtual void AttachToData(void * ptr) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting AttachToData() on attached DataArray1D");
		}

		m_data = reinterpret_cast<T *>(ptr);
		m_fOwnsData = false;
	}

	///	<summary>
	///		Detach data from this DataChunk.
	///	</summary>
	virtual void Detach() {
		if ((m_fOwnsData) && (m_data != NULL)) {
			delete[] m_data;
		}
		m_fOwnsData = true;
		m_data = NULL;
	}

	///	<summary>
	///		Deallocate data from this DataChunk.
	///	</summary>
	void Deallocate() {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting to Deallocate an attached DataArray1D");
		}

		Detach();
	}

public:
	///	<summary>
	///		Get the size of the data.
	///	</summary>
	size_t GetTotalSize() const {
		return (m_sSize);
	}

	///	<summary>
	///		Get the number of elements in this DataArray1D.
	///	</summary>
	inline size_t GetRows() const {
		return m_sSize;
	}

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	void Assign(const DataArray1D<T> & da) {

		// Verify source array existence
		if (!da.IsAttached()) {
			if (!IsAttached()) {
				m_sSize = da.m_sSize;
				return;
			}

			_EXCEPTIONT("Attempting to assign unattached DataArray1D\n"
				"to attached DataArray1D (undefined behavior)");
		}

		// Allocate if necessary
		if (!IsAttached()) {
			Allocate(da.m_sSize);
		}
		if (IsAttached() && m_fOwnsData) {
			if (m_sSize != da.m_sSize) {
				Deallocate();
				Allocate(da.m_sSize);
			}
		}

		// Verify array consistency
		if (da.GetRows() != GetRows()) {
			_EXCEPTIONT("Size mismatch in assignment of DataArray1D");
		}

		// Copy data
		memcpy(m_data, da.m_data, GetByteSize());
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	DataArray1D<T> & operator= (const DataArray1D<T> & da) {
		Assign(da);
		return (*this);
	}

public:
	///	<summary>
	///		Zero the data content of this object.
	///	</summary>
	void Zero() {

		// Check that this DataArray1D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on uninitialized DataArray1D");
		}

		// Set content to zero
		memset(m_data, 0, m_sSize * sizeof(T));
	}

	///	<summary>
	///		Scale data by a given constant.
	///	</summary>
	void Scale(const T & x) {

		// Check that this DataArray1D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray1D");
		}

		// Scale data values
		for (size_t i = 0; i < m_sSize; i++) {
			m_data[i] *= x;
		}
	}

	///	<summary>
	///		Add a factor of the given DataArray1D to this DataArray1D.
	///	</summary>
	void AddProduct(
		const DataArray1D<T> & da,
		const T & x
	) {
		// Check that this DataArray1D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray1D");
		}
		if (!da.IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray1D");
		}
		if (da.GetRows() != GetRows()) {
			_EXCEPTIONT("Size mismatch in DataArray1D");
		}

		// Scale data values
		for (size_t i = 0; i < m_sSize; i++) {
			m_data[i] += x * da.m_data[i];
		}
	}

public:
	///	<summary>
	///		Implicit converstion to a pointer.
	///	</summary>
	inline operator T const*() const {
		return m_data;
	}
	///	<summary>
	///		Implicit converstion to a pointer.
	///	</summary>
	inline operator T*() {
		return m_data;
	}

	///	<summary>
	///		Parenthetical array accessor.
	///	</summary>
	inline T & operator()(size_t i) {
#if defined(DEBUG_ARRAYOUTOFBOUNDS)
		if (i > m_sSize) {
			_EXCEPTIONT("Array access out of bounds");
		}
#endif
		return m_data[i];
	}

	///	<summary>
	///		Parenthetical array accessor.
	///	</summary>
	inline const T & operator()(size_t i) const {
#if defined(DEBUG_ARRAYOUTOFBOUNDS)
		if (i > m_sSize) {
			_EXCEPTIONT("Array access out of bounds");
		}
#endif
		return m_data[i];
	}

private:
	///	<summary>
	///		A flag indicating this array owns its data.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		The number of rows in this DataArray1D.
	///	</summary>
	size_t m_sSize;

	///	<summary>
	///		A pointer to the data for this DataArray1D.
	///	</summary>
	T * m_data;

};

///////////////////////////////////////////////////////////////////////////////

#endif

