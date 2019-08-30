///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataArray4D.h
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

#ifndef _DATAARRAY4D_H_
#define _DATAARRAY4D_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"
#include "Subscript.h"

#include <cstdlib>
#include <cstring>

template <typename T>
class DataArray4D {

public:
	typedef T ValueType;

	///	<summary>
	///		Constructor.
	///	</summary>
	DataArray4D() :
		m_fOwnsData(true),
		m_data1D(NULL)
	{
		m_sSize[0] = 0;
		m_sSize[1] = 0;
		m_sSize[2] = 0;
		m_sSize[3] = 0;
	}

	///	<summary>
	///		Constructor allowing specification of size.
	///	</summary>
	DataArray4D(
		size_t sSize0,
		size_t sSize1,
		size_t sSize2,
		size_t sSize3,
		bool fAllocate = true
	) :
		m_fOwnsData(true),
		m_data1D(NULL)
	{
		m_sSize[0] = sSize0;
		m_sSize[1] = sSize1;
		m_sSize[2] = sSize2;
		m_sSize[3] = sSize3;

		if (fAllocate) {
			Allocate();
		}
	}

	///	<summary>
	///		Copy constructor.
	///	</summary>
	DataArray4D(const DataArray4D<T> & da) :
		m_fOwnsData(true),
		m_data1D(NULL)
	{
		if (da.IsAttached()) {
			m_sSize[0] = 0;
			m_sSize[1] = 0;
			m_sSize[2] = 0;
			m_sSize[3] = 0;

			Assign(da);

		} else {
			m_sSize[0] = da.m_sSize[0];
			m_sSize[1] = da.m_sSize[1];
			m_sSize[2] = da.m_sSize[2];
			m_sSize[3] = da.m_sSize[3];

			m_fOwnsData = true;

			m_data1D = NULL;
		}
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	virtual ~DataArray4D() {
		Detach();
	}

	///	<summary>
	///		Get the size of this DataChunk, in bytes.
	///	</summary>
	virtual size_t GetByteSize() const {
		size_t sSize =
			(m_sSize[0] * m_sSize[1] * m_sSize[2] * m_sSize[3] * sizeof(T));
		
		if (sSize % sizeof(size_t) == 0) {
			return sSize;
		} else {
			return (sSize / sizeof(size_t) + 1) * sizeof(size_t);
		}
	}

	///	<summary>
	///		Allocate data in this DataArray4D.
	///	</summary>
	void Allocate(
		size_t sSize0 = 0,
		size_t sSize1 = 0,
		size_t sSize2 = 0,
		size_t sSize3 = 0
	) {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting to Allocate() on attached DataArray4D");
		}

		if (sSize0 == 0) {
			sSize0 = m_sSize[0];
		}
		if (sSize1 == 0) {
			sSize1 = m_sSize[1];
		}
		if (sSize2 == 0) {
			sSize2 = m_sSize[2];
		}
		if (sSize3 == 0) {
			sSize3 = m_sSize[3];
		}
		if ((sSize0 == 0) ||
		    (sSize1 == 0) ||
		    (sSize2 == 0) ||
			(sSize3 == 0)
		) {
			_EXCEPTIONT("Attempting to Allocate() zero-size DataArray4D");
		}
		if ((m_data1D == NULL) ||
		    (m_sSize[0] != sSize0) ||
		    (m_sSize[1] != sSize1) ||
		    (m_sSize[2] != sSize2) ||
		    (m_sSize[3] != sSize3)
		) {
			Detach();

			m_sSize[0] = sSize0;
			m_sSize[1] = sSize1;
			m_sSize[2] = sSize2;
			m_sSize[3] = sSize3;

			m_data1D = reinterpret_cast<T *>(malloc(GetByteSize()));
		}

		Zero();

		m_fOwnsData = true;
	}

	///	<summary>
	///		Set the dimension sizes of this DataArray4D.
	///	</summary>
	inline void SetSize(
		size_t sSize0,
		size_t sSize1,
		size_t sSize2,
		size_t sSize3
	) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting SetSize() on attached DataArray4D");
		}

		m_sSize[0] = sSize0;
		m_sSize[1] = sSize1;
		m_sSize[2] = sSize2;
		m_sSize[3] = sSize3;
	}

public:
	///	<summary>
	///		Determine if this DataChunk is attached to a data array.
	///	</summary>
	virtual bool IsAttached() const {
		return (m_data1D != NULL);
	}

	///	<summary>
	///		Attach this DataChunk to an array of pre-allocated data.
	///	</summary>
	virtual void AttachToData(void * ptr) {
		if (IsAttached()) {
			_EXCEPTIONT("Attempting AttachToData() on attached DataArray4D");
		}

		m_data1D = reinterpret_cast<T *>(ptr);
		m_fOwnsData = false;
	}

	///	<summary>
	///		Detach data from this DataChunk.
	///	</summary>
	virtual void Detach() {
		if ((m_fOwnsData) && (m_data1D != NULL)) {
			delete[] m_data1D;
		}
		m_fOwnsData = true;
		m_data1D = NULL;
	}

	///	<summary>
	///		Deallocate data from this DataChunk.
	///	</summary>
	void Deallocate() {
		if (!m_fOwnsData) {
			_EXCEPTIONT("Attempting Deallocate() on attached DataArray4D");
		}

		Detach();
	}

public:
	///	<summary>
	///		Get the size of the data.
	///	</summary>
	size_t GetTotalSize() const {
		return (m_sSize[0] * m_sSize[1] * m_sSize[2] * m_sSize[3]);
	}


	///	<summary>
	///		Get the size of the specified dimension.
	///	</summary>
	inline size_t GetSize(int dim) const {
		return m_sSize[dim];
	}

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	void Assign(const DataArray4D<T> & da) {

		// Verify source array existence
		if (!da.IsAttached()) {
			if (!IsAttached()) {
				m_sSize[0] = da.m_sSize[0];
				m_sSize[1] = da.m_sSize[1];
				m_sSize[2] = da.m_sSize[2];
				m_sSize[3] = da.m_sSize[3];

				return;
			}

			_EXCEPTIONT("Attempting to assign unattached DataArray4D\n"
				"to attached DataArray4D (undefined behavior)");
		}

		// Allocate if necessary
		if (!IsAttached()) {
			Allocate(
				da.m_sSize[0],
				da.m_sSize[1],
				da.m_sSize[2],
				da.m_sSize[3]);
		}
		if (IsAttached() && m_fOwnsData) {
			if ((m_sSize[0] != da.m_sSize[0]) ||
			    (m_sSize[1] != da.m_sSize[1]) ||
			    (m_sSize[2] != da.m_sSize[2]) ||
			    (m_sSize[3] != da.m_sSize[3])
			) {
				Deallocate();
				Allocate(
					da.m_sSize[0],
					da.m_sSize[1],
					da.m_sSize[2],
					da.m_sSize[3]);
			}
		}

		// Check initialization status
		if (da.GetSize(0) != GetSize(0)) {
			_EXCEPTIONT("Dimension 0 mismatch in DataArray4D");
		}
		if (da.GetSize(1) != GetSize(1)) {
			_EXCEPTIONT("Dimension 1 mismatch in DataArray4D");
		}
		if (da.GetSize(2) != GetSize(2)) {
			_EXCEPTIONT("Dimension 2 mismatch in DataArray4D");
		}
		if (da.GetSize(3) != GetSize(3)) {
			_EXCEPTIONT("Dimension 3 mismatch in DataArray4D");
		}

		// Copy data
		memcpy(m_data1D, da.m_data1D, GetByteSize());
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	DataArray4D<T> & operator= (const DataArray4D<T> & da) {
		Assign(da);
		return (*this);
	}

public:
	///	<summary>
	///		Zero the data content of this object.
	///	</summary>
	void Zero() {

		// Check that this DataArray4D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray4D");
		}

		// Set content to zero
		memset(m_data1D, 0, GetByteSize());
	}

	///	<summary>
	///		Scale data by a given constant.
	///	</summary>
	void Scale(const T & x) {

		// Check that this DataArray4D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray4D");
		}

		// Scale data values
		size_t sTotalSize = GetTotalSize();

		for (size_t i = 0; i < sTotalSize; i++) {
			m_data1D[i] *= x;
		}
	}

	///	<summary>
	///		Add a factor of the given DataArray4D to this DataArray4D.
	///	</summary>
	void AddProduct(
		const DataArray4D<T> & da,
		const T & x
	) {
		// Check that this DataArray4D is attached to a data object
		if (!IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray4D");
		}
		if (!da.IsAttached()) {
			_EXCEPTIONT("Attempted operation on unattached DataArray4D");
		}
		if (da.GetSize(0) != GetSize(0)) {
			_EXCEPTIONT("Dimension 0 mismatch in DataArray4D");
		}
		if (da.GetSize(1) != GetSize(1)) {
			_EXCEPTIONT("Dimension 1 mismatch in DataArray4D");
		}
		if (da.GetSize(2) != GetSize(2)) {
			_EXCEPTIONT("Dimension 2 mismatch in DataArray4D");
		}
		if (da.GetSize(3) != GetSize(3)) {
			_EXCEPTIONT("Dimension 3 mismatch in DataArray4D");
		}

		// Scale data values
		size_t sTotalSize = GetTotalSize();

		for (size_t i = 0; i < sTotalSize; i++) {
			m_data1D[i] += x * da.m_data1D[i];
		}
	}

public:
	///	<summary>
	///		Subscript DSEL operator.
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline Subscript<DataArray4D<T> const, 3, 4>
	operator[](std::ptrdiff_t idx) const
#else
	inline Subscript<DataArray4D<T> const, 3, 4>
	operator[](std::ptrdiff_t idx) const noexcept
#endif
	{
		Subscript<DataArray4D<T> const, 4, 4> s(*this);
		return s[idx];
	}

	///	<summary>
	///		Subscript DSEL operator.
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline Subscript<DataArray4D<T>, 3, 4>
	operator[](std::ptrdiff_t idx)
#else
	inline Subscript<DataArray4D<T>, 3, 4>
	operator[](std::ptrdiff_t idx) noexcept
#endif
	{
		Subscript<DataArray4D<T>, 4, 4> s(*this);
		return s[idx];
	}

	///	<summary>
	///		Parenthetical array accessor.
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline T const&
	operator()(index_array<std::ptrdiff_t, 4> indices) const
#else
	inline T const&
	operator()(index_array<std::ptrdiff_t, 4> indices) const noexcept
#endif
	{
		return (*this)(indices[0], indices[1], indices[2], indices[3]);
	}

	///	<summary>
	///		Parenthetical array accessor.
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline T&
	operator()(index_array<std::ptrdiff_t, 4> indices)
#else
	inline T&
	operator()(index_array<std::ptrdiff_t, 4> indices) noexcept
#endif
	{
		return (*this)(indices[0], indices[1], indices[2], indices[3]);
	}

	///	<summary>
	///		Parenthetical array accessor.
	///	</summary>
	inline const T & operator()(size_t i, size_t j, size_t k, size_t l) const {
#if defined(DEBUG_ARRAYOUTOFBOUNDS)
		if ((i > m_sSize[0]) || (j > m_sSize[1]) || (k > m_sSize[2]) || (l > m_sSize[3])) {
			_EXCEPTIONT("Array access out of bounds");
		}
#endif
		return (*(m_data1D + i * m_sSize[1] * m_sSize[2] * m_sSize[3]
				+ j * m_sSize[2] * m_sSize[3] + k * m_sSize[3] + l));
	}
	///	<summary>
	///		Parenthetical array accessor.
	///	</summary>
	inline T & operator()(size_t i, size_t j, size_t k, size_t l) {
#if defined(DEBUG_ARRAYOUTOFBOUNDS)
		if ((i > m_sSize[0]) || (j > m_sSize[1]) || (k > m_sSize[2]) || (l > m_sSize[3])) {
			_EXCEPTIONT("Array access out of bounds");
		}
#endif
		return (*(m_data1D + i * m_sSize[1] * m_sSize[2] * m_sSize[3]
				+ j * m_sSize[2] * m_sSize[3] + k * m_sSize[3] + l));
	}

	///	<summary>
	///		Parenthetical array accessor (unit-stride slicer).
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline T const*
	operator()(index_array<std::ptrdiff_t, 3> indices) const
#else
	inline T const*
	operator()(index_array<std::ptrdiff_t, 3> indices) const noexcept
#endif
	{
		return (*this)(indices[0], indices[1], indices[2]);
	}
	///	<summary>
	///		Parenthetical array accessor (unit-stride slicer).
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline T*
	operator()(index_array<std::ptrdiff_t, 3> indices)
#else
	inline T*
	operator()(index_array<std::ptrdiff_t, 3> indices) noexcept
#endif
	{
		return (*this)(indices[0], indices[1], indices[2]);
	}

	///	<summary>
	///		Parenthetical array accessor (unit-stride slicer).
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline T const* operator()(size_t i, size_t j, size_t k) const {
#else
	inline T const* operator()(size_t i, size_t j, size_t k) const noexcept {
#endif
#if defined(DEBUG_ARRAYOUTOFBOUNDS)
		if ((i > m_sSize[0]) || (j > m_sSize[1]) || (k > m_sSize[2])) {
			_EXCEPTIONT("Array access out of bounds");
		}
#endif
		return m_data1D + i * m_sSize[1] * m_sSize[2] * m_sSize[3]
				+ j * m_sSize[2] * m_sSize[3] + k * m_sSize[3];
	}

	///	<summary>
	///		Parenthetical array accessor (unit-stride slicer).
	///	</summary>
#if defined(__INTEL_COMPILER)
	inline T* operator()(size_t i, size_t j, size_t k) {
#else
	inline T* operator()(size_t i, size_t j, size_t k) noexcept {
#endif
#if defined(DEBUG_ARRAYOUTOFBOUNDS)
		if ((i > m_sSize[0]) || (j > m_sSize[1]) || (k > m_sSize[2])) {
			_EXCEPTIONT("Array access out of bounds");
		}
#endif
		return m_data1D + i * m_sSize[1] * m_sSize[2] * m_sSize[3]
				+ j * m_sSize[2] * m_sSize[3] + k * m_sSize[3];
	}

private:
	///	<summary>
	///		A flag indicating this array owns its data.
	///	</summary>
	bool m_fOwnsData;

	///	<summary>
	///		The size of each dimension of this DataArray4D.
	///	</summary>
	size_t m_sSize[4];

	///	<summary>
	///		A pointer to the data for this DataArray4D.
	///	</summary>
	T * m_data1D;
};

///////////////////////////////////////////////////////////////////////////////

#endif

