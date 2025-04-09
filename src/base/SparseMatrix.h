///////////////////////////////////////////////////////////////////////////////
///
///	\file    SparseMatrix.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include "DataArray1D.h"
#include <map>

///////////////////////////////////////////////////////////////////////////////

template <typename DataType>
class SparseMatrix {

public:
	///	<summary>
	///		Sparse matrix map.
	///	</summary>
	typedef typename std::pair<int, int> IndexType;
	typedef typename std::map<IndexType, DataType> SparseMap;
	typedef typename SparseMap::value_type SparseMapPair;
	typedef typename SparseMap::iterator SparseMapIterator;
	typedef typename SparseMap::const_iterator SparseMapConstIterator;
	typedef typename std::pair<SparseMapIterator, bool> SparseMapInsertResult;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	SparseMatrix() :
		m_nRows(0),
		m_nCols(0)
	{ }

public:
	///	<summary>
	///		Accessor.
	///	</summary>
	DataType & operator()(int iRow, int iCol) {
		SparseMapIterator iter = m_mapEntries.find(IndexType(iRow, iCol));
		if (iter == m_mapEntries.end()) {
			if (iRow >= m_nRows) {
				m_nRows = iRow + 1;
			}
			if (iCol >= m_nCols) {
				m_nCols = iCol + 1;
			}

			SparseMapInsertResult result =
				m_mapEntries.insert(
					SparseMapPair(
						IndexType(iRow, iCol), (DataType)(0)));

			return result.first->second;
		} else {
			return iter->second;
		}
	}

	///	<summary>
	///		Get the number of rows in the SparseMatrix.
	///	</summary>
	int GetRows() const {
		return m_nRows;
	}

	///	<summary>
	///		Get the number of columns in the SparseMatrix.
	///	</summary>
	int GetColumns() const {
		return m_nCols;
	}

	///	<summary>
	///		Get the entries of the SparseMatrix.
	///	</summary>
	void GetEntries(
		DataArray1D<int> & dataRows,
		DataArray1D<int> & dataCols,
		DataArray1D<DataType> & dataEntries
	) const {
		dataRows.Allocate(m_mapEntries.size());
		dataCols.Allocate(m_mapEntries.size());
		dataEntries.Allocate(m_mapEntries.size());

		int ix = 0;

		SparseMapConstIterator iter = m_mapEntries.begin();
		for (; iter != m_mapEntries.end(); iter++) {
			dataRows[ix] = iter->first.first;
			dataCols[ix] = iter->first.second;
			dataEntries[ix] = iter->second;

			ix++;
		}
	}

	///	<summary>
	///		Set the entries of the SparseMatrix in bulk.
	///	</summary>
	void SetEntries(
		const DataArray1D<int> & dataRows,
		const DataArray1D<int> & dataCols,
		const DataArray1D<DataType> & dataEntries
	) {
		if (dataRows.GetRows() != dataCols.GetRows()) {
			_EXCEPTIONT("Mismatch between size of dataRows and dataCols");
		}
		if (dataRows.GetRows() != dataEntries.GetRows()) {
			_EXCEPTIONT("Mismatch between size of dataRows and dataEntries");
		}

		m_nRows = 0;
		m_nCols = 0;

		m_mapEntries.clear();

		for (unsigned i = 0; i < dataRows.GetRows(); i++) {
			if (dataRows[i] >= m_nRows) {
				m_nRows = dataRows[i] + 1;
			}
			if (dataCols[i] >= m_nCols) {
				m_nCols = dataCols[i] + 1;
			}

			m_mapEntries.insert(
				SparseMapPair(
					IndexType(dataRows[i], dataCols[i]),
						dataEntries[i]));
		}
	}

	///	<summary>
	///		Clear the operator.
	///	</summary>
	void Clear() {
		m_nRows = 0;
		m_nCols = 0;
		m_mapEntries.clear();
	}

	///	<summary>
	///		Iterator to beginning of SparseMap.
	///	</summary>
	SparseMapIterator begin() {
		return m_mapEntries.begin();
	}

	///	<summary>
	///		Iterator to end of SparseMap.
	///	</summary>
	SparseMapIterator end() {
		return m_mapEntries.end();
	}

public:
	///	<summary>
	///		Apply the sparse matrix to a DataArray1D.
	///	</summary>
	void Apply(
		const DataArray1D<DataType> & dataVectorIn,
		DataArray1D<DataType> & dataVectorOut,
		bool fZeroOutputArray = true
	) const {
/*
		if (dataVectorIn.GetRows() != m_nCols) {
			_EXCEPTION1("dataVectorIn has incorrect row count (%i)", m_nCols);
		}
		if (dataVectorOut.GetRows() != m_nRows) {
			_EXCEPTION1("dataVectorOut has incorrect row count (%i)", m_nRows);
		}
*/
		if (fZeroOutputArray) {
			dataVectorOut.Zero();
		}

		SparseMapConstIterator iter = m_mapEntries.begin();
		for (; iter != m_mapEntries.end(); iter++) {
			if (dataVectorIn[iter->first.second]>1e19){
				// dataVectorOut[iter->first.first] += std::numeric_limits<double>::quiet_NaN();
				continue;
			} else {
				dataVectorOut[iter->first.first] +=
				iter->second * dataVectorIn[iter->first.second];
			}
		}
	}

protected:
	///	<summary>
	///		Number of rows in the sparse matrix.
	///	</summary>
	int m_nRows;

	///	<summary>
	///		Number of columns in the sparse matrix.
	///	</summary>
	int m_nCols;

	///	<summary>
	///		Entries of the sparse matrix.
	///	</summary>
	SparseMap m_mapEntries;
};

///////////////////////////////////////////////////////////////////////////////

#endif

