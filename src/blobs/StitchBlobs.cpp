///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchBlobs.cpp
///	\author  Paul Ullrich Hongyu Chen
///	\version July 2, 2023
///
///	<remarks>
///		Copyright 2023 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

#include "Constants.h"
#include "CoordTransforms.h"
#include "BlobUtilities.h"
#include <cmath>
#include <memory>

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "FilenameList.h"
#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"
#include "SimpleGrid.h"
#include "Variable.h"
#include "kdtree.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <queue>
#include <exception>

///////////////////////////////////////////////////////////////////////////////

struct Tag {

	///	<summary>
	///		Identifier associated with this tag.
	///	</summary>
	int id;

	///	<summary>
	///		Time index associated with this tag.
	///	</summary>
	int time;

	///	<summary>
	///		Global id associated with this tag (minimum value 1).
	///	</summary>
	int global_id;

	///	<summary>
	///		Default constructor.
	///	</summary>
	Tag() :
		id(0),
		time(0),
		global_id(0)
	{ }

	///	<summary>
	///		Value-based constructor.
	///	</summary>
	Tag(int a_id, int a_time) :
		id(a_id),
		time(a_time),
		global_id(0)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	Tag(const Tag & tag) :
		id(tag.id),
		time(tag.time),
		global_id(tag.global_id)
	{ }

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const Tag & tag) const {
		if (time < tag.time) {
			return true;
		}
		if (time > tag.time) {
			return false;
		}
		if (id < tag.id) {
			return true;
		}
		return false;
	}
};



#if defined(TEMPEST_MPIOMP)
///	<summary>
///		A class used for MPI collective communication for vecAllBlobTags
///		(gather and scatter) between processors)
///		This class will serialize the vecAllBlobTags into an 1-D array
///		and deserialize the received 1-D vecAllBlobTags into the original
///		2-D vectors after the communication. 
///		This class can also generate and return the setAllTags based on the local vecAllBlobTags.
///	</summary>
class TagCollectiveOP {
	private:
		///	<summary>
		///		Tag for MPI communication for sending/receiving Tag
		///	</summary>
		int gather_tag = 107;

		///	<summary>
		///		Tag for MPI communication for sending/receiving Tag index
		///	</summary>
		int gather_tag_index = 108;

		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		The MPI Datatype for Tag
		///	</summary>
		MPI_Datatype MPI_Tag_type;

		///	<summary>
		///		The Flag that marks whether the vecAllTags are in serialized state. 
		///		(default with 0; 1 after calling Serialize(); False after calling Deserailize()).	
		///	</summary>
		int serializedFlag;

		///	<summary>
		///		Serialize the std::vector< std::vector<Tag>> vecAllBlobTags into a 1-D array 
		///		and generate the index array serialVecAllBlobTags_index
		///	</summary>
		void Serialize() {
			serialVecAllBlobTags.clear();
			serialVecAllBlobTags_index.clear();
			int curIndx = 0;//Point to the next empty slot for inserting a new vector<Tag>
			serialVecAllBlobTags_index.push_back(curIndx);
			for (int i = 0; i < _vecAllBlobTags.size(); i++) {
				for (int j = 0; j < _vecAllBlobTags[i].size(); j++) {
					serialVecAllBlobTags.push_back(_vecAllBlobTags[i][j]);
					curIndx++;
				}
				serialVecAllBlobTags_index.push_back(curIndx);
			}
			serializedFlag = 1;
		}

		///	<summary>
		///		Deserialize the the local std::vector<Tag> serialVecAllBlobTags array 
		///		to generate the deserialSetAllTags and desirialVecAllBlobTags
		///		the this.serialVecAllBlobTags and this.serialVecAllBlobTags_index will be cleared
		///		after the deserialization
		///	</summary>
		void Deserialize() {
			for (int i = 0; i < serialVecAllBlobTags_index.size() - 1; i++) {
				int startIndx = serialVecAllBlobTags_index[i];
				int endIndx = std::min(serialVecAllBlobTags_index[i+1],int(serialVecAllBlobTags.size()));
				std::vector<Tag> curVecBlobTags;
				for (int i = startIndx; i < endIndx; i++) {
					curVecBlobTags.push_back(serialVecAllBlobTags[i]);
					deserialSetAllTags.insert(serialVecAllBlobTags[i]);
				}
				desirialVecAllBlobTags.push_back(curVecBlobTags);
			}
			serialVecAllBlobTags.clear();
			serialVecAllBlobTags_index.clear();
			serializedFlag = 0;
		}
		


	protected:
		///	<summary>
		///		The unexchanged original local 2-D vecAllBlobTags that needs to be serialize before communication.		
		///	</summary>
		std::vector< std::vector<Tag>> _vecAllBlobTags;

		///	<summary>
		///		The serialized 1D vecAllBlobTags that needs to be sent.	
		///	</summary>
		std::vector<Tag> serialVecAllBlobTags;

		///	<summary>
		///		The serialized 1D vecAllBlobTags index array that needs to be sent.	
		///	</summary>
		std::vector<int> serialVecAllBlobTags_index;

		///	<summary>
		///		Output setAllTags based on the local desirialVecAllBlobTags. (Only used for processor 0)	
		///	</summary>	
		std::set<Tag> deserialSetAllTags;

		///	<summary>
		///		The deserial vecAllBlobTags after communication (Only used for processor 0)	
		///	</summary>
		std::vector< std::vector<Tag>> desirialVecAllBlobTags;

		///	<summary>
		///		Vector that records the size information of vecAllBlobTags on each processors	
		///	</summary>
		std::vector<int> vecScatterCounts;

		///	<summary>
		///		Vector that records the index information of the serial vecAllBlobTags on each processors	
		///	</summary>
		std::vector<int> vecScatterCounts_index;

	public:

		///	<summary>
		///		Constructor that will read in std::vector< std::vector<Tag>> vecAllBlobTags and MPI communicator
		///		It will also create the derived MPI_Datatype for Tag and commit it.
		///	</summary>	
		TagCollectiveOP(
			MPI_Comm communicator, 
			const std::vector< std::vector<Tag>> & vecAllBlobTags
		) {
			this->_vecAllBlobTags = vecAllBlobTags;
			this->m_comm = communicator;
			this->serializedFlag = 0;
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			// Create an MPI datatype for the Tag:
			struct Tag sampleTag;
			int tagFieldsCount = 3;	
			MPI_Datatype Tag_typesig[3] = {MPI_INT,MPI_INT,MPI_INT};
			int Tag_block_lengths[3] = {1,1,1};
			MPI_Aint Tag_displacements[3];
	
			MPI_Aint base_address;
			MPI_Get_address(&sampleTag, &base_address);
			MPI_Get_address(&sampleTag.id, &Tag_displacements[0]);
			MPI_Get_address(&sampleTag.time, &Tag_displacements[1]);
			MPI_Get_address(&sampleTag.global_id, &Tag_displacements[2]);

			Tag_displacements[0] = MPI_Aint_diff(Tag_displacements[0], base_address);
			Tag_displacements[1] = MPI_Aint_diff(Tag_displacements[1], base_address);
			Tag_displacements[2] = MPI_Aint_diff(Tag_displacements[2], base_address);

			MPI_Type_create_struct(tagFieldsCount, Tag_block_lengths, Tag_displacements, Tag_typesig, &MPI_Tag_type);

			int result = MPI_Type_commit(&MPI_Tag_type);
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Type_commit(&MPI_Tag_type) failed (code %i)", result);
			}
		}

		///	<summary>
		///		Destructor.
		///	</summary>
		~TagCollectiveOP(){
			MPI_Type_free(&MPI_Tag_type);			
		}

		///	<summary>
		///		The MPI gather process that will gather each local vecAllBlobTags to the processor 0.
		///	</summary>
		void Gather() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			this->Serialize();
			int d = std::ceil(std::log2(size)); // here we use floor.
			
			for (int j = 0; j < d; j++) {

				// First send the serialVecAllBlobTags
				if ((rank & (int)std::round(std::pow(2,j))) != 0) {
					
					// Send to (world_rank ^ pow(2,j)
					int destRank = rank ^ (int)round(std::pow(2,j));
					if (destRank > size - 1) {
						continue;
					}

					MPI_Send (serialVecAllBlobTags.data(), serialVecAllBlobTags.size(), MPI_Tag_type, destRank, gather_tag, m_comm);
					MPI_Send (serialVecAllBlobTags_index.data(), serialVecAllBlobTags_index.size(), MPI_INT, destRank, gather_tag_index, m_comm);
					// Simply need to break the algorithm here (juest return, not Finalize())
					return;


				} else {
					// Receive from (world_rank ^ pow(2,j))
					MPI_Status status;
					int recvCount;
					int sourceRank = rank ^ (int)std::round(std::pow(2,j));

					if (sourceRank > size - 1) {
						continue;
					}

					MPI_Probe(sourceRank, gather_tag,  m_comm, &status);

					MPI_Get_count( &status, MPI_Tag_type, &recvCount);					 
					std::vector<Tag> recvTags;
					recvTags.resize(recvCount);
					MPI_Recv(recvTags.data(), recvCount, MPI_Tag_type, sourceRank, gather_tag, m_comm, &status);

					// Pack the receive Tag into the local serialVecAllBlobTags.
					for (auto recvTag : recvTags) {
						serialVecAllBlobTags.push_back(recvTag);
					}

					MPI_Status status_index;
					int recvCount_index;
					MPI_Probe(sourceRank, gather_tag_index,  m_comm, &status_index);

					MPI_Get_count( &status_index, MPI_INT, &recvCount_index);					 
					std::vector<int> recvTagsIndx;
					recvTagsIndx.resize(recvCount_index);
					MPI_Recv(recvTagsIndx.data(), recvCount_index, MPI_INT, sourceRank, gather_tag_index, m_comm, &status_index);

					
					// Update the received index and then Pack the receive Tag index into the
					// local serialVecAllBlobTags_index.
					// Example:
					// Initial:
					// P0 serialVecAllBlobTags: 0, 3, 5, 7;   P1 serialVecAllBlobTags: 0, 3, 5, 7
					// After Gather:
					// P0 serialVecAllBlobTags: 0, 3, 5, 7, 9, 11, 13
					int serialVecAllBlobTags_index_size = serialVecAllBlobTags_index.size();
					int curLocalTagSize = serialVecAllBlobTags_index[serialVecAllBlobTags_index_size - 1];
					for (int i = 1; i < recvTagsIndx.size(); i++) {
						// Update the index
						int index = recvTagsIndx[i];
						index += curLocalTagSize;
						serialVecAllBlobTags_index.push_back(index);
					}
				}
			}
		}

		///	<summary>
		///		Return the gathered setAllTags (only called by the processor 0)
		///	</summary>
		std::set<Tag> GetGatheredSetAllTags() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			if (rank == 0) {
				if (serializedFlag == 1) {
					this->Deserialize();
				}

				return deserialSetAllTags;
			} else {
				_EXCEPTIONT("Only processor 0 should call GetGatheredSetAllTags().");
			}
		}

		///	<summary>
		///		Return the gathered/scattered vecAllTags (only called by the processor 0)
		///		GatheredFlag = 1 indicates that returning the vecAllBlobTags after gathering to the processor 0
		///		(in this case, vecAllBlobTags is the global vecAllBlobTags and only processor 0 can call the function)
		///		GatheredFalg = 0 indicates that returning the vecAllBlobTags after scattering to each processot.
		///		(in this case, vecAllBlobTags is the local vecAllBlobTags and all valid processor can call the function)
		///	</summary>
		std::vector< std::vector<Tag> > GetUnserialVecAllTags(int GatheredFlag) {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			if ((GatheredFlag == 1 && rank == 0) || GatheredFlag== 0) {
				if (serializedFlag == 1) {
					this->Deserialize();
				}
				return desirialVecAllBlobTags;
			} else {
				_EXCEPTIONT("Only processor 0 should call GetUnserialVecAllTags().");
			}
		}

		///	<summary>
		///		Gather the size info of the original unexchanged vecAllBlobTags to processor 0
		///		On processor 0, the _vecAllBlobTags will be the gathered global vecAllBlobTags;
		///		On other processors, the _vecAllBlobTags will be updated to the input vecAllBlobTags
		///	</summary>
		void GatherTagCounts(const std::vector< std::vector<Tag> > & vecAllBlobTags) {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			int curSize = 0;		
			for ( auto vecBlobTags : vecAllBlobTags) {
				curSize += vecBlobTags.size();			
			}
			int localSize = vecAllBlobTags.size();
			vecScatterCounts.resize(size);
			vecScatterCounts_index.resize(size);

			if (rank == 0) {				
				MPI_Gather(&curSize, 1, MPI_INT, vecScatterCounts.data(), 1, MPI_INT, 0, m_comm);
				MPI_Gather(&localSize, 1, MPI_INT, vecScatterCounts_index.data(), 1, MPI_INT, 0, m_comm);

			} else {
				this->_vecAllBlobTags = vecAllBlobTags;
				MPI_Gather(&curSize, 1, MPI_INT, NULL, 0, MPI_INT, 0,m_comm);
				MPI_Gather(&localSize, 1, MPI_INT, NULL, 0, MPI_INT, 0,m_comm);

			}
		}


		///	<summary>
		///		Scatter the vecAllBlobs to each processor
		///		The displacement is calculated based on the vecGlobalTimes
		///		On processor 0, the _vecAllBlobTags will be the reduced global vecAllBlobTags;
		///		On other processors, the _vecAllBlobTags will be the original unreduced _vecAllBlobTags
		///	</summary>
		void Scatter(){
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			this->Serialize();

			if (rank == 0) {
				// For Tags
				int arrayScatterCounts[size];
				int arrayScatterDisplacements[size];
				_ASSERT(vecScatterCounts.size() > 1);
				arrayScatterCounts[0] = vecScatterCounts[0];
				arrayScatterDisplacements[0] = 0;		
				for (int i = 1; i < vecScatterCounts.size(); i++) {
					arrayScatterCounts[i] = vecScatterCounts[i];				
					arrayScatterDisplacements[i] = vecScatterCounts[i - 1] + arrayScatterDisplacements[i - 1];
				}

				// For Tags Index
				int arrayScatterCounts_index[size];
				int arrayScatterDisplacements_index[size];
				_ASSERT(vecScatterCounts_index.size() > 1);
				arrayScatterCounts_index[0] = vecScatterCounts_index[0];
				arrayScatterDisplacements_index[0] = 0;		
				for (int i = 1; i < vecScatterCounts.size(); i++) {
					arrayScatterCounts_index[i] = vecScatterCounts_index[i];				
					arrayScatterDisplacements_index[i] = vecScatterCounts_index[i - 1] + arrayScatterDisplacements_index[i - 1];
				}
				
				//-------------------Scatter---------------------------
				// For Tags
				auto scatterBuffer = this->serialVecAllBlobTags;
				this->serialVecAllBlobTags.clear();
				this->serialVecAllBlobTags.resize(arrayScatterCounts[rank]);
				MPI_Scatterv(scatterBuffer.data(), arrayScatterCounts, arrayScatterDisplacements, MPI_Tag_type, 
							serialVecAllBlobTags.data(), arrayScatterCounts[rank], MPI_Tag_type, 0, m_comm);

				// For index
				auto scatterBuffer_index = this->serialVecAllBlobTags_index;
				this->serialVecAllBlobTags_index.clear();
				this->serialVecAllBlobTags_index.resize(arrayScatterCounts_index[rank]);
				MPI_Scatterv(scatterBuffer_index.data(), arrayScatterCounts_index, arrayScatterDisplacements_index, 
							MPI_INT, serialVecAllBlobTags_index.data(), arrayScatterCounts_index[rank], MPI_INT, 0, m_comm);
	
			} else {
				int localTagSize = serialVecAllBlobTags.size();
				this->serialVecAllBlobTags.clear();
				this->serialVecAllBlobTags.resize(localTagSize);
				MPI_Scatterv(NULL, NULL, NULL, MPI_Tag_type, serialVecAllBlobTags.data(), localTagSize, MPI_Tag_type, 0, m_comm);

				int localTagSize_index = serialVecAllBlobTags_index.size();
				this->serialVecAllBlobTags_index.clear();
				this->serialVecAllBlobTags_index.resize(localTagSize_index);
				MPI_Scatterv(NULL, NULL, NULL, MPI_INT, serialVecAllBlobTags_index.data(), localTagSize_index, MPI_INT, 0, m_comm);
			}

			// Now modify the received serialVecAllBlobTags_index for deserailization call
			int prevCount = serialVecAllBlobTags_index[0];
			serialVecAllBlobTags_index[0] = 0;				
			for (int i = 1; i < serialVecAllBlobTags_index.size() - 1; i++ ) {
				serialVecAllBlobTags_index[i] = serialVecAllBlobTags_index[i] - prevCount;
			}
			if (rank == 0) {
				serialVecAllBlobTags_index.push_back(serialVecAllBlobTags.size());

			} else {
				serialVecAllBlobTags_index[serialVecAllBlobTags_index.size() - 1] = serialVecAllBlobTags.size();
			}
		}
};

///	<summary>
///		Enumerator of exchange directions.
///	</summary>
typedef enum {
	DIR_LEFT = 0, 
	DIR_RIGHT = 1
} CommDirection;

///	<summary>
///		A Class used for exchanging vecAllBlobsTag between processors
///		The exchange process will first update the exchanged Tags.time to the actual global time
///		And then start the exchange process.
///	</summary>
class TagExchangeOP {
	private:

		///	<summary>
		///		Tag for MPI communication for blob Tag
		///	</summary>
		int tag = 100;

		///	<summary>
		///		The MPI Datatype for Tag
		///	</summary>
		MPI_Datatype MPI_Tag_type;
		
		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;

		///	<summary>
		///		The tool function that uses the prefix sum algorithmn to assign global time for each Tag
		///	</summary>		
		void UpdateTime(){
			int err, rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			this->gloablTimeIndx.resize(size);
			int localEndTime = _vecAllBlobTags.size();
			int preFixTime;
			
			MPI_Scan(&localEndTime, &preFixTime, 1, MPI_INT, MPI_SUM, m_comm);
			this->gloablTimeIndx[rank].resize(2);
			this->gloablTimeIndx[rank][0] = preFixTime - localEndTime;
			this->gloablTimeIndx[rank][1] = preFixTime;

			//Update Tag.time
			int globalTime = this->gloablTimeIndx[rank][0];
			for (int i = 0; i < _vecAllBlobTags.size(); i++) {
				_ASSERT(globalTime < gloablTimeIndx[rank][1]);
				
				for (int j = 0; j < _vecAllBlobTags[i].size(); j++) {
					_vecAllBlobTags[i][j].time = globalTime;
				}
				globalTime++;

			}
		}

	protected:
		///	<summary>
		///		The initial vecAllBlobTags that need to be exchanged
		///	</summary>
		std::vector< std::vector<Tag>> _vecAllBlobTags;

		///	<summary>
		///		The vecAllBlobTags after the exchange. it's column number is _vecAllBlobTags' column number plus 2 (left and right) 
		///		except for p0 and pn-1 (for these two processors, the column number is _vecAllBlobTags' column number plus 1).
		///	</summary>
		std::vector< std::vector<Tag>> exchangedvecAllBlobTags;

		///	<summary>
		///		The buffer for vecAllBlobTags that will be sent
		///		sendTags[0] is the left vector and sendTags[1] is the right vector
		///	</summary>
		std::vector<std::vector<Tag>> sendTags;

		///	<summary>
		///		The buffer for vecAllBlobTags that will be received
		///		recvTags[0] is the left vector and recvTags[1] is the right vector
		///	</summary>
		std::vector<std::vector<Tag>> recvTags;

		///	<summary>
		///		The array (size is nMPISize) records the start(inclusive) and end (exclusive) global time index in each processor
		///		start: gloablTimeIndx[p_i][0]
		///		end: gloablTimeIndx[p_i][1]
		///	</summary>	
		std::vector<std::vector<int>> gloablTimeIndx;	

	public:

		///	<summary>
		///		Construct the Operator with vecAllBlobTags
		///		It will contruct the this->m_comm and this->_vecAllBlobTags based on the input communicator and vecAllBlobTags
		///		And also construct the derived MPI_Datatype for Tag and commit it.
		///	</summary>
		TagExchangeOP(MPI_Comm communicator, 
					  const std::vector< std::vector<Tag> > & vecAllBlobTags){
			this->_vecAllBlobTags = vecAllBlobTags;
			this->m_comm = communicator;
			//Initialize the size for the sendTags:
			sendTags.resize(2);
			sendTags[0].resize(_vecAllBlobTags[0].size());
			sendTags[1].resize(_vecAllBlobTags[_vecAllBlobTags.size()-1].size());

			//Initialize the size for the recvTags:
			recvTags.resize(2);

			//Create an MPI datatype for the Tag:
			struct Tag sampleTag;
			int tagFieldsCount = 3;	
			MPI_Datatype Tag_typesig[3] = {MPI_INT,MPI_INT,MPI_INT};
			int Tag_block_lengths[3] = {1,1,1};
			MPI_Aint Tag_displacements[3];
	
			MPI_Aint base_address;
			MPI_Get_address(&sampleTag, &base_address);
			MPI_Get_address(&sampleTag.id, &Tag_displacements[0]);
			MPI_Get_address(&sampleTag.time, &Tag_displacements[1]);
			MPI_Get_address(&sampleTag.global_id, &Tag_displacements[2]);
			Tag_displacements[0] = MPI_Aint_diff(Tag_displacements[0], base_address);
			Tag_displacements[1] = MPI_Aint_diff(Tag_displacements[1], base_address);
			Tag_displacements[2] = MPI_Aint_diff(Tag_displacements[2], base_address);
			MPI_Type_create_struct(tagFieldsCount, Tag_block_lengths, Tag_displacements, Tag_typesig, &MPI_Tag_type);

			int result = MPI_Type_commit(&MPI_Tag_type);
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Type_commit(&MPI_Tag_type) failed (code %i)", result);
			}
		}


		///	<summary>
		///		Destructor
		///	</summary>
		~TagExchangeOP(){
			MPI_Type_free(&MPI_Tag_type);
			MPIrequests.clear();
			MPIstatuses.clear();
		}

		///	<summary>
		///		Return the original unexchanged vecAllBlobTags
		///	</summary>
		std::vector< std::vector<Tag> > GetOriginalVecAllBlobTags(){
			return _vecAllBlobTags;
		}

		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the TagExchangeOP should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {
			
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			
			//First update all processors' Tag.time to the actual global time
			this->UpdateTime();

			//----------------------Send sendTags data first----------------------
			// Pack data into the send buffer
			sendTags[0] = _vecAllBlobTags[0];
			sendTags[1] = _vecAllBlobTags[_vecAllBlobTags.size()-1];

			// Send data
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						//Rank 0 Do Nothing
						destRank = -1;
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						//Rank n-1 Do Nothing
						destRank = -1;
						continue;
					} else {
						destRank = rank + 1;
					}
				}
				if (destRank > size - 1) {
						continue;
				}

				//----------------------Send sendBlobs----------------------

				// Only the odd number processors will send out the data
				if (rank % 2 != 0) {

					MPI_Request request;
					int result = MPI_Isend(sendTags[dir].data(), sendTags[dir].size(), MPI_Tag_type,
							destRank, tag, m_comm, &request);
					if (result != MPI_SUCCESS) {
						_EXCEPTION1("The MPI routine MPI_Isend failed (code %i)", result);
					}
				}
			}

			//----------------------Then Receive data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank = -1;				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					// Rank 0 will not receive from the left		 
					if (rank == 0) {			
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					// rank n-1 will not receive from the right
					if (rank == size - 1) {
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}
				if (sourceRank > size - 1) {
						continue;
				}



				//----------------------Receive----------------------	
				// Only the prime number processors will receive data
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount;

					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, tag, m_comm, &flag, &status );
					}
					MPI_Get_count( &status, MPI_Tag_type, &recvCount );
					recvTags[dir].resize(recvCount);

					int result =
						MPI_Irecv(recvTags[dir].data(), recvTags[dir].size(), MPI_Tag_type,
							sourceRank, tag, m_comm, &request);
					if (result != MPI_SUCCESS) {
						_EXCEPTION1("The MPI routine MPI_Isend failed (code %i)", result);
					}
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());
				}
			}
		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the TagExchangeOP without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {
			// Wait for all Irecv to complete
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}

			MPIrequests.clear();
			MPIstatuses.clear();

			// Pack the data into the vecAllBlobTags
			// Only need to pack data for the prime number processors 
			if (rank % 2 == 0) {
				if (rank == 0) {
					exchangedvecAllBlobTags.resize(_vecAllBlobTags.size() + 1);					
					for (int i = 0; i < exchangedvecAllBlobTags.size() - 1; i++) {
						exchangedvecAllBlobTags[i] = _vecAllBlobTags[i];
					}
					exchangedvecAllBlobTags[exchangedvecAllBlobTags.size() - 1] = recvTags[1];

				} else if (rank == size - 1) {
					exchangedvecAllBlobTags.resize(_vecAllBlobTags.size() + 1);
					exchangedvecAllBlobTags[0] = recvTags[0];
					for (int i = 1; i < exchangedvecAllBlobTags.size(); i++) {
						exchangedvecAllBlobTags[i] = _vecAllBlobTags[i -1 ];
					}

				} else {	

					exchangedvecAllBlobTags.resize(_vecAllBlobTags.size() + 2);
					exchangedvecAllBlobTags[0] = recvTags[0];
					for (int i = 1; i < exchangedvecAllBlobTags.size() - 1; i++) {
						exchangedvecAllBlobTags[i] = _vecAllBlobTags[i-1];
					}
					exchangedvecAllBlobTags[exchangedvecAllBlobTags.size() - 1] = recvTags[1];

				}

			} else {
				exchangedvecAllBlobTags = _vecAllBlobTags;
			}
		}

		///	<summary>
		///		Return the exchanged vecAllBlobTags
		///	</summary>
		std::vector< std::vector<Tag>> GetExchangedVecAllBlobTags(){
			return this->exchangedvecAllBlobTags;
		}
};


///	<summary>
///		Class for exchanging vecAllBlobBoxesDeg among processors
///	</summary>
class BlobBoxesDegExchangeOP {
	private:
		///	<summary>
		///		Tag for MPI communication for LatlonBox<double>
		///	</summary>
		int tag = 103;

		// ///	<summary>
		// ///		The MPI Datatype for LatlonBox<double>
		// ///	</summary>
		// MPI_Datatype MPI_LatonBox_double_type;

		// ///	<summary>
		// ///		The MPI Datatype for double[2]
		// ///	</summary>
		// MPI_Datatype MPI_doubleArray;

		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;

	protected:
		///	<summary>
		///		The initial vecAllBlobBoxesDeg that will be sent
		///	</summary>
		std::vector< std::vector< LatLonBox<double> > > _vecAllBlobBoxesDeg;

		///	<summary>
		///		The vecAllBlobBoxesDeg that is after the exchange. it's column number is _vecAllBlobBoxesDegs' column number plus 2 (left and right) 
		///		except for p0 and pn-1 (for these two processors, the column number is _vecAllBBlobBoxesDegs' column number plus 1).
		///	</summary>
		std::vector< std::vector< LatLonBox<double> > > exchangedvecAllBlobBoxesDeg;

		///	<summary>
		///		The buffer for vecAllBlobBoxesDeg that will be sent
		///		sendBlobBoxesDeg[0] is the left vector and sendBlobBoxesDeg[1] is the right vector
		///	</summary>
		std::vector< std::vector< LatLonBox<double> > > sendBlobBoxesDeg;

		///	<summary>
		///		The buffer for vecAllBlobBoxesDeg that will be received
		///		sendBlobBoxesDeg[0] is the left vector and sendBlobBoxesDeg[1] is the right vector
		///	</summary>
		std::vector< std::vector< LatLonBox<double> > > recvBlobBoxesDeg;

	public:

		///	<summary>
		///		Construct the Operator with vecAllBlobBoxesDeg
		///		It will contruct the this->m_comm and this->_vecAllBlobBoxesDeg based on the input communicator and vecAllBlobBoxesDeg
		///	</summary>
		BlobBoxesDegExchangeOP(MPI_Comm communicator, 
							   const std::vector< std::vector< LatLonBox<double> > > & vecAllBlobBoxesDeg){
			this->_vecAllBlobBoxesDeg = vecAllBlobBoxesDeg;
			this->m_comm = communicator;

		//########################### Notes for  derived MPI datatype of the LatLonBox<double>(Hongyu Chen) ############################################################################ 
		//1.  "Because vector<bool> holds bits instead of bools, it can't return a bool& from its indexing operator or iterator dereference" (src: https://isocpp.org/blog/2012/11/on-vectorbool)
		//    Therefore, creating the userdefined datatype for LatLonBox<double> and then use vector.data() to send/recv like TagExhangeOP is not working here
		//2.  Now we're using the MPI_BYTE to manually calculate the send/recv buffer size in byte here, which is working currently
		//3.  If the program breaks again here, please consider going to the BlobUtilities.h and modify the constructer at line 232 and line 236 according to the description there.

			// //Create an MPI datatype for the LatLonBox<double>:
			// //First create the datatype for double[2]
			// MPI_Type_contiguous	(2,	MPI_DOUBLE,&MPI_doubleArray);	
			// MPI_Type_commit	(&MPI_doubleArray);
			// //Then use this doubleArray to construct LatlonBox<double>
			// LatLonBox<double> sampleBox;
			// int LatlonBoxFieldCount = 5;
			// MPI_Datatype LatlonBox_typesig[5] = {MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_DOUBLE, MPI_doubleArray, MPI_doubleArray};
			// int LatlonBox_block_lengths[5] = {1,1,1,1,1};
			// MPI_Aint LatlongBox_displacements[5];
			// MPI_Aint base_address;
			// MPI_Get_address(&sampleBox.is_null, &LatlongBox_displacements[0]);
			// MPI_Get_address(&sampleBox.lon_periodic, &LatlongBox_displacements[1]);
			// MPI_Get_address(&sampleBox.lon_width, &LatlongBox_displacements[2]);
			// MPI_Get_address(&sampleBox.lon, &LatlongBox_displacements[3]);
			// MPI_Get_address(&sampleBox.lat, &LatlongBox_displacements[4]);
			// LatlongBox_displacements[0] = MPI_Aint_diff(LatlongBox_displacements[0], base_address);
			// LatlongBox_displacements[1] = MPI_Aint_diff(LatlongBox_displacements[1], base_address);
			// LatlongBox_displacements[2] = MPI_Aint_diff(LatlongBox_displacements[2], base_address);
			// LatlongBox_displacements[3] = MPI_Aint_diff(LatlongBox_displacements[3], base_address);
			// LatlongBox_displacements[4] = MPI_Aint_diff(LatlongBox_displacements[4], base_address);
			// MPI_Type_create_struct(LatlonBoxFieldCount, LatlonBox_block_lengths, LatlongBox_displacements, LatlonBox_typesig, &MPI_LatonBox_double_type);
			// MPI_Type_commit(&MPI_LatonBox_double_type);

			//########################### End Notes for  derived MPI datatype of the LatLonBox<double>(Hongyu Chen) ############################################################################ 

			sendBlobBoxesDeg.resize(2);
			recvBlobBoxesDeg.resize(2);
		}


		///	<summary>
		///		Destructor.
		///	</summary>
		~BlobBoxesDegExchangeOP(){
			// MPI_Type_free(&MPI_LatonBox_double_type);
			// MPI_Type_free(&MPI_doubleArray);
			MPIrequests.clear();
			MPIstatuses.clear();
		}

		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the BlobBoxesDegExchangeOP should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			//----------------------Send data first----------------------
			// Pack data into the send buffer
			sendBlobBoxesDeg[0] = _vecAllBlobBoxesDeg[0];
			sendBlobBoxesDeg[1] = _vecAllBlobBoxesDeg[_vecAllBlobBoxesDeg.size()-1];

			// Send data
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						// Rank 0 Do Nothing
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						//Rank n-1 Do Nothing
						continue;
					} else {
						destRank = rank + 1;
					}
				}

				if (destRank > size - 1) {
					continue;
				}

				//----------------------Send sendBlobs----------------------
				// Only the odd number processors will send out data
				if (rank % 2 != 0) {
					MPI_Request request;
					MPI_Isend(sendBlobBoxesDeg[dir].data(), sendBlobBoxesDeg[dir].size() * sizeof(LatLonBox<double>), MPI_BYTE,
					destRank,tag , m_comm, &request);
				}
			}

			//----------------------Then Receive data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank;				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					if (rank == 0) {// Rank 0 will not receive from the left
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					if (rank == size - 1) {// rank n-1 will not receive from the right
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}

				if (sourceRank > size - 1) {
					continue;
				}

				// Only the even number processors will receive data 
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount;
					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, tag, m_comm, &flag, &status );
					}
					MPI_Get_count( &status, MPI_BYTE, &recvCount );	
					_ASSERT(recvCount % sizeof(LatLonBox<double>) == 0);			
					recvBlobBoxesDeg[dir].resize(recvCount / sizeof(LatLonBox<double>));			
					MPI_Irecv(recvBlobBoxesDeg[dir].data(), recvCount,MPI_BYTE,
							sourceRank, tag, m_comm, &request);
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());
				}
			}
		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the BlobBoxesDegExchangeOP data without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {
			// Wait for all Irecv to complete

			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}

			MPIrequests.clear();
			MPIstatuses.clear();
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			// Pack the data into the vecAllBlobBoxesDeg
			// The std::move() here is to avoid the possible errors in LatLonBox<double> default copy constructor
			// Only the prime number processors need to pack the data
			if (rank % 2 == 0) {
				if (rank == 0) {
					for (int i = 0; i < _vecAllBlobBoxesDeg.size();i++) {
						exchangedvecAllBlobBoxesDeg.emplace_back(std::move(_vecAllBlobBoxesDeg[i]));
					}
					exchangedvecAllBlobBoxesDeg.emplace_back(std::move(recvBlobBoxesDeg[1]));

				} else if (rank == size - 1) {
					exchangedvecAllBlobBoxesDeg.emplace_back(std::move(recvBlobBoxesDeg[0]));
					for (int i = 0; i < _vecAllBlobBoxesDeg.size(); i++) {
						exchangedvecAllBlobBoxesDeg.emplace_back(std::move(_vecAllBlobBoxesDeg[i]));
					}

				} else {
					exchangedvecAllBlobBoxesDeg.emplace_back(std::move(recvBlobBoxesDeg[0]));
					for (int i = 0; i < _vecAllBlobBoxesDeg.size(); i++) {
						exchangedvecAllBlobBoxesDeg.emplace_back(std::move(_vecAllBlobBoxesDeg[i]));
					}
					exchangedvecAllBlobBoxesDeg.emplace_back(std::move(recvBlobBoxesDeg[1]));
				}

			} else {
				exchangedvecAllBlobBoxesDeg = _vecAllBlobBoxesDeg;
			}


		}

		///	<summary>
		///		Return the exchanged vecAllBlobTags
		///	</summary>
		std::vector< std::vector< LatLonBox<double> > > GetExchangedVecAllBlobBoxesDeg(){
			return this->exchangedvecAllBlobBoxesDeg;
		}
};

#endif 

///////////////////////////////////////////////////////////////////////////////

// Set of indicator locations stored as grid indices
typedef std::set<int> IndicatorSet;
typedef IndicatorSet::iterator IndicatorSetIterator;
typedef IndicatorSet::const_iterator IndicatorSetConstIterator;


#if defined(TEMPEST_MPIOMP) 

///	<summary>
///		A Class used for exchanging vecAllBlobs between processors
///	</summary>
class BlobsExchangeOp {
	private:

		///	<summary>
		///		Tag for MPI communication for sending/receiving blobs
		///	</summary>
		int blob_tag = 101;

		///	<summary>
		///		Tag for MPI communication for sending/receiving serialized index
		///	</summary>
		int indx_tag = 102;


		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;


		/// <summary>
		/// 	Serialize the boundary core column from _vecAllBlobsWithGhosts and generate the sendBlobsIndx array.
		/// 	For DIR_LEFT, we serialize the first core column; for DIR_RIGHT, we serialize the last core column.
		/// 	Note: This function is only used on odd-numbered processors, where _vecAllBlobsWithGhosts holds only core data.
		/// </summary>
		void Serialize(){
			// Get MPI rank and size.
			int rank, size;
			MPI_Comm_rank(m_comm, &rank);
			MPI_Comm_size(m_comm, &size);

			// For odd-numbered processors, no ghost columns exist.
			int ghostLeft, ghostRight;
			if (rank % 2 != 0) {
				ghostLeft = 0;
				ghostRight = 0;
			} else {
				// For even-numbered processors (which have been resized to include ghost cells):
				// Include this implementation here just in case
				ghostLeft = (rank == 0) ? 0 : 1;
				ghostRight = (rank == size - 1) ? 0 : 1;
			}
			
			// _vecAllBlobsWithGhosts is allocated to have:
			// totalColumns = ghostLeft + coreColumns + ghostRight.
			int totalColumns = static_cast<int>(_vecAllBlobsWithGhosts.size());
			// Determine the number of core columns.
			int coreColumns = totalColumns - ghostLeft - ghostRight;
			
			// Calculate indices of the core boundary columns.
			int leftCoreIndex = ghostLeft;                // first core column.
			int rightCoreIndex = ghostLeft + coreColumns - 1; // last core column.
			
			// Clear previous serialization data.
			sendBlobs.clear();
			sendBlobsIndx.clear();
			sendBlobs.resize(2);
			sendBlobsIndx.resize(2);
			
			// Iterate over both directions.
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int curIndx = 0; // Next empty slot for serialized data.
				
				// For odd processors, the container holds core data only.
				std::vector<IndicatorSet> sendVecBlobs = (dir == DIR_LEFT) ?
					_vecAllBlobsWithGhosts[leftCoreIndex] :
					_vecAllBlobsWithGhosts[rightCoreIndex];
				
				// Record the starting index.
				sendBlobsIndx[dir].push_back(curIndx);
				
				// Serialize each IndicatorSet from the selected column.
				for (size_t i = 0; i < sendVecBlobs.size(); i++) {
					IndicatorSet curSet = sendVecBlobs[i];
					for (auto it = curSet.begin(); it != curSet.end(); ++it) {
						sendBlobs[dir].push_back(*it);
						curIndx++;
					}
					// Record the boundary index after each set.
					sendBlobsIndx[dir].push_back(curIndx);
				}
			}
		}


		///	<summary>
		///		Deserialize the received vector<int>recvBlobs into vector<IndicatorSet> and clear the recvBlobsIndx array
		///	</summary>
		void Deserialize(){
			recvBlobsUnserial.clear();
			recvBlobsUnserial.resize(2);
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);


			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				//Rank 0 will not receive from its left, hence no deserialization.
				if (rank == 0 && dir == DIR_LEFT) {
					continue;
				}

				//the last processor will now receive from its right, hence no deserialization
				if ((rank == size - 1) && dir == DIR_RIGHT) {
					continue;
				}
				for (int i = 0; i < recvBlobsIndx[dir].size()-1; i++){
					IndicatorSet curSet;
					int startIndx = recvBlobsIndx[dir][i];
					int endIndx = std::min(recvBlobsIndx[dir][i+1],int(recvBlobs[dir].size()));

					//Deserialize each set
					for (int i = startIndx; i < endIndx; i++){
						curSet.insert(recvBlobs[dir][i]);
					}

					//push each set into the vector
					recvBlobsUnserial[dir].push_back(curSet);
				}
			}
			recvBlobsIndx.clear();
			sendBlobsIndx.clear();
		}
	protected:

		/// <summary>
		/// 	The local vecAllBlobs including ghost cells.
		/// 	This container stores both the original (core) data and the extra columns for ghost cells.
		/// 	The core region corresponds to the original unexchanged data, while the ghost cells (left and right)
		/// 	are reserved for the exchanged boundary data.
		/// 	For interior processors, the total number of columns is coreColumns + 2 (one ghost column on each side).
		/// 	For boundary processors, it may be coreColumns + 1.
		/// </summary>
		std::vector<std::vector<IndicatorSet>> _vecAllBlobsWithGhosts;

		///	<summary>
		///		The buffer for vecBlobs that is serialized and will be sent
		///		sendBlobs[0] is the left vector and sendBlobs[1] is the right vector
		///	</summary>
		std::vector<std::vector<int>> sendBlobs;


		///	<summary>
		///		The Array recording the starting index for each set that is serialized
		///		sendBlobsIndx[0] is for the left vector and sendBlobsIndx[1] is for the right vector
		///	</summary>
		std::vector<std::vector<int>> sendBlobsIndx;

		///	<summary>
		///		The buffer for the received serailized blobs array and need to be unserialized
		///		recvBlobs[0] is for the left vector and recvBlobs[1] is for the right vector
		///	</summary>
		std::vector<std::vector<int>> recvBlobs;

		///	<summary>
		///		The Buffer  recording the starting index for each set that is in the recvBlobs
		///		recvBlobsIndx[0] is for the left vector and recvBlobsIndx[1] is for the right vector
		///	</summary>
		std::vector<std::vector<int>> recvBlobsIndx;

		///	<summary>
		///		The array recording the unserialize recvBlobs
		///		recvBlobsUnserial[0] is for the left vector and recvBlobsUnserial[1] is for the right vector	
		///	</summary>
		std::vector<std::vector<IndicatorSet>> recvBlobsUnserial;


	public:
		/// <summary>
		/// 	Construct the Operator with BlobsExchangeOp.
		/// 	Simply save the original core data (vecAllBlobs) in _vecAllBlobsWithGhosts.
		/// 	No extra ghost cell allocation is done at this stage.
		/// </summary>
		BlobsExchangeOp(MPI_Comm communicator, 
			const std::vector< std::vector<IndicatorSet> > & vecAllBlobs) {
			this->m_comm = communicator;
			// Save the core data as-is.
			_vecAllBlobsWithGhosts = vecAllBlobs;
		}


		///	<summary>
		///		Destructor for BlobsExchangeOp
		///	</summary>
		~BlobsExchangeOp(){
			MPIrequests.clear();
			MPIstatuses.clear();
			
		}


		///	<summary>
		/// 	Start the exchange process (non-blocking).
		/// 	This function initiates the halo exchange by sending out the serialized boundary (core) data
		/// 	from _vecAllBlobsWithGhosts and posting receives for ghost data from neighbors.
		/// 	The received data (still serialized) will later be integrated directly into the ghost cell regions
		/// 	of _vecAllBlobsWithGhosts in EndExchange().
		///	</summary>
		void StartExchange() {			
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			recvBlobs.resize(2);
			recvBlobsIndx.resize(2);

			// Only even-numbered (prime) processors RECIEVE ghost data and need to expand the container.
			if (rank % 2 == 0) {
				// Determine ghost offsets.
				int ghostLeft  = (rank == 0) ? 0 : 1;
				int ghostRight = (rank == size - 1) ? 0 : 1;
				// Core data columns count (the current size of _vecAllBlobsWithGhosts).
				int coreColumns = static_cast<int>(_vecAllBlobsWithGhosts.size());
				// Total columns after adding ghost columns.
				int newTotalColumns = ghostLeft + coreColumns + ghostRight;

				// Allocate a temporary container for the expanded data.
				std::vector<std::vector<IndicatorSet>> ghostContainer;
				ghostContainer.resize(newTotalColumns);

				// Copy the core data into the new container at the proper offset.
				for (int i = 0; i < coreColumns; i++) {
					ghostContainer[i + ghostLeft] = _vecAllBlobsWithGhosts[i];
				}
				
				// Replace the in-place container with the expanded one.
				_vecAllBlobsWithGhosts = std::move(ghostContainer);
			}
			// For odd-number processors, no resizing is done; _vecAllBlobsWithGhosts remains as the original core data.


			//----------------------Send sendBlobs/sendBlobsIndx data----------------------
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						//Rank 0 Do Nothing
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						// Rank n-1 Do Nothing
						continue;
					} else {
						destRank = rank + 1;
					}
				}
				if (destRank > size - 1) {
					continue;

				}

				// only the odd number processors will send out data
				if (rank % 2 != 0) {
					// First Serialize the Sending Buffer.
					this->Serialize();

					
					//----------------------Send sendBlobs----------------------
					MPI_Request request;
					MPI_Isend(sendBlobs[dir].data(), sendBlobs[dir].size(), MPI_INT,
					destRank, blob_tag, m_comm, &request);

					//----------------------Send sendBlobsIndx----------------------
					MPI_Request indx_Request;
					MPI_Isend(sendBlobsIndx[dir].data(), sendBlobsIndx[dir].size(), MPI_INT,
					destRank, indx_tag, m_comm, &indx_Request);

				}



			}

			//----------------------Receive sendBlobs/recvBlobsIndx Data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank;

				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					if (rank == 0) {// Rank 0 will not receive from the left
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					if (rank == size - 1) {// rank n-1 will not receive from the right
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}
				if (sourceRank > size - 1) {
					continue;

				}

				//----------------------Receive the serialized Blobs----------------------

				//only the prime number processor will receive the blobs
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount;

					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, blob_tag, m_comm, &flag, &status );
					}
					MPI_Get_count( &status, MPI_INT, &recvCount );
					recvBlobs[dir].resize(recvCount);
					MPI_Irecv(recvBlobs[dir].data(), recvBlobs[dir].size(), MPI_INT,
							sourceRank, blob_tag, m_comm, &request);
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());

					//----------------------Receive the index info for the Blobs----------------------
					MPI_Status indxStatus;
					MPI_Request indxRequest;
					int indxRecvCount;

					// Use a non-blocking probe to know the incoming data size
					int indxFlag = 0;
					while(!indxFlag)
					{
						MPI_Iprobe( sourceRank, indx_tag, m_comm, &indxFlag, &indxStatus);
					}
					MPI_Get_count( &indxStatus, MPI_INT, &indxRecvCount);
					recvBlobsIndx[dir].resize(indxRecvCount);
					MPI_Irecv(recvBlobsIndx[dir].data(), recvBlobsIndx[dir].size(), MPI_INT,
							sourceRank, indx_tag, m_comm, &indxRequest);
					MPIrequests.emplace_back(std::move(indxRequest));
					MPIstatuses.push_back(MPI_Status());

				}

			}


		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the BlobsExchangeOp data without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {

			// Wait for all Irecv to complete
			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}

			MPIrequests.clear();
			MPIstatuses.clear();
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			// Compute ghost offsets based on processor rank.
			int ghostLeft  = (rank == 0) ? 0 : 1;
			int ghostRight = (rank == size - 1) ? 0 : 1;
			
			// Total number of columns in our in-place container.
			int totalColumns = static_cast<int>(_vecAllBlobsWithGhosts.size());
			// Core columns are the original columns (without ghost cells).
			int coreColumns = totalColumns - ghostLeft - ghostRight;

			// For odd processors, _vecAllBlobsWithGhosts remains at core size.
				

			// Only the even number processors need deserialize
			if (rank % 2 == 0) {
				this->Deserialize();

				if (rank == 0) {
					// Rank 0 has no left ghost cell; update right ghost column (if exists).
					if (ghostRight == 1) {
						// Right ghost column is at index ghostLeft + coreColumns.
						_vecAllBlobsWithGhosts[ghostLeft + coreColumns] = recvBlobsUnserial[1];
					}

				} else if (rank == size - 1) {
					// The last processor has no right ghost cell; update left ghost column.
					if (ghostLeft == 1) {
						_vecAllBlobsWithGhosts[0] = recvBlobsUnserial[0];
					}
				} else {
					// Interior processors update both ghost columns.
					_vecAllBlobsWithGhosts[0] = recvBlobsUnserial[0];
					_vecAllBlobsWithGhosts[ghostLeft + coreColumns] = recvBlobsUnserial[1];
				}
			} 
			// Odd-numbered processors did not expand their container; no modifications are needed.
			// They already hold exactly the original core data.
		}

		/// <summary>
		/// 	Return the exchanged vecAllBlobs.
		/// 	For even-numbered processors, this is the expanded container with ghost cells;
		/// 	for odd-numbered processors, it is the same as the original core data.
		/// </summary>
		std::vector<std::vector<IndicatorSet>> GetExchangedVecAllBlobs() {
			// Simply return the container as it stands after EndExchange.
			// Even processors have ghost cells; odd processors have core data.
			return _vecAllBlobsWithGhosts;
		}

		/// <summary>
		/// 	Return the original (core) vecAllBlobs without ghost cells.
		/// 	For even-numbered processors, we extract the core region from the expanded container;
		/// 	for odd-numbered processors, the container already holds exactly the core data.
		/// </summary>
		std::vector<std::vector<IndicatorSet>> GetOriginalVecAllBlobs() {
			int rank, size;
			MPI_Comm_rank(m_comm, &rank);
			MPI_Comm_size(m_comm, &size);
			
			int ghostLeft  = (rank == 0) ? 0 : 1;
			int ghostRight = (rank == size - 1) ? 0 : 1;
			
			// For even processors, _vecAllBlobsWithGhosts is expanded:
			// total columns = coreColumns + ghostLeft + ghostRight.
			// For odd processors, _vecAllBlobsWithGhosts is exactly the core data.
			std::vector<std::vector<IndicatorSet>> coreData;
			
			if (rank % 2 == 0) {
				// Even processor: extract the core region.
				int totalColumns = static_cast<int>(_vecAllBlobsWithGhosts.size());
				int coreColumns = totalColumns - ghostLeft - ghostRight;
				coreData.resize(coreColumns);
				for (int i = 0; i < coreColumns; i++) {
					coreData[i] = _vecAllBlobsWithGhosts[i + ghostLeft];
				}
			}
			else {
				// Odd processors: the container is already core data.
				coreData = _vecAllBlobsWithGhosts;
			}
			
			return coreData;
		}

};

///	<summary>
///		A Class used for exchanging vecGlobalTimes between processors
///		This class will update the local global time accordingly after the exchange process.
///	</summary>	
class GlobalTimesExchangeOp {
	private:
		///	<summary>
		///		Tag for MPI communication for vecGlobalTimes
		///	</summary>
		int tag = 104;

		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;

	protected:
	///	<summary>
	///		The initial vecGlobalTimes before exchange
	///	</summary>
	std::vector< std::vector<Time> > _vecGlobalTimes;

	///	<summary>
	///		The vecGlobalTimes that is after the exchange. exchangedVecGlobalTimes[pi]'s size will be  is _vecGlobalTimes[pi]'s size plus 2 (left and right) 
	///		except for p0 and pn-1 (for these two processors, exchangedVecGlobalTimes[pi]'s size is _vecGlobalTimes[pi]'s size plus 1).
	///	</summary>
	std::vector< std::vector<Time> > exchangedVecGlobalTimes;

	///	<summary>
	///		The buffer for vecGlobalTimes that will be sent
	///		sendTimes[0] is the left vector and sendTimes[1] is the right vector
	///	</summary>
	std::vector<Time> sendTimes;

	///	<summary>
	///		The buffer for vecGlobalTimes that will be received
	///		recvTimes[0] is the left vector and sendTimes[1] is the right vector
	///	</summary>
	std::vector<Time> recvTimes;	

	///	<summary>
	///		Lowerbound of the file number that this processor is reading (starts with index 0 and inclusive)
	///	</summary>	
	int fileLowerBound;

	///	<summary>
	///		Upperbound of the file number that this processor is reading (starts with index 0 and exclusive)
	///	</summary>
	int fileUpperBound;		

	public:

		///	<summary>
		///		It will contruct the this->m_comm this->_vecGlobalTimes this->fileLowerBound and this->fileUpperBound based on the input.
		///	</summary>
		GlobalTimesExchangeOp(
			MPI_Comm communicator, 
			const std::vector< std::vector<Time> > & vecGlobalTimes, 
			const int & processorResponsibalForFile_LB, 
			const int & processorResponsibalForFile_UB
		) {
			this->_vecGlobalTimes = vecGlobalTimes;
			this->m_comm = communicator;
			this->fileLowerBound = processorResponsibalForFile_LB;
			this->fileUpperBound = processorResponsibalForFile_UB;
			sendTimes.resize(2);
			recvTimes.resize(2);
		}

		///	<summary>
		///		Destructor
		///	</summary>
		~GlobalTimesExchangeOp(){
			MPIrequests.clear();
			MPIstatuses.clear();
		}

		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the GlobalTimesExchangeOp should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);


			//----------------------Send Data----------------------
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						// Rank 0 Do Nothing
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						// Rank n-1 Do Nothing
						continue;
					} else {
						destRank = rank + 1;
					}
				}
				if (destRank > size - 1) {
					continue;
				}


				// Only the odd number processors will send out data
				if (rank % 2 != 0) {
					// Pack data into the send buffer
					sendTimes[0] = _vecGlobalTimes[fileLowerBound][0];
					sendTimes[1] = _vecGlobalTimes[fileUpperBound-1][_vecGlobalTimes[fileUpperBound-1].size()-1];
					MPI_Request request;
					MPI_Isend(&sendTimes[dir], sizeof(Time), MPI_BYTE,
					destRank,tag , m_comm, &request);
				}
			}

			//----------------------Then Receive data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank;				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					if (rank == 0) {// Rank 0 will not receive from the left
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					if (rank == size - 1) {// rank n-1 will not receive from the right
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}

				if (sourceRank > size - 1) {
					continue;
				}

				// Only the peven number processors will receive data 
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount;
					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, tag, m_comm, &flag, &status );
					}
					MPI_Get_count( &status, MPI_BYTE, &recvCount );						
					MPI_Irecv(&recvTimes[dir], recvCount, MPI_BYTE,
							sourceRank, tag, m_comm, &request);
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());
				}
			}
		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the GlobalTimesExchangeOp data without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {
			// Wait for all Irecv to complete
			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}

			MPIrequests.clear();
			MPIstatuses.clear();
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			exchangedVecGlobalTimes.resize(_vecGlobalTimes.size());

			// Pack the data into the vecGlobalTimes 
			// Only the even number processors need to pack the data
			if (rank % 2 == 0) {
				exchangedVecGlobalTimes = _vecGlobalTimes;
								
				if (rank == 0) {
					// Put the Time received from the left at the [fileUpperBound][end] position
										
					exchangedVecGlobalTimes[fileUpperBound - 1].push_back(recvTimes[1]);;

				} else if (rank == size - 1) {
					// Put the Time received from the left at the [0][0] position
					exchangedVecGlobalTimes[fileLowerBound].clear();
					exchangedVecGlobalTimes[fileLowerBound].push_back(recvTimes[0]);
					for (int i = 0; i < _vecGlobalTimes[fileLowerBound].size(); i++) {
						exchangedVecGlobalTimes[fileLowerBound].push_back(_vecGlobalTimes[fileLowerBound][i]);
					}

				} else {
					// Put the Time received from the left at the [0][0] position
					exchangedVecGlobalTimes[fileLowerBound].clear();
					
					exchangedVecGlobalTimes[fileLowerBound].push_back(recvTimes[0]);
					for (int i = 0; i < _vecGlobalTimes[fileLowerBound].size(); i++) {
						exchangedVecGlobalTimes[fileLowerBound].push_back(_vecGlobalTimes[fileLowerBound][i]);
					}
					// Put the Time received from the right at the [fileUpperBound][end] position
					exchangedVecGlobalTimes[fileUpperBound - 1].push_back(recvTimes[1]);;
				}
			} else {
				// for odd number processors, nothing is modified.
				exchangedVecGlobalTimes = _vecGlobalTimes;
			}

		}


		///	<summary>
		///		Return the exchanged vecGlobalTimes
		///	</summary>
		std::vector< std::vector<Time> >  GetExchangedVecGlobalTimes(){
			return this->exchangedVecGlobalTimes;
		}

		///	<summary>
		///		Return the unexchanged vecGlobalTimes
		///	</summary>
		std::vector< std::vector<Time> >  GetUnExchangedVecGlobalTimes(){
			return this->_vecGlobalTimes;
		}
};
#endif 

// Array of equivalent tags
typedef std::multimap<Tag, Tag> MapGraph;
typedef MapGraph::const_iterator MapGraphConstIterator;
typedef MapGraph::iterator MapGraphIterator;


#if defined(TEMPEST_MPIOMP) 


///	<summary>
///		A Class for MPI communication for sending/receiving Tag pairs (MapGraph)
///	</summary>
class MapGraphGatherOp {
	private:
		///	<summary>
		///		Tag for MPI communication for sending/receiving Tag pairs
		///	</summary>
		int reduce_tag = 106;

		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
		std::vector<MPI_Status> MPIstatuses;


		///	<summary>
		///		Serialize the MapGraph and generate the local std::pair<Tag, Tag> array
		///	</summary>
		void Serialize(){
			localPairs.clear();
			for (auto it = _multimapTagGraph.begin(); it != _multimapTagGraph.end(); ++it) {
				localPairs.push_back(std::pair<Tag, Tag>(it->first, it->second));
			}

		}

		///	<summary>
		///		Deserialize the local std::pair<Tag, Tag> array and generate back
		///		the MapGraph (Only Processor 0 will call it)
		///	</summary>
		void Deserialize() {
			for (std::pair<Tag, Tag> tagPair : localPairs) {
				outputTagGraph.insert(tagPair);
			}

		}

	protected:
		///	<summary>
		///		Original multimapTagGraph.	
		///	</summary>
		MapGraph _multimapTagGraph;


		///	<summary>
		///		Serialized local multimapTagGraph.	
		///	</summary>
		std::vector<std::pair<Tag, Tag>> localPairs;

		///	<summary>
		///		Output multimapTagGraph. (Only used for processor 0)	
		///	</summary>	
		MapGraph outputTagGraph;
	
	public:
		///	<summary>
		///		Construct the Operator with multimapTagGraph
		///		It will contruct the this->m_comm and this->_multimapTagGraph based on
		///		the input communicator and multimapTagGraph
		///		And also construct the derived MPI_Datatype for Tag and commit it.
		///	</summary>
		MapGraphGatherOp(MPI_Comm communicator, 
						 const MapGraph & multimapTagGraph) {

			this->_multimapTagGraph = multimapTagGraph;
			this->m_comm = communicator;
		}

		///	<summary>
		///		Destructor
		///	</summary>
		~MapGraphGatherOp(){
			MPIrequests.clear();
			MPIstatuses.clear();
		}



		///	<summary>
		///		The MPI gather that will gather each local MapGraph to the processor 0.
		///	</summary>
		void Gather() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			this->Serialize();
			int d = std::ceil(std::log2(size));//here we use ceil.
			for (int j = 0; j < d; j++) {
				if ((rank & (int)std::round(std::pow(2,j))) != 0) {
					
					//send to (world_rank ^ pow(2,j)
					int destRank = rank ^ (int)round(std::pow(2,j));
					if (destRank > size - 1) {
						continue;
					}

					MPI_Send (localPairs.data(), localPairs.size() * sizeof(std::pair<Tag, Tag>), MPI_BYTE, destRank, reduce_tag, m_comm);
					//  We just simply need to break the algorithm here (juest return, not Finalize())
					return;
				} else {
					//receive from (world_rank ^ pow(2,j))
					MPI_Status status;
					int byteCount;
					int sourceRank = rank ^ (int)std::round(std::pow(2,j));

					if (sourceRank > size - 1) {
						continue;
					}

					MPI_Probe(sourceRank, reduce_tag,  m_comm, &status);

					MPI_Get_count( &status, MPI_BYTE, &byteCount);					 
					std::vector<std::pair<Tag, Tag>> recvTagPairs;
					recvTagPairs.resize(byteCount / sizeof(std::pair<Tag, Tag>));
					MPI_Recv(recvTagPairs.data(), byteCount, MPI_BYTE, sourceRank, reduce_tag, m_comm, &status);

					//Pack the receive Pairs into the localPairs.
					for (auto recvPair : recvTagPairs) {
						localPairs.push_back(recvPair);
					}



				}

			}

		}

		///	<summary>
		///	Return the gathered multiGraph (only called by the processor 0)
		///	</summary>
		MapGraph GetGatheredTagGraph() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			if (rank == 0) {
				this->Deserialize();
				return outputTagGraph;
			} else {
				_EXCEPTIONT("Only processor 0 should call GetGatheredTagGraph().");
			}
			

		}

};

#endif 
///////////////////////////////////////////////////////////////////////////////

class BlobThresholdOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum ThresholdQuantity {
		MinArea,
		MaxArea,
		MinArealFraction,
		MaxArealFraction,
		//EastwardOrientation,
		//WestwardOrientation
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Quantity,
			ReadMode_Value,
			//ReadMode_MinCount,
			ReadMode_Invalid
		} eReadMode = ReadMode_Quantity;

		// Loop through string
		int iLast = 0;
		for (int i = 0; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in threshold quantity
				if (eReadMode == ReadMode_Quantity) {
					if (strSubStr == "minarea") {
						m_eQuantity = MinArea;
						eReadMode = ReadMode_Value;

					} else if (strSubStr == "maxarea") {
						m_eQuantity = MaxArea;
						eReadMode = ReadMode_Value;

					} else if (strSubStr == "minarealfraction") {
						m_eQuantity = MinArealFraction;
						eReadMode = ReadMode_Value;

					} else if (strSubStr == "maxarealfraction") {
						m_eQuantity = MaxArealFraction;
						eReadMode = ReadMode_Value;
/*
					} else if (strSubStr == "eastwardorientation") {
						m_eQuantity = EastwardOrientation;
						eReadMode = ReadMode_Invalid;

					} else if (strSubStr == "westwardorientation") {
						m_eQuantity = WestwardOrientation;
						eReadMode = ReadMode_Invalid;
*/
					} else {
						_EXCEPTION1("Threshold invalid quantity \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					m_dValue = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;
/*
				// Read in minimum count
				} else if (eReadMode == ReadMode_MinCount) {
					if (strSubStr == "all") {
						m_nMinimumCount = (-1);
					} else {
						m_nMinimumCount = atoi(strSubStr.c_str());
					}

					if (m_nMinimumCount < -1) {
						_EXCEPTION1("Invalid minimum count \"%i\"",
							m_nMinimumCount);
					}

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;
*/
				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("Too many entries in threshold string \"%s\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("Insufficient entries in threshold string \"%s\"",
					strOp.c_str());
		}

		// Output announcement
		std::string strDescription;

		char szValue[128];
		snprintf(szValue, 128, "%f", m_dValue);

		if (m_eQuantity == MinArea) {
			strDescription += "Minimum area ";
			strDescription += szValue;
		} else if (m_eQuantity == MaxArea) {
			strDescription += "Maximum area ";
			strDescription += szValue;
		} else if (m_eQuantity == MinArealFraction) {
			strDescription += "Minimum areal fraction ";
			strDescription += szValue;
		} else if (m_eQuantity == MaxArealFraction) {
			strDescription += "Maximum areal fraction ";
			strDescription += szValue;
		//} else if (m_eQuantity == EastwardOrientation) {
		//	strDescription += "Eastward orientation ";
		//} else if (m_eQuantity == WestwardOrientation) {
		//	strDescription += "Westward orientation ";
		}

		Announce("%s", strDescription.c_str());
	}

	///	<summary>
	///		Verify that the specified path satisfies the threshold op.
	///	</summary>
	bool Apply(
		const SimpleGrid & grid,
		const IndicatorSet & setBlobPoints,
		const LatLonBox<double> & boxBlobDeg
	) {
		// Thresholds related to area
		if ((m_eQuantity == MinArea) ||
		    (m_eQuantity == MaxArea) ||
		    (m_eQuantity == MinArealFraction) ||
		    (m_eQuantity == MaxArealFraction)
		) {
			// Calculate the area of each blob box
			double dBoxArea = DegToRad(boxBlobDeg.lon[1]) - DegToRad(boxBlobDeg.lon[0]);
			if (boxBlobDeg.lon[0] > boxBlobDeg.lon[1]) {
				dBoxArea += 2.0 * M_PI;
			}
			dBoxArea *=
				fabs(sin(DegToRad(boxBlobDeg.lat[1])) - sin(DegToRad(boxBlobDeg.lat[0])));

			// Calculate the area and mean lat/lon of each blob
			double dBlobArea = 0.0;
			IndicatorSetIterator iterBlob = setBlobPoints.begin();
			for (; iterBlob != setBlobPoints.end(); iterBlob++) {
				dBlobArea += grid.m_dArea[*iterBlob];
			}

			// Minimum area
			if (m_eQuantity == MinArea) {
				if (dBlobArea < m_dValue) {
					return false;
				}

			// Maximum area
			} else if (m_eQuantity == MaxArea) {
				if (dBlobArea > m_dValue) {
					return false;
				}

			// Minimum areal fraction
			} else if (m_eQuantity == MinArealFraction) {
				if (dBlobArea < m_dValue * dBoxArea) {
					return false;
				}

			// Maximum areal fraction
			} else if (m_eQuantity == MaxArealFraction) {
				if (dBlobArea > m_dValue * dBoxArea) {
					return false;
				}
			}
/*
		// Thresholds related to orientation
		} else if ((m_eQuantity == EastwardOrientation) ||
		           (m_eQuantity == WestwardOrientation)
		) {
			double dNorthHemiMeanLat = 0.0;
			double dNorthHemiMeanLon = 0.0;
			double dNorthHemiMeanLon2 = 0.0;
			double dNorthHemiCoLatLon = 0.0;

			double dSouthHemiMeanLat = 0.0;
			double dSouthHemiMeanLon = 0.0;
			double dSouthHemiMeanLon2 = 0.0;
			double dSouthHemiCoLatLon = 0.0;

			// Calculate regression coefficients for this blob
			IndicatorSetIterator iterBlob = setBlobPoints.begin();
			for (; iterBlob != setBlobPoints.end(); iterBlob++) {

				double dAltLon = 0.0;
				if (dLatDeg[iterBlob->lat] > 0.0) {
					if (iterBlob->lon < boxBlob.lon[0]) {
						dAltLon = dLonDeg[iterBlob->lon] + 360.0;
					} else {
						dAltLon = dLonDeg[iterBlob->lon];
					}

					dNorthHemiMeanLat += dLatDeg[iterBlob->lat];
					dNorthHemiMeanLon += dAltLon;
					dNorthHemiMeanLon2 += dAltLon * dAltLon;
					dNorthHemiCoLatLon += dLatDeg[iterBlob->lat] * dAltLon;

				} else if (dLatDeg[iterBlob->lat] < 0.0) {
					if (iterBlob->lon < boxBlob.lon[0]) {
						dAltLon = dLonDeg[iterBlob->lon] + 360.0;
					} else {
						dAltLon = dLonDeg[iterBlob->lon];
					}

					dSouthHemiMeanLat += dLatDeg[iterBlob->lat];
					dSouthHemiMeanLon += dAltLon;
					dSouthHemiMeanLon2 += dAltLon * dAltLon;
					dSouthHemiCoLatLon += dLatDeg[iterBlob->lat] * dAltLon;
				}
			}

			double dBlobCount = static_cast<double>(setBlobPoints.size());

			dNorthHemiMeanLat /= dBlobCount;
			dNorthHemiMeanLon /= dBlobCount;
			dNorthHemiMeanLon2 /= dBlobCount;
			dNorthHemiCoLatLon /= dBlobCount;

			dSouthHemiMeanLat /= dBlobCount;
			dSouthHemiMeanLon /= dBlobCount;
			dSouthHemiMeanLon2 /= dBlobCount;
			dSouthHemiCoLatLon /= dBlobCount;

			// Calculate the slope of the regression line
			double dNorthSlopeNum =
				dNorthHemiCoLatLon
					- dNorthHemiMeanLat * dNorthHemiMeanLon;

			double dNorthSlopeDen =
				dNorthHemiMeanLon2
					- dNorthHemiMeanLon * dNorthHemiMeanLon;

			double dSouthSlopeNum =
				dSouthHemiCoLatLon
					- dSouthHemiMeanLat * dSouthHemiMeanLon;

			double dSouthSlopeDen =
				dSouthHemiMeanLon2
					- dSouthHemiMeanLon * dSouthHemiMeanLon;

			// Check orientation
			if (m_eQuantity == EastwardOrientation) {

				if (dNorthSlopeNum * dNorthSlopeDen < 0.0) {
					return false;
				}
				if (dSouthSlopeNum * dSouthSlopeDen > 0.0) {
					return false;
				}

			} else if (m_eQuantity == WestwardOrientation) {
				if (dNorthSlopeNum * dNorthSlopeDen > 0.0) {
					return false;
				}
				if (dSouthSlopeNum * dSouthSlopeDen < 0.0) {
					return false;
				}
			}
*/
		// Invalid quantity
		} else {
			_EXCEPTIONT("Invalid quantity");
		}

		return true;
	}

protected:
	///	<summary>
	///		Threshold quantity.
	///	</summary>
	ThresholdQuantity m_eQuantity;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;
};

///////////////////////////////////////////////////////////////////////////////

class RestrictRegion {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RestrictRegion() :
		m_fActive(false),
		m_dBox(),
		m_nCount(0)
	{ }

	///	<summary>
	///		Check if this RestrictRegion is active.
	///	</summary>
	bool IsActive() {
		return m_fActive;
	}

	///	<summary>
	///		Get the minimum count associated with this threshold.
	///	</summary>
	int GetMinimumCount() const {
		return m_nCount;
	}

	///	<summary>
	///		Parse a restrict region argument string, consisting of a list of
	///		5 elements.
	///	</summary>
	void Parse(
		const std::string & strArg
	) {
		int iMode = 0;

		int iLast = 0;
		for (int i = 0; i <= strArg.length(); i++) {

			if ((i == strArg.length()) || (strArg[i] == ',')) {

				std::string strSubStr =
					strArg.substr(iLast, i - iLast);

				if (iMode == 0) {
					m_dBox.lat[0] = atof(strSubStr.c_str());

				} else if (iMode == 1) {
					m_dBox.lat[1] = atof(strSubStr.c_str());

				} else if (iMode == 2) {
					m_dBox.lon[0] = atof(strSubStr.c_str());

				} else if (iMode == 3) {
					m_dBox.lon[1] = atof(strSubStr.c_str());

				} else if (iMode == 4) {
					m_nCount = atoi(strSubStr.c_str());
				}

				iMode++;
				iLast = i + 1;
			}
		}

		if (iMode != 5) {
			_EXCEPTION1("Incorrect number of arguments in --restrict_region: "
				"5 arguments expected, %i found", iMode);
		}

		if (m_nCount < 1) {
			_EXCEPTION1("Invalid value of count in --restrict_region (%i)", m_nCount);
		}

		m_fActive = true;
		m_dBox.is_null = false;

		Announce("Blob must enter region [%1.5f %1.5f] x [%1.5f %1.5f] for at least %i times",
			m_dBox.lat[0], m_dBox.lat[1], m_dBox.lon[0], m_dBox.lon[1], m_nCount);
	}

	///	<summary>
	///		Apply the operator.
	///	</summary>
	bool ContainsPoint(
		double dLatDeg,
		double dLonDeg
	) {
		return m_dBox.contains(dLatDeg, dLonDeg);
	}

protected:
	///	<summary>
	///		A flag indicating this operator is active
	///	</summary>
	bool m_fActive;

	///	<summary>
	///		A latitude-longitude box indicating the region of interest.
	///	</summary>
	LatLonBox<double> m_dBox;

	///	<summary>
	///		The count of the number of time points the blob must be in
	///		this region.
	///	</summary>
	int m_nCount;
};

///////////////////////////////////////////////////////////////////////////////

struct Node3 {
	double dX;
	double dY;
	double dZ;
};

///////////////////////////////////////////////////////////////////////////////
			
	//########################### HPC Notes (Hongyu Chen) #####################
	// Each Processor will read in the commandline information.
	// 1. Each processor will have the information of the entire input file
	//    lists.
	// 2. Then each processor will read in a number of input files accordingly
	//    and will have a fragment of the global time series
	//    But each processor will have the entire spatial dimension. In other
	//    words, the input files is a 4*6 Matrix, it means time * spatial
	//    Matrix. Each processor has a n_time/p*6 local matrix.
	// 3. When read in the benchmark file, each processor will still read in
	//    the vecInputfiles[0](the global input file array)
	// 4. Then each processor will generate the local vecAllBlobs;
	//    vecAllBlobTags; vecPrevBlobBoxesDeg; vecGlobalTimes according to
	//    their reponsible input files and do the exchange with the neighbor
	//    processors.
	// 6. The each processor will build the multigraph locally and then gather
	//    it to the root processor. (I should only gather the multigraph to P0, all others remains at local)
	// 7. The root processor will build the connectivity graph based on the
	//    gathered multigraph and then reassign tag numbers
	// 8. Then root processor will scatter out the updated vecAllBlobTags to
	//    each processor. (should The time.localTagId only update the id, not entire tag)
	// 9. Each processor will write their local vecAllBlobTags to the output
	//    file individually. The time of the output file will be identical to
	//    the input file.
	// For example, if processor 1 reads in the 1979/01/01/00~1979/01/31/23 time
	// from input file, it will also write the result of
	// 1979/01/01/00~1979/01/31/23 to the output file.
	//########################### End HPC Notes (Hongyu Chen) ##################

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
		// Initialize MPI
		int result = MPI_Init(&argc, &argv);
		if (result != MPI_SUCCESS) {
			_EXCEPTION1("The MPI routine MPI_Init failed (code %i)", result);
		}

#endif

	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {



	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Output file
	std::string strOutputFile;

	// Output file list
	std::string strOutputFileList;

	// Variable name
	std::string strVariable;

	// Output variable name
	std::string strOutputVariable;

	// Minimum blob size (in grid points)
	int nMinBlobSize;

	// Minimum duration of blob
	std::string strMinTime;

	// Tag only
	bool fTagOnly;

	// Minimum percentage of the earlier blob that needs to be covered by the later blob
	double dMinPercentOverlapPrev;

	// Maximum percentage of the earlier blob that needs to be covered by the later blob
	double dMaxPercentOverlapPrev;

	// Minimum percentage of the later blob that needs to be covered by the earlier blob
	double dMinPercentOverlapNext;

	// Maximum percentage of the later blob that needs to be covered by the earlier blob
	double dMaxPercentOverlapNext;

	// Restrict blobs to those passing through a given region
	std::string strRestrictRegion;

	// Minimum latitude for detections
	double dMinLatDeg;

	// Maximum latitude for detections
	double dMaxLatDeg;

	// Minimum longitude for detections
	double dMinLonDeg;

	// Maximum longitude for detections
	double dMaxLonDeg;

	// Merge distance (maximum distance between two blobs for them to be considered the same)
	double dMergeDistDeg;

	// Flatten data
	bool fFlatten;

	// Operate on a regional area (no periodic boundaries)
	bool fRegional;

	// Threshold commands
	std::string strThresholdCmd;

	// Latitude variable dimension name
	std::string strLatitudeName;

	// Longitude variable dimension name
	std::string strLongitudeName;

	// Time variable dimension name
	//std::string strTimeName;

	// Time variable units
	std::string strOutTimeUnits;

	// Verbose output
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "in_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputFileList, "out_list", "");
		CommandLineString(strVariable, "var", "binary_tag");
		CommandLineString(strOutputVariable, "outvar", "object_id");
		CommandLineInt(nMinBlobSize, "minsize", 1);
		CommandLineString(strMinTime, "mintime", "1");
		CommandLineBool(fTagOnly, "tagonly");
		CommandLineDoubleD(dMinPercentOverlapPrev, "min_overlap_prev", 0.0, "(%)")
		CommandLineDoubleD(dMaxPercentOverlapPrev, "max_overlap_prev", 100.0, "(%)")
		CommandLineDoubleD(dMinPercentOverlapNext, "min_overlap_next", 0.0, "(%)")
		CommandLineDoubleD(dMaxPercentOverlapNext, "max_overlap_next", 100.0, "(%)")
		CommandLineDouble(dMergeDistDeg, "merge_dist", 0.0); 
		CommandLineStringD(strRestrictRegion, "restrict_region", "", "(lat0,lat1,lon0,lon1,count)");
		CommandLineBool(fRegional, "regional");
		CommandLineDouble(dMinLatDeg, "minlat", -90.0);
		CommandLineDouble(dMaxLatDeg, "maxlat", 90.0);
		CommandLineDouble(dMinLonDeg, "minlon", 0.0);
		CommandLineDouble(dMaxLonDeg, "maxlon", 360.0);
		CommandLineBool(fFlatten, "flatten");
		CommandLineString(strLatitudeName, "latname","lat");
		CommandLineString(strLongitudeName, "lonname","lon");
		CommandLineString(strOutTimeUnits,"outtimeunits","");
		CommandLineString(strThresholdCmd, "thresholdcmd", "");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Check input
	if ((strInputFile == "") && (strInputFileList == "")) {
		_EXCEPTIONT("No input file (--in) or (--in_list) specified");
	}
	if ((strInputFile != "") && (strInputFileList != "")) {
		_EXCEPTIONT("Only one of input file (--in) or (--in_list) allowed");
	}

	// Check output
	if ((strOutputFile == "") && (strOutputFileList == "")) {
		_EXCEPTIONT("No output file (--out) or (--out_list) specified");
	}
	if ((strOutputFile != "") && (strOutputFileList != "")) {
		_EXCEPTIONT("Only one of input file (--in) or (--in_list) allowed");
	}

	// Check variable
	if (strVariable == "") {
		_EXCEPTIONT("No variable name (--var) specified");
	}

	// Register variable
	int varix = varreg.FindOrRegister(strVariable);

	// Check output variable
	if (strOutputVariable.length() == 0) {
		strOutputVariable = strVariable + "tag";
	}

	// Input file list
	FilenameList vecInputFiles;

	if (strInputFile != "") {
		vecInputFiles.push_back(strInputFile);
	}
	if (strInputFileList != "") {
		vecInputFiles.FromFile(strInputFileList, false);
	}

	int nFiles = vecInputFiles.size();

	// Output file list
	FilenameList vecOutputFiles;

	if (strOutputFile != "") {
		vecOutputFiles.push_back(strOutputFile);
	}
	if (strOutputFileList != "") {
		vecOutputFiles.FromFile(strOutputFileList, false);

		if (vecOutputFiles.size() != vecInputFiles.size()) {
			_EXCEPTION2("Mismatch in number of rows of --in_list (%lu) and --out_list (%lu)",
				vecInputFiles.size(), vecOutputFiles.size());
		}
	}

	// Parse --mintime
	if (fTagOnly && (strMinTime != "1")) {
		_EXCEPTIONT("--tagonly and --mintime cannot be combined");
	}
	int nMinTime = 1;
	double dMinTimeSeconds = 0.0;
	if (STLStringHelper::IsIntegerIndex(strMinTime)) {
		nMinTime = std::stoi(strMinTime);
		if (nMinTime < 1) {
			_EXCEPTIONT("Invalid value of --mintime; expected positive integer or time delta");
		}

	} else {
		Time timeMinTime;
		timeMinTime.FromFormattedString(strMinTime);
		dMinTimeSeconds = timeMinTime.AsSeconds();
		if (dMinTimeSeconds <= 0.0) {
			_EXCEPTIONT("Invalid value for --mintime; expected positive integer or positive time delta");
		}
	}

	// Convert percent overlap into a decimal value
	if ((dMinPercentOverlapPrev < 0.0) || (dMinPercentOverlapPrev > 100.0)) {
		_EXCEPTIONT("--min_overlap_prev must take on values between 0 and 100");
	}
	if ((dMaxPercentOverlapPrev < 0.0) || (dMaxPercentOverlapPrev > 100.0)) {
		_EXCEPTIONT("--max_overlap_prev must take on values between 0 and 100");
	}
	if ((dMinPercentOverlapNext < 0.0) || (dMinPercentOverlapNext > 100.0)) {
		_EXCEPTIONT("--min_overlap_next must take on values between 0 and 100");
	}
	if ((dMaxPercentOverlapNext < 0.0) || (dMaxPercentOverlapNext > 100.0)) {
		_EXCEPTIONT("--max_overlap_next must take on values between 0 and 100");
	}

	dMinPercentOverlapPrev /= 100.0;
	dMaxPercentOverlapPrev /= 100.0;
	dMinPercentOverlapNext /= 100.0;
	dMaxPercentOverlapNext /= 100.0;

	// Merge nearby blobs
	if ((dMergeDistDeg < 0.0) || (dMergeDistDeg > 180.0)) {
		_EXCEPTION1("--merge_dist must be between 0.0 and 180.0 (given %1.5f)", dMergeDistDeg);
	}

	// Ensure --minlat and --maxlat values are in range
	if ((dMinLatDeg < -90.0) || (dMinLatDeg > 90.0)) {
		_EXCEPTION1("--minlat must be in range [-90,90] (given %1.5f)", dMinLatDeg);
	}
	if ((dMaxLatDeg < -90.0) || (dMaxLatDeg > 90.0)) {
		_EXCEPTION1("--maxlat must be in range [-90,90] (given %1.5f)", dMaxLatDeg);
	}

	// Convert longitude values to standard range
	if (dMinLonDeg != 0.0) {
		dMinLonDeg = LonDegToStandardRange(dMinLonDeg);
	}
	if (dMaxLonDeg != 360.0) {
		dMaxLonDeg = LonDegToStandardRange(dMaxLonDeg);
	}
	if (dMinLonDeg == dMaxLonDeg) {
		_EXCEPTION2("--minlon (%1.5f) and --maxlon (%1.5f) must not be equal",
			dMinLonDeg, dMaxLonDeg);
	}
	if (dMinLatDeg >= dMaxLonDeg) {
		_EXCEPTION2("--minlat (%1.5f) must be strictly less than --maxlat (%1.5f)",
			dMinLatDeg, dMaxLatDeg);
	}

	// Parse the restrict region string
	RestrictRegion opRestrictRegion;
	if (strRestrictRegion != "") {
		opRestrictRegion.Parse(strRestrictRegion);
	}

	// Parse the threshold string
	std::vector<BlobThresholdOp> vecThresholdOp;

	if (strThresholdCmd != "") {
		AnnounceStartBlock("Parsing thresholds");

		int iLast = 0;
		for (int i = 0; i <= strThresholdCmd.length(); i++) {

			if ((i == strThresholdCmd.length()) ||
			    (strThresholdCmd[i] == ';')
			) {
				std::string strSubStr =
					strThresholdCmd.substr(iLast, i - iLast);

				int iNextOp = (int)(vecThresholdOp.size());
				vecThresholdOp.resize(iNextOp + 1);
				vecThresholdOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}


	// Define the SimpleGrid
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");

		// Load in the benchmark file
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[0]);//[HC] What is stored at vecInputFiles[0]?

		_ASSERT(vecNcFiles.size() > 0);

		if (vecNcFiles.size() < 1) {
			_EXCEPTIONT("No data files specified; unable to generate grid");
		}

		grid.GenerateLatitudeLongitude(
			vecNcFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		_ASSERT(grid.m_nGridDim.size() == 2);

		AnnounceEndBlock("Done");
	}

	// Check for area
	if (!grid.HasAreas()) {
		_EXCEPTIONT("SimpleGrid has no area information (needed for StitchBlobs)");
	}

	// Get time dimension over all files
	// strOutTimeUnits is either predetermined or set at the command line
	AnnounceStartBlock("Concatenating times");
	NcType nctypeTime;
	std::vector< std::pair<int, int> > vecGlobalTimeIxToFileTimeIx;

	std::vector< std::vector<Time> > vecGlobalTimes;
	vecGlobalTimes.resize(vecOutputFiles.size());

	
	#if defined(TEMPEST_MPIOMP)
		//============================= Spread files across nodes=================================
		// Note: if vecInputFiles.size() < total processor numbers, only <vecInputFiles.size()>
		//   number of processor will be used.
		// Calculate how many files each processor should process
		int processorResponsibalForFile_UB;
		int processorResponsibalForFile_LB;
		MPI_Comm MPI_REAL_COMM;
		int tempMPIRank;
		int tempMPISize;
		int nMPIRank;
		int nMPISize;
		int valid_flag;//Used to indicate the current rank is valid or not.

		MPI_Comm_rank(MPI_COMM_WORLD, &tempMPIRank);	
		MPI_Comm_size(MPI_COMM_WORLD, &tempMPISize);

		if (tempMPISize > 1) {
			// Assign each processor with corresponding file index. The remainder will be
			// evenly spread across the first few processors
			
			int avgNumFiles;
			avgNumFiles = vecInputFiles.size() / tempMPISize; 
			int localFileNumber;
			int remainder = vecInputFiles.size() % tempMPISize;
			if (remainder != 0) {
				if (tempMPIRank < remainder ) {
					localFileNumber = avgNumFiles + 1;
					processorResponsibalForFile_LB = localFileNumber * tempMPIRank;
					processorResponsibalForFile_UB = localFileNumber * (tempMPIRank + 1);
				} else {
					processorResponsibalForFile_LB = (avgNumFiles + 1) * remainder + (tempMPIRank - remainder) * avgNumFiles;
					processorResponsibalForFile_UB = (avgNumFiles + 1) * remainder + (tempMPIRank - remainder + 1) * avgNumFiles;					
				}
				
			} else {
				processorResponsibalForFile_LB = avgNumFiles * tempMPIRank;
				processorResponsibalForFile_UB = avgNumFiles * (tempMPIRank + 1);

			}

			if (vecInputFiles.size() < tempMPISize) {
				if (tempMPIRank >= vecInputFiles.size() ) {
					valid_flag = 0;
				} else {
					// Create a new communicator with only at most vecInputFiles.size() size
					valid_flag = 1;					
				}

			} else {
				valid_flag = 1;		
			}
			MPI_Comm_split(MPI_COMM_WORLD, valid_flag, tempMPIRank, &MPI_REAL_COMM);
			MPI_Comm_rank(MPI_REAL_COMM, &nMPIRank);	
			MPI_Comm_size(MPI_REAL_COMM, &nMPISize);

		} else {
			processorResponsibalForFile_LB = 0;
			processorResponsibalForFile_UB = vecInputFiles.size();
			MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);	
			MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

		}
	#endif

	for (int f = 0; f < vecInputFiles.size(); f++){
	
		#if defined(TEMPEST_MPIOMP)
			if ((f >= processorResponsibalForFile_UB) || f < processorResponsibalForFile_LB) {
				continue;
			}
		#endif 

		// Load in the time variable from all files
		NcFileVector vecNcFiles;//also known as the local vecNcFiles when MPI is enabled.
		vecNcFiles.ParseFromString(vecInputFiles[f]);
		_ASSERT(vecNcFiles.size() > 0);

		// Get the time variable
		NcVar * varTime = NcGetTimeVariable(*(vecNcFiles[0]));
		if (varTime == NULL) {
			_EXCEPTION1("File \"%s\" does not contain \"time\" variable",
				vecNcFiles.GetFilename(0).c_str());
		}

		nctypeTime = varTime->type();

		// Get time units (if not specified on command line)
		if (strOutTimeUnits == "") {
			NcAtt * attTime = varTime->get_att("units");
			if (attTime == NULL) {
				_EXCEPTION1("File \"%s\" missing \"time:units\" attribute",
					vecNcFiles.GetFilename(0).c_str());
			}

			strOutTimeUnits = attTime->as_string(0);
		}

		// Load in CF-compliant time data
		const NcTimeDimension & vecTimes = vecNcFiles.GetNcTimeDimension(0);
		if (vecTimes.size() == 0) {
			_EXCEPTION1("WARNING: File group does not contain any time data (%s)",
				vecInputFiles[f].c_str());
		}
		//ReadCFTimeDataFromNcFile(
		//	vecNcFiles[0],
		//	vecInputFiles[f].c_str(),
		//	vecTimes,
		//	true);

		//std::cout << vecTimes[0].ToString() << " " << vecTimes[0].GetCalendarType() << std::endl;

		if (vecOutputFiles.size() == 1) {
			for (int t = 0; t < vecTimes.size(); t++) {
				int iGlobalTime = vecGlobalTimes[0].size();
				vecGlobalTimes[0].push_back(vecTimes[t]);
				vecGlobalTimeIxToFileTimeIx.push_back( std::pair<int,int>(0,iGlobalTime) );
			}
		} else {
			for (int t = 0; t < vecTimes.size(); t++) {
				vecGlobalTimes[f].push_back(vecTimes[t]);
				vecGlobalTimeIxToFileTimeIx.push_back( std::pair<int,int>(f,t) );
			}
		}
	}

	int nGlobalTimes = 0;
	for (int f = 0; f < vecGlobalTimes.size(); f++) {
		#if defined(TEMPEST_MPIOMP)
			if ((f >= processorResponsibalForFile_UB) || f < processorResponsibalForFile_LB) {
				continue;
			}
		#endif 
		nGlobalTimes += vecGlobalTimes[f].size();
	}
	_ASSERT(nGlobalTimes > 0);
	_ASSERT(nGlobalTimes == vecGlobalTimeIxToFileTimeIx.size());

	AnnounceEndBlock("Done");

	///////////////////////////////////////////////////////////////////////////
	// Build the set of nodes at each time contained in each blob
	///////////////////////////////////////////////////////////////////////////

	// Build blobs at each time level
	AnnounceStartBlock("Building blob set at each time level");

	// Set of nodes at each time contained in each blob
	std::vector< std::vector<IndicatorSet> > vecAllBlobs;//Sending and Receiving Blobs to nearby processors [Halo Var]
	vecAllBlobs.resize(nGlobalTimes);

	// Bounding boxes at each time for each blob
	std::vector< std::vector< LatLonBox<double> > > vecAllBlobBoxesDeg;//[Halo Var]
	vecAllBlobBoxesDeg.resize(nGlobalTimes);

	// Time index across all (local) files
	int iTime = 0;

	// Loop through all files
	int startIndx = 0;//The starting index for looping through all files
	#if defined(TEMPEST_MPIOMP) 
		//If MPI is enabled, then modify the nFiles to the local file numbers
		startIndx = processorResponsibalForFile_LB;
		nFiles = processorResponsibalForFile_UB;		
	#endif 

	for (int f = startIndx; f < nFiles; f++) {
		// Clear existing data in the register
		varreg.UnloadAllGridData();

		// Load in the benchmark file
		NcFileVector vecNcFiles; //[HC] Why it is called "vecNcFiles" what information does it have? Coz u can have multiple files at one line
		vecNcFiles.ParseFromString(vecInputFiles[f]);
		_ASSERT(vecNcFiles.size() > 0);
		


		// Number of times in this input file
		NcDim * dimTimeInput = vecNcFiles[0]->get_dim("time");
		if (dimTimeInput == NULL) {
			_EXCEPTION1("No dimension \"time\" in file \"%s\"",
				vecNcFiles.GetFilename(0).c_str());
		}
		int nLocalTimes = dimTimeInput->size();

		// Loop through all loacal time
		for (int t = 0; t < nLocalTimes; t++, iTime++) {
			
			// Get the current patch vector
			std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[iTime];
			

			std::vector< LatLonBox<double> > & vecBlobBoxesDeg = vecAllBlobBoxesDeg[iTime];


			// KD-trees for points in blob at this time step (used for --merge_dist)
			std::vector<kdtree*> vecBlobTrees;


			std::vector<IndicatorSet> vecBlobPerimeters;


			// New announcement block for timestep
			if (vecGlobalTimes.size() == 1) {
				_ASSERT((iTime >= 0) && (iTime < vecGlobalTimes[0].size()));
				AnnounceStartBlock("Time %i (%s)", iTime, 
					vecGlobalTimes[0][iTime].ToString().c_str());
			} else {
				_ASSERT((t >= 0) && (t < vecGlobalTimes[f].size()));//allow t = vecGlobalTimes[f].size())
				AnnounceStartBlock("Time %i (%s)", iTime,
					vecGlobalTimes[f][t].ToString().c_str());
			}

			// Load the search variable data [HC] What does the following block do exactly, how does it read through the time slice? What does the following i mean?
			Variable & var = varreg.Get(varix);
			vecNcFiles.SetConstantTimeIx(t);
			// [HC] Check in every loop how does thedata change
			// Memeory scale by the number of the threads used
			var.LoadGridData(varreg, vecNcFiles, grid); //
			const DataArray1D<float> & dataIndicator = var.GetData();
/*
			float dChecksum = 0.0;
			for (int i = 0; i < dataState.GetRows(); i++) {
				dChecksum += dataState[i];
			}
			std::cout << dChecksum << std::endl;
*/

			// Buffer variables used for looping through indicators
			IndicatorSet setIndicators;

			// Insert all detected locations into set
			// (accounting for points out of range)
			if ((dMinLatDeg != -90.0) || (dMaxLatDeg != 90.0) ||
			    (dMinLonDeg != 0.0) || (dMaxLonDeg != 360.0)
			) {
				LatLonBox<double> boxBoundsDeg(!fRegional);

				if (fRegional) {
					boxBoundsDeg.lon[0] = dMinLonDeg;
					boxBoundsDeg.lon[1] = dMaxLonDeg;

					if (dMinLonDeg > dMaxLonDeg) {
						_EXCEPTIONT("On regional grids, --minlon must be less than --maxlon");
					}

				} else {
					boxBoundsDeg.lon[0] = LonDegToStandardRange(dMinLonDeg);
					boxBoundsDeg.lon[1] = LonDegToStandardRange(dMaxLonDeg);
				}
				boxBoundsDeg.lat[0] = dMinLatDeg;
				boxBoundsDeg.lat[1] = dMaxLatDeg;

				for (int i = 0; i < grid.GetSize(); i++) {
					if (dataIndicator[i] == 0.0f) {
						continue;
					}

					double dLonDeg = RadToDeg(grid.m_dLon[i]);
					double dLatDeg = RadToDeg(grid.m_dLat[i]);

					dLonDeg = LonDegToStandardRange(dLonDeg);

					if (!boxBoundsDeg.contains(dLatDeg, dLonDeg)) {
						continue;
					}

					setIndicators.insert(i);
				}

			// Insert all detected locations into set
			// (no bounds checking)
			} else {
				for (int i = 0; i < grid.GetSize(); i++) {
					if (dataIndicator[i] != 0.0f) {
						setIndicators.insert(i);
					}
				}
			}

			Announce("Finding blobs (%i tagged points)", setIndicators.size());

			// Backup of original indicator set
			IndicatorSet setIndicatorsBackup = setIndicators;

			// Rejections due to insufficient node count
			int nRejectedMinSize = 0;

			// Rejections due to a given threshold
			DataArray1D<int> nRejectedThreshold(vecThresholdOp.size());

			// Find all contiguous patches
			for (; setIndicators.size() != 0;) {

				// Next starting location
				int ixNode = *(setIndicators.begin());

				// Current patch
				int ixBlob = vecBlobs.size();
				vecBlobs.resize(ixBlob+1);
				vecBlobBoxesDeg.resize(ixBlob+1, LatLonBox<double>(fRegional));

				if (dMergeDistDeg > 0.0) {
					vecBlobPerimeters.resize(ixBlob+1);
					vecBlobTrees.resize(ixBlob+1);
					vecBlobTrees[ixBlob] = kd_create(3);
				}

				// Initialize bounding box
				LatLonBox<double> & boxDeg = vecBlobBoxesDeg[ixBlob];
				boxDeg.lat[0] = RadToDeg(grid.m_dLat[ixNode]);
				boxDeg.lat[1] = RadToDeg(grid.m_dLat[ixNode]);
				boxDeg.lon[0] = LonDegToStandardRange(RadToDeg(grid.m_dLon[ixNode]));
				boxDeg.lon[1] = LonDegToStandardRange(RadToDeg(grid.m_dLon[ixNode]));

				//printf("=================================== BLOB %i\n", ixBlob);
				// Find all connecting nodes in patch
				IndicatorSet setNeighbors;
				setNeighbors.insert(ixNode);
				while (setNeighbors.size() != 0) {
					ixNode = *(setNeighbors.begin());
					setNeighbors.erase(setNeighbors.begin());

					// This node is already included in the blob
					if (vecBlobs[ixBlob].find(ixNode) != vecBlobs[ixBlob].end()) {
						//printf("..%i already in set\n", ixNode);
						continue;
					}

					// This node has not been tagged
					IndicatorSetIterator iterIndicator = setIndicators.find(ixNode);
					if (iterIndicator == setIndicators.end()) {
						//printf("..%i has not been tagged\n", ixNode);
						continue;
					}

					// Remove this from the set of available indicators
					setIndicators.erase(iterIndicator);

					// Insert the node into the blob
					vecBlobs[ixBlob].insert(ixNode);

					// Update bounding box
					boxDeg.insert(
						RadToDeg(grid.m_dLat[ixNode]),
						LonDegToStandardRange(RadToDeg(grid.m_dLon[ixNode])));

					// Insert neighbors
					bool fPerimeter = false;
					//printf("..%i has %lu neighbors\n", ixNode, grid.m_vecConnectivity[ixNode].size());
					for (int i = 0; i < grid.m_vecConnectivity[ixNode].size(); i++) {
						int ixNeighbor = grid.m_vecConnectivity[ixNode][i];
						//printf("....neighbor %i %1.5f ", ixNeighbor, dataIndicator[ixNeighbor]);
						//if (setIndicatorsBackup.find(ixNeighbor) != setIndicatorsBackup.end()) {
						//	printf("IN INDICATORSET\n");
						//} else {
						//	printf("\n");
						//}

						// Perimeter point
						if (setIndicatorsBackup.find(ixNeighbor) == setIndicatorsBackup.end()) {
							fPerimeter = true;
						} else {
							setNeighbors.insert(ixNeighbor);
						}
					}
					if (fPerimeter && (vecBlobPerimeters.size() != 0)) {
						vecBlobPerimeters[ixBlob].insert(ixNode);
						double dX, dY, dZ;
						RLLtoXYZ_Rad(
							grid.m_dLon[ixNode],
							grid.m_dLat[ixNode],
							dX, dY, dZ);
						kd_insert3(vecBlobTrees[ixBlob], dX, dY, dZ, NULL);
					}
				}
			}

			// setIndicatorsBackup no longer needed
			setIndicatorsBackup.clear();

			// Merge blobs
			if (vecBlobTrees.size() != 0) {

				AnnounceStartBlock("Merging blobs (from %lu blobs)", vecBlobs.size());

				Announce("Building merge graph");

				// Use squared Cartesian distance for comparisons
				double dMergeDistChordLength2 = ChordLengthFromGreatCircleDistance_Deg(dMergeDistDeg);
				dMergeDistChordLength2 *= dMergeDistChordLength2;

				// Build a set of all pairs that need to be merged
				std::vector< std::set<int> > vecMergeGraph(vecBlobs.size());
				for (int ixBlob = 0; ixBlob < vecBlobs.size(); ixBlob++) {
					//Announce("Blob %i has %lu total points, %lu perimeter points",
					//	ixBlob,
					//	vecBlobs[ixBlob].size(),
					//	vecBlobPerimeters[ixBlob].size());

					// Convert perimeter points to xyz
					std::vector<Node3> vecPerimeterXYZ;
					vecPerimeterXYZ.reserve(vecBlobPerimeters[ixBlob].size());
					for (auto it = vecBlobPerimeters[ixBlob].begin(); it != vecBlobPerimeters[ixBlob].end(); it++) {
						Node3 xyzPerim;
						RLLtoXYZ_Rad(
							grid.m_dLon[*it],
							grid.m_dLat[*it],
							xyzPerim.dX, xyzPerim.dY, xyzPerim.dZ);
/*
						if ((RadToDeg(grid.m_dLon[*it]) < 187.25) &&
							(RadToDeg(grid.m_dLon[*it]) > 178.00) &&
							(RadToDeg(grid.m_dLat[*it]) < -24.0) &&
							(RadToDeg(grid.m_dLat[*it]) > -50.0)
							//(RadToDeg(grid.m_dLat[*it]) < -50.0) &&
							//(RadToDeg(grid.m_dLat[*it]) > -57.0)
						) {
							printf("%i %i\n", ixBlob, *it);
						}
*/
						vecPerimeterXYZ.push_back(xyzPerim);
					}

					// Search all perimeter points of this blob against perimeter points of other blobs
					for (int ixBlobNext = ixBlob+1; ixBlobNext != vecBlobs.size(); ixBlobNext++) {
						for (int i = 0; i < vecPerimeterXYZ.size(); i++) {

							const Node3 & xyzPerim = vecPerimeterXYZ[i];
							kdres * kdr = kd_nearest3(vecBlobTrees[ixBlobNext], xyzPerim.dX, xyzPerim.dY, xyzPerim.dZ);
							if (kdr == NULL) {
								_EXCEPTION4("Logic error: Cannot find nearest point to blob %i from (%1.5e %1.5e %1.5e)",
									ixBlobNext, xyzPerim.dX, xyzPerim.dY, xyzPerim.dZ);
							}
							if (kd_res_size(kdr) == 0) {
								_EXCEPTION4("Logic error: Cannot find nearest point to blob %i from (%1.5e %1.5e %1.5e)",
									ixBlobNext, xyzPerim.dX, xyzPerim.dY, xyzPerim.dZ);
							}
							double dX1, dY1, dZ1;
							kd_res_item3(kdr, &dX1, &dY1, &dZ1);
							kd_res_free(kdr);

							double dDX = xyzPerim.dX - dX1;
							double dDY = xyzPerim.dY - dY1;
							double dDZ = xyzPerim.dZ - dZ1;
							double dDist = dDX * dDX + dDY * dDY + dDZ * dDZ;
/*
							if ((ixBlob == 10) && (ixBlobNext == 15)) {
								printf("%1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.10f\n",
									xyzPerim.dX, xyzPerim.dY, xyzPerim.dZ,
									dX1, dY1, dZ1, sqrt(dDist));
							}
*/
							if (dDist < dMergeDistChordLength2) {
								//printf("Merge %i %i\n", ixBlob, ixBlobNext);
								vecMergeGraph[ixBlob].insert(ixBlobNext);
								vecMergeGraph[ixBlobNext].insert(ixBlob);
								break;
							}
						}
					}
				}

				// Actually merge the blobs
				Announce("Retagging merged blobs");

				for (int ixBlob = 0; ixBlob < vecBlobs.size(); ixBlob++) {

					// Find the ultimate blob id of each blob via graph search for minimum id
					int ixFinalBlobId = ixBlob;

					std::set<int> setVisited;
					std::queue<int> queueToVisit;
					queueToVisit.push(ixBlob);

					while (queueToVisit.size() != 0) {
						int ixCurrentBlobId = queueToVisit.front();
						queueToVisit.pop();

						if (setVisited.find(ixCurrentBlobId) != setVisited.end()) {
							continue;
						}
						setVisited.insert(ixCurrentBlobId);
						//std::cout << ixBlob << " " << ixCurrentBlobId << std::endl;

						if (ixCurrentBlobId < ixFinalBlobId) {
							ixFinalBlobId = ixCurrentBlobId;
						}

						for (auto it = vecMergeGraph[ixCurrentBlobId].begin(); it != vecMergeGraph[ixCurrentBlobId].end(); it++) {
							queueToVisit.push(*it);
						}
					}

					// Move points into final blob tag
					if (ixFinalBlobId != ixBlob) {
						//printf("Merging blob %i into blob %i\n", ixBlob, ixFinalBlobId);
						for (auto it = vecBlobs[ixBlob].begin(); it != vecBlobs[ixBlob].end(); it++) {
							vecBlobs[ixFinalBlobId].insert(*it);
							vecBlobBoxesDeg[ixFinalBlobId].insert(
								RadToDeg(grid.m_dLat[*it]),
								LonDegToStandardRange(RadToDeg(grid.m_dLon[*it])));
						}
						vecBlobs[ixBlob].clear();
					}
				}

				// Clean up
				for (int ixBlob = 0; ixBlob < vecBlobs.size(); ixBlob++) {
					if (vecBlobs[ixBlob].size() == 0) {
						vecBlobs.erase(vecBlobs.begin() + ixBlob);
						vecBlobBoxesDeg.erase(vecBlobBoxesDeg.begin() + ixBlob);
						ixBlob--;
					}
				}
				for (int ixBlob = 0; ixBlob < vecBlobTrees.size(); ixBlob++) {
					kd_free(vecBlobTrees[ixBlob]);
				}
				vecBlobTrees.clear();
				vecBlobPerimeters.clear();

				AnnounceEndBlock(NULL);
			}

			// Filter blobs
			if ((nMinBlobSize > 1) || (vecThresholdOp.size() != 0)) {
				AnnounceStartBlock("Applying filters (from %lu blobs)", vecBlobs.size());

				for (int ixBlob = 0; ixBlob < vecBlobs.size(); ixBlob++) {

					// Check blob size
					if (vecBlobs[ixBlob].size() < nMinBlobSize) {
						nRejectedMinSize++;
						vecBlobs.erase(vecBlobs.begin() + ixBlob);
						vecBlobBoxesDeg.erase(vecBlobBoxesDeg.begin() + ixBlob);
						ixBlob--;

					// Check other thresholds
					} else {
						for (int x = 0; x < vecThresholdOp.size(); x++) {

							bool fSatisfies =
								vecThresholdOp[x].Apply(
									grid,
									vecBlobs[ixBlob],
									vecBlobBoxesDeg[ixBlob]);

							if (!fSatisfies) {
								nRejectedThreshold[x]++;
								vecBlobs.erase(vecBlobs.begin() + ixBlob);
								vecBlobBoxesDeg.erase(vecBlobBoxesDeg.begin() + ixBlob);
								ixBlob--;
								break;
							}
						}
					}
				}

				Announce("Rejected (min size): %i", nRejectedMinSize);
				for (int x = 0; x < vecThresholdOp.size(); x++) {
					Announce("Rejected (threshold %i): %i",
						x, nRejectedThreshold[x]);
				}

				AnnounceEndBlock(NULL);
			}

			AnnounceStartBlock("Blobs detected: %i", vecBlobs.size());

			if (fVerbose) {
				for (int p = 0; p < vecBlobBoxesDeg.size(); p++) {
					Announce("Blob %i [%1.5f, %1.5f] x [%1.5f, %1.5f]",
						p+1,
						vecBlobBoxesDeg[p].lat[0],
						vecBlobBoxesDeg[p].lat[1],
						vecBlobBoxesDeg[p].lon[0],
						vecBlobBoxesDeg[p].lon[1]);
				}
			}

			AnnounceEndBlock(NULL);

			AnnounceEndBlock("Done");
		}
	}

	AnnounceEndBlock("Done");

	///////////////////////////////////////////////////////////////////////////
	// Stitch blobs together in time using graph search
	///////////////////////////////////////////////////////////////////////////




	AnnounceStartBlock("Assign local tags to each blob");

	// Tags for each blob at each time slice
	std::vector< std::vector<Tag> > vecAllBlobTags;
	vecAllBlobTags.resize(nGlobalTimes);

	// Next available patch tag
	Tag tagNextBlob(1, 0);

	// Give blob tags to the initial set of blobs
	std::set<Tag> setAllTags;

	for (int t = 0; t < nGlobalTimes; t++) {
		vecAllBlobTags[t].resize(vecAllBlobs[t].size());

		tagNextBlob.id = 0;
		for (int p = 0; p < vecAllBlobTags[t].size(); p++) {
			vecAllBlobTags[t][p] = tagNextBlob;
			setAllTags.insert(tagNextBlob);

			tagNextBlob.id++;
		}
		tagNextBlob.time++;
	}
	AnnounceEndBlock("Done");

	//================================Actual Parallization Starts===================================
	//1. Exchang the vecAllBlobs; vecAllBlobTags; vecPrevBlobBoxesDeg; vecGlobalTimes
	//==============================================================================================

#if defined(TEMPEST_MPIOMP)

	//We still need the original unexchanged data for these two variables later
	std::unique_ptr<GlobalTimesExchangeOp> MPI_exchangedGlobalTimes;
	std::unique_ptr<TagExchangeOP>       MPI_exchangedTags;
	std::unique_ptr<BlobsExchangeOp>       MPI_exchangedBlobs;
	if (nMPISize > 1 && valid_flag) {

    // Now allocate the objects.
		MPI_exchangedGlobalTimes.reset(new GlobalTimesExchangeOp(MPI_REAL_COMM, vecGlobalTimes,
				processorResponsibalForFile_LB,
				processorResponsibalForFile_UB));
		MPI_exchangedTags.reset(new TagExchangeOP(MPI_REAL_COMM, vecAllBlobTags));
		MPI_exchangedBlobs.reset(new BlobsExchangeOp(MPI_REAL_COMM, vecAllBlobs));

		// Exchange vecGlobalTimes.
		MPI_exchangedGlobalTimes->StartExchange();
		MPI_exchangedGlobalTimes->EndExchange();
		vecGlobalTimes = MPI_exchangedGlobalTimes->GetExchangedVecGlobalTimes();

		// Exchange vecAllBlobTags.
		MPI_exchangedTags->StartExchange();
		MPI_exchangedTags->EndExchange();
		vecAllBlobTags = MPI_exchangedTags->GetExchangedVecAllBlobTags();

		// Exchange vecAllBlobs.
		MPI_exchangedBlobs->StartExchange();
		MPI_exchangedBlobs->EndExchange();
		vecAllBlobs = MPI_exchangedBlobs->GetExchangedVecAllBlobs();

		//Exchange vecPrevBlobBoxesDeg
		BlobBoxesDegExchangeOP MPI_exchangedBlobBoxesDeg(MPI_REAL_COMM, vecAllBlobBoxesDeg);
		MPI_exchangedBlobBoxesDeg.StartExchange();
		MPI_exchangedBlobBoxesDeg.EndExchange();
		vecAllBlobBoxesDeg = MPI_exchangedBlobBoxesDeg.GetExchangedVecAllBlobBoxesDeg();
		
		//Make sure all processors finish the exchange
		MPI_Barrier(MPI_REAL_COMM);

	}
#endif



	MapGraph multimapTagGraph;

	// Set of Tags that are within the restrict region
	std::set<Tag> setRestrictRegion;

	AnnounceStartBlock("Building connectivity graph");

	// Loop through all remaining time steps
	int iFileLocal = 0;
	int iTimeLocal = 0;
	
	#if defined(TEMPEST_MPIOMP)
		//recalculate the nGlobalTimes here based on the updated vecGlobalTimes file (will be inverted after the connectivity graph is built)
		//[HC] Isn't the original one the same as the following one?[Answer: it should be, but need double check]
		int original_nGlobalTimes = nGlobalTimes;
			
		if (nMPISize > 1 && valid_flag) {
			nGlobalTimes = 0;
			iFileLocal = processorResponsibalForFile_LB;
			for (int f = processorResponsibalForFile_LB; f < processorResponsibalForFile_UB; f++) {
				nGlobalTimes += vecGlobalTimes[f].size();
			}
		}
	
	#endif
	for (int t = 1; t < nGlobalTimes; t++) {

		// New announcement block for timestep
		if (vecGlobalTimes.size() == 1) {
			_ASSERT((t >= 0) && (t < vecGlobalTimes[0].size()));
			Announce("Time %i (%s)", t,
				vecGlobalTimes[0][t].ToString().c_str());
		} else {
			iTimeLocal++;
			_ASSERT(iFileLocal < vecGlobalTimes.size());
			if (iTimeLocal >= vecGlobalTimes[iFileLocal].size()) {
				iTimeLocal = 0;
				iFileLocal++;
				#if defined(TEMPEST_MPIOMP)
					_ASSERT(iFileLocal < processorResponsibalForFile_UB);
				#else
					_ASSERT(iFileLocal < vecGlobalTimes.size());
				#endif
				
			}
			_ASSERT(iTimeLocal < vecGlobalTimes[iFileLocal].size());
			Announce("Time %i (%s)", t,
				vecGlobalTimes[iFileLocal][iTimeLocal].ToString().c_str());
		}

		// Get the current blob vector
		const std::vector<Tag> & vecPrevBlobTags = vecAllBlobTags[t-1];

		std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		const std::vector<IndicatorSet> & vecPrevBlobs = vecAllBlobs[t-1];

		const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

		const std::vector< LatLonBox<double> > & vecPrevBlobBoxesDeg
			= vecAllBlobBoxesDeg[t-1];

		const std::vector< LatLonBox<double> > & vecBlobBoxesDeg
			= vecAllBlobBoxesDeg[t];

		// Determine overlaps between these blobs and previous blobs
		vecBlobTags.resize(vecBlobs.size());
		int nCountRemove = 0;
		for (int p = 0; p < vecBlobTags.size(); p++) {

			const LatLonBox<double> & boxP = vecBlobBoxesDeg[p];

			// Find overlap with bounding boxes at previous time
			for (int q = 0; q < vecPrevBlobTags.size(); q++) {

				const LatLonBox<double> & boxQ = vecPrevBlobBoxesDeg[q];

				// Check if tagonly
				if (fTagOnly) {
					continue;
				}

				// Check if bounding boxes overlap
				if (!boxP.overlaps(boxQ)) {
					continue;
				}

				// Verify that at least one node overlaps between blobs
				{
					bool fHasOverlapNode = false;
					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						if (vecPrevBlobs[q].find(*iter) !=
							vecPrevBlobs[q].end()
						) {
							fHasOverlapNode = true;
							break;
						}
					}

					if (!fHasOverlapNode) {
						continue;
					}
				}

				

				// Verify that blobs meet percentage overlap criteria
				// (as a percentage of the current blob)
				if ((dMinPercentOverlapNext != 0.0) ||
				    (dMaxPercentOverlapNext != 1.0)
				) {
					double dCurrentArea = 0.0;
					double dOverlapArea = 0.0;

					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						dCurrentArea += grid.m_dArea[*iter];
						if (vecPrevBlobs[q].find(*iter) != vecPrevBlobs[q].end()) {
							dOverlapArea += grid.m_dArea[*iter];
						}
					}
					if (dCurrentArea == 0.0) {
						_EXCEPTIONT("Logic error (zero area blob)");
					}
					if (dOverlapArea == 0.0) {
						_EXCEPTIONT("Logic error (zero overlap area)");
					}
					if (dOverlapArea < dCurrentArea * dMinPercentOverlapNext) {
						continue;
					}
					if (dOverlapArea > dCurrentArea * dMaxPercentOverlapNext) {
						continue;
					}
				}

				// Verify the blobs meet percentage overlap criteria
				// (as a percentage of the earlier blob)
				if ((dMinPercentOverlapPrev != 0.0) ||
				    (dMaxPercentOverlapPrev != 1.0)
				) {
					double dCurrentArea = 0.0;
					double dOverlapArea = 0.0;

					IndicatorSetConstIterator iter = vecPrevBlobs[q].begin();
					for (; iter != vecPrevBlobs[q].end(); iter++) {
						dCurrentArea += grid.m_dArea[*iter];
						if (vecBlobs[p].find(*iter) != vecBlobs[p].end()) {
							dOverlapArea += grid.m_dArea[*iter];
						}
					}
					if (dCurrentArea == 0.0) {
						_EXCEPTIONT("Logic error (zero area blob)");
					}
					if (dOverlapArea == 0.0) {
						_EXCEPTIONT("Logic error (zero overlap area)");
					}
					if (dOverlapArea < dCurrentArea * dMinPercentOverlapPrev) {
						continue;
					}
					if (dOverlapArea > dCurrentArea * dMaxPercentOverlapPrev) {
						continue;
					}
				}

				// Check restrict_region criteria
				if (opRestrictRegion.IsActive()) {
					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					for (; iter != vecBlobs[p].end(); iter++) {
						bool fInRestrictRegion =
							opRestrictRegion.ContainsPoint(
								RadToDeg(grid.m_dLat[*iter]),
								RadToDeg(grid.m_dLon[*iter]));

						if (fInRestrictRegion) {
							setRestrictRegion.insert(vecBlobTags[p]);
							break;
						}
					}
				}


				// Insert bidirectional edge in graph
				multimapTagGraph.insert(
					std::pair<Tag, Tag>(
						vecBlobTags[p], vecPrevBlobTags[q]));

				multimapTagGraph.insert(
					std::pair<Tag, Tag>(
						vecPrevBlobTags[q], vecBlobTags[p]));
			}
		}
	}

	AnnounceEndBlock("Done");
	


	
	#if defined(TEMPEST_MPIOMP)
		if (nMPISize > 1 && valid_flag) {
			//Gather the connectivity graph (multimapTagGraph) to P0
			MapGraphGatherOp MPI_MapGraph(MPI_REAL_COMM,multimapTagGraph);
			MPI_MapGraph.Gather();
			if (nMPIRank == 0){
				multimapTagGraph = MPI_MapGraph.GetGatheredTagGraph();
			}

			//Gather the setAllTags to P0 (Set All Tags can remain at local)
			TagCollectiveOP MPI_TagsGather(MPI_REAL_COMM, MPI_exchangedTags->GetOriginalVecAllBlobTags());
			MPI_TagsGather.Gather();


			//And then reduced the original global time (At this point, the "local" Global Times) to P0 for next step:
			int reducedNGlobalTimes = 0;
			MPI_Reduce(&original_nGlobalTimes, &reducedNGlobalTimes, 1, MPI_INT, MPI_SUM, 0, MPI_REAL_COMM);
			if (nMPIRank == 0) {
				nGlobalTimes = reducedNGlobalTimes;
			}

			MPI_Barrier(MPI_REAL_COMM);
			if (nMPIRank == 0){
				setAllTags = MPI_TagsGather.GetGatheredSetAllTags();
				vecAllBlobTags = MPI_TagsGather.GetUnserialVecAllTags(1);
			}
			
		} 


	#endif

	

	



	AnnounceStartBlock("Identify components in connectivity graph");//Only In processor 0;

	// Total number of blobs
	int nTotalBlobCount = 0;

	// Find all components using a bidirectional graph search
	std::map<Tag, Tag> mapEquivalentTags;

	std::set<Tag>::const_iterator iterTag = setAllTags.begin();

	for (; iterTag != setAllTags.end(); iterTag++) {
		#if defined(TEMPEST_MPIOMP)
			if (nMPIRank != 0) {
				continue;
			}

		#endif



		std::set<Tag> setTagsVisited;

		std::set<Tag> setTagsToVisit;
		setTagsToVisit.insert(*iterTag);

		Tag tagMinimum = *iterTag;

		// Check if this tag is already part of an explored component
		if (mapEquivalentTags.find(*iterTag) != mapEquivalentTags.end()) {
			continue;
		}

		// New blob
		nTotalBlobCount++;

		// All time indices associated with this blob
		std::set<int> setBlobTimes;

		// All time indices where this blob is in the restrict region
		std::set<int> setBlobTimesInRestrictRegion;

		// Find component containing this node
		for (;;) {

			// Out of tags to visit; done
			if (setTagsToVisit.size() == 0) {
				break;
			}

			// Get next tag to visit
			Tag tagNext = *(setTagsToVisit.begin());
			setTagsToVisit.erase(setTagsToVisit.begin());

			// Verify we haven't visited this tag already
			if (setTagsVisited.find(tagNext) != setTagsVisited.end()) {
				continue;
			}
			setTagsVisited.insert(tagNext);

			// Check minimum tag
			if (tagNext < tagMinimum) {
				tagMinimum = tagNext;
			}

			setBlobTimes.insert(tagNext.time);

			if (setRestrictRegion.find(tagNext) != setRestrictRegion.end()) {
				setBlobTimesInRestrictRegion.insert(tagNext.time);
			}

			// Get edges from this node
			std::pair<MapGraphIterator, MapGraphIterator> iterGraphEdges
				= multimapTagGraph.equal_range(tagNext);

			MapGraphIterator iter = iterGraphEdges.first;
			for (; iter != iterGraphEdges.second; iter++) {
				setTagsToVisit.insert(iter->second);
			}
		}

		// Apply filters to the blob
		bool fAcceptBlob = true;

		// Filter on RestrictRegion count for this global_id
		if (opRestrictRegion.IsActive()) {
			if (setBlobTimesInRestrictRegion.size() < opRestrictRegion.GetMinimumCount()) {
				fAcceptBlob = false;
			}
		}

		// Filter on min_time
		if (setBlobTimes.size() < nMinTime) {
			fAcceptBlob = false;
		}
		if (dMinTimeSeconds > 0.0) {
			int iFirstTime = *(setBlobTimes.begin());
			int iLastTime = *(setBlobTimes.rbegin());

			_ASSERT(iFirstTime >= 0);
			_ASSERT(iLastTime < vecGlobalTimeIxToFileTimeIx.size());

			std::pair<int,int> iFileTimeIxFirst = vecGlobalTimeIxToFileTimeIx[iFirstTime];
			std::pair<int,int> iFileTimeIxLast = vecGlobalTimeIxToFileTimeIx[iLastTime];

			const Time & timeFirst = vecGlobalTimes[iFileTimeIxFirst.first][iFileTimeIxFirst.second];
			const Time & timeLast = vecGlobalTimes[iFileTimeIxLast.first][iFileTimeIxLast.second];

			double dBlockDurationSeconds = timeFirst.DeltaSeconds(timeLast);
			if (dBlockDurationSeconds < dMinTimeSeconds) {
				fAcceptBlob = false;
			}
		}

		// Set global id
		if (fAcceptBlob) {
			tagMinimum.global_id = nTotalBlobCount;
		} else {
			tagMinimum.global_id = 0;
			nTotalBlobCount--;
		}

		// Refer all tags in component to minimum tag
		std::set<Tag>::const_iterator iterTagsVisited = setTagsVisited.begin();
		for (; iterTagsVisited != setTagsVisited.end(); iterTagsVisited++) {
			mapEquivalentTags.insert(
				std::pair<Tag,Tag>(*iterTagsVisited, tagMinimum));
		}
	}

	AnnounceEndBlock("Done");

	// Merge blobs at each time step with equivalent tags
	AnnounceStartBlock("Reassign blob tags");
	for (int t = 0; t < nGlobalTimes; t++) {
		#if defined(TEMPEST_MPIOMP)
			if (nMPIRank != 0) {
				break;
			}
		#endif
		_ASSERT(nGlobalTimes == vecAllBlobTags.size());//[HC] This might need to update for the MPI optimization purpose
		std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];//[HC] So we actually use the vecAllBlobTags at P0

		for (int p = 0; p < vecBlobTags.size(); p++) {
			std::map<Tag, Tag>::const_iterator iterTagPair =
				mapEquivalentTags.find(vecBlobTags[p]);

			if (iterTagPair != mapEquivalentTags.end()) {
					vecBlobTags[p] = iterTagPair->second;
			}
		}
	}







	Announce("Unique tags found: %i", nTotalBlobCount);

	AnnounceEndBlock("Done");



/*
	// Apply post-hoc threshold operators
	std::vector<bool> fRejectedBlob;
	fRejectedBlob.resize(nTotalBlobCount+1, false);

	if (vecThresholdOp.size() != 0) {
		AnnounceStartBlock("Applying threshold commands");

		// Whether or not each blob (global_id) satisfies each threshold
		// at each time
		std::vector< std::vector< std::vector<bool> > >
			vecBlobSatisfiesThreshold;

		vecBlobSatisfiesThreshold.resize(nTotalBlobCount+1);
		for (int gid = 1; gid <= nTotalBlobCount; gid++) {
			vecBlobSatisfiesThreshold[gid].resize(
				vecThresholdOp.size());
		}

		// Determine if blobs satisfy the threshold at each time
		for (int t = 0; t < nTime; t++) {

			std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

			std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

			std::vector<LatLonBox> & vecBlobBoxesDeg = vecAllBlobBoxesDeg[t];

			for (int x = 0; x < vecThresholdOp.size(); x++) {

				std::map<int, bool> mapSatisfiesThreshold;

				vecThresholdOp[x].Apply(
					dCellArea,
					dataLatDeg,
					dataLonDeg,
					vecBlobTags,
					vecBlobs,
					vecBlobBoxesDeg,
					mapSatisfiesThreshold);

				std::map<int,bool>::const_iterator iter =
					mapSatisfiesThreshold.begin();

				for (; iter != mapSatisfiesThreshold.end(); iter++) {
					int iGlobalBlobIndex = iter->first;

					if ((iGlobalBlobIndex < 1) ||
					    (iGlobalBlobIndex > nTotalBlobCount)
					) {
						_EXCEPTION1("global_id out of range (%i)",
							iGlobalBlobIndex);
					}
					vecBlobSatisfiesThreshold[iGlobalBlobIndex][x].
						push_back(iter->second);
				}
			}
		}

		// Number of blobs rejected for each reason
		std::vector<int> nBlobsRejected;
		nBlobsRejected.resize(vecThresholdOp.size());

		for (int gid = 1; gid <= nTotalBlobCount; gid++) {
		for (int x = 0; x < vecThresholdOp.size(); x++) {
			int nCount = 0;
			for (int t = 0; t < vecBlobSatisfiesThreshold[gid][x].size(); t++) {
				if (vecBlobSatisfiesThreshold[gid][x][t]) {
					nCount++;
				}
			}

			// Criteria must be satisfied at all times
			if (vecThresholdOp[x].GetMinimumCount() == (-1)) {
				if (nCount != vecBlobSatisfiesThreshold[gid][x].size()) {
					fRejectedBlob[gid] = true;
					nBlobsRejected[x]++;
					break;
				}

			// Criteria must be satisfied at a minimum number of times
			} else if (nCount < vecThresholdOp[x].GetMinimumCount()) {
				fRejectedBlob[gid] = true;
				nBlobsRejected[x]++;
				break;
			}
		}
		}

		// Announce rejections
		for (int x = 0; x < vecThresholdOp.size(); x++) {
			Announce("Rejected (threshold %i): %i", x, nBlobsRejected[x]);
		}

		AnnounceEndBlock("Done");
	}
*/


	///////////////////////////////////////////////////////////////////////////
	// Output results
	///////////////////////////////////////////////////////////////////////////
	AnnounceStartBlock("Output blobs");
	{	
		#if defined(TEMPEST_MPIOMP)
			// Processor 0 scatter the updated vecAllBlobsTag to other processors
			if (nMPISize > 1 && valid_flag) {
				// Restore the original vecGlobalTimes and nGlobalTimes values for correct output,
				// undoing any modifications made during MPI data exchange.
				vecGlobalTimes = MPI_exchangedGlobalTimes->GetUnExchangedVecGlobalTimes();
				nGlobalTimes = original_nGlobalTimes;

				// Create a TagCollectiveOP object for MPI collective operations on vecAllBlobTags.
				// Note: The vecAllBlobTags is only meaningful on processor 0 in this constructor context.
				TagCollectiveOP MPI_TagScatter(MPI_REAL_COMM, vecAllBlobTags);

				// Each processor sends its original vecAllBlobTags size/count information to processor 0.
				// The GatherTagCounts() function gathers this metadata from all processors,
				// updating the local _vecAllBlobTags within the TagCollectiveOP instance for every processor
				// except processor 0, which uses the combined information to build the global vecAllBlobTags.
				MPI_TagScatter.GatherTagCounts(MPI_exchangedTags->GetOriginalVecAllBlobTags());

				// Synchronize all processors to ensure that the gathering process is complete.
				MPI_Barrier(MPI_REAL_COMM);

				// Scatter the combined global vecAllBlobTags data from processor 0 to all processors.
				MPI_TagScatter.Scatter();

				// Synchronize again to ensure the scattering process is complete.
				MPI_Barrier(MPI_REAL_COMM);

				// Retrieve the deserialized vecAllBlobTags on processor 0 (global view) after scattering.
				vecAllBlobTags = MPI_TagScatter.GetUnserialVecAllTags(0);

				// Also restore the original vecAllBlobs data from the exchanged blobs.
				vecAllBlobs = MPI_exchangedBlobs->GetOriginalVecAllBlobs();

			}
		#endif

		// Load in the benchmark file
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[0]);
		_ASSERT(vecNcFiles.size() > 0);

		// Loop through all output files
		_ASSERT(vecOutputFiles.size() == vecGlobalTimes.size());

		int iGlobalTimeIx = 0;

		for (int f = 0; f < vecOutputFiles.size(); f++) {
			#if defined(TEMPEST_MPIOMP) //[Commented out for auto-complete, need to uncomment later]
				if (nMPISize > 1 && valid_flag) {
					
					// Assign each processor with their responsible chunks of files
					if ((f >= processorResponsibalForFile_UB) || f < processorResponsibalForFile_LB) {
						continue;
					}
				}


			#endif //[Commented out for auto-complete, need to uncomment later]

			Announce("Writing file \"%s\"", vecOutputFiles[f].c_str());

			// Open output file
			NcFile ncOutput(vecOutputFiles[f].c_str(), NcFile::Replace);
			if (!ncOutput.is_valid()) {
				_EXCEPTION1("Unable to open output file \"%s\"",
					vecOutputFiles[f].c_str());
			}

			// Output time dimension
			int nLocalTimes = vecGlobalTimes[f].size();

			NcDim * dimOutputTime = ncOutput.add_dim("time", nLocalTimes);
			if (dimOutputTime == NULL) {
				_EXCEPTIONT("Unable to create dimension \"time\" in output");
			}
			NcVar * varOutputTime =
				ncOutput.add_var("time", ncDouble, dimOutputTime);

			DataArray1D<double> dOutputTimes(nLocalTimes);
			for (int t = 0; t < vecGlobalTimes[f].size(); t++) {
				dOutputTimes[t] =
					vecGlobalTimes[f][t].GetCFCompliantUnitsOffsetDouble(strOutTimeUnits);
			}

			varOutputTime->add_att("long_name","time");
			varOutputTime->add_att("units",strOutTimeUnits.c_str());
			varOutputTime->add_att("calendar",vecGlobalTimes[f][0].GetCalendarName().c_str());

			varOutputTime->put(&(dOutputTimes[0]), nLocalTimes);

			// Create output variable
			NcDim * dimOut0 = NULL;
			NcDim * dimOut1 = NULL;
			NcVar * varTagOut = NULL;

			PrepareBlobOutputVar(
				*(vecNcFiles[0]),
				ncOutput,
				vecOutputFiles[f],
				grid,
				strOutputVariable,
				strLatitudeName,
				strLongitudeName,
				ncInt,
				dimOutputTime,
				&dimOut0,
				&dimOut1,
				&varTagOut);

			_ASSERT(varTagOut != NULL);

			int nDimOutSize0 = 0;
			int nDimOutSize1 = 0;

			if (dimOut0 != NULL) {
				nDimOutSize0 = dimOut0->size();
			}
			if (dimOut1 != NULL) {
				nDimOutSize1 = dimOut1->size();
			}

			// Write all time steps
			DataArray1D<int> dataBlobTag(grid.GetSize());

			for (int t = 0; t < nLocalTimes; t++) {

				dataBlobTag.Zero();

				_ASSERT(iGlobalTimeIx + t < vecAllBlobTags.size());
				_ASSERT(iGlobalTimeIx + t < vecAllBlobs.size());
	
				// Get the current blob vectors
				const std::vector<Tag> & vecBlobTags = vecAllBlobTags[iGlobalTimeIx + t];
				const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[iGlobalTimeIx + t];

				_ASSERT(vecBlobTags.size() == vecBlobs.size());

				// Put blob information into dataBlobTag
				for (int p = 0; p < vecBlobTags.size(); p++) {

					if (vecBlobTags[p].global_id == 0) {
						continue;
					}
		
					IndicatorSetConstIterator iter = vecBlobs[p].begin();
					if (!fFlatten) {
						for (; iter != vecBlobs[p].end(); iter++) {
							dataBlobTag[*iter] =
								vecBlobTags[p].global_id;
						}
					} else {
						for (; iter != vecBlobs[p].end(); iter++) {
							dataBlobTag[*iter] = 1.0;
						}
					}
				}

				// Write to file
				if (grid.DimCount() == 1) {
					varTagOut->set_cur(t, 0);
					varTagOut->put(&(dataBlobTag[0]), 1, nDimOutSize0);
				} else {
					varTagOut->set_cur(t, 0);
					varTagOut->put(&(dataBlobTag[0]), 1, nDimOutSize0, nDimOutSize1);
				}
			}

			// Update global time index
			iGlobalTimeIx += nLocalTimes;

			// Close the output file
			ncOutput.close();

		}

		AnnounceEndBlock("Done");
	}
/*
	// Copy variable attributes from first input file
	{
		NcFile ncInput(vecInputFiles[0].c_str());

		NcVar * varLat = ncInput.get_var(strLatitudeName.c_str());
		NcVar * varLon = ncInput.get_var(strLongitudeName.c_str());

		CopyNcVarAttributes(varLat, varOutputLat);
		CopyNcVarAttributes(varLon, varOutputLon);

		NcVar * varTime = ncInput.get_var(strTimeName.c_str());
		if (varTime != NULL) {
			CopyNcVarAttributes(varTime, varOutputTime);
		}
	}
*/
/*
	// Change the time units to the preset
	NcAtt * oldUnits = varOutputTime->get_att("units");
	oldUnits->remove();
	varOutputTime->add_att("units",strOutTimeUnits.c_str());
	NcVar * varData =
		ncOutput.add_var(
			strOutputVariable.c_str(),
			ncInt,
			dimOutputTime,
			dimOutputLat,
			dimOutputLon);

	// Loop through all time steps
	DataArray2D<int> dataBlobTag(nLat, nLon);

	for (int t = 0; t < nTime; t++) {

		dataBlobTag.Zero();

		// Get the current blob vectors
		const std::vector<Tag> & vecBlobTags = vecAllBlobTags[t];

		const std::vector<IndicatorSet> & vecBlobs = vecAllBlobs[t];

		// Put blob information into matrix
		for (int p = 0; p < vecBlobTags.size(); p++) {

			if (vecBlobTags[p].global_id == 0) {
				continue;
			}

			//if (fRejectedBlob[vecBlobTags[p].global_id]) {
			//	continue;
			//}

			IndicatorSetConstIterator iter = vecBlobs[p].begin();
			for (; iter != vecBlobs[p].end(); iter++) {
				dataBlobTag[iter->lat][iter->lon] =
					vecBlobTags[p].global_id;
			}
		}

		// Write to file
		varData->set_cur(t, 0, 0);
		varData->put(&(dataBlobTag[0][0]), 1, nLat, nLon);
	}
*/

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());

#if defined(TEMPEST_MPIOMP)
	MPI_Abort(MPI_COMM_WORLD, 1);
#endif
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////