/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors:
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 *
 ******************************************************************************/



/** \file StdMpi.h */

#ifndef AMDIS_STDMPI_H
#define AMDIS_STDMPI_H

#include <map>
#include <stdint.h>
#include <mpi.h>
#include "MeshStructure.h"
#include "parallel/InteriorBoundary.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    /** \brief
     * The class StdMpiHelper defines for a type a set of variables, types and
     * functions  that makes it able to transfer objects of the original type to
     * a buffer which is suitable to be send and received using MPI
     * communication.
     *
     *    mpiDataType     Specifies the MPI data type that should be used for
     *                    communication.
     *    cppDataType     Specifies the C++ data type of the buffer.
     *    getBufferSize   This functions calculates the size of the buffer for a
     *                    specific object.
     *    createBuffer    Create the buffer for a specific object.
     *    makeFromBuffer  Create an object of a buffer.
     */
    template<typename T>
    struct StdMpiHelper {};

    template<>
    struct StdMpiHelper<int>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(int& data);

      static void createBuffer(int& data, int* buf);

      static void makeFromBuffer(int& data, int* buf, int bufSize);
    };


    template<>
    struct StdMpiHelper<std::vector<int>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<int>& data);

      static int getBufferSize(std::vector<const int*>& data);

      static void createBuffer(std::vector<int>& data, int* buf);

      static void makeFromBuffer(std::vector<int>& data, int* buf, int bufSize);
    };


    template<>
    struct StdMpiHelper<std::set<int>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::set<int>& data);

      static void createBuffer(std::set<int>& data, int* buf);

      static void makeFromBuffer(std::set<int>& data, int* buf, int bufSize);
    };


    template<>
    struct StdMpiHelper<std::vector<std::set<int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<std::set<int>>& data);

      static void createBuffer(std::vector<std::set<int>>& data, int* buf);

      static void makeFromBuffer(std::vector<std::set<int>>& data,
                                 int* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::set<std::pair<int, int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::set<std::pair<int, int>>& data);

      static void createBuffer(std::set<std::pair<int, int>>& data, int* buf);

      static void makeFromBuffer(std::set<std::pair<int, int>>& data,
                                 int* buf, int bufSize);
    };


    template<>
    struct StdMpiHelper<std::vector<double>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::vector<double>& data);

      static void createBuffer(std::vector<double>& data, double* buf);

      static void makeFromBuffer(std::vector<double>& data, double* buf, int bufSize);
    };


    template<>
    struct StdMpiHelper<std::vector<std::vector<double>>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::vector<std::vector<double>>& data);

      static void createBuffer(std::vector<std::vector<double>>& data, double* buf);

      static void makeFromBuffer(std::vector<std::vector<double>>& data, double* buf, int bufSize);
    };


    // MeshDistributor::MeshCodeVec
    template<>
    struct StdMpiHelper<std::vector<MeshStructure>>
    {
      static MPI_Datatype mpiDataType;

      typedef uint64_t cppDataType;

      static int getBufferSize(std::vector<MeshStructure>& data);

      static void createBuffer(std::vector<MeshStructure>& data, uint64_t* buf);

      static void makeFromBuffer(std::vector<MeshStructure>& data,
                                 uint64_t* buf, int bufSize);
    };


    template<>
    struct StdMpiHelper<std::vector<AtomicBoundary>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<AtomicBoundary>& data);

      static void createBuffer(std::vector<AtomicBoundary>& data, int* buf);

      static void makeFromBuffer(std::vector<AtomicBoundary>& data, int* buf, int bufSize);
    };


    // PeriodicMap::PeriodicDofMap
    template<>
    struct StdMpiHelper<std::map<BoundaryType, std::map<DegreeOfFreedom, DegreeOfFreedom>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::map<BoundaryType, std::map<DegreeOfFreedom, DegreeOfFreedom>>& data);

      static void createBuffer(std::map<BoundaryType, std::map<DegreeOfFreedom, DegreeOfFreedom>>& data, int* buf);

      static void makeFromBuffer(std::map<BoundaryType, std::map<DegreeOfFreedom, DegreeOfFreedom>>& data, int* buf, int bufSize);
    };


    // PetscSolver::createPetscNnzStructure::MatrixNnzEntry
    template<>
    struct StdMpiHelper<std::vector<std::pair<int, int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<std::pair<int, int>>& data);

      static void createBuffer(std::vector<std::pair<int, int>>& data, int* buf);

      static void makeFromBuffer(std::vector<std::pair<int, int>>& data,
                                 int* buf, int bufSize);
    };


    // ParallelDebug::CoordsVec
    template<>
    struct StdMpiHelper<std::vector<WorldVector<double>>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::vector<WorldVector<double>>& data);

      static void createBuffer(std::vector<WorldVector<double>>& data, double* buf);

      static void makeFromBuffer(std::vector<WorldVector<double>>& data, double* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::vector<WorldVector<int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<WorldVector<int>>& data);

      static void createBuffer(std::vector<WorldVector<int>>& data, int* buf);

      static void makeFromBuffer(std::vector<WorldVector<int>>& data, int* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::vector<Vector<int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<Vector<int>>& data);

      static void createBuffer(std::vector<Vector<int>>& data, int* buf);

      static void makeFromBuffer(std::vector<Vector<int>>& data, int* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::vector<WorldMatrix<double>>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::vector<WorldMatrix<double>>& data);

      static void createBuffer(std::vector<WorldMatrix<double>>& data, double* buf);

      static void makeFromBuffer(std::vector<WorldMatrix<double>>& data, double* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::vector<WorldMatrix<int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<WorldMatrix<int>>& data);

      static void createBuffer(std::vector<WorldMatrix<int>>& data, int* buf);

      static void makeFromBuffer(std::vector<WorldMatrix<int>>& data, int* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::vector<WorldVector<WorldVector<double>>>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::vector<WorldVector<WorldVector<double>>>& data);

      static void createBuffer(std::vector<WorldVector<WorldVector<double>>>& data, double* buf);

      static void makeFromBuffer(std::vector<WorldVector<WorldVector<double>>>& data, double* buf, int bufSize);
    };

    template<>
    struct StdMpiHelper<std::vector<Matrix<int>>>
    {
      static MPI_Datatype mpiDataType;

      typedef int cppDataType;

      static int getBufferSize(std::vector<Matrix<int>>& data);

      static void createBuffer(std::vector<Matrix<int>>& data, int* buf);

      static void makeFromBuffer(std::vector<Matrix<int>>& data, int* buf, int bufSize);
    };


    // ParallelDebug::testGlobalIndexByCoords::CoordsIndexMap
    template<>
    struct StdMpiHelper<std::map<WorldVector<double>, int>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::map<WorldVector<double>, int>& data);

      static void createBuffer(std::map<WorldVector<double>, int>& data, double* buf);

      static void makeFromBuffer(std::map<WorldVector<double>, int>& data, double* buf, int bufSize);
    };
    template<>
    struct StdMpiHelper<std::map<int, WorldVector<double>>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      static int getBufferSize(std::map<int, WorldVector<double>>& data);

      static void createBuffer(std::map<int, WorldVector<double>>& data, double* buf);

      static void makeFromBuffer(std::map<int, WorldVector<double>>& data, double* buf, int bufSize);
    };



    template<>
    struct StdMpiHelper<std::vector<std::vector<WorldVector<double>>>>
    {
      static MPI_Datatype mpiDataType;

      typedef double cppDataType;

      typedef std::vector<std::vector<WorldVector<double>>>  DataType;

      static int getBufferSize(DataType& data);

      static void createBuffer(DataType& data, double* buf);

      static void makeFromBuffer(DataType& data, double* buf, int bufSize);
    };


    /** \brief
     * This class is used to easily send and receive STL containers using MPI.
     */
    template<typename SendT, typename RecvT=SendT>
    class StdMpi
    {
    public:
      /** \brief
       * Constructor of the class.
       *
       * \param[in]  comm   Reference to the MPI communicator that should be used.
       * \param[in]  b      If true, the size of the data that will be communicated
       *                    between the ranks will be communicated before to achieve
       *                    correct communication. If the value is false, the data size
       *                    must be set by hand and is assumed to be correct on all
       *                    sending and receiving ranks.
       */
      StdMpi(MPI::Intracomm& comm, bool b = true) :
        mpiComm(comm),
        commPrepared(true),
        exchangeDataSize(b)
      {}

      /// Reset all variables.
      void clear()
      {
        sendData.clear();
        recvData.clear();
        sendDataSize.clear();
        recvDataSize.clear();

        commPrepared = true;
      }


      /** \brief
       * Send data to rank.
       *
       * \param[in]  toRank   Rank number to which the data is send to.
       * \param[in]  data     Data to be send.
       */
      void send(int toRank, SendT& data)
      {
        FUNCNAME("StdMpi::send()");

        TEST_EXIT_DBG(commPrepared)("Communication is not prepared!\n");

        sendData[toRank] = data;
        sendDataSize[toRank] = StdMpiHelper<SendT>::getBufferSize(data);
      }


      /** \brief
       * Send data to rank.
       *
       * \param[in]  data    Maps rank number to the data that should be send to
       *                     these ranks.
       */
      void send(std::map<int, SendT>& data)
      {
        FUNCNAME("StdMpi::send()");

        TEST_EXIT_DBG(commPrepared)("Communication is not prepared!\n");

        for (typename std::map<int, SendT>::iterator it = data.begin();
             it != data.end(); ++it)
        {
          sendData[it->first] = it->second;
          sendDataSize[it->first] = StdMpiHelper<SendT>::getBufferSize(it->second);
        }
      }


      /// Returns sending data, see \ref sendData.
      inline std::map<int, SendT>& getSendData()
      {
        return sendData;
      }


      /// Returns the data that should be send to a specific rank, see \ref sendData.
      inline SendT& getSendData(int rank)
      {
        return sendData[rank];
      }


      /// Updates the buffer sizes for all sending data.
      void updateSendDataSize()
      {
        for (typename std::map<int, SendT>::iterator it = sendData.begin();
             it != sendData.end(); ++it)
          sendDataSize[it->first] = StdMpiHelper<SendT>::getBufferSize(it->second);
      }


      /** \brief
       * Is used to specify from which rank data should be received.
       *
       * \param[in]  fromRank   Rank number.
       * \param[in]  size
       */
      void recv(int fromRank, int size = -1)
      {
        FUNCNAME("StdMpi::recv()");

        TEST_EXIT_DBG(commPrepared)("Communication is not prepared!\n");

        recvDataSize[fromRank] = size;
      }


      /// Is used to specify that data is received from all ranks.
      void recvFromAll(int size = -1)
      {
        FUNCNAME("StdMpi::recvFromAll()");

        TEST_EXIT_DBG(commPrepared)("Communication is no prepared!\n");

        for (int i = 0; i < mpiComm.Get_size(); i++)
          recvDataSize[i] = size;
      }


      /// Is used to specify from which rank data will be received.
      template<class T>
      void recv(std::map<int, T>& fromRanks)
      {
        FUNCNAME("StdMpi::recv()");

        TEST_EXIT_DBG(commPrepared)("Communication is not prepared!\n");

        for (typename std::map<int, T>::iterator it = fromRanks.begin();
             it != fromRanks.end(); ++it)
          recvDataSize[it->first] =
            (exchangeDataSize ? -1 : StdMpiHelper<SendT>::getBufferSize(it->second));
      }


      /// Returns received data, see \ref recvData.
      std::map<int, RecvT>& getRecvData()
      {
        return recvData;
      }


      /// Returns received data from a specific rank, see \ref recvData.
      RecvT& getRecvData(int rank)
      {
        FUNCNAME("StdMpi::getRecvData()");
        TEST_EXIT_DBG(recvData.count(rank))("No recv data from rank %d\n", rank);
        return recvData[rank];
      }


      /// If data sizes should be exchanged before the data itself is communicated,
      /// this will be done in this function.
      void commDataSize()
      {
        FUNCNAME("StdMpi::commDataSize()");

        MPI::Request* request= new MPI::Request[sendData.size() + recvDataSize.size()];
        std::vector<int> sendBuffers;
        sendBuffers.resize(sendDataSize.size());

        int requestCounter = 0;

        for (std::map<int, int>::iterator sendIt = sendDataSize.begin();
             sendIt != sendDataSize.end(); ++sendIt)
        {
          sendBuffers[requestCounter] = sendIt->second;
          request[requestCounter] =
            mpiComm.Isend(&(sendBuffers[requestCounter]), 1,
                          MPI_INT, sendIt->first, 0);
          requestCounter++;
        }

        for (std::map<int, int>::iterator recvIt = recvDataSize.begin();
             recvIt != recvDataSize.end(); ++recvIt)
        {
          request[requestCounter++] =
            mpiComm.Irecv(&(recvIt->second), 1, MPI_INT, recvIt->first, 0);
        }

        MPI::Request::Waitall(requestCounter, request);

#if (DEBUG != 0)
        bool testall = MPI::Request::Testall(requestCounter, request);
        TEST_EXIT(testall)("Should not happen!\n");
#endif

        delete[] request;
      }


      /// Main function of the class, makes the communication. When the function
      /// returns, all data is send and received.
      void startCommunication()
      {
        FUNCNAME("StdMpi::startCommunication()");

        TEST_EXIT_DBG(commPrepared)("Communication is not prepared!\n");

        typedef typename StdMpiHelper<SendT>::cppDataType cppDataType;
        MPI_Datatype mpiDataType = StdMpiHelper<SendT>::mpiDataType;

        if (exchangeDataSize)
          commDataSize();

        // === Remove empty data communication. ===

        {
          std::map<int, int>::iterator it = sendDataSize.begin();
          while (it != sendDataSize.end())
          {
            TEST_EXIT_DBG(it->second >= 0)("Should not happen!\n");

            if (it->second == 0)
            {
              sendData.erase(it->first);
              sendDataSize.erase(it++);
            }
            else
              ++it;
          }
        }

        {
          std::map<int, int>::iterator it = recvDataSize.begin();
          while (it != recvDataSize.end())
          {
            TEST_EXIT_DBG(it->second >= 0)("Should not happen!\n");

            if (it->second == 0)
            {
              recvData.erase(it->first);
              recvDataSize.erase(it++);
            }
            else
              ++it;
          }
        }


        // === Start communication. ===

        MPI::Request* request = new MPI::Request[sendData.size() + recvDataSize.size()];
        int requestCounter = 0;
        std::vector<cppDataType*> sendBuffers, recvBuffers;

        for (typename std::map<int, SendT>::iterator sendIt = sendData.begin();
             sendIt != sendData.end(); ++sendIt)
        {
          int bufferSize = sendDataSize[sendIt->first];

          // Ommit sending empty buffers
          if (bufferSize == 0)
            continue;

          cppDataType* buf = new cppDataType[bufferSize];
          StdMpiHelper<SendT>::createBuffer(sendIt->second, buf);

          request[requestCounter++] =
            mpiComm.Isend(buf, bufferSize, mpiDataType, sendIt->first, 0);

          sendBuffers.push_back(buf);
        }

        for (std::map<int, int>::iterator recvIt = recvDataSize.begin();
             recvIt != recvDataSize.end(); ++recvIt)
        {
          // Ommit receiving empty buffers
          if (recvIt->second == 0)
            continue;

          cppDataType* buf = new cppDataType[recvIt->second];

          request[requestCounter++] =
            mpiComm.Irecv(buf, recvIt->second, mpiDataType, recvIt->first, 0);

          recvBuffers.push_back(buf);
        }

        MPI::Request::Waitall(requestCounter, request);

        for (unsigned int i = 0; i < sendBuffers.size(); i++)
          delete [] sendBuffers[i];
        sendBuffers.clear();

        int i = 0;
        for (std::map<int, int>::iterator recvIt = recvDataSize.begin();
             recvIt != recvDataSize.end(); ++recvIt)
        {
          if (recvIt->second == 0)
            continue;

          StdMpiHelper<SendT>::makeFromBuffer(recvData[recvIt->first],
                                              recvBuffers[i],
                                              recvIt->second);
          delete [] recvBuffers[i];
          i++;
        }
        delete[] request;
        commPrepared = false;
      }

    protected:
      ///
      MPI::Intracomm mpiComm;

      ///
      std::map<int, SendT> sendData;

      ///
      std::map<int, RecvT> recvData;

      std::map<int, int> sendDataSize;

      std::map<int, int> recvDataSize;

      bool commPrepared;

      bool exchangeDataSize;

      static int ccc;
    };

  }
}

#endif
