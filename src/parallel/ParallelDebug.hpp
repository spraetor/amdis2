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



/** \file ParallelDebug.h */

#ifndef AMDIS_PARALLELDEBUG_H
#define AMDIS_PARALLELDEBUG_H

#include "parallel/ParallelTypes.hpp"
#include "parallel/MeshDistributor.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    class ParallelDebug
    {
    protected:
      typedef std::vector<WorldVector<double>> CoordsVec;

      /// Defines a mapping type from rank numbers to sets of coordinates.
      typedef std::map<int, CoordsVec> RankToCoords;

    public:
      /** \brief
      * Tests the interior and the periodic boundaries on all ranks if their order
      * fits together.
      *
      * \param[in]  pdb   Parallel problem definition used for debugging.
      */
      static void testInteriorBoundary(MeshDistributor& pdb);

      /** \brief
      * Test if all periodic boundaries are set in a consistent way on all ranks.
      *
      * \param[in]  pdb   Parallel problem definition used for debugging.
      */
      static void testPeriodicBoundary(MeshDistributor& pdb);

      /** \brief
      * Test if all periodic boundaries are set in a consistent way on all ranks.
      *
      * \param[in]  pdb       Parallel problem definition used for debugging.
      * \oaram[in]  feSpace   FE space for which the DOFs are tested.
      */
      static void testPeriodicBoundary(MeshDistributor& pdb,
                                       const FiniteElemSpace* feSpace);

      /** \brief
      * This function is used for debugging only. It traverses all interior
      * boundaries and compares the DOF indices on them with the DOF indices of
      * the boundarys neighbours. The function fails, when DOF indices on an
      * interior boundary do not fit together.
      *
      * \param[in]  pdb           Parallel problem definition used for debugging.
      * \param[in]  printCoords   If true, the coords of all common dofs are
      *                           printed to the screen.
      */
      static void testCommonDofs(MeshDistributor& pdb, Mesh* mesh, bool printCoords = false);

      /** \brief
      * This function is used for debugging only. It checks if on all ranks DOFs
      * with the same coordinates have the same global index. For this, the
      * function generates on all ranks a list of all DOFs with their coordinates
      * and global indices and  sends the list to all neighbouring ranks.
      *
      * \param[in]  pdb           Parallel problem definition used for debugging.
      */
      static void testGlobalIndexByCoords(MeshDistributor& pdb, Mesh* mesh);

      /** \brief
      * Tests if all elements in the macro mesh are memeber of exactly one rank.
      *
      * \param[in]  pdb           Parallel problem definition used for debugging.
      */
      static void testAllElements(MeshDistributor& pdb);

      /** \brief
      * Tests for all ranks simultaneously, if the number of all send and received
      * DOFs fits together for all communication partners.
      *
      * \param[in]  pdb        Parallel problem definition used for debugging.
      */
      static void testDofContainerCommunication(MeshDistributor& pdb, Mesh* mesh);

      /// Tests if there are multiple DOFs in mesh with the same coords.
      static void testDoubleDofs(Mesh* mesh);

      /** \brief
      * This function is used for debugging only. It prints all information from
      * the local to global dof mapping, see \ref mapLocalGlobalDofs.
      *
      * \param[in]  pdb     Parallel problem definition used for debugging.
      * \param[in]  rank    If specified, only the information from the given rank
      *                     is printed.
      */
      static void printMapLocalGlobal(MeshDistributor& pdb, int rank = -1);

      /** \brief
      * This function is used for debugging only. It prints all information about
      * the periodic mapping of dofs, that are on periodic boundaries.
      *
      * \param[in] pdb    Parallel problem definition used for debugging.
      * \param[in] rank   If specified, only the information from the given rank
      *                   is printed.
      */
      static void printMapPeriodic(MeshDistributor& pdb, int rank = -1);

      /** \brief
      * This function is used for debugging only. It prints information about DOFs
      * in rank's partition.
      *
      * \param[in]  pdb           Parallel problem definition used for debugging.
      * \param[in]  rank          If specified, only the information from the
      *                           given rank is printed.
      * \param[in]  rankDofs      List of all dofs in ranks partition that are
      *                           owned by rank.
      * \param[in]  rankAllDofs   List of all dofs in ranks partition.
      */
      static void printRankDofs(MeshDistributor& pdb,
                                int rank,
                                DofContainer& rankDofs,
                                DofContainer& rankAllDofs);

      /** \brief
      * This functions prints all information about all interior boundaries on
      * all ranks.
      *
      * \param[in]  intBoundary   The boundary object to be printed.
      * \param[in]  force         If true, the information is always printed to
      *                           screen. Otherwise, this is done only if AMDiS
      *                           is compiled in debug mode or if the init file
      *                           parameter "parallel->debug->print boundary info"
      *                           is set.
      */
      static void printBoundaryInfo(InteriorBoundary& intBoundary,
                                    bool force = false);


      static void writeDebugFile(const FiniteElemSpace* feSpace,
                                 ParallelDofMapping& dofMap,
                                 std::string prefix,
                                 std::string postfix);

      static void writeDebugFile(MeshToFeSpaces& meshToFeSpaces,
                                 ParallelDofMapping& dofMap,
                                 std::string debugOutputDir);

      /** \brief
      * This functions create a Paraview file with the macro mesh where the
      * elements are colored by the partition they are part of.
      */
      static void writePartitioning(MeshDistributor& pdb, std::string filename);

      /** \brief
      * The mesh is written to a value and all values are assigned by the rank
      * number where the vertex is contained in.
      *
      * \param[in]  filename    Name of the output file without extension (.vtu).
      * \param[in]  counter     Counter index. If not negative, this number is
      *                         added to the filename.
      * \param[in]  feSpace
      */
      static void writePartitioningFile(std::string filename,
                                        int counter,
                                        const FiniteElemSpace* feSpace);


      static bool followThisBound(int rankElIndex, int neighElIndex);

      static void followBoundary(MeshDistributor& pdb);

      static void followBoundary(Mesh* mesh,
                                 AtomicBoundary& bound,
                                 MeshStructure& code);

      /** \brief
      * Writes Element Map of local Rank
      * Map containes for each DOF in each Element (resulting in massive doubling of DOFs):
      * localElementNumber	elementLevel	localDOFNumber	dofType dofCoords(0-2) elementLocalDOFNumber	typeOfElement (0-2)
      */
      static void writeCsvElementMap(const FiniteElemSpace* feSpace,
                                     ParallelDofMapping& dofMap,
                                     std::string prefix,
                                     std::string postfix);

      static void writeDofMap(ParallelDofMapping& dofMap,
                              int component,
                              std::string filename,
                              std::string postfix);

      static void writePeriodicElObjInfo(MeshDistributor& pdb, std::string debugOutputDir);

      static void writeInterchangeVector(MeshDistributor& pdb, std::string debugOutputDir);
    };
  } // end namespace Parallel
} // end namespace AMDiS

#endif
