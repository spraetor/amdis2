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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <stdint.h>
#include <boost/filesystem.hpp>

#include "parallel/MeshDistributor.h"
#include "parallel/MeshManipulation.h"
#include "parallel/ParallelDebug.h"
#include "parallel/StdMpi.h"
#include "parallel/MeshPartitioner.h"
#include "parallel/ParMetisPartitioner.h"
#include "parallel/ZoltanPartitioner.h"
#include "parallel/SimplePartitioner.h"
#include "parallel/CheckerPartitioner.h"
#include "parallel/MpiHelper.h"
#include "parallel/DofComm.h"
#include "parallel/ParallelProblemStat.h"
#include "io/ElementFileWriter.h"
#include "io/MacroInfo.h"
#include "io/MacroWriter.h"
#include "io/VtkWriter.h"
#include "io/ArhReader.h"
#include "io/Arh2Reader.h"
#include "Mesh.h"
#include "Traverse.h"
#include "ElInfo.h"
#include "Element.h"
#include "MacroElement.h"
#include "DOFMatrix.h"
#include "DOFVector.h"
#include "SystemVector.h"
#include "ElementDofIterator.h"
#include "ProblemStatBase.h"
#include "StandardProblemIteration.h"
#include "VertexVector.h"
#include "MeshStructure.h"
#include "ProblemStat.h"
#include "ProblemInstat.h"
#include "RefinementManager3d.h"
#include "Debug.h"
#include "Timer.h"
#include "io/MacroReader.h"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace boost::filesystem;
    using namespace std;

    MeshDistributor* MeshDistributor::globalMeshDistributor = NULL;

    const Flag MeshDistributor::BOUNDARY_SUBOBJ_SORTED              = 0X01L;
    const Flag MeshDistributor::BOUNDARY_FILL_INFO_SEND_DOFS        = 0X02L;
    const Flag MeshDistributor::BOUNDARY_FILL_INFO_RECV_DOFS        = 0X04L;

    inline bool cmpDofsByValue(const DegreeOfFreedom* dof1, const DegreeOfFreedom* dof2)
    {
      return (*dof1 < *dof2);
    }

    MeshDistributor::MeshDistributor()
      : problemStat(0),
        initialized(false),
        name("parallel"),
        macroMesh(NULL),
        refineManager(NULL),
        partitioner(NULL),
        initialPartitioner(NULL),
        deserialized(false),
        writeSerializationFile(false),
        repartitioningAllowed(false),
        repartitionOnlyOnce(false),
        repartitionIthChange(20),
        repartitioningWaitAfterFail(20),
        nMeshChangesAfterLastRepartitioning(0),
        repartitioningCounter(0),
        repartitioningFailed(0),
        debugOutputDir(""),
        createBoundaryDofFlag(0),
        boundaryDofInfo(1),
        meshAdaptivity(true),
        hasPeriodicBoundary(false),
        printTimings(false),
        printMemoryUsage(false)
    {
      FUNCNAME("MeshDistributor::MeshDistributor()");

      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;
      mpiRank = mpiComm.Get_rank();

      Parameters::get(name + "->repartitioning", repartitioningAllowed);
      Parameters::get(name + "->debug output dir", debugOutputDir);
      Parameters::get(name + "->repartition only once", repartitionOnlyOnce);
      Parameters::get(name + "->repartition ith change", repartitionIthChange);
      Parameters::get(name + "->repartition wait after fail", repartitioningWaitAfterFail);
      Parameters::get(name + "->mesh adaptivity", meshAdaptivity);

      nMeshChangesAfterLastRepartitioning = repartitionIthChange - 1;

      // === Create partitioner object. ===

      string partStr = "parmetis";
      Parameters::get(name + "->partitioner", partStr);

      if (partStr == "parmetis")
        partitioner = new ParMetisPartitioner("parallel->partitioner", &mpiComm);

      if (partStr == "zoltan")
      {
#ifdef HAVE_ZOLTAN
        partitioner = new ZoltanPartitioner("parallel->partitioner", &mpiComm);
#else
        ERROR_EXIT("AMDiS was compiled without Zoltan support. Therefore you cannot make use of it!\n");
#endif
      }

      if (partStr == "checker")
        partitioner = new CheckerPartitioner("parallel->partitioner", &mpiComm);

      if (partStr == "simple")
        partitioner = new SimplePartitioner("parallel->partitioner", &mpiComm);

      if (!partitioner)
      {
        ERROR_EXIT("Unknown partitioner or no partitioner specified!\n");
      }

      // === Create initial partitioner object. ===

      partStr = "";
      Parameters::get(name + "->initial partitioner", partStr);
      if (partStr == "")
      {
        initialPartitioner = partitioner;
      }
      else
      {
        if (partStr == "checker")
        {
          initialPartitioner =
            new CheckerPartitioner("parallel->initial partitioner", &mpiComm);
        }
        else
        {
          ERROR_EXIT("Not yet supported, but very easy to implement!\n");
        }
      }


      // === And read some more parameters. ===

      int tmp = 0;
      Parameters::get(name + "->box partitioning", tmp);
      partitioner->setBoxPartitioning(static_cast<bool>(tmp));
      initialPartitioner->setBoxPartitioning(static_cast<bool>(tmp));

      Parameters::get(name + "->print timings", printTimings);
      Parameters::get(name + "->print memory usage", printMemoryUsage);

      // If required, create hierarchical mesh level structure.
      createMeshLevelStructure();
    }


    MeshDistributor::~MeshDistributor()
    {
      if (partitioner)
      {
        delete partitioner;
        partitioner = NULL;
      }
    }

    void MeshDistributor::addInterchangeVector(SystemVector* vec)
    {
      for (int i = 0; i < vec->getSize(); i++)
        interchangeVectors.push_back(vec->getDOFVector(i));
    }

    void MeshDistributor::removeInterchangeVector(SystemVector* vec)
    {
      for (int i = 0; i < vec->getSize(); i++)
        removeInterchangeVector(vec->getDOFVector(i));
    }

    void MeshDistributor::initParallelization()
    {
      FUNCNAME("MeshDistributor::initParallelization()");

      if (initialized)
        return;

      double first = MPI::Wtime();
      MSG("Initialization phase 1 needed %.5f seconds\n",
          first - ParallelProblemStat::initTimeStamp);

      TEST_EXIT(MPI::COMM_WORLD.Get_size() > 1)
      ("Parallelization does not work with only one process!\n");
      TEST_EXIT(feSpaces.size() > 0)
      ("No FE space has been defined for the mesh distributor!\n");
      TEST_EXIT(meshes.size() > 0)("No mesh has been defined for the mesh distributor!\n");


      // === Sort FE spaces with respect to the degree of the basis ===
      // === functions. Use stuiped bubble sort for this.           ===

      //TODO Maybe unnecessary anymore
      bool doNext = false;
      do
      {
        doNext = false;
        for (unsigned int i = 0; i < feSpaces.size() - 1; i++)
        {
          if (feSpaces[i]->getBasisFcts()->getDegree() >
              feSpaces[i + 1]->getBasisFcts()->getDegree())
          {
            const FiniteElemSpace* tmp = feSpaces[i + 1];
            feSpaces[i + 1] = feSpaces[i];
            feSpaces[i] = tmp;
            doNext = true;
          }
        }
      }
      while (doNext);

      // Sort FE spaces seperately for each mesh.
      map<Mesh*, vector<const FiniteElemSpace*>>::iterator iter = meshToFeSpaces.begin();
      while(iter != meshToFeSpaces.end())
      {
        vector<const FiniteElemSpace*>& spaces = iter->second;

        bool doNext = false;
        do
        {
          doNext = false;
          for (size_t i = 0; i < spaces.size() - 1; i++)
          {
            if (spaces[i]->getBasisFcts()->getDegree() >
                spaces[i + 1]->getBasisFcts()->getDegree())
            {
              const FiniteElemSpace* tmp = spaces[i + 1];
              spaces[i + 1] = spaces[i];
              spaces[i] = tmp;
              doNext = true;
            }
          }
        }
        while (doNext);
        ++iter;
      }

      // Test, if the mesh is the macro mesh only! Paritioning of the mesh is
      // supported only for macro meshes, so it will not work yet if the mesh is
      // already refined in some way.
      testForMacroMesh();

#if (DEBUG != 0)
      // Check whether meshes come from the same macro mesh. The way is to compare
      // the node coords of each macro element in the meshes.
      debug::ElementIdxToCoords macroCoords;
      debug::createNodeCoords(macroMesh, macroCoords);

      for (size_t i = 1; i < meshes.size(); i++)
        debug::testNodeCoords(meshes[i], macroCoords);

#endif

      // Initialize dof communicators which have been created in addProblemStat.
      for (size_t i = 0; i < meshes.size(); i++)
        dofComms[meshes[i]].init(levelData, meshToFeSpaces[meshes[i]]);


      // Initialize element object DB
      elObjDb.init(meshes, meshToFeSpaces);

      // If the problem has been already read from a file, we need only to set
      // isRankDofs to all matrices and rhs vector and to remove periodic
      // boundary conditions (if there are some).
      if (deserialized)
      {
        updateMacroElementInfo();

        removePeriodicBoundaryConditions();

        elObjDb.createMacroElementInfo(allMacroElements);

        updateDofRelatedStruct();

        elObjDb.setData(partitionMap, levelData);

#if (DEBUG != 0)
        TEST_EXIT_DBG(dofMaps.size())("No DOF mapping defined!\n");
        ParallelDebug::writeDebugFile(feSpaces[feSpaces.size() - 1],
                                      *(dofMaps[0]),
                                      debugOutputDir + "mpi-dbg", "dat");
#endif

        initialized = true;
        return;
      }


      // Create a very first pariitioning of the mesh.
      createInitialPartitioning();


#if (DEBUG != 0)
      std::vector<debug::ElementIdxToDofs> elMap(meshes.size());
      for (size_t i = 0; i < meshes.size(); i++)
      {
        debug::createSortedDofs(meshes[i], elMap[i]);
      }
#endif

      if (mpiRank == 0)
      {
#if (DEBUG != 0)
        int writePartMesh = 1;
#else
        int writePartMesh = 0;
#endif
        Parameters::get("parallel->debug->write mesh partitioning", writePartMesh);

        if (writePartMesh > 0)
        {
          debug::writeElementIndexMesh(macroMesh , debugOutputDir + "elementIndex");
          ParallelDebug::writePartitioning(*this, debugOutputDir + "part");
        }
      }

      // Create interior boundary information.
      createInteriorBoundary(true);

      // === Remove neighbourhood relations due to periodic bounday conditions. ===

      for (size_t i = 0; i < meshes.size(); i++)
      {
        for (deque<MacroElement*>::iterator it = meshes[i]->firstMacroElement();
             it != meshes[i]->endOfMacroElements(); ++it)
        {
          for (int j = 0; j < macroMesh->getGeo(NEIGH); j++)
          {
            if ((*it)->getNeighbour(j) &&
                meshes[i]->isPeriodicAssociation((*it)->getBoundary(j)))
            {

              int neighIndex = (*it)->getNeighbour(j)->getIndex();

              (*it)->getNeighbour(j)->setNeighbour((*it)->getOppVertex(j), NULL);
              (*it)->setNeighbour(j, NULL);
              (*it)->setBoundary(j, 0);

              if (i == 0)
              {
                macroElementNeighbours[(*it)->getIndex()][j] = -1;
                macroElementNeighbours[neighIndex][(*it)->getOppVertex(j)] = -1;
              }
            }
          }
        }
      }

      for (map<Mesh*, std::vector<MacroElement*>>::iterator iter = allMacroElements.begin();
           iter != allMacroElements.end(); iter++)
      {
        Mesh* mesh = iter->first;
        vector<MacroElement*>& allMacros = iter->second;

        for (vector<MacroElement*>::iterator it = allMacros.begin();
             it != allMacros.end(); ++it)
        {
          for (int i = 0; i < macroMesh->getGeo(NEIGH); i++)
          {
            if ((*it)->getNeighbour(i) &&
                mesh->isPeriodicAssociation((*it)->getBoundary(i)))
            {

              int neighIndex = (*it)->getNeighbour(i)->getIndex();

              (*it)->getNeighbour(i)->setNeighbour((*it)->getOppVertex(i), NULL);
              (*it)->setNeighbour(i, NULL);
              (*it)->setBoundary(i, 0);

              if (i == 0)
              {
                macroElementNeighbours[(*it)->getIndex()][i] = -1;
                macroElementNeighbours[neighIndex][(*it)->getOppVertex(i)] = -1;
              }
            }
          }
        }
      }

      // === Remove all macro elements that are not part of the rank partition. ===

      removeMacroElements();

      // === Create new global and local DOF numbering. ===

      // We have to remove the VertexVectors, which contain periodic assoiciations,
      // because they are not valid anymore after some macro elements have been
      // removed and the corresponding DOFs were deleted.
      for (size_t i = 0; i < meshes.size(); i++)
        for (map<BoundaryType, VertexVector*>::iterator it = meshes[i]->getPeriodicAssociations().begin();
             it != meshes[i]->getPeriodicAssociations().end(); ++it)
          const_cast<DOFAdmin&>(meshes[i]->getDofAdmin(0)).removeDOFContainer(dynamic_cast<DOFContainer*>(it->second));



      updateDofRelatedStruct();

      // === In 3D we have to fix the mesh to allow local refinements. ===

      fix3dMeshRefinement();

      // === If in debug mode, make some tests. ===

#if (DEBUG != 0)
      MSG("AMDiS runs in debug mode, so make some test ...\n");

      ParallelDebug::testAllElements(*this);
      for (size_t i = 0; i < meshes.size(); i++)
        debug::testSortedDofs(meshes[i], elMap[i]);
      ParallelDebug::testInteriorBoundary(*this);
      ParallelDebug::followBoundary(*this);

      debug::writeMesh(feSpaces[0], -1, debugOutputDir + "macro_mesh");

      MSG("Debug mode tests finished!\n");
#endif

      // Remove periodic boundary conditions in sequential problem definition.
      removePeriodicBoundaryConditions();

#if (DEBUG != 0)
      ParallelDebug::testPeriodicBoundary(*this);
#endif

      // === Global refinements. ===
      bool oneMeshRefined = false;
      for (size_t i = 0; i < meshes.size(); i++)
      {

        int globalRefinement = 0;
        Parameters::get(meshes[i]->getName() + "->global refinements", globalRefinement);

        if (globalRefinement > 0)
        {
          oneMeshRefined = true;
          bool doRefineInter = true;
          Parameters::get(meshes[i]->getName() + "->refinement interpol", doRefineInter);
          std::map<DOFVector<double>*, RefineCoarsenOperation> rememberOp;
          if (!doRefineInter)
          {
            // No refinement during initial global refinement
            for (int iadmin = 0; iadmin < meshes[i]->getNumberOfDOFAdmin(); iadmin++)
            {
              std::list<DOFIndexedBase*>::iterator it;
              DOFAdmin* admin = const_cast<DOFAdmin*>(&meshes[i]->getDofAdmin(iadmin));
              std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
              for (it = admin->beginDOFIndexed(); it != end; it++)
              {
                DOFVector<double>* vec = dynamic_cast<DOFVector<double>*>(*it);
                if (vec)
                {
                  rememberOp[vec] = vec->getRefineOperation();
                  vec->setRefineOperation(NO_OPERATION);
                }
              }
            }
          }
          refineManager->globalRefine(meshes[i], globalRefinement);
          if (!doRefineInter)
          {
            // No refinement during initial global refinement
            for (int iadmin = 0; iadmin < meshes[i]->getNumberOfDOFAdmin(); iadmin++)
            {
              std::list<DOFIndexedBase*>::iterator it;
              DOFAdmin* admin = const_cast<DOFAdmin*>(&meshes[i]->getDofAdmin(iadmin));
              std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
              for (it = admin->beginDOFIndexed(); it != end; it++)
              {
                DOFVector<double>* vec = dynamic_cast<DOFVector<double>*>(*it);
                if (vec)
                  vec->setRefineOperation(rememberOp[vec]);
              }
            }
          }
          updateDofRelatedStruct(meshes[i]);

#if (DEBUG != 0)
          ParallelDebug::testPeriodicBoundary(*this);
#endif
        }
      }
      if (oneMeshRefined)
        updateLocalGlobalNumbering();

      // And delete some data, we there is no mesh adaptivty.
      if (!meshAdaptivity)
        elObjDb.clear();

      initialized = true;
      MSG("Init parallelization needed %.5f seconds\n", MPI::Wtime() - first);
    }


    void MeshDistributor::createInitialPartitioning()
    {
      FUNCNAME("MeshDistributor::createInitialPartitioning()");

      // For later mesh repartitioning, we need to store some information about
      // the macro mesh.
      createMacroElementInfo();

      // Create an initial partitioning of the mesh
      bool useInitialPartitioning =
        initialPartitioner->createInitialPartitioning();

      // Set the element weights, which are 1 at the very first begin
      setInitialElementWeights();

      if (!useInitialPartitioning)
      {
        // And now partition the mesh
        bool partitioningSucceed =
          initialPartitioner->partition(elemWeights, INITIAL);
        TEST_EXIT(partitioningSucceed)("Initial partitioning does not work!\n");
      }

      initialPartitioner->createPartitionMap(partitionMap);

      if (initialPartitioner != partitioner)
      {
        *partitioner = *initialPartitioner;
      }
    }


    void MeshDistributor::setInitialElementWeights()
    {
      FUNCNAME("MeshDistributor::setInitialElementWeights()");

      elemWeights.clear();

      string filename = "";
      Parameters::get(macroMesh->getName() + "->macro weights", filename);
      if (filename != "")
      {
        MSG("Read macro weights from %s\n", filename.c_str());

        ifstream infile;
        infile.open(filename.c_str(), ifstream::in);
        while (!infile.eof())
        {
          int elNum, elWeight;
          infile >> elNum;
          if (infile.eof())
            break;
          infile >> elWeight;

          elemWeights[elNum] = elWeight;
        }

        infile.close();

        return;
      }

      filename = "";
      Parameters::get("parallel->partitioner->read meta arh", filename);
      if (filename != "")
      {
        map<int, int> arhElInRank;
        map<int, int> arhElCodeSize;

        /*int nProc = */ io::Arh2Reader::readMetaData(filename, arhElInRank, arhElCodeSize);
        for (map<int, int>::iterator it = arhElCodeSize.begin();
             it != arhElCodeSize.end(); ++it)
          elemWeights[it->first] = it->second;

        return;
      }

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(macroMesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        elemWeights[elInfo->getElement()->getIndex()] = 1.0;
        elInfo = stack.traverseNext(elInfo);
      }
    }


    void MeshDistributor::addProblemStat(ProblemStatSeq* probStat)
    {
      FUNCNAME("MeshDistributor::addProblemStat()");

      TEST_EXIT_DBG(probStat->getFeSpaces().size())
      ("No FE spaces in stationary problem!\n");


      // === Add all FE spaces from stationary problem. ===

      vector<const FiniteElemSpace*> newFeSpaces = probStat->getFeSpaces();

      for (size_t i = 0; i < newFeSpaces.size(); i++)
        if (find(feSpaces.begin(), feSpaces.end(), newFeSpaces[i]) ==
            feSpaces.end())
        {
          // Find new fespace, then add mesh if it's new
          Mesh* mesh = newFeSpaces[i]->getMesh();
          if (find(meshes.begin(), meshes.end(), mesh) == meshes.end())
          {
            TEST_EXIT(meshes.size() < 2)
            ("Currently max two meshes supported on parallel mode.\n");

            meshes.push_back(mesh);
            dofComms[mesh] = MultiLevelDofComm();
            lastMeshChangeIndexs[mesh] = 0;
          }

          feSpaces.push_back(newFeSpaces[i]);
          meshToFeSpaces[mesh].push_back(newFeSpaces[i]);
        }

      // === Set the first mesh to be the macro mesh and pass it to partitioner.    ===
      // === Create a corresponding refinement manager object.                      ===

      if (problemStat.empty())
      {
        macroMesh = meshes[0];

        switch (macroMesh->getDim())
        {
        case 2:
          refineManager = new RefinementManager2d();
          break;
        case 3:
          refineManager = new RefinementManager3d();
          break;
        default:
          ERROR_EXIT("This should not happen for dim = %d!\n", macroMesh->getDim());
        }

        partitioner->setMesh(macroMesh);
        initialPartitioner->setMesh(macroMesh);
      }


      // === Check whether the stationary problem should be serialized. ===


      // Create parallel serialization file writer, if needed.
      int writeSerialization = 0;
      Parameters::get(probStat->getName() + "->output->write serialization",
                      writeSerialization);
      if (writeSerialization && !writeSerializationFile)
      {
        string filename = "";
        Parameters::get(name + "->output->serialization filename", filename);

        TEST_EXIT(filename != "")
        ("No filename defined for parallel serialization file!\n");

        int tsModulo = -1;
        Parameters::get(probStat->getName() + "->output->write every i-th timestep",
                        tsModulo);

        probStat->getFileWriterList().push_back(new Serializer<MeshDistributor>(this, filename, tsModulo));
        writeSerializationFile = true;
      }


      // === Check whether the stationary problem should be deserialized. ===

      int readSerialization = 0;
      Parameters::get(probStat->getName() + "->input->read serialization",
                      readSerialization);
      if (readSerialization)
      {
        string filename = "";
        Parameters::get(probStat->getName() + "->input->serialization filename",
                        filename);
        filename += ".p" + std::to_string(mpiRank);
        MSG("Start deserialization with %s\n", filename.c_str());
        ifstream in(filename.c_str());

        TEST_EXIT(!in.fail())("Could not open deserialization file: %s\n",
                              filename.c_str());

        // Read the revision number of the AMDiS version which was used to create
        // the serialization file.
        int revNumber = -1;
        SerUtil::deserialize(in, revNumber);

        probStat->deserialize(in);
        in.close();
        MSG("Deserialization from file: %s\n", filename.c_str());

        filename = "";
        Parameters::get(name + "->input->serialization filename", filename);

        TEST_EXIT(filename != "")
        ("No filename defined for parallel deserialization file!\n");

        string rankFilename = filename + ".p" + std::to_string(mpiRank);
        in.open(rankFilename.c_str());

        TEST_EXIT(!in.fail())("Could not open parallel deserialization file: %s\n",
                              filename.c_str());

        // Read the revision number of the AMDiS version which was used to create
        // the serialization file.
        revNumber = -1;
        SerUtil::deserialize(in, revNumber);

        deserialize(in);
        in.close();
        MSG("Deserializtion of mesh distributor from file: %s\n", rankFilename.c_str());
        deserialized = true;
      }

      problemStat.push_back(probStat);


      // === If the mesh distributor is already initialized, don't forget to  ===
      // === remove the periodic boundary conditions on these objects.        ===

      if (initialized)
        removePeriodicBoundaryConditions(probStat);
    }


    void MeshDistributor::addProblemStatGlobal(ProblemStatSeq* probStat)
    {
      if (globalMeshDistributor == NULL)
        globalMeshDistributor = new MeshDistributor();

      globalMeshDistributor->addProblemStat(probStat);
    }


    void MeshDistributor::exitParallelization()
    {}


    void MeshDistributor::registerDofMap(ParallelDofMapping& dofMap)
    {
      FUNCNAME("MeshDistributor::registerDofMap()");

      TEST_EXIT(find(dofMaps.begin(), dofMaps.end(), &dofMap) ==
                dofMaps.end())
      ("Parallel DOF mapping already registerd in mesh distributor object!\n");

      dofMaps.push_back(&dofMap);
    }


    void MeshDistributor::removeDofMap(ParallelDofMapping& dofMap)
    {
      vector<ParallelDofMapping*>::iterator it =
        find(dofMaps.begin(), dofMaps.end(), &dofMap);

      if (it != dofMaps.end())
        dofMaps.erase(it);
    }


    void MeshDistributor::testForMacroMesh()
    {
      FUNCNAME("MeshDistributor::testForMacroMesh()");

      for (size_t i = 0; i < meshes.size(); i++)
      {

        int nMacroElements = 0;

        TraverseStack stack;
        ElInfo* elInfo = stack.traverseFirst(meshes[i], -1, Mesh::CALL_LEAF_EL);
        while (elInfo)
        {
          TEST_EXIT(elInfo->getLevel() == 0)
          ("Mesh is already refined! This does not work with parallelization!\n");

          TEST_EXIT(elInfo->getType() == 0)
          ("Only macro elements with level 0 are supported!\n");

          nMacroElements++;

          elInfo = stack.traverseNext(elInfo);
        }

        TEST_EXIT(nMacroElements >= MPI::COMM_WORLD.Get_size())
        ("The mesh has less macro elements than number of mpi processes!\n");
      }
    }


    void MeshDistributor::synchVector(SystemVector& vec)
    {
      FUNCNAME_DBG("MeshDistributor::synchVector()");

      int nLevels = levelData.getNumberOfLevels();
      for (int level = nLevels - 1; level >= 0; level--)
      {
        MPI::Intracomm& mpiComm = levelData.getMpiComm(level);

        if (mpiComm == MPI::COMM_SELF)
          continue;

        vector<Mesh*> sysMeshes;

        // std::set of mesh pointer cannot be used for parallel
        for (size_t i = 0; i < meshes.size(); i++)
          for (int j = 0; j < vec.getSize(); j++)
            if (meshes[i] == vec.getFeSpace(j)->getMesh())
            {
              sysMeshes.push_back(meshes[i]);
              break;
            }

        TEST_EXIT_DBG(sysMeshes.size() <= 2)
        ("This really should not happen.\n");

        vector<Mesh*>::iterator meshIt = sysMeshes.begin();
        while(meshIt != sysMeshes.end())
        {
          StdMpi<vector<double>> stdMpi(mpiComm);

          MultiLevelDofComm& dofComm = dofComms[*meshIt];

          for (DofComm::Iterator it(dofComm[level].getSendDofs());
               !it.end(); it.nextRank())
          {
            vector<double> dofs;

            for (int j = 0; j < vec.getSize(); j++)
            {
              if (vec.getFeSpace(j)->getMesh() != *meshIt)
                continue;

              DOFVector<double>& dofVec = *(vec.getDOFVector(j));
              for (it.beginDofIter(vec.getFeSpace(j)); !it.endDofIter(); it.nextDof())
                dofs.push_back(dofVec[it.getDofIndex()]);
            }

            int rank = it.getRank();
            // if (level > 0)
            //   rank = levelData.mapRank(rank, 0, level);
            stdMpi.send(rank, dofs);
          }

          for (DofComm::Iterator it(dofComm[level].getRecvDofs());
               !it.end(); it.nextRank())
          {
            int rank = it.getRank();
            // if (level > 0)
            //   rank = levelData.mapRank(rank, 0, level);
            stdMpi.recv(rank);
          }

          stdMpi.startCommunication();

          for (DofComm::Iterator it(dofComm[level].getRecvDofs());
               !it.end(); it.nextRank())
          {
            int rank = it.getRank();
            // if (level > 0)
            //   rank = levelData.mapRank(rank, 0, level);

            size_t counter = 0;
            vector<double>& recvData =  stdMpi.getRecvData(rank);

            for (int j = 0; j < vec.getSize(); j++)
            {
              if (vec.getFeSpace(j)->getMesh() != *meshIt)
                continue;

              DOFVector<double>& dofVec = *(vec.getDOFVector(j));

              for (it.beginDofIter(vec.getFeSpace(j)); !it.endDofIter(); it.nextDof())
              {
                TEST_EXIT_DBG(counter < recvData.size())
                ("Recv data from rank %d has only %d entries!\n", rank, recvData.size());
                dofVec[it.getDofIndex()] = recvData[counter++];
              }
            }
          }
          meshIt++;
        }	// End while
      }
    }


    void MeshDistributor::fix3dMeshRefinement()
    {
      FUNCNAME("MeshDistributor::fix3dMeshRefinement()");

      if (macroMesh->getDim() != 3)
        return;

      Timer t;

      // === Create set of all edges and all macro element indices in rank's ===
      // === subdomain.                                                      ===

      std::set<DofEdge> allEdges;
      std::set<int> rankMacroEls;

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(macroMesh, 0, Mesh::CALL_EL_LEVEL);
      while (elInfo)
      {
        for (int i = 0; i < macroMesh->getGeo(EDGE); i++)
        {
          ElementObjectData elData(elInfo->getElement()->getIndex(), i);
          DofEdge e = elObjDb.getEdgeLocalMap(elData);
          allEdges.insert(e);
        }

        rankMacroEls.insert(elInfo->getElement()->getIndex());

        elInfo = stack.traverseNext(elInfo);
      }


      // === Reset fixed edges, and start new search for edges that must  ===
      // === be fixed.                                                    ===

      FixRefinementPatch::connectedEdges.clear();

      // Traverse all edges
      for (std::set<DofEdge>::iterator it = allEdges.begin(); it != allEdges.end(); ++it)
      {
        // For each edge get all macro elements that contain this edge.
        vector<ElementObjectData>& edgeEls = elObjDb.getElements(*it);

        TEST_EXIT_DBG(edgeEls.size() > 0)
        ("No edge %d/%d in elObjDB!\n", it->first, it->second);


        // We have now to check, if the current edge connects two macro elements,
        // which are otherwise not connected.

        // All elements on this edge.
        std::set<int> allEls;
        // Maps from each element index to the local edge index of the common edge.
        map<int, int> edgeNoInEl;

        for (unsigned int i = 0; i < edgeEls.size(); i++)
        {
          if (rankMacroEls.count(edgeEls[i].elIndex))
          {
            allEls.insert(edgeEls[i].elIndex);
            edgeNoInEl[edgeEls[i].elIndex] = edgeEls[i].ithObject;
          }
        }

        TEST_EXIT_DBG(!allEls.empty())("Should not happen!\n");


        // If there is only one element in rank on this edge, there is nothing
        // to do.
        if (allEls.size() == 1)
          continue;

        // Create a set of element index sets. All element within one set are
        // neighbours of each other, but can not reach any other element of the
        // other index sets by neighbourhood relations.
        vector<std::set<int>> disconnectedEls;
        helpToFix(allEls, disconnectedEls);

        // Fix the edges
        for (int dc0 = 0; dc0 < static_cast<int>(disconnectedEls.size()); dc0++)
        {
          for (int dc1 = 0; dc1 < static_cast<int>(disconnectedEls.size()); dc1++)
          {
            if (dc0 == dc1)
              continue;

            for (size_t i = 0; i < meshes.size(); i++)
            {

              int elIdx = *(disconnectedEls[dc1].begin());

              pair<Element*, int> edge0 =
                make_pair(elObjDb.getElementPtr(elIdx, meshes[i]), edgeNoInEl[elIdx]);

              for (std::set<int>::iterator elIt = disconnectedEls[dc0].begin();
                   elIt != disconnectedEls[dc0].end(); ++elIt)
              {
                pair<Element*, int> edge1 =
                  make_pair(elObjDb.getElementPtr(*elIt, meshes[i]), edgeNoInEl[*elIt]);

#if (DEBUG != 0)
                DofEdge dofEdge0 = edge0.first->getEdge(edge0.second);
                DofEdge dofEdge1 = edge1.first->getEdge(edge1.second);

                WorldVector<double> c0, c1;
                meshes[i]->getDofIndexCoords(dofEdge0.first, meshToFeSpaces[meshes[i]][0], c0);
                meshes[i]->getDofIndexCoords(dofEdge0.second, meshToFeSpaces[meshes[i]][0], c1);

                // 	      MSG("FOUND EDGE %d/%d <-> %d/%d\n",
                // 		  edge0.first->getIndex(), edge0.second,
                // 		  edge1.first->getIndex(), edge1.second);

                // 	      MSG("FOUND EDGE: %d %d with coords %f %f %f   %f %f %f\n",
                // 		  dofEdge0.first, dofEdge0.second,
                // 		  c0[0], c0[1], c0[2],
                // 		  c1[0], c1[1], c1[2]);

                TEST_EXIT_DBG(dofEdge0 == dofEdge1)("Should noth happen!\n");
#endif

                FixRefinementPatch::connectedEdges[meshes[i]].push_back(make_pair(edge1, edge0));
              }
            }
          }
        }
      }


      MSG("Fix 3D mesh refinement needed %.5f seconds\n", t.elapsed());
    }


    void MeshDistributor::helpToFix(std::set<int>& elems,
                                    vector<std::set<int>>& disconnectedEls)
    {
      std::set<int> firstElem;
      firstElem.insert(*(elems.begin()));

      std::set<int> otherElems(++(elems.begin()), elems.end());

      bool found = false;
      do
      {
        found = false;
        for (std::set<int>::iterator elIt = firstElem.begin(); elIt != firstElem.end(); ++elIt)
        {
          for (int i = 0; i < macroMesh->getGeo(NEIGH); i++)
          {
            int neighEl = macroElementNeighbours[*elIt][i];
            if (neighEl != -1 && otherElems.count(neighEl))
            {
              firstElem.insert(neighEl);
              otherElems.erase(neighEl);
              found = true;
            }
          }

          if (found)
            break;
        }
      }
      while (found);

      disconnectedEls.push_back(firstElem);
      if (otherElems.size())
        helpToFix(otherElems, disconnectedEls);
    }


    void MeshDistributor::getAllBoundaryDofs(const FiniteElemSpace* feSpace,
        int level,
        DofContainer& dofs)
    {
      DofContainerSet dofSet;
      MultiLevelDofComm& dofComm = dofComms[feSpace->getMesh()];
      for (DofComm::Iterator it(dofComm[level].getSendDofs(), feSpace);
           !it.end(); it.nextRank())
        dofSet.insert(it.getDofs().begin(), it.getDofs().end());
      for (DofComm::Iterator it(dofComm[level].getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
        dofSet.insert(it.getDofs().begin(), it.getDofs().end());

      dofs.clear();
      dofs.insert(dofs.begin(), dofSet.begin(), dofSet.end());
    }


    void MeshDistributor::removePeriodicBoundaryConditions()
    {
      // Remove periodic boundaries in boundary manager on matrices and vectors.
      for (size_t i = 0; i < problemStat.size(); i++)
        removePeriodicBoundaryConditions(problemStat[i]);

      // Remove periodic boundaries on elements in mesh.
      for (size_t i = 0; i < meshes.size(); i++)
      {

        TraverseStack stack;
        ElInfo* elInfo = stack.traverseFirst(meshes[i], -1, Mesh::CALL_EVERY_EL_PREORDER);
        while (elInfo)
        {
          elInfo->getElement()->deleteElementData(PERIODIC);
          elInfo = stack.traverseNext(elInfo);
        }

        // Remove periodic vertex associations
        meshes[i]->getPeriodicAssociations().clear();
      }
    }


    void MeshDistributor::removePeriodicBoundaryConditions(ProblemStatSeq* probStat)
    {
      int nComponents = probStat->getNumComponents();

      for (int j = 0; j < nComponents; j++)
      {
        for (int k = 0; k < nComponents; k++)
        {
          DOFMatrix* mat = probStat->getSystemMatrix(j, k);
          if (mat && mat->getBoundaryManager())
            removePeriodicBoundaryConditions(const_cast<BoundaryIndexMap&>(mat->getBoundaryManager()->getBoundaryConditionMap()));
        }

        if (probStat->getSolution()->getDOFVector(j)->getBoundaryManager())
          removePeriodicBoundaryConditions(const_cast<BoundaryIndexMap&>(probStat->getSolution()->getDOFVector(j)->getBoundaryManager()->getBoundaryConditionMap()));

        if (probStat->getRhs()->getDOFVector(j)->getBoundaryManager())
          removePeriodicBoundaryConditions(const_cast<BoundaryIndexMap&>(probStat->getRhs()->getDOFVector(j)->getBoundaryManager()->getBoundaryConditionMap()));
      }
    }


    void MeshDistributor::removePeriodicBoundaryConditions(BoundaryIndexMap& boundaryMap)
    {
      BoundaryIndexMap::iterator it = boundaryMap.begin();
      while (it != boundaryMap.end())
      {
        if (it->second->isPeriodic())
          boundaryMap.erase(it++);
        else
          ++it;
      }
    }


    void MeshDistributor::createMeshLevelStructure()
    {
      FUNCNAME("MeshDistributor::createMeshLevelStructure()");

      int levelMode = 0;
      Parameters::get("parallel->level mode", levelMode);

      TEST_EXIT(levelMode >= 0)("Wrong level mode %d!\n", levelMode);

      if (levelMode == 0)
      {
        std::set<int> neighbours;
        levelData.init(neighbours);
        return;
      }

      if (levelMode == 1)
      {
        std::set<int> neighbours;
        levelData.init(neighbours);
        levelData.addLevelMode1();
        return;
      }

      TEST_EXIT(MPI::COMM_WORLD.Get_size() == 16)
      ("Only mpiSize == 16 supported!\n");

      if (levelMode == 2)
      {

        std::set<int> neighbours;
        switch (mpiRank)
        {
        case 0:
          neighbours.insert(1);
          neighbours.insert(2);
          neighbours.insert(3);
          break;
        case 1:
          neighbours.insert(0);
          neighbours.insert(4);
          neighbours.insert(2);
          neighbours.insert(3);
          neighbours.insert(6);
          break;
        case 2:
          neighbours.insert(0);
          neighbours.insert(1);
          neighbours.insert(3);
          neighbours.insert(8);
          neighbours.insert(9);
          break;
        case 3:
          neighbours.insert(0);
          neighbours.insert(1);
          neighbours.insert(4);
          neighbours.insert(2);
          neighbours.insert(6);
          neighbours.insert(8);
          neighbours.insert(9);
          neighbours.insert(12);
          break;
        case 4:
          neighbours.insert(1);
          neighbours.insert(3);
          neighbours.insert(6);
          neighbours.insert(7);
          neighbours.insert(5);
          break;
        case 5:
          neighbours.insert(4);
          neighbours.insert(6);
          neighbours.insert(7);
          break;
        case 6:
          neighbours.insert(1);
          neighbours.insert(4);
          neighbours.insert(5);
          neighbours.insert(3);
          neighbours.insert(7);
          neighbours.insert(9);
          neighbours.insert(12);
          neighbours.insert(13);
          break;
        case 7:
          neighbours.insert(4);
          neighbours.insert(5);
          neighbours.insert(6);
          neighbours.insert(12);
          neighbours.insert(13);
          break;
        case 8:
          neighbours.insert(2);
          neighbours.insert(3);
          neighbours.insert(9);
          neighbours.insert(10);
          neighbours.insert(11);
          break;
        case 9:
          neighbours.insert(2);
          neighbours.insert(3);
          neighbours.insert(6);
          neighbours.insert(8);
          neighbours.insert(12);
          neighbours.insert(10);
          neighbours.insert(11);
          neighbours.insert(14);
          break;
        case 10:
          neighbours.insert(8);
          neighbours.insert(9);
          neighbours.insert(11);
          break;
        case 11:
          neighbours.insert(8);
          neighbours.insert(9);
          neighbours.insert(12);
          neighbours.insert(10);
          neighbours.insert(14);
          break;
        case 12:
          neighbours.insert(3);
          neighbours.insert(6);
          neighbours.insert(7);
          neighbours.insert(9);
          neighbours.insert(13);
          neighbours.insert(11);
          neighbours.insert(14);
          neighbours.insert(15);
          break;
        case 13:
          neighbours.insert(6);
          neighbours.insert(7);
          neighbours.insert(12);
          neighbours.insert(14);
          neighbours.insert(15);
          break;
        case 14:
          neighbours.insert(9);
          neighbours.insert(12);
          neighbours.insert(13);
          neighbours.insert(11);
          neighbours.insert(15);
          break;
        case 15:
          neighbours.insert(12);
          neighbours.insert(13);
          neighbours.insert(14);
          break;
        }

        TEST_EXIT(neighbours.size() > 0)("Should not happen!\n");

        levelData.init(neighbours);

        switch (mpiRank)
        {
        case 0:
        case 1:
        case 2:
        case 3:
        {
          std::set<int> m;
          m.insert(0);
          m.insert(1);
          m.insert(2);
          m.insert(3);
          levelData.addLevel(m, 0);
        }
        break;
        case 4:
        case 5:
        case 6:
        case 7:
        {
          std::set<int> m;
          m.insert(4);
          m.insert(5);
          m.insert(6);
          m.insert(7);
          levelData.addLevel(m, 1);
        }
        break;
        case 8:
        case 9:
        case 10:
        case 11:
        {
          std::set<int> m;
          m.insert(8);
          m.insert(9);
          m.insert(10);
          m.insert(11);
          levelData.addLevel(m, 2);
        }
        break;
        case 12:
        case 13:
        case 14:
        case 15:
        {
          std::set<int> m;
          m.insert(12);
          m.insert(13);
          m.insert(14);
          m.insert(15);
          levelData.addLevel(m, 3);
        }
        break;
        }

        levelData.addLevelMode1();

        return;
      }
    }


    void MeshDistributor::checkMeshChange(bool tryRepartition)
    {
      FUNCNAME("MeshDistributor::checkMeshChange()");

      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;
      mpiComm.Barrier();
      Timer t;

      int skip = 0;
      Parameters::get("parallel->debug->skip check mesh change", skip);

      // === If mesh has not been changed on all ranks, return. ===
      vector<int> meshAllValues(meshes.size(), 0);

      for (size_t i = 0; i < meshes.size(); i++)
      {
#if (DEBUG != 0)
        MSG("mesh[%d] change index = %d, stored last index = %d.\n",
            i, meshes[i]->getChangeIndex(), lastMeshChangeIndexs[meshes[i]]);
#endif
        int sendValue =
          static_cast<int>(meshes[i]->getChangeIndex() != lastMeshChangeIndexs[meshes[i]]);
        mpiComm.Allreduce(&sendValue, &meshAllValues[i], 1, MPI_INT, MPI_SUM);
      }

      if (meshAllValues.size() == count(meshAllValues.begin(), meshAllValues.end(), 0))
        return;

      TEST_EXIT(levelData.getNumberOfLevels() == 1)
      ("Not yet implemented for multi level!\n");

      // === At least one rank mesh has been changed, so the boundaries must be ===
      // === adapted to the new mesh structure.                                 ===

      if (skip != 0)
      {
        for (size_t i = 0; i < meshes.size(); i++)
          if (meshAllValues[i] != 0)
            updateDofRelatedStruct(meshes[i]);
      }
      else
      {
        // To check the interior boundaries, the ownership of the boundaries is not
        // important. Therefore, we add all boundaries to one boundary container.
        RankToBoundMap allBound;

        for (InteriorBoundary::iterator it(intBoundary[0].getOwn()); !it.end(); ++it)
          if ((macroMesh->getDim() == 2 && it->rankObj.subObj == EDGE) ||
              (macroMesh->getDim() == 3 && it->rankObj.subObj == FACE))
            allBound[it.getRank()].push_back(*it);

        for (InteriorBoundary::iterator it(intBoundary[0].getOther());
             !it.end(); ++it)
          if ((macroMesh->getDim() == 2 && it->rankObj.subObj == EDGE) ||
              (macroMesh->getDim() == 3 && it->rankObj.subObj == FACE))
            allBound[it.getRank()].push_back(*it);

        for (InteriorBoundary::iterator it(intBoundary[0].getPeriodic());
             !it.end(); ++it)
          if (it.getRank() != mpiRank)
            if ((macroMesh->getDim() == 2 && it->rankObj.subObj == EDGE) ||
                (macroMesh->getDim() == 3 && it->rankObj.subObj == FACE))
              allBound[it.getRank()].push_back(*it);

        for (size_t i = 0; i < meshes.size(); i++)
        {
          if (meshAllValues[i] == 0)
            continue;

          int iterationCounter = 0;
          do
          {
            bool meshChanged = false;

            // === Check the boundaries and adapt mesh if necessary. ===
#if (DEBUG != 0)
            MSG("Run checkAndAdaptBoundary for mesh[%d]...\n", i);
#endif

            // Check for periodic boundaries within rank's subdomain.
            for (InteriorBoundary::iterator it(intBoundary[0].getPeriodic());
                 !it.end(); ++it)
            {
              if (it.getRank() == mpiRank)
              {
                if ((macroMesh->getDim() == 2 && it->rankObj.subObj == EDGE) ||
                    (macroMesh->getDim() == 3 && it->rankObj.subObj == FACE))
                {
                  MeshStructure elCode;
                  elCode.init(it->rankObj, elObjDb.getElementPtr(it->rankObj.elIndex, meshes[i]));

                  MeshManipulation mm(meshes[i]);
                  it->neighObj.el = elObjDb.getElementPtr(it->neighObj.elIndex, meshes[i]);
                  meshChanged |= mm.fitElementToMeshCode(elCode, it->neighObj);
                }
              }
            }

            meshChanged |= checkAndAdaptBoundary(allBound, meshes[i]);


            // === Check on all ranks if at least one rank's mesh has changed. ===

            int sendValue = static_cast<int>(meshChanged);
            meshAllValues[i] = 0;
            mpiComm.Allreduce(&sendValue, &meshAllValues[i], 1, MPI_INT, MPI_SUM);

            MSG("Mesh changed on %d ranks!\n", meshAllValues[i]);
            iterationCounter++;
          }
          while (meshAllValues[i] != 0);

          MSG("Number of iteration to adapt mesh: %d\n", iterationCounter);
          updateDofRelatedStruct(meshes[i]);
        }
      }
      mpiComm.Barrier();
      MSG("Parallel mesh adaption needed %.5f seconds\n", t.elapsed());


      // === Update the DOF numbering and mappings. ===

      updateLocalGlobalNumbering();

#if (DEBUG != 0)
      debug::writeMesh(feSpaces[0], -1, debugOutputDir + "mesh");
      ParallelDebug::testPeriodicBoundary(*this);
#endif

      // === The mesh has changed, so check if it is required to repartition ===
      // === the mesh.                                                       ===

      nMeshChangesAfterLastRepartitioning++;


      if (repartitioningFailed > 0)
      {
#if (DEBUG != 0)
        MSG("Repartitioning not tried because it has failed in the past!\n");
#endif

        repartitioningFailed--;
      }
      else if (tryRepartition &&
               repartitioningAllowed &&
               (repartitionOnlyOnce || nMeshChangesAfterLastRepartitioning >= repartitionIthChange))
      {
        repartitionMesh();
        nMeshChangesAfterLastRepartitioning = 0;
        if (repartitionOnlyOnce)
          repartitioningAllowed = false;
      }
      else
      {
        MSG("Repartitioning not tried because tryRepartitioning = %d   repartitioningAllowed = %d   nMeshChange = %d   repartitionIthChange = %d\n",
            tryRepartition, repartitioningAllowed,
            nMeshChangesAfterLastRepartitioning, repartitionIthChange);
      }

      // === Print imbalance factor. ===

      printImbalanceFactor();
    }


    void MeshDistributor::getImbalanceFactor(double& imbalance,
        int& minDofs,
        int& maxDofs,
        int& sumDofs)
    {
      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;
      int mpiSize = mpiComm.Get_size();

      vector<int> nDofsInRank(mpiSize);
      int nDofs = 0;

      for (size_t i = 0; i < meshes.size(); i++)
        nDofs += meshes[i]->getDofAdmin(0).getUsedDofs();

      mpiComm.Gather(&nDofs, 1, MPI_INT, &(nDofsInRank[0]), 1, MPI_INT, 0);

      if (mpiRank == 0)
      {
        sumDofs = 0;
        minDofs = numeric_limits<int>::max();
        maxDofs = numeric_limits<int>::min();

        for (int i = 0; i < mpiSize; i++)
        {
          sumDofs += nDofsInRank[i];
          minDofs = std::min(minDofs, nDofsInRank[i]);
          maxDofs = std::max(maxDofs, nDofsInRank[i]);
        }
        int avrgDofs = sumDofs / mpiSize;
        imbalance = ((static_cast<double>(maxDofs) / avrgDofs) - 1.0);
      }
    }


    double MeshDistributor::getImbalanceFactor()
    {
      double factor;
      int a = 0;
      int b = 0;
      int c = 0;
      getImbalanceFactor(factor, a, b, c);
      return factor;
    }


    void MeshDistributor::printImbalanceFactor()
    {
      FUNCNAME("MeshDistributor::printImbalanceFactor()");

      double imbalanceFactor = 0.0;
      int minDofs = 0;
      int maxDofs = 0;
      int sumDofs = 0;

      getImbalanceFactor(imbalanceFactor, minDofs, maxDofs, sumDofs);
      if (mpiRank == 0)
        MSG("Imbalancing factor: %.2f %%  [ minDofs = %d, maxDofs = %d, sumDofs = %d ]\n",
            imbalanceFactor * 100.0, minDofs, maxDofs, sumDofs);
    }


    bool MeshDistributor::checkAndAdaptBoundary(RankToBoundMap& allBound, Mesh* mesh)
    {
      FUNCNAME_DBG("MeshDistributor::checkAndAdaptBoundary()");

      // === Create mesh structure codes for all ranks boundary elements. ===

      map<int, MeshCodeVec> sendCodes;

      for (RankToBoundMap::iterator it = allBound.begin();
           it != allBound.end(); ++it)
      {
        for (vector<AtomicBoundary>::iterator boundIt = it->second.begin();
             boundIt != it->second.end(); ++boundIt)
        {
          MeshStructure elCode;
          elCode.init(boundIt->rankObj, elObjDb.getElementPtr(boundIt->rankObj.elIndex, mesh));
          sendCodes[it->first].push_back(elCode);
        }
      }

      StdMpi<MeshCodeVec> stdMpi(MPI::COMM_WORLD, true);
      stdMpi.send(sendCodes);
      for (RankToBoundMap::iterator it = allBound.begin();
           it != allBound.end(); ++it)
        stdMpi.recv(it->first);

      stdMpi.startCommunication();


      // === Compare received mesh structure codes. ===

      bool meshChanged = false;

      for (RankToBoundMap::iterator it = allBound.begin();
           it != allBound.end(); ++it)
      {

        MeshCodeVec& recvCodes = stdMpi.getRecvData()[it->first];
        int i = 0;

        for (vector<AtomicBoundary>::iterator boundIt = it->second.begin();
             boundIt != it->second.end(); ++boundIt, i++)
        {

          MeshStructure elCode;
          elCode.init(boundIt->rankObj, elObjDb.getElementPtr(boundIt->rankObj.elIndex, mesh));

#if (DEBUG != 0)
          ParallelDebug::followBoundary(mesh, *boundIt, elCode);
#endif

          if (elCode.getCode() != recvCodes[i].getCode())
          {
            // 	  MSG("MACRO EL %d %d %d DOES NOT FIT WITH MACRO EL %d %d %d ON RANK %d\n",
            // 	      boundIt->rankObj.elIndex, boundIt->rankObj.subObj, boundIt->rankObj.ithObj,
            // 	      boundIt->neighObj.elIndex, boundIt->neighObj.subObj, boundIt->neighObj.ithObj,
            // 	      it->first);
            TEST_EXIT_DBG(refineManager)("Refinement manager is not set correctly!\n");

            MeshManipulation mm(mesh);
            boundIt->rankObj.el = elObjDb.getElementPtr(boundIt->rankObj.elIndex, mesh);
            meshChanged |= mm.fitElementToMeshCode(recvCodes[i], boundIt->rankObj);
            boundIt->rankObj.el = NULL;
          }
        }
      }

      return meshChanged;
    }


    void MeshDistributor::serialize(ostream& out, DofContainer& data)
    {
      int vecSize = data.size();
      SerUtil::serialize(out, vecSize);
      for (int i = 0; i < vecSize; i++)
      {
        int dofIndex = *(data[i]);
        SerUtil::serialize(out, dofIndex);
      }
    }


    void MeshDistributor::serialize(ostream& out,
                                    map<int, map<const FiniteElemSpace*, DofContainer>>& data)
    {
      int mapSize = data.size();
      SerUtil::serialize(out, mapSize);
      for (map<int, map<const FiniteElemSpace*, DofContainer>>::iterator it = data.begin();
           it != data.end(); ++it)
      {
        int rank = it->first;
        SerUtil::serialize(out, rank);

        for (map<const FiniteElemSpace*, DofContainer>::iterator dcIt = it->second.begin();
             dcIt != it->second.end(); ++dcIt)
          serialize(out, dcIt->second);
      }
    }


    void MeshDistributor::deserialize(istream& in, DofContainer& data,
                                      map<int, const DegreeOfFreedom*>& dofIndexMap)
    {
      FUNCNAME("MeshDistributor::deserialize()");

      int vecSize = 0;
      SerUtil::deserialize(in, vecSize);
      data.clear();
      data.resize(vecSize);
      for (int i = 0; i < vecSize; i++)
      {
        int dofIndex = 0;
        SerUtil::deserialize(in, dofIndex);

        TEST_EXIT_DBG(dofIndexMap.count(dofIndex) != 0)
        ("Dof index could not be deserialized correctly!\n");

        data[i] = dofIndexMap[dofIndex];
      }
    }


    void MeshDistributor::deserialize(istream& in,
                                      map<int, map<const FiniteElemSpace*, DofContainer>>& data,
                                      map<const FiniteElemSpace*, map<int, const DegreeOfFreedom*>>& dofIndexMap)
    {
      FUNCNAME("MeshDistributor::deserialize()");

      ERROR_EXIT("Must be reimplemented!\n");
#if 0
      data.clear();

      int mapSize = 0;
      SerUtil::deserialize(in, mapSize);
      for (int i = 0; i < mapSize; i++)
      {
        int rank = 0;
        SerUtil::deserialize(in, rank);

        for (unsigned int j = 0; j < componentSpaces.size(); j++)
          deserialize(in,
                      data[rank][componentSpaces[j]],
                      dofIndexMap[componentSpaces[j]]);
      }
#endif
    }

    bool MeshDistributor::isRepartitionNecessary()
    {
      FUNCNAME("MeshDistributor::isRepartitionNecessary()");

      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;

      // === First, check if the load is unbalanced on the ranks. ===

      int repartitioning = 0;
      double imbalanceFactor = getImbalanceFactor();
      double imbalanceRepartitionBound = 0.2;
      Parameters::get("parallel->repartitioning->imbalance",
                      imbalanceRepartitionBound);

      if (mpiRank == 0)
      {
        if (imbalanceFactor > imbalanceRepartitionBound)
          repartitioning = 1;

        mpiComm.Bcast(&repartitioning, 1, MPI_INT, 0);
      }
      else
      {
        mpiComm.Bcast(&repartitioning, 1, MPI_INT, 0);
      }
#if (DEBUG != 0)
      if (repartitioning == 0)
      {
        MSG("imbalanceFactor = %f < %f = imbalanceRepartitionBound\n", imbalanceFactor, imbalanceRepartitionBound);
      }
#endif
      return (repartitioning != 0);
    }

    void MeshDistributor::calculateElemWeights()
    {
      FUNCNAME("MeshDistributor::calculateElemWeights()");

      elemWeights.clear();
      {
        for (size_t i = 0; i < meshes.size(); i++)
        {
          TraverseStack stack;
          ElInfo* elInfo = stack.traverseFirst(meshes[i], -1, Mesh::CALL_LEAF_EL);
          while (elInfo)
          {
            elemWeights[elInfo->getMacroElement()->getIndex()]++;
            elInfo = stack.traverseNext(elInfo);
          }
        }
      }

      double maxWeight = -1.0;
      double sumWeight = 0.0;
      for (map<int, double>::iterator it = elemWeights.begin();
           it != elemWeights.end(); ++it)
      {
        maxWeight = std::max(maxWeight, it->second);
        sumWeight += it->second;
      }

      mpi::globalMax(maxWeight);
      mpi::globalAdd(sumWeight);

      MSG("Partition weight: sum = %e max = %e\n", sumWeight, maxWeight);
    }


    bool MeshDistributor::repartitionMesh()
    {
      FUNCNAME("MeshDistributor::repartitionMesh()");

      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;

      // === First, check if the load is unbalanced on the ranks. ===

      if (!isRepartitionNecessary())
        return false;

      Timer t;

#if (DEBUG != 0)
      for (size_t i = 0; i < meshes.size(); i++)
        ParallelDebug::testDoubleDofs(meshes[i]);
      int writePartMesh = 1;
#else
      int writePartMesh = 0;
#endif
      Parameters::get("parallel->debug->write part mesh", writePartMesh);
      if (writePartMesh > 0 && repartitioningCounter == 0)
        ParallelDebug::writePartitioningFile(debugOutputDir + "partitioning",
                                             repartitioningCounter,
                                             feSpaces[0]);

      repartitioningCounter++;

      // === Create new element weights. ===

      int externalWeight = 0;
      Parameters::get("parallel->use external weights", externalWeight);
      if (!externalWeight)
        calculateElemWeights();

      // === Run mesh partitioner to calculate a new mesh partitioning.  ===

      TEST_EXIT(dofMaps.size())("No DOF mapping defined!\n");
      ParallelDofMapping& dofMap = *(dofMaps[0]);
      //     partitioner->setLocalGlobalDofMap(&(dofMap[feSpaces[0]].getMap())); // TODO: hier kommt es zu Problemen!!!
      partitioner->setLocalGlobalDofMap(&(dofMap[0].getMap()));
      bool partitioningSucceed =
        partitioner->partition(elemWeights, ADAPTIVE_REPART);
      if (!partitioningSucceed)
      {
        mpiComm.Barrier();
        repartitioningFailed = repartitioningWaitAfterFail;
        MSG("Mesh partitioner created empty partition!\n");
        MSG("Mesh repartitioning needed %.5f seconds\n", t.elapsed());
        return false;
      }


      // In the case the partitioner does not create a new mesh partition, return
      // without and changes.
      if (!partitioner->meshChanged())
      {
        MSG("Mesh partition does not create a new partition!\n");
        MSG("Try to refine partitioning!\n");
        partitioningSucceed = partitioner->partition(elemWeights, REFINE_PART);
        if (!partitioningSucceed || !partitioner->meshChanged())
        {
          mpiComm.Barrier();
          repartitioningFailed = repartitioningWaitAfterFail;;
          MSG("Mesh repartitioning needed %.5f seconds\n", t.elapsed());
          return false;
        }
      }

      TEST_EXIT_DBG(!(partitioner->getSendElements().size() == macroMesh->getMacroElements().size() &&
                      partitioner->getRecvElements().size() == 0))
      ("Partition is empty, should not happen!\n");

      int strategy = 0;
      Parameters::get("parallel->repartitioning->strategy", strategy);

      switch (strategy)
      {
      case 0:
        // Always starts from meshes[0], the macro mesh
        for (size_t i = 0; i < meshes.size(); i++)
          quickRepartition(meshes[i]);
        break;
      case 1:
        for (size_t i = 0; i < meshes.size(); i++)
          fullRepartition(meshes[i]);
        break;
      default:
        ERROR_EXIT("Unknown repartition strategy = %d!\n", strategy);
      }
      // This is done outside repartition, don't forget!
      updateLocalGlobalNumbering();


#if (DEBUG != 0)
      MSG("AMDiS runs in debug mode, so make some test ...\n");

      ParallelDebug::writePartitioningFile(debugOutputDir + "partitioning",
                                           repartitioningCounter,
                                           feSpaces[0]);
      ParallelDebug::testAllElements(*this);
      for (size_t i = 0; i < meshes.size(); i++)
        ParallelDebug::testDoubleDofs(meshes[i]);
      ParallelDebug::testInteriorBoundary(*this);
      ParallelDebug::testPeriodicBoundary(*this);

      MSG("Debug mode tests finished!\n");
#endif


      mpiComm.Barrier();
      MSG("Mesh repartitioning needed %.5f seconds\n", t.elapsed());
      return true;
    }

    void MeshDistributor::quickRepartition(Mesh* mesh)
    {
      FUNCNAME("MeshDistributor::quickRepartition()");

#if (DEBUG != 0)
      MSG("... Run quickRepartition ...\n");
#endif

      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;
      // === Create map that maps macro element indices to pointers to the ===
      // === macro elements.                                               ===

      map<int, MacroElement*> elIndexMap;
      for (vector<MacroElement*>::iterator it = allMacroElements[mesh].begin();
           it != allMacroElements[mesh].end(); ++it)
        elIndexMap[(*it)->getIndex()] = *it;


      // === Create set of all new macro elements this rank will receive from ===
      // === other ranks.                                                     ===

      vector<MacroElement*> newMacroEl;
      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
      {
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {
          TEST_EXIT_DBG(elIndexMap.count(*elIt) == 1)
          ("Could not find macro element %d\n", *elIt);

          TEST_EXIT_DBG(find(newMacroEl.begin(), newMacroEl.end(),
                             elIndexMap[*elIt]) == newMacroEl.end())
          ("The mesh received the same macro element %d twice! Something is wrong...\n", *elIt);

          newMacroEl.push_back(elIndexMap[*elIt]);
        }
      }

      std::set<MacroElement*> allMacroEl;
      for (vector<MacroElement*>::iterator it = newMacroEl.begin();
           it != newMacroEl.end(); ++it)
        allMacroEl.insert(*it);

      for (deque<MacroElement*>::iterator it = mesh->firstMacroElement();
           it != mesh->endOfMacroElements(); ++it)
        allMacroEl.insert(*it);


      // === Add new macro elements to mesh. ===
      for (vector<MacroElement*>::iterator it = newMacroEl.begin();
           it != newMacroEl.end(); ++it)
      {
        MacroElement* mel = *it;

        // First, reset all neighbour relations. The correct neighbours will be
        // set later.
        for (int i = 0; i < mesh->getGeo(NEIGH); i++)
          mel->setNeighbour(i, NULL);

        // Create new DOFs for the macro element.
        mel->getElement()->createNewDofPtrs(true);

        // Push the macro element to all macro elements in mesh.
        mesh->getMacroElements().push_back(mel);

      }


      // === Send and receive mesh structure codes. ===

      map<int, MeshCodeVec> sendCodes;
      map<int, vector<vector<double>>> sendValues;

      TEST_EXIT(interchangeVectors.size() > 0)
      ("There are no interchange vectors defined!\n");

      for (map<int, vector<int>>::iterator it = partitioner->getSendElements().begin();
           it != partitioner->getSendElements().end(); ++it)
      {
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {
          MeshStructure elCode;
          elCode.init(mesh, *elIt);
          sendCodes[it->first].push_back(elCode);

          for (size_t i = 0; i < interchangeVectors.size(); i++)
          {
            //only the vectors defined on the corresponding mesh
            if (interchangeVectors[i]->getFeSpace()->getMesh() != mesh)
              continue;

            vector<double> valVec;
            elCode.getMeshStructureValues(*elIt, interchangeVectors[i], valVec);
            sendValues[it->first].push_back(valVec);
          }
          TEST_EXIT_DBG(!sendValues[it->first].empty())
          ("There is no DOFVector defined on this mesh.\n");
        }
      }

      StdMpi<MeshCodeVec> stdMpi(mpiComm, true);
      stdMpi.send(sendCodes);
      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
        stdMpi.recv(it->first);
      stdMpi.startCommunication();

      StdMpi<vector<vector<double>>> stdMpi2(mpiComm, true);
      stdMpi2.send(sendValues);
      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
        stdMpi2.recv(it->first);
      stdMpi2.startCommunication();


      // === Adapte the new macro elements due to the received mesh   ===
      // === structure codes.                                         ===

      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
      {
        MeshCodeVec& recvCodes = stdMpi.getRecvData()[it->first];
        int i = 0;

        TEST_EXIT_DBG(recvCodes.size() == it->second.size())
        ("Should not happen!\n");

        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
          recvCodes[i++].fitMeshToStructure(mesh, refineManager, false, *elIt);
      }


      // === Set correct neighbour information on macro elements. ===

      for (std::set<MacroElement*, bool(*)(MacroElement*, MacroElement*)>::iterator it = allMacroEl.begin();
           it != allMacroEl.end(); ++it)
      {
        int elIndex = (*it)->getIndex();
        for (int i = 0; i < mesh->getGeo(NEIGH); i++)
        {

          TEST_EXIT_DBG(macroElementNeighbours.count(elIndex) == 1)
          ("Should not happen!\n");
          TEST_EXIT_DBG(static_cast<int>(macroElementNeighbours[elIndex].size()) ==
                        mesh->getGeo(NEIGH))
          ("Should not happen!\n");

          int neighIndex = macroElementNeighbours[elIndex][i];

          if (neighIndex == -1 ||
              partitioner->getElementInRank()[neighIndex] == false)
          {
            (*it)->setNeighbour(i, NULL);
          }
          else
          {
            TEST_EXIT_DBG(elIndexMap.count(neighIndex) == 1)
            ("Should not happen!\n");

            (*it)->setNeighbour(i, elIndexMap[neighIndex]);
          }
        }
      }


      // === Delete macros from mesh. ===

      std::set<MacroElement*> deleteMacroElements;

      for (map<int, vector<int>>::iterator it = partitioner->getSendElements().begin();
           it != partitioner->getSendElements().end(); ++it)
      {
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {
          TEST_EXIT_DBG(elIndexMap.count(*elIt) == 1)
          ("Could not find macro element %d\n", *elIt);

          deleteMacroElements.insert(elIndexMap[*elIt]);
        }
      }


      // Note that also if there are no macros to be deleted, this function will
      // update the number of elements, vertices, etc. of the mesh.
      mesh->removeMacroElements(deleteMacroElements, meshToFeSpaces[mesh]);


      // === Remove double DOFs. ===

      MeshManipulation meshManipulation(mesh);
      meshManipulation.deleteDoubleDofs(meshToFeSpaces[mesh], newMacroEl, elObjDb);
      mesh->dofCompress();


      // === After mesh adaption, set values of interchange vectors. ===

      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
      {

        MeshCodeVec& recvCodes = stdMpi.getRecvData()[it->first];
        vector<vector<double>>& recvValues = stdMpi2.getRecvData()[it->first];

        int i = 0, j = 0;
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {

          for (size_t k = 0; k < interchangeVectors.size(); k++)
          {
            if (interchangeVectors[k]->getFeSpace()->getMesh() != mesh)
              continue;

            recvCodes[i].setMeshStructureValues(*elIt,
                                                interchangeVectors[k],
                                                recvValues[j++]);
          }
          i++;
        }
      }

      if (mesh == macroMesh)
      {
        // In 3D we have to make some test, if the resulting mesh is valid.
        fix3dMeshRefinement();
        // PartitionMap is needed for createInteriorBoundary
        partitioner->createPartitionMap(partitionMap);
        // InteriorBoundary is needed for updateDofRelatedStruct
        createInteriorBoundary(false);
      }
      updateDofRelatedStruct(mesh);

    }

    void MeshDistributor::fullRepartition(Mesh* mesh)
    {
      FUNCNAME("MeshDistributor::fullRepartition()");

#if (DEBUG != 0)
      MSG("... Run fullRepartition ...\n");
#endif

      TEST_EXIT(interchangeVectors.size() > 0)
      ("There are no interchange vectors defined!\n");

      TEST_EXIT(mesh->getMacroFileInfo())
      ("Please turn on \"preserve macroFileInfo\". Should not happen.\n");

      MPI::Intracomm& mpiComm = MPI::COMM_WORLD;

#if (DEBUG != 0)
      int nOldLeaves = mesh->getNumberOfLeaves();
      mpi::globalAdd(mpiComm, nOldLeaves);
#endif

      // === Create map that maps macro element indices to pointers to the ===
      // === macro elements.                                               ===

      map<int, MacroElement*> elIndexMap;
      for (vector<MacroElement*>::iterator it = allMacroElements[mesh].begin();
           it != allMacroElements[mesh].end(); ++it)
        elIndexMap[(*it)->getIndex()] = *it;


      // === Create set of all new macro elements this rank will receive from ===
      // === other ranks.                                                     ===

      std::set<MacroElement*> newMacroEl;
      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
      {
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {
          TEST_EXIT_DBG(elIndexMap.count(*elIt) == 1)
          ("Could not find macro element %d\n", *elIt);

          newMacroEl.insert(elIndexMap[*elIt]);
        }
      }

      // === Store mesh structure code and value for domain macro elements. ===

      std::set<MacroElement*> allMacroEl;
      for (deque<MacroElement*>::iterator it = mesh->firstMacroElement();
           it != mesh->endOfMacroElements(); ++it)
        allMacroEl.insert(*it);


      map<int, MeshStructure> domainMacroCodes;
      map<int, vector<vector<double>>> domainMacroValues;

      for (std::set<MacroElement*>::iterator it = allMacroEl.begin();
           it != allMacroEl.end(); ++it)
      {

        MeshStructure elCode;
        vector<vector<double>> elValue;

        int macroId = (*it)->getElement()->getIndex();

        elCode.init(mesh, macroId);
        for (size_t i = 0; i < interchangeVectors.size(); i++)
        {
          if (interchangeVectors[i]->getFeSpace()->getMesh() != mesh)
            continue;

          vector<double> valVec;
          elCode.getMeshStructureValues(macroId, interchangeVectors[i], valVec);
          elValue.push_back(valVec);
        }

        domainMacroCodes[macroId] = elCode;
        domainMacroValues[macroId] = elValue;
      }

      // === Rebuild all macros in order to clean and renumber the dofs. ===

      mesh->removeAllMacroElements();

      TEST_EXIT_DBG(!mesh->getNumberOfMacros())
      ("Still %d macros remain?\n", mesh->getNumberOfMacros());

      // Add all macros again.
      for (vector<MacroElement*>::iterator it = allMacroElements[mesh].begin();
           it != allMacroElements[mesh].end(); ++it)
      {

        MacroElement* mel = *it;
        int elIndex = mel->getIndex();

        TEST_EXIT_DBG(macroElementNeighbours.count(elIndex) == 1)
        ("Macro (%d) is not recorded for its neighbour information. Should not happen.\n", elIndex);
        TEST_EXIT_DBG(static_cast<int>(macroElementNeighbours[elIndex].size()) ==
                      mesh->getGeo(NEIGH))
        ("Should not happen!\n");

        // Set initial neighbour information.
        for (int i = 0; i < mesh->getGeo(NEIGH); i++)
        {

          int neighIndex = macroElementNeighbours[elIndex][i];

          if (neighIndex == -1)
            mel->setNeighbour(i, NULL);
          else
          {
            TEST_EXIT_DBG(elIndexMap.count(neighIndex) == 1)
            ("Should not happen!\n");

            (*it)->setNeighbour(i, elIndexMap[neighIndex]);
          }
        }

        // Allocate dof pointers without initialization. Since neighbour info
        // is set, neighbours share common dofs
        mel->getElement()->createNewDofPtrs();

        // Push the macro element to all macro elements in mesh.
        mesh->getMacroElements().push_back(mel);
      }

      // New dof indices for all macro elements.
      io::MacroReader::restoreMacroDofs(*(mesh->getMacroFileInfo()), 0);

      // === Update the current mesh partition. ===

      for (map<int, vector<int>>::iterator it = partitioner->getSendElements().begin();
           it != partitioner->getSendElements().end(); ++it)
      {
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {
          TEST_EXIT_DBG(elIndexMap.count(*elIt) == 1)
          ("Could not find macro element %d\n", *elIt);

          allMacroEl.erase(elIndexMap[*elIt]);
        }
      }

      std::set<MacroElement*> localMacroEl = allMacroEl;

      for (std::set<MacroElement*>::iterator it = newMacroEl.begin();
           it != newMacroEl.end(); ++it)
        allMacroEl.insert(*it);


      // === Delete macros which don't belong to the domain in this ===
      // === round of partition.                                    ===

      std::set<MacroElement*> deleteMacroElements;

      for (vector<MacroElement*>::iterator it = allMacroElements[mesh].begin();
           it != allMacroElements[mesh].end(); ++it)
        if (allMacroEl.find(*it) == allMacroEl.end())
          deleteMacroElements.insert(*it);

      // Note that also if there are no macros to be deleted, this function will
      // update the number of elements, vertices, etc. of the mesh.
      mesh->removeMacroElements(deleteMacroElements, meshToFeSpaces[mesh]);
      mesh->dofCompress();
      // In fullRepartition, we need to call fix3d here instead of in the end.
      if (mesh == macroMesh)
        fix3dMeshRefinement();

      // === Adapte all domain macro elements. ===

      // Send and receive mesh structure codes.
      map<int, MeshCodeVec> sendCodes;
      map<int, vector<vector<double>>> sendValues;

      for (map<int, vector<int>>::iterator it = partitioner->getSendElements().begin();
           it != partitioner->getSendElements().end(); ++it)
      {
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {

          sendCodes[it->first].push_back(domainMacroCodes[*elIt]);
          sendValues[it->first].reserve(sendValues[it->first].size() + domainMacroValues.size());
          sendValues[it->first].insert(
            sendValues[it->first].end(), domainMacroValues[*elIt].begin(), domainMacroValues[*elIt].end());
        }
      }

      StdMpi<MeshCodeVec> stdMpi(mpiComm, true);
      stdMpi.send(sendCodes);
      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
        stdMpi.recv(it->first);
      stdMpi.startCommunication();

      StdMpi<vector<vector<double>>> stdMpi2(mpiComm, true);
      stdMpi2.send(sendValues);
      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
        stdMpi2.recv(it->first);
      stdMpi2.startCommunication();


      // === Adapte the new macro elements due to the received mesh   ===
      // === structure codes, reserved elements due to the local
      // === mesh structure codes.                                     ===

      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
      {
        MeshCodeVec& recvCodes = stdMpi.getRecvData()[it->first];
        int i = 0;

        TEST_EXIT_DBG(recvCodes.size() == it->second.size())
        ("Should not happen!\n");

        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {

          TEST_EXIT_DBG(allMacroEl.find(elIndexMap[*elIt]) != allMacroEl.end())
          ("Received a macro element[%d] which doesn't belong to the domain. Something wrong.\n", *elIt);


          recvCodes[i++].fitMeshToStructure(mesh, refineManager, false, *elIt);
        }
      }

      // === After mesh adaption, set values of interchange vectors. ===

      for (std::set<MacroElement*>::iterator it = localMacroEl.begin();
           it != localMacroEl.end(); ++it)
      {

        int macroId = (*it)->getIndex();

        TEST_EXIT_DBG(domainMacroCodes.count(macroId))
        ("Reserved macro element[%d] has no stored structure data. Something wrong.\n", macroId);


        domainMacroCodes[macroId].fitMeshToStructure(mesh, refineManager, false, macroId);

        int j = 0;
        for (size_t i = 0; i < interchangeVectors.size(); i++)
        {
          if (interchangeVectors[i]->getFeSpace()->getMesh() != mesh)
            continue;

          domainMacroCodes[macroId].setMeshStructureValues(macroId,
              interchangeVectors[i],
              domainMacroValues[macroId][j++]);

#if (DEBUG != 0)
          MeshStructure code;
          code.init(mesh, macroId);
          TEST_EXIT(code.getCode() == domainMacroCodes[macroId].getCode())
          ("Local mesh structure code not equal. Due to parallel adaption?\n");
#endif
        }
      }

      for (map<int, vector<int>>::iterator it = partitioner->getRecvElements().begin();
           it != partitioner->getRecvElements().end(); ++it)
      {
        MeshCodeVec& recvCodes = stdMpi.getRecvData()[it->first];
        vector<vector<double>>& recvValues = stdMpi2.getRecvData()[it->first];
        int i = 0, j = 0;
        for (vector<int>::iterator elIt = it->second.begin();
             elIt != it->second.end(); ++elIt)
        {
          for (size_t k = 0; k < interchangeVectors.size(); k++)
          {
            if (interchangeVectors[k]->getFeSpace()->getMesh() != mesh)
              continue;

            recvCodes[i].setMeshStructureValues(*elIt,
                                                interchangeVectors[k],
                                                recvValues[j++]);
          }
          i++;
        }
      }

      if (mesh == macroMesh)
      {
        partitioner->createPartitionMap(partitionMap);
        createInteriorBoundary(false);
      }
      updateDofRelatedStruct(mesh);

#if (DEBUG != 0)
      int nNewLeaves = mesh->getNumberOfLeaves();
      mpi::globalAdd(mpiComm, nNewLeaves);

      TEST_EXIT(nOldLeaves == nNewLeaves)
      ("Overall number of leaves change from %d to %d\n", nOldLeaves, nNewLeaves);
#endif
    }


    void MeshDistributor::createInteriorBoundary(bool firstCall)
    {
      FUNCNAME("MeshDistributor::createInteriorBoundary()");

      if (firstCall)
        elObjDb.create(partitionMap, levelData);
      elObjDb.updateRankData();

#if (DEBUG != 0)
      if (mpiRank == 0)
        ParallelDebug::writePeriodicElObjInfo(*this, debugOutputDir);
#endif

      // === If requested, calculate and print memory usage of element  ===
      // === object database.                                           ===

      //TODO To be reimplemented
      if (printMemoryUsage && firstCall)
      {
        unsigned long memsize = elObjDb.calculateMemoryUsage();
        MSG("Memory usage of element object database = %5.f KByte\n",
            static_cast<double>(memsize / 1024.0));
      }


      intBoundary.create(levelData, elObjDb);

      for (int level = 0; level < levelData.getNumberOfLevels(); level++)
        ParallelDebug::printBoundaryInfo(intBoundary[level]);

      if (firstCall)
      {
        int tmpSend = static_cast<int>(intBoundary[0].hasPeriodic());
        mpi::globalAdd(MPI::COMM_WORLD, tmpSend);
        hasPeriodicBoundary = static_cast<bool>(tmpSend);
      }
    }


    void MeshDistributor::createBoundaryDofs(Mesh* mesh)
    {
      FUNCNAME("MeshDistributor::createBoundaryDofs()");

      Timer* t = NULL;
      if (printTimings)
        t = new Timer();

      // === Create all DOF communicators, if mesh is NULL. ===
      // === Otherwise create only special one.                ===
      if (!mesh)
      {
        std::map<Mesh*, MultiLevelDofComm>::iterator it = dofComms.begin();
        for (; it != dofComms.end(); it++)
          it->second.create(it->first, intBoundary);
      }
      else
      {
        TEST_EXIT_DBG(dofComms.find(mesh) != dofComms.end())("No such mesh.\n");
        dofComms[mesh].create(mesh, intBoundary);
      }

      // === If requested, create more information on communication DOFs. ===

      if (createBoundaryDofFlag.isSet(BOUNDARY_SUBOBJ_SORTED))
      {

        int nLevels = levelData.getNumberOfLevels();
        boundaryDofInfo.resize(nLevels);

        for (unsigned int i = 0; i < feSpaces.size(); i++)
        {
          const FiniteElemSpace* feSpace = feSpaces[i];

          for (int level = 0; level < nLevels; level++)
          {

            // === Clear data. ===
            for (int geo = FACE; geo >= VERTEX; geo--)
              boundaryDofInfo[level][feSpace].geoDofs[static_cast<GeoIndex>(geo)].clear();

            // === Create send DOFs. ===
            for (int geo = FACE; geo >= VERTEX; geo--)
            {
              for (InteriorBoundary::iterator it(intBoundary[level].getOwn());
                   !it.end(); ++it)
              {
                if (it->rankObj.subObj == geo)
                {
                  DofContainer dofs;

                  intBoundary[level].getElementPtr(it->rankObj.elIndex, feSpace->getMesh())
                  ->getAllDofs(feSpace, it->rankObj, dofs);

                  if (createBoundaryDofFlag.isSet(BOUNDARY_FILL_INFO_SEND_DOFS))
                    boundaryDofInfo[level][feSpace].
                    geoDofs[static_cast<GeoIndex>(geo)].insert(dofs.begin(), dofs.end());
                }
              }
            }

            // === Create recv DOFs. ===
            for (int geo = FACE; geo >= VERTEX; geo--)
            {
              for (InteriorBoundary::iterator it(intBoundary[level].getOther());
                   !it.end(); ++it)
              {
                if (it->rankObj.subObj == geo)
                {
                  DofContainer dofs;

                  intBoundary[level].getElementPtr(it->rankObj.elIndex, feSpace->getMesh())
                  ->getAllDofs(feSpace, it->rankObj, dofs);

                  if (createBoundaryDofFlag.isSet(BOUNDARY_FILL_INFO_RECV_DOFS))
                    boundaryDofInfo[level][feSpace].
                    geoDofs[static_cast<GeoIndex>(geo)].insert(dofs.begin(), dofs.end());
                }
              }
            }
          }
        }
      }

      if (printTimings && t)
      {
        MSG("MESHDIST 01 (createBoundaryDofs) needed %.5f seconds\n", t->elapsed());
        delete t;
      }
    }


    void MeshDistributor::removeMacroElements()
    {
      for (size_t i = 0; i < meshes.size(); i++)
      {
        std::set<MacroElement*> macrosToRemove;
        for (deque<MacroElement*>::iterator it = meshes[i]->firstMacroElement();
             it != meshes[i]->endOfMacroElements(); ++it)
          if (partitioner->getElementInRank()[(*it)->getIndex()] == false)
            macrosToRemove.insert(*it);

        meshes[i]->removeMacroElements(macrosToRemove, meshToFeSpaces[meshes[i]]);
      }
    }

    void MeshDistributor::updateDofRelatedStruct()
    {
      for (size_t i = 0; i < meshes.size(); i++)
        updateDofRelatedStruct(meshes[i]);

      updateLocalGlobalNumbering();
    }

    void MeshDistributor::updateDofRelatedStruct(Mesh* mesh)
    {
      FUNCNAME("MeshDistributor::updateDofRelatedStruct()");

      TEST_EXIT(mesh)("No mesh.\n");

      Timer t;

      mesh->dofCompress();

#if (DEBUG != 0)
      debug::ElementIdxToDofs elMap;
      debug::createSortedDofs(mesh, elMap);
#endif

      // === Update DOF communicator objects. ===
      createBoundaryDofs(mesh);


      // === Update all registered DOF mapping objects. ===
      updateDofsToDofMapping(mesh);


      // === Update DOF admins due to new number of DOFs. ===
      lastMeshChangeIndexs[mesh] = mesh->getChangeIndex();


#if (DEBUG != 0)
      static int fileNumber(0);		//improvised counter for adapt Iteration
      stringstream ss;
      ss << debugOutputDir << "elementMaps." << fileNumber ;

      // Write local Element Maps to csv file
      ParallelDebug::writeCsvElementMap(feSpaces[0],
                                        *(dofMaps[0]),
                                        ss.str(), "csv");
      fileNumber++;


      ParallelDebug::testDofContainerCommunication(*this, mesh);

      //     debug::writeElementIndexMesh(mesh, debugOutputDir + "elementIndex-" +
      // 				 std::to_string(mpiRank) + ".vtu");

      debug::testSortedDofs(mesh, elMap);

      int test = 0;
      Parameters::get("parallel->remove periodic boundary", test);
      if (test == 0)
        ParallelDebug::testCommonDofs(*this, mesh, true);
#endif

      MPI::COMM_WORLD.Barrier();
      MSG("Update parallel data structures needed %.5f seconds\n", t.elapsed());
    }

    void MeshDistributor::updateLocalGlobalNumbering()
    {
      FUNCNAME("MeshDistributor::updateLocalGlobalNumbering()");

      Timer t;
#if (DEBUG != 0)
      bool printInfo = true;
#else
      bool printInfo = false;
#endif
      Parameters::get("parallel->print dofmap info", printInfo);

      for (size_t i = 0; i < dofMaps.size(); i++)
      {
        dofMaps[i]->update();
        if (printInfo)
          dofMaps[i]->printInfo();
      }

      // === Create periodic DOF maps, if there are periodic boundaries.  ===

      if (hasPeriodicBoundary)
      {
        createPeriodicMap();

        for (size_t i = 0; i < dofMaps.size(); i++)
          dofMaps[i]->updateMatIndex();
      }

      if (printInfo)
      {
        int test = 0;
        Parameters::get("parallel->remove periodic boundary", test);

        MSG("------------- Debug information -------------\n");
        for (size_t k = 0; k < meshes.size(); k++)
        {
          Mesh* mesh = meshes[k];

          MSG("|  mesh pointer adr:         %p\n", mesh);
          MSG("|  number of levels:         %d\n", levelData.getNumberOfLevels());
          MSG("|  number of FE spaces:      %d\n", meshToFeSpaces[mesh].size());


          for (size_t i = 0; i < dofMaps.size(); i++)
          {
            vector<const FiniteElemSpace*> dofMapSpaces = dofMaps[i]->getFeSpaces(mesh);

            if (dofMapSpaces.empty())
              continue;

            for (size_t j = 0; j < dofMapSpaces.size(); j++)
            {
              const FiniteElemSpace* feSpace = dofMapSpaces[j];
              MSG("|  FE space = %d  (pointer adr %p):\n", j, feSpace);
              MSG("|      nRankDofs    = %d\n", (*(dofMaps[i]))[feSpace].nRankDofs);
              MSG("|      nOverallDofs = %d\n", (*(dofMaps[i]))[feSpace].nOverallDofs);
              MSG("|      rStartDofs   = %d\n", (*(dofMaps[i]))[feSpace].rStartDofs);
            }
          }

          //        ParallelDebug::writeDebugFile(meshToFeSpaces[mesh][meshToFeSpaces[mesh].size() - 1],
          //   				  *(dofMaps[0]),
          //   				  debugOutputDir + "mpi-dbg", "dat");

          if (test == 0)
            ParallelDebug::testGlobalIndexByCoords(*this, mesh);
        }

      }
      else
      {

        //       int tmp = 0;
        //       Parameters::get(name + "->write parallel debug file", tmp);

        for (size_t k = 0; k < meshes.size(); k++)
        {
          Mesh* mesh = meshes[k];

          for (size_t i = 0; i < dofMaps.size(); i++)
          {

            vector<const FiniteElemSpace*> dofMapSpaces = dofMaps[i]->getFeSpaces(mesh);

            if (dofMapSpaces.empty())
              continue;

            for (size_t j = 0; j < dofMapSpaces.size(); j++)
              MSG("FE space %d (mesh ptr %p):  nRankDofs    = %d   nOverallDofs = %d\n", j, dofMapSpaces[j]->getMesh(),
                  (*(dofMaps[i]))[dofMapSpaces[j]].nRankDofs,
                  (*(dofMaps[i]))[dofMapSpaces[j]].nOverallDofs);
          }

          // Write debug file for each mesh and each dofmap?
          //      if (tmp)
          //        ParallelDebug::writeDebugFile(meshToFeSpaces[mesh][meshToFeSpaces[mesh].size() - 1],
          //				      *(dofMaps[0]),
          //				      debugOutputDir + "mpi-dbg", "dat");
        }
      }

      MPI::COMM_WORLD.Barrier();
      MSG("Update dof mapping needed %.5f seconds\n", t.elapsed());
    }


    void MeshDistributor::updateDofsToDofMapping(Mesh* mesh)
    {
      FUNCNAME("MeshDistributor::updateDofsToDofMapping()");

      Timer* t = NULL;
      if (printTimings)
        t = new Timer();

      TEST_EXIT(dofMaps.size())("No DOF mapping defined!\n");

      for (size_t i = 0; i < dofMaps.size(); i++)
      {
        vector<const FiniteElemSpace*> dofMapSpaces = dofMaps[i]->getFeSpaces(mesh);

        if (dofMapSpaces.empty())
          continue;

        dofMaps[i]->clear(mesh);

        for (size_t j = 0; j < dofMapSpaces.size(); j++)
          updateDofsToDofMapping(*(dofMaps[i]), dofMapSpaces[j]);

      }

      if (printTimings && t)
      {
        MSG("MESHDIST 02 (update DOF mappings) needed %.5f seconds\n", t->elapsed());
        delete t;
      }
    }


    void MeshDistributor::updateDofsToDofMapping(ParallelDofMapping& dofMap,
        const FiniteElemSpace* feSpace)
    {

      DofComm& dcom = dofMap[feSpace].getDofComm();
      Mesh* mesh = feSpace->getMesh();


      // === Get all DOFs in ranks partition. ===

      std::set<const DegreeOfFreedom*> rankDofSet;
      mesh->getAllDofs(feSpace, rankDofSet);

      DofContainer rankDofs(rankDofSet.begin(), rankDofSet.end());
      sort(rankDofs.begin(), rankDofs.end(), cmpDofsByValue);


      // === Traverse interior boundaries and get all DOFs on them. ===

      DofContainerSet nonRankDofs;
      for (DofComm::Iterator it(dcom.getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
        for (; !it.endDofIter(); it.nextDof())
          nonRankDofs.insert(it.getDof());

      for (DofContainer::iterator it = rankDofs.begin();
           it != rankDofs.end(); ++it)
      {
        if (nonRankDofs.count(*it))
          dofMap[feSpace].insertNonRankDof(**it);
        else
          dofMap[feSpace].insertRankDof(**it);
      }
    }


    void MeshDistributor::createPeriodicMap()
    {
      FUNCNAME("MeshDistributor::createPeriodicMap()");

      // Clear all periodic DOF mappings calculated before. We do it from scratch.
      periodicMap.clear();

      // If there are no periodic boundaries, return.
      if (!intBoundary[0].hasPeriodic())
        return;

      TEST_EXIT(levelData.getNumberOfLevels() == 1)
      ("Periodic DOF map does not support multi level domain decomposition!\n");

      //    MPI::COMM_WORLD.Barrier();   [TODO: CHANGE BECAUSE NOT ALL RANKS HAVE PERIODIC MAP!!!]
      double first = MPI::Wtime();

      for (size_t i = 0; i < feSpaces.size(); i++)
        createPeriodicMap(feSpaces[i]);

      //    MPI::COMM_WORLD.Barrier();
      MSG("Creation of periodic mapping needed %.5f seconds\n",  MPI::Wtime() - first);
    }


    void MeshDistributor::createPeriodicMap(const FiniteElemSpace* feSpace)
    {
      FUNCNAME("MeshDistributor::createPeriodicMap()");

      TEST_EXIT(dofMaps.size())("No DOF mapping defined!\n");

      ParallelDofMapping* parDofMap = NULL;
      for (size_t i = 0; i < dofMaps.size(); i++)
        if (dofMaps[i]->hasFeSpace(feSpace))
        {
          parDofMap = dofMaps[i];
          break;
        }

      TEST_EXIT_DBG(parDofMap)("FE space not found in parallel dof map.\n");
      ComponentDofMap& dofMap = (*parDofMap)[feSpace];
      Mesh* mesh = feSpace->getMesh();
      DofComm::DataType& periodicDofs = dofComms[mesh][0].getPeriodicDofs();

      StdMpi<vector<int>> stdMpi(MPI::COMM_WORLD, false);


      // === Each rank traverse its periodic boundaries and sends the DOF      ===
      // === indices to the rank "on the other side" of the periodic boundary. ===

      map<int, vector<int>> rankToDofType;

      for (RankToBoundMap::iterator it = intBoundary[0].getPeriodic().begin();
           it != intBoundary[0].getPeriodic().end(); ++it)
      {

        if (it->first == mpiRank)
        {
          // Here we have a periodic boundary within rank's subdomain. So we can
          // compute the periodic mappings without communication.

          TEST_EXIT_DBG(it->second.size() % 2 == 0)("Should not happen!\n");

          for (size_t i = 0; i < it->second.size(); i++)
          {
            AtomicBoundary& bound = it->second[i];

            DofContainer dofs0, dofs1;

            intBoundary[0].getElementPtr(bound.rankObj.elIndex, mesh)
            ->getAllDofs(feSpace, bound.rankObj, dofs0);
            intBoundary[0].getElementPtr(bound.neighObj.elIndex, mesh)
            ->getAllDofs(feSpace, bound.neighObj, dofs1);

            TEST_EXIT_DBG(dofs0.size() == dofs1.size())
            ("Number of DOFs does not fit together: %d %d\n",
             dofs0.size(), dofs1.size());

            BoundaryType type = bound.type;

            for (size_t j = 0; j < dofs0.size(); j++)
            {
              DegreeOfFreedom globalDof0 = dofMap[*(dofs0[j])].global;
              DegreeOfFreedom globalDof1 = dofMap[*(dofs1[j])].global;

              if (!periodicMap.isPeriodicOnBound(feSpace, type, globalDof0))
                periodicMap.add(feSpace, type, globalDof0, globalDof1);
            }
          }

        }
        else
        {
          // Here we have a periodic boundary between two ranks.

          // Create DOF indices on the boundary.
          DofContainer& dofs = periodicDofs[it->first][feSpace];
          for (vector<AtomicBoundary>::iterator boundIt = it->second.begin();
               boundIt != it->second.end(); ++boundIt)
          {

            int nDofs = dofs.size();
            intBoundary[0].getElementPtr(boundIt->rankObj.elIndex, mesh)
            ->getAllDofs(feSpace, boundIt->rankObj, dofs);

            for (unsigned int i = 0; i < (dofs.size() - nDofs); i++)
              rankToDofType[it->first].push_back(boundIt->type);
          }

          // Send the global indices to the rank on the other side.
          stdMpi.getSendData(it->first).reserve(dofs.size());
          for (unsigned int i = 0; i < dofs.size(); i++)
            stdMpi.getSendData(it->first).push_back(dofMap[*(dofs[i])].global);

          // Receive from this rank the same number of dofs.
          stdMpi.recv(it->first, dofs.size());
        }
      }

      stdMpi.updateSendDataSize();
      stdMpi.startCommunication();


      // === The rank has received the DOFs from the rank on the other side of ===
      // === the boundary. Now it can use them to create the mapping between   ===
      // === the periodic DOFs in this rank and the corresponding periodic     ===
      // === DOFs from the other ranks.                                        ===


      for (RankToBoundMap::iterator it = intBoundary[0].getPeriodic().begin();
           it != intBoundary[0].getPeriodic().end(); ++it)
      {
        DofContainer& dofs = periodicDofs[it->first][feSpace];
        vector<int>& types = rankToDofType[it->first];

        TEST_EXIT_DBG(dofs.size() == types.size())("Should not happen!\n");

        // Added the received DOFs to the mapping.
        for (unsigned int i = 0; i < dofs.size(); i++)
        {
          int globalDofIndex = dofMap[*(dofs[i])].global;
          int mapGlobalDofIndex = stdMpi.getRecvData(it->first)[i];
          BoundaryType type = types[i];

          // Check if this global DOF with the corresponding boundary type was
          // not added before by another periodic boundary from other rank.
          if (!periodicMap.isPeriodicOnBound(feSpace, type, globalDofIndex))
            periodicMap.add(feSpace, type, globalDofIndex, mapGlobalDofIndex);
        }
      }


      StdMpi<PeriodicDofMap> stdMpi2(MPI::COMM_WORLD);


      for (RankToBoundMap::iterator it = intBoundary[0].getPeriodic().begin();
           it != intBoundary[0].getPeriodic().end(); ++it)
      {
        if (it->first == mpiRank)
          continue;

        PeriodicDofMap perObjMap;

        for (vector<AtomicBoundary>::iterator boundIt = it->second.begin();
             boundIt != it->second.end(); ++boundIt)
        {
          // If the boundary type is a normal periodic BC, continue. We just
          // operate on indirectly connected periodic boundaries.
          if (elObjDb.isValidPeriodicType(boundIt->type))
            continue;

          DofContainer dofs;
          intBoundary[0].getElementPtr(boundIt->rankObj.elIndex, mesh)
          ->getAllDofs(feSpace, boundIt->rankObj, dofs);

          for (size_t i = 0; i < dofs.size(); i++)
          {
            DegreeOfFreedom globalDof = dofMap[*dofs[i]].global;

            std::set<BoundaryType>& assoc =
              periodicMap.getAssociations(feSpace, globalDof);
            TEST_EXIT_DBG(assoc.size() > 0)("Should not happen!\n");

            for (std::set<BoundaryType>::iterator perAscIt = assoc.begin();
                 perAscIt != assoc.end(); ++perAscIt)
              if (elObjDb.isValidPeriodicType(*perAscIt))
                perObjMap[*perAscIt][globalDof] =
                  periodicMap.map(feSpace, *perAscIt, globalDof);
          }
        }

        stdMpi2.send(it->first, perObjMap);
        stdMpi2.recv(it->first);
      }

      stdMpi2.startCommunication();

      // NOTE: Before changing the data structure of periodic boundary DOFS,
      //       at this place no periodicDofAssociations are set for the received
      //       DOFs. I'm note sure what is correct here.
      for (map<int, PeriodicDofMap>::iterator it = stdMpi2.getRecvData().begin();
           it != stdMpi2.getRecvData().end(); ++it)
        periodicMap.add(feSpace, it->second);
    }


    void MeshDistributor::createMacroElementInfo()
    {
      // Store neighbours info.
      for (deque<MacroElement*>::iterator it = macroMesh->firstMacroElement();
           it != macroMesh->endOfMacroElements(); ++it)
      {

        int elIndex = (*it)->getIndex();

        macroElementNeighbours[elIndex].resize(macroMesh->getGeo(NEIGH));
        for (int i = 0; i < macroMesh->getGeo(NEIGH); i++)
          macroElementNeighbours[elIndex][i] =
            (*it)->getNeighbour(i) ? (*it)->getNeighbour(i)->getIndex() : -1;
      }

      // Store macro elements
      for (size_t i = 0; i < meshes.size(); i++)
        for (deque<MacroElement*>::iterator it = meshes[i]->firstMacroElement();
             it != meshes[i]->endOfMacroElements(); ++it)
          allMacroElements[meshes[i]].push_back(*it);

    }

    void MeshDistributor::updateMacroElementInfo()
    {
      for (size_t i = 0; i < meshes.size(); i++)
      {
        deque<MacroElement*>& meshMacros = meshes[i]->getMacroElements();

        vector<MacroElement*>& allMacros = allMacroElements[meshes[i]];

        for (size_t j = 0; j < allMacros.size(); j++)
        {
          for (deque<MacroElement*>::iterator it = meshMacros.begin();
               it != meshMacros.end(); ++it)
          {
            if ((*it)->getIndex() == allMacros[j]->getIndex() &&
                *it != allMacros[j])
            {
              delete allMacros[j];
              allMacros[j] = *it;
            }
          }
        }
      }
    }

  }
}
