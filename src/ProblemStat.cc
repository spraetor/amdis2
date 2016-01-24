#include <sstream>
#include <string>

#include "ProblemStat.h"
#include "Operator.h"
#include "SystemVector.h"
#include "DOFMatrix.h"
#include "FiniteElemSpace.h"
#include "Marker.h"
#include "AdaptInfo.h"
#include "io/FileWriter.h"
#include "CoarseningManager.h"
#include "RefinementManager.h"
// #include "DualTraverse.h"
#include "Mesh.h"
#include "solver/LinearSolverInterface.h"
#include "DirichletBC.h"
#include "RobinBC.h"
#include "PeriodicBC.h"
#include "Lagrange.h"
#include "Flag.h"
#include "Timer.h"
#include "est/Estimator.h"
#include "io/VtkWriter.h"
#include "io/ValueReader.h"
#include "ProblemStatDbg.h"
// #include "Debug.h"

namespace AMDiS
{
  using namespace std;

  ProblemStatSeq::ProblemStatSeq(string nameStr)
    : name(nameStr),
      nComponents(-1),
      nAddComponents(0),
      nMeshes(0),
      traverseInfo(0),
      solver(NULL),
      solution(NULL),
      rhs(NULL),
      systemMatrix(NULL),
      refinementManager(NULL),
      coarseningManager(NULL),
      info(10),
      boundaryConditionSet(false),
      writeAsmInfo(false),
      solutionTime(0.0),
      buildTime(0.0)
  {
    Parameters::get(name + "->components", nComponents);
    Parameters::get(name + "->additional components", nAddComponents);
    TEST_EXIT(nComponents > 0)("No value set for parameter \"%s->components\"!\n",
                               name.c_str());
    TEST_EXIT(nAddComponents >= 0)("Wrong parameter \"%s->additional components\"!\n",
                                   name.c_str());
    estimator.resize(nComponents, NULL);
    marker.resize(nComponents, NULL);

    assembleMatrixOnlyOnce.resize(nComponents);
    assembledMatrix.resize(nComponents);
    for (int i = 0; i < nComponents; i++)
    {
      assembleMatrixOnlyOnce[i].resize(nComponents);
      assembledMatrix[i].resize(nComponents);
      for (int j = 0; j < nComponents; j++)
      {
        assembleMatrixOnlyOnce[i][j] = false;
        assembledMatrix[i][j] = false;
      }
    }


    // === Initialize name of components. ===

    componentNames.resize(nComponents, "");
    for (int i = 0; i < nComponents; i++)
      componentNames[i] = "solution[" + std::to_string(i) + "]";

    // Parameters::get(name + "->names", componentNames);
    // componentNames.resize(nComponents);
    // TODO: need to be corrected

    for (int i = 0; i < nComponents; i++)
      Parameters::get(name + "->name[" + std::to_string(i) + "]",
                      componentNames[i]);

    Parameters::get("debug->write asm info", writeAsmInfo);
  }


  ProblemStatSeq::~ProblemStatSeq()
  {
    if (rhs) {
      delete rhs;
      rhs = NULL;
    }
    if (solution) {
      delete solution;
      solution = NULL;
    }

    if (solver) {
      delete solver;
      solver = NULL;
    }

    if (systemMatrix)
    {
      for (int i = 0; i < nComponents; i++)
        for (int j = 0; j < nComponents; j++)
          if ((*systemMatrix)[i][j])
          {
            delete (*systemMatrix)[i][j];
            (*systemMatrix)[i][j] = NULL;
          }

      delete systemMatrix;
      systemMatrix = NULL;
    }

    for (unsigned int i = 0; i < meshes.size(); i++)
      if (meshes[i])
      {
	// TODO: Why are meshes not deleted????
        // 	delete meshes[i];
        // 	meshes[i] = NULL;
      }

    for (unsigned int i = 0; i < estimator.size(); i++)
      if (estimator[i])
      {
        delete estimator[i];
        estimator[i] = NULL;
      }

    for (unsigned int i = 0; i < marker.size(); i++)
      if (marker[i])
      {
        delete marker[i];
        marker[i] = NULL;
      }
  }

  void ProblemStatSeq::initialize(Flag initFlag,
                                  ProblemStatSeq* adoptProblem,
                                  Flag adoptFlag)
  {
    FUNCNAME("ProblemStat::initialize()");

    // === create meshes ===
    if (meshes.size() != 0)
    {
      WARNING("meshes already created\n");
    }
    else
    {
      if (initFlag.isSet(CREATE_MESH) ||
          (!adoptFlag.isSet(INIT_MESH) &&
           (initFlag.isSet(INIT_SYSTEM) || initFlag.isSet(INIT_FE_SPACE))))
        createMesh();

      if (adoptProblem &&
          (adoptFlag.isSet(INIT_MESH) ||
           adoptFlag.isSet(INIT_SYSTEM) ||
           adoptFlag.isSet(INIT_FE_SPACE)))
      {
        meshes = adoptProblem->getMeshes();
        if (meshes.size() == 1)
          componentMeshes.resize(nComponents, meshes[0]);
        else if (adoptProblem->getNumComponents() >= nComponents)
        {
          componentMeshes.resize(nComponents);
          std::copy(adoptProblem->componentMeshes.begin(),
                    adoptProblem->componentMeshes.begin() + nComponents,
                    componentMeshes.begin());
        }
        else
        {
          componentMeshes.resize(nComponents, meshes[0]);
          WARNING("componentMeshes may not be derived correctly from the adoptProblem. You have to do this manually!\n");
        }

        if (nAddComponents > 0)
        {
          WARNING("Additional meshed can not be adopted from adoptProblem. You have to initialize these meshes manually!\n");
        }
      }
    }

    if (meshes.size() == 0)
      WARNING("no mesh created\n");

    // === create refinement/corasening-manager ===
    if (refinementManager != NULL && coarseningManager != NULL)
    {
      WARNING("refinement-/coarseningmanager already created\n");
    }
    else
    {
      if (initFlag.isSet(CREATE_MESH) ||
          (!adoptFlag.isSet(INIT_MESH) &&
           (initFlag.isSet(INIT_SYSTEM) || initFlag.isSet(INIT_FE_SPACE))))
        createRefCoarseManager();

      if (adoptProblem &&
          (adoptFlag.isSet(INIT_MESH) ||
           adoptFlag.isSet(INIT_SYSTEM) ||
           adoptFlag.isSet(INIT_FE_SPACE)))
      {
        refinementManager = adoptProblem->refinementManager;
        coarseningManager = adoptProblem->coarseningManager;
      }
    }

    if (refinementManager == NULL || coarseningManager == NULL)
      WARNING("no refinement-/coarseningmanager created\n");

    // === create fespace ===
    if (feSpaces.size() != 0)
    {
      WARNING("feSpaces already created\n");
    }
    else
    {
      if (initFlag.isSet(INIT_FE_SPACE) ||
          (initFlag.isSet(INIT_SYSTEM) && !adoptFlag.isSet(INIT_FE_SPACE)))
        createFeSpace(NULL);

      if (adoptProblem &&
          (adoptFlag.isSet(INIT_FE_SPACE) || adoptFlag.isSet(INIT_SYSTEM)))
      {
        feSpaces = adoptProblem->getFeSpaces();
        traverseInfo = adoptProblem->traverseInfo;

        if (feSpaces.size() == 1)
          componentSpaces.resize(nComponents, feSpaces[0]);
        else if (adoptProblem->getNumComponents() >= nComponents)
        {
          componentSpaces.resize(nComponents);
          std::copy(adoptProblem->componentSpaces.begin(),
                    adoptProblem->componentSpaces.begin() + nComponents,
                    componentSpaces.begin());
        }
        else
        {
          componentSpaces.resize(nComponents, feSpaces[0]);
          WARNING("componentSpaces may not be derived correctly from the adoptProblem. You have to do this manually!\n");
        }

        if (nAddComponents > 0)
        {
          WARNING("Additional feSpaces can not be adopted from adoptProblem. You have to initialize these fespaces manually!\n");
        }
      }
    }

    if (feSpaces.size() == 0)
      WARNING("no feSpace created\n");

    // === create system ===
    if (initFlag.isSet(INIT_SYSTEM))
      createMatricesAndVectors();

    if (adoptProblem && adoptFlag.isSet(INIT_SYSTEM))
    {
      solution = adoptProblem->getSolution();
      rhs = adoptProblem->getRhs();
      systemMatrix = adoptProblem->getSystemMatrix();
    }

    // === create solver ===
    if (solver)
    {
      WARNING("solver already created\n");
    }
    else
    {
      if (initFlag.isSet(INIT_SOLVER))
        createSolver();

      if (adoptProblem && adoptFlag.isSet(INIT_SOLVER))
      {
        TEST_EXIT(!solver)("solver already created\n");
        solver = adoptProblem->getSolver();
      }
    }

    if (!solver)
      WARNING("no solver created\n");

    // === create estimator ===
    if (initFlag.isSet(INIT_ESTIMATOR))
      createEstimator();

    if (adoptProblem && adoptFlag.isSet(INIT_ESTIMATOR))
      estimator = adoptProblem->getEstimators();

    // === create marker ===
    if (initFlag.isSet(INIT_MARKER))
      createMarker();

    if (adoptProblem && adoptFlag.isSet(INIT_MARKER))
      marker = adoptProblem->getMarkers();

    // === create file writer ===
    if (initFlag.isSet(INIT_FILEWRITER))
      createFileWriter();

    int globalRefinements = 0;

    // If AMDiS is compiled for parallel computations, the global refinements are
    // ignored here. Later, each rank will add the global refinements to its
    // private mesh.
#ifndef HAVE_PARALLEL_DOMAIN_AMDIS
    Parameters::get(meshes[0]->getName() + "->global refinements",
                    globalRefinements);
#endif

    bool initMesh = initFlag.isSet(INIT_MESH);

    // Initialize the meshes
    for (int i = 0; i < static_cast<int>(meshes.size()); i++)
      if (initMesh && meshes[i] && !(meshes[i]->isInitialized()))
        meshes[i]->initialize();

    // === read value file and use it for the mesh values ===
    string valueFilename("");
    Parameters::get(meshes[0]->getName() + "->value file name", valueFilename);
    if (valueFilename.length())
      io::ValueReader::readValue(valueFilename,
                                 meshes[0],
                                 solution->getDOFVector(0),
                                 meshes[0]->getMacroFileInfo());

    // === do global refinements ===

    for (unsigned int i = 0; i < meshes.size(); i++)
      if (initMesh && meshes[i])
        refinementManager->globalRefine(meshes[i], globalRefinements);
  }


  void ProblemStatSeq::createMesh()
  {
    FUNCNAME("ProblemStat::createMesh()");

    map<int, Mesh*> meshForRefinementSet;

    string meshName("");
    Parameters::get(name + "->mesh", meshName);
    TEST_EXIT(meshName != "")("No mesh name specified for \"%s->mesh\"!\n",
                              name.c_str());

    int dim = 0;
    Parameters::get(name + "->dim", dim);
    TEST_EXIT(dim)("No problem dimension specified for \"%s->dim\"!\n",
                   name.c_str());

    componentMeshes.resize(nComponents + nAddComponents);

    for (int i = 0; i < nComponents + nAddComponents; i++)
    {
      int refSet = -1;
      Parameters::get(name + "->refinement set[" + std::to_string(i) + "]",
                      refSet);
      if (refSet < 0)
        refSet = 0;

      if (meshForRefinementSet[refSet] == NULL)
      {
        Mesh* newMesh = new Mesh(meshName, dim);
        meshForRefinementSet[refSet] = newMesh;
        meshes.push_back(newMesh);
        nMeshes++;
      }
      componentMeshes[i] = meshForRefinementSet[refSet];
    }
  }


  void ProblemStatSeq::createRefCoarseManager()
  {
    FUNCNAME("ProblemStat::createRefCoarseManager()");

    int dim = Parameters::get<int>(name + "->dim", 0);
    TEST_EXIT(dim)("No problem dimension specified for \"%s->dim\"!\n",
                   name.c_str());

    switch (dim)
    {
    case 1:
      refinementManager = new RefinementManager1d();
      coarseningManager = new CoarseningManager1d();
      break;
    case 2:
      refinementManager = new RefinementManager2d();
      coarseningManager = new CoarseningManager2d();
      break;
    case 3:
      refinementManager = new RefinementManager3d();
      coarseningManager = new CoarseningManager3d();
      break;
    default:
      ERROR_EXIT("invalid dim!\n");
    }
  }


  void ProblemStatSeq::createFeSpace(DOFAdmin* admin)
  {
    FUNCNAME("ProblemStat::createFeSpace()");

    map<pair<Mesh*, string>, FiniteElemSpace*> feSpaceMap;
    int dim = Parameters::get<int>(name + "->dim", -1);
    TEST_EXIT(dim != -1)("no problem dimension specified!\n");

    componentSpaces.resize(nComponents + nAddComponents, NULL);
    traverseInfo.resize(nComponents);

    for (int i = 0; i < nComponents + nAddComponents; i++)
    {
      if (componentSpaces[i] != NULL)
      {
        WARNING("feSpace already created\n");
        continue;
      }
      string componentString = "[" + std::to_string(i) + "]";

      string feSpaceName = "";
      string initFileStr = name + "->feSpace" + componentString;
      Parameters::get(initFileStr, feSpaceName);

      // synonym for "feSpace"
      if (feSpaceName.size() == 0)
      {
        initFileStr = name + "->finite element space" + componentString;
        Parameters::get(initFileStr, feSpaceName);
      }

      // for backward compatibility also accept the old syntax
      if (feSpaceName.size() == 0)
      {
        initFileStr = name + "->polynomial degree" + componentString;
        int degree = Parameters::get<int>(initFileStr, -1);
        TEST_EXIT(degree > 0)
        ("Poynomial degree in component %d must be larger than zero!\n", i);

        feSpaceName = "Lagrange" + std::to_string(degree);
      }

      if (feSpaceMap[pair<Mesh*, string>(componentMeshes[i], feSpaceName)] == NULL)
      {
        BasisFunctionCreator* basisFctCreator =
          dynamic_cast<BasisFunctionCreator*>(CreatorMap<BasisFunction>::getCreator(feSpaceName, initFileStr));
        TEST_EXIT(basisFctCreator)("No valid basisfunction type found in parameter \"%s\"\n", 
				   initFileStr.c_str());
        basisFctCreator->setDim(dim);

        FiniteElemSpace* newFeSpace =
          FiniteElemSpace::provideFeSpace(admin, basisFctCreator->create(),
                                          componentMeshes[i], "FeSpace" + componentString + " (" + feSpaceName + ")");

        feSpaceMap[pair<Mesh*, string>(componentMeshes[i], feSpaceName)] = newFeSpace;
        feSpaces.push_back(newFeSpace);
      }

      componentSpaces[i] = feSpaceMap[pair<Mesh*, string>(componentMeshes[i], feSpaceName)];
    }

    for (int i = 0; i < nComponents; i++)
    {
      for (int j = 0; j < nComponents; j++)
        traverseInfo.getMatrix(i, j).setFeSpace(componentSpaces[i], componentSpaces[j]);

      traverseInfo.getVector(i).setFeSpace(componentSpaces[i]);
    }

    // create dof admin for vertex dofs if neccessary
    for (int i = 0; i < static_cast<int>(meshes.size()); i++)
    {
      if (meshes[i]->getNumberOfDofs(VERTEX) == 0)
      {
        DimVec<int> ln_dof(meshes[i]->getDim(), 0);
        ln_dof[VERTEX] = 1;
        meshes[i]->createDOFAdmin("vertex dofs", ln_dof);
      }
    }
  }


  void ProblemStatSeq::createMatricesAndVectors()
  {
    // === create vectors and system matrix ===

    systemMatrix = new Matrix<DOFMatrix*>(nComponents, nComponents);
    systemMatrix->set((DOFMatrix*)(NULL));
    rhs = new SystemVector("rhs", componentSpaces, nComponents);
    solution = new SystemVector("solution", componentSpaces, nComponents);

    for (int i = 0; i < nComponents; i++)
    {
      (*systemMatrix)[i][i] = new DOFMatrix(componentSpaces[i],
                                            componentSpaces[i], "A_ii");
      (*systemMatrix)[i][i]->setCoupleMatrix(false);

      //set parallel synchronization later in ParalleProblemStat
      rhs->setDOFVector(i, new DOFVector<double>(componentSpaces[i],
                        "rhs[" + std::to_string(i) + "]", false));

      //set parallel synchronization later in ParalleProblemStat
      solution->setDOFVector(i, new DOFVector<double>(componentSpaces[i],
                             componentNames[i], false));
      solution->getDOFVector(i)->setCoarsenOperation(COARSE_INTERPOL);
      solution->getDOFVector(i)->set(0.0);
    }
  }


  void ProblemStatSeq::createSolver()
  {
    FUNCNAME("ProblemStat::createSolver()");

    // definition of standard-backends
#if defined HAVE_PARALLEL_PETSC
    string backend("p_petsc");
#elif defined HAVE_PARALLEL_MTL
    string backend("p_mtl");
#elif defined HAVE_SEQ_PETSC
    string backend("petsc");
#else
    string backend("mtl");
#endif

    // === read backend-name ===
    string initFileStr = name + "->solver";
    Parameters::get(initFileStr + "->backend", backend);

    // === read solver-name ===
    string solverType = Parameters::get<string>(initFileStr, "0");

    if (backend != "0" && backend != "no" && backend != "")
      solverType = backend + "_" + solverType;

    // === create solver ===
    LinearSolverCreator* solverCreator =
      dynamic_cast<LinearSolverCreator*>(CreatorMap<LinearSolverInterface>::getCreator(solverType, initFileStr));
    TEST_EXIT(solverCreator)("No valid solver type found in parameter \"%s\"\n", 
			     initFileStr.c_str());
    solverCreator->setName(initFileStr);
    solver = solverCreator->create();
  }


  void ProblemStatSeq::createEstimator()
  {
    FUNCNAME("ProblemStat::createEstimator()");

    // create and set leaf data prototype
    for (Mesh* mesh : meshes)
      mesh->setElementDataPrototype(new LeafDataEstimatableVec(new LeafDataCoarsenableVec));

    for (int i = 0; i < nComponents; i++)
    {
      TEST_EXIT(estimator[i] == NULL)("estimator already created\n");
      string estName = name + "->estimator[" + std::to_string(i) + "]";

      // === create estimator ===
      string estimatorType = Parameters::get<string>(estName, "0");
      if (estimatorType == "0")
        Parameters::get(estName + "->name", estimatorType);

      EstimatorCreator* estimatorCreator =
        dynamic_cast<EstimatorCreator*>(CreatorMap<Estimator>::getCreator(estimatorType, estName));
      if (estimatorCreator)
      {
        estimatorCreator->setName(estName);
        estimatorCreator->setRow(i);
        estimatorCreator->setSolution(solution->getDOFVector(i));
        estimator[i] = estimatorCreator->create();
      }


      if (estimator[i])
        for (int j = 0; j < nComponents; j++)
          estimator[i]->addSystem((*systemMatrix)[i][j],
                                  solution->getDOFVector(j),
                                  rhs->getDOFVector(i));
    }
  }


  void ProblemStatSeq::createMarker()
  {
    int nMarkersCreated = 0;

    for (int i = 0; i < nComponents; i++)
    {
      marker[i] = Marker::createMarker(name + "->marker[" + std::to_string(i) + "]", i);

      if (marker[i])
      {
        nMarkersCreated++;

        // If there is more than one marker, and all components are defined
        // on the same mesh, the maximum marking has to be enabled.
        if (nMarkersCreated > 1 && nMeshes == 1)
          marker[i]->setMaximumMarking(true);
      }
    }
  }


  void ProblemStatSeq::createFileWriter()
  {
    FUNCNAME("ProblemStat::createFileWriter()");

    // Create one filewriter for all components of the problem
    string numberedName  = name + "->output";
    string filename = Parameters::get<string>(numberedName + "->filename", "");

    if (filename != "")
    {
      vector<DOFVector<double>*> solutionList(nComponents);

      for (int i = 0; i < nComponents; i++)
      {
        TEST_EXIT(componentMeshes[0] == componentMeshes[i])
        ("All Meshes have to be equal to write a vector file.\n");

        solutionList[i] = solution->getDOFVector(i);
      }

      fileWriters.push_back(new FileWriter(numberedName,
                                           componentMeshes[0],
                                           solutionList));
    }

    // Create own filewriters for each component of the problem
    for (int i = 0; i < nComponents; i++)
    {
      numberedName = name + "->output[" + std::to_string(i) + "]";
      filename = "";
      Parameters::get(numberedName + "->filename", filename);

      if (filename != "")
        fileWriters.push_back(new FileWriter(numberedName,
                                             componentMeshes[i],
                                             solution->getDOFVector(i)));
    }

    // create a filewrite for groups of components to write vector-valued output
    int nVectors = Parameters::get<int>(name + "->output->num vectors", 0);
    if (nVectors > 0)
    {
      for (int j = 0; j < nVectors; j++)
      {
        numberedName = name + "->output->vector[" + std::to_string(j) + "]";

        filename = Parameters::get<string>(numberedName + "->filename", "");
        string componentName = Parameters::get<string>(numberedName + "->name", "");
        vector<string> names;
        if (componentName != "")
          names.push_back(componentName);
        vector<int> comp;
        Parameters::get(numberedName + "->components", comp);

        if (filename != "" && comp.size() > 0)
        {
          // Create own filewriters for each component of the problem
          vector<DOFVector<double>*> vectors;
          for (size_t i = 0; i < comp.size(); i++)
            vectors.push_back(solution->getDOFVector(comp[i]));

          fileWriters.push_back(new FileWriter(numberedName,
                                               componentMeshes[comp[0]],
                                               vectors, names
                                              ));

        }
      }
    }
  }


  void ProblemStatSeq::solve(AdaptInfo& adaptInfo,
                             bool createMatrixData,
                             bool storeMatrixData)
  {
    FUNCNAME("Problem::solve()");

    if (!solver)
    {
      WARNING("no solver\n");
      return;
    }

    Timer t;
    solver->solveSystem(solverMatrix, *solution, *rhs,
                        createMatrixData, storeMatrixData);

    INFO(info, 8)("solution of discrete system needed %.5f seconds\n",
                  t.elapsed());
    solutionTime = t.elapsed();

    adaptInfo.setSolverIterations(solver->getIterations());
    adaptInfo.setMaxSolverIterations(solver->getMaxIterations());
    adaptInfo.setSolverTolerance(solver->getTolerance());
    adaptInfo.setSolverResidual(solver->getResidual());
  }


  void ProblemStatSeq::estimate(AdaptInfo& adaptInfo)
  {
    FUNCNAME("ProblemStat::estimate()");

    Timer t;

    for (int i = 0; i < nComponents; i++)
    {
      Estimator* scalEstimator = estimator[i];

      if (scalEstimator)
      {
        traverseInfo.updateStatus();
        scalEstimator->setTraverseInfo(traverseInfo);
        scalEstimator->estimate(adaptInfo.getTimestep());

        adaptInfo.setEstSum(scalEstimator->getErrorSum(), i);
        adaptInfo.setEstMax(scalEstimator->getErrorMax(), i);

        if (adaptInfo.getRosenbrockMode() == false) // TODO: remove dependence on rosenbrock method
        {
          adaptInfo.setTimeEstSum(scalEstimator->getTimeEst(), i);
          adaptInfo.setTimeEstMax(scalEstimator->getTimeEstMax(), i);
        }
      }
    }


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    MPI::COMM_WORLD.Barrier();
#endif
    INFO(info, 8)("estimation of the error needed %.5f seconds\n",
                  t.elapsed());
  }


  Flag ProblemStatSeq::markElements(AdaptInfo& adaptInfo)
  {
    FUNCNAME_DBG("ProblemStat::markElements()");

    TEST_EXIT_DBG(static_cast<unsigned int>(nComponents) == marker.size())
    ("Wrong number of markers!\n");

    Flag markFlag = 0;
    for (int i = 0; i < nComponents; i++)
      if (marker[i])
        markFlag |= marker[i]->markMesh(adaptInfo, componentMeshes[i]);

    return markFlag;
  }


  Flag ProblemStatSeq::refineMesh(AdaptInfo& adaptInfo)
  {
    int nMeshes = static_cast<int>(meshes.size());
    Flag refineFlag = 0;
    for (int i = 0; i < nMeshes; i++)
      if (i >= nComponents || adaptInfo.isRefinementAllowed(i))
        refineFlag |= refinementManager->refineMesh(meshes[i]);

    return refineFlag;
  }


  Flag ProblemStatSeq::coarsenMesh(AdaptInfo& adaptInfo)
  {
    int nMeshes = static_cast<int>(meshes.size());
    Flag coarsenFlag = 0;
    for (int i = 0; i < nMeshes; i++)
      if (i >= nComponents || adaptInfo.isCoarseningAllowed(i))
        coarsenFlag |= coarseningManager->coarsenMesh(meshes[i]);

    return coarsenFlag;
  }


  void ProblemStatSeq::buildAfterCoarsen(AdaptInfo& adaptInfo, Flag flag,
                                         bool asmMatrix, bool asmVector)
  {
    FUNCNAME("ProblemStat::buildAfterCoarsen()");

    Timer t;

    for (Mesh* mesh : meshes)
      mesh->dofCompress();


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    MPI::COMM_WORLD.Barrier();
    INFO(info, 8)("dof compression needed %.5f seconds\n", t.elapsed());
#endif
    t.reset();

    Flag assembleFlag =
      flag |
      (*systemMatrix)[0][0]->getAssembleFlag() |
      rhs->getDOFVector(0)->getAssembleFlag()  |
      Mesh::CALL_LEAF_EL                        |
      Mesh::FILL_COORDS                         |
      Mesh::FILL_DET                            |
      Mesh::FILL_GRD_LAMBDA |
      Mesh::FILL_NEIGH;

    assembleFlag |= Mesh::FILL_BOUND;


    traverseInfo.updateStatus();
    // Used to calculate the overall number of non zero entries.
    int nnz = 0;

    for (int rowComponent = 0; rowComponent < nComponents; rowComponent++)
    {

      MSG("%d DOFs for %s\n",
          componentSpaces[rowComponent]->getAdmin()->getUsedDofs(),
          componentSpaces[rowComponent]->getName().c_str());

      if (asmVector)
        rhs->getDOFVector(rowComponent)->set(0.0);

      for (int colComponent = 0; colComponent < nComponents; colComponent++)
      {

        if (writeAsmInfo)
          MSG("--------- %d %d -------------\n", rowComponent, colComponent);

        // Only if this variable is true, the current matrix will be assembled.
        bool assembleMatrix = true;
        // The DOFMatrix which should be assembled (or not, if assembleMatrix
        // will be set to false).
        DOFMatrix* matrix =
          (asmMatrix ? (*systemMatrix)[rowComponent][colComponent] : NULL);

        if (writeAsmInfo && matrix)
        {
          MSG("  -> matrix has %d operators\n", matrix->getOperators().size());
          for (vector<Operator*>::iterator it = matrix->getOperatorsBegin();
               it != matrix->getOperatorsEnd(); ++it)
          {
            Assembler* assembler = (*it)->getAssembler();
            if (assembler->getZeroOrderAssembler())
              cout << "ZOA: " << assembler->getZeroOrderAssembler()->getName() << endl;
            if (assembler->getFirstOrderAssembler(GRD_PSI))
              cout << "FOA GRD_PSI: " << assembler->getFirstOrderAssembler(GRD_PSI)->getName() << endl;
            if (assembler->getFirstOrderAssembler(GRD_PHI))
              cout << "FOA GRD_PHI: " << assembler->getFirstOrderAssembler(GRD_PHI)->getName() << endl;
            if (assembler->getSecondOrderAssembler())
              cout << "SOA: " << assembler->getSecondOrderAssembler()->getName() << endl;
          }
        }

        if (matrix)
          matrix->calculateNnz();

        // If the matrix was assembled before and it is marked to be assembled
        // only once, it will not be assembled.
        if (assembleMatrixOnlyOnce[rowComponent][colComponent] &&
            assembledMatrix[rowComponent][colComponent])
        {
          assembleMatrix = false;
        }
        else if (matrix)
        {
          matrix->getBaseMatrix().
          change_dim(componentSpaces[rowComponent]->getAdmin()->getUsedSize(),
                     componentSpaces[colComponent]->getAdmin()->getUsedSize());

          set_to_zero(matrix->getBaseMatrix());
        }

        // If there is no DOFMatrix, e.g., if it is completly 0, do not assemble.
        if (!matrix || !assembleMatrix)
          assembleMatrix = false;

        // If the matrix should not be assembled, the rhs vector has to be considered.
        // This will be only done, if i == j. So, if both is not true, we can jump
        // to the next matrix.
        if (!assembleMatrix && rowComponent != colComponent)
        {
          if (matrix)
            nnz += matrix->getBaseMatrix().nnz();

          continue;
        }

        if (assembleMatrix && matrix->getBoundaryManager())
          matrix->getBoundaryManager()->initMatrix(matrix);

        // This OpenMP barrier is required for the case of periodc boundary
        // conditions. In this case, we may have a data race between
        // fillBoundaryCondition and exitMatrix, where both make use of the
        // vertex vector associated.
        // #pragma omp barrier

        // The simplest case: either the right hand side has no operaters, no aux
        // fe spaces, or all aux fe spaces are equal to the row and col fe space.
        assembleOnOneMesh(componentSpaces[rowComponent],
                          assembleFlag,
                          assembleMatrix ? matrix : NULL,
                          ((rowComponent == colComponent) && asmVector) ?
                          rhs->getDOFVector(rowComponent) :
                          NULL);

        // #pragma omp barrier

        assembledMatrix[rowComponent][colComponent] = true;

        if (assembleMatrix)
          matrix->finishInsertion();

        if (assembleMatrix && matrix->getBoundaryManager())
          matrix->getBoundaryManager()->exitMatrix(matrix);

        if (matrix)
          nnz += matrix->getBaseMatrix().nnz();
      }


      // And now assemble boundary conditions on the vectors
      assembleBoundaryConditions(rhs->getDOFVector(rowComponent),
                                 solution->getDOFVector(rowComponent),
                                 componentMeshes[rowComponent],
                                 assembleFlag);
    }

    if (asmMatrix)
    {
      solverMatrix.setMatrix(*systemMatrix);

      INFO(info, 8)("fillin of assembled matrix: %d\n", nnz);
    }


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    MPI::COMM_WORLD.Barrier();
#endif
    INFO(info, 8)("buildAfterCoarsen needed %.5f seconds\n", t.elapsed());
    buildTime = t.elapsed();
  }

  void ProblemStatSeq::writeFiles(AdaptInfo& adaptInfo, bool force)
  {
    FUNCNAME("ProblemStat::writeFiles()");

    Timer t;

    for (auto fileWriter : fileWriters)
      fileWriter->writeFiles(adaptInfo, force);

    INFO(info, 8)("writeFiles needed %.5f seconds\n", t.elapsed());
  }


  void ProblemStatSeq::interpolInitialSolution(vector<std::function<double(WorldVector<double>)>>& fct)
  {
    solution->interpol(fct);
  }


  void ProblemStatSeq::addMatrixOperator(Operator& op, int i, int j,
                                         double* factor, double* estFactor)
  {
    FUNCNAME("ProblemStat::addMatrixOperator()");

    TEST_EXIT(i < nComponents && j < nComponents)
    ("Cannot add matrix operator at position %d/%d. The stationary problem has only %d components!\n",
     i, j, nComponents);

    TEST_EXIT(!boundaryConditionSet)
    ("Do not add operators after boundary conditions were set!\n");

    if (!(*systemMatrix)[i][j])
    {
      TEST_EXIT(i != j)("should have been created already\n");
      (*systemMatrix)[i][j] = new DOFMatrix(componentSpaces[i], componentSpaces[j], "");
      (*systemMatrix)[i][j]->setCoupleMatrix(true);
      (*systemMatrix)[i][j]->getBoundaryManager()->
      setBoundaryConditionMap((*systemMatrix)[i][i]->getBoundaryManager()->
                              getBoundaryConditionMap());

      if (estimator[i])
        estimator[i]->setNewMatrix(j, (*systemMatrix)[i][j]);
    }

    (*systemMatrix)[i][j]->addOperator(&op, factor, estFactor);

    traverseInfo.getMatrix(i, j).setAuxFeSpaces(op.getAuxFeSpaces());

    OperatorPos opPos = {i, j, factor, estFactor};
    operators[&op].push_back(opPos);
  }


  void ProblemStatSeq::addVectorOperator(Operator& op, int i,
                                         double* factor, double* estFactor)
  {
    FUNCNAME("ProblemStat::addVectorOperator()");

    TEST_EXIT(i < nComponents)
    ("Cannot add vector operator at position %d. The stationary problem has only %d components!\n",
     i, nComponents);

    TEST_EXIT(!boundaryConditionSet)
    ("Do not add operators after boundary conditions were set!\n");

    rhs->getDOFVector(i)->addOperator(&op, factor, estFactor);

    traverseInfo.getVector(i).setAuxFeSpaces(op.getAuxFeSpaces());

    OperatorPos opPos = {i, -1, factor, estFactor};
    operators[&op].push_back(opPos);
  }


  void ProblemStatSeq::addBoundaryMatrixOperator(BoundaryType type,
      Operator& op, int row, int col)
  {
    addRobinBC(type, row, col, (Operator*)NULL, &op);
  }


  void ProblemStatSeq::addBoundaryVectorOperator(BoundaryType type,
      Operator& op, int row)
  {
    addRobinBC(type, row, row, &op, (Operator*)NULL);
  }


  void ProblemStatSeq::addRobinBC(BoundaryType type, int row, int col,
                                  Operator* n, Operator* r)
  {
    // TODO: not yet implemented

    //     boundaryConditionSet = true;
    //
    //     RobinBC *robin =
    //       new RobinBC(type, n, r, componentSpaces[row], componentSpaces[col]);
    //
    //     if (systemMatrix && (*systemMatrix)[row][col])
    //       (*systemMatrix)[row][col]->getBoundaryManager()->addBoundaryCondition(robin);
    //     if (rhs)
    //       rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(robin);
  }


  void ProblemStatSeq::addPeriodicBC(BoundaryType type, int row, int col)
  {
    boundaryConditionSet = true;
    PeriodicBC* periodic = new PeriodicBC(type, componentSpaces[row]);

    if (systemMatrix && (*systemMatrix)[row][col])
      (*systemMatrix)[row][col]->getBoundaryManager()->addBoundaryCondition(periodic);

    if (rhs && row == col)
      rhs->getDOFVector(row)->getBoundaryManager()->addBoundaryCondition(periodic);
  }


  void ProblemStatSeq::assembleOnOneMesh(const FiniteElemSpace* feSpace,
                                         Flag assembleFlag,
                                         DOFMatrix* matrix,
                                         DOFVector<double>* vector)
  {
    Mesh* mesh = feSpace->getMesh();
    const BasisFunction* basisFcts = feSpace->getBasisFcts();

    TraverseStack stack;

    BoundaryType* bound = new BoundaryType[basisFcts->getNumber()];

    if (matrix)
      matrix->startInsertion(matrix->getNnz());

    if (vector)
      vector->set(0.0);

    // == Traverse mesh and assemble. ==
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, assembleFlag);
    while (elInfo)
    {
      basisFcts->getBound(elInfo, bound);

      if (matrix)
      {
        matrix->assemble(1.0, elInfo, bound);

        // assemble the boundary conditions on the matrix.
        if (matrix->getBoundaryManager())
          matrix->getBoundaryManager()->fillBoundaryConditions(elInfo, matrix);
      }

      if (vector)
        vector->assemble(1.0, elInfo, bound, NULL);

      elInfo = stack.traverseNext(elInfo);
    }

    if (matrix)
    {
      matrix->clearDirichletRows();
      matrix->finishAssembling();
    }

    if (vector)
      vector->finishAssembling();

    delete [] bound;
  }


  void ProblemStatSeq::assembleBoundaryConditions(DOFVector<double>* rhs,
						  DOFVector<double>* solution,
						  Mesh* mesh,
						  Flag assembleFlag)
  {
    // === Initialization of vectors ===

    if (rhs->getBoundaryManager())
      rhs->getBoundaryManager()->initVector(rhs);
    if (solution->getBoundaryManager())
      solution->getBoundaryManager()->initVector(solution);

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, assembleFlag);
    while (elInfo)
    {
      if (rhs->getBoundaryManager())
        rhs->getBoundaryManager()->fillBoundaryConditions(elInfo, rhs);

      if (solution->getBoundaryManager())
        solution->getBoundaryManager()->fillBoundaryConditions(elInfo, solution);

      elInfo = stack.traverseNext(elInfo);
    }


    // === Finalize vectors ===

    if (rhs->getBoundaryManager())
      rhs->getBoundaryManager()->exitVector(rhs);
    if (solution->getBoundaryManager())
      solution->getBoundaryManager()->exitVector(solution);
  }


  void ProblemStatSeq::writeResidualMesh(int comp,
                                         AdaptInfo& adaptInfo,
                                         string name)
  {
    map<int, double> vec;
    TraverseStack stack;
    ElInfo* elInfo =
      stack.traverseFirst(this->getMesh(comp),  -1, Mesh::CALL_LEAF_EL);

    while (elInfo)
    {
      vec[elInfo->getElement()->getIndex()] =
        elInfo->getElement()->getEstimation(comp);
      elInfo = stack.traverseNext(elInfo);
    }

    io::ElementFileWriter::writeFile(vec, this->getMesh(comp), name);
  }

} // end namespace AMDiS
