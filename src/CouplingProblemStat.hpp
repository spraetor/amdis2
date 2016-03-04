#pragma once

// std c++ headers
#include <list>
#include <set>
#include <vector>

// AMDiS includes
#include "AMDiS_fwd.hpp"
#include "Initfile.hpp"
#include "ProblemStat.hpp"

namespace AMDiS
{
  namespace detail
  {
    /** \brief
    * This class defines a coupled stationary problem definition in sequential
    * computations.
    */
    template <class ProblemStatType>
    class CouplingProblemStat : public ProblemStatSeq
    {
    protected:
      typedef ProblemStatSeq  super;

      using super::nComponents;
      using super::meshes;
      using super::nMeshes;
      using super::feSpaces;
      using super::name;
      using super::refinementManager;
      using super::coarseningManager;

    public:
      /// Constructor
      CouplingProblemStat(std::string name_)
        : super(name_),
          dim(-1)
      { }

      /// add problem by number
      virtual void addProblem(ProblemStatType* prob)
      {
        problems.push_back(prob);
        nComponents += prob->getNumComponents();
      };

      /// Initialisation of the problem.
      virtual void initialize(Flag initFlag,
                              ProblemStatSeq* adoptProblem = NULL,
                              Flag adoptFlag = INIT_NOTHING) override
      {
        FUNCNAME("CouplingProblemStat::initialize()");

        super::initialize(initFlag - INIT_MESH);

        const Flag DEFAULT_INIT = (INIT_FE_SPACE | INIT_MESH | CREATE_MESH | INIT_SYSTEM |
                                   INIT_SOLVER | INIT_ESTIMATOR | INIT_MARKER | INIT_FILEWRITER);
        for (size_t p = 0; p < problems.size(); ++p)
        {
          problems[p]->initialize(initFlag - DEFAULT_INIT);
        }

        for (size_t i = 0; i < meshes.size(); i++)
        {
          int globalRefinements = 0;

          // If AMDiS is compiled for parallel computations, the global refinements are
          // ignored here. Later, each rank will add the global refinements to its
          // private mesh.
#ifndef HAVE_PARALLEL_DOMAIN_AMDIS
          Parameters::get(meshes[i]->getName() + "->global refinements",
                          globalRefinements);
#endif

          bool initMesh = initFlag.isSet(INIT_MESH);

          // Initialize the meshes if there is no serialization file.
          if (initMesh && meshes[i] && !(meshes[i]->isInitialized()))
          {
            meshes[i]->initialize();
            refinementManager->globalRefine(meshes[i], globalRefinements);
          }
        }
      }


      /// Used in \ref initialize().
      virtual void createMesh() override
      {
        // all problems must have the same dimension (?)
        dim = 0;
        Parameters::get(name + "->dim", dim);
        TEST_EXIT(dim)("No problem dimension specified for \"%s->dim\"!\n",
                       name.c_str());

        std::map<std::pair<std::string, int>, Mesh*> meshByName;  // (name, refSet) --> Mesh*
        typedef std::map<std::pair<std::string, int>, Mesh*>::iterator MeshIterator;

        std::vector<std::set<Mesh*>> meshesForProblems(problems.size());

        for (size_t i = 0; i < problems.size(); ++i)
        {
          TEST_EXIT(problems[i])("problem[%d] does not exist!\n",i);

          int nComponents = problems[i]->getNumComponents();

          int nAddComponents = 0;
          Parameters::get(problems[i]->getName() + "->additional components", nAddComponents);

          problems[i]->componentMeshes.resize(nComponents + nAddComponents);

          for (size_t j = 0; j < nComponents + nAddComponents; j++)
          {
            // name of the mesh
            std::string meshName("");
            Parameters::get(problems[i]->getName() + "->mesh", meshName);
            TEST_EXIT(meshName != "")("No mesh name specified for \"%s->mesh\"!\n",
                                      problems[i]->getName().c_str());

            // dimension of the mesh
            int mesh_dim = 0;
            Parameters::get(problems[i]->getName() + "->dim", mesh_dim);
            TEST_EXIT(dim == mesh_dim)("Mesh-dimension must be the same for all problems!\n");

            // refinement set (optional)
            int refSet = 0;
            Parameters::get(problems[i]->getName() + "->refinement set[" + to_string(j) + "]", refSet);

            // create a new Mesh only if not already created for other problem
            Mesh* componentMesh;
            MeshIterator meshIt = meshByName.find(std::make_pair(meshName, refSet));
            if (meshIt == meshByName.end())
            {
              Mesh* newMesh = new Mesh(meshName, dim);
              meshes.push_back(newMesh);
              meshByName[std::make_pair(meshName, refSet)] = newMesh;
              componentMesh = newMesh;
              nMeshes++;
            }
            else
            {
              componentMesh = meshIt->second;
            }
            problems[i]->componentMeshes[j] = componentMesh;
          }

          // copy unqiue set of meshes to problem[i]->meshes
          std::set<Mesh*> uniqueMeshes;
          for (size_t j = 0; j < problems[i]->componentMeshes.size(); ++j)
            uniqueMeshes.insert(problems[i]->componentMeshes[j]);
          problems[i]->meshes.clear();
          problems[i]->meshes.insert(problems[i]->meshes.begin(), uniqueMeshes.begin(), uniqueMeshes.end());
        }
      }

      /// Used in \ref initialize().
      virtual void createFeSpace(DOFAdmin* admin) override
      {
        std::vector<std::set<FiniteElemSpace const*>> feSpacesForProblems(problems.size());
        std::map<std::pair<Mesh*, std::string>, FiniteElemSpace*> feSpaceMap;

        for (size_t p = 0; p < problems.size(); ++p)
        {
          TEST_EXIT(problems[p])("problem[%d] does not exist!\n",p);

          int nComponents = problems[p]->getNumComponents();

          int nAddComponents = 0;
          Parameters::get(problems[p]->getName() + "->additional components", nAddComponents);
          problems[p]->componentSpaces.resize(nComponents + nAddComponents, NULL);
          problems[p]->traverseInfo.resize(nComponents);

          for (size_t i = 0; i < nComponents + nAddComponents; i++)
          {

            std::string componentString = "[" + to_string(i) + "]";

            std::string feSpaceName = "";
            std::string initFileStr = problems[p]->getName() + "->feSpace" + componentString;
            Parameters::get(initFileStr, feSpaceName);

            // synonym for "feSpace"
            if (feSpaceName.size() == 0)
            {
              initFileStr = problems[p]->getName() + "->finite element space" + componentString;
              Parameters::get(initFileStr, feSpaceName);
            }

            // for backward compatibility also accept the old syntax
            if (feSpaceName.size() == 0)
            {
              int degree = 1;
              initFileStr = problems[p]->getName() + "->polynomial degree" + componentString;
              Parameters::get(initFileStr, degree);
              TEST_EXIT(degree > 0)
              ("Poynomial degree in component %d must be larger than zero!\n", i);

              feSpaceName = "Lagrange" + to_string(degree);
            }

            if (feSpaceName.size() == 0)
              feSpaceName = "Lagrange1";

            if (feSpaceMap[std::make_pair(problems[p]->componentMeshes[i], feSpaceName)] == NULL)
            {
              BasisFunctionCreator* basisFctCreator =
                dynamic_cast<BasisFunctionCreator*>(CreatorMap<BasisFunction>::getCreator(feSpaceName, initFileStr));
              TEST_EXIT(basisFctCreator)
              ("No valid basisfunction type found in parameter \"%s\"\n", initFileStr.c_str());
              basisFctCreator->setDim(dim);

              FiniteElemSpace* newFeSpace =
                FiniteElemSpace::provideFeSpace(admin, basisFctCreator->create(),
                                                problems[p]->componentMeshes[i],
                                                "FeSpace" + componentString + " (" + feSpaceName + ")");

              feSpaceMap[std::make_pair(problems[p]->componentMeshes[i], feSpaceName)] = newFeSpace;
              feSpaces.push_back(newFeSpace);
            }

            problems[p]->componentSpaces[i] = feSpaceMap[std::make_pair(problems[p]->componentMeshes[i], feSpaceName)];
          }

          // copy unqiue set of meshes to problem[p]->meshes
          std::set<FiniteElemSpace const*> uniqueFeSpaces;
          for (size_t i = 0; i < problems[p]->componentSpaces.size(); ++i)
            uniqueFeSpaces.insert(problems[p]->componentSpaces[i]);
          problems[p]->feSpaces.clear();
          problems[p]->feSpaces.insert(problems[p]->feSpaces.begin(), uniqueFeSpaces.begin(), uniqueFeSpaces.end());

          // create traverseInfo
          for (int i = 0; i < nComponents; i++)
          {
            for (int j = 0; j < nComponents; j++)
              problems[p]->traverseInfo.getMatrix(i, j).setFeSpace(problems[p]->componentSpaces[i], problems[p]->componentSpaces[j]);

            problems[p]->traverseInfo.getVector(i).setFeSpace(problems[p]->componentSpaces[i]);
          }
        }

        // create dof admin for vertex dofs if neccessary
        for (size_t i = 0; i < meshes.size(); i++)
        {
          if (meshes[i]->getNumberOfDofs(VERTEX) == 0)
          {
            DimVec<int> ln_dof(meshes[i]->getDim(), 0);
            ln_dof[VERTEX] = 1;
            meshes[i]->createDOFAdmin("vertex dofs", ln_dof);
          }
        }
      }

      /// Used in \ref initialize().
      virtual void createRefCoarseManager() override
      {
        FUNCNAME("CouplingProblemStat::createRefCoarseManager()");
        assert( dim > 0);

        switch (dim)
        {
        case 1:
          coarseningManager = new CoarseningManager1d();
          refinementManager = new RefinementManager1d();
          break;
        case 2:
          coarseningManager = new CoarseningManager2d();
          refinementManager = new RefinementManager2d();
          break;
        case 3:
          coarseningManager = new CoarseningManager3d();
          refinementManager = new RefinementManager3d();
          break;
        default:
          ERROR_EXIT("invalid dim!\n");
        }

        for (size_t p = 0; p < problems.size(); p++)
        {
          problems[p]->setRefinementManager(refinementManager);
          problems[p]->setCoarseningManager(coarseningManager);
        }
      }

      /// Used in \ref initialize().
      virtual void createMatricesAndVectors() override
      {
        for (size_t p = 0; p < problems.size(); ++p)
        {
          assert( problems[p] );
          problems[p]->createMatricesAndVectors();
        }
      }

      /// Used in \ref initialize().
      virtual void createSolver() override
      {
        for (size_t p = 0; p < problems.size(); ++p)
        {
          assert( problems[p] );
          problems[p]->createSolver();
        }
      }

      /// Used in \ref initialize().
      virtual void createEstimator() override
      {
        for (size_t p = 0; p < problems.size(); ++p)
        {
          assert( problems[p] );
          problems[p]->createEstimator();
        }
      }

      /// Used in \ref initialize().
      virtual void createMarker() override
      {
        for (size_t p = 0; p < problems.size(); ++p)
        {
          assert( problems[p] );
          problems[p]->createMarker();
        }
      }

      /// Used in \ref initialize().
      virtual void createFileWriter() override
      {
        for (size_t p = 0; p < problems.size(); ++p)
        {
          assert( problems[p] );
          problems[p]->createFileWriter();
        }
      }


      /// Returns number of managed problems
      virtual int getNumProblems() const
      {
        return problems.size();
      }

      /** \brief
      * Returns the problem with the given number. If only one problem
      * is managed by this master problem, the number hasn't to be given.
      */
      virtual ProblemStatType& getProblem(int number = 0)
      {
        return *problems[number];
      }

      /// Returns \ref meshes[i]
      inline Mesh* getMesh(int number = 0)
      {
        return meshes[number];
      }

      using super::getNumComponents;
      using super::getRefinementManager;
      using super::getCoarseningManager;
      using super::getName;

    protected:
      /// unqiue mesh-dimension for all problems
      int dim;

      std::vector<ProblemStatType*> problems;
    };

  } // end namespace detail

  typedef detail::CouplingProblemStat<ProblemStat> CouplingProblemStat;

} // end namespace AMDiS
