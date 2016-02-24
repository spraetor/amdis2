#pragma once

#include <functional>
#include <list>
#include <vector>

#include "AMDiS_fwd.hpp"
#include "Boundary.hpp"
#include "ComponentTraverseInfo.hpp"
#include "Initfile.hpp"
#include "MatrixVector.hpp"
#include "ProblemStatBase.hpp"
#include "StandardProblemIteration.hpp"
#include "SystemVector.hpp"
#include "io/ElementFileWriter.hpp"
#include "solver/SolverMatrix.hpp"

namespace AMDiS
{
  struct OperatorPos
  {
    int row, col;
    double* factor, *estFactor;
  };


  /**
   * \ingroup Problem
   *
   * \brief
   * This class defines the stationary problem definition in sequential
   * computations. For parallel computations, see \ref ParallelProblemStatBase.
   **/
  class ProblemStatSeq : public ProblemStatBase
  {
  public:
    /// Constructor
    /**
     * Parameters read by ProblemStatSeq, with name 'PROB'
     *   PROB->components:             \ref nComponents
     *   PROB->additional components:  \ref nAddComponents
     *   PROB->names:                  \ref componentNames
     *   PROB->name[i]:                \ref componentNames[i] i=0...nComponents
     *   debug->write asm info:        \ref writeAsmInfo
     **/
    explicit ProblemStatSeq(std::string name);

    /// Destructor
    virtual ~ProblemStatSeq();

    /// Initialisation of the problem.
    /**
     * Parameters read in initialize() for problem with name 'PROB'
     *   MESH[0]->global refinements:  nr of initial global refinements
     *   MESH[0]->value file name:     read value file and use it for the mesh values
     **/
    virtual void initialize(Flag initFlag,
                            ProblemStatSeq* adoptProblem = NULL,
                            Flag adoptFlag = INIT_NOTHING);

    /// Used in \ref initialize().
    /**
     * Parameters read in createMesh() for problem with name 'PROB'
     *   PROB->mesh:              Name of the mesh of this problem, creates one
     *                            or more meshes (w.r.t. different refinement sets)
     *                            and store it in \ref meshes and \ref componentsMeshes.
     *   PROB->dim:               Problem dimension, i.e. dim of the elements
     *   PROB->refinement set[i]: Index for an individual refinement of mesh
     *                            of component i, for i=0...(nComponents+nAddComponents)
     **/
    virtual void createMesh();

    /// Used in \ref initialize().
    virtual void createRefCoarseManager();

    /// Used in \ref initialize().
    /**
     * Parameters read in createFeSpace() for problem with name 'PROB', 
     * for i=0...(nComponents+nAddComponents).
     *   PROB->feSpace[i]:               Name of the i-th feSpace
     *   PROB->finite element space[i]:  Name of the i-th feSpace
     *   PROB->polynomial degree[i]:     Same as feSpace[i] = "Lagrange" + DEGREE
     * 
     * FeSpaces are stored in \ref feSpaces and \ref componentSpaces.
     **/
    virtual void createFeSpace(DOFAdmin* admin);

    /// Used in \ref initialize().
    virtual void createMatricesAndVectors();

    /// Used in \ref initialize().
    /**
     * Parameters read in createMatricesAndVectors() for problem with name 'PROB'
     *   PROB->solver:              Name of the iterative solver
     *   PROB->solver->backend:     Name of the backend the iterative solver belongs to
     **/
    virtual void createSolver();

    /// Used in \ref initialize().
    /**
     * Parameters read in createEstimator() for problem with name 'PROB'
     *   PROB->estimator[i]->name:  Name of the i-th estimator
     * 
     * Estimators are stored in \ref estimator.
     **/
    virtual void createEstimator();

    /// Used in \ref initialize().
    /**
     * Creates a marker corresponding to parameter
     *   PROB->marker[i]
     * 
     * Markers are stored in \ref marker.
     **/
    virtual void createMarker();

    /// Used in \ref initialize().
    /**
     * Creates filewriters corresponding to parameter
     *   PROB->output:              Filewriter for all components
     *   PROB->output->filename:    Filename for the filewriter
     *   PROB->output[i]:           Filewriter for component i
     *   PROB->output[i]->filename: Filename for the filewriter i
     *   PROB->output->num vectors: Combines some components to a vector filewriter
     *   PROB->output->vectors[j]:  Vector-Filewriter with index j=0...num vectors
     *   PROB->output->vectors[j]->components:   The Components corresponding to Vector-Filewriter j
     *   PROB->output->vectors[j]->filename:     The filename of Vector-Filewriter j
     *   PROB->output->vectors[j]->name:         The corresponding component-name
     * 
     * Filewriters are stored in \ref fileWriters.
     **/
    virtual void createFileWriter();

    /// Implementation of ProblemStatBase::solve(). Deligates the solving
    /// of problems system to \ref solver.
    virtual void solve(AdaptInfo& adaptInfo,
                       bool createMatrixData = true,
                       bool storeMatrixData = false) override;

    /// Implementation of ProblemStatBase::estimate(). Deligates the estimation
    /// to \ref estimator.
    virtual void estimate(AdaptInfo& adaptInfo) override;

    /// Implementation of ProblemStatBase::markElements().
    /// Deligated to \ref adapt.
    virtual Flag markElements(AdaptInfo& adaptInfo) override;

    /// Implementation of ProblemStatBase::refineMesh(). Deligated to the
    /// RefinementManager of \ref adapt.
    virtual Flag refineMesh(AdaptInfo& adaptInfo) override;

    /// Implementation of ProblemStatBase::coarsenMesh(). Deligated to the
    /// CoarseningManager of \ref adapt.
    virtual Flag coarsenMesh(AdaptInfo& adaptInfo) override;

    /// Implementation of ProblemStatBase::buildBeforeRefine().
    /// Does nothing here.
    virtual void buildBeforeRefine(AdaptInfo& adaptInfo, Flag) override {}

    /// Implementation of ProblemStatBase::buildBeforeCoarsen().
    /// Does nothing here.
    virtual void buildBeforeCoarsen(AdaptInfo& adaptInfo, Flag) override {}

    /// Implementation of ProblemStatBase::buildAfterCoarsen().
    /// Assembles \ref A and \ref rhs. With the last two parameters, assembling
    /// can be restricted to matrices or vectors only.
    virtual void buildAfterCoarsen(AdaptInfo& adaptInfo, Flag flag,
                                   bool assembleMatrix = true,
                                   bool assembleVector = true) override;

    void assemble(AdaptInfo& adaptInfo)
    {
      buildAfterCoarsen(adaptInfo, Flag(0));
    }

    /// Returns nr of components \ref nComponents
    virtual int getNumComponents() const
    {
      return nComponents;
    }

    /// Returns nr of additional components \ref nAddComponents
    virtual int getNumAddComponents() const
    {
      return nAddComponents;
    }
    
    /// Writes output files.
    void writeFiles(AdaptInfo& adaptInfo, bool force);

    /// Interpolates fct to \ref solution.
    void interpolInitialSolution(std::vector<std::function<double(WorldVector<double>)>>& fct);

    /// Adds an operator to \ref A.
    void addMatrixOperator(Operator& op, int i, int j,
                           double* factor = NULL, double* estFactor = NULL);

    /// Adds an operator to \ref rhs.
    void addVectorOperator(Operator& op, int i,
                           double* factor = NULL, double* estFactor = NULL);

    /// add boundary operator to matrix side
    virtual void addBoundaryMatrixOperator(BoundaryType type,
                                           Operator& op, int row, int col);

    /// add boundary operator to rhs (vector) side
    virtual void addBoundaryVectorOperator(BoundaryType type,
                                           Operator& op, int row);

    /// Adds a Dirichlet boundary condition, where the rhs is given by an
    /// abstract function.
    template <class Term>
    void addDirichletBC(BoundaryType type, int row, int col, Term&& term);

    /// Adds a Neumann boundary condition, where the rhs is given by an
    /// abstract function.
    template <class Term>
    void addNeumannBC(BoundaryType type, int row, int col, Term&& term);

    /// Adds Robin boundary condition.
    template <class ExprRhs, class ExprLhs>
    void addRobinBC(BoundaryType type, int row, int col,
                    ExprRhs&& exprRhs, ExprLhs&& exprLhs);

    /// Adds Robin boundary condition.
    // TODO: replace by generic expression
    virtual void addRobinBC(BoundaryType type, int row, int col,
                            Operator* n, Operator* r);

    /// Adds a periodic boundary condition.
    virtual void addPeriodicBC(BoundaryType type, int row, int col);

    /// This function assembles a DOFMatrix and a DOFVector for the case,
    /// the meshes from row and col FE-space are equal.
    void assembleOnOneMesh(FiniteElemSpace const* feSpace,
                           Flag assembleFlag,
                           DOFMatrix* matrix, DOFVector<double>* vector);


    ///
    void assembleBoundaryConditions(DOFVector<double>* rhs,
                                    DOFVector<double>* solution,
                                    Mesh* mesh,
                                    Flag assembleFlag);


    /** \name getting methods
     * \{
     */

    /// Returns \ref solution.
    SystemVector* getSolution()
    {
      return solution;
    }

    DOFVector<double>* getSolution(int i)
    {
      return solution->getDOFVector(i);
    }

    /// Returns \ref rhs.
    SystemVector* getRhs()
    {
      return rhs;
    }

    DOFVector<double>* getRhsVector(int i = 0)
    {
      return rhs->getDOFVector(i);
    }

    /// Returns \ref systemMatrix.
    Matrix<DOFMatrix*>* getSystemMatrix()
    {
      return systemMatrix;
    }

    /// Returns a pointer to the corresponding DOFMatrix.
    DOFMatrix* getSystemMatrix(int row, int col)
    {
      return (*systemMatrix)[row][col];
    }

    /// Returns mesh of given component
    Mesh* getMesh(int comp = 0)
    {
      FUNCNAME("ProblemStatSeq::getMesh()");
      TEST_EXIT(comp < static_cast<int>(componentMeshes.size()) && comp >= 0)
      ("invalid component number\n");
      return componentMeshes[comp];
    }

    /// Returns \ref meshes
    std::vector<Mesh*>& getMeshes()
    {
      return meshes;
    }

    /// Returns \ref feSpace_.
    FiniteElemSpace const* getFeSpace(int comp = 0) const
    {
      FUNCNAME("ProblemStatSeq::getFeSpace()");
      TEST_EXIT(comp < static_cast<int>(componentSpaces.size()) && comp >= 0)
      ("invalid component number: %d\n", comp);
      return componentSpaces[comp];
    }

    /// Returns \ref feSpaces.
    std::vector<FiniteElemSpace const*>& getFeSpaces()
    {
      return feSpaces;
    }

    /// Returns \ref componentSpaces;
    std::vector<FiniteElemSpace const*>& getComponentSpaces()
    {
      return componentSpaces;
    }

    /// Returns \ref estimator.
    std::vector<Estimator*>& getEstimators()
    {
      return estimator;
    }

    /// Returns \ref estimator.
    Estimator* getEstimator(int comp = 0)
    {
      return estimator[comp];
    }

    /// Returns \ref refinementManager.
    RefinementManager* getRefinementManager(int comp = 0)
    {
      return refinementManager;
    }

    /// Returns \ref refinementManager.
    CoarseningManager* getCoarseningManager(int comp = 0)
    {
      return coarseningManager;
    }

    /// Returns \ref solver.
    LinearSolverInterface* getSolver()
    {
      return solver;
    }

    /// Returns \ref marker.
    Marker* getMarker(int comp = 0)
    {
      return marker[comp];
    }

    /// Returns \ref marker.
    std::vector<Marker*>& getMarkers()
    {
      return marker;
    }

    /// Returns the name of the problem
    virtual std::string getName() const override
    {
      return name;
    }

    /// Returns the name of the problem
    std::string getComponentName(int comp = 0)
    {
      TEST_EXIT(comp < static_cast<int>(componentNames.size()) && comp >= 0)
      ("invalid component number\n");
      return componentNames[comp];
    }

    /// Returns \ref info.
    int getInfo() const
    {
      return info;
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref estimator.
    void setEstimator(std::vector<Estimator*> est)
    {
      estimator = est;
    }

    /// Sets the FE space for the given component.
    void setFeSpace(FiniteElemSpace const* feSpace, int comp = 0)
    {
      feSpaces[comp] = feSpace;
    }

    void setFeSpaces(std::vector<FiniteElemSpace const*> feSpaces_)
    {
      feSpaces = feSpaces_;
    }

    void setComponentSpace(int comp, FiniteElemSpace const* feSpace)
    {
      if (static_cast<int>(componentSpaces.size()) < nComponents)
        componentSpaces.resize(nComponents);
      TEST_EXIT(comp >= 0 && comp < nComponents + nAddComponents)
      ("Component number not in feasable range!");

      componentSpaces[comp] = feSpace;
    }

    /// Sets \ref estimator.
    void setEstimator(Estimator* est, int comp = 0)
    {
      estimator[comp] = est;
    }

    /// Sets \ref marker.
    void setMarker(Marker* mark, int comp = 0)
    {
      marker[comp] = mark;
    }

    /// Sets \ref solver.
    void setSolver(LinearSolverInterface* sol)
    {
      solver = sol;
    }

    void setSolverMatrix(SolverMatrix<Matrix<DOFMatrix*>>& solverMatrix_)
    {
      solverMatrix.setMatrix(*solverMatrix_.getOriginalMat());
    }

    ///
    void setAssembleMatrixOnlyOnce(int i = 0, int j = 0, bool value = true)
    {
      assembleMatrixOnlyOnce[i][j] = value;
    }

    /// Sets \ref writeAsmInfo;
    void setWriteAsmInfo(bool b)
    {
      writeAsmInfo = b;
    }

    void setMeshes(std::vector<Mesh*> meshes_)
    {
      meshes = meshes_;
      nMeshes = static_cast<int>(meshes.size());
    }

    void setComponentMesh(int comp, Mesh* mesh)
    {
      if (static_cast<int>(componentMeshes.size()) < nComponents)
        componentMeshes.resize(nComponents);
      TEST_EXIT(comp >= 0 && comp < nComponents + nAddComponents)
      ("Component number not in feasable range!");

      componentMeshes[comp] = mesh;
    }

    void setRefinementManager(RefinementManager* ref)
    {
      refinementManager = ref;
    }

    void setCoarseningManager(CoarseningManager* coarse)
    {
      coarseningManager = coarse;
    }
    /** \} */

    /// Outputs the mesh of the given component, but the values are taken from
    /// the residual error estimator.
    void writeResidualMesh(int comp, AdaptInfo& adaptInfo, std::string name);

    /// Returns \ref fileWriters.
    std::vector<FileWriterInterface*>& getFileWriterList()
    {
      return fileWriters;
    }

    /// Returns \ref solutionTime.
    double getSolutionTime() const
    {
      return solutionTime;
    }

    /// Returns \ref buildTime.
    double getBuildTime() const
    {
      return buildTime;
    }

  protected:

    /// Name of this problem.
    std::string name;

    /// Number of problem components
    int nComponents = -1;

    /// Number of additional components
    int nAddComponents = 0;

    /// Stores the names for all components. Is used for naming the solution
    /// vectors, \ref solution.
    std::vector<std::string> componentNames;

    /// Number of problem meshes. If all components are defined on the same mesh,
    /// this number is 1. Otherwise, this variable is the number of different meshes
    /// within the problem.
    int nMeshes = 0;

    /// FE spaces of this problem.
    std::vector<FiniteElemSpace const*> feSpaces;

    /// Meshes of this problem.
    std::vector<Mesh*> meshes;

    /// Pointer to the fe spaces for the different problem components
    std::vector<FiniteElemSpace const*> componentSpaces;

    /// Pointer to the meshes for the different problem components
    std::vector<Mesh*> componentMeshes;

    /// Stores information about which meshes must be traversed to assemble the
    /// specific components. I.e., it was implemented to make use of different
    /// meshes for different components.
    ComponentTraverseInfo traverseInfo = 0;

    /// Responsible for element marking.
    std::vector<Marker*> marker;

    /// Estimator of this problem. Used in \ref estimate().
    std::vector<Estimator*> estimator;

    /// Linear solver of this problem. Used in \ref solve().
    LinearSolverInterface* solver = NULL;

    /// System vector  storing the calculated solution of the problem.
    SystemVector* solution = NULL;

    /// System vector for the right hand side
    SystemVector* rhs = NULL;

    /// System matrix
    Matrix<DOFMatrix*>* systemMatrix = NULL;

    /// Composed system matrix
    SolverMatrix<Matrix<DOFMatrix*>> solverMatrix;

    /// Some DOFMatrices of the systemMatrix may be assembled only once (for
    /// example if they are independent of the time or older solutions). If
    /// [i][j] of this field is set to true, the corresponding DOFMatrix will
    /// be assembled only once. All other matrices will be assembled at every
    /// time step.
    std::vector<std::vector<bool>> assembleMatrixOnlyOnce;

    /// If [i][j] of this field is set to true, the corresponding DOFMatrix of
    /// the systemMatrix has been assembled at least once. This field is used
    /// to determine, if assembling of a matrix can be ommitted, if it is set
    /// to be assembled only once.
    std::vector<std::vector<bool>> assembledMatrix;

    /// Writes the meshes and solution after the adaption loop.
    std::vector<FileWriterInterface*> fileWriters;

    /// All actions of mesh refinement are performed by refinementManager.
    /// If new refinement algorithms should be realized, one has to override
    /// RefinementManager and give one instance of it to AdaptStationary.
    RefinementManager* refinementManager = NULL;

    /// All actions of mesh coarsening are performed by coarseningManager.
    /// If new coarsening algorithms should be realized, one has to override
    /// CoarseningManager and give one instance of it to AdaptStationary.
    CoarseningManager* coarseningManager = NULL;

    /// Info level.
    int info = 10;

    /// If at least on boundary condition is set, this variable is true. It is
    /// used to ensure that no operators are added after boundary condition were
    /// set. If this would happen, boundary conditions could set wrong on off
    /// diagonal matrices.
    bool boundaryConditionSet = false;

    /// If true, AMDiS prints information about the assembling procedure to
    /// the screen.
    bool writeAsmInfo = false;

    std::map<Operator*, std::vector<OperatorPos>> operators;

    /// time needed to solve the linear system
    double solutionTime = 0.0;

    /// time needed to assemble the linear system
    double buildTime = 0.0;
  };

  namespace detail
  {
    template <class ProblemStatType>
    struct ProblemStat : public ProblemStatType,
                         public StandardProblemIteration
    {
      using ProblemStatType::getName;

      /// Constructor
      explicit ProblemStat(std::string nameStr)
        : ProblemStatType(nameStr),
          StandardProblemIteration(dynamic_cast<ProblemStatBase&>(*this))
      {}

      /// Determines the execution order of the single adaption steps. If adapt is
      /// true, mesh adaption will be performed. This allows to avoid mesh adaption,
      /// e.g. in timestep adaption loops of timestep adaptive strategies.
      // implements StandardProblemIteration::oneIteration(AdaptInfo&, Flag)
      virtual Flag oneIteration(AdaptInfo& adaptInfo, Flag toDo = FULL_ITERATION) override
      {
        for (int i = 0; i < ProblemStatType::getNumComponents(); i++)
          if (adaptInfo.spaceToleranceReached(i))
            adaptInfo.allowRefinement(false, i);
          else
            adaptInfo.allowRefinement(true, i);

        return StandardProblemIteration::oneIteration(adaptInfo, toDo);
      }
    };
  }

#ifndef HAVE_PARALLEL_DOMAIN_AMDIS
  using ProblemStat = detail::ProblemStat<ProblemStatSeq>;
#endif

} // end namespace AMDiS

#include "ProblemStat.hh"
