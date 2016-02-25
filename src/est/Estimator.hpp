#pragma once

// AMDiS includes
#include "AMDiS_fwd.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"
#include "Global.hpp"
#include "SystemVector.hpp"
#include "CreatorInterface.hpp"
#include "ComponentTraverseInfo.hpp"
#include "DualTraverse.hpp"

namespace AMDiS
{

  /**
   * \ingroup Estimator
   *
   * \brief
   * Estimator for scalar problems.
   */
  class Estimator
  {
  public:
    Estimator() {}

    /// Constructor.
    Estimator(std::string name_, int r);

    /// destructor
    virtual ~Estimator() {}

    /// Returns \ref name of the Estimator
    inline std::string getName() const
    {
      return name;
    }

    /// Performs the estimation and returns the final \ref est_sum
    virtual double estimate(double timestep = 0.0);

    ///
    virtual void init(double timestep) =0;

    /** \brief
     * Estimates the error on an element. If there is more than one mesh used in the
     * problem description, it may be necessary to used the dual mesh traverse. In this
     * case elInfo is the element of the main mesh, i.e., the mesh of the row FE space,
     * and dualElInfo contains all elInfo informations about the main mesh element and
     * the col (or aux) mesh element.
     */
    virtual void estimateElement(ElInfo* elInfo, DualElInfo* dualElInfo = NULL) =0;

    ///
    virtual void exit(bool output = true) =0;

    /// Returns \ref est_sum of the Estimator
    inline double getErrorSum() const
    {
      return est_sum;
    }

    /// Sets \ref est_sum of the Estimator
    inline void setErrorSum(double sum)
    {
      est_sum = sum;
    }

    /// Returns \ref est_max of the Estimator
    inline double getErrorMax() const
    {
      return est_max;
    }

    /// Returns the estimated time error.
    virtual double getTimeEst() const
    {
      return est_t_sum;
    }

    /// Returns the maximal time estimation.
    virtual double getTimeEstMax() const
    {
      return est_t_max;
    }

    /// Sets \ref est_max of the Estimator
    inline void setErrorMax(double m)
    {
      est_max = m;
    }

    /// Returns \ref norm.
    inline Norm getErrorNorm()
    {
      return norm;
    }

    /// Adds one system to the estimator.
    virtual void addSystem(DOFMatrix* matrix_,
                           DOFVector<double>* uh_,
                           DOFVector<double>* fh_,
                           DOFVector<double>* uhOld_ = NULL)
    {
      matrix.push_back(matrix_);
      uh.push_back(uh_);
      fh.push_back(fh_);
      uhOld.push_back(uhOld_);
    }

    void setNewMatrix(int i, DOFMatrix* m)
    {
      matrix[i] = m;
    }

    /// Adds pointer to old solution to the given system.
    virtual void addUhOldToSystem(int system, DOFVector<double>* uhOld_)
    {
      FUNCNAME("Estimator::addUhOldToSystem()");

      TEST_EXIT(static_cast<int>(uhOld.size()) > system)("Invalid system!\n");
      TEST_EXIT(uhOld[system] == NULL)("There is already an uhOld!\n");

      uhOld[system] = uhOld_;
    }

    /// Returns number of systems.
    inline int getNumSystems()
    {
      return static_cast<int>(matrix.size());
    }

    inline Flag getTraverseFlag()
    {
      return traverseFlag;
    }

    inline Mesh* getMesh()
    {
      return mesh;
    }

    inline int getRow()
    {
      return row;
    }

    /// Sets \ref traverseInfo.
    void setTraverseInfo(const ComponentTraverseInfo& ti)
    {
      traverseInfo = ti;
    }

  protected:
    /// Traverse one mesh to estimate the error.
    void singleMeshTraverse();

    /// Traverses two meshes to estimate the error.
    void dualMeshTraverse();

  protected:
    /// Name of the Estimator
    std::string name;

    /// Used norm
    Norm norm = NO_NORM;

    /// Sum of all error estimates
    double est_sum = 0.0;

    /// Maximal error estimate
    double est_max = 0.0;

    /// Sum of all time error estimates
    double est_t_sum = 0.0;

    /// Max of all time error estimates
    double est_t_max = 0.0;

    /** \brief
     * Vector of DOFMatrix pointers. There can be more than one
     * DOFMatrix, if the Estimator is part of a vector valued
     * estimator. Then it contains also coupling matrices
     * of the different vector components.
     */
    std::vector<DOFMatrix*> matrix;

    /// Vector of solution vectors for the different systems.
    std::vector<DOFVector<double>*> uh;

    /** \brief
     * Vector of old solutions vectors for the different systems.
     * Used for instationary problems.
     */
    std::vector<DOFVector<double>*> uhOld;

    /// Vector of RHS vectors for the different systems.
    std::vector<DOFVector<double>*> fh;

    /** \brief
     * Used, if the scalar estimator builds one row of a vector valued estimator.
     * Then row gives the position in the vector valued estimator, used
     * for calculation of the time derivative.
     */
    int row;

    Flag traverseFlag;

    /** \brief
     * The mesh on which the error must be estimated. If there is more than one mesh
     * used, here the main, i.e., the row mesh, is stored.
     */
    Mesh* mesh = NULL;

    /** \brief
     * If there is only one mesh used at all, this variable is not used. In the case
     * that the error must be estimated on a system row with more than one mesh, here
     * either the column mesh or the auxiliary mesh is stored.
     */
    Mesh* auxMesh = NULL;

    double timestep = 0.0;

    /** \brief
     * Stores information about which mesh(es) must be traversed to estimate
     * the error on the component matrices.
     */
    ComponentTraverseInfo traverseInfo = 0;
  };


  /**
   * \ingroup Estimator
   *
   * \brief
   * Interface for creators of concrete estimators.
   */
  class EstimatorCreator : public CreatorInterface<Estimator>
  {
  public:
    /// constructor
    EstimatorCreator() : row(-1), uh(NULL) {}

    /// destructor
    virtual ~EstimatorCreator() {}

    /// Sets \ref name
    inline void setName(std::string name_)
    {
      name = name_;
    }

    /// Sets \ref row
    inline void setRow(int r)
    {
      row = r;
    }

    inline void setSolution(DOFVector<double>* uh_)
    {
      uh = uh_;
    }

  protected:
    /// Name of the estimator to be created.
    std::string name;

    /// Row of the estimator.
    int row;

    /// Pointer to solution vector
    DOFVector<double>* uh;
  };

}
