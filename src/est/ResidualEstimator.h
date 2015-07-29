/** \file ResidualEstimator.h */

/** \defgroup Estimator Estimator module
 * @{ <img src="estimator.png"> @}
 */

#pragma once

#include "Estimator.h"
#include "FixVec.h"

namespace AMDiS 
{

  /** \brief
   * Returns residual square at quadrature point. Not Member of
   * Estimator to avoid multiple instantiation.
   */
  void r(const ElInfo *elInfo,
	 int nPoints,
	 const DenseVector<double>& uhIq,
	 const DenseVector<WorldVector<double> > &grdUhIq,
	 const DenseVector<WorldMatrix<double> > &D2UhIq,
	 const DenseVector<double>& uhOldIq,
	 const DenseVector<WorldVector<double> > &grdUhOldIq,
	 const DenseVector<WorldMatrix<double> > &D2UhOldIq,
	 DOFMatrix *A, 
	 DOFVector<double> *fh,
	 Quadrature *quad,
	 DenseVector<double>& result);
 
  /// Returns pow(det,2.0/dim). Not Member of Estimator to avoid multiple instantiation.
  inline double h2_from_det(double det, int dim) 
  {
    return std::pow(det, 2.0 / dim);
  }

  /**
   * \ingroup Estimator
   * 
   * \brief
   * Residual estimator.
   */
  class ResidualEstimator : public Estimator
  {
  public:
    class Creator : public EstimatorCreator
    {
    public:
      Creator() 
	: EstimatorCreator() 
      {}

      virtual ~Creator() 
      {}

      /// Returns a new ODirSolver object.
      Estimator* create() 
      { 
	return new ResidualEstimator(name, row);
      }
    };
  
    /// Constructor.
    ResidualEstimator(std::string name, int r);

    void init(double timestep);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    void initParallel();
#endif

    /// Estimates the error on an element. For more information about the
    /// parameter, see the description \ref Estimator::estimateElement.
    void estimateElement(ElInfo *elInfo, DualElInfo *dualElInfo = NULL);

    void exit(bool output = true);

  protected:
    /// Computes the element residual for a given element.
    double computeElementResidual(ElInfo *elInfo, DualElInfo *dualElInfo);

    /// Computes the jump residual for a given element.
    double computeJumpResidual(ElInfo *elInfo, DualElInfo *dualElInfo);

  protected:
    /// Constant in front of element residual
    double C0;
 
    /// Constant in front of edge/face residual
    double C1;

    /// Not used! Was thought to be the constant in front of coarsening term.
    double C2;

    /// Constant in front of the time
    double C3;
    
    /// Is true, if C1 != 0, and C0 and C3 = 0, hence only the jump residual must be
    /// calculated. In this case, some optimizations can be done, because the jump
    /// residual is calculated only on second order terms.
    bool jumpResidualOnly;

    /// Number of systems, e.g., number of variables in the equation.
    int nSystems;

    /// Number of quadrature points.
    int nPoints;

    int dim;

    int degree;

    Quadrature *quad;

    FastQuadrature **quadFast;

    const BasisFunction **basFcts;

    std::vector<ElementVector> uhEl;

    std::vector<ElementVector> uhOldEl;

    std::vector<ElementVector> uhNeigh;

    DenseVector<double> uhQP;

    DenseVector<double> uhOldQP;

    /// Stores the element residual computed at the quadrature points of the element.
    ElementVector riq;

    DenseVector<WorldVector<double> > grdUhQp;

    DenseVector<WorldMatrix<double> > D2UhQp;

    /// not used:
    DenseVector<WorldVector<double> > grdUhOldQp;
    DenseVector<WorldMatrix<double> > D2UhOldQp;
    
    ElInfo *neighInfo;

    Quadrature *surfaceQuad;

    int nPointsSurface;

    std::vector<WorldVector<double> > grdUhEl;

    std::vector<WorldVector<double> > grdUhNeigh;

    std::vector<WorldVector<double> > jump;

    std::vector<WorldVector<double> > localJump;

    WorldVector<int> faceIndEl;
    
    WorldVector<int> faceIndNeigh;

    DimVec<WorldVector<double> > *lambdaNeigh;

    DimVec<double> *lambda;

    /// Maximal number of neighbours an element may have in the used dimension.
    int nNeighbours;

    /// Defines for every system if there are second order terms. These values
    /// are used to ommit computations of the jump residual that is defined
    /// only on second order terms.
    std::vector<bool> secondOrderTerms;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    std::map<BoundaryObject, double> elBoundDet;

    std::map<BoundaryObject, std::vector<WorldVector<double> > > elBoundGrdUhNeigh;
#endif
  };
  
} // end namespace AMDiS
