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



/** \file CFE_NormAndErrorFcts.h */

#ifndef AMDIS_CFE_NORMANDERRORFCTS_H
#define AMDIS_CFE_NORMANDERRORFCTS_H

#include "AbstractFunction.h"
#include "DOFVector.h"
#include "FixVec.h"
#include "Quadrature.h"

#include "ElementLevelSet.h"
#include "ScalableQuadrature.h"
#include "SubPolytope.h"


namespace compositeFEM {
  
  using namespace AMDiS;
  
  class ElementNorm
  {
  public:
    /// Constructor.
    ElementNorm(Quadrature *q_)
      : q(q_),
	nQPts(0)
    {
      if (q)
	nQPts = q->getNumPoints();
    }

    /// Destructor.
    virtual ~ElementNorm() {}

    /// Calculates element norm on elInfo.
    virtual double calcElNorm(ElInfo *elInfo, 
			      const double &det, 
			      const double &fac = 1.0) = 0;

    /// Sets quadrature to q_.
    inline void setQuadrature(Quadrature *q_) 
    {
      q = q_;
      nQPts = q->getNumPoints();
    }

  protected:
    /// Quadrature formula.
    Quadrature *q;

    /// Number of quadrature points.
    int nQPts;
  };

  class ElementL1Norm_Analyt : public ElementNorm
  {
  public:
    /// Constructor.
    ElementL1Norm_Analyt(Quadrature *q_, 
			 AbstractFunction<double, WorldVector<double> > *f_)
      : ElementNorm(q_),
	f(f_)
    {}

    /// Calculates element norm on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

  protected:
    /// Abstract function for which norm is calculated.
    AbstractFunction<double, WorldVector<double> > *f;
  };

  class ElementL2Norm_Analyt : public ElementNorm
  {
  public:
    /// Constructor.
    ElementL2Norm_Analyt(Quadrature *q_, 
			 AbstractFunction<double, WorldVector<double> > *f_)
      : ElementNorm(q_),
	f(f_)
    {}

    /// Calculates element norm on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

  protected:
    /// Abstract function for which norm is calculated.
    AbstractFunction<double, WorldVector<double> > *f;
  };

  class ElementH1Norm_Analyt : public ElementNorm
  {
  public:
    /// Constructor.
    ElementH1Norm_Analyt(Quadrature *q_, 
			 AbstractFunction<WorldVector<double>, WorldVector<double> > *grd_,
			 int dim_)
      : ElementNorm(q_),
	grd(grd_),
	dim(dim_)
    {}

    /// Calculates element norm on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

  protected:
    /// Abstract function for which norm is calculated.
    AbstractFunction<WorldVector<double>, WorldVector<double> > *grd; 

    /// Mesh dimension.
    int dim;
  };

  class ElementL1Norm_DOF : public ElementNorm
  {
  public:
    /// Constructor.
    ElementL1Norm_DOF(Quadrature *q_, DOFVector<double> *dofVec_)
      : ElementNorm(q_),
	dofVec(dofVec_)
    {}

    /// Calculates element norm on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

  protected:
    /// DOF vector for which norm is calculated.
    DOFVector<double> *dofVec;
  };

  class ElementL2Norm_DOF : public ElementNorm
  {
  public:
    /// Constructor.
    ElementL2Norm_DOF(Quadrature *q_, DOFVector<double> *dofVec_)
      : ElementNorm(q_),
	dofVec(dofVec_)
    {}

    /// Calculates element norm on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

  protected:
    /// DOF vector for which norm is calculated.
    DOFVector<double> *dofVec;
  };

  class ElementH1Norm_DOF : public ElementNorm
  {
  public:
    /// Constructor.
    ElementH1Norm_DOF(Quadrature *q_, DOFVector<double> *dofVec_, int dim_)
      : ElementNorm(q_),
	dofVec(dofVec_),
	dim(dim_)
    {}

    /// Calculates element norm on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

  protected:
    /// DOF vector for which norm is calculated.
    DOFVector<double> *dofVec;

    /// Mesh dimension.
    int dim;
  };

  class ElementL2Err : public ElementNorm
  {
  public:
    /// Constructor.
    ElementL2Err(Quadrature *q_, 
		 AbstractFunction<double, WorldVector<double> > *u_,
		 DOFVector<double> *uh_,
		 int relErr_)
      : ElementNorm(q_),
	u(u_),
	uh(uh_),
	relErr(relErr_),
	nrmU(0.0)
    {}

    /// Calculates element error on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

    /// Get norm of u.
    inline double getNormU() const 
    {
      return nrmU;
    }

  protected:
    /// Abstract function for which error is calculated.
    AbstractFunction<double, WorldVector<double> > *u;

    /// DOF vector for which error is calculated.
    DOFVector<double> *uh;

    /// Indicates whether relative (1) or absolute error (0) is calculated.
    int relErr;

    /// Norm of u in case relative error is calculated.
    double nrmU;
  };

  class ElementH1Err : public ElementNorm
  {
  public:
    /// Constructor.
    ElementH1Err(Quadrature *q_, 
		 AbstractFunction<WorldVector<double>, WorldVector<double> > *grdu_,
		 DOFVector<double> *uh_,
		 int relErr_,
		 int dim_)
      : ElementNorm(q_),
	grdu(grdu_),
	uh(uh_),
	relErr(relErr_),
	nrmGrdU(0.0),
	dim(dim_)
    {}

    /// Calculates element error on elInfo.
    double calcElNorm(ElInfo *elInfo, 
		      const double &det,
		      const double &fac = 1.0);

    /// Get norm of grdu.
    inline double getNormGrdU() const 
    {
      return nrmGrdU;
    }

  protected:
    /// Abstract function for which norm is calculated.
    AbstractFunction<WorldVector<double>, WorldVector<double> > *grdu; 

    /// DOF vector for which error is calculated.
    DOFVector<double> *uh;

    /// Indicates whether relative (1) or absolute error (0) is calculated.
    int relErr;

    /// Norm of grdu in case relative error is calculated.
    double nrmGrdU;

    /// Mesh dimension.
    int dim;
  };

  class CFE_NormAndErrorFcts
  {
  public:
    // ========================================================================
    //  Calculation of L1-Norm, L2-Norm and H1-Norm on domain depending
    //  on the level set function.
    //
    //  Here, the H1-norm is the H1-seminorm, i.e. you get the "common"
    //  H1-norm by
    //     H1_Norm() = sqrt(L2_Norm^2() + H1_Seminorm^2()) .
    //
    //  Parameters for domainFlag:
    //      -3:  all elements in the domain with negative level set function
    //           values which do not intersect the zero level set
    //      -2 : all elements which intersect the domain with negative level set
    //           function values, in particular all (complete) elements which
    //           cut the zero level set
    //      -1 : domain with negative level set function values
    //       0 : all elements cut by the zero level set
    //       1 : complete mesh
    // ========================================================================
    static double L1Norm_Analyt(AbstractFunction<double, WorldVector<double> > *f,
				ElementLevelSet *elLS,
				int domainFlag,
				int deg = 1,
				Quadrature* q = NULL);
    static double L2Norm_Analyt(AbstractFunction<double, WorldVector<double> > *f,
				ElementLevelSet *elLS,
				int domainFlag,
				int deg = 2,
				Quadrature* q = NULL);
    static double L2NormSquare_Analyt(AbstractFunction<double, WorldVector<double> > *f,
				      ElementLevelSet *elLS,
				      int domainFlag,
				      int deg = 2,
				      Quadrature* q = NULL);
    static double H1Norm_Analyt(AbstractFunction<WorldVector<double>, WorldVector<double> > *grd,
				ElementLevelSet *elLS,
				int domainFlag,
				int deg = 0,
				Quadrature* q = NULL);
    static double H1NormSquare_Analyt(AbstractFunction<WorldVector<double>, WorldVector<double> > *grd,
				      ElementLevelSet *elLS,
				      int domainFlag,
				      int deg = 0,
				      Quadrature* q = NULL);
    static double L1Norm_DOF(DOFVector<double> *dof,
			     ElementLevelSet *elLS,
			     int domainFlag,
			     int deg = 1,
			     Quadrature* q = NULL);
    static double L2Norm_DOF(DOFVector<double> *dof,
			     ElementLevelSet *elLS,
			     int domainFlag,
			     int deg = 2,
			     Quadrature* q = NULL);
    static double L2NormSquare_DOF(DOFVector<double> *dof,
				   ElementLevelSet *elLS,
				   int domainFlag,
				   int deg = 2,
				   Quadrature* q = NULL);
    static double H1Norm_DOF(DOFVector<double> *dof,
			     ElementLevelSet *elLS,
			     int domainFlag,
			     int deg = 0,
			     Quadrature* q = NULL);
    static double H1NormSquare_DOF(DOFVector<double> *dof,
				   ElementLevelSet *elLS,
				   int domainFlag,
				   int deg = 0,
				   Quadrature* q = NULL);

    // ========================================================================
    //  Calculation of error between
    //     -> true solution (analytically given) and
    //     -> FE solution (dof vector).
    //
    //  Here, the H1-norm is the H1-seminorm, i.e. you get the "common"
    //  H1-norm by
    //    H1_Norm() = sqrt(L2_Norm^2() + H1_Seminorm^2()) .
    //
    //  Parameter \ref domainFlag:
    //      -3:  all elements in the domain with negative level set function
    //           values which do not intersect the zero level set
    //      -2 : all elements which intersect the domain with negative level set
    //           function values, in particular all (complete) elements which
    //           cut the zero level set
    //      -1 : domain with negative level set function values
    //       0 : all elements cut by the zero level set
    //       1 : complete mesh
    // ========================================================================
    static double L2Err(AbstractFunction<double, WorldVector<double> > *u,
			DOFVector<double> *uh,
			ElementLevelSet *elLS,
			int domainFlag,
			int relErr = 0,
			int deg = 2,
			Quadrature *q = NULL);
    static double H1Err(
			AbstractFunction<WorldVector<double>, WorldVector<double> > *grdU,
			DOFVector<double> *uh,
			ElementLevelSet *elLS,
			int domainFlag,
			int relErr = 0,
			int deg = 0,
			Quadrature *q = NULL);

    /**
     * Get absolute L2 error.
     * (If relative L2 error is calculated, the absolute L2 error is also
     *  calculated and stored in L2_err_abs).
     */
    inline static double getL2ErrAbs()
    {
      return L2_err_abs;
    }

    /**
     * Get absolute H1 error.
     * (If relative H1 error is calculated, the absolute H1 error is also
     *  calculated and stored in H1_err_abs).
     */
    inline static double getH1ErrAbs()
    {
      return H1_err_abs;
    }

    /**
     * Get L2 norm of solution u.
     * (If relative L2 error is calculated, the L2 norm of the solution u is 
     *  also calculated and stored in L2_u_norm).
     */
    inline static double getL2_U_Norm()
    {
      return L2_u_norm;
    }

    /**
     * Get H1 norm of solution u.
     * (If relative H1 error is calculated, the H1 norm of the solution u is 
     *  also calculated and stored in H1_u_norm).
     */
    inline static double getH1_U_Norm()
    {
      return H1_u_norm;
    }

  protected:
    static double Norm_IntNoBound(ElementNorm *elNorm,
				  ElementLevelSet *elLS,
				  Flag fillFlag,
				  int deg,
				  Quadrature* q);
    static double Norm_IntBound(ElementNorm *elNorm,
				ElementLevelSet *elLS,
				Flag fillFlag,
				int deg,
				Quadrature* q);
    static double Norm_Int(ElementNorm *elNorm,
			   ElementLevelSet *elLS,
			   Flag fillFlag,
			   int deg,
			   Quadrature* q);
    static double Norm_Bound(ElementNorm *elNorm,
			     ElementLevelSet *elLS,
			     Flag fillFlag,
			     int deg,
			     Quadrature* q);
    static double Norm_Complete(ElementNorm *elNorm,
				ElementLevelSet *elLS,
				Flag fillFlag,
				int deg,
				Quadrature* q);

    /// Calculate norm on subpolytope.
    static double calcSubPolNorm(ElInfo *elInfo,
				 SubPolytope *subPolytope,
				 ElementNorm *elNorm,
				 ScalableQuadrature *scalQuad,
				 const double &subPolFac = 1.0);

  protected:
    /// Absolute L2 error (last L2 error calculation !).
    static double L2_err_abs;

    /// L2 norm of correct solution u (last L2 error calculation !).
    static double L2_u_norm;

    /// Absolute H1 error (last H1 error calculation !).
    static double H1_err_abs;

    /// H1 norm of correct solution u (last H1 error calculation !).
    static double H1_u_norm;
  };

}

#endif  // AMDIS_CFE_NORMANDERRORFCTS_H
