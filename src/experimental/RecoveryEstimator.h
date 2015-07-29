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



/** \file RecoveryEstimator.h */

#ifndef AMDIS_MYRECOVERYESTIMATOR_H
#define AMDIS_MYRECOVERYESTIMATOR_H

#include "Estimator.h"
#include "Recovery.h"

namespace AMDiS {

  /// Error estimator using the recovery gradient of the finite element solution.
  class RecoveryEstimator : public Estimator
  {
  public:
    /// Creator class.
    class Creator : public EstimatorCreator
    {
    public:
      Creator() : EstimatorCreator() {}

      virtual ~Creator() {}

      /// Returns a new Estimator object.
      virtual Estimator* create()
      {
	return new RecoveryEstimator(name, uh, row);
      }
    };

    /// constructor
    RecoveryEstimator(std::string name, DOFVector<double> *uh_, int r = -1);

    /// destructor.
    virtual ~RecoveryEstimator() {}

    /// implements \ref Estimator::init(double).
    virtual void init(double ts);

    /// implements \ref Estimator::estimateElement(ElInfo*, DualElInfo*).
    virtual void estimateElement(ElInfo *elInfo, DualElInfo *dualElInfo = NULL);

    /// implements \ref Estimator::exit(bool).
    virtual void exit(bool output = true);

    /// Sets uh.
    inline void setUh(DOFVector<double> *uh_)
    {
      uh = uh_;
    }

    /// Sets f.
    inline void setFct(AbstractFunction<double, WorldVector<double> > *fct)
    {
      f_vec = fct;
    }

    ///
    inline void setFct(AbstractFunction<double, double> *fct)
    {
      f_scal = fct;
    }

    /// Sets auxiliar vector.
    inline void setAuxVec(DOFVector<double> *uh)
    {
      aux_vec = uh;
    }

    /// Gets recovery gradient.
    inline DOFVector<WorldVector<double> >* getRecGrd()
    {
      return rec_grd;
    }

    /// Gets higher-order recovery solution.
    inline DOFVector<double>* getRecUh()
    {
      return rec_uh;
    }


  protected:   
    /// finite element solution
    DOFVector<double> *uh;

    /// absolute or relative error?
    int relative;

    /// constant for scaling the estimator
    double C;

    /// recovery method
    int method;

    /// Working finite element space
    const FiniteElemSpace *feSpace;

    /// Degree of corresponding basic functions
    int degree;

    /// Basis functions for recovery vector.
    const BasisFunction *rec_basFcts;

    /// Recovery gradient
    DOFVector<WorldVector<double> > *rec_grd;

    /// Higher-order recovery solution
    DOFVector<double> *rec_uh;

    /// Diffusion coefficient (for flux recovery)
    AbstractFunction<double, WorldVector<double> > *f_vec;
    AbstractFunction<double, double> *f_scal;

    /// auxiliar vector
    DOFVector<double> *aux_vec;

    /// Recovery structure.
    Recovery *rec_struct;

    /// Number of quadrature points.
    int nPoints;

    Quadrature *quad;
    FastQuadrature *quadFast, *rec_quadFast;

    /// Basis functions 
    const BasisFunction *basFcts;
    
    double h1Norm2;

    WorldVector<double> quad_pt;
    mtl::dense_vector<WorldVector<double> > grdAtQP;
    mtl::dense_vector<WorldVector<double> > recoveryGrdAtQP;
    mtl::dense_vector<double> uhAtQP;
    mtl::dense_vector<double> recoveryUhAtQP;
  };

}

#endif
