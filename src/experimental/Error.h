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



/** \file Error.h */

#ifndef AMDIS_ERROR_H
#define AMDIS_ERROR_H

#include "Global.h"
#include "Mesh.h"
#include "AMDiS_fwd.h"


// NOTE: Evneutell rausschmeissen
namespace AMDiS {

  /** \ingroup Common
   * \brief
   * True error calculator
   */
  template<typename T>
  class Error
  {
  public:  
    static void setWriteLeafData() 
    { 
      writeInLeafData = true; 
    }
    
    static void unsetWriteLeafData() 
    { 
      writeInLeafData = false; 
    }
    
    static bool writeLeafData() 
    { 
      return writeInLeafData; 
    }
  
    static double maxErrAtQp(const AbstractFunction<T, WorldVector<double> >& u,
			     const DOFVector<T>& uh,
			     const Quadrature* q);

    static double H1Err(const AbstractFunction<WorldVector<T>, WorldVector<double> >& grdU,
			const DOFVector<T>& uh,
			int relErr,
			double* max,
			bool writeLeafData = false,
			int comp = 0);

    static double L2Err(const AbstractFunction<T, WorldVector<double> >& u,
			const DOFVector<T>& uh,
			int relErr,
			double* max,
			bool writeLeafData = false,
			int comp = 0);

    static double L2Err_ElementWise(const AbstractFunction<T, WorldVector<double> >& u,
				    const DOFVector<T>& uh,
				    double* max,
				    bool writeLeafData,
				    int comp,
				    int level,
				    std::map<int, double>& estMap);
  public:
    static T errUFct(const DimVec<double>& lambda);

    static WorldVector<T> grdErrUFct(const DimVec<double>& lambda);

    class AbstrFctErrU  : public AbstractFunction<T, DimVec<double> >
    { 
    public:
      AbstrFctErrU() : AbstractFunction<T, DimVec<double> >(0) {}

      inline T operator()(const DimVec<double> & x) const 
      { 
	return Error<T>::errUFct(x);
      }
    };

    static AbstrFctErrU errU;

    class AbstrFctGrdErrU  : public AbstractFunction<WorldVector<T>, DimVec<double> >
    { 
    public:
      AbstrFctGrdErrU() : AbstractFunction<WorldVector<T>, DimVec<double> >(0) {}

      inline WorldVector<T> operator()(const DimVec<double> & x) const 
      { 
	return Error<T>::grdErrUFct(x);
      }
    };

    static AbstrFctGrdErrU   grdErrU;

  private:
    static ElInfo* elinfo;
    /*   static const Parametric* el_parametric; */
    static const FastQuadrature* quadFast;
    static const AbstractFunction<T, WorldVector<double> >* pU;
    static const AbstractFunction<WorldVector<T>, WorldVector<double> >* pGrdU;
    static const BasisFunction* basFct;
    static const DOFVector<T>* errUh;
    static bool writeInLeafData;
    static int component;
  };

  template<typename T> ElInfo* Error<T>::elinfo = NULL;
  template<typename T> const FastQuadrature* Error<T>::quadFast = NULL;
  template<typename T> const AbstractFunction<T, WorldVector<double> >* Error<T>::pU = NULL;
  template<typename T> const AbstractFunction<WorldVector<T>, WorldVector<double> >* Error<T>::pGrdU = NULL;
  template<typename T> const BasisFunction* Error<T>::basFct = NULL;
  template<typename T> const DOFVector<T>* Error<T>::errUh = NULL;
  template<typename T> typename Error<T>::AbstrFctErrU Error<T>::errU;
  template<typename T> typename Error<T>::AbstrFctGrdErrU Error<T>::grdErrU;
  template<typename T> bool Error<T>::writeInLeafData = false;
  template<typename T> int Error<T>::component = 0;

}

#include "Error.hh"

#endif  // AMDIS_ERROR_H



