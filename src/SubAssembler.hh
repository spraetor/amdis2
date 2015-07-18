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


#include "Quadrature.h"
#include "FiniteElemSpace.h"
#include "ElInfo.h"

namespace AMDiS {

  template <class T>
  void SubAssembler::getVectorAtQPs(DOFVectorBase<T>* vec, 
				    const ElInfo* elInfo,
				    Quadrature *quad,
				    DenseVector<T>& vecAtQPs)
  {
    FUNCNAME_DBG("SubAssembler::getVectorAtQPs()");

    TEST_EXIT_DBG(vec)("No DOF vector!\n");
    TEST_EXIT_DBG(elInfo->getMesh() == vec->getFeSpace()->getMesh())
      ("Vector and FE space do not fit together!\n");

    Quadrature *localQuad = quad ? quad : quadrature;
    vecAtQPs.change_dim(localQuad->getNumPoints());

    if (cachedValuesAtQPs[vec] && cachedValuesAtQPs[vec]->valid) {
      vecAtQPs = boost::any_cast<DenseVector<T>& >(cachedValuesAtQPs[vec]->values);
      return;
    }
    
    boost::any swap;
    if (!cachedValuesAtQPs[vec]) {
      swap = DenseVector<T>();
      cachedValuesAtQPs[vec] = new ValuesAtQPs(swap, false);
    }
    boost::any_cast<DenseVector<T>& >(cachedValuesAtQPs[vec]->values).change_dim(localQuad->getNumPoints());

    DenseVector<T>& values = boost::any_cast<DenseVector<T>& >(cachedValuesAtQPs[vec]->values);
    bool sameFeSpaces = 
      vec->getFeSpace() == rowFeSpace || vec->getFeSpace() == colFeSpace;

    if (opt && !quad && sameFeSpaces) {
      const BasisFunction *psi = rowFeSpace->getBasisFcts();
      const BasisFunction *phi = colFeSpace->getBasisFcts();
      if (vec->getFeSpace()->getBasisFcts() == psi)
	psiFast = updateFastQuadrature(psiFast, psi, INIT_PHI);
      else if(vec->getFeSpace()->getBasisFcts() == phi)
	phiFast = updateFastQuadrature(phiFast, phi, INIT_PHI);      
    }

    // calculate new values
    const BasisFunction *basFcts = vec->getFeSpace()->getBasisFcts();

    if (opt && !quad && sameFeSpaces) {
      if (psiFast->getBasisFunctions() == basFcts) {
	vec->getVecAtQPs(elInfo, NULL, psiFast, values);
      } else if (phiFast->getBasisFunctions() == basFcts) {
	vec->getVecAtQPs(elInfo, NULL, phiFast, values);
      } else {
	vec->getVecAtQPs(elInfo, localQuad, NULL, values);
      }
    } else {
      vec->getVecAtQPs(elInfo, localQuad, NULL, values);
    }

    cachedValuesAtQPs[vec]->valid = true;
    cachedValuesAtQPs[vec]->quad = localQuad;
    vecAtQPs = values;
  }


  template <class T>
  void SubAssembler::getGradientsAtQPs(DOFVectorBase<T>* vec,
				       const ElInfo* elInfo,
				       Quadrature *quad,
				       DenseVector<typename GradientType<T>::type>& grdAtQPs)
  {
    FUNCNAME_DBG("SubAssembler::getGradientsAtQPs()");

    TEST_EXIT_DBG(vec)("No DOF vector!\n");
    Quadrature *localQuad = quad ? quad : quadrature;
    grdAtQPs.change_dim(localQuad->getNumPoints());

    if (cachedGradientsAtQPs[vec] && cachedGradientsAtQPs[vec]->valid) {
      grdAtQPs = 
	boost::any_cast<DenseVector<typename GradientType<T>::type>& >(cachedGradientsAtQPs[vec]->values);
      return;
    }

    boost::any swap;
    if (!cachedGradientsAtQPs[vec]) {
      swap = DenseVector<typename GradientType<T>::type>();
      cachedGradientsAtQPs[vec] = new ValuesAtQPs(swap, false);
    }

    boost::any_cast<DenseVector<typename GradientType<T>::type>& >(cachedGradientsAtQPs[vec]->values).change_dim(localQuad->getNumPoints());

    DenseVector<typename GradientType<T>::type>& values = 
      boost::any_cast<DenseVector<typename GradientType<T>::type>& >(cachedGradientsAtQPs[vec]->values);

    const BasisFunction *psi = rowFeSpace->getBasisFcts();
    const BasisFunction *phi = colFeSpace->getBasisFcts();

    bool sameFeSpaces =
      vec->getFeSpace() == rowFeSpace || vec->getFeSpace() == colFeSpace;

    if (opt && !quad && sameFeSpaces) {
      if (vec->getFeSpace()->getBasisFcts() == psi)
	psiFast = updateFastQuadrature(psiFast, psi, INIT_GRD_PHI);
      else if(vec->getFeSpace()->getBasisFcts() == phi)
	phiFast = updateFastQuadrature(phiFast, phi, INIT_GRD_PHI);
    }

    // calculate new values
    const BasisFunction *basFcts = vec->getFeSpace()->getBasisFcts();

    if (opt && !quad && sameFeSpaces) {
      if (psiFast->getBasisFunctions() == basFcts)
	vec->getGrdAtQPs(elInfo, NULL, psiFast, values);
      else if (phiFast->getBasisFunctions() == basFcts)
	vec->getGrdAtQPs(elInfo, NULL, phiFast, values);
      else
	vec->getGrdAtQPs(elInfo, NULL, phiFast, values);
    } else {
      vec->getGrdAtQPs(elInfo, localQuad, NULL, values);
    }

    cachedGradientsAtQPs[vec]->valid = true;
    grdAtQPs = values;
  }

  template <class T>
  void SubAssembler::getDerivativeAtQPs(DOFVectorBase<T>* vec,
					const ElInfo* elInfo,
					Quadrature *quad,
					int comp,
					DenseVector<T>& grdAtQPs)
  {
    FUNCNAME_DBG("SubAssembler::getGradientsAtQPs()");

    TEST_EXIT_DBG(vec)("No DOF vector!\n");

    Quadrature *localQuad = quad ? quad : quadrature;
    grdAtQPs.change_dim(localQuad->getNumPoints());

    if (cachedGradientsAtQPs[vec] && cachedGradientsAtQPs[vec]->valid) {      
      DenseVector<typename GradientType<T>::type> tmp = boost::any_cast<DenseVector<typename GradientType<T>::type>& >(cachedGradientsAtQPs[vec]->values);
      for (size_t iq = 0; iq < num_rows(tmp); iq++)
	grdAtQPs[iq] = tmp[iq][comp];
      return;
    }

    const BasisFunction *psi = rowFeSpace->getBasisFcts();
    const BasisFunction *phi = colFeSpace->getBasisFcts();

    bool sameFeSpaces =
      vec->getFeSpace() == rowFeSpace || vec->getFeSpace() == colFeSpace;

    if (opt && !quad && sameFeSpaces) {
      if (vec->getFeSpace()->getBasisFcts() == psi)
	psiFast = updateFastQuadrature(psiFast, psi, INIT_GRD_PHI);
      else if(vec->getFeSpace()->getBasisFcts() == phi)
	phiFast = updateFastQuadrature(phiFast, phi, INIT_GRD_PHI);
    }

    // calculate new values
    const BasisFunction *basFcts = vec->getFeSpace()->getBasisFcts();

    if (opt && !quad && sameFeSpaces) {
      if (psiFast->getBasisFunctions() == basFcts)
	vec->getDerivativeAtQPs(elInfo, NULL, psiFast, comp, grdAtQPs);
      else if (phiFast->getBasisFunctions() == basFcts)
	vec->getDerivativeAtQPs(elInfo, NULL, phiFast, comp, grdAtQPs);
      else
	vec->getDerivativeAtQPs(elInfo, NULL, phiFast, comp, grdAtQPs);
    } else {
      vec->getDerivativeAtQPs(elInfo, localQuad, NULL, comp, grdAtQPs);
    }
  }

} // end namespace AMDiS
