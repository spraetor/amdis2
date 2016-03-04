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


#include "NormEps.h"
#include "MatrixVectorOperations.h"

namespace reinit
{

  double NormEps::eps = 0.0;

  double
  NormEps::calcNormEps(const WorldVector<double>& x)
  {
    FUNCNAME("NormEps::calcNormEps()");

    static double result;

    TEST_EXIT(eps > 1.e-15)("illegal eps for norm regularization !\n");
    result = x*x + eps*eps;

    return(sqrt(result));
  }

  double
  NormEps::calcNormEps2(const WorldVector<double>& x)
  {
    FUNCNAME("NormEps::calcNormEps2()");

    static double result;

    TEST_EXIT(eps > 1.e-15)("illegal eps for norm regularization !\n");
    result = x*x + eps*eps;

    return(result);
  }

}
