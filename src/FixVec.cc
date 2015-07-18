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


#include <stdarg.h>
#include "FixVec.h"

namespace AMDiS {

  VectorOfFixVecs<DimVec<double> > *createAndInit(int dim, int size, ...)
  {
    va_list arg;
    va_start(arg, size);
    VectorOfFixVecs<DimVec<double> > *result = 
      new VectorOfFixVecs<DimVec<double> >(dim, size, NO_INIT);
    for (int i = 0; i < size; i++)
      for (int j = 0; j < dim + 1; j++)
	(*result)[i][j] = va_arg(arg, double);
    va_end(arg);
    return result;
  }

  double* createAndInitArray(int size, ...)
  {
    va_list arg;
    va_start(arg, size);
    double *result = new double[size];
  
    for (int i = 0; i < size; i++)
      result[i] = va_arg(arg, double);

    va_end(arg);
    return result;
  }

}
