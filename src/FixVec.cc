#include "FixVec.h"

#include <cstdarg>

namespace AMDiS
{

  VectorOfFixVecs<DimVec<double>>* createAndInit(int dim, int size, ...)
  {
    va_list arg;
    va_start(arg, size);
    VectorOfFixVecs<DimVec<double>>* result =
                                   new VectorOfFixVecs<DimVec<double>>(dim, size, NO_INIT);
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
    double* result = new double[size];

    for (int i = 0; i < size; i++)
      result[i] = va_arg(arg, double);

    va_end(arg);
    return result;
  }

} // end namespace AMDiS
