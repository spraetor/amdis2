#include "AMDiS_base.hpp"

namespace AMDiS
{
  template <>
  EnumParser<MatrixTranspose>::EnumParser()
  {
    enumMap["NoTranspose"] = NoTranspose;
    enumMap["Transpose"] = Transpose;
    enumMap["ConjugateTranspose"] = ConjugateTranspose;
  }

  template <>
  EnumParser<Norm>::EnumParser()
  {
    enumMap["NO_NORM"] = NO_NORM;
    enumMap["H1_NORM"] = H1_NORM;
    enumMap["L2_NORM"] = L2_NORM;
  }

  template <>
  EnumParser<RefineCoarsenOperation>::EnumParser()
  {
    enumMap["NO_OPERATION"] = NO_OPERATION;
    enumMap["COARSE_RESTRICT"] = COARSE_RESTRICT;
    enumMap["COARSE_INTERPOL"] = COARSE_INTERPOL;
    enumMap["REFINE_INTERPOL"] = REFINE_INTERPOL;
  }

  template <>
  EnumParser<FirstOrderType>::EnumParser()
  {
    enumMap["GRD_PSI"] = GRD_PSI;
    enumMap["GRD_PHI"] = GRD_PHI;
  }

} // end namespace AMDiS
