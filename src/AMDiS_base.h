/** \file AMDiS_base.h */

#pragma once

#include <string>		// std::string
#include <map>
#include <array>
#include <vector>
#include <set>

namespace AMDiS
{
  static constexpr int amdisRevisionNumber = 1700; // TODO: how to update this value
  
  template <class T>
  struct EnumParser
  {
    EnumParser() {}
    void operator()(std::string const& valStr, T& value)
    {
      auto it = enumMap.find(valStr);
      if (it != enumMap.end())
	value = it->second;
    }
    
    T operator()(std::string const& valStr)
    {
      T result; nullify(result);
      this->operator()(valStr, result);
      return result;
    }
    
  private:
    std::map<std::string, T> enumMap;
  };

  /// Used by matrix vector multiplication
  enum MatrixTranspose 
  { 
    NoTranspose,
    Transpose,
    ConjugateTranspose
  };
  template <> EnumParser<MatrixTranspose>::EnumParser();

  /// Speciefies the norm used by Estimator.
  enum Norm
  { 
    NO_NORM = 0,
    H1_NORM = 1,
    L2_NORM = 2
  };
  template <> EnumParser<Norm>::EnumParser();

  /// Specifies which operation should be done after coarsening
  enum RefineCoarsenOperation
  { 
    NO_OPERATION = 0,
    COARSE_RESTRICT = 1,
    COARSE_INTERPOL = 2,
    REFINE_INTERPOL = 4
  };
  template <> EnumParser<RefineCoarsenOperation>::EnumParser();

  /// Specifies the type of a FirstOrderTerm
  enum FirstOrderType 
  { 
    GRD_PSI,
    GRD_PHI
  };
  template <> EnumParser<FirstOrderType>::EnumParser();

  struct DofIndex
  {
    using size_type = signed long;
  };

  /// Datatype for degrees of freedom
  using DegreeOfFreedom = DofIndex::size_type;

  /// Defines type for a vector of DOF pointers.
  using DofContainer = std::vector<DegreeOfFreedom const*>;

  using DofContainerSet = std::set<DegreeOfFreedom const*>;

  /// Defines a type for global edge identification via its DOFs.
  // TODO: replace by std::array<DegreeOfFreedom, 2>
  using DofEdge = std::pair<DegreeOfFreedom, DegreeOfFreedom>;

  /// Defines a tzpe for global face identiication via its DOFs.
  using DofFace = std::array<DegreeOfFreedom, 3>;
}
