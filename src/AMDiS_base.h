/** \file AMDiS_base.h */

#pragma once

#include <string>		// std::string
#include <functional>		// std::binary_function (deprecated)
#include <array>

#ifdef _MSC_VER
#include <io.h>			// _access
#else
#include <unistd.h>
#endif


#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>		// boost::algorithm::trim

namespace AMDiS 
{
  static constexpr int amdisRevisionNumber = 1700; // TODO: how to update this value

  /// Used by matrix vector multiplication
  typedef enum { NoTranspose,
            		 Transpose,
            		 ConjugateTranspose } MatrixTranspose;

  /// Speciefies the norm used by Estimator.
  typedef enum { NO_NORM = 0, 
            		 H1_NORM = 1, 
            		 L2_NORM = 2 } Norm;
		 

  /// Specifies the type of a FirstOrderTerm 
  enum FirstOrderType {
    GRD_PSI,
    GRD_PHI
  };
  
  struct DofIndex 
  {
    using size_type = signed long;
  };
  
//   std::ostream& operator<<(std::ostream& os, const DofIndex& di);  
//   std::istream& operator>>(std::istream& is, DofIndex& di);

  /// Datatype for degrees of freedom 
  using DegreeOfFreedom = DofIndex::size_type;

  /// Defines type for a vector of DOF pointers.
  using DofContainer = std::vector<const DegreeOfFreedom*>;

  using DofContainerSet = std::set<const DegreeOfFreedom*>;

  /// Defines a type for global edge identification via its DOFs. 
  // TODO: replace by std::array<DegreeOfFreedom, 2>
  using DofEdge = std::pair<DegreeOfFreedom, DegreeOfFreedom>;

  /// Defines a tzpe for global face identiication via its DOFs.
  using DofFace = std::array<DegreeOfFreedom, 3>;

  // ===== some simple template functions ====================================== 
  
  /// Content comparision of two pointers. Used e.g. for find_if
  template <class T>
  struct comparePtrContents : public std::binary_function<T*, T*, bool>
  {
    bool operator()(T* a, T* b) const 
    {
      return (*a == *b);
    }
  };

  
  /// check for file existence
  inline bool file_exists(const std::string filename)
  {
#ifdef _MSC_VER
    return _access(filename.c_str(), 0) == 0;
#else
    return access(filename.c_str(), F_OK) == 0;
#endif
  }

  
  /// trim std::string
  inline std::string trim(const std::string& oldStr)
  {
    std::string swap(oldStr);
    boost::algorithm::trim(swap);
    return swap;
  }
}
