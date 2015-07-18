#include "OperatorTerm.h"
#include "ElInfo.h"
#include "DOFVector.h"

namespace AMDiS 
{
  const Flag OperatorTerm::PW_CONST = 1;
  const Flag OperatorTerm::SYMMETRIC = 2;


  void OperatorTerm::setSymmetric(bool symm)
  {
    if (symm)
      properties.setFlag(SYMMETRIC);
    else
      properties.unsetFlag(SYMMETRIC);       
  }


  bool OperatorTerm::isSymmetric()
  {
    return properties.isSet(SYMMETRIC);
  }

} // end namespace AMDiS
