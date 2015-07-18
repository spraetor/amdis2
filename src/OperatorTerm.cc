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


#include "OperatorTerm.h"
#include "ElInfo.h"
#include "DOFVector.h"

namespace AMDiS {

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

}
