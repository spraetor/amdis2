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


#include "ElInfoStack.h"
#include "Mesh.h"

namespace AMDiS {
  
  ElInfoStack::ElInfoStack(Mesh *mesh) 
    : mesh_(mesh),
      stackPosition_(-1)
  {
    elInfoStack_.resize(20);
    for (int i = 0; i < 20; i++) {
      elInfoStack_[i] = mesh_->createNewElInfo();
    }
  }

  
  ElInfoStack::~ElInfoStack()
  {
    for (int i = 0; i < static_cast<int>(elInfoStack_.size()); i++) {
      delete elInfoStack_[i];
    }     
  }
  

  ElInfo* ElInfoStack::getNextElement()
  {
    // Check if the stack if large enough. If not, the stack is enlarged
    // and new elements are created.
    if (stackPosition_ + 1 >= static_cast<int>(elInfoStack_.size())) {
      int oldSize = elInfoStack_.size();
      elInfoStack_.resize(oldSize + oldSize / 2);
      for (int i = oldSize; i < static_cast<int>(elInfoStack_.size()); i++) {
	elInfoStack_[i] = mesh_->createNewElInfo();
      }
    }

    stackPosition_++;
    
    return elInfoStack_[stackPosition_];
  }

  
  void ElInfoStack::getBackElement()
  {
    TEST_EXIT_DBG(stackPosition_ >= 0)("Invalid stack position!\n");

    stackPosition_--;
  }

  
  ElInfo* ElInfoStack::getCurrentElement()
  {
    TEST_EXIT_DBG(stackPosition_ >= 0)("Invalid stack position!\n");
    
    return elInfoStack_[stackPosition_];
  }
  
} // end namespace AMDiS
