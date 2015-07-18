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


#include <typeinfo>

#include "ElInfo3d.h"
#include "BasisFunction.h"
#include "Element.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "FiniteElemSpace.h"
#include "Flag.h"
#include "MacroElement.h"
#include "Mesh.h"
#include "Global.h"
#include "FixVec.h"
#include "DOFVector.h"

namespace AMDiS {

  double ElInfo3d::mat_d1_left_val[4][4] = {{1.0, 0.0, 0.0, 0.5}, 
					    {0.0, 0.0, 0.0, 0.5},
					    {0.0, 1.0, 0.0, 0.0},
					    {0.0, 0.0, 1.0, 0.0}};
  mtl::dense2D<double> ElInfo3d::mat_d1_left(mat_d1_left_val);

  
  double ElInfo3d::mat_d1_l0_right_val[4][4] = {{0.0, 0.0, 0.0, 0.5}, 
						{1.0, 0.0, 0.0, 0.5},
						{0.0, 0.0, 1.0, 0.0},
						{0.0, 1.0, 0.0, 0.0}};
  mtl::dense2D<double> ElInfo3d::mat_d1_l0_right(mat_d1_l0_right_val);

  double ElInfo3d::mat_d1_l12_right_val[4][4] = {{0.0, 0.0, 0.0, 0.5}, 
						 {1.0, 0.0, 0.0, 0.5},
						 {0.0, 1.0, 0.0, 0.0},
						 {0.0, 0.0, 1.0, 0.0}};
  mtl::dense2D<double> ElInfo3d::mat_d1_l12_right(mat_d1_l12_right_val);



  double ElInfo3d::mat_d4_left_val[35][35] = 
    {
      {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.273437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.0, 0.0, 0.023438},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.093750, 1.0, 0.468750, 0.0, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.031250, 0.0, 0.156250, 0.0, 0.0, 0.156250, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.703125, 0.0, 0.0, 0.0, 0.015625, 0.0, 0.140625, 0.015625, 0.0, 0.140625, 0.015625, 0.015625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.218750, 0.0, 0.0, 0.0, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.031250, 0.0, 0.156250, 0.093750, 0.0, 0.156250, 0.093750, 0.0, 0.0, 0.0, 0.0, 0.093750},
      {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.3125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.3125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0625},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0625},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.375},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.9375, 0.375, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9375, 0.375, 1.0, 0.0, 0.0, 0.0, 0.1875},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75}
    };
  mtl::dense2D<double> ElInfo3d::mat_d4_left(mat_d4_left_val);

  double ElInfo3d::mat_d4_l0_right_val[35][35] = 
    {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.273437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.023438},
      {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.218750, 0.0, 0.0, 0.0, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.031250, 0.0, 0.156250, 0.093750, 0.0, 0.156250, 0.093750, 0.0, 0.0, 0.0, 0.0, 0.093750},
      {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.703125, 0.0, 0.0, 0.0, 0.015625, 0.0, 0.140625, 0.015625, 0.0, 0.140625, 0.015625, 0.015625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.093750, 1.0, 0.468750, 0.0, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.031250, 0.0, 0.156250, 0.0, 0.0, 0.156250, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0625},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0625},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.3125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.3125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.375},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9375, 0.375, 1.0, 0.0, 0.0, 0.0, 0.1875},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.9375, 0.375, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75}
    };
  mtl::dense2D<double> ElInfo3d::mat_d4_l0_right(mat_d4_l0_right_val);

  double ElInfo3d::mat_d4_l12_right_val[35][35] = 
    {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.023437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.273438, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.023438, 0.0, 0.0, 0.0, 0.0, 0.023438},
      {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.218750, 0.0, 0.0, 0.0, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.031250, 0.0, 0.156250, 0.093750, 0.0, 0.156250, 0.093750, 0.0, 0.0, 0.0, 0.0, 0.093750},
      {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.703125, 0.0, 0.0, 0.0, 0.015625, 0.0, 0.140625, 0.015625, 0.0, 0.140625, 0.015625, 0.015625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.093750, 1.0, 0.468750, 0.0, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.0, 0.0, 0.031250, 0.031250, 0.0, 0.156250, 0.0, 0.0, 0.156250, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0625},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.0625, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0625},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.3125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875, 0.0, 0.0, 0.0625, 0.125, 0.0, 0.3125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.375},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.9375, 0.375, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1875},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9375, 0.375, 1.0, 0.0, 0.0, 0.0, 0.1875},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75},      
    };
  mtl::dense2D<double> ElInfo3d::mat_d4_l12_right(mat_d4_l12_right_val);



  void ElInfo3d::fillMacroInfo(const MacroElement * mel)
  {
    FUNCNAME_DBG("ElInfo3d::fillMacroInfo()");

    TEST_EXIT_DBG(mel)("No macro element given!\n");

    Element *nb;
    MacroElement *mnb;
    Flag fill_opp_coords;

    macroElement = const_cast<MacroElement*>(mel);
    element = const_cast<Element*>(mel->getElement());
    parent = NULL;
    level = 0;
    elType = const_cast<MacroElement*>(mel)->getElType();

    int vertices = mesh->getGeo(VERTEX);

    if (fillFlag.isSet(Mesh::FILL_COORDS) || 
	fillFlag.isSet(Mesh::FILL_DET) ||
	fillFlag.isSet(Mesh::FILL_GRD_LAMBDA))
      for (int i = 0; i < vertices; i++)
	coord[i] = mel->coord[i];

    int neighbours = mesh->getGeo(NEIGH);

    if (fillFlag.isSet(Mesh::FILL_OPP_COORDS) || 
        fillFlag.isSet(Mesh::FILL_NEIGH)) {

      fill_opp_coords.setFlags(fillFlag & Mesh::FILL_OPP_COORDS);
      for (int i = 0; i < neighbours; i++) {
	if ((mnb = const_cast<MacroElement*>(mel->getNeighbour(i)))) {
	  neighbour[i] = const_cast<Element*>(mel->getNeighbour(i)->getElement());
	  nb = const_cast<Element*>(neighbour[i]);
	  int k;
	  k = oppVertex[i] = mel->getOppVertex(i);

	  if (nb->getChild(0) && (k < 2)) {   /*make nb nearest element.*/
	    if (k == 1) {
	      neighbour[i] = const_cast<Element*>(nb->getChild(0));
	      nb = const_cast<Element*>(neighbour[i]);
	    } else {
	      neighbour[i] = const_cast<Element*>(nb->getChild(1));
	      nb = const_cast<Element*>(neighbour[i]);
	    }
	    k = oppVertex[i] = 3;
	    if (fill_opp_coords.isAnySet()) {
	      /* always edge between vertices 0 and 1 is bisected! */
	      if (mnb->getElement()->isNewCoordSet())
		oppCoord[i] = *(mnb->getElement()->getNewCoord());
	      else
		oppCoord[i] = (mnb->coord[0] + mnb->coord[1]) * 0.5;
	    }
	  } else {
	    if  (fill_opp_coords.isAnySet()) {
	      oppCoord[i] = mnb->coord[k];
	    }
	  }
	} else {
	  neighbour[i] = NULL;
	}
      }
    }

    if (fillFlag.isSet(Mesh::FILL_BOUND)) {
      for (int i = 0; i < element->getGeo(BOUNDARY); i++)
	boundary[i] = mel->getBoundary(i);     

      for (int i = 0; i < element->getGeo(PROJECTION); i++)
	projection[i] = mel->getProjection(i);
    }

    if (fillFlag.isSet(Mesh::FILL_ORIENTATION)) {
      WorldVector<WorldVector<double> > a;
      double s;

      for (int i = 0; i < 3; i++) {
	a[i] = mel->coord[i + 1];
	a[i] -= mel->coord[0];
      }

      s = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * a[2][0]
	+ (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * a[2][1]
	+ (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * a[2][2];

      if (s >= 0)
	orientation = 1;
      else
	orientation = -1;
    }
  }


  double ElInfo3d::calcGrdLambda(DimVec<WorldVector<double> >& grd_lam)
  {
    FUNCNAME("ElInfo3d::calcGrdLambda()");

    TEST_EXIT_DBG(dimOfWorld == 3)
      ("dim != dim_of_world ! use parametric elements!\n");

    std::vector<double> &e1 = tmpWorldVecs[0];
    std::vector<double> &e2 = tmpWorldVecs[1];
    std::vector<double> &e3 = tmpWorldVecs[2];

    testFlag(Mesh::FILL_COORDS);

    for (int i = 0; i < 3; i++) {
      e1[i] = coord[1][i] - coord[0][i];
      e2[i] = coord[2][i] - coord[0][i];
      e3[i] = coord[3][i] - coord[0][i];
    }

    double det = 
      e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) -
      e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) +
      e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);

    double adet = abs(det);

    if (adet < 1.0E-25) {
      MSG("abs(det) = %f\n",adet);
      for (int i = 0; i < 4; i++)
	for (int j = 0; j < 3; j++)
	  grd_lam[i][j] = 0.0;
    } else {
      det = 1.0 / det;
      /* (a_ij) = A^{-T} */

      grd_lam[1][0] = (e2[1] * e3[2] - e2[2] * e3[1]) * det;
      grd_lam[1][1] = (e2[2] * e3[0] - e2[0] * e3[2]) * det;
      grd_lam[1][2] = (e2[0] * e3[1] - e2[1] * e3[0]) * det;
      grd_lam[2][0] = (e1[2] * e3[1] - e1[1] * e3[2]) * det;
      grd_lam[2][1] = (e1[0] * e3[2] - e1[2] * e3[0]) * det;
      grd_lam[2][2] = (e1[1] * e3[0] - e1[0] * e3[1]) * det;
      grd_lam[3][0] = (e1[1] * e2[2] - e1[2] * e2[1]) * det;
      grd_lam[3][1] = (e1[2] * e2[0] - e1[0] * e2[2]) * det;
      grd_lam[3][2] = (e1[0] * e2[1] - e1[1] * e2[0]) * det;

      grd_lam[0][0] = -grd_lam[1][0] - grd_lam[2][0] - grd_lam[3][0];
      grd_lam[0][1] = -grd_lam[1][1] - grd_lam[2][1] - grd_lam[3][1];
      grd_lam[0][2] = -grd_lam[1][2] - grd_lam[2][2] - grd_lam[3][2];
    }

    return adet;
  }


  int ElInfo3d::worldToCoord(const WorldVector<double>& xy,
				   DimVec<double>* lambda) const
  {
    FUNCNAME("ElInfo::worldToCoord()");

    DimVec<WorldVector<double> > edge(mesh->getDim(), NO_INIT);
    WorldVector<double> x;
    double  x0, det, det0, det1, det2;
  
    static DimVec<double> vec(mesh->getDim(), NO_INIT);

    TEST_EXIT_DBG(lambda)("lambda must not be NULL\n");

    int dim = mesh->getDim();

    TEST_EXIT_DBG(dim == dimOfWorld)("dim!=dimOfWorld not yet implemented\n");
    
    /*  wir haben das gleichungssystem zu loesen: */
    /*       ( q1x q2x q3x)  (lambda1)     (qx)      */
    /*       ( q1y q2y q3y)  (lambda2)  =  (qy)      */
    /*       ( q1z q2z q3z)  (lambda3)     (qz)      */
    /*      mit qi=pi-p3, q=xy-p3                 */

    for (int j = 0; j < dimOfWorld; j++) {
      x0 = coord[dim][j];
      x[j] = xy[j] - x0;

      for (int i = 0; i < dim; i++)
	edge[i][j] = coord[i][j] - x0;
    }

    det =  edge[0][0] * edge[1][1] * edge[2][2]
      + edge[0][1] * edge[1][2] * edge[2][0]
      + edge[0][2] * edge[1][0] * edge[2][1]
      - edge[0][2] * edge[1][1] * edge[2][0]
      - edge[0][0] * edge[1][2] * edge[2][1]
      - edge[0][1] * edge[1][0] * edge[2][2];
    det0 =       x[0] * edge[1][1] * edge[2][2]
      +       x[1] * edge[1][2] * edge[2][0]
      +       x[2] * edge[1][0] * edge[2][1]
      -       x[2] * edge[1][1] * edge[2][0]
      -       x[0] * edge[1][2] * edge[2][1]
      -       x[1] * edge[1][0] * edge[2][2];
    det1 = edge[0][0] *       x[1] * edge[2][2]
      + edge[0][1] *       x[2] * edge[2][0]
      + edge[0][2] *       x[0] * edge[2][1]
      - edge[0][2] *       x[1] * edge[2][0]
      - edge[0][0] *       x[2] * edge[2][1]
      - edge[0][1] *       x[0] * edge[2][2];
    det2 = edge[0][0] * edge[1][1] *       x[2]
      + edge[0][1] * edge[1][2] *       x[0]
      + edge[0][2] * edge[1][0] *       x[1]
      - edge[0][2] * edge[1][1] *       x[0]
      - edge[0][0] * edge[1][2] *       x[1]
      - edge[0][1] * edge[1][0] *       x[2];
  
    if (abs(det) < DBL_TOL) {
      ERROR("det = %le; abort\n", det);

      for (int i = 0; i <= dim; i++)
	(*lambda)[i] = 1.0 / dim;

      return 0;
    }

    (*lambda)[0] = det0 / det;
    (*lambda)[1] = det1 / det;
    (*lambda)[2] = det2 / det;
    (*lambda)[3] = 1.0 - (*lambda)[0] - (*lambda)[1] - (*lambda)[2];
  
    int k = -1;
    double lmin = 0.0;

    for (int i = 0; i <= dim; i++) {
      if ((*lambda)[i] < -1.E-5) {
	if ((*lambda)[i] < lmin) {
	  k = i;
	  lmin = (*lambda)[i];
	}
      }
    }

    return k;
  }


  /****************************************************************************/
  /*   update EL_INFO structure after refinement (of some neighbours)	    */
  /****************************************************************************/

  void ElInfo3d::update()
  {
    FUNCNAME_DBG("ElInfo::update()");

    int neighbours = mesh->getGeo(NEIGH);
    int vertices = mesh->getGeo(VERTEX);
  
    if (fillFlag.isSet(Mesh::FILL_NEIGH) || fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {
      Tetrahedron *nb;
      Flag fill_opp_coords = fillFlag & Mesh::FILL_OPP_COORDS;
      
      for (int ineigh = 0; ineigh < neighbours; ineigh++) {
	if ((nb = dynamic_cast<Tetrahedron*>(const_cast<Element*>(neighbour[ineigh])))) {
	  int ov = oppVertex[ineigh];
	  if (ov < 2 && nb->getFirstChild()) {
	    if (fill_opp_coords != Flag(0)) {
	      int k = -1;
	      for (int j = 0; j < vertices; j++)
		if (element->getDof(j) == nb->getDof(1 - ov)) 
		  k = j;
	      
	      if (k == -1) {
		for (int j = 0; j < vertices; j++)
		  if (mesh->associated(element->getDof(j, 0), nb->getDof(1 - ov, 0)))
		    k = j;
	      }
	      TEST_EXIT_DBG(k >= 0)("neighbour dof not found\n");
	      
	      if (nb->isNewCoordSet())
		oppCoord[ineigh] = *(nb->getNewCoord());
	      else
		for (int j = 0; j < dimOfWorld; j++)
		  oppCoord[ineigh][j] = (oppCoord[ineigh][j] + coord[k][j]) / 2;
	    }
	    neighbour[ineigh] = dynamic_cast<Tetrahedron*>(const_cast<Element*>(nb->getChild(1-ov)));
	    oppVertex[ineigh] = 3;
	  }
	}
      }
    }
  }


  double ElInfo3d::getNormal(int face, WorldVector<double> &normal) const
  {
    FUNCNAME("ElInfo3d::getNormal()");

    double det = 0.0;

    WorldVector<double> e0, e1, e2;

    if (dimOfWorld == 3) {
      int i0 = (face + 1) % 4;
      int i1 = (face + 2) % 4;
      int i2 = (face + 3) % 4;

      for (int i = 0; i < dimOfWorld; i++) {
	e0[i] = coord[i1][i] - coord[i0][i];
	e1[i] = coord[i2][i] - coord[i0][i];
	e2[i] = coord[face][i] - coord[i0][i];
      }

      vectorProduct(e0, e1, normal);

      if ((e2 * normal) < 0.0)
	for (int i = 0; i < dimOfWorld; i++)
	  normal[i] = -normal[i];

      det = norm(&normal);
      TEST_EXIT_DBG(det > 1.e-30)("det = 0 on face %d\n", face);

      normal[0] /= det;
      normal[1] /= det;
      normal[2] /= det;
    } else {
      MSG("Not implemented for DIM_OF_WORLD = %d in 3D!\n", dimOfWorld);
    }

    return det;
  }


  void ElInfo3d::fillElInfo(int ichild, const ElInfo *elInfoOld)
  {
    FUNCNAME("ElInfo3d::fillElInfo()");

    TEST_EXIT_DBG(elInfoOld)("Missing old elInfo!\n");

    int ochild = 0;             /* index of other child = 1-ichild */
    int *cv = NULL;             /* cv = child_vertex[el_type][ichild] */
    const int (*cvg)[4] = NULL;     /* cvg = child_vertex[el_type] */
    int *ce;                    /* ce = child_edge[el_type][ichild] */
    Element *nb, *nbk;
    Element *elOld = elInfoOld->element;
    Flag fillFlag_local = elInfoOld->fillFlag;
    DegreeOfFreedom *dof;
    int ov = -1;
    Mesh *mesh = elInfoOld->getMesh();

    TEST_EXIT_DBG(elOld->getChild(0))("missing child?\n"); 

    element = const_cast<Element*>(elOld->getChild(ichild));
    macroElement = elInfoOld->macroElement;
    fillFlag = fillFlag_local;
    parent = elOld;
    level = elInfoOld->level + 1;
    iChild = ichild;
    int el_type_local = 0;
    try {
      el_type_local = (dynamic_cast<const ElInfo3d*>(elInfoOld))->getType();
    } catch (const std::bad_cast& e) {
      ERROR_EXIT("ElInfo is not of type ElInfo3D but %s!\n", typeid(*elInfoOld).name());
    }

    elType = (el_type_local + 1) % 3;

    TEST_EXIT_DBG(element)("missing child %d?\n", ichild);

    if (fillFlag_local.isAnySet()) {
      cvg = Tetrahedron::childVertex[el_type_local];
      cv = const_cast<int*>(cvg[ichild]);
      ochild = 1 - ichild;
    }

    if (fillFlag_local.isSet(Mesh::FILL_COORDS) || 
	fillFlag.isSet(Mesh::FILL_DET) ||
	fillFlag.isSet(Mesh::FILL_GRD_LAMBDA)) {
      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < dimOfWorld; j++)
	  coord[i][j] = elInfoOld->coord[cv[i]][j];
      }
      if (elOld->getNewCoord()) {
	coord[3] = *(elOld->getNewCoord());
      } else {
	for (int j = 0; j < dimOfWorld; j++)
	  coord[3][j] = (elInfoOld->coord[0][j] + elInfoOld->coord[1][j]) / 2;
      }
    }

    if (fillFlag_local.isSet(Mesh::FILL_NEIGH) || 
	fillFlag.isSet(Mesh::FILL_OPP_COORDS)) {

      FixVec<Element*, NEIGH> *neigh_local = &neighbour;
      const FixVec<Element*, NEIGH> *neigh_old = &elInfoOld->neighbour;

      Flag fill_opp_coords;
      fill_opp_coords.setFlags(fillFlag_local & Mesh::FILL_OPP_COORDS);
      
      // === nb[0] is other child ===
      
      if (elOld->getChild(0) &&  
	  (nb = const_cast<Element*>(elOld->getChild(ochild)))) {
	
	if (nb->getChild(0)) {         /* go down one level for direct neighbour */
	  if (fill_opp_coords.isAnySet()) {
	    if (nb->getNewCoord()) {
	      oppCoord[0]= *(nb->getNewCoord());
	    } else {
	      int k = cvg[ochild][1];
	      for (int j = 0; j < dimOfWorld; j++)
		oppCoord[0][j] = (elInfoOld->coord[ochild][j] + elInfoOld->coord[k][j]) / 2;
	    }
	  }
	  (*neigh_local)[0] = const_cast<Element*>(nb->getChild(1));
	  oppVertex[0] = 3;
	} else {
	  if (fill_opp_coords.isAnySet())
	    for (int j = 0; j < dimOfWorld; j++)
	      oppCoord[0][j] = elInfoOld->coord[ochild][j];

	  (*neigh_local)[0] = nb;
	  oppVertex[0] = 0;
	}
      } else {
	ERROR_EXIT("no other child");
	(*neigh_local)[0] = NULL;
      }


      /*----- nb[1], nb[2] are childs of old neighbours nb_old[cv[i]] ----------*/
      
      for (int i = 1; i < 3; i++) {
	if ((nb = const_cast<Element*>((*neigh_old)[cv[i]]))) {
	  TEST_EXIT_DBG(nb->getChild(0))("nonconforming triangulation\n");

	  int k;
	  for (k = 0; k < 2; k++) { /* look at both childs of old neighbour */
	    nbk = const_cast<Element*>(nb->getChild(k));

	    if (nbk->getDof(0) == elOld->getDof(ichild)) {
	      /* opp. vertex */
	      dof = const_cast<DegreeOfFreedom*>(nb->getDof(elInfoOld->oppVertex[cv[i]])); 
	      
	      if (dof == nbk->getDof(1)) {
		ov = 1;
		if (nbk->getChild(0)) {
		  if (fill_opp_coords.isAnySet()) {
		    if (nbk->getNewCoord()) {
		      oppCoord[i] = *(nbk->getNewCoord());
		    } else {
		      for (int j = 0; j < dimOfWorld; j++)
			oppCoord[i][j] = (elInfoOld->oppCoord[cv[i]][j] + 
					   elInfoOld->coord[ichild][j]) / 2;		      
		    }
		  }
		  (*neigh_local)[i] = nbk->getChild(0);
		  oppVertex[i] = 3;
		  break;
		}
	      } else {
		if (dof != nbk->getDof(2)) { 
		  ov = -1; 
		  break; 
		}
		ov = 2;
	      }

	      if (fill_opp_coords.isAnySet())
		for (int j = 0; j < dimOfWorld; j++)
		  oppCoord[i][j] = elInfoOld->oppCoord[cv[i]][j];

	      (*neigh_local)[i] = nbk;
	      oppVertex[i] = ov;

	      break;
	    }
	    
	  } /* end for k */


	  // === Test if periodic. ===

	  if (k == 2 || ov == -1) {

	    // Look at both childs of old neighbour.
	    for (k = 0; k < 2; k++) {
	      nbk = const_cast<Element*>(nb->getChild(k));
	      if (nbk->getDof(0) == elOld->getDof(ichild) ||
		  mesh->associated(nbk->getDof(0, 0), elOld->getDof(ichild, 0))) {

		// opp. vertex 
		dof = const_cast<DegreeOfFreedom*>(nb->getDof(elInfoOld->oppVertex[cv[i]])); 
		
		if (dof == nbk->getDof(1) || 
		    mesh->associated(dof[0], nbk->getDof(1, 0))) {
		  ov = 1;
		  if (nbk->getChild(0)) {
		    if (fill_opp_coords.isAnySet()) {
		      if (nbk->getNewCoord()) {
			oppCoord[i] = *(nbk->getNewCoord());
		      } else {
			for (int j = 0; j < dimOfWorld; j++)
			  oppCoord[i][j] = (elInfoOld->oppCoord[cv[i]][j] + 
					    elInfoOld->coord[ichild][j]) / 2;
		      }
		    }
		    (*neigh_local)[i] = nbk->getChild(0);
		    oppVertex[i] = 3;
		    break;
		  }
		} else {
		  TEST_EXIT_DBG(dof == nbk->getDof(2) || 
				mesh->associated(dof[0], nbk->getDof(2, 0)))
		    ("opp_vertex not found\n");
		  ov = 2;
		}

		if (fill_opp_coords.isAnySet())
		  for (int j = 0; j < dimOfWorld; j++)
		    oppCoord[i][j] = elInfoOld->oppCoord[cv[i]][j];

		(*neigh_local)[i] = nbk;
		oppVertex[i] = ov;
		break;
	      }
	      
	    } /* end for k */

	    TEST_EXIT_DBG(k < 2)
	      ("Fill %d child of element %d (0-Level: %d): Child not found with vertex!\n",
	       ichild, elOld->getIndex(), elInfoOld->getMacroElement()->getIndex());
	  }
	} else {
	  (*neigh_local)[i] = NULL;
	}
      }  /* end for i */
      
      
      /*----- nb[3] is old neighbour neigh_old[ochild] ------------------------*/
      
      if (((*neigh_local)[3] = (*neigh_old)[ochild])) {
	oppVertex[3] = elInfoOld->oppVertex[ochild];

	if (fill_opp_coords.isAnySet())
	  for (int j = 0; j < dimOfWorld; j++)
	    oppCoord[3][j] = elInfoOld->oppCoord[ochild][j];
      }
    }

    if (fillFlag_local.isSet(Mesh::FILL_BOUND)) {
      for (int i = 0; i < 3; i++)
	boundary[10 + i] = elInfoOld->getBoundary(10 + cv[i]);
      
      boundary[13] = elInfoOld->getBoundary(4);
      
      boundary[0] = INTERIOR;
      boundary[1] = elInfoOld->getBoundary(cv[1]);
      boundary[2] = elInfoOld->getBoundary(cv[2]);
      boundary[3] = elInfoOld->getBoundary(ochild);
      
      int geoFace = mesh->getGeo(FACE);

      ce = const_cast<int*>(Tetrahedron::childEdge[el_type_local][ichild]);
      for (int iedge = 0; iedge < 4; iedge++)
	boundary[geoFace + iedge] = elInfoOld->getBoundary(geoFace + ce[iedge]);      
      for (int iedge = 4; iedge < 6; iedge++)
	boundary[geoFace + iedge] = elInfoOld->getBoundary(5 - cv[iedge - 3]);

      if (elInfoOld->getProjection(0) &&
	  elInfoOld->getProjection(0)->getType() == VOLUME_PROJECTION) {
	
	projection[0] = elInfoOld->getProjection(0);      
      } else { // boundary projection
	projection[0] = NULL;
	projection[1] = elInfoOld->getProjection(cv[1]);
	projection[2] = elInfoOld->getProjection(cv[2]);
	projection[3] = elInfoOld->getProjection(ochild);
	
	for (int iedge = 0; iedge < 4; iedge++)
	  projection[geoFace + iedge] = elInfoOld->getProjection(geoFace + ce[iedge]);
	for (int iedge = 4; iedge < 6; iedge++)
	  projection[geoFace + iedge] = elInfoOld->getProjection(5 - cv[iedge - 3]);
      }
    }

    
    if (fillFlag.isSet(Mesh::FILL_ORIENTATION)) {
      orientation = 
	(dynamic_cast<ElInfo3d*>(const_cast<ElInfo*>(elInfoOld)))->orientation 
	* Tetrahedron::childOrientation[el_type_local][ichild];
    }
  }


  mtl::dense2D<double>& ElInfo3d::getSubElemCoordsMat(int degree) const
  {
    FUNCNAME("ElInfo3d::getSubElemCoordsMat()");

    using namespace mtl;

    if (subElemMatrices[degree].count(std::make_pair(refinementPathLength, refinementPath)) == 0) {
      switch (degree) {
      case 1:
	{
	  dense2D<double> mat(4, 4), tmpMat(4, 4);
	  mat = 1;
	  
	  for (int i = 0; i < refinementPathLength; i++) {
	    if (refinementPath & (1 << i)) {
	      if ((level + i) % 3 == 0)
		tmpMat = mat * mat_d1_l0_right;
	      else
		tmpMat = mat * mat_d1_l12_right;

	      mat = tmpMat;
	    } else  {
	      tmpMat = mat * mat_d1_left;
	      mat = tmpMat;
	    }
	  }

	  subElemMatrices[degree][std::make_pair(refinementPathLength, refinementPath)] = mat;  
	}
	break;
      case 4:
	{
	  dense2D<double> mat(35, 35), tmpMat(35, 35);
	  mat = 1;
	  
	  for (int i = 0; i < refinementPathLength; i++) {
	    if (refinementPath & (1 << i)) {
	      if ((level + i) % 3 == 0)
		tmpMat = mat * mat_d4_l0_right;
 	      else
 		tmpMat = mat * mat_d4_l12_right;

	      mat = tmpMat;
	    } else  {
	      tmpMat = mat * mat_d4_left;
	      mat = tmpMat;
	    }
	  }

	  subElemMatrices[degree][std::make_pair(refinementPathLength, refinementPath)] = mat;  
	}
	break;	
      default:
	ERROR_EXIT("Not supported for basis function degree: %d\n", degree);
      }
    }

    return subElemMatrices[degree][std::make_pair(refinementPathLength, refinementPath)];
  }

}
