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



/** \file expressions.h */

#ifndef AMDIS_EXPRESSIONS_H
#define AMDIS_EXPRESSIONS_H

#include "expr_traits.hpp"	// collection of traits necessary to construct expressions
#include "value_expr.hpp"	// constants
#include "coords_expr.hpp"	// coordinates / normal vectors
#include "valueOf.hpp"		// value of DOFVector at QPs
#include "gradientOf.hpp"	// gradient of DOFVector at QPs
#include "hessianOf.hpp"	// second derivatives of DOFVector at QPs

#include "add_expr.hpp"		// add two expressions
#include "mult_expr.hpp"	// multiply two expressions
#include "functor_expr.hpp"	// apply a functor with 1/2/3 arguments to expressions
#include "cmath_expr.hpp"	// apply a cmath function to expressions
#include "vec_functors.hpp"	// apply a vector function to expressions

#endif // AMDIS_EXPRESSIONS_H
