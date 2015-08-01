/** \file expressions.h */

#pragma once

#include "expr_traits.hpp"	// collection of traits necessary to construct expressions
#include "value_expr.hpp"	// constants
#include "coords_expr.hpp"	// coordinates / normal vectors
#include "valueOf.hpp"		// value of DOFVector at QPs
#include "gradientOf.hpp"	// gradient of DOFVector at QPs
#include "hessianOf.hpp"	// second derivatives of DOFVector at QPs

#include "add_expr.hpp"		// add two expressions
#include "mult_expr.hpp"	// multiply two expressions
#include "functorN_expr.hpp"	// apply a functor with 1/2/3 arguments to expressions
#include "cmath_expr.hpp"	// apply a cmath function to expressions
#include "vec_functors.hpp"	// apply a vector function to expressions
