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



/** \file ValueWriter.h */

#ifndef AMDIS_VALUEWRITER_H
#define AMDIS_VALUEWRITER_H

#include "Global.h"
#include "Flag.h"
#include "Mesh.h"
#include "DataCollector.h"
#include "AMDiS_fwd.h"

namespace AMDiS { namespace io {

  /** \ingroup Output
   * \brief
   * ValueWriter is a namespace which provied writers to write the values of a DOFVector
   * to an ascii file named values->name.'dat'. This output is done
   * via two leaf-traversals of values->feSpace->mesh. In the first traversal
   * the values at the vertices are printed, in the second these at the
   * interpolation points of each element. For a closer disription of the
   * output format see (...link fehlt noch)
   */
  namespace ValueWriter
  {
    /// Writes DOFVector values to a file.
    void writeValues(DataCollector<> *dc, std::string filename,
		     double time = 0.0, int level = -1,
		     Flag traverseFlag = Mesh::CALL_LEAF_EL,
		     bool (*writeElem)(ElInfo*) = NULL);

  } // end namespace ValueWriter
} } // end namespace io, AMDiS

#endif // AMDIS_VALUEWRITER_H
