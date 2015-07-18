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



/** \file AdaptStationary.h */

/** \defgroup Adaption Adaption module
 * @{ <img src="adaption.png"> @}
 *
 * \brief
 * Contains all classes needed for adaption.
 */

#ifndef AMDIS_ADPATSTATIONARY_H
#define AMDIS_ADPATSTATIONARY_H

#include <string>
#include "Flag.h"
#include "AdaptInfo.h"
#include "AdaptBase.h"
#include "AMDiS_fwd.h"

namespace AMDiS {

  /** \ingroup Adaption 
   * \brief
   * AdaptStationary contains information about the adaptive procedure and the
   * adapt procedure itself
   */
  class AdaptStationary : public AdaptBase
  {
  public:
    /// Creates a AdaptStationary object with given name.
    AdaptStationary(std::string name,
		    ProblemIterationInterface &prob,
		    AdaptInfo &info);

    /// Destructor
    virtual ~AdaptStationary() {}

    /// Implementation of AdaptBase::adapt()
    virtual int adapt() override;
  };

} // end namespace AMDiS

#endif // AMDIS_ADPATSTATIONARY_H
