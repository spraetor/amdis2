/** \file AdaptStationary.h */

/** \defgroup Adaption Adaption module
 * @{ <img src="adaption.png"> @}
 *
 * \brief
 * Contains all classes needed for adaption.
 */

#pragma once

#include <string>

#include "AMDiS_fwd.h"
#include "AdaptBase.h"

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

    /// Implementation of AdaptBase::adapt()
    virtual int adapt() override;
  };

} // end namespace AMDiS
