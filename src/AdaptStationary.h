/** \file AdaptStationary.h */

/** \defgroup Adaption Adaption module
 * @{ <img src="adaption.png"> @}
 *
 * \brief
 * Contains all classes needed for adaption.
 */

#pragma once

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
