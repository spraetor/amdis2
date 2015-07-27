/** \file Global.h */

/** \mainpage AMDiS
 * @{ <img src="vis.png"> @}
 */

/** \defgroup Common Common
 */

#pragma once

#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <functional>
#include <float.h>
#include <time.h>

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include "AMDiS_fwd.h"
#include "Config.h"
#include "GeoIndex.h"
#include "Log.h"
#include "traits/meta_basic.hpp"

namespace AMDiS 
{
  /** \ingroup Common
   * \brief
   * Static class wich holds global information about the world and all types of
   * elements.
   */
  class Global
  {
  public:
    /// returns a pointer to \ref referenceElement [dim]. With this pointer you
    /// can get information about the element via Element's getGeo method.
    static const Element *getReferenceElement(int dim) 
    {
      FUNCNAME("Global::getReferenceElement()");
      TEST_EXIT(dim > 0 && dim < 4)("invalid dim: %d\n", dim);
      return referenceElement[dim];
    }

    /// returns geometrical information. Currently this is only dimOfWorld.
    static inline int getGeo(GeoIndex p) 
    {
      if (WORLD == p) 
        return dimOfWorld; 

      ERROR_EXIT("Illegal request for geometry data: part = %d!\n", p);
      return 0;
    }

    /// Returns geometrical information about elements of the dimension dim.
    /// getGeo(VERTEX, 3) returns 4 because a Tetrahedron has 4 vertices.
    static inline int getGeo(GeoIndex p, int dim) 
    {
      TEST_EXIT_DBG(p >= MINPART && p <= MAXPART)
        ("Calling for invalid geometry value %d\n",p);
      TEST_EXIT_DBG(dim >= 0 && dim < 4)
        ("invalid dim: %d\n", dim);
      TEST_EXIT_DBG((dim != 0) || (p == PARTS || p == VERTEX || p == EDGE || p == FACE))
        ("dim = 0\n");

      return geoIndexTable[dim][p - MINPART];
    }

    /// Inits the Global class with the help of Parameters.
    static void init();

    static void clear();

  private:
    /// Global is a pure static class. So the constructor is private to avoid
    /// instantiation.
    Global();

  private:
    /// Dimension of the simulated world
    static int dimOfWorld;

    /// contains a pointer to a Line, a Triangle, and a Tetrahedron.
    /// This allows the access to information of the concrete elements via
    /// the dimension index.
    /// referenceElement[3]->getGeo(VERTEX) gives the number of vertices of a
    /// Tetrahedron wich is 4 => no switch statement is necessary.
    static Element *referenceElement[4];

    /// Stores the precalculated results that should be returned by Global::getGeo.
    static std::vector<std::vector<int> > geoIndexTable;
  };
  
  
  /// Size-policy used in definition of WorldVector and FixVec
  template <GeoIndex G>
  struct FixedSize
  {
    /// return argument \param s
    static constexpr size_t eval(size_t s)
    {
      return aux(_geo<G>(), s);
    }

  private:
    /// calculate size for DimVector
    template <GeoIndex H>
    static constexpr size_t aux(_geo<H>, size_t dim)
    {
      return Global::getGeo(H, dim);
    }
    
    /// calculate size for WorldVector
    static constexpr size_t aux(_geo<WORLD>, size_t = 0)
    {
      return Global::getGeo(WORLD);
    }
  };
  
  /// maximal size to allocate for container types, based on GeoIndex
  template <GeoIndex> struct MaxSize : int_<-1> {};
  
  /// \cond HIDDEN_SYMBOLS
  template <> struct MaxSize<CENTER> : int_< 1> {};
  
#ifdef DOW
  template <> struct MaxSize<WORLD>  : int_<DOW> {};
#else
  template <> struct MaxSize<WORLD>  : int_< 3> {}; // upper bound
#endif
    
#ifdef DIM
  template <> struct MaxSize<DIMEN>  : int_<DIM> {};
  template <> struct MaxSize<VERTEX> : int_<DIM+1> {};
  template <> struct MaxSize<PARTS>  : int_<DIM+1> {};
  template <> struct MaxSize<NEIGH>  : int_<DIM+1> {};
  template <> struct MaxSize<EDGE>   : int_<(DIM==1?1:(DIM==2?3:6))> {};
  template <> struct MaxSize<FACE>   : int_<(DIM==1?0:(DIM==2?1:4))> {};
#else
  // upper bounds
  template <> struct MaxSize<DIMEN>  : int_< 3> {};
  template <> struct MaxSize<VERTEX> : int_< 4> {};
  template <> struct MaxSize<PARTS>  : int_< 4> {};
  template <> struct MaxSize<NEIGH>  : int_< 4> {};
  template <> struct MaxSize<EDGE>   : int_< 6> {};
  template <> struct MaxSize<FACE>   : int_< 4> {};
#endif
  /// \endcond

} // end namespace AMDiS

