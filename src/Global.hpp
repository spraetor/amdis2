/** \file Global.h */

/** \mainpage AMDiS
 * @{ <img src="vis.png"> @}
 */

/** \defgroup Common Common
 */

#pragma once

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include "AMDiS_fwd.hpp"
#include "Config.hpp"
#include "GeoIndex.hpp"
#include "Log.hpp"
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
    static Element const* getReferenceElement(int dim)
    {
      FUNCNAME("Global::getReferenceElement()");
      TEST_EXIT(dim > 0 && dim < 4)("invalid dim: %d\n", dim);
      return referenceElement[dim];
    }

    /// returns geometrical information. Currently this is only dimOfWorld.
#if AMDIS_FIXED_SIZE && defined(DOW)
    static constexpr int getGeo(GeoIndex p)
    {
      return DOW;
      //       return p == WORLD ? DOW : throw std::runtime_error("Illegal request for geometry data!");
    }
#else
    static int getGeo(GeoIndex p)
    {
      if (WORLD == p)
        return dimOfWorld;

      ERROR_EXIT("Illegal request for geometry data: part = %d!\n", p);
      return 0;
    }
#endif

    /// Returns geometrical information about elements of the dimension dim.
    /// getGeo(VERTEX, 3) returns 4 because a Tetrahedron has 4 vertices.
    static int getGeo(GeoIndex p, int dim)
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
#if AMDIS_FIXED_SIZE && defined(DOW)
    static constexpr int dimOfWorld = DOW;
#else
    static int dimOfWorld;
#endif

    /// contains a pointer to a Line, a Triangle, and a Tetrahedron.
    /// This allows the access to information of the concrete elements via
    /// the dimension index.
    /// referenceElement[3]->getGeo(VERTEX) gives the number of vertices of a
    /// Tetrahedron wich is 4 => no switch statement is necessary.
    static Element* referenceElement[4];

    /// Stores the precalculated results that should be returned by Global::getGeo.
    static std::vector<std::vector<int>> geoIndexTable;
  };



  /// Size-policy used in definition of WorldVector and FixVec
  template <GeoIndex G>
  struct FixedSize
  {
    /// return argument \param s
    template <class size_type>
    static size_type eval(size_type s)
    {
      return evalImpl(_geo<G>(), s);
    }

  private:
    /// calculate size for DimVector
    template <GeoIndex H, class size_type>
    static size_type evalImpl(_geo<H>, size_type dim)
    {
      return size_type(Global::getGeo(H, dim));
    }

    /// calculate size for WorldVector
    template <class size_type>
#if AMDIS_FIXED_SIZE && defined(DOW)
    static constexpr size_t evalImpl(_geo<WORLD>, size_type /* dim */ = 0)
#else
    static size_type evalImpl(_geo<WORLD>, size_type /* dim */ = 0)
    {
      return size_type(Global::getGeo(WORLD));
    }
#endif
  };

} // end namespace AMDiS

