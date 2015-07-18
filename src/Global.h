/** \file Global.h */

/** \mainpage AMDiS
 * @{ <img src="vis.png"> @}
 */

/** \defgroup Common Common
 */

#pragma once

#include "Config.h"

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

#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "boost/tuple/tuple.hpp"
#include "AMDiS_fwd.h"
#include "OpenMP.h"

namespace AMDiS {

  extern const char *funcName;

  /// Speciefies the norm used by Estimator.
  typedef enum { NO_NORM = 0, 
		 H1_NORM = 1, 
		 L2_NORM = 2 } Norm;
  
  // may be used in the future
  // TODO: No good idea to use DofIndex class as index type, better: use typedef
  struct DofIndex 
  {
    typedef signed int size_type;
  };

  /// Datatype for degrees of freedom 
//   typedef signed int DegreeOfFreedom;
  typedef DofIndex::size_type DegreeOfFreedom;

  /// Defines type for a vector of DOF pointers.
  typedef std::vector<const DegreeOfFreedom*> DofContainer;

  typedef std::set<const DegreeOfFreedom*> DofContainerSet;

  /// Defines a type for global edge identification via its DOFs.
  typedef std::pair<DegreeOfFreedom, DegreeOfFreedom> DofEdge;

  /// Defines a tzpe for global face identiication via its DOFs.
  typedef boost::tuple<DegreeOfFreedom, DegreeOfFreedom, DegreeOfFreedom> DofFace;


  /// Returns the GeoIndex of d for dimension dim.
#define INDEX_OF_DIM(d, dim) (static_cast<GeoIndex>((d == dim) ? CENTER : d + 1))

  /// Returns the dimension of GeoIndex ind for dimension dim
#define DIM_OF_INDEX(ind, dim) ((static_cast<int>(ind) == 0) ? dim : static_cast<int>(ind) - 1)


#if SUPPRESS_OUTPUT
#define PRINT_LINE(stream, line)
#else
#if HAVE_PARALLEL_DOMAIN_AMDIS
#define PRINT_LINE(stream, line) stream << "[" << MPI::COMM_WORLD.Get_rank() << "] " << line
#else
#define PRINT_LINE(stream, line) stream << line
#endif
#endif

  void waitSec(int seconds);

  void processMemUsage(double& vm_usage, double& resident_set, bool inMegaByte = true);

  /// Content comparision of two pointers. Used e.g. for find_if
  template<typename T>
  struct comparePtrContents : public std::binary_function<T*, T*, bool>
  {
    /// Overrides binary_function::operator()
    bool operator()(T* a, T* b) const 
    {
      return (*a == *b);
    }
  };

  
  /// check for file existence
  inline bool file_exists(const std::string filename)
  {
#ifdef _MSC_VER
    return _access(filename.c_str(), 0) == 0;
#else
    return access(filename.c_str(), F_OK) == 0;
#endif
  }

  
  /// trim std::string
  inline std::string trim(const std::string& oldStr)
  {
    std::string swap(oldStr);
    boost::algorithm::trim(swap);
    return swap;
  }


  // TODO: move to Logging-class
  /** \brief
   * Manages the output of messages, warnings, errors, ...
   * Used by the macros FUNCNAME, ERROR, ERROR_EXIT, WARNING, TEST, MSG, INFO,
   * PRINT_INFO, WAIT, WAIT_REALLY.
   * Don't use this class directly but only via these macros!
   */
  class Msg
  {
  public:
    /// Prints a formated message to the message stream
    static void print(const char *format, ...);

    /// Prints a formated message to the error stream
    static void print_error(const char *format, ...);

    /// Prints a formated message to the error stream and exits 
    static void print_error_exit(const char *format, ...);

    ///
    static void catch_error_exit(const char *format, ...) {}

    /// Prints an error message with funcname, file, and line to the error stream
    static void print_error_funcname(const char *funcname,
				     const char *file, 
				     int line);

    /// Prints a warning to the message stream
    static void print_warn(const char *format, ...);

    /// Prints a warning with funcname, file, and line to the message stream
    static void print_warn_funcname(const char *funcname,
				    const char *file, 
				    int line);

    /// Prints the funcname to the message stream
    static void print_funcname(const char *funcname);

    /// Changes the message stream
    static void change_out(std::ostream*);

    /// Changes the error stream 
    static void change_error_out(std::ofstream *fp);

    /// Creates a filestream and sets the error stream to this filestream
    static void open_error_file(const char *filename, OPENMODE);

    /// Sets \ref msgInfo
    static void setMsgInfo(int info) 
    { 
      msgInfo = info; 
    }

    /// Returns \ref msgInfo
    static int  getMsgInfo() 
    { 
      return msgInfo; 
    }

    /// Sets \ref msgWait
    static void setMsgWait(bool wait) 
    { 
      msgWait = wait; 
    }

    /// Returns \ref msgWait
    static bool getMsgWait() 
    { 
      return msgWait; 
    }

    /// Waits for enter if w is true
    static void wait(bool w);

    /// Returns \ref out
    static std::ostream *getOutStream() 
    { 
      return out; 
    }

    /// Returns \ref error
    static std::ostream *getErrorStream() 
    { 
      return error; 
    }

  public:
#if HAVE_PARALLEL_DOMAIN_AMDIS
    /// In parallel computations, when this variable is true, only the 0 rank will
    /// print messages to the output stream. Error messages and warnings are always
    /// printed from all ranks.
    static bool outputMainRank;
#endif

  protected:
    /// Message stram
    static std::ostream *out;

    /// Error stream
    static std::ostream *error;

    /// Remember funcName to avoid multiple output of funcName within the same
    /// function call
    static ThreadPrivate<const char*> oldFuncName;

    /// Global info level
    static int msgInfo;

    /// Spezifies whether to wait when WAIT is called
    static bool msgWait;
  };

  // ===========================================================================
  // ===== message macros ======================================================
  // ===========================================================================

  /// Should be the first call in every functions. It defines the current 
  /// function name nn for message output via MSG, WARNING, ...
#define FUNCNAME(nn) const char *funcName; funcName = nn;
  
#ifdef NDEBUG
  #define FUNCNAME_DBG(nn)
#else
  #define FUNCNAME_DBG(nn) const char *funcName; funcName = nn;
#endif

  /// prints an error message 
#define ERROR Msg::print_error_funcname(funcName,__FILE__, __LINE__),	\
    Msg::print_error

  /// prints an error message and exits 
#define ERROR_EXIT Msg::print_error_funcname(funcName,__FILE__, __LINE__), \
    Msg::print_error_exit

  /// prints a warning
#define WARNING Msg::print_warn_funcname(funcName,__FILE__, __LINE__),	\
    Msg::print_warn

  /// if test is false, an error message is printed
#define TEST(test) if ((test));else ERROR

  /// if test is false, an error message is printed and the program exits
#define TEST_EXIT(test) if ((test));else ERROR_EXIT

  /// In debug mode, it corresponds to ERROR_EXIT, otherwise it is noop.
#ifdef NDEBUG
  #define TEST_EXIT_DBG(test) if (false) Msg::catch_error_exit
  #define DBG_VAR(var)
#else
  #define TEST_EXIT_DBG(test) if ((test));else ERROR_EXIT
  #define DBG_VAR(var) var
#endif

  /// prints a message
#define MSG Msg::print_funcname(funcName), Msg::print

#ifdef NDEBUG
  #define MSG_DBG
#else
  #define MSG_DBG Msg::print_funcname(funcName), Msg::print
#endif

  /// prints a message, if min(Msg::msgInfo, info) >= noinfo
#define INFO(info,noinfo)						\
  if (Msg::getMsgInfo() && (std::min(Msg::getMsgInfo(), (info)) >= (noinfo))) MSG

  /// prints a message, if min(Msg::msgInfo, info) >= noinfo
#define PRINT_INFO(info,noinfo)						\
  if (Msg::getMsgInfo() && (std::min(Msg::getMsgInfo(), (info)) >= (noinfo))) Msg::print


  /// If the value of Msg::wait is not zero the macro will produce the message 
  /// 'wait for <enter> ...' and will continue after pressing the return or enter
  /// key. Otherwise the program continues without a message.
#define WAIT Msg::wait(Msg::getMsgWait())

  /// produces the message 'wait for <enter> ...' and will continue after
  /// pressing the return or enter key.
#define WAIT_REALLY Msg::wait(true)

  /// internal used indices to represent the different geometrical objects.
  /// Used as parameter for getGeo() and as template parameter for FixVec. 
  typedef enum
    {
      CENTER   = 0, /**< in 1d the center is at the edge, in 2d at the face, in 3d 
		     * at the interior of an element. So a Line has no edge but
		     * only a center, a Triangle has no face but only a center.
		     */
      VERTEX   = 1, /**< index for element vertices.
		     * number of vertices is equal to number of parts and 
		     * neighbours.
		     */
      EDGE     = 2, /**< index for element edges */
      FACE     = 3, /**< index for element faces */
      DIMEN    =-1, /**< index for problem dimension */
      PARTS    =-2, /**< index for parts of an element (vertices in 1d, edges in 2d
		     * , faces in 3d). Number of element parts is equal to number
		     * of vertices and neighbours. 
		     */
      NEIGH    =-3, /**< index for neighbours of an element.
		     * Number of element neighbours is equal to number of 
		     * vertices and parts.
		     */
      WORLD    =-4, /**< index for world dimension */
      BOUNDARY =-5, /**< index for boundary nodes of an element. This could be
		     * vertices, edges or faces.
		     */
      PROJECTION=-6, /**< index for element and boundary projections */
      
      NO_INDEX =-127
    } GeoIndex;

#define MAXPART FACE
#define MINPART PROJECTION


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

#define COMMA ,

  /**
   * \ingroup Assembler
   * \brief
   * Specifies the type of a FirstOrderTerm 
   */
  enum FirstOrderType {
    GRD_PSI,
    GRD_PHI
  };

} // end namespace AMDiS
