#define AMDIS_NO_EXTERN_INITFILE
#include "Initfile.hpp"
#undef AMDIS_NO_EXTERN_INITFILE

#include <string>
#ifdef _MSC_VER
#include <io.h>			// _access
#else
#include <unistd.h>
#endif

#include <boost/property_tree/json_parser.hpp>

// a parser for arithmetic expressions
#include <muParser.h>

#include "AMDiS_base.hpp"
#include "Log.hpp"

namespace AMDiS
{
  /// check for file existence
  inline bool file_exists(std::string const& filename)
  {
#ifdef _MSC_VER
    return _access(filename.c_str(), 0) == 0;
#else
    return access(filename.c_str(), F_OK) == 0;
#endif
  }


  namespace detail
  {
    double mu_parser_eval(std::string const& valStr)
    {
      mu::Parser parser;
      parser.DefineConst(_T("M_PI"), m_pi);
      parser.DefineConst(_T("M_E"), m_e);

      parser.SetExpr(valStr);
      return parser.Eval();
    }
  }

  Initfile* Initfile::singlett = NULL;

  /// initialize singleton object an global parameters
  void Initfile::init(std::string in)
  {
    initIntern();
    singlett->read(in);
    singlett->getInternalParameters();
  }


  /// Fill an initfile from a file with filename fn
  void Initfile::read(std::string fn, bool /*force*/)
  {
    TEST_EXIT( file_exists(fn) )
    ("init-file '%s' cannot be opened for reading", fn.c_str());

    using namespace boost::property_tree;
    json_parser::read_json(fn, pt);
  }


  /// read standard values for output and information of parameter-values
  void Initfile::getInternalParameters()
  {
    int val = 0;
    get("level of information", val, 0);
    msgInfo = val;

    val = 1;
    get("WAIT", val, 0);
    msgWait = val;

    val = 1;
    get("parameter information", val, 0);
    paramInfo = val;

    val = 0;
    get("break on missing tag", val, 0);
    breakOnMissingTag = val;

    if (msgInfo == 0)
      paramInfo = 0;
  }


  /// print all parameters to std::cout
  void Initfile::printParameters()
  {
    initIntern();
    // TODO: implement printing of all parameters
  }

  // explicit template instatiation
  template void Initfile::get(std::string, int&, int);
  template void Initfile::get(std::string, double&, int);
  template void Initfile::get(std::string, std::string&, int);


} // end namespace AMDiS
