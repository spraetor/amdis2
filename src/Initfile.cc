#define AMDIS_NO_EXTERN_INITFILE
#include "Initfile.h"
#undef AMDIS_NO_EXTERN_INITFILE

#include "AMDiS_base.h"
#include "Log.h"

#include <string>

namespace AMDiS
{
  Initfile* Initfile::singlett = NULL;

  /// initialize singleton object an global parameters
  void Initfile::init(std::string in)
  {
    initIntern();
    singlett->read(in);
    singlett->getInternalParameters();
  }


  /// Fill an initfile from a file with filename fn
  void Initfile::read(std::string fn, bool force)
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
  template void Initfile::get(const std::string, int&, int);
  template void Initfile::get(const std::string, double&, int);
  template void Initfile::get(const std::string, std::string&, int);


} // end namespace AMDiS
