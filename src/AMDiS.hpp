#pragma once

// std c++ headers
#include <string>

namespace AMDiS
{
  void init(int argc, char** argv, std::string initFileName = "");

  void init(std::string initFileName);

  void finalize();

} // end namespace AMDiS
