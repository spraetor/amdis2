#include "AMDiS.h"

using AMDiS::Msg;

int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  MSG("Hello World!\n");

  AMDiS::finalize();
}


