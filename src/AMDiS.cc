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


#include "AMDiS.h"
#include <string>
#ifdef HAVE_ZOLTAN
#include <zoltan_cpp.h>
#endif
#ifdef MTL_HAS_HYPRE
#include "mpi.h"
#endif
#include "boost/program_options.hpp"

#if defined HAVE_PETSC || defined HAVE_SEQ_PETSC || defined HAVE_PARALLEL_PETSC
#include <petsc.h>
#endif

namespace AMDiS {

  using namespace std;

#if defined(HAVE_PARALLEL_MTL4)
  mtl::par::environment* mtl_environment = NULL;
#endif

  void init(int argc, char **argv, std::string initFileName)
  {
  #ifdef HAVE_PARALLEL_MTL4
    mtl_environment = new mtl::par::environment(argc, argv);
  #elif defined HAVE_PETSC || defined HAVE_SEQ_PETSC || defined HAVE_PARALLEL_PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);
  #elif defined MTL_HAS_HYPRE
    MPI_Init(&argc, &argv);
  #endif
    
  #if defined(HAVE_PARALLEL_DOMAIN_AMDIS)
    Parallel::mpi::startRand();
  #else
    srand(time(0));
  #endif    
    
  #ifdef HAVE_ZOLTAN
    float zoltanVersion = 0.0;
    Zoltan_Initialize(argc, argv, &zoltanVersion);
  #endif
    
    Parameters::clearData();
    
    // read commandline arguments
    namespace po = boost::program_options;
    
    // Declare the supported options.
    po::options_description desc("Usage: " + std::string(argv[0]) + " init-file [options]\nAllowed options");
    desc.add_options()
	("help", "produce help message")
	("init-file", po::value<std::string>(), "set init file")
	("parameters", po::value<std::string>(), "set parameter in init file\nsyntax: \"key1: value1; key2: value2...\"");

    po::options_description hidden("Hidden options");
    hidden.add_options()
	("unknown", po::value<std::vector<std::string> >(), "unknown options");
	
    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);
	
    // first argument is init-filename
    po::positional_options_description p;
    p.add("init-file", 1);
    p.add("unknown", -1);

    // parse comandline
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).allow_unregistered().run(), vm);
    po::notify(vm);    

    // print help message
    if (vm.count("help")) {
      cout << desc << "\n";
      exit(1);
    }
      
    // set parameters before reading the initfile
    if (vm.count("parameters"))
      Parameters::readArgv(vm["parameters"].as<std::string>());
    
    if (initFileName == "") {
      if (vm.count("init-file"))
	Parameters::init(vm["init-file"].as<std::string>());
      else
	throw(std::runtime_error("No init file specified!"));
    } else {
      Parameters::init(initFileName);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parameters::get("parallel->log main rank", Msg::outputMainRank);
#endif

    // reset parameters from command line
    bool ignoreCommandline = false;
    Parameters::get("ignore commandline options", ignoreCommandline);
    if (vm.count("parameters") && !ignoreCommandline)
      Parameters::readArgv(vm["parameters"].as<std::string>(),0);
    
    
    // initialize global strcutures using parameters
    Global::init();
  }


  void init(std::string initFileName)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    FUNCNAME("AMDiS::init()");
    ERROR_EXIT("Does not work in parallel!\n");
#endif

    Parameters::init(initFileName);
  }


  void finalize()
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->exitParallelization();
    delete Parallel::MeshDistributor::globalMeshDistributor;
#endif
    
#ifdef HAVE_PARALLEL_MTL4
    if (mtl_environment)
      delete mtl_environment;
#elif defined HAVE_PETSC || defined HAVE_SEQ_PETSC || defined HAVE_PARALLEL_PETSC
    PetscFinalize();
#elif defined MTL_HAS_HYPRE
    MPI_Finalize();
#endif
    
  }

}
