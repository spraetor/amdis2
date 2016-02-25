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


/** \file Reader.h */

#ifndef AMDIS_READER_H
#define AMDIS_READER_H

#include <cstring>
#include "DOFVector.h"
#include "SystemVector.h"

#include "ArhReader.h"
#include "Arh2Reader.h"
#include "Arh3Reader.h"
#include "MacroReader.h"
#include "ValueReader.h"
#include "XYZReader.h"

#ifdef HAVE_PNG
#include "PngReader.h"
#endif

#ifdef HAVE_EXTENSIONS
#include "VtkReader.h"
#endif

#include "detail/ReaderWriter.h"


namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     * \brief reader interface that selects a reader depending on the file extension of the filename.
     *
     * @param filename Filename of the input file, used to extract extension.
     * @param container One of DOFVector&, std::vector<DOFVector*>,
     * 			SystemVector&, Mesh&.
     **/
    template<typename Container>
    void readFile(std::string filename,
                  Container& container)
    {
      std::string ext = filename.substr(filename.find_last_of("."));

      Mesh* mesh = detail::getMesh(container);
      if (ext == ".1d" || ext == ".2d" || ext == ".3d")
      {
        std::string periodicFilename = "";
        int check = 1;
        Parameters::get(mesh->getName() + "->periodic file", periodicFilename);
        Parameters::get(mesh->getName() + "->check", check);
        MacroReader::readMacro(filename, mesh, periodicFilename, check);
      }
      else if (ext == ".arh")
      {
        if (Arh2Reader::isReadable(filename))
          Arh2Reader::readFile(filename, container);
        else if (Arh3Reader::isReadable(filename))
          Arh3Reader::readFile(filename, container);
        else
          ArhReader::readFile(filename, container);
      }
      else if (ext == ".dat")
      {
        bool preserveMacroFileInfo = false;
        Parameters::get(mesh->getName() + "->preserve macroFileInfo",
                        preserveMacroFileInfo);
        TEST_EXIT(mesh->getMacroFileInfo() != NULL && preserveMacroFileInfo)
        ("Yout have to set the flag 'mesh_name->preserve macroFileInfo' "
         "to 'true' (%d), in order to read .dat-files",
         int(preserveMacroFileInfo));

        ValueReader::readValue(filename, mesh, container,
                               mesh->getMacroFileInfo());
      }
#ifdef HAVE_PNG
      else if (ext == ".png")
      {
        PngReader::readFile(filename, container);
      }
#endif

#ifdef HAVE_EXTENSIONS
      else if (ext == ".vtu")
      {
        VtkReader::readFile(filename, container);
      }
#endif
      else
      {
        ERROR_EXIT("File-extensions %s can not be assigned to a reader!\n", ext.c_str());
      }
    }


    /** \ingroup Output
     * \brief Wrapper for pointers to container types
     **/
    template<typename Container>
    void readFile(std::string filename,
                  Container* container)
    {
      readFile(filename, *container);
    }


    /** \ingroup Output
     * \brief specialization for read meshes
     **/
    inline
    void readFile(std::string filename,
                  Mesh* mesh)
    {

      std::string ext = filename.substr(filename.find_last_of("."));
      if (ext == ".1d" || ext == ".2d" || ext == ".3d")
      {
        std::string periodicFilename = "";
        int check = 1;
        Parameters::get(mesh->getName() + "->periodic file", periodicFilename);
        Parameters::get(mesh->getName() + "->check", check);
        MacroReader::readMacro(filename, mesh, periodicFilename, check);
      }
      else if (ext == ".arh")
      {
        if (Arh2Reader::isReadable(filename))
          Arh2Reader::readFile(filename, mesh);
        else if (Arh3Reader::isReadable(filename))
          Arh3Reader::readFile(filename, mesh);
        else
          ArhReader::readFile(filename, mesh);
      }
      else
      {
        ERROR_EXIT("File-extensions %s can not be assigned to a reader!\n", ext.c_str());
      }
    }
  } // end namespace io
} // end namespace AMDiS

#endif // AMDIS_READER_H
