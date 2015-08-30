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



/** \file FileWriterInterface.h */

/** \defgroup Output Output module
 * @{ <img src="output.png"> @}
 */

#ifndef AMDIS_FILEWRITER_INTERFACE_H
#define AMDIS_FILEWRITER_INTERFACE_H

#include <vector>
#include <string>
#include "AMDiS_fwd.h"
#include "Mesh.h"

namespace AMDiS
{

  class FileWriterInterface
  {
  public:
    FileWriterInterface()
      : filename(""),
        appendIndex(0),
        indexLength(5),
        indexDecimals(3),
        createSubDir(-1),
        tsModulo(1),
        timeModulo(-1.0),
        lastWriteTime(-1.0),
        traverseLevel(-1),
        traverseFlag(Mesh::CALL_LEAF_EL),
        writeElement(NULL)
    {}

    virtual ~FileWriterInterface() {}

    /** \brief
     * Interface. Must be overridden in subclasses.
     * \param time time index of solution std::vector.
     * \param force enforces the output operation for the last timestep.
     */
    virtual void writeFiles(AdaptInfo* adaptInfo, bool force,
                            int level = -1,
                            Flag traverseFlag = Mesh::CALL_LEAF_EL,
                            bool (*writeElem)(ElInfo*) = NULL) = 0;

    /// Test whether timestep should be written
    virtual bool doWriteTimestep(AdaptInfo* adaptInfo, bool force);

    std::string getFilename()
    {
      return filename;
    }

    void setFilename(std::string n)
    {
      filename = n;
    }

    void setWriteModulo(int tsModulo_ = 1, double timeModulo_ = -1.0)
    {
      tsModulo = tsModulo_;
      timeModulo = timeModulo_;
    }

    void setTraverseProperties(int level,
                               Flag flag,
                               bool (*writeElem)(ElInfo*))
    {
      traverseLevel = level;
      traverseFlag |= flag;
      writeElement = writeElem;
    }

  protected:
    /// Reads all file writer dependend parameters from the init file.
    virtual void readParameters(std::string name);

    /// create a filename that includes the timestep and possibly a processor ID in parallel mode
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    void getFilename(AdaptInfo* adaptInfo, std::string& fn, std::string& paraFilename, std::string& postfix);
#else
    void getFilename(AdaptInfo* adaptInfo, std::string& fn);
#endif

    /// Used filename prefix.
    std::string filename;

    /// 0: Don't append time index to filename prefix.
    /// 1: Append time index to filename prefix.
    int appendIndex;

    /// Total length of appended time index.
    int indexLength;

    /// Number of decimals in time index.
    int indexDecimals;

    /// create a subdirectory where to put the files
    int createSubDir;

    /// Timestep modulo: write only every tsModulo-th timestep!
    int tsModulo;

    /// Time modulo: write at first iteration after lastWriteTime + timeModulo
    double timeModulo;
    double lastWriteTime;

    /// Traverse properties
    int traverseLevel;
    Flag traverseFlag;
    bool (*writeElement)(ElInfo*);
  };

}

#endif
