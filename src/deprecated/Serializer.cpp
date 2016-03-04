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


#include "Serializer.h"

namespace AMDiS
{

  namespace SerUtil
  {

    void serialize(std::ostream& out, DofEdge& data)
    {
      out.write(reinterpret_cast<const char*>(&(data.first)), sizeof(int));
      out.write(reinterpret_cast<const char*>(&(data.second)), sizeof(int));
    }

    void deserialize(std::istream& in, DofEdge& data)
    {
      in.read(reinterpret_cast<char*>(&(data.first)), sizeof(int));
      in.read(reinterpret_cast<char*>(&(data.second)), sizeof(int));
    }

    void serialize(std::ostream& out, DofFace& data)
    {
      out.write(reinterpret_cast<const char*>(&(data.get<0>())), sizeof(int));
      out.write(reinterpret_cast<const char*>(&(data.get<1>())), sizeof(int));
      out.write(reinterpret_cast<const char*>(&(data.get<2>())), sizeof(int));
    }

    void deserialize(std::istream& in, DofFace& data)
    {
      in.read(reinterpret_cast<char*>(&(data.get<0>())), sizeof(int));
      in.read(reinterpret_cast<char*>(&(data.get<1>())), sizeof(int));
      in.read(reinterpret_cast<char*>(&(data.get<2>())), sizeof(int));
    }

  }

}
