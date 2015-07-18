#include <string>

#include "AdaptInfo.h"

namespace AMDiS 
{
  void AdaptInfo::setScalContents(int newSize) 
  {
    int oldSize = static_cast<int>(scalContents.size());

    if (newSize > oldSize) { 
      scalContents.resize(newSize);

      for (int i = oldSize; i < newSize; i++)
	scalContents[i] = 
	  new ScalContent(name + "[" + std::to_string(i) + "]"); 
    }
  }

} // end namespace AMDiS
