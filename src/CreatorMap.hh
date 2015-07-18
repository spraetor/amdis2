#include <map>

namespace AMDiS {

  template <class BaseClass>
  void CreatorMap<BaseClass>::clear() 
  {
    typename std::map< std::string, CreatorInterface<BaseClass>* >::iterator it;

    for (it = creatorMap.begin(); it != creatorMap.end(); ++it)
      delete it->second;
  }
  
} // end namespace AMDiS

