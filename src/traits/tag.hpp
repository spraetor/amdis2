/** \file tag.hpp */

#pragma once

namespace AMDiS 
{
  namespace tag
  {
    struct scalar {};
    struct vector {};
    struct matrix {};
    struct unknown {};
    struct dummy {};
    
    struct expression {};
  }
   
  namespace traits
  {
    template <class T, class = void>
    struct Tag {};
    
    template <class T>
    struct Tag<T, typename T::tag>
    {
      using type = typename T::tag;
    };
    
    template <class T>
    using Tag_t = typename Tag<T>::type;
    
    template <class T, class Tag_>
    using HasTag = std::is_same<Tag_t<T>, Tag_>;
    
  } // end namespace traits   
    
} // end namespace AMDiS
