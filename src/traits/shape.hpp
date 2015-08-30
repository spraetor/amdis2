/** \file tag.hpp */

#pragma once

namespace AMDiS
{
  namespace shape
  {
    struct scalar {};

    struct vector {};
    struct row_vector : vector {};
    struct col_vector : vector {};

    struct matrix {};
    struct sqr_matrix : matrix {};

    struct unknown {};

  } // end namespace shape

  namespace traits
  {
    template <class T, class = void>
    struct Shape {};

    template <class T>
    struct Shape<T, typename T::shape>
    {
      using type = typename T::shape;
    };

    template <class T>
    using Shape_t = typename Shape<T>::type;

    template <class T, class Shape_>
    using HasShape = std::is_base_of<Shape_, Shape_t<T>>;

    template <class T>
    using BaseShape = if_then_else<
                      HasShape<T, shape::scalar>, shape::scalar, if_then_else<
                      HasShape<T, shape::vector>, shape::vector, if_then_else<
                      HasShape<T, shape::matrix>, shape::matrix, shape::unknown>>>;

  } // end namespace traits
} // end namespace AMDiS
