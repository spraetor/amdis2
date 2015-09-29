/** \file Traits.h */

#pragma once

// AMDiS headers
#include <traits/basic.hpp>
#include <traits/shape.hpp>


namespace AMDiS
{
  namespace traits
  {
    // multiplication shapes
    // _________________________________________________________________________
    namespace detail
    {
      template <class Shape1, class Shape2>
      struct mult_shape
      {
        typedef no_valid_type type;
      };

      /// Scalar*Scalar => Scalar
      template <>
      struct mult_shape<shape::scalar, shape::scalar> : Id<shape::scalar> {};

      /// Vec*Vec => Scalar (dot-product) - if no orientation of the vector is given
      template <>
      struct mult_shape<shape::vector, shape::vector> : Id<shape::scalar> {};

      // inner product
      template <>
      struct mult_shape<shape::row_vector, shape::col_vector> : Id<shape::scalar> {};

      // outer product
      template <>
      struct mult_shape<shape::col_vector, shape::row_vector> : Id<shape::matrix> {};

      /// Mat*Mat => Mat
      template <>
      struct mult_shape<shape::matrix, shape::matrix> : Id<shape::matrix> {};

      /// Vec*Scalar => Vector
      template <>
      struct mult_shape<shape::vector, shape::scalar> : Id<shape::vector> {};

      /// Scalar*Vector => Vector
      template <>
      struct mult_shape<shape::scalar, shape::vector> : Id<shape::vector> {};

      /// Matrix*Scalar => Matrix
      template <>
      struct mult_shape<shape::matrix, shape::scalar> : Id<shape::matrix> {};

      /// Scalar*Matrix => Matrix
      template <>
      struct mult_shape<shape::scalar, shape::matrix> : Id<shape::matrix> {};

      /// Matrix*Vector => Vector
      template <>
      struct mult_shape<shape::matrix, shape::vector> : Id<shape::vector> {};

    } // end namespace detail

    /// determines the type of the product T1*T2
    template <class T1, class T2>
    struct MultShape
    {
    private:
      class S1 = typename detail::mult_shape<Shape_t<T1>,   Shape_t<T2>>::type;
      class S2 = typename detail::mult_shape<BaseShape<T1>, BaseShape<T2>>::type;
    public:
      using type = if_then_else<std::is_same<S1, no_valid_type>, S2, S1>;
    };

    template <class T1, class T2>
    using MultShape_t = typename MultShape<T1, T2>::type;


    // addition shapes
    // _________________________________________________________________________
    namespace detail
    {
      template <class Shape1, class Shape2>
      struct add_shape
      {
        typedef no_valid_type type;
      };

      /// add types of the same shape
      template <class Shape>
      struct add_shape<Shape, Shape> : Id<Shape> {};

      /// shape + scalar => shape
      template <class Shape>
      struct add_shape<Shape, shape::scalar> : Id<Shape> {};

      /// scalar + shape => shape
      template <class Shape>
      struct add_shape<shape::scalar, Shape> : Id<Shape> {};

    } // end namespace detail

    /// determines the type of the product T1*T2
    template <class T1, class T2>
    using AddShape = detail::add_shape<BaseShape<T1>, BaseShape<T2>>;

    template <class T1, class T2>
    using AddShape_t = typename AddShape<T1, T2>::type;

  } // end namespace traits

} // end namespace AMDiS
