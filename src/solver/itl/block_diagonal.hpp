#ifndef ITL_PC_BLOCK_DIAGONAL_INCLUDE
#define ITL_PC_BLOCK_DIAGONAL_INCLUDE

#include <boost/numeric/itl/pc/diagonal.hpp>

#include "solver/BlockMTLMatrix.h"

namespace itl
{
  namespace pc
  {

    /// Diagonal Preconditioner
    template <>
    class diagonal<AMDiS::BlockMTLMatrix, mtl::Collection<AMDiS::BlockMTLMatrix>::value_type>
    {
    public:
      typedef mtl::Collection<AMDiS::BlockMTLMatrix>::value_type   Value;
      typedef Value                                                value_type;
      typedef mtl::Collection<AMDiS::BlockMTLMatrix>::size_type    size_type;
      typedef diagonal                                             self;

      /// Constructor takes matrix reference
      explicit diagonal(const AMDiS::BlockMTLMatrix& A) : inv_diag(num_rows(A))
      {
        MTL_THROW_IF(num_rows(A) != num_cols(A), mtl::matrix_not_square());
        using math::reciprocal;

        for (size_t i = 0; i < A.n_rows; i++)
          inv_diag[A.getRowRange(i)] = mtl::matrix::diagonal(A.getSubMatrix(i,i));

        for (size_type i = 0; i < num_rows(A); ++i)
          inv_diag[i]= reciprocal(inv_diag[i]);
      }

      /// Member function solve, better use free function solve
      template <typename Vector>
      Vector solve(const Vector& x) const
      {
        Vector y(resource(x));
        solve(x, y);
        return y;
      }

      template <typename VectorIn, typename VectorOut>
      void solve(const VectorIn& x, VectorOut& y) const
      {
        y.checked_change_resource(x);
        MTL_THROW_IF(size(x) != size(inv_diag), mtl::incompatible_size());
        for (size_type i= 0; i < size(inv_diag); ++i)
          y[i]= inv_diag[i] * x[i];
      }

      /// Member function for solving adjoint problem, better use free function adjoint_solve
      template <typename Vector>
      Vector adjoint_solve(const Vector& x) const
      {
        Vector y(resource(x));
        adjoint_solve(x, y);
        return y;
      }

      template <typename VectorIn, typename VectorOut>
      void adjoint_solve(const VectorIn& x, VectorOut& y) const
      {
        using mtl::conj;
        y.checked_change_resource(x);
        MTL_THROW_IF(size(x) != size(inv_diag), mtl::incompatible_size());
        for (size_type i= 0; i < size(inv_diag); ++i)
          y[i]= conj(inv_diag[i]) * x[i];
      }

    protected:
      mtl::vector::dense_vector<value_type>    inv_diag;
    };


  }
} // namespace itl::pc

#endif // ITL_PC_DIAGONAL_INCLUDE
