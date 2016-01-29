/** \file ITL_Solver.h */

#pragma once

// AMDiS includes
#include <solver/LinearSolver.h>
#include <solver/ITL_Runner.h>
#include <MTL4Types.h>

// MTL4 includes
#include <boost/numeric/itl/krylov/bicg.hpp>
#include <boost/numeric/itl/krylov/bicgstab_2.hpp>
#include <boost/numeric/itl/krylov/bicgstab_ell.hpp>
#include <boost/numeric/itl/krylov/bicgstab.hpp>
#include <boost/numeric/itl/krylov/cg.hpp>
#include <boost/numeric/itl/krylov/cgs.hpp>
#include <boost/numeric/itl/krylov/gmres.hpp>
#include <boost/numeric/itl/krylov/idr_s.hpp>
#include <boost/numeric/itl/krylov/qmr.hpp>
#include <boost/numeric/itl/krylov/tfqmr.hpp>

// more solvers defined in AMDiS
#include <solver/itl/minres.hpp>
#include <solver/itl/gcr.hpp>
#include <solver/itl/fgmres.hpp>
#include <solver/itl/fgmres_householder.hpp>
#include <solver/itl/gmres2.hpp>
#include <solver/itl/gmres_householder.hpp>
#include <solver/itl/preonly.hpp>


namespace AMDiS
{
  /**
   * \ingroup Solver
   *
   * \brief
   * Wrapper for MTL4 itl-solvers.
   *
   * One of the following solvers can be chosen:
   * - @ref CGSolver "cg" (conjugate gradient method)
   * - @ref CGSSolver "cgs" (squared conjugate gradient method)
   * - @ref BiCGSolver "bicg" (biconjugate gradient method)
   * - @ref BiCGStabSolver "bicgstab" (stabilized BiCG method)
   * - @ref BiCGStab2Solver "bicgstab2" (stabilized BiCG(l) method with l=2)
   * - @ref QMRSolver "qmr" (Quasi-Minimal Residual method)
   * - @ref TFQMRSolver "tfqmr" (Transposed-Free Quasi-Minimal Residual method)
   * - @ref BiCGStabEllSolver "bicgstab_ell" (stabilized BiCG(l) method)
   * - @ref GMResSolver "gmres" (generalized minimal residual method)
   * - @ref IDRsSolver "idr_s" (Induced Dimension Reduction method)
   * - @ref MinResSolver "minres" (minimal residual method)
   * - @ref GcrSolver "gcr" (generalized conjugate residual method)
   * - @ref FGMResSolver "fgmres" (flexible GMRes method)
   * - @ref PreOnly "preonly" (solver that implements pure preconditioning applied to the rhs)
   */
  template <class SolverType>
  using ITL_Solver = LinearSolver<
    MTLTypes::MTLMatrix, 
    MTLTypes::MTLVector, 
    ITL_Runner<SolverType, MTLTypes::MTLMatrix, MTLTypes::MTLVector>>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::CGSolver
   * \brief ITL_Solver <\ref cg_solver_type> implementation of conjugate gradient 
   * method \implements ITL_Solver
   *
   * Solves a linear system \f$ Ax=b \f$ by the conjugate gradient method (CG) 
   * and can be used for symmetric positive definite system matrices.
   * Right preconditioner is ignored.
   */

  class cg_solver_type
  {
  public:
    cg_solver_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::cg(A, x, b, l, r, iter);
    }
  };
  using CGSolver = ITL_Solver<cg_solver_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::CGSSolver
   * \brief ITL_Solver <\ref cgs_solver_type> implementation of squared conjugate 
   * gradient method \implements ITL_Solver
   *
   * Solves a linear system \f$ Ax=b \f$ by the squared conjugate gradient method 
   * (CGS). Right preconditioner is ignored.
   */

  class cgs_solver_type
  {
  public:
    cgs_solver_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const&, I& iter)
    {
      return itl::cgs(A, x, b, l, iter);
    }
  };
  using CGSSolver = ITL_Solver<cgs_solver_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::BiCGSolver
   * \brief ITL_Solver <\ref bicg_solver_type> implementation of bi-conjugate 
   * gradient method \implements ITL_Solver
   *
   * Solves a linear system \f$ Ax=b \f$ by a BiCG method and can be used for
   * system matrices.
   */

  class bicg_solver_type
  {
  public:
    bicg_solver_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const&, I& iter)
    {
      return itl::bicg(A, x, b, l, iter);
    }
  };
  using BiCGSolver = ITL_Solver<bicg_solver_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::BiCGStabSolver
   * \brief ITL_Solver <\ref bicgstab_type> implementation of stabilized 
   * bi-conjugate gradient method \implements ITL_Solver
   *
   * Solves a linear system \f$ Ax=b \f$ by a stabilized BiCG method and can be 
   * used for system matrices.
   */

  class bicgstab_type
  {
  public:
    bicgstab_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const&, I& iter)
    {
      return itl::bicgstab(A, x, b, l, iter);
    }
  };
  using BiCGStabSolver = ITL_Solver<bicgstab_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::BiCGStab2Solver
   * \brief ITL_Solver <\ref bicgstab2_type> implementation of BiCGStab(l) method 
   * with l=2 \implements ITL_Solver
   *
   * Solves a linear system \f$ Ax=b \f$ by a stabilized BiCG(2) method and can 
   * be used for system matrices.
   */

  class bicgstab2_type
  {
  public:
    bicgstab2_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const&, I& iter)
    {
      return itl::bicgstab_2(A, x, b, l, iter);
    }
  };
  using BiCGStab2Solver = ITL_Solver<bicgstab2_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::QMRSolver
   * \brief ITL_Solver <\ref qmr_solver_type> implementation of Quasi-Minimal 
   * Residual method \implements ITL_Solver
   *
   * Solves a linear system \f$ Ax=b \f$ by the Quasi-Minimal Residual method (QMR).
   */

  class qmr_solver_type
  {
  public:
    qmr_solver_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::qmr(A, x, b, l, r, iter);
    }
  };
  using QMRSolver = ITL_Solver<qmr_solver_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::TFQMRSolver
   * \brief ITL_Solver <\ref tfqmr_solver_type> implementation of Transposed-Free 
   * Quasi-Minimal Residual method \implements ITL_Solver
   *
   * Solves a linear system by the Transposed-Free Quasi-Minimal Residual method 
   * (TFQMR). Does not use preconditioning currently.
   */

  class tfqmr_solver_type
  {
  public:
    tfqmr_solver_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::tfqmr(A, x, b, l, r, iter);
    }
  };
  using TFQMRSolver = ITL_Solver<tfqmr_solver_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::BiCGStabEllSolver
   * \brief ITL_Solver <\ref bicgstab_ell_type> implementation of stabilized 
   * BiCG(ell) method \implements ITL_Solver
   *
   * Solves a linear system by a stabilized BiCG(ell) method and can be used for
   * system matrices. The parameter ell [3] can be specified.
   */

  class bicgstab_ell_type
  {
    int ell;
  public:
    bicgstab_ell_type(std::string name) : ell(3)
    {
      Parameters::get(name + "->ell", ell);
    }
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::bicgstab_ell(A, x, b, l, r, iter, ell);
    }
  };
  using BiCGStabEllSolver = ITL_Solver<bicgstab_ell_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::GMResSolver
   * \brief ITL_Solver <\ref gmres_type> implementation of generalized minimal 
   * residual method \implements ITL_Solver
   *
   * Solves a linear system by the GMRES method.
   * The parameter restart [30] is the maximal number of orthogonalized vectors.
   * The method is not preconditioned
   */

  enum ORTHOGONALIZATION
  {
    GRAM_SCHMIDT = 1,
    HOUSEHOLDER = 2
  };

  class gmres_type
  {
    int restart;
    int ortho;

  public:
    gmres_type(std::string name) : restart(30), ortho(GRAM_SCHMIDT)
    {
      Parameters::get(name + "->restart", restart);
      Parameters::get(name + "->orthogonalization", ortho);
    }
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      switch ((ORTHOGONALIZATION)ortho)
      {
      default:
      case GRAM_SCHMIDT:
        return itl::gmres2(A, x, b, l, r, iter, restart);
        break;
#ifndef HAVE_PARALLEL_MTL4
      case HOUSEHOLDER:
        return itl::gmres_householder(A, x, b, l, iter, restart);
        break;
#endif
      }
    }
  };
  using GMResSolver = ITL_Solver<gmres_type>;


  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::IDRsSolver
   * \brief ITL_Solver <\ref idr_s_type> implementation of Induced Dimension 
   * Reduction method \implements ITL_Solver
   *
   * Solves a linear system by an Induced Dimension Reduction method and can be 
   * used for system matrices.  The parameter s [30] can be specified.
   * 
   * Peter Sonneveld and Martin B. van Gijzen, IDR(s): a family of simple and fast 
   * algorithms for solving large nonsymmetric linear systems.
   * SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035-1062 (2008). (copyright SIAM)
   */

  class idr_s_type
  {
    int s;
  public:
    idr_s_type(std::string name) : s(30)
    {
      Parameters::get(name + "->s", s);
    }
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::idr_s(A, x, b, l, r, iter, s);
    }
  };
  using IDRsSolver = ITL_Solver<idr_s_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::MinResSolver
   * \brief ITL_Solver <\ref minres_solver_type> implementation of minimal 
   * residual method \implements ITL_Solver
   *
   * Solves a linear system by the Minres method. Can be used for symmetric
   * indefinite systems.
   */

  class minres_solver_type
  {
  public:
    minres_solver_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::minres(A, x, b, l, r, iter);
    }
  };
  using MinResSolver = ITL_Solver<minres_solver_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::GcrSolver
   * \brief ITL_Solver <\ref gcr_type> implementation of generalized conjugate 
   * residual method \implements ITL_Solver
   *
   * Solves a linear system by the GCR method - generalized conjugate residual 
   * method. The parameter restart [30] is the maximal number of orthogonalized 
   * vectors.
   */

  class gcr_type
  {
    int restart;

  public:
    gcr_type(std::string name) : restart(30)
    {
      Parameters::get(name + "->restart", restart);
    }
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::gcr(A, x, b, l, r, iter, restart);
    }
  };
  using GcrSolver = ITL_Solver<gcr_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::FGMResSolver
   * \brief ITL_Solver <\ref fgmres_type> implementation of flexible GMRes method 
   * \implements ITL_Solver
   *
   * Solves a linear system by the FGMRES method.
   * The parameter restart [30] is the maximal number of orthogonalized vectors.
   * See reference "A Flexible Inner-Outer Preconditiones GMRES Algorithm", 
   * Youcef Saad, (1993)
   */

  class fgmres_type
  {
    int restart;
    int orthogonalization;

  public:
    fgmres_type(std::string name) : restart(30), orthogonalization(GRAM_SCHMIDT)
    {
      Parameters::get(name + "->restart", restart);
      Parameters::get(name + "->orthogonalization", orthogonalization);
    }
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      switch ((ORTHOGONALIZATION)orthogonalization)
      {
      default:
      case GRAM_SCHMIDT:
        return itl::fgmres(A, x, b, l, r, iter, restart);
        break;
#ifndef HAVE_PARALLEL_MTL4
      case HOUSEHOLDER:
        return itl::fgmres_householder(A, x, b, r, iter, restart);
        break;
#endif
      }
    }
  };
  using FGMResSolver = ITL_Solver<fgmres_type>;

  // ===========================================================================

  /**
   * \ingroup Solver
   * \class AMDiS::PreOnly
   * \brief ITL_Solver <\ref preonly_type> implementation of preconditioner as 
   * \implements ITL_Solver
   *
   * Solves a linear system by applying a preconditioner only.
   */
  class preonly_type
  {
  public:
    preonly_type(std::string name) {}
    template <class LinOp, class X, class B, class L, class R, class I>
    int operator()(LinOp const& A, X& x, B const& b, L const& l, R const& r, I& iter)
    {
      return itl::preonly(A, x, b, l, iter);
    }
  };
  using PreOnly = ITL_Solver<preonly_type>;

  // ===========================================================================

} // end namespace AMDiS
