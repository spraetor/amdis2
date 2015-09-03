/** \file diff_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "simplify_expr.hpp"

namespace AMDiS
{
  // by default the derivative is zero
  namespace expressions
  {
    struct DummyType
    {
      typedef _unknown id;
      typedef double value_type;
    };

    template <class Identifier, class Term, class Direction = DummyType>
    struct Diff
    {
      typedef CValue<0> type;
      typedef CValue<0> dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&&, Direction_&&)
      {
        return {};
      }
    };

  } // end namespace expressions

  template <class Id, class Term>
  inline typename expressions::Diff<Id,Term>::type
  diff(Term&& t)
  {
    return expressions::Diff<Id,Term>::eval(std::forward<Term>(t));
  }

  template <class Id, class Term, class Direction>
  inline typename expressions::Diff<Id,Term,Direction>::dir_type
  diff(Term&& t, Direction&& d)
  {
    return expressions::Diff<Id,Term,Direction>::eval(std::forward<Term>(t), std::forward<Direction>(d));
  }


  // _____________________________________________________________________________
  // value/gradient expressions
  namespace expressions
  {
    // ValueOf
    template <class Id, class Vector, class Name, class Direction>
    struct Diff<Id, ValueOf<Vector, Name>, Direction>
    {
      typedef ValueOf<Vector, Name> Term;
      typedef CValue<0> type;
      typedef CValue<0> dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&&, Direction_&&)
      {
        return {};
      }
    };

    template <class Id, class Vector, class Direction>
    struct Diff<Id, ValueOf<Vector, Id>, Direction>
    {
      typedef expressions::ValueOf<Vector, Id> Term;
      typedef CValue<1> type;
      typedef Direction dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&&, Direction_&& d)
      {
        return d;
      }
    };


    // GradientOf
    template <class Id, class Vector, class Name, class Direction>
    struct Diff<Id, GradientOf<Vector, Name>, Direction>
    {
      typedef GradientOf<Vector, Name> Term;
      typedef CValue<0> type;
      typedef CValue<0> dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&& t, Direction_&& d)
      {
        return {};
      }
    };

    template <class Id, class Vector, class Direction>
    struct Diff<Id, GradientOf<Vector, Id>, Direction>
    {
      typedef GradientOf<Vector, Id> Term;
      typedef CValue<0> type;
      typedef GradientOf<Vector, typename Direction::id> dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&&, Direction_&& d)
      {
        return gradientOf<typename Direction::id>(d.vecDV);
      }
    };


    // DerivativeOf
    template <class Id, int I, class Vector, class Name, class Direction>
    struct Diff<Id, DerivativeOf<I, Vector, Name>, Direction>
    {
      typedef DerivativeOf<I, Vector, Name> Term;
      typedef CValue<0> type;
      typedef CValue<0> dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&& t, Direction_&& d)
      {
        return {};
      }
    };

    template <class Id, int I, class Vector, class Direction>
    struct Diff<Id, DerivativeOf<I, Vector, Id>, Direction>
    {
      typedef DerivativeOf<I, Vector, Id>                      original_type;
      typedef CValue<0>                                        type;
      typedef DerivativeOf<I, Vector, typename Direction::id>  dir_type;

      template <class Term_>
      static type eval(Term_&&)
      {
        return {};
      }

      template <class Term_, class Direction_>
      static dir_type eval(Term_&&, Direction_&& d)
      {
        return derivativeOf<typename Direction::id, I>(d.vecDV);
      }
    };

  } // end namespace expressions


  // _____________________________________________________________________________
  // multiplication expressions
  namespace expressions
  {
    template<class Id, class Term1, class Term2, class Direction>
    class Diff<Id, Mult<Term1, Term2>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term1>::type>::type D1;
      typedef typename Simplify<typename Diff<Id, Term2>::type>::type D2;
      typedef typename Simplify<typename Diff<Id, Term1, Direction>::dir_type>::type D1_;
      typedef typename Simplify<typename Diff<Id, Term2, Direction>::dir_type>::type D2_;

    public:
      typedef Mult<Term1, Term2> original_type;
      typedef typename Simplify<Add<Mult<D1, Term2>,
              Mult<Term1, D2>>>::type type;
      typedef typename Simplify<Add<Mult<D1_, Term2>,
              Mult<Term1, D2_>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * t.getTerm(_1) + t.getTerm(_0) * simplify(diff<Id>(t.getTerm(_1))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0),d)) * t.getTerm(_1) + t.getTerm(_0) * simplify(diff<Id>(t.getTerm(_1),d)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, MultInverse<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef MultInverse<Term> original_type;
      typedef typename Simplify<Mult<Negative<D>,
              MultInverse<Pow<2, Term>>>>::type type;
      typedef typename Simplify<Mult<Negative<D_>,
              MultInverse<Pow<2, Term>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify((-simplify(diff<Id>(t.getTerm(_0)))) / pow<2>(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify((-simplify(diff<Id>(t.getTerm(_0),d))) / pow<2>(t.getTerm(_0)));
      }
    };

  } // end namespace expressions


  // _____________________________________________________________________________
  // add expressions
  namespace expressions
  {

    template<typename Id, typename Term1, typename Term2, typename Direction>
    class Diff<Id, Add<Term1, Term2>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term1>::type>::type D1;
      typedef typename Simplify<typename Diff<Id, Term2>::type>::type D2;
      typedef typename Simplify<typename Diff<Id, Term1, Direction>::dir_type>::type D1_;
      typedef typename Simplify<typename Diff<Id, Term2, Direction>::dir_type>::type D2_;

    public:
      typedef Add<Term1, Term2> original_type;
      typedef typename Simplify<Add<D1, D2>>::type type;
      typedef typename Simplify<Add<D1_, D2_>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) + simplify(diff<Id>(t.getTerm(_1))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0),d)) + simplify(diff<Id>(t.getTerm(_1),d)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    struct Diff<Id, Negative<Term>, Direction>
    {
      typedef Negative<Term> original_type;
      typedef typename Simplify<Negative<typename Diff<Id, Term>::type>>::type type;
      typedef typename Simplify<Negative<typename Diff<Id, Term, Direction>::dir_type>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(-diff<Id>(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(-diff<Id>(t.getTerm(_0),d));
      }
    };

  } // end namespace expressions


  // _____________________________________________________________________________
  // absolute value and signum
  namespace expressions
  {
    template<typename Term, typename Id>
    struct Diff<Id, Abs<Term>> {}; // not yet implemented

  } // end namespace expressions


  // _____________________________________________________________________________
  // powers and roots
  namespace expressions
  {
    template<typename Id, int I, typename Term, typename Direction>
    class Diff<Id, Pow<I, Term>,Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Pow<I, Term> original_type;
      typedef typename Simplify< Mult< Mult<D,  CValue<I>>, Pow<I-1, Term> > >::type type;
      typedef typename Simplify< Mult< Mult<D_,  CValue<I>>, Pow<I-1, Term> > >::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * CValue<I>() * pow<I-1>(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * CValue<I>() * pow<I-1>(t.getTerm(_0)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Pow<1, Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Pow<1, Term> original_type;
      typedef D type;
      typedef D_ dir_type;

      static type eval(original_type const& t)
      {
        return simplify(diff<Id>(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(diff<Id>(t.getTerm(_0), d));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Sqrt<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Sqrt<Term> original_type;
      typedef typename Simplify<Mult<D, MultInverse<Mult<CValue<2>, Sqrt<Term>>>>>::type type;
      typedef typename Simplify<Mult<D_, MultInverse<Mult<CValue<2>, Sqrt<Term>>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) / (CValue<2>() * sqrt(t.getTerm(_0))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) / (CValue<2>() * sqrt(t.getTerm(_0))));
      }
    };

  } // end namespace expressions


  // _____________________________________________________________________________
  // exponential function and logarithms
  namespace expressions
  {
    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Exp<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Exp<Term> original_type;
      typedef typename Simplify<Mult<D, Exp<Term>>>::type type;
      typedef typename Simplify<Mult<D_, Exp<Term>>>::type dir_type;

      static type eval(Exp<Term> const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * exp(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * exp(t.getTerm(_0)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Log<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Log<Term> original_type;
      typedef typename Simplify<Mult<D, MultInverse<Term>>>::type type;
      typedef typename Simplify<Mult<D_, MultInverse<Term>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) / t.getTerm(_0));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) / t.getTerm(_0));
      }
    };

  } // end namespace expressions

  // _____________________________________________________________________________
  // cosine, sine and tangence
  namespace expressions
  {
    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Cos<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Cos<Term> original_type;
      typedef typename Simplify<Negative<Mult<D, Sin<Term>>>>::type type;
      typedef typename Simplify<Negative<Mult<D_, Sin<Term>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(-(simplify(diff<Id>(t.getTerm(_0))) * sin(t.getTerm(_0))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(-(simplify(diff<Id>(t.getTerm(_0), d)) * sin(t.getTerm(_0))));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Sin<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Sin<Term> original_type;
      typedef typename Simplify<Mult<D, Cos<Term>>>::type type;
      typedef typename Simplify<Mult<D_, Cos<Term>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * cos(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * cos(t.getTerm(_0)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Tan<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Tan<Term> original_type;
      typedef typename Simplify<Mult<D, Add<CValue<1>, Negative<Pow<2, Tan<Term>>>>>>::type type;
      typedef typename Simplify<Mult<D_, Add<CValue<1>, Negative<Pow<2, Tan<Term>>>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * (CValue<1>() - pow<2>(tan(t.getTerm(_0)))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * (CValue<1>() - pow<2>(tan(t.getTerm(_0)))));
      }
    };

  } // end namespace expressions

  // _____________________________________________________________________________
  // arcus-(cosine, sine, tangence)
  namespace expressions
  {
    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Acos<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Acos<Term> original_type;
      typedef typename Simplify<Negative<Mult<D, MultInverse<Add<CValue<1>, Negative<Pow<2, Term>>>>>>>::type type;
      typedef typename Simplify<Negative<Mult<D_, MultInverse<Add<CValue<1>, Negative<Pow<2, Term>>>>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(-(simplify(diff<Id>(t.getTerm(_0))) / (CValue<1>() - pow<2>(t.getTerm(_0)))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(-(simplify(diff<Id>(t.getTerm(_0), d)) / (CValue<1>() - pow<2>(t.getTerm(_0)))));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Asin<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Asin<Term> original_type;
      typedef typename Simplify<Mult<D, MultInverse<Sqrt<Add<CValue<1>, Negative<Pow<2, Term>>>>>>>::type type;
      typedef typename Simplify<Mult<D_, MultInverse<Sqrt<Add<CValue<1>, Negative<Pow<2, Term>>>>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) / sqrt(CValue<1>() - pow<2>(t.getTerm(_0))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) / sqrt(CValue<1>() - pow<2>(t.getTerm(_0))));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Atan<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Atan<Term> original_type;
      typedef typename Simplify<Mult<D, MultInverse<Add<CValue<1>, Pow<2, Term>>>>>::type type;
      typedef typename Simplify<Mult<D_, MultInverse<Add<CValue<1>, Pow<2, Term>>>>>::type dir_type;

      static type eval(original_type const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) / (CValue<1>() + pow<2>(t.getTerm(_0))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) / (CValue<1>() + pow<2>(t.getTerm(_0))));
      }
    };

    template<typename Id, typename Term1, typename Term2>
    struct Diff<Id, Atan2<Term1, Term2>> {};

  } // end namespace expressions

  // _____________________________________________________________________________
  // (cose, sine ,tangence)-hyperbolicus
  namespace expressions
  {
    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Cosh<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Cosh<Term> original_type;
      typedef typename Simplify<Mult<D, Sinh<Term>>>::type type;
      typedef typename Simplify<Mult<D_, Sinh<Term>>>::type dir_type;

      static type eval(Cosh<Term> const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * sinh(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * sinh(t.getTerm(_0)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Sinh<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Sinh<Term> original_type;
      typedef typename Simplify<Mult<D, Cosh<Term>>>::type type;
      typedef typename Simplify<Mult<D_, Cosh<Term>>>::type dir_type;

      static type eval(Sinh<Term> const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * cosh(t.getTerm(_0)));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * cosh(t.getTerm(_0)));
      }
    };

    template<typename Id, typename Term, typename Direction>
    class Diff<Id, Tanh<Term>, Direction>
    {
      typedef typename Simplify<typename Diff<Id, Term>::type>::type D;
      typedef typename Simplify<typename Diff<Id, Term, Direction>::dir_type>::type D_;

    public:
      typedef Tanh<Term> original_type;
      typedef typename Simplify<Mult<D, Add<CValue<1>, Negative<Pow<2, Tanh<Term>>>>>>::type type;
      typedef typename Simplify<Mult<D_, Add<CValue<1>, Negative<Pow<2, Tanh<Term>>>>>>::type dir_type;

      static type eval(Tanh<Term> const& t)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0))) * (CValue<1>() - pow<2>(tanh(t.getTerm(_0)))));
      }

      static dir_type eval(original_type const& t, Direction const& d)
      {
        return simplify(simplify(diff<Id>(t.getTerm(_0), d)) * (CValue<1>() - pow<2>(tanh(t.getTerm(_0)))));
      }
    };

  } // end namespace expressions


#if 0
  // _____________________________________________________________________________
  // arcus-(cosine, sine, tangence)-hyperbolicus
  namespace expressions
  {

    template<typename Id, typename Term>
    struct Diff<Id, Acosh<Term>> {}; // not yet implemented

    template<typename Id, typename Term>
    struct Diff<Id, Asinh<Term>> {}; // not yet implemented

    template<typename Id, typename Term>
    struct Diff<Id, Atanh<Term>> {}; // not yet implemented

  } // end namespace expressions


  // _____________________________________________________________________________
  // maximum and mininmum
  namespace expressions
  {

    template<typename Id, typename Term1, typename Term2>
    struct Diff<Id, Max<Term1, Term2>> {};  // not yet implemented

    template<typename Id, typename Term1, typename Term2>
    struct Diff<Id, Min<Term1, Term2>> {};  // not yet implemented

  } // end namespace expressions
#endif

  // =============================================================================
  // higher order derivatives

  // forward declaration
  namespace expressions
  {
    template<int N, class Id, class Term> struct DiffN;

  } // end namespace expressions

  template<int N, typename Id, typename Term>
  inline typename expressions::DiffN<N,Id,Term>::type
  diff(Term&& t)
  {
    return expressions::DiffN<N,Id,Term>::eval(std::forward<Term>(t));
  }

  namespace expressions
  {
    template <int N, class Id, class Term>
    struct DiffN
    {
      typedef typename DiffN<N-1, Id, typename Diff<Id, Term>::type>::type type;

      template <class Term_>
      static type eval(Term_&& t)
      {
        return diff<N-1, Id>(diff<Id>(std::forward<Term_>(t)));
      }
    };

    template <class Id, class Term>
    struct DiffN<1, Id, Term>
    {
      typedef typename Diff<Id, Term>::type type;

      template <class Term_>
      static type eval(Term_&& t)
      {
        return diff<Id>(std::forward<Term_>(t));
      }
    };

    template <class Id, class Term>
    struct DiffN<0, Id, Term>
    {
      typedef Term type;

      template <class Term_>
      static type eval(Term_&& t)
      {
        return t;
      }
    };

  } // end namespace expressions


} // end namespace AMDiS
