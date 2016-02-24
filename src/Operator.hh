// #include "Expressions.h"

namespace AMDiS
{
  template <class C>
  void Operator::addZOTImpl(C&& c)
  {
    using C_ = Decay_t<C>;
    OperatorTerm* term = new GenericZeroOrderTerm<C_>(std::forward<C>(c));
    zeroOrder.push_back(term);
    term->setOperator(this);
    c.insertFeSpaces(auxFeSpaces);
  }


  template <class B>
  void Operator::addFOTImpl(tag::scalar, B&& b, FirstOrderType type, int i)
  {
    using B_ = Decay_t<B>;
    OperatorTerm* term;
    if (i >= 0)
      term = new GenericFirstOrderTerm_i<-1,B_>(std::forward<B>(b), i);
    else
      term = new GenericFirstOrderTerm_1<B_>(std::forward<B>(b));

    if (type == GRD_PSI)
      firstOrderGrdPsi.push_back(term);
    else
      firstOrderGrdPhi.push_back(term);

    term->setOperator(this);
    b.insertFeSpaces(auxFeSpaces);
  }


  template <class B>
  void Operator::addFOTImpl(tag::vector, B&& b, FirstOrderType type, int /*i*/)
  {
    using B_ = Decay_t<B>;
    OperatorTerm* term = new GenericFirstOrderTerm_b<B_>(std::forward<B>(b));

    if (type == GRD_PSI)
      firstOrderGrdPsi.push_back(term);
    else
      firstOrderGrdPhi.push_back(term);

    term->setOperator(this);
    b.insertFeSpaces(auxFeSpaces);
  }


  template <class A>
  void Operator::addSOTImpl(tag::scalar, A&& a, int i, int j, bool /*sym*/)
  {
    using A_ = Decay_t<A>;
    OperatorTerm* term;
    if (i >= 0 && j >= 0)
      term = new GenericSecondOrderTerm_ij<-1,-1,A_>(std::forward<A>(a), i, j);
    else
      term = new GenericSecondOrderTerm_1<A_>(std::forward<A>(a));

    secondOrder.push_back(term);
    term->setOperator(this);
    a.insertFeSpaces(auxFeSpaces);
  }


  template <class A>
  void Operator::addSOTImpl(tag::matrix, A&& a, int /*i*/, int /*j*/, bool sym)
  {
    using A_ = Decay_t<A>;
    OperatorTerm* term;
    if (sym)
      term = new GenericSecondOrderTerm_A<A_, true>(std::forward<A>(a));
    else
      term = new GenericSecondOrderTerm_A<A_, false>(std::forward<A>(a));

    secondOrder.push_back(term);
    term->setOperator(this);
    a.insertFeSpaces(auxFeSpaces);
  }

} // end namespace AMDiS
