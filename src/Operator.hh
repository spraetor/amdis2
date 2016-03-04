// #include "Expressions.h"

namespace AMDiS
{
  template <class C>
  void Operator::addZOTImpl(C const& c)
  {
    GenericOperatorTerm<C, 0>* operatorTerm = new GenericZeroOrderTerm<C>(c);
    
    zeroOrder.push_back(operatorTerm);
    operatorTerm->setOperator(this);
    operatorTerm->term.insertFeSpaces(auxFeSpaces);
  }


  template <class B>
  void Operator::addFOTImpl(tag::scalar, B const& b, FirstOrderType type, int i)
  {
    GenericOperatorTerm<B, 1>* operatorTerm;
    if (i >= 0)
      operatorTerm = new GenericFirstOrderTerm_i<-1,B>(b, i);
    else
      operatorTerm = new GenericFirstOrderTerm_1<B>(b);

    if (type == GRD_PSI)
      firstOrderGrdPsi.push_back(operatorTerm);
    else
      firstOrderGrdPhi.push_back(operatorTerm);

    operatorTerm->setOperator(this);
    operatorTerm->term.insertFeSpaces(auxFeSpaces);
  }


  template <class B>
  void Operator::addFOTImpl(tag::vector, B const& b, FirstOrderType type, int /*i*/)
  {
    GenericOperatorTerm<B, 1>* operatorTerm = new GenericFirstOrderTerm_b<B>(b);

    if (type == GRD_PSI)
      firstOrderGrdPsi.push_back(operatorTerm);
    else
      firstOrderGrdPhi.push_back(operatorTerm);

    operatorTerm->setOperator(this);
    operatorTerm->term.insertFeSpaces(auxFeSpaces);
  }


  template <class A>
  void Operator::addSOTImpl(tag::scalar, A const& a, int i, int j, bool /*sym*/)
  {
    GenericOperatorTerm<A, 2>* operatorTerm;
    if (i >= 0 && j >= 0)
      operatorTerm = new GenericSecondOrderTerm_ij<-1,-1,A>(a, i, j);
    else
      operatorTerm = new GenericSecondOrderTerm_1<A>(a);

    secondOrder.push_back(operatorTerm);
    operatorTerm->setOperator(this);
    operatorTerm->term.insertFeSpaces(auxFeSpaces);
  }


  template <class A>
  void Operator::addSOTImpl(tag::matrix, A const& a, int /*i*/, int /*j*/, bool sym)
  {
    GenericOperatorTerm<A, 2>* operatorTerm;
    if (sym)
      operatorTerm = new GenericSecondOrderTerm_A<A, true>(a);
    else
      operatorTerm = new GenericSecondOrderTerm_A<A, false>(a);

    secondOrder.push_back(operatorTerm);
    operatorTerm->setOperator(this);
    operatorTerm->term.insertFeSpaces(auxFeSpaces);
  }

} // end namespace AMDiS
