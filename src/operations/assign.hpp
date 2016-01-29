#pragma once

namespace AMDiS
{
  namespace assign
  {

    // -------------------------------------------------------------------------
    // unary assigners

    /// assign::value(v) --> v = val
    template <class T, class S=T>
    struct value
    {
      using result_type = T;
      value(S val = 0) : val(val) {}

      T& operator()(T& v) const
      {
        return (v = val);
      }

    private:
      S val;
    };

    /// assign::ct_value(v) --> v = Val
    template <class T, class S, S Val>
    struct ct_value : value<T, S>
    {
      ct_value() : value<T, S>(Val) {}
    };

    /// assign::min_value(v) --> v = min<S>
    template <class T, class S=T>
    struct min_value : value<T, S>
    {
      min_value() : value<T, S>(std::numeric_limits<S>::min()) {}
    };

    /// assign::max_value(v) --> v = max<S>
    template <class T, class S=T>
    struct max_value : value<T, S>
    {
      max_value() : value<T, S>(std::numeric_limits<S>::max()) {}
    };

    /// assign::plus_value(v) --> v += val
    template <class T, class S=T>
    struct plus_value
    {
      using result_type = T;
      plus_value(S val) : val(val) {}

      T& operator()(T& v) const
      {
        return (v += val);
      }

    private:
      S val;
    };

    /// assign::minus_value(v) --> v -= val
    template <class T, class S=T>
    struct minus_value
    {
      using result_type = T;
      minus_value(S val) : val(val) {}

      T& operator()(T& v) const
      {
        return (v -= val);
      }

    private:
      S val;
    };

    /// assign::mult_value(v) --> v *= val
    template <class T, class S=T>
    struct mult_value
    {
      using result_type = T;
      mult_value(S val) : val(val) {}

      T& operator()(T& v) const
      {
        return (v *= val);
      }

    private:
      S val;
    };

    /// assign::div_value(v) --> v /= val
    template <class T, class S=T>
    struct div_value
    {
      using result_type = T;
      div_value(S val) : val(val) {}

      T& operator()(T& v)
      {
        return (v /= val);
      }

    private:
      S val;
    };

    
    // -------------------------------------------------------------------------
    /// binary assigners


    /// functor for operator=
    template <class T>
    struct assign
    {
      using result_type = T;

      constexpr static T& apply(T& v, T const& v0)
      {
        return (v = v0);
      }
      T& operator()(T& v, T const& v0) const
      {
        return (v = v0);
      }
    };

    /// functor for operator+=
    template <class T>
    struct plus
    {
      using result_type = T;

      constexpr static T& apply(T& v, T const& v0)
      {
        return (v += v0);
      }
      T& operator()(T& v, T const& v0) const
      {
        return (v += v0);
      }
    };

    /// functor for operator*=
    template <class T>
    struct multiplies
    {
      using result_type = T;

      constexpr static T& apply(T& v, T const& v0)
      {
        return (v *= v0);
      }
      T& operator()(T& v, T const& v0)
      {
        return (v *= v0);
      }
    };

    /// functor for operator/=
    template <class T>
    struct divides
    {
      using result_type = T;

      constexpr static T& apply(T& v, T const& v0)
      {
        return (v /= v0);
      }
      T& operator()(T& v, T const& v0)
      {
        return (v /= v0);
      }
    };

    /// functor for v = max(v, v0)
    template <class T>
    struct max
    {
      using result_type = T;

      constexpr static T& apply(T& v, T const& v0)
      {
        return v = std::max(v, v0);
      }
      T& operator()(T& v, T const& v0)
      {
        return v = std::max(v, v0);
      }
    };

    /// functor for v = min(v, v0)
    template <class T>
    struct min
    {
      using result_type = T;

      constexpr static T& apply(T& v, T const& v0)
      {
        return v = std::min(v, v0);
      }
      T& operator()(T& v, T const& v0)
      {
        return v = std::min(v, v0);
      }
    };

  } // end namespace assign
} // end namespace AMDiS
