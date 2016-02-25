#pragma once

#ifdef _OPENMP
#include <omp.h>
#include <vector>
#endif

namespace AMDiS
{

#ifdef _OPENMP
  template <class T>
  class ThreadPrivate
  {
  public:
    ThreadPrivate()
      : data(omp_get_max_threads())
    {}

    explicit ThreadPrivate(T val)
      : data(omp_get_max_threads(), val)
    {}

    T& get()
    {
#ifndef NDEBUG
      if (omp_get_thread_num() >= data.size())
      {
        std::cout << "Error in ThreadPrivate::get()!\n";
        exit(0);
      }
#endif
      return data[omp_get_thread_num()];
    }

    void set(T& val)
    {
#ifndef NDEBUG
      if (omp_get_thread_num() >= data.size())
      {
        std::cout << "Error in ThreadPrivate::set()!\n";
        exit(0);
      }
#endif
      data[omp_get_thread_num()] = val;
    }

  private:
    std::vector<T> data;
  };

#else

  template <class T>
  class ThreadPrivate
  {
  public:
    ThreadPrivate() {}

    explicit ThreadPrivate(T val)
      : data(val)
    {}

    T& get()
    {
      return data;
    }

    void set(T& val)
    {
      data = val;
    }

  private:
    T data;
  };
#endif

} // end namespace AMDiS
