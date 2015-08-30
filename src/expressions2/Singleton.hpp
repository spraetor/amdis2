/** \file Singleton.hpp */

#pragma once

#include <mutex>
#include <memory>

namespace AMDiS
{
  template <class T>
  class Singleton
  {
  private:
    static std::shared_ptr<T> instance;
    static std::mutex         mutex_holder;

    Singleton() = delete;
    Singleton(Singleton const& rhs) = delete;
    Singleton& operator=(Singleton const& rhs) = delete;

  public:
    ~Singleton() {}

    template <class... Args>
    static T& getInstance(Args&& ... args)
    {
      if (!Singleton::instance)
      {
        std::lock_guard<std::mutex> lock(mutex_holder);

        if (!instance)
          instance.reset(new T(std::forward<Args>(args)...));
      }
      return *Singleton::instance;
    }
  };

  std::mutex Singleton::mutex_holder;
  template <class T> template* Singleton<T>::instance = nullptr;

} // end namespace AMDiS
