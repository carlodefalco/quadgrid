#ifndef MY_ATOMIC_ADD
#define MY_ATOMIC_ADD

//requires c++20

#include <atomic>

template <typename T>
void atomicAdd (T *x, const T y)
{
  std::atomic_ref<T> x_ref(*x);
  x_ref += y;
}

#endif
