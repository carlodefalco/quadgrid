#ifndef HAVE_TBASISFUN_H
#define HAVE_TBASISFUN_H

#include <algorithm>
#include <iterator>


template<typename INT, typename FLT, typename ITERATOR>
FLT
onebasisfun (FLT const u, INT const p, ITERATOR const Ubegin, ITERATOR const Uend) {
  FLT N{0.0};
  const FLT Umin = *(std::min_element (Ubegin, Uend));
  const FLT Umax = *(std::max_element (Ubegin, Uend));

  if ((u <= Umin) || ( u > Umax)) {
    N = 0.0;
    return N;
  } else if (p == 0) {
    N = 1.0;
    return N;
  } else if (p == 1) {
    if (u < *std::next(Ubegin, 1)) {
      N = (u - *Ubegin) / (*std::next(Ubegin, 1) - *Ubegin);
      return N;
    } else {
      N = (*std::next(Ubegin, 2) - u) / (*std::next(Ubegin, 2) - *std::next(Ubegin, 1));
      return N;
    }
  } else if (p == 2) {
    const FLT ln = u - *Ubegin;
    const FLT dn = *std::next(Ubegin, 3) - u;
    const FLT ld = *std::next(Ubegin, 2) - *Ubegin; 
    const FLT dd = *std::next(Ubegin, 3) - *std::next(Ubegin, 1);
    if (u < *std::next(Ubegin, 1)) {
      N = ln*ln / (ld * (*std::next(Ubegin, 1) - *Ubegin));
      return N;
    } else if (u > *std::next(Ubegin, 2)) {
      N = dn*dn / (dd * (*std::next(Ubegin, 3) - *std::next(Ubegin, 2)));
      return N;
    } else {
      if (ld != FLT(0.0)) {
	N += ln * (*std::next(Ubegin, 2) - u) / ((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * ld);
      }
      if (dd != FLT(0.0)) {
	N += dn * (u - *std::next(Ubegin, 1)) / ((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * dd);
      }
      return N;
    }
  }

  const FLT ln = u - *Ubegin;
  const FLT ld = *std::next(Uend, - 2) - *Ubegin;
  if (ld != FLT(0.0)) {
    N += ln * onebasisfun (u, p-1, Ubegin, std::next(Uend, - 1)) / ld; 
  }
  
  const FLT dn = *std::next (Uend, - 1) - u;
  const FLT dd = *std::next (Uend, - 1) - *std::next (Ubegin , 1);
  if (dd != 0) {
    N += dn * onebasisfun (u, p-1, std::next (Ubegin, 1), Uend) / dd;
  }
  
  return N;
};

#endif
  

