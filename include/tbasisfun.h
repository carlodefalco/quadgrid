#ifndef HAVE_TBASISFUN_H
#define HAVE_TBASISFUN_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>

namespace bspline {
  
  template<typename FLT>
  FLT
  ratio (const FLT num, const FLT denom) {
    if (denom != FLT(0.0)) {
      return num / denom;
    } else if (num != FLT(0.0)) {
      throw ("division by zero!");
    } else {
      return FLT(1.0);
    }
  };
 
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
	N = ratio ((u - *Ubegin),
		   (*std::next(Ubegin, 1) - *Ubegin));
	return N;
      } else {
	N = ratio ((*std::next(Ubegin, 2) - u),
		   (*std::next(Ubegin, 2) - *std::next(Ubegin, 1)));
	return N;
      }
    } else if (p == 2) {
      const FLT ln = u - *Ubegin;
      const FLT dn = *std::next(Ubegin, 3) - u;
      const FLT ld = *std::next(Ubegin, 2) - *Ubegin; 
      const FLT dd = *std::next(Ubegin, 3) - *std::next(Ubegin, 1);
      if (u < *std::next(Ubegin, 1)) {
	N = ratio (ln*ln, (ld * (*std::next(Ubegin, 1) - *Ubegin)));
	return N;
      } else if (u > *std::next(Ubegin, 2)) {
	N = ratio (dn*dn,
		   (dd * (*std::next(Ubegin, 3) - *std::next(Ubegin, 2))));
	return N;
      } else {
	if (ld > FLT(0.0)) {
	  N += ratio (ln * (*std::next(Ubegin, 2) - u),
		      ((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * ld));
	}
	if (dd > FLT(0.0)) {
	  N += ratio (dn * (u - *std::next(Ubegin, 1)),
		      ((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * dd));
	}
	return N;
      }
    }

    const FLT ln = u - *Ubegin;
    const FLT ld = *std::next(Uend, - 2) - *Ubegin;
    if (ld > FLT(0.0)) {
      N += ln * onebasisfun (u, p-1, Ubegin, std::next(Uend, - 1)) / ld; 
    }
  
    const FLT dn = *std::next (Uend, - 1) - u;
    const FLT dd = *std::next (Uend, - 1) - *std::next (Ubegin , 1);
    if (dd > FLT(0.0)) {
      N += dn * onebasisfun (u, p-1, std::next (Ubegin, 1), Uend) / dd;
    }
  
    return N;
  };

  template<typename INT, typename FLT, typename ITERATOR>
  FLT
  onebasisfunder (FLT u, INT p, ITERATOR Ubegin, ITERATOR Uend)
  {

    FLT Nder{0.0};
    const FLT Umin = *(std::min_element (Ubegin, Uend));
    const FLT Umax = *(std::max_element (Ubegin, Uend));
    
    if ((u <= Umin) || ( u > Umax)) {
      Nder = FLT(0.0);
      return Nder;
    } else if (p == 0) {
      Nder = FLT(0.0);
      return Nder;
    } else {    

      const FLT ld = *std::next (Uend, -2) - *Ubegin;
      if (ld != FLT(0.0)) {
	Nder += p * onebasisfun (u, p-1, Ubegin, std::next (Uend, -1)) / ld;
      }
    
      const FLT dd = *std::next (Uend, -1) - *std::next (Ubegin, 1);
      if (dd != FLT(0.0)) { 
	Nder -= p * onebasisfun (u, p-1, std::next (Ubegin, 1), Uend) / dd;
      }
      
    }

    return Nder;
  };
  
}

#endif
  

