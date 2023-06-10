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

  //! \brief function to compute the value of a BSpline basis
  //! function at a given point. 

  //! @param u point where the derivative is to be evaluated.
  //! @param p basis function degree.
  //! @param Ubegin iterator to the beginning of the local knot vector.
  //! @param Uend iterator past the end of the local knot vector.

  template<typename INT, typename FLT, typename ITERATOR>
  FLT
  onebasisfun (FLT const u, INT const p, ITERATOR const Ubegin, ITERATOR const Uend) {
    FLT N{0.0};
    const FLT Umin = *(std::min_element (Ubegin, Uend));
    const FLT Umax = *(std::max_element (Ubegin, Uend));

    if (u >= Umin && u <= Umax) {

      if (p == 0) {
	N = 1.0;
      }
      
      else if (p == 1) {
	if (u < *std::next(Ubegin, 1)) {
	  N = ratio ((u - *Ubegin),
		     (*std::next(Ubegin) - *Ubegin));
	}
	else {
	  N = ratio ((*std::next(Ubegin, 2) - u),
		     (*std::next(Ubegin, 2) - *std::next(Ubegin)));
	}
      }
      
      else if (p == 2) {
	const FLT ln = u - *Ubegin;
	const FLT dn = *std::next(Ubegin, 3) - u;
	const FLT ld = *std::next(Ubegin, 2) - *Ubegin; 
	const FLT dd = *std::next(Ubegin, 3) - *std::next(Ubegin, 1);
	if (u < *std::next(Ubegin, 1)) {
	  N = ratio (ln*ln, (ld * (*std::next(Ubegin, 1) - *Ubegin)));
	}
	else if (u > *std::next(Ubegin, 2)) {
	  N = ratio (dn*dn,
		     (dd * (*std::next(Ubegin, 3) - *std::next(Ubegin, 2))));
	}
	else {
	  if (ld > FLT(0.0)) {
	    N += ratio (ln * (*std::next(Ubegin, 2) - u),
			((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * ld));
	  }
	  if (dd > FLT(0.0)) {
	    N += ratio (dn * (u - *std::next(Ubegin, 1)),
			((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * dd));
	  }
	}
      }

      else {
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
      }
    }
    return N;
  };

  //! \brief function to compute the derivative of a BSpline basis function.

  //! @param u point where the derivative is to be evaluated.
  //! @param p basis function degree.
  //! @param Ubegin iterator to the beginning of the local knot vector.
  //! @param Uend iterator past the end of the local knot vector.
  template<typename INT, typename FLT, typename ITERATOR>
  FLT
  onebasisfunder (FLT u, INT p, ITERATOR Ubegin, ITERATOR Uend)
  {

    FLT Nder{0.0};
    const FLT Umin = *(std::min_element (Ubegin, Uend));
    const FLT Umax = *(std::max_element (Ubegin, Uend));
    
    if ((u >= Umin) && ( u <= Umax)) {
   
      if (p == 0) {
	Nder = FLT(0.0);
      }
      else {    

	const FLT ld = *std::next (Uend, -2) - *Ubegin;
	if (std::abs (ld) > FLT(0.0)) {
	  Nder += ratio (p * onebasisfun (u, p-1, Ubegin, std::next (Uend, -1)), ld);
	}
	
	const FLT dd = *std::next (Uend, -1) - *std::next (Ubegin, 1);
	if (std::abs (dd) > FLT(0.0)) { 
	  Nder -= ratio (p * onebasisfun (u, p-1, std::next (Ubegin, 1), Uend), dd);
	}
	
      }
    }

    return Nder;
  };

   //! \brief create an uniform open knot vector given the list
  //! of breaks, a (scalar) degree p and a (scalar) regularity r.

  //! @param bb begin of list of breaks.
  //! @param be end of list of breaks. 
  //! @param p basis function degree.
  //! @param r basis functions regularity.
  //! @return k knoct vector.

  template<typename INT, typename ITERATOR, typename FLT=double>
  std::vector<FLT>
  open_knot_vector (ITERATOR const bb, ITERATOR const be, INT const p, INT const r) {

    const INT Nb = std::distance (bb, be);
    const INT Nk = (Nb - 2) + 2 * (p + 1);

    std::vector<FLT> k (Nk, *bb);

    auto bi = std::next (bb);
    auto ki = std::copy (bi, be, std::next (k.begin (), p + 1));
    std::fill (ki, k.end (), *(std::next (be, -1)));

    return k;
  };

  
  template<typename INT, typename FLT, typename ITERATOR>
  FLT
  onebasisfun2d (FLT const u, FLT const v, INT const pu, INT const pv,
		 ITERATOR const Ubegin, ITERATOR const Uend,
		 ITERATOR const Vbegin, ITERATOR const Vend) {
    FLT N = onebasisfun (u, pu, Ubegin, Uend);
    N *= onebasisfun (v, pv, Vbegin, Vend);
    return N;
  }
  
}

#endif
  

