#ifndef HAVE_TBASISFUN_H
#define HAVE_TBASISFUN_H

#include <algorithm>
#include <iomanip>
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
      return FLT(0.0);
    }
  };

  //! \brief function to compute the value of a BSpline basis
  //! function at a given point. 

  //! @param u point where the function is to be evaluated.
  //! @param p basis function degree.
  //! @param Ubegin iterator to the beginning of the local knot vector.
  //! @param Uend iterator past the end of the local knot vector.

  template<typename INT, typename FLT, typename ITERATOR>
  FLT
  onebasisfun (FLT const u, INT const p, ITERATOR const Ubegin, ITERATOR const Uend) {
    FLT N{0.0};
    const FLT Umin = *(std::min_element (Ubegin, Uend));
    const FLT Umax = *(std::max_element (Ubegin, Uend));

    //std::cout << std::setprecision (16) <<"Umin = " << Umin << " u = " << u << " Umax = " << Umax << std::endl;
    if (u >= Umin && u < Umax) {

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
	const FLT dd = *std::next(Ubegin, 3) - *std::next(Ubegin);
	  
	if (u < *std::next (Ubegin)) {
	  N = ratio (ln*ln, ld * (*std::next (Ubegin) - *Ubegin));
	}
	else if (u >= *std::next (Ubegin, 2)) {
	  N = ratio (dn*dn,
		     (dd * (*std::next (Ubegin, 3) - *std::next (Ubegin, 2))));
	}
	else {

	  if (u>=*std::next(Ubegin,1) && u < *std::next(Ubegin,2)) {
	    N += ratio (ln * (*std::next(Ubegin, 2) - u),
			((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * ld));
	  
	    N += ratio (dn * (u - *std::next(Ubegin, 1)),
			((*std::next(Ubegin, 2) - *std::next(Ubegin, 1)) * dd));
	  }
	}
      }

      else {
	
	const FLT ld = *std::next (Uend, - 2) - *Ubegin;
	const FLT dd = *std::prev (Uend) - *std::next (Ubegin);
	
	if (u<*std::next(Ubegin,p)) {
	  const FLT ln = u - *Ubegin;
	  N += ratio(ln * onebasisfun (u, p-1, Ubegin, std::prev (Uend)) , ld); 
	}
	
	if (u>=*std::next(Ubegin,1)) {
	  const FLT dn = *std::prev (Uend) - u;
	  N += ratio(dn * onebasisfun (u, p-1, std::next (Ubegin), Uend) , dd);
	}

	//std::cout << "p = " << p << std::endl;
	//std::cout << "u - *Ubegin = " << u - *Ubegin << " ld = " << ld << " ";

	//std::cout << " *std::next (Uend, - 1) - *std::next (Ubegin , 1) = "
	//	  << *std::next (Uend, - 1) - *std::next (Ubegin , 1)
	//	  << " dd = " << dd << " N = " << N << std::endl;
      }
    }

    //std::cout << " N = " << N << std::endl;
    return N;
  };


 //! \brief function to compute the value of a BSpline basis
  //! function at a given point. 

  //! @param u point where the function is to be evaluated.
  //! @param p basis function degree.
  //! @param U knot-vector
  //! @param i index of the BSpline basis funct


  template<typename INT, typename FLT>
  FLT
  onebasisfun (FLT const u, INT const p, std::vector<FLT> const & U, INT i) {
    FLT N{0.0};
  
  if((u==U.back() && i==(static_cast<INT>(U.size())-p-2)) || // added static_cast on size 
  (u==U[0] && i==0)){// 2nd case not needed

  return FLT{1.0};
  }

  //std::cout << std::setprecision (16) <<"Umin = " << Umin << " u = " << u << " Umax = " << Umax << std::endl;
  if (u >= U[i] && u < U[i+p+1]) {


    if (p == 0) {
        N = 1.0;
        }
        
    else if (p == 1) {
      if (u < U[i+1]) {
        N = ratio ((u - U[i]),
            (U[i+1] - U[i]));
      }
      else {
        N = ratio ((U[i+2] - u),
            (U[i+2] - U[i+1]));
      }
        }
        
    else if (p == 2) {
    
      const FLT ln = u - U[i];
      const FLT dn = U[i+3] - u;
      const FLT ld = U[i+2] - U[i]; 
      const FLT dd = U[i+3] - U[i+1];
        
      if (u <U[i+1]) {
        N = ratio (ln*ln, ld * (U[i+1] - U[i]));
      }
      else if (u >= U[i+2]) {
        N = ratio (dn*dn,
            (dd * (U[i+3] - U[i+2])));
      }
      else {

        if (u>=U[i+1] && u <U[i+2]) {
          N += ratio (ln * (U[i+2] - u),
          ((U[i+2] - U[i+1]) * ld));
        
          N += ratio (dn * (u - U[i+1]),
          ((U[i+2] - U[i+1]) * dd));
        }
      }
    }

    else {
    
      const FLT ld = U[i+p] - U[i];
      const FLT dd = U[i+p+1] - U[i+1];
      
      if (u<U[i+p]) {
        const FLT ln = u - U[i];
        N += ratio(ln * onebasisfun (u, p-1, U , i ) , ld); 
      }
      
      if (u>=U[i+1]) {
        const FLT dn = U[i+p+1] - u;
        N += ratio(dn * onebasisfun (u, p-1, U, i+1 ) , dd);
      }

    }
}

  return N;
};

// Recursive version


  template<typename INT, typename FLT, int P>
  FLT onebasisfun (FLT const u, std::vector<FLT> const & U, INT i) {
    FLT N{0.0};
  
    if((u==U.back() && i==(static_cast<INT>(U.size())-P-2)) || // added static_cast on size 
      (u==U[0] && i==0)){// 2nd case not needed

      return FLT{1.0};
    }

    if (u >= U[i] && u < U[i+P+1]) {

      if constexpr (P == 0) {
        return FLT{1.0};
        }

      else{
      const FLT ld = U[i+P] - U[i];
      const FLT dd = U[i+P+1] - U[i+1];
      
      if (u<U[i+P]) {
        const FLT ln = u - U[i];
        N += ratio(ln * onebasisfun<INT,FLT,P-1> (u, U , i ) , ld); 
      }
      
      if (u>=U[i+1]) {
        const FLT dn = U[i+P+1] - u;
        N += ratio(dn * onebasisfun<INT,FLT,P-1> (u, U, i+1 ) , dd);
      }
    }
    }
      return N;

    }





  //! \brief function to compute the derivative of a BSpline basis function.

  //! @todo make shure that evaluation at the endopoints works correctly.
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
  // std::cout<<  onebasisfun (u, p-1, std::next (Ubegin, 1), Uend)<<std::endl;
	// 	std::cout<<"deriv in u: "<<u<<"è "<<Nder<<std::endl;
      }
    }

    return Nder;
  };


// New version

 //! \brief function to compute the value of the derivative of a BSpline basis
  //! function at a given point. 

  //! @param u point where the derivative has to be evaluated.
  //! @param p basis function degree.
  //! @param U knot-vector
  //! @param i index of the BSpline basis funct to be derived

  template<typename INT, typename FLT>
  FLT
  onebasisfunder (FLT u, INT p, std::vector<FLT> const & U, INT i)
  {

    FLT Nder{0.0};

    if (u==U[i+p+1] && i==(static_cast<INT>(U.size())-p-2)){
        double tol=1e-10;
         u=u-tol;
    }
  
    if ((u >=U[i] ) && ( u < U[i+p+1])) {
   
      if (p == 0) {
	Nder = FLT(0.0);
      }
      else {    

	const FLT ld = U[i+p] - U[i];
	if (u<U[i+p]) {
	  Nder += ratio (p * onebasisfun (u, p-1, U, i), ld);
	}
	const FLT dd = U[i+p+1] - U[i+1];
	if (u>=U[i+1]) { 
	  Nder -= ratio (p * onebasisfun (u, p-1, U, i+1), dd);
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
    const INT Nk = (p - r) * (Nb - 2) + 2 * (p + 1); // (Nb - 2 + 2) *p 

    std::vector<FLT> k (Nk, *bb);

    
    auto ki = std::next (k.begin (), p + 1);
    for (auto bi = std::next (bb); bi != be; bi = std::next (bi)) {
      ki = std::fill_n (ki, p - r, *bi);
    }    
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

  template<typename INT, typename FLT>
  FLT
  onebasisfun2d (FLT const u, FLT const v, INT const pu, INT const pv,
           std::vector<FLT> const & U, INT iu,
           std::vector<FLT> const & V, INT iv) {
    FLT Nu = onebasisfun<INT, FLT>(u, pu, U, iu);
    FLT Nv = onebasisfun<INT, FLT>(v, pv, V, iv);
    return Nu * Nv;
  };
  
}

#endif
  

