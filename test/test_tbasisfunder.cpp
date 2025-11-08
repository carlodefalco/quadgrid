#include <tbasisfun.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace bspline;

int main(int argc, char *argv[]) {

int p=3;
int N=20;
std::vector<double> breakpts{ 1.0, 2.0, 3.0};
std::vector<double> U = open_knot_vector (breakpts.begin (), breakpts.end (), p,0 );

std::cout << "Knot vector: " << std::endl;
for (double const & ii : U)
 std::cout << ii << " ";

std::cout << std::endl;
std::vector<double> x(N+1, 0.0);
double dx=(breakpts.back()-breakpts.front())/double(N);
for (int j=0;j<x.size();j++){
    x[j]=breakpts.front()+j*dx;
}

std::cout<<dx<<std::endl;
std::vector<double> y1(N+1, 0.0);
std::vector<double> y2(N+1, 0.0);
std::cout<<static_cast<int>(U.size())-(p+1)<<std::endl;
std::cout<<" Starting evaluation"<<std::endl;

  for (int ii = 0; ii < static_cast<int>(U.size())-(p+1); ++ii) {
    
    auto kb = std::next (U.begin (), ii);
    auto ke = std::next (kb, p + 2);

    for (int j=0; j < N+1; j++){
      y1[j] = onebasisfunder (x[j], p, kb, ke);
      y2[j]= onebasisfunder<int,double>(x[j],p,U,ii);
      std::cout << std::setprecision (19) << x[j] << "   " << y1[j] <<"   "<<y2[j]<< std::endl;
    }
    std::cout<<"-----------------------"<<ii<<"--------------------"<<std::endl;

  }
  return 0;  
}

