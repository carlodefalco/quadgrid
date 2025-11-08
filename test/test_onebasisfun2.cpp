#include <tbasisfun.h>
#include <iostream>
#include <vector>




using namespace bspline;

int main(int argc, char *argv[]) {

int p=3;
int N=20;
std::vector<double> breakpts{ 1.0, 2.0, 3.0};
std::vector<double> k = open_knot_vector (breakpts.begin (), breakpts.end (), p, p-3);

std::cout << "Knot vector: " << std::endl;
for (double const & ii : k)
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
std::vector<double> y3(N+1, 0.0);
std::cout<<static_cast<int>(k.size())-(p+1)<<std::endl;
std::cout<<" Starting evaluation"<<std::endl;

for (int i=0; i< static_cast<int>(k.size())-(p+1);i++){

    auto kb = std::next (k.begin (), i);
    auto ke = std::next (kb, p + 2);// iterator past the end of the local knot vector.
    
    for (int j=0; j < N+1; j++){
        y1[j]=onebasisfun(x[j],p,kb,ke);
        y2[j]= onebasisfun<int,double>(x[j],p,k,i);
        y3[j]= onebasisfun<int,double,3>(x[j],k,i);
        std::cout << std::setprecision (19) << x[j] << "   " << y1[j] <<"   "<<y2[j]<<"    "<<y3[j]<< std::endl;
    }    
std::cout<<"-----------------------"<<i<<"--------------------"<<std::endl;
}
return 0;
}