#include <Rcpp.h>
using namespace Rcpp;

/* Calculates extinction */
// [[Rcpp::export]]
int is_extinct(double c1, double c2, double c3, double maxtime) {
  double X=100, Y=100;
  double time = 0;
  double u, total;
  
  while(time < maxtime) {
    total = c1*X + c2*X*Y + c3*Y;
    time = time + R::rexp(1/total);
    u = R::runif(0, 1);
    if(c1*X/total > u) {
      X = X + 1;
    } else if((c1*X+c2*X*Y)/total > u) {
      X = X - 1;
      Y = Y + 1;
    } else {
      Y = Y - 1;
    }
    
    if(X < 0.5) {return(1);}
    if(Y < 0.5) {return(2);}
  }
  return(3);
}