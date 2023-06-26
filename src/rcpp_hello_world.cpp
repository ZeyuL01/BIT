
#include <Rcpp.h>
using namespace Rcpp;

//' Name
//'
//' Description
//' @param variable Parameter_Description
//' @return Return_Description
// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}
