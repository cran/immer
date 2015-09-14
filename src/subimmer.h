

// includes from the plugin
#include <RcppArmadillo.h>
// #include <Rcpp.h>


using namespace Rcpp;



    
//****************************************************************
// sample from multinomial distribution
Rcpp::NumericVector subimmer_sample_prob_index( Rcpp::NumericMatrix probs ,        
	Rcpp::NumericVector rn ){

int N = rn.size();
int K = probs.ncol();

Rcpp::NumericVector xi(N);
double t1=0;
for (int nn=0;nn<N;nn++){
	t1=0;
	for (int kk=0;kk<K;kk++){
		t1 = t1 + probs(nn,kk);
//		Rcpp::Rcout << "nn=" << nn << " kk=" << kk <<
//		      " t1=" << t1 << std::endl;			
		if ( rn[nn] < t1 ){
			xi[nn] = kk ;
			break ;
				}
		}
	}

//*************************************************    
// OUTPUT            

return( wrap( xi ) );
	}

 
//**********************************************************************
//*********** probabilities GPCM ***************************************
Rcpp::NumericVector subimmer_probs_gpcm_rcpp( Rcpp::NumericVector x ,        
    Rcpp::NumericVector theta ,  Rcpp::NumericVector b ,
    Rcpp::NumericVector a , int K , Rcpp::NumericVector x_ind
 ){

int N = x.size();
Rcpp::NumericVector l1(K+1);

//    N <- length(theta)
//    KM <- matrix( 0:K , nrow = N , ncol= K+1 , byrow=TRUE)
//    b0 <- c( 0 , b[1:K] )
Rcpp::NumericVector b0(K+1);
b0[0] = 0 ;
for (int kk=1;kk<K+1;kk++){
	b0[kk] = b[kk-1] ;
		}	

//    bM <- matrix( b0 , nrow = N , ncol= K+1 , byrow=TRUE)
//    probs <- exp( a * KM *  theta - bM )
//    probs <- probs / rowSums(probs , na.rm=TRUE)


Rcpp::NumericVector probs(N);
// Rcpp::NumericMatrix p1(N,K+1);

double t1 = 0;

for (int nn=0;nn<N;nn++){	
	if ( x_ind[nn] > 0 ){
		t1 = 0 ;
		for (int kk=0; kk<K+1;kk++){
			l1[kk] = exp( a[0] * kk * theta[nn] - b0[kk] ) ;
			t1 = t1 + l1[kk] ;
				}		
		probs[nn] = l1[ x[nn]  ] / t1 ;
			} else {
		probs[nn] = 1 ;
			}
		}
		
		
//*************************************************    
// OUTPUT            

return( wrap( probs) );
// return( wrap( l1) );

}

//***************************************************************
//****** probabilities HRM ************************************

Rcpp::NumericVector subimmer_probs_hrm_rcpp( Rcpp::NumericVector x ,        
    Rcpp::NumericVector xi ,  Rcpp::NumericVector phi ,
    Rcpp::NumericVector psi ,  int K  ,
    Rcpp::NumericVector x_ind ) {

int N = x.size();
Rcpp::NumericVector l1(K+1);
Rcpp::NumericVector probs(N);

double t1 = 0;

for (int nn=0; nn < N ; nn++){
	//	KM <- matrix( 0:K , nrow=N , ncol=K+1 , byrow=TRUE )
	//    p1 <- exp( - ( KM - ( xi + phi ) )^2 / ( 2 * psi ) )	
	//    probs <- p1 / rowSums(p1 , na.rm=TRUE)
if ( x_ind[nn] > 0 ){	
	t1 = 0 ;
	for (int kk=0;kk<K+1;kk++){
		l1[kk] = exp( - pow( kk - xi[nn] - phi[nn] , 2.0 ) / 2 / psi[nn] ) ;
		t1 += l1[kk] ;
				}
	probs[nn] = l1[ x[nn] ] / t1 ;
  } else {
 probs[nn] = 1 ;
 	}
 }
		
//*************************************************    
// OUTPUT            

return( wrap( probs) );
    }




//**********************************************************************
//*********** probabilities GPCM testlet *******************************
Rcpp::NumericVector subimmer_probs_gpcm_testlet_rcpp( Rcpp::NumericVector x ,        
    Rcpp::NumericVector theta ,  Rcpp::NumericVector u , Rcpp::NumericVector b ,
    Rcpp::NumericVector a , int K , Rcpp::NumericVector x_ind
 ){

int N = x.size();
Rcpp::NumericVector l1(K+1);

//    N <- length(theta)
//    KM <- matrix( 0:K , nrow = N , ncol= K+1 , byrow=TRUE)
//    b0 <- c( 0 , b[1:K] )
Rcpp::NumericVector b0(K+1);
b0[0] = 0 ;
for (int kk=1;kk<K+1;kk++){
	b0[kk] = b[kk-1] ;
		}	

//    bM <- matrix( b0 , nrow = N , ncol= K+1 , byrow=TRUE)
//    probs <- exp( a * KM *  theta - bM )
//    probs <- probs / rowSums(probs , na.rm=TRUE)


Rcpp::NumericVector probs(N);
// Rcpp::NumericMatrix p1(N,K+1);

     double t1 = 0;         
     for (int nn=0;nn<N;nn++){	  
     	if ( x_ind[nn] > 0 ){  
     		t1 = 0 ;  
     		for (int kk=0; kk<K+1;kk++){  
     			l1[kk] = exp( a[0] * kk * theta[nn] + kk*u[nn] - b0[kk] ) ;  
     			t1 = t1 + l1[kk] ;  
     				}		  
     		probs[nn] = l1[ x[nn]  ] / t1 ;  
     			} else {  
     		probs[nn] = 1 ;  
     			}  
     		}  
     		  
		
		
//*************************************************    
// OUTPUT            

return( wrap( probs) );
// return( wrap( l1) );

}    

