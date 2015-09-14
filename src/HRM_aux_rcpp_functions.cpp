

// includes from the plugin

#include <RcppArmadillo.h>


// user includes
#include "subimmer.h"


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// declarations
extern "C" {
SEXP immer_sampling_xi( SEXP x_, SEXP theta_, SEXP b_, SEXP a_, SEXP K_, 
	SEXP x_ind_, SEXP phi_, SEXP psi_, SEXP eps_, SEXP pid_, 
	SEXP rater_, SEXP N_) ;
}

// definition

SEXP immer_sampling_xi( SEXP x_, SEXP theta_, SEXP b_, SEXP a_, SEXP K_,
	SEXP x_ind_, SEXP phi_, SEXP psi_, SEXP eps_, SEXP pid_, 
	SEXP rater_, SEXP N_ ){
BEGIN_RCPP
  
       
     // probs_gpcm <- function( x , theta , b , a , K , x_ind = NULL )  
       
     Rcpp::NumericVector x(x_);          
     Rcpp::NumericVector theta(theta_);  
     Rcpp::NumericVector b(b_);  
     Rcpp::NumericVector a(a_);  
     int K = as<int>(K_);  
     Rcpp::NumericVector x_ind(x_ind_);  
     Rcpp::NumericVector phi(phi_);  
     Rcpp::NumericVector psi(psi_);  
     // Rcpp::NumericVector xi(xi_);  
     double eps=as<double>(eps_);  
     Rcpp::NumericVector pid(pid_);  
     Rcpp::NumericVector rater(rater_);  
     int N = as<int>(N_) ;  
     // int R = as<int>(R_) ;  
       
       
       
     // create matrix of probabilities  
     Rcpp::NumericMatrix probs(N,K+1);  
     int ND = x.size();   
       
       
     // vector for cumulated probabilities  
     Rcpp::NumericVector p3(N);  
     Rcpp::NumericVector x_kk(N);  
     // vector of phi for computation of HRM probabilities  
     Rcpp::NumericVector phi_rater(ND);  
     Rcpp::NumericVector psi_rater(ND);  
     Rcpp::NumericVector xi_kk(ND);  
       
     // int ii=0;  // item ii  
       
       
     for( int kk=0; kk < K+1 ; kk++){ // category kk  
       
     for (int nn=0;nn<N;nn++){  
     	x_kk[nn] = kk ;  
     		}  
       
       
     // probabilities GPCM  
     Rcpp::NumericVector p1 = subimmer_probs_gpcm_rcpp( x_kk , theta ,   b ,  
          a ,  K ,  x_ind ) ;  
     // probabilities HRM  
     for (int nn=0;nn<ND;nn++){  
     	phi_rater[nn] = phi[ rater[nn] - 1 ] ;  
     	psi_rater[nn] = psi[ rater[nn] - 1 ] ;  
     	xi_kk[nn] = kk ;  
     			}	  
     Rcpp::NumericVector p2 = subimmer_probs_hrm_rcpp(  x ,          
          xi_kk ,  phi_rater ,  psi_rater ,  K  ,   x_ind ) ;   
       
     for (int nn=0;nn<N;nn++){  
       p3[nn] = 0 ;  
       		}  
     		  
     for (int nn=0;nn<ND;nn++){  
        p2[nn] = log( p2[nn] + eps );     
        p3[ pid[nn] - 1 ] += p2[nn] ;     
                     }  
       
     for (int nn=0; nn<N ; nn++){ 		  
     	probs(nn,kk) = exp( p3[nn] ) * p1[nn] ;                  
                     	}  
     } // end category kk 	  
       
     // standardize probabilities  
     double t1=0;  
       
     for (int nn=0;nn<N;nn++){  
     t1 = 0;  
     for (int kk=0;kk<K+1;kk++){  
        t1 += probs(nn,kk);  
        		}  
     for (int kk=0;kk<K+1;kk++){  
        probs(nn,kk) = probs(nn,kk) / t1 ;  
        		}  
        	}  
       
     // draw a random number   
     Rcpp::NumericVector rn(N);  
     for (int nn=0;nn<N;nn++){  
        rn[nn] = Rf_runif(  0.0, 1.0 ) ;  
        			}  
        	  
     // draw xi vector   			  
     Rcpp::NumericVector xi_samp = subimmer_sample_prob_index( probs ,   rn ) ;  
       
       
     //*************************************************      
     // OUTPUT              
       
     return( wrap( xi_samp ) );  
     // return( wrap( p3 ) );  
       
       
     // return Rcpp::List::create(  
     // 	 _["rn"] = rn , _["xi_samp"] = xi_samp ,   
     // 	 _["probs"] = probs   		) ;
END_RCPP
}



// user includes


// declarations
extern "C" {
SEXP probs_gpcm_rcpp( SEXP x_, SEXP theta_, SEXP b_, SEXP a_, SEXP K_, SEXP x_ind_) ;
}

// definition

SEXP probs_gpcm_rcpp( SEXP x_, SEXP theta_, SEXP b_, SEXP a_, SEXP K_, SEXP x_ind_ ){
BEGIN_RCPP
  
       
     // probs_gpcm <- function( x , theta , b , a , K , x_ind = NULL )  
       
     Rcpp::NumericVector x(x_);          
     Rcpp::NumericVector theta(theta_);  
     Rcpp::NumericVector b(b_);  
     Rcpp::NumericVector a(a_);  
     int K = as<int>(K_);  
     Rcpp::NumericVector x_ind(x_ind_);  
       
       
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
END_RCPP
}




// declarations
extern "C" {
SEXP probs_hrm_rcpp( SEXP x_, SEXP xi_, SEXP phi_, SEXP psi_, SEXP K_, SEXP x_ind_) ;
}

// definition

SEXP probs_hrm_rcpp( SEXP x_, SEXP xi_, SEXP phi_, SEXP psi_, SEXP K_, SEXP x_ind_ ){
BEGIN_RCPP
  
       
     // # probs_hrm <- function( x , xi , phi , psi , K , x_ind = NULL ){  
       
     Rcpp::NumericVector x(x_);          
     Rcpp::NumericVector xi(xi_);  
     Rcpp::NumericVector phi(phi_);  
     Rcpp::NumericVector psi(psi_);  
     int K = as<int>(K_);  
     Rcpp::NumericVector x_ind(x_ind_);  
       
       
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
     // return( wrap( l1) );
END_RCPP
}




// declarations
extern "C" {
SEXP sample_prob_index( SEXP probs_, SEXP rn_) ;
}

// definition

SEXP sample_prob_index( SEXP probs_, SEXP rn_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix probs(probs_);          
     Rcpp::NumericVector rn(rn_);  
       
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
END_RCPP
}






// declarations
extern "C" {
SEXP probs_gpcm_testlet_rcpp( SEXP x_, SEXP theta_, SEXP u_, 
			SEXP b_,  SEXP a_, SEXP K_, SEXP x_ind_) ;
}

// definition

SEXP probs_gpcm_testlet_rcpp( SEXP x_, SEXP theta_, SEXP u_, 
			SEXP b_, SEXP a_, SEXP K_, SEXP x_ind_ ){
BEGIN_RCPP
  
       
     // probs_gpcm <- function( x , theta , b , a , K , x_ind = NULL )  
       
     Rcpp::NumericVector x(x_);          
     Rcpp::NumericVector theta(theta_);
     Rcpp::NumericVector u(u_);
     Rcpp::NumericVector b(b_);  
     Rcpp::NumericVector a(a_);  
     int K = as<int>(K_);  
     Rcpp::NumericVector x_ind(x_ind_);  
       
       
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
END_RCPP
}



