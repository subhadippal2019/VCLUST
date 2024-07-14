
fun<-function(x){    return(Marginal_post_kappa(kappa=x, nu=nu, norm_Y_bar=norm_Y_bar, n=N,ifLog=FALSE))}
dMinus_Log_Fun<-function(x){  return(-Marginal_log_post_kappa_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar, n=N))   }
d2Minus_Log_Fun<-function(x){ return(-Marginal_log_post_kappa_prime_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar, n=N))  }

#x <- rMARS(100,"exp(-(4-x^2)^2)",-Inf,Inf, c(-2.5,0,2.5),c(-2/sqrt(3),2/sqrt(3)))
N=10;norm_Y_bar=.9;nu=1

vals=find_inflex(nu=1,norm_Y_bar=.9, n=N, xUpper=20,yLim=c(-2,20))
x1 <- rMARSs(n = 10,fun=fun,dMinus_Log_Fun=dMinus_Log_Fun, d2Minus_Log_Fun=d2Minus_Log_Fun ,0.00001,Inf, c(vals[1]/2, 2*vals[2]),infp = c(vals[1], 100))



