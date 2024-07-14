#a="f<-function(x){x^2};plot(f, 0, 2)"
#library(Bessel)

############################################################################################
#############################################################################################
#log_bessel_I <- function(x, nu){ return( log(BesselI(z=x,nu=nu, expon.scaled = TRUE))+x)  }


############################################################################################
#############################################################################################

 ############################################################################################
 #############################################################################################

############################################################################################
#############################################################################################

#plot(function(x){f(x,nu = .5,norm_Y_bar = .5,N = 100)}, xlim=c(0,4),col="red")


approx_f<-function(x){
  nu = .5;norm_Y_bar = .5;N = 100
  log_val=-x*N*(1-norm_Y_bar)+log(x)*(N-1)*(nu+.5)
}
# term1=-.5*log(2*pi*x)
#
# seq=
#
#
#   cum_term<-function(k, x, nu){
#     browser()
#     seq_k=( 2*(1:k)-1)^2
#     terms_seq=sum(log(4*nu^2-seq_k))-k*log(8*x)-log(gamma(k+1))
#     return(exp(terms_seq))
#
#   }
#




sample_iid<-function(sample_size=10, kappa_true=.9){


  f<-function(x, nu, N,norm_Y_bar){
    #browser()
    val=nu*(N-1)*log(x)+ log_Bessel_I(x=N*norm_Y_bar*x, nu = nu)-N*log_Bessel_I(x=x, nu=nu)
    return(val)
  }
  ############################################################################################
  #############################################################################################
  fprima<-function(x,nu, N,norm_Y_bar){

    a=N*norm_Y_bar
    val=  a*ratio_Bessel_I(x=a*x,nu=nu)-N*ratio_Bessel_I(x=x, nu=nu)

    return(val)
  }
  N=2000
theta=rvmf(n=N, mu=c(1,0,0), k=kappa_true)
theta_bar=apply(theta,2, mean);
nu=length(theta_bar)/2-1
norm_Y_bar=norm_vec(theta_bar)
print(norm_Y_bar)
########################################################################
########################################################################




#f<-function(x,mu=0,sigma=1){-1/(2*sigma^2)*(x-mu)^2}
#fprima<-function(x,mu=0,sigma=1){-1/sigma^2*(x-mu)}
N=2000
mysample<-ars(n=sample_size,f,fprima,lb = TRUE,x=c( .76, 1, 2),xlb=.75, ub = TRUE,xub=30, nu=nu, N=N, norm_Y_bar=norm_Y_bar   )
mysample
plot(density(mysample))
hist(mysample)
return(mysample)
}


text<-function(){
  library(tictoc)
  library(ars)
  tic()
x1=sample_iid(sample_size =10000, kappa_true = 10 )
toc()

tic()
x_2=gen_sample(N=100000, ifPlot = FALSE, n=2000, nu = .5, norm_Y_bar =0.8987719, initial_num_of_tangents = 50 )
toc()
tic(); x_3=gen_sample(N=100000, ifPlot = FALSE, n=2000, nu = .5, norm_Y_bar =0.8987719, method = "aa" , initial_num_of_tangents = 50); toc()
tic(); x_3=rGenSample(N=100000, ifPlot = FALSE, n=2000, nu = .5, norm_Y_bar =0.8987719, method = "aa" , initial_num_of_points_get_inflex  = 50, BugTest1 = FALSE); toc()

library(profvis)
profvis(gen_sample(N=10000, ifPlot = FALSE, n=2000, nu = .5, norm_Y_bar =0.8987719, initial_num_of_tangents = 50 ))
profvis(rGenSample(N=10000, ifPlot = FALSE, n=2000, nu = .5, norm_Y_bar =0.8987719 , initial_num_of_points_get_inflex  = 50, BugTest1 = FALSE))
profvis(rGenSample1(N=10000, ifPlot = FALSE, n=2000, nu = .5, norm_Y_bar =0.8987719 , initial_num_of_points_get_inflex  = 50, BugTest1 = FALSE))
}











