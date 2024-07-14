
#' calculates log(K(x, nu)) where K is the modified Bessel function of the Second kind.
#'
#' @param x A complex number.
#' @param nu degree of the Bessel function
#' @return The log(K(\code{x},\code{nu}))
#' @examples
#' log_Bessel_K(10,.5)
#' log_Bessel_K(10+3*1i,0)
#'  log_Bessel_K(10000,0)
#' @export
log_Bessel_K<-function(x, nu){ return( log( suppressWarnings(BesselK(x,nu =nu, expon.scaled = TRUE  )))-x ) }


#' calculates log(K(x, nu)) where K is the modified Bessel function of the Second kind.
#'
#' @param x A complex number.
#' @param nu degree of the Bessel function
#' @return The log(K(\code{x},\code{nu}))
#' @examples
#' log_Bessel_K(10,.5)
#' log_Bessel_K(10+3*1i,0)
#'  log_Bessel_K(10000,0)
#' @export
log_Bessel_I<-function(x, nu){
  val=log(suppressWarnings(BesselI(x,nu =nu, expon.scaled = TRUE  )))+x
  return( val)
  }

#' @export
ratio_Bessel_I<-function(x,nu, Iflog=FALSE){

  #val= log_bessel_I(x, nu+1)-log_bessel_I(x, nu)
  val= log_Bessel_I(x, nu+1)-log_Bessel_I(x, nu)
  if(!Iflog){ val=exp(val)}
  return(val)
}


norm_vec<-function(x){  sqrt(sum(x^2)) }
#' calculates log(sinh(x)) for small and large arguments.
#'
#' @param x A real number, a vector of real numbers.
#' @return The log(sinh(\code{x}))
#' @examples
#' log_sinh(10)
#' log_sinh(c(1000,1,10,10000))
#' @export
log_sinh<-function(x){
  log_sinh_single<-function(y){
    if(y<200){return(log(sinh(y)))}
    if(y>=200){return(y-log(2))}
  }
  val=apply(X = matrix(x, ncol=1),MARGIN = 1,FUN = log_sinh_single)
  return(val)
}




inner_prod<-function(x,y){(sum(x*y))}
#### Data Generation #####
alt_besselI<-function(x, nu){
  if(nu==0){
    val=gsl::bessel_I0(x)
  }
  if(nu>0){
    val=besselI(10,nu = nu)
  }
  return(val)
}



############################################


#' @export
#'
f1<-function(x){
  return(log_Bessel_I(x=x, nu=0))
}




#' @export
plot_bessel_I_0_approx<-function(){
              f1_L<-function(x){
                return(.5*(log(sinh(2*x))-log(pi*x) ) )
              }

              f1_U<-function(x){
                return(.5*(log_sinh(2*x)-log(2*x) ) )
              }
              #
              #
              # f1_L<-function(x){
              #   return(.5*(log_sinh(2*x)-log(pi*x) ) )
              # }

              f2_U<-function(x){
                    return(.5*(log_sinh(x)-log(x) )+.5*log(cosh(p*x))/p )
              }

              f2_L<-function(x){
                    return(.5*(log_sinh(x)-log(x) )+.5*log(cosh(q*x))/q )
              }

        p=2/3 # or greater
        q=log(2)/log(pi) # or less



        Ylim=c(0,50)

        plot(f1, 0.0000001, 50, ylim=Ylim)
        par(new=T)
        plot(f1_L, 0.0000001, 50, ylim=Ylim, col="red")
        par(new=T)
        plot(f1_U, 0.0000001, 50, ylim=Ylim, col="red")

        par(new=T)
        plot(f2_L, 0.0000001, 50, ylim=Ylim, col="blue")
        par(new=T)
        plot(f2_U, 0.0000001, 50, ylim=Ylim, col="blue")
return(NULL)
}


#########################################################################################
#########################################################################################
#########################################################################################

#' @export
Marginal_post_kappa_old<-function(kappa, nu, norm_Y_bar, n,ifLog=FALSE){
  # if(kappa==0){
  #   val=nu*( log(n)+log(norm_Y_bar) )  +(n-1)*(  nu*log(2)+lgamma(nu+1)    )
  # }
  # else if(kappa>0){
        val=log_Bessel_I(x =n*norm_Y_bar* kappa, nu = nu) -n*log_Bessel_I(x = kappa, nu = nu)+nu*(n-1)*log(kappa)
  #}

  if(!ifLog){  val=exp(val)   }
  return(val)
}


#' @export
Marginal_post_kappa_test<-function(kappa, nu, norm_Y_bar, n,ifLog=FALSE, method=2){

  if(method==1){
    val=ifelse(kappa==0, nu*( log(n)+log(norm_Y_bar) )  +(n-1)*(  nu*log(2)+lgamma(nu+1)    ), log_Bessel_I(x =n*norm_Y_bar* kappa, nu = nu) -n*log_Bessel_I(x = kappa, nu = nu)+nu*(n-1)*log(kappa) )
  }
  if(method!=1){
    val=0*kappa-1; non_zero_ind<-(kappa!=0)
              #if(sum(!non_zero_ind)){
                val[!non_zero_ind]=nu*( log(n)+log(norm_Y_bar) )  +(n-1)*(  nu*log(2)+lgamma(nu+1)    )
                non_zero_kappa=kappa[non_zero_ind]
                val[non_zero_ind]=log_Bessel_I(x =n*norm_Y_bar* non_zero_kappa, nu = nu) -n*log_Bessel_I(x = non_zero_kappa, nu = nu)+nu*(n-1)*log(non_zero_kappa)
              #}
  }
  if(!ifLog){  val=exp(val)   }
  return(val)
}


#' @export
Marginal_post_kappa<-function(kappa, nu, norm_Y_bar, n,ifLog=FALSE){

    val=kappa+NA; #benchmark({kappa+NA}, replicate(n=length(kappa), NA), replications = 10000)
    non_zero_ind<-(kappa!=0); non_zero_kappa=kappa[non_zero_ind]

    val[!non_zero_ind]=nu*( log(n)+log(norm_Y_bar) )  +(n-1)*(  nu*log(2)+lgamma(nu+1)    )
    val[non_zero_ind]=log_Bessel_I(x =n*norm_Y_bar* non_zero_kappa, nu = nu) -n*log_Bessel_I(x = non_zero_kappa, nu = nu)+nu*(n-1)*log(non_zero_kappa)

  if(!ifLog){  val=exp(val)   }
  return(val)
}

#########################################################################################
#########################################################################################
#########################################################################################




#' @export
Marginal_log_post_kappa_prime_old<-function(kappa,nu,norm_Y_bar, n){
  a=n*norm_Y_bar
  val=  a*ratio_Bessel_I(x=a*kappa,nu=nu)-n*ratio_Bessel_I(x=kappa, nu=nu)
  return(val)
}




#' Marginal_log_post_kappa_prime(kappa,nu,norm_Y_bar, n)
#'
#' @examples
#' Marginal_log_post_kappa_prime(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10)
#' @export
Marginal_log_post_kappa_prime<-function(kappa,nu,norm_Y_bar, n){
  a=n*norm_Y_bar

  val=kappa+NA; #benchmark({kappa+NA}, replicate(n=length(kappa), NA), replications = 10000)
  non_zero_ind<-(kappa!=0); non_zero_kappa=kappa[non_zero_ind]
  val[!non_zero_ind]=0

  val[non_zero_ind]=  a*ratio_Bessel_I(x=a*non_zero_kappa,nu=nu)-n*ratio_Bessel_I(x=non_zero_kappa, nu=nu)
  return(val)
}



#' @export
Marginal_log_post_kappa_prime_prime_old<-function(kappa,nu,norm_Y_bar, n){

  a=n*norm_Y_bar
  (a^2-n)/(2*(nu+1))
  ratio_bes1=ratio_Bessel_I(x=a*kappa,nu=nu)
  ratio_bes2=ratio_Bessel_I(x=kappa, nu=nu)
  val=  a^2*(1-  (2*nu+1)*ratio_bes1/(a*kappa) - (ratio_bes1)^2  ) -  n*(1-(2*nu+1)*ratio_bes2/kappa-ratio_bes2^2  )
  return(val)
}




#'Marginal_log_post_kappa_prime_prime
#'
#' @examples
#' Marginal_log_post_kappa_prime_prime(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10)
#' benchmark(Marginal_log_post_kappa_prime_prime(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10),
#' Marginal_log_post_kappa_prime_prime_old(kappa=c(0, 0, 1, 2, 3), nu = 1, norm_Y_bar = .8, n = 10),
#' replications = 10000 )
#' @export
Marginal_log_post_kappa_prime_prime<-function(kappa,nu,norm_Y_bar, n){

  a=n*norm_Y_bar

  val=kappa+NA; #benchmark({kappa+NA}, replicate(n=length(kappa), NA), replications = 10000)
  non_zero_ind<-(kappa!=0); non_zero_kappa=kappa[non_zero_ind]
  val[!non_zero_ind]=(a^2-n)/(2*(nu+1))


  ratio_bes1=ratio_Bessel_I(x=a*non_zero_kappa,nu=nu)
  ratio_bes2=ratio_Bessel_I(x=non_zero_kappa, nu=nu)
  val[non_zero_ind]=  a^2*(1-  (2*nu+1)*ratio_bes1/(a*non_zero_kappa) - (ratio_bes1)^2  ) -  n*(1-(2*nu+1)*ratio_bes2/non_zero_kappa-ratio_bes2^2  )
  return(val)
}


#' plot_marginl_post
#'
#' @examples
#' plot_marginl_post(nu=.5,norm_Y_bar = .3, x_lim = c(.00001,10), ifLog = TRUE)
#' VCLUST::plot_marginl_post(nu=.5,norm_Y_bar = .8, x_lim = c(0,10), ifLog = TRUE)
#' @export
plot_marginl_post<-function(nu=0,norm_Y_bar=.9, n=10,x_lim=c(0.0000001,5),ifLog=FALSE){
  base <- ggplot() + xlim( x_lim[1], x_lim[2])

  plt=base +
    geom_function(aes(colour = "blue"), fun = Marginal_post_kappa, args = list(nu=nu,norm_Y_bar=norm_Y_bar, n=n,ifLog=ifLog))
  return(plt)
}


#' plot_marginl_post
#'
#' @examples
#' par_list=list(nu=1,norm_Y_bar=.9, n=10,ifLog=TRUE)
#' plot_marginl_post(nu=.5,norm_Y_bar = .8, x_lim = c(0,10), ifLog = TRUE)
#' plot_g_function(fun=Marginal_post_kappa,par_list=par_list )
#' par_list1=list(nu=1,norm_Y_bar=.9, n=10)
#' plot_g_function(fun=Marginal_log_post_kappa_prime_prime,par_list=par_list1 )
#'
#' @export
plot_g_function<-function(fun, par, x_lim=c(0.0000001,5), par_list){
  #par_list=list(nu=nu,norm_Y_bar=norm_Y_bar, n=n,ifLog=ifLog)
  base <- ggplot() + xlim( x_lim[1], x_lim[2])

  plt=base +
    geom_function( fun = fun, args =par_list )
  return(plt)
}



#' @examples
#' find_inflex(nu=1,norm_Y_bar=.9, n=10, xUpper=20,yLim=c(-3,350))
#' @export
find_inflex<-function(nu=1,norm_Y_bar=.9, n=10, xUpper=400,yLim=c(-2,20), ifPlot=FALSE){
#browser()
  error_floating_point<-0#1e-16
  if( (n*norm_Y_bar*norm_Y_bar)<=(1+error_floating_point)){
   return( c(inflex=0, mode=0))
  }
  fun <- function(x) {Marginal_log_post_kappa_prime_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n)}
  fun1<-function(x) {Marginal_log_post_kappa_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n)}
 #if(fun1(xUpper)>0){xUpper=3*xUpper }  #### need to adjust here

  uni <- uniroot(fun, c(1e-10, xUpper), extendInt="yes", trace=4)$root
  uni1 <- uniroot(fun1, c(1e-10, xUpper), extendInt="yes", trace=4)$root
        if(ifPlot){
          fun2<-function(x) {Marginal_post_kappa(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n, ifLog = TRUE)}
               curve(fun(x), 0.000001, xUpper, ylim=yLim)
               abline(h = 0, lty = 3)

               points(uni, 0, pch = 16, cex = 1)
               abline(v =uni , lty = 3)
               curve(fun1(x), 0.000001, xUpper,add = TRUE)

               points(uni1, 0, pch = 16, cex = 1)
               curve(fun2(x), 0.000001, xUpper,add = TRUE)
         }
   return(c(inflex=uni, mode=uni1))
}


#' @export
find_inflex_1<-function(nu=1,norm_Y_bar=.9, n=10, xUpper=20){
    fun <- function(x) {Marginal_log_post_kappa_prime_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n)}
    xstart <- c(2,0.5)
    #browser()
    val=nleqslv(xstart, fun, control=list(btol=.01))
    curve(fun(x), 0.000001, xUpper, ylim=c(-2, 20))
    abline(h = 0, lty = 3)
    points(min(val$x), 0, pch = 16, cex = 2)
    return(val$x)
}


#' @export
Plot_inflex<-function(nu=1,norm_Y_bar=.9, n=10, xUpper=20,yLim=c(-2,20), ifPlot=FALSE){
  val=find_inflex(nu=nu,norm_Y_bar=norm_Y_bar, n=n, xUpper=xUpper,yLim=yLim, ifPlot = TRUE)
  return(val)
}




#' @export
Find_approx_range_of_marginal_post_kappa<-function(norm_Y_bar=.8, n=10, nu=1, eps=.01, xUpper=400, piv_val=NULL){
#browser()
  if(is.null(piv_val)){
      piv_val=find_inflex(norm_Y_bar=norm_Y_bar, n=n, nu=nu)
  }
  mode=piv_val[2]; inflex<-piv_val[1]
  a= norm_Y_bar*n
  log_max_val=Marginal_post_kappa(kappa=mode, norm_Y_bar=norm_Y_bar, n=n, nu=nu, ifLog = TRUE)

  #fun <- function(x) {Marginal_log_post_kappa_prime_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n,)-eps*}
  #fun1<-function(x) {Marginal_log_post_kappa_prime(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n)}
  fun2<-function(x) {Marginal_post_kappa(kappa=x,nu=nu,norm_Y_bar=norm_Y_bar , n=n, ifLog = TRUE)- log(eps)-(log_max_val)}
  #browser()
  uni_right <- uniroot(fun2, c(mode, xUpper), extendInt="yes", trace=2)$root
  if(inflex==0){uni_left = 0}
  else if(fun2(inflex) >= 0){ uni_left = inflex}
  else{
    #print("here")
    uni_left <- uniroot(fun2, c(inflex, mode),extendInt="yes", trace=2 )$root
  }
  #uni1 <- uniroot(fun1, c(1e-10, xUpper))$root
return(range_par=data.frame(inflex=inflex,  mode=mode, left_point=uni_left,  right_point=uni_right))
}




#' @export
Compute_the_initial_grid_around_mode<-function(norm_Y_bar=.8, n=10, nu=1, eps=.01, xUpper=400, piv_val=NULL,   initial_num_of_points=10,curv_a=1.5){
  range_par<-Find_approx_range_of_marginal_post_kappa(norm_Y_bar=norm_Y_bar, n=n, nu=nu, eps=eps, xUpper=xUpper, piv_val=piv_val)
  mode=range_par$mode; inflex=range_par$inflex; left_point=range_par$left_point; right_point=range_par$right_point
#browser()

    curv_a=max(1, curv_a)

      #m=10; lll=5; uuu=20; n_left=5; n_right=8; curv_a=1.5
    if(inflex>0){
      n_left<-floor(initial_num_of_points/2);n_right<-ceiling(initial_num_of_points/2)
          left_seq= mode  - (mode-left_point)*seq(1,0.1, length.out =n_left )^curv_a ;
          right_seq=  mode  + (right_point-mode)*seq(0,1, length.out =n_right )^curv_a ;
          all_seq=c(inflex, left_seq, right_seq)
    }
    if(inflex==0){
      #left_seq= mode  - (mode-left_point)*seq(1,0.1, length.out =n_left )^curv_a ;
      n_right<-(initial_num_of_points)
      right_seq=   (right_point)*seq(0,1, length.out =n_right )^curv_a ;
      all_seq= right_seq
    }
     #plot(all_seq, all_seq)
      return(all_seq)
}

#par<-Find_approx_range_of_marginal_post_kappa()
#seq_all <- Compute_the_initial_grid_around_mode(mode = par$mode, left_point = par$left_point, right_point = par$right_point )


