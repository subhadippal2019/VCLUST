

#' @author Dr. Subhadip Pal, UAE University
#' @export
LB_add_tangent_lines<-function(norm_Y_bar=.8, n=10, nu=1, X_val=NULL, left_point=NULL, right_point=NULL, inflex=NULL, left_most_point=1e-16){

  h<-function(x){ Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )}
  h_p<-function(x){Marginal_log_post_kappa_prime(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n  )}


  #if(is.infinite(right_point)){right_point=NULL}
  #left_eps

  if(is.null(X_val)){
    piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
    inflex=piv[1]; mode=piv[2]
    X_val<-seq( left_most_point, inflex, length.out=4)
    left_point_x=inflex
  }
  if(is.null(inflex)) {piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n);  inflex=piv[1]; }
  if(is.null(left_point)){left_point_x=left_most_point}
  if(is.null(right_point)){right_point_x=inflex;right_point_z=inflex}
  if(!is.null(left_point)){left_point_x=left_point$x}
  if(!is.null(right_point)){right_point_x=right_point$x;right_point_z=right_point$end_line}
  X_val=unique(c(left_point_x, X_val, right_point_x))
  h_x=h(X_val)
  h_p_x=h_p(X_val)
  x_i=X_val[-length(X_val)]; x_i_p_1=X_val[-1];

  h_x_i=h_x[-length(X_val)]; h_x_i_p_1=h_x[-1]
  h_p_x_i=h_p_x[-length(X_val)]; h_p_x_i_p_1=h_p_x[-1]

  #Z_val= x_i + ( h(x_i)-h(x_i_p_1) + h_p(x_i_p_1)*(x_i_p_1-x_i)  )/( h_p(x_i_p_1)-h_p(x_i) )
  Z_val= x_i + ( h_x_i-h_x_i_p_1 + h_p_x_i_p_1*(x_i_p_1-x_i)  )/( h_p_x_i_p_1-h_p_x_i )
  #browser()
  if(is.null(left_point)){
    Z_val=c(left_point_x,Z_val )
  }
  if(!is.null(left_point)){
    Z_val=c(left_point$start_line,Z_val )
  }

  intercepts<- h_x-h_p_x*X_val;
  slopes=h_p_x

  start_line=Z_val
  end_line=c(Z_val[-1],right_point_z)
  lower_bound=data.frame( intercepts=intercepts,slopes=slopes, start_line=start_line,end_line=end_line, x=X_val  )
  return(lower_bound)
}




#'
#' @author Dr. Subhadip Pal, UAE University
#' @examples
#' LB_add_cord_lines(orm_Y_bar=.8, n=10, nu=1)
#' @export
LB_add_cord_lines<-function(norm_Y_bar=.8, n=10, nu=1, X_val=NULL, left_point=NULL, right_point=NULL, inflex=NULL){

  h<-function(x){ Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )}
  #h_p<-function(x){Marginal_log_post_kappa_prime(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n  )}
  #if(is.infinite(right_point)){right_point=NULL}
#browser()
  if(is.null(X_val)){
      piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
      inflex=piv[1]; mode=piv[2]
      #X_val<-c( mean(piv), mode, sum(piv), 2*mode)
      X_val<-Compute_the_initial_grid_around_mode(norm_Y_bar=norm_Y_bar, n=n, nu=nu, eps=.001, xUpper=400, piv_val=piv,   initial_num_of_points=5,curv_a=2)

      left_point_x=inflex
      right_point_x=NULL;
      left_point_x=inflex # we want to replace it by 0
  }

  if(is.null(inflex)) {piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n);  inflex=piv[1]; }
  if(is.null(left_point)){left_point_x=inflex} # we want to replace it by 0
  if(is.null(right_point)){right_point_x=NULL; right_point_z=Inf}
  if(!is.null(left_point)){left_point_x=left_point$x}
  if(!is.null(right_point)){right_point_x=right_point$x;right_point_z=right_point$end_line}
      X_val=unique(c(left_point_x, X_val, right_point_x))
      h_x=h(X_val)

  # Z_val= x_i + ( h(x_i)-h(x_i_p_1) + h_p(x_i_p_1)*(x_i_p_1-x_i)  )/( h_p(x_i_p_1)-h_p(x_i) )
  #browser()
  #Z_val=X_val

  # if(is.null(right_point)){
  #   Z_val=c(Z_val,inflex )
  # }
  # if(!is.null(right_point)){
  #   Z_val=c(left_point$start_line,Z_val )
  # }

  x_i   =  X_val[-length(X_val)];   x_i_p_1   = X_val[-1];
  h_x_i =  h_x[-length(X_val)];     h_x_i_p_1 = h_x[-1]

  slopes    =    ( h_x_i_p_1-h_x_i ) / (x_i_p_1-x_i )
  intercepts<-   h_x_i-slopes*x_i;


  start_line=x_i
  end_line=x_i_p_1
  ### list artificial line


  #browser()
  lower_bound=data.frame( intercepts=intercepts,slopes=slopes, start_line=start_line,end_line=end_line, x=X_val[-length(X_val)]  )
  start_last_line=X_val[length(X_val)]
  lower_bound=rbind(lower_bound, c(-Inf, -Inf,start_last_line, Inf, start_last_line ))
  return(lower_bound)
}





#'
#' @author Dr. Subhadip Pal, UAE University
#' @export
plot_h_and_h_l_seg<-function(norm_Y_bar=.8, n=10, nu=1,Ulim=20, lower_bound, left_most_point=1e-16){

  h<-function(x){
    Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )
  }

  #piv_quantities<-compute_h_u(norm_Y_bar=norm_Y_bar, n=n, nu=nu)
  #browser()
  piv_quantities=lower_bound
  if(is.na(piv_quantities$intercepts[nrow(piv_quantities)])){piv_quantities=piv_quantities[-nrow(piv_quantities),]}

  intercepts=piv_quantities$intercepts; intercepts=c(intercepts, last_val(intercepts))
  slopes=piv_quantities$slopes; slopes=c(slopes, last_val(slopes))

  k=length(intercepts)
  #if(Z_val[k-1]>=Ulim){Ulim=Z_val[k-1]+2}
  if(last_val(piv_quantities$start_line)>=Ulim){Ulim=last_val(piv_quantities$start_line)}
  Z_val=c(piv_quantities$start_line)

  #y_0=intercepts[-length(Z_val)]+slopes[-length(Z_val)]*Z_val[-length(Z_val)]
  # y_1=intercepts[-1]+slopes[-1]*Z_val[-1]

  y_0=intercepts[-k]+slopes[-k]*piv_quantities$start_line
  y_1=intercepts[-1]+slopes[-1]*piv_quantities$end_line
  plot(h, left_most_point, Ulim)
  #segments(x0=Z_val[-length(Z_val)], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue')
  segments(x0=piv_quantities$start_line, y0=y_0, x1=piv_quantities$end_line, y1=y_1, col='red', lwd = .5)

}


#'
#' @author Dr. Subhadip Pal, UAE University
#' @export
LB_Initiate_tangent_cord_generation<-function(norm_Y_bar=.8, n=10, nu=1, leftmost_point=1e-16, initial_num_of_cords=20){
  #browser()
  piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
  inflex=piv[1]; mode=piv[2]
  X_val_R<-  seq(from=inflex, to =3*mode, length.out = initial_num_of_cords)   #c( mean(piv), mode, sum(piv), 2*mode)
  X_val_L<-  seq(from=leftmost_point, to =inflex, length.out = 4)

  lines_R=LB_add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
  lines_L=LB_add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
  lines=rbind(lines_L,lines_R )
  return(lines)
}


#' @author Dr. Subhadip Pal, UAE University
#' @export
LB_Initiate_tangent_cord_generation_V1<-function(norm_Y_bar=.8, n=10, nu=1, leftmost_point=1e-16, initial_num_of_cords=20, xUpper=400, eps=.01){
  #browser()
  piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
  inflex=piv[1]; mode=piv[2]

  #X_val_R<-  seq(from=inflex, to =3*mode, length.out = initial_num_of_cords)   #c( mean(piv), mode, sum(piv), 2*mode)
  if(inflex>0){
        X_val_R<-Compute_the_initial_grid_around_mode(norm_Y_bar=norm_Y_bar, n=n, nu=nu, eps=eps, xUpper=xUpper, piv_val=piv,   initial_num_of_points=initial_num_of_cords,curv_a=2)

        X_val_L<-  seq(from=leftmost_point, to =inflex, length.out = 4)

        lines_R=LB_add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
        lines_L=LB_add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
        lines=rbind(lines_L,lines_R )
  }
  if(inflex==0){
    X_val_R<-Compute_the_initial_grid_around_mode(norm_Y_bar=norm_Y_bar, n=n, nu=nu, eps=eps, xUpper=xUpper, piv_val=piv,   initial_num_of_points=initial_num_of_cords,curv_a=2)

    #X_val_L<-  seq(from=leftmost_point, to =inflex, length.out = 4)

    lines_R=LB_add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
    #lines_L=LB_add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
    lines=lines_R
  }
  return(lines)
}



#' @author Dr. Subhadip Pal, UAE University
#' @export
LB_Update_tangent_and_cord_lines_recalculate_all<-function(norm_Y_bar=.8, n=10, nu=1, old_lines=NULL, xnew=xnew, inflex=NULL){
  if(is.null(inflex)){
    piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
    inflex=piv[1]; mode=piv[2]
  }
  #browser()
  x_new=sort(c(old_lines$x,xnew ))
  X_val_L= x_new[x_new<=inflex]
  X_val_R= x_new[x_new>inflex]

  lines_R=LB_add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
  lines_L=LB_add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
  lines=rbind(lines_L,lines_R )
  return(lines)
}


#' @author Dr. Subhadip Pal, UAE University
#' @export
plot_log_density_and_Upper_lower_bound<-function(norm_Y_bar=.8, n=10, nu=1,Ulim=20, upper_bound, lower_bound, left_most_point=1e-16){

  plot_h_and_h_u_seg(norm_Y_bar=norm_Y_bar, n=n, nu=nu,Ulim=Ulim,upper_bound = upper_bound )
  par(new=TRUE)
  plot_h_and_h_l_seg(norm_Y_bar=norm_Y_bar, n=n, nu=nu,Ulim=Ulim,lower_bound = lower_bound )
}






#'rGenSample(N=10, norm_Y_bar = .7,nu = 1,n = 100)
#' @param N : Number of sample to be Generated
#' @param norm_Y_bar : Norm of the Response Directional Vector. Note that norm_Y_bar is a number between 0 and 1.
#' @param nu : nu=d/2-1 where the directional responses belongs to S^{d-1}. For circular data d=2 and \nu= \frac{1}{2}
#' @param n : Sample size for the collected data.
#' @param ifPlot=TRUE : If we want to plot/visualize the log marginal densityand the log of the proposal bounds.
#' @param newPlot=FALSE : There are two types of plot. newPlot=TRUE provides more fancy plot but requires more time and the package ggplot2
#' @param initial_num_of_points_get_inflex=4 : The role of initial_num_of_points_ greater than the inflex point  is number of points to start creating the proposal. We may increase efficiency by choosing appropriate number. If we want to generate a large number of samples for example,  10000, we may choose this parameter to be 50 or 100 for a more time efficient sampling.
#' @param sleep_time=4 : When ifPlot=TRUE, the sleep_time is the pause amount to display the adaptive change in the corresponding display.
#' @author Dr. Subhadip Pal, UAE University
#' @examples
#' rGenSample(N=10, norm_Y_bar = .7,nu = 1,n = 100, ifPlot = FALSE, BugTest1 = FALSE)
#' rGenSample(N=10, norm_Y_bar = .7,nu = 1,n = 100, ifPlot = TRUE,newPlot=FALSE , BugTest1 = FALSE, Ulim=15, sleep_time=4)
#' rGenSample(N=10, norm_Y_bar = .7,nu = 1,n = 100, ifPlot = TRUE,newPlot=TRUE , BugTest1 = FALSE, Ulim=15, sleep_time=4)
#'
#' @export
rGenSample<-function(N=10, norm_Y_bar = .7,nu = 1,n = 100, ifPlot=TRUE, initial_num_of_points_get_inflex=4, method="UpdateALL", BugTest1=TRUE, Ulim = 20, sleep_time=4, newPlot=FALSE, Version_Init_Gen=1 ){
  #browser()
   # if(norm_Y_bar ^2*n<=1){
   #        print(paste0("The value of norm_Y_bar ^2*n=", norm_Y_bar ^2*n))
   #        print("This function is only applicable for the case when norm_Y_bar ^2*n>1 ")
   #        #return(NULL)
   #  }

  if(Version_Init_Gen==1){
      lines   =    Initiate_tangent_cord_generation(initial_num_of_tangents = initial_num_of_points_get_inflex , n = n, norm_Y_bar =norm_Y_bar, nu = nu)
      lines_l = LB_Initiate_tangent_cord_generation(initial_num_of_cords    = initial_num_of_points_get_inflex , n = n, norm_Y_bar =norm_Y_bar, nu = nu)
  }
#browser()
  if(Version_Init_Gen==2){
      lines   =    Initiate_tangent_cord_generation_V1(initial_num_of_tangents = initial_num_of_points_get_inflex , n = n, norm_Y_bar =norm_Y_bar, nu = nu, eps = .001)
      lines_l = LB_Initiate_tangent_cord_generation_V1(initial_num_of_cords    = initial_num_of_points_get_inflex , n = n, norm_Y_bar =norm_Y_bar, nu = nu, eps = .001)
  }
  if(ifPlot){
    #plot_log_density_and_Upper_lower_bound(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound = lines, lower_bound = lines_l)
    if(!newPlot){
      plot_h_and_h_u_seg(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound =lines, Ulim = Ulim)
    }
      if(newPlot){
    p=plot_h_and_bounds(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound = lines, lower_bound = lines_l, Ulim=Ulim)
    print(p)
    }
    Sys.sleep(sleep_time)
  }

  #browser()

  x_prop=sample_from_proposal(n = N, proposal =lines )
  #Num<-Marginal_post_kappa(kappa=x_prop, nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
  Line_ind<-  apply( matrix(x_prop, ncol=1), MARGIN=1,FUN=function(xx){ sum(lines$start_line<=xx)})
  Line_ind_lb<-  apply( matrix(x_prop, ncol=1), MARGIN=1,FUN=function(xx){ sum(lines_l$start_line<=xx)})

  u=runif(n = N);

  log_prop_kernel<- (lines$slopes[Line_ind])*x_prop  + lines$intercepts[Line_ind]
  log_prop_kernel_lb<- (lines_l$slopes[Line_ind_lb])*x_prop  + lines_l$intercepts[Line_ind_lb]

  accepted_ind_init<-(log(u)<=log_prop_kernel_lb-log_prop_kernel)
  x_accept_init=x_prop[accepted_ind_init]; x_rej_init=x_prop[!accepted_ind_init]

  if(BugTest1){
    Num1<-Marginal_post_kappa(kappa=x_prop, nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
    print(paste0("Number points violiting the test=",sum((log_prop_kernel-Num1)*(Num1-log_prop_kernel_lb)<0)))

  }

  Num<-Marginal_post_kappa(kappa=x_prop[!accepted_ind_init], nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
  accepted_ind_2_stage<-(log(u[!accepted_ind_init])<=Num-log_prop_kernel[!accepted_ind_init])

  x_accept_2nd_stage=x_rej_init[accepted_ind_2_stage]
  x_sample=c( x_accept_init,x_accept_2nd_stage )


  size=N-length(x_sample)


  while(size>0){
    # print(size)
    #add_nodes<- x_prop[((log_prop_kernel_lb-log_prop_kernel)<log(.9))]
    add_nodes<-x_rej_init[(Num-log_prop_kernel[!accepted_ind_init])<.9]
    # lines=Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL  )
    if(length(add_nodes)>0){
      if(method=='UpdateALL'){
        lines   = Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL  )
        lines_l = LB_Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL )
        }
      if(method!='UpdateALL'){
        lines=Update_tangent_and_cord_lines(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL  )
        lines_l = LB_Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL )

        }
    }
    if(ifPlot){
      #plot_log_density_and_Upper_lower_bound(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound = lines, lower_bound = lines_l)
      if(!newPlot){
        plot_h_and_h_u_seg(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound =lines, Ulim = Ulim)
      }
      if(newPlot){
      #plot_h_and_h_u_seg(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound =lines, Ulim = 20)
      p=plot_h_and_bounds(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound = lines, lower_bound = lines_l,Ulim=Ulim)
     print(p)
      }
     Sys.sleep(sleep_time)
    }

    # need to change inflex in the above line.
    x_prop=sample_from_proposal(n = size, proposal =lines )

    Line_ind<-  apply( matrix(x_prop, ncol=1), MARGIN=1,FUN=function(xx){ sum(lines$start_line<=xx)})
    Line_ind_lb<-  apply( matrix(x_prop, ncol=1), MARGIN=1,FUN=function(xx){ sum(lines_l$start_line<=xx)})

    log_prop_kernel<- (lines$slopes[Line_ind])*x_prop  + lines$intercepts[Line_ind]
    log_prop_kernel_lb<- (lines_l$slopes[Line_ind_lb])*x_prop  + lines_l$intercepts[Line_ind_lb]
    u=runif(n = size);
    accepted_ind_init<-(log(u)<=log_prop_kernel_lb-log_prop_kernel)
    x_accept_init=x_prop[accepted_ind_init]; x_rej_init=x_prop[!accepted_ind_init]

    if(BugTest1){
      Num1<-Marginal_post_kappa(kappa=x_prop, nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
      print(paste0("Number points violiting the test=",sum((log_prop_kernel-Num1)*(Num1-log_prop_kernel_lb)<0)))
    }


    Num<-Marginal_post_kappa(kappa=x_prop[!accepted_ind_init], nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
    accepted_ind_2_stage<-(log(u[!accepted_ind_init])<=Num-log_prop_kernel[!accepted_ind_init])

    x_accept_2nd_stage=x_rej_init[accepted_ind_2_stage]
    x_sample=c( x_sample, x_accept_init,x_accept_2nd_stage )
    size=N-length(x_sample)
  }

  if(ifPlot){
    Sys.sleep(sleep_time)
    #plot(density(x_sample), xlab =c(.00000001, max(x_sample)) )
  }
  return(x_sample)

}



