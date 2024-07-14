

#norm_Y_bar=.8; n=10; nu=1
#' @export
last_val<-function(x){return(x[length(x)])}

#' compute_h_u
#'
#' @examples
#' compute_h_u(orm_Y_bar=.8, n=10, nu=1)
#' @export
compute_h_u<-function(norm_Y_bar=.8, n=10, nu=1){

h<-function(x){ Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )}

h_p<-function(x){Marginal_log_post_kappa_prime(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n  )}

piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
inflex=piv[1]; mode=piv[2]


X_val<-c(inflex, mean(piv), mode, sum(piv), 2*mode)


x_i=X_val[-length(X_val)]; x_i_p_1=X_val[-1];
Z_val= x_i + ( h(x_i)-h(x_i_p_1) + h_p(x_i_p_1)*(x_i_p_1-x_i)  )/( h_p(x_i_p_1)-h_p(x_i) )
Z_val=c(inflex,Z_val )


intercepts<- h(X_val)-h_p(X_val)*X_val;
slopes=h_p(X_val)



return(data.frame(intercepts=intercepts,slopes=slopes, Z_val=Z_val ))
}







#' plot_h_and_h_u
#'
#' @examples
#' plot_h_and_h_u(norm_Y_bar=.8, n=10, nu=1,Ulim=20)
#' @export
plot_h_and_h_u<-function(norm_Y_bar=.8, n=10, nu=1,Ulim=20){
  h<-function(x){
    Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )
  }

  piv_quantities<-compute_h_u(norm_Y_bar=norm_Y_bar, n=n, nu=nu)
  #browser()
  intercepts=piv_quantities$intercepts; intercepts=c(intercepts, last_val(intercepts))
  slopes=piv_quantities$slopes; slopes=c(slopes, last_val(slopes))
  Z_val=c(piv_quantities$Z_val, Ulim)
  plot(h, 0.000001, Ulim)
 #y_0=intercepts[-length(Z_val)]+slopes[-length(Z_val)]*Z_val[-length(Z_val)]
# y_1=intercepts[-1]+slopes[-1]*Z_val[-1]
  k=length(intercepts)
  y_0=intercepts[-k]+slopes[-k]*Z_val[-k]
  y_1=intercepts[-1]+slopes[-1]*Z_val[-1]

 #segments(x0=Z_val[-length(Z_val)], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue')
  segments(x0=Z_val[-k], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue')

}




##############################################################################################################



#' compute_h_u
#'
#' @examples
#' compute_h_u(orm_Y_bar=.8, n=10, nu=1)
#' @export
add_tangent_lines<-function(norm_Y_bar=.8, n=10, nu=1, X_val=NULL, left_point=NULL, right_point=NULL, inflex=NULL){

        h<-function(x){ Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )}
        h_p<-function(x){Marginal_log_post_kappa_prime(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n  )}


        #if(is.infinite(right_point)){right_point=NULL}


        if(is.null(X_val)){
          piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
          inflex=piv[1]; mode=piv[2]
          X_val<-c( mean(piv), mode, sum(piv), 2*mode)
          left_point_x=inflex
        }
        if(is.null(inflex)) {piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n);  inflex=piv[1]; }
        if(is.null(left_point)){left_point_x=inflex}
        if(is.null(right_point)){right_point_x=NULL;right_point_z=Inf}
        if(!is.null(left_point)){left_point_x=left_point$x}
        if(!is.null(right_point)){right_point_x=right_point$x;right_point_z=right_point$end_line}
        X_val=unique(c(left_point_x, X_val, right_point_x))
        h_x=h(X_val); h_p_x=h_p(X_val)

        x_i=X_val[-length(X_val)]; x_i_p_1=X_val[-1];
        h_x_i=h_x[-length(X_val)]; h_x_i_p_1=h_x[-1]
        h_p_x_i=h_p_x[-length(X_val)]; h_p_x_i_p_1=h_p_x[-1]

        Z_val= x_i + ( h_x_i-h_x_i_p_1 + h_p_x_i_p_1*(x_i_p_1-x_i)  )/( h_p_x_i_p_1-h_p_x_i )
        #Z_val= x_i + ( h(x_i)-h(x_i_p_1) + h_p(x_i_p_1)*(x_i_p_1-x_i)  )/( h_p(x_i_p_1)-h_p(x_i) )
#browser()
        if(is.null(left_point)){
            Z_val=c(inflex,Z_val )
        }
        if(!is.null(left_point)){
          Z_val=c(left_point$start_line,Z_val )
        }
        intercepts<- h_x-h_p_x*X_val;
        slopes=h_p_x
        #intercepts<- h(X_val)-h_p(X_val)*X_val;
        #slopes=h_p(X_val)

      start_line=Z_val
      end_line=c(Z_val[-1],right_point_z)
      upper_bound=data.frame( intercepts=intercepts,slopes=slopes, start_line=start_line,end_line=end_line, x=X_val  )
  return(upper_bound)
}



#' @export
plot_h_and_h_u_seg<-function(norm_Y_bar=.8, n=10, nu=1,Ulim=20, upper_bound,left_most_point=1e-16){
  h_u=upper_bound
  h<-function(x){
    Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )
  }

  #piv_quantities<-compute_h_u(norm_Y_bar=norm_Y_bar, n=n, nu=nu)
  #browser()
  piv_quantities=h_u
  intercepts=piv_quantities$intercepts; intercepts=c(intercepts, last_val(intercepts))
  slopes=piv_quantities$slopes; slopes=c(slopes, last_val(slopes))

    k=length(intercepts)
    #if(Z_val[k-1]>=Ulim){Ulim=Z_val[k-1]+2}
    if(last_val(piv_quantities$start_line)>=Ulim){Ulim=last_val(piv_quantities$start_line)+2}
  Z_val=c(piv_quantities$start_line, Ulim)

  #y_0=intercepts[-length(Z_val)]+slopes[-length(Z_val)]*Z_val[-length(Z_val)]
  # y_1=intercepts[-1]+slopes[-1]*Z_val[-1]

  y_0=intercepts[-k]+slopes[-k]*Z_val[-k]
  y_1=intercepts[-1]+slopes[-1]*Z_val[-1]
  plot(h, left_most_point, Ulim)
  #segments(x0=Z_val[-length(Z_val)], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue')
  segments(x0=Z_val[-k], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue', lwd = .5)

}





append_array<-function(base_arr, ind, arr_to_add ){
  base_len=nrow(base_arr)
  add_len=nrow(arr_to_add)

  if(ind==1){
    new_arr=rbind( arr_to_add, base_arr[(ind+2):base_len,])
  }

  if((ind>1)*(ind<base_len)){
    new_arr=rbind(base_arr[1:(ind-1),], arr_to_add, base_arr[(ind+2):base_len,])
  }

  if(ind==base_len){
    new_arr=rbind(base_arr[1:(base_len-1),],arr_to_add )
  }
  return(new_arr)
}



##################################################################
##################################################################


#' compute_h_u
#'
#' @examples
#' compute_h_u(orm_Y_bar=.8, n=10, nu=1)
#' @export
add_cord_lines<-function(norm_Y_bar=.8, n=10, nu=1, X_val=NULL, left_point=NULL, right_point=NULL, inflex=NULL, left_most_x=1e-16){

  h<-function(x){ Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )}
  #h_p<-function(x){Marginal_log_post_kappa_prime(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n  )}


  #if(is.infinite(right_point)){right_point=NULL}


  if(is.null(X_val)){
    piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
    inflex=piv[1]; mode=piv[2]
    X_val<-c(inflex/4, inflex/2, inflex*3/4 )
    right_point_x=inflex;
    left_point_x=left_most_x # we want to replace it by 0
  }
  if(is.null(inflex)) {piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n);  inflex=piv[1]; }
  if(is.null(left_point)){left_point_x=left_most_x} # we want to replace it by 0
  if(is.null(right_point)){right_point_x=inflex; right_point_z=inflex}
  if(!is.null(left_point)){left_point_x=left_point$x}
  if(!is.null(right_point)){right_point_x=right_point$x;right_point_z=right_point$end_line}
  X_val=unique(c(left_point_x, X_val, right_point_x))


 # Z_val= x_i + ( h(x_i)-h(x_i_p_1) + h_p(x_i_p_1)*(x_i_p_1-x_i)  )/( h_p(x_i_p_1)-h_p(x_i) )
  #browser()
  #Z_val=X_val

  # if(is.null(right_point)){
  #   Z_val=c(Z_val,inflex )
  # }
  # if(!is.null(right_point)){
  #   Z_val=c(left_point$start_line,Z_val )
  # }
  #browser()
  h_x=h(X_val)
  x_i      =  X_val[-length(X_val)];    x_i_p_1 = X_val[-1];
  h_x_i=  h_x[-length(X_val)];      h_x_i_p_1   = h_x[-1]
  slopes= ( h_x_i_p_1-h_x_i) /(x_i_p_1-x_i )
  #slopes1=(h(x_i_p_1)-h(x_i))/( (x_i_p_1-x_i ))
  intercepts<-   h_x_i-slopes*x_i;


  start_line=x_i
  end_line=x_i_p_1
  #browser()
  upper_bound=data.frame( intercepts=intercepts,slopes=slopes, start_line=start_line,end_line=end_line, x=X_val[-length(X_val)]  )
  return(upper_bound)
}




#' @export
sample_from_proposal<-function(n=1, proposal){
#proposal=lines2

Uper=proposal$end_line
Lwer=proposal$start_line
m=proposal$slopes
inter<-proposal$intercepts
#browser()
common_factor=max(Lwer*m+inter, Uper*m+inter)


Prob_size=(exp(Uper*m+inter-common_factor)-exp(Lwer*m+inter-common_factor))/m
if(any(m==0)){
  M_0_indicator=which(m==0)
  Prob_size[M_0_indicator]=(Uper[M_0_indicator]-Lwer[M_0_indicator])*exp(inter[M_0_indicator]-common_factor)
}

total_size=sum(Prob_size)
prob=Prob_size/total_size
cum_prob<-cumsum(prob)

 u=runif(n);

  piece_selector<-function(u){1+(sum(u >cum_prob)) }

  piece_to_gen_rv<-apply(matrix(u, ncol=1), MARGIN = 1, FUN = piece_selector)
  ss=piece_to_gen_rv
  cum_sum_0<-c(0,cum_prob)
  var_1=-(inter[ss]-common_factor)+ log( (u-cum_sum_0[ss]) * (m[ss])*total_size + exp(  m[ss]*Lwer[ss]+inter[ss]- common_factor  ) )
  gen_x=var_1/m[ss]

  if(any(m[ss]==0)){
    sel_ind<-which(m[ss]==0)
    #gen_x[sel_ind]=runif(Lwer[ss], )
    gen_x[sel_ind]=(u[sel_ind]-cum_sum_0[ss[sel_ind]]) * total_size *exp(  -inter[ss[sel_ind]]+ common_factor  ) + Lwer[ss[sel_ind]]
    if(any(gen_x[sel_ind]>Uper[ss[sel_ind]])){print("Generation from the proposal distribution error. Inside function 'sample_from_proposal'")}
    }


return(gen_x)
}




##################################################

#' @export
where_to_append<-function(xval,xnew ){

  k_0=length(xval);k_1=length(xnew);k=k_0+k_1
  indicator=c( 1:k_0,replicate(k_1, -1)); row_index=1:(k_0+k_1);

  x_tot=c(xval, xnew)
  ar1=cbind(x_tot,indicator ); ar=cbind(ar1[order(x_tot),],row_index)

  #start=data.frame(ar[which(((ar[1:(k-1),2]>0)*(ar[2:k,2])<0)),], ncol=3)
  #end=data.frame(ar[which(((ar[1:(k-1),2]<0)*(ar[2:k,2])>0))+1,], ncol=3)

  start=  matrix(ar[which(((ar[1:(k-1),2]>0)*(ar[2:k,2])<0)),], ncol=3)
  end=matrix(ar[which(((ar[1:(k-1),2]<0)*(ar[2:k,2])>0))+1,], ncol=3)

  if(nrow(start)>nrow(end)){end=rbind(end, c(Inf, NA, k+1))}

  # if(nrow(start)<nrow(end)){
  #   start=rbind( c(Inf, NA, k+1))
  # }

  val_df<-data.frame(x_start     =start[, 1],
                     old_indi_start=start[,2],
                     new_indi_start= start[,3],
                     x_end         =end[,1],
                     old_indi_end=  end[,2],
                     new_indi_end_pls_1= end[,3]
                     )
  return(list(val_df=val_df,all_xval=ar[, 1] ))
}


#######################################


#' @export
Update_tangent_lines<-function(norm_Y_bar=.8, n=10, nu=1,old_lines=NULL, xnew=xnew){
  xval=old_lines$x;
  lst<-where_to_append(xval = xval,xnew=xnew )
  df_append=lst$val_df; x_all=lst$all_xval
  #browser()
  num_segments<-length(df_append$x_start)
  new_lines=old_lines[1:df_append$old_indi_start[1], ]
  for( ii in 1:num_segments){


    new_line_last_index=length(new_lines$intercepts)
    l_point <-new_lines[new_line_last_index,]
    new_lines=new_lines[ -new_line_last_index,]


    if(is.na(df_append$old_indi_end[ii])){  r_point<-NULL}
    if(!is.na(df_append$old_indi_end[ii])){r_point= old_lines[df_append$old_indi_end[ii],]}
    x_add=x_all[( df_append$new_indi_start[ii]+1):(df_append$new_indi_end_pls_1[ii]-1)]

    lines1=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu,X_val =x_add, left_point=l_point, right_point =  r_point, inflex = NULL  )
    new_lines=rbind(new_lines,lines1 )

    if(ii<num_segments){
      if(df_append$old_indi_end[ii]+1<=df_append$old_indi_start[ii+1]){
        append_old_unchanged_lines<-old_lines[(df_append$old_indi_end[ii]+1):df_append$old_indi_start[ii+1], ]
        new_lines=rbind(new_lines,append_old_unchanged_lines )
      }
    }
    if(ii==num_segments){
      end_old_lines_index=length(old_lines$x)
      if(!is.na(df_append$old_indi_end[ii])){
        if(df_append$old_indi_end[ii]<end_old_lines_index){
          new_lines=rbind(new_lines, old_lines[(df_append$old_indi_end[ii]+1):end_old_lines_index ,] )
        }
      }
    }
  }

  return(new_lines)
}


###############################################################

#' Update_tangent_lines_recalculate_all  returns the same output as Update_tangent_lines
#' Update_tangent_lines_recalculate_all calculates the lines all over again
#' whereas Update_tangent_lines adds the new one keeping the old lines
#' @export
Update_tangent_lines_recalculate_all<-function(norm_Y_bar=.8, n=10, nu=1, old_lines=NULL, xnew=xnew){
  x_new=sort(c(old_lines$x,xnew ))
  lines1=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =x_new, left_point=NULL, right_point =  NULL  )
  return(lines1)
}


##################################################################
#'
#' @export

Initiate_tangent_cord_generation<-function(norm_Y_bar=.8, n=10, nu=1, leftmost_point=1e-16, initial_num_of_tangents=20){
  #browser()
  piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
  inflex=piv[1]; mode=piv[2]
  X_val_R<-  seq(from=inflex, to =3*mode, length.out = initial_num_of_tangents)   #c( mean(piv), mode, sum(piv), 2*mode)
  X_val_L<-  seq(from=leftmost_point, to =inflex, length.out = 4)

  lines_R=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
  lines_L=add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
  lines=rbind(lines_L,lines_R )
  return(lines)
}

##################################################################
#'
#' @export

Initiate_tangent_cord_generation_V1<-function(norm_Y_bar=.8, n=10, nu=1, leftmost_point=1e-16, initial_num_of_tangents=20, xUpper=400, eps=.01){
  #browser()
  piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
  inflex=piv[1]; mode=piv[2]
  #browser()
  #X_val_R<-  seq(from=inflex, to =3*mode, length.out = initial_num_of_tangents)   #c( mean(piv), mode, sum(piv), 2*mode)
  if(inflex>0){
        X_val_R<-Compute_the_initial_grid_around_mode(norm_Y_bar=norm_Y_bar, n=n, nu=nu, eps=eps, xUpper=xUpper, piv_val=piv,   initial_num_of_points=initial_num_of_tangents,curv_a=2)
        X_val_L<-  seq(from=leftmost_point, to =inflex, length.out = 4)

        lines_R=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
        lines_L=add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
        lines=rbind(lines_L,lines_R )
  }

  if(inflex==0){
    X_val_R<-Compute_the_initial_grid_around_mode(norm_Y_bar=norm_Y_bar, n=n, nu=nu, eps=eps, xUpper=xUpper, piv_val=piv,   initial_num_of_points=initial_num_of_tangents,curv_a=2)

    #X_val_L<-  seq(from=leftmost_point, to =inflex, length.out = 4)

    lines_R=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
    #lines_L=LB_add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
    lines=lines_R
  }
  return(lines)
}







##################################################################

#' @export
Update_tangent_and_cord_lines<-function(norm_Y_bar=.8, n=10, nu=1,old_lines=NULL, xnew=xnew, inflex=NULL){

  if(is.null(inflex)){
    piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
    inflex=piv[1]; mode=piv[2]}
  xval=old_lines$x;
  #browser()
  lst<-where_to_append(xval = xval,xnew=xnew )
  df_append=lst$val_df; x_all=lst$all_xval
  #browser()
  num_segments<-length(df_append$x_start)
  new_lines=old_lines[1:df_append$old_indi_start[1], ]
  for( ii in 1:num_segments){


    new_line_last_index=length(new_lines$intercepts)
    l_point <-new_lines[new_line_last_index,]
    new_lines=new_lines[ -new_line_last_index,]


    if(is.na(df_append$old_indi_end[ii])){  r_point<-NULL}
    if(!is.na(df_append$old_indi_end[ii])){r_point= old_lines[df_append$old_indi_end[ii],]}
    x_add=x_all[( df_append$new_indi_start[ii]+1):(df_append$new_indi_end_pls_1[ii]-1)]
    if(max(x_add)<=inflex){
      lines1=add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu,X_val =x_add, left_point=l_point, right_point =  r_point, inflex = NULL  )

      new_lines=rbind(new_lines,lines1 )

      if(ii<num_segments){
        if(df_append$old_indi_end[ii]+1<=df_append$old_indi_start[ii+1]){
          append_old_unchanged_lines<-old_lines[(df_append$old_indi_end[ii]):df_append$old_indi_start[ii+1], ]
          new_lines=rbind(new_lines,append_old_unchanged_lines )
        }
      }
      if(ii==num_segments){
        end_old_lines_index=length(old_lines$x)
        if(df_append$old_indi_end[ii]<end_old_lines_index){
          new_lines=rbind(new_lines, old_lines[(df_append$old_indi_end[ii]):end_old_lines_index ,] )
        }
      }


    }
    if(min(x_add)>=inflex){
      lines1=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu,X_val =x_add, left_point=l_point, right_point =  r_point, inflex = NULL  )

      new_lines=rbind(new_lines,lines1 )

      if(ii<num_segments){
        if(df_append$old_indi_end[ii]+1<=df_append$old_indi_start[ii+1]){
          append_old_unchanged_lines<-old_lines[(df_append$old_indi_end[ii]+1):df_append$old_indi_start[ii+1], ]
          new_lines=rbind(new_lines,append_old_unchanged_lines )
        }
      }
      if(ii==num_segments){
        end_old_lines_index=length(old_lines$x)
        if(!is.na(df_append$old_indi_end[ii])){
          if(df_append$old_indi_end[ii]<end_old_lines_index){
            new_lines=rbind(new_lines, old_lines[(df_append$old_indi_end[ii]+1):end_old_lines_index ,] )
          }
        }
      }


    }

  }

  return(new_lines)
}


###############################################################

#' Update_tangent_lines_recalculate_all  returns the same output as Update_tangent_lines
#' Update_tangent_lines_recalculate_all calculates the lines all over again
#' whereas Update_tangent_lines adds the new one keeping the old lines
#' @export
Update_tangent_and_cord_lines_recalculate_all<-function(norm_Y_bar=.8, n=10, nu=1, old_lines=NULL, xnew=NULL, inflex=NULL){
  if(is.null(inflex)){
    piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
    inflex=piv[1]; mode=piv[2]
  }
#browser()
  x_new=sort(c(old_lines$x,xnew ))
  X_val_L= x_new[x_new<=inflex]
  X_val_R= x_new[x_new>inflex]

  lines_R=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
  lines_L=add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
  lines=rbind(lines_L,lines_R )
  return(lines)
}





#' Update_tangent_lines_recalculate_all  returns the same output as Update_tangent_lines
#' Update_tangent_lines_recalculate_all calculates the lines all over again
#' whereas Update_tangent_lines adds the new one keeping the old lines
#' @export
Update_tangent_and_cord_lines_recalculate_all_V1<-function(norm_Y_bar=.8, n=10, nu=1, old_lines=NULL, xnew=xnew, inflex=NULL){
  if(is.null(inflex)){
    piv<-find_inflex(nu = nu, norm_Y_bar = norm_Y_bar, n = n)
    inflex=piv[1]; mode=piv[2]
  }
  #browser()
  x_new=sort(c(old_lines$x,xnew ))
  X_val_R= x_new[x_new>inflex]

  lines_R=add_tangent_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_R,  inflex = inflex)
  lines=(lines_R )
  if(inflex>0){
  X_val_L= x_new[x_new<=inflex]
  lines_L=add_cord_lines(norm_Y_bar=norm_Y_bar, n=n, nu=nu, X_val =X_val_L,  inflex = inflex)
  lines=rbind(lines_L,lines_R )
  }

  return(lines)
}











gen_sample<-function(N=10, norm_Y_bar = .7,nu = 1,n = 100, ifPlot=TRUE, initial_num_of_tangents=4, method="UpdateALL" ){
  #browser()
  lines=Initiate_tangent_cord_generation(initial_num_of_tangents =initial_num_of_tangents , n = n, norm_Y_bar =norm_Y_bar, nu = nu )

  if(ifPlot){
    plot_h_and_h_u_seg(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound =lines, Ulim = 20)
    Sys.sleep(1)
  }


  x_prop=sample_from_proposal(n = N, proposal =lines )
  Num<-Marginal_post_kappa(kappa=x_prop, nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
  Line_ind<-  apply( matrix(x_prop, ncol=1), MARGIN=1,FUN=function(xx){ sum(lines$start_line<=xx)})
  log_prop_kernel<- (lines$slopes[Line_ind])*x_prop  + lines$intercepts[Line_ind]
  u=runif(n = N);
  accepted_ind<-(log(u)<=Num-log_prop_kernel)
  x_sample=x_prop[accepted_ind]
  size=N-length(x_sample)


  while(size>0){
   # print(size)
    add_nodes<- x_prop[((Num-log_prop_kernel)<log(.95))]
   # lines=Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL  )
    if(length(add_nodes)>0){
       if(method=='UpdateALL'){
            lines=Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL  )
       }
      if(method!='UpdateALL'){
          lines=Update_tangent_and_cord_lines(old_lines = lines,xnew = add_nodes, norm_Y_bar =norm_Y_bar, n = n, nu = nu, inflex=NULL  )
      }
    }
    if(ifPlot){
      plot_h_and_h_u_seg(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound =lines, Ulim = 20)
      #plot_h_and_bounds(norm_Y_bar =norm_Y_bar ,nu = nu,n = n,upper_bound =lines, Ulim = 20, lower_boun)
      Sys.sleep(1)
    }

    # need to change inflex in the above line.
    x_prop=sample_from_proposal(n = size, proposal =lines )
    Num<-Marginal_post_kappa(kappa=x_prop, nu = nu, norm_Y_bar =norm_Y_bar, n = n,ifLog = TRUE )
    Line_ind<-  apply( matrix(x_prop, ncol=1), MARGIN=1,FUN=function(xx){ sum(lines$start_line<=xx)})
    log_prop_kernel<- (lines$slopes[Line_ind])*x_prop  + lines$intercepts[Line_ind]
    u=runif(n = size);
    accepted_ind<-(log(u)<=Num-log_prop_kernel)
    x_sample=c( x_sample, x_prop[accepted_ind])
    size=N-length(x_sample)
  }

  if(ifPlot){
    Sys.sleep(1)
    plot(density(x_sample), xlab =c(.00000001, max(x_sample)) )
  }
  return(x_sample)

}

