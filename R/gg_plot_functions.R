

#' @export
plot_h_and_bounds<-function(norm_Y_bar=.8, n=10, nu=1,Ulim=20, upper_bound, lower_bound,left_most_point=1e-16, PlotLb=TRUE, PlotUb=TRUE,markPoints=TRUE ){
  h_u=upper_bound
  h<-function(x){
    Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )
  }
  #browser()
  #+ geom_area(stat = "function", fun = f_unscaled, fill = "green",color='black', xlim = c(0.001, 3), alpha=.4,args = list(alpha=alpha, beta=beta, gamma=gamma))


  #piv_quantities<-compute_h_u(norm_Y_bar=norm_Y_bar, n=n, nu=nu)
  #browser()
  piv_quantities=h_u
  intercepts=piv_quantities$intercepts;
  slopes=piv_quantities$slopes;

  # k=length(intercepts)
  #if(Z_val[k-1]>=Ulim){Ulim=Z_val[k-1]+2}
  if(last_val(piv_quantities$x)>=Ulim){Ulim=last_val(piv_quantities$x)+.1}
  #Z_val=(piv_quantities$start_line)

  #y_0=intercepts[-length(Z_val)]+slopes[-length(Z_val)]*Z_val[-length(Z_val)]
  # y_1=intercepts[-1]+slopes[-1]*Z_val[-1]
  ub<-function(xx ){
    Line_ind<-  apply( matrix(xx, ncol=1), MARGIN=1,FUN=function(xy){ sum(piv_quantities$start_line<=xy)})
    #browser()
    val=piv_quantities$intercepts[Line_ind]+ piv_quantities$slopes[Line_ind]* xx
    return(val)
  }
  #y=c( 0, intercepts+slopes*Z_val, 0)
  #Z_val=c(0,Z_val, last_val(Z_val) )
  #browser()
  y=c( intercepts+slopes*upper_bound$x)
  x=(upper_bound$x )
  #browser()
  lb<-function(xx ){
    Line_ind<-  apply( matrix(xx, ncol=1), MARGIN=1,FUN=function(xy){ sum(lower_bound$start_line<=xy)})

    val=lower_bound$intercepts[Line_ind]+ lower_bound$slopes[Line_ind]* xx
    val[is.na(val)]=0
    return(val)
  }

  #df=data.frame(x=Z_val, y=y)
  #p<-ggplot(df, aes( x=x, y=y )) +
  #geom_polygon(alpha=.5, col="black", fill="green", size=.1)+
  #geom_area(stat = "function", fun = h, fill = "yellow",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)
  #browser()
  p<-ggplot()+xlim(left_most_point, Ulim)
  if(PlotUb){
    #p<-p+  geom_area(stat = "function", fun = ub, fill = "orange",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)
    p<-p+stat_function(fun = ub, geom = "area",fill="orange",  color='black',  alpha=.5, size=.2, n=10000)
  }
  #p<-p+ geom_area(stat = "function", fun = h, fill = "brown",color='white', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)
  p<-p+stat_function(fun = h, geom = "area",fill="brown",  color='white',  alpha=.5, size=.2, n=10000)
  if(PlotLb){
    p<-p+stat_function(fun = lb, geom = "area",fill="white",  color='black',  alpha=.5, size=.2, n=10000)
    # p<-p+geom_area(stat = "function", fun = lb, fill = "white",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)

  }
  # geom_area(stat = "function", fun = ub, fill = "white",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.5)
  p=p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p=p+xlab("Support of Distribution")+ylab("Log Density")
  if(markPoints){
    p<-p+geom_point(size=5, aes(x =x, y = y),alpha=.6,  color="#9933FF")
    p<-p+geom_point(size=2, aes(x = x, y = y),alpha=.95,  color="#660066")
  }

  #segments(x0=Z_val[-length(Z_val)], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue')
  #segments(x0=Z_val[-k], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue', lwd = .5)
  return(p)
}


#' @export
plot_h_and_upper<-function(norm_Y_bar=.8, n=10, nu=1,Ulim=20, upper_bound, lower_bound,left_most_point=1e-16, PlotLb=TRUE, PlotUb=TRUE,markPoints=TRUE ){
  h_u=upper_bound
  h<-function(x){
    Marginal_post_kappa(kappa = x, nu = nu, norm_Y_bar =norm_Y_bar, n = n , ifLog = TRUE )
  }
  #browser()
  #+ geom_area(stat = "function", fun = f_unscaled, fill = "green",color='black', xlim = c(0.001, 3), alpha=.4,args = list(alpha=alpha, beta=beta, gamma=gamma))


  #piv_quantities<-compute_h_u(norm_Y_bar=norm_Y_bar, n=n, nu=nu)
  #browser()
  piv_quantities=h_u
 intercepts=piv_quantities$intercepts;
  slopes=piv_quantities$slopes;

 # k=length(intercepts)
  #if(Z_val[k-1]>=Ulim){Ulim=Z_val[k-1]+2}
 if(last_val(piv_quantities$x)>=Ulim){Ulim=last_val(piv_quantities$x)+.1}
  #Z_val=(piv_quantities$start_line)

  #y_0=intercepts[-length(Z_val)]+slopes[-length(Z_val)]*Z_val[-length(Z_val)]
  # y_1=intercepts[-1]+slopes[-1]*Z_val[-1]
  ub<-function(xx ){
    Line_ind<-  apply( matrix(xx, ncol=1), MARGIN=1,FUN=function(xy){ sum(piv_quantities$start_line<=xy)})
    #browser()
    val=piv_quantities$intercepts[Line_ind]+ piv_quantities$slopes[Line_ind]* xx
    return(val)
  }
  #y=c( 0, intercepts+slopes*Z_val, 0)
  #Z_val=c(0,Z_val, last_val(Z_val) )
  #browser()
  y=c( intercepts+slopes*upper_bound$x)
  x=(upper_bound$x )
  #browser()
  lb<-function(xx ){
    Line_ind<-  apply( matrix(xx, ncol=1), MARGIN=1,FUN=function(xy){ sum(lower_bound$start_line<=xy)})

    val=lower_bound$intercepts[Line_ind]+ lower_bound$slopes[Line_ind]* xx
    val[is.na(val)]=0
    return(val)
  }

  #df=data.frame(x=Z_val, y=y)
  #p<-ggplot(df, aes( x=x, y=y )) +
    #geom_polygon(alpha=.5, col="black", fill="green", size=.1)+
    #geom_area(stat = "function", fun = h, fill = "yellow",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)
   #browser()
  p<-ggplot(NULL, aes( seq(0, Ulim, by=.001) ))
  if(PlotUb){
  p<-p+  geom_area(stat = "function", fun = ub, fill = "orange",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)
  #p<-p+stat_function(fun = ub, geom = "area",fill="orange",  color='black',  alpha=.5, size=.2)
  }
  p<-p+ geom_area(stat = "function", fun = h, fill = "brown",color='white', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)
  if(PlotLb){
  p<-p+geom_area(stat = "function", fun = lb, fill = "white",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.2)

  }
   # geom_area(stat = "function", fun = ub, fill = "white",color='black', xlim = c(left_most_point, Ulim), alpha=.5, size=.5)
  p=p+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p=p+xlab("Support of Distribution")+ylab("Log Density")
  if(markPoints){
      p<-p+geom_point(size=5, aes(x =x, y = y),alpha=.6,  color="#9933FF")
      p<-p+geom_point(size=2, aes(x = x, y = y),alpha=.95,  color="#660066")
  }

  #segments(x0=Z_val[-length(Z_val)], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue')
  #segments(x0=Z_val[-k], y0=y_0, x1=Z_val[-1], y1=y_1, col='blue', lwd = .5)
  return(p)
}




