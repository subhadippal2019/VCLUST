





lines=add_tangent_lines()
k=length(lines$x)
xval=seq(from=lines$x[k], to=2*lines$x[k], length.out = 6 )[2:5]
lines1=add_tangent_lines(X_val =xval, left_point = lines[k, ], right_point  = NULL )

lines_add= lines


lines_l=LB_Initiate_tangent_cord_generation(norm_Y_bar = 0.8, n = 10, nu = 1, initial_num_of_cords = 5, leftmost_point = 0.000001)
lines_u=Initiate_tangent_cord_generation(norm_Y_bar = 0.8, n = 10, nu = 1, initial_num_of_tangents = 5)

#p=plot_h_and_upper(upper_bound = lines_u, lower_bound =lines_l )
p=plot_h_and_bounds(norm_Y_bar = 0.8, n = 10, nu = 1,Ulim = 10, upper_bound = lines_u, lower_bound =lines_l )
p

plot_log_density_and_Upper_lower_bound(norm_Y_bar = 0.8, n = 10, nu = 1, upper_bound = lines_u, lower_bound = lines_l )



plot_h_and_h_u_seg(upper_bound = lines_u )
#if(last_val(piv_quantities$start_line)>=Ulim){Ulim=last_val(piv_quantities$start_line)+2}
plot_h_and_h_u_seg(upper_bound = lines)


plot_h_and_h_u_seg(upper_bound = lines)


lines2=append_array(lines, k, lines1)
plot_h_and_h_u_seg(h_u=lines2)



library(profvis)
profvis({
  rGenSample(N=100, BugTest1=FALSE, ifPlot = FALSE, n = 1000000, initial_num_of_points_get_inflex = 200)
}
)

profvis({
  gen_sample(N=100, ifPlot = FALSE, n = 1000000, initial_num_of_tangents = 200 )
})


library(rbenchmark)

benchmark(gen_sample(N=10000, ifPlot = FALSE, n = 1000, initial_num_of_tangents = 30 ),
          rGenSample(N=10000, BugTest1=FALSE, ifPlot = FALSE, n = 1000, initial_num_of_points_get_inflex = 30),
          columns=c('test', 'elapsed', 'replications'), replications=20)


rGenSample(N=10000, BugTest1=FALSE, ifPlot = TRUE, n = 100, initial_num_of_points_get_inflex = 6, )


lines1_u=add_cord_lines(X_val = c(.01, .04))
lines2_u=add_tangent_lines(X_val = c(1, 5, 7, 11, 15))

lines=rbind(lines2, lines1)
lines_u=rbind(lines1_u, lines2_u)
plot_h_and_h_u_seg(upper_bound = lines_u, Ulim = 20)





###############
lines=add_tangent_lines(norm_Y_bar = .7,nu = 1,n = 100)
lines1=add_cord_lines(norm_Y_bar = .7,nu = 1,n = 100)
lines2=rbind(lines1, lines)

add_cord_lines()



profvis({x_sample=gen_sample(N=70000, initial_num_of_tangents = 50, norm_Y_bar = .7,n = 10,  ifPlot = FALSE, method = "aa")})

benchmark(gen_sample(N=70000, initial_num_of_tangents = 200, norm_Y_bar = .7,n = 10,  ifPlot = FALSE),
          gen_sample(N=70000, initial_num_of_tangents = 20, norm_Y_bar = .7,n = 10,  ifPlot = FALSE),
          columns=c('test', 'elapsed', 'replications'))








x_sample=gen_sample(N=1000)
plt=plot_marginl_post(norm_Y_bar = .7,nu = 1,n = 100, x_lim=c(.00000001, max(x_sample)))


lines=Initiate_tangent_cord_generation(initial_num_of_tangents = 4)
xval=lines$x
xnew=c(8, 10)
xnew=c( .10, .2)
new_lines1=Update_tangent_and_cord_lines(old_lines = lines,xnew = xnew )
new_lines=Update_tangent_lines(old_lines = lines,xnew = xnew )

new_lines2<-Update_tangent_and_cord_lines_recalculate_all(old_lines = lines,xnew = xnew)

plot_h_and_h_u_seg(upper_bound = lines)
plot_h_and_h_u_seg(upper_bound = new_lines)

new_lines11=Update_tangent_lines_recalculate_all(norm_Y_bar=.8, n=10, nu=1,old_lines = lines,xnew = xnew)





benchmark(a={rGenSample(N=10000, ifPlot = FALSE, n = 1000, initial_num_of_points_get_inflex = 30, Version_Init_Gen = 1 )},
          b={rGenSample(N=10000, BugTest1=FALSE, ifPlot = FALSE, n = 1000, initial_num_of_points_get_inflex = 30, Version_Init_Gen = 2)},
          replications=10)













Update_tangent_lines_old<-function(old_lines=NULL, xnew=xnew){
  xval=old_lines$x;
lst<-where_to_append(xval = xval,xnew=xnew )
df_append=lst$val_df; x_all=lst$all_xval
browser()

new_lines=old_lines[1:df_append$old_indi_start[1], ]
for( ii in 1:length(df_append$x_start)){
  new_line_last_index=length(new_lines$intercepts)
  l_point <-new_lines[new_line_last_index,]
  new_lines=new_lines[ -new_line_last_index,]


if(is.na(df_append$old_indi_end[ii])){  r_point<-NULL}
if(!is.na(df_append$old_indi_end[ii])){r_point= old_lines[df_append$old_indi_end[ii],]}
x_add=x_all[( df_append$new_indi_start[ii]+1):(df_append$new_indi_end_pls_1[ii]-1)]

lines1=add_tangent_lines(X_val =x_add, left_point=l_point, right_point =  r_point  )
new_lines=rbind(new_lines,lines1 )
}

return(new_lines)
}

#
# k_0=length(xval);k_1=length(xnew)
# indicator=c( 1:k_0,replicate(k_1, -1))
# row_index=1:(k_0+k_1)
# ar=ar1[order(x_tot),]
# k=nrow(ar)
# start=ar[which(((ar[1:(k-1),2]>0)*(ar[2:k,2])<0)),]
# end=ar[which(((ar[1:(k-1),2]<0)*(ar[2:k,2])>0))+1,]

