
#######################################################################################
###### Data Generation#############


#######################################################################################
#' export
#Norm_vec <- function(x) sqrt(sum(x^2))
#######################################################################################

#' Generate_random_direction
#'
#' @examples
#' #' Norm_vec(Generate_random_direction(3)) #checking the accuracy
#' @export
Generate_random_direction<-function(space_dim){x <- runif(space_dim, -1, 1); return(x/norm_vec(x))}
#######################################################################################



#'Generate_mix_of_VMF
#'
#' @examples
#' data_par=Generate_mix_of_VMF(100, K=1)
#' @export
Generate_mix_of_VMF<-function(n,K, mean_direction_list=NULL, concentration_list=NULL ){
  #set.seed(666)
  if(is.null(mean_direction_list)){
    mean_direction_list=list()
    #muv<-list();
    for(j in 1:K){ mean_direction_list[[j]]<-Generate_random_direction(space_dim = 3) }
  }
  if(is.null(concentration_list)){
    #muv<-list();
    concentration_list=list()
    for(j in 1:K){
      concentration_list[[j]]<-10+rgamma(1, shape = 10, rate=.1)
    }
  }
  #n<-10000 # Number of samples ; #K<-10 # Number of Clusters.
  mixture_prob <- runif(K); mixture_prob= mixture_prob/sum(mixture_prob)
  Z<-sample(1:K , n , replace=TRUE,prob=mixture_prob) # Generating Sample for cluster
  Y<-list() # Sample of Von Mises Distribution
  simulation_par_true<-list() # Noting down the simulated parameters
  for(i in 1:n){ Y[[i]]<-rvmf(1,mu= mean_direction_list[[Z[i]]],k=concentration_list[[Z[i]]]) }
   Simulation_data_and_parameter= list(data=Y, True_parameter=list(mean_direction_list=mean_direction_list,concentration_list=concentration_list, ClusterIndicator=Z , mixture_prob=mixture_prob))
   return(Simulation_data_and_parameter)
}



#######################################################################################


#' dVonMIsesFishar
#'
#' @export
dVonMIsesFishar<-function(x, mu,kappa , IfLog=FALSE ){
  #if(sum(abs(x))!=1){print("make sure that the data is on sphere")}
  #if(sum(abs(mu))!=1){print("make sure that the mu is on sphere")}
  #if(kappa<0){print("concentration, k, should be non negative. ")}
  p=length(mu); nu=p/2-1
  #log_density=(kappa*(x%*%mu))+(nu)*log(kappa)-( log(besselI( kappa,nu, expon.scaled = TRUE))+kappa)- (p/2)*log(2*pi)
  log_density=(kappa*(sum(x*mu)))+(nu)*log(kappa)-( log(besselI( kappa,nu, expon.scaled = TRUE))+kappa)- (p/2)*log(2*pi)



  if(IfLog){return(log_density)}
  if(!IfLog){return(exp(log_density))}

}

#######################################################################################

#' initial_cluster_generation
#'
#' @export
initial_cluster_generation<-function(Y, cluster_size){

          k<-cluster_size
          cdM<-matrix(nrow=length(Y),ncol=3)
            for(i in 1:length(Y))
          {
            vec<-as.vector(unlist(Y[[i]],use.names = "FALSE"))
            cdM[i,]<-vec

            }

        ######## Procedure for K means clusterin####
            #cd<-rdist(cdM,metric="angular")
            #clust1<-hclust(cd)
            #labels<-cutree(Spherical_Cluster,k)

        ###### Procedure for Spherical Cluestering###########
          Spherical_Cluster<-skmeans(cdM,k)
          labels<-Spherical_Cluster$cluster

        ##################################################
          mu0<-aggregate(cdM,list(labels),mean)  # calculating mean

            mu0<-mu0[,-1]  #removing constant coloumn

            mu0NC<-apply(mu0, 1,norm_vec) # calculating the norm

            mu0NM<-apply(mu0,1,function(x) x/norm_vec(x)) # Normalizing the means
            apply(mu0NM,2,norm_vec)                        ## Norm is not coming out to be perfect 1
            mu0M<-as.matrix(mu0)
            mu0Mt<-t(mu0M) #Transposing the raw mean matrix
            mu0N<-as.matrix(mu0NM) # matrix Conversion of Normalized mean vector
            kappa0<-c() #initialzing the kappa

            for(i in 1:k)
            {
              v1<-as.vector(mu0Mt[,i])
              v2<-as.vector(mu0N[,i])
                #kappa0<-(3*t(mu0M[,-1])%*%mu0N)/(1-t(mu0M[,-1])%*%mu0N) #Kappa
              q<-t(v1)%*%(v2)
              kappa0[i]<-(3*q)/(1-q) #Kappa
            }

            #kappa0<-3*data.frame(mapply(`*`,mu0[,-1],mu0N,SIMPLIFY=FALSE))/(1-data.frame(mapply(`*`,mu0[,-1],mu0N,SIMPLIFY=FALSE)))
            #pi0<-matrix(nrow=length(Y),ncol=cluster_size)

            pi0<-matrix(nrow=length(Y),ncol=k)

            for( i in 1:length(Y))
            {
              for(j in 1:k)
              {
                pi0[i,j]<-dVonMIsesFishar(Y[[i]],mu=mu0NM[,j],kappa=kappa0[j],IfLog=FALSE)

                #pi0[i,j]<-dvm(Y[[i]],mu=mu0NM[,j],kappa=kappa0[j],log=FALSE)

                #pi0[j]<-dvm(Y[[1]],kappa=kappa0[1],mu=mu0N[1,],l=FALSE)

              }

            }

            pi0N<-apply(pi0,1,function(x) x/sum(x))
             mu0L<-list()
            kappa0L<-list()
            for(i in 1:k)
            { mu0L[[i]]<-mu0N[,i]
              kappa0L[[i]]<-kappa0[i]
            }

            Pi_initial<-table(labels)/length(Y)  # check table comments for the label accuracy

            param<-list(Mu0=mu0L, Z=labels,Pi0=pi0N,Kappa0=kappa0L,Pi_initial= Pi_initial)

            #initial_val_in_list<- convert_matrix_to_list(param)
            # return(as.data.frame(mu0))  # Calculating means of clusters
            return(param)
  }





# Updating Mu and Kappa ####

#' Update_parameter_for_single_cluster
#'
#' @examples
#' data_var=Generate_mix_of_VMF(n=100, K = 1)
#' data=data_var$data
#' param=list(mu=unlist(data_var$True_parameter$mean_direction_list),
#'            kappa=unlist(data_var$True_parameter$concentration_list))
#' Update_parameter_for_single_cluster(sel_data=data,curr_param_j_th=param )

Update_parameter_for_single_cluster<-function(sel_data, curr_param_j_th , hyper=NULL ){
#browser()

  if(is.null(hyper)){hyper$beta_prior_kappa=0;hyper$alpha_prior_kappa=1}
  #datasummary
  if(length(sel_data)==0){
    beta_Kappa_post=beta_prior_kappa
    alpha_kappa_post<-alpha_prior_kappa
    curr_param_j_th$kappa<-rgamma(n = 1,shape = alpha_kappa_post,rate = beta_Kappa_post)
    post_mean_dir<-delta_j/norm_vec(delta_j); post_concentration<-norm_vec(delta_j);
    #browser()
    mu<-as.vector(rvmf(n=1, mu=curr_param_j_th$mu, k=curr_param_j_th$kappa))
    curr_param_j_th$mu=mu
    return(curr_param_j_th)
    }
  else{
    #datasummary
    #K<-length()
    beta_prior_kappa<-hyper$beta_prior_kappa;
    alpha_prior_kappa<-hyper$alpha_prior_kappa

    Y<-sel_data

    Y_SUM<-as.vector(Reduce("+", Y));nj=length(Y);Y_bar=Y_SUM/nj;  #Reduce here is working as cumsum
    #initialization

    kappa<-curr_param_j_th$kappa;
    mu<-curr_param_j_th$mu

    #McSample_kappa=replicate(MCSamplerSize,-1);McSample_mu=replicate(MCSamplerSize,0*Y_bar)
    #Augmenttion Step
    Vj<-rnbinom(1,size=nj,p=1-exp(-2*kappa))  # Data Augumentation Step for Jth Cluster

    #Z_NB=tryCatch(rnbinom(n=1, size=n,  prob=1-exp(-2*kappa)), error=function(e){browser()})
    #rgeom(n=1, prob=(1-exp(-2*kappa)))

    #### kappa: concentration sampler posterior #####################
    #beta_Kappa_post= calculate_beta_kappa_post_GEO(Z_GEO,mu,Y_SUM)

    beta_Kappa_post=(2*Vj+nj)-(mu)%*%(Y_SUM)+beta_prior_kappa

    alpha_kappa_post<-nj+alpha_prior_kappa

    kappa<-rgamma(n = 1,shape = alpha_kappa_post,rate = beta_Kappa_post)

    #### mu: mean direction sampler posterior #####################

    theta_prior_mu<-curr_param_j_th$mu
    zeta_prior_mu<-curr_param_j_th$kappa

    delta_j<-theta_prior_mu*zeta_prior_mu+Y_SUM*kappa

    post_mean_dir<-delta_j/norm_vec(delta_j); post_concentration<-norm_vec(delta_j);
    #browser()
    mu<-as.vector(rvmf(n=1, mu=post_mean_dir, k=post_concentration))
    curr_param_j_th$kappa=kappa
    curr_param_j_th$mu=mu
    return(curr_param_j_th)
  }
}






#' Updating Z (Cluster membership) posterior#####
#'
#' @export
updateZ<-function(Y, mu_all_cluster, kappa_all_cluster,pi_vec ){
 #browser()
     K=length(mu_all_cluster); n=length(Y)
  Z=NULL
  for(i in 1:n){
    log_post_prob=rep(-1, K)
    for(j in 1:K){
      log_post_prob[j]=dVonMIsesFishar(as.vector(Y[[i]]),mu=mu_all_cluster[[j]],kappa=kappa_all_cluster[[j]],IfLog=TRUE)
      #post_prob[j]=dvm(Y[[i]],kappa=kappa_all_cluster[[j]],mu=mu_all_cluster[[j]])*pi_all_cluster[j]
    }
    log_post_prob=log_post_prob-max(log_post_prob)
    post_prob=exp(log_post_prob)*pi_vec/sum(exp(log_post_prob)*pi_vec)
    Z[i]=sample(x = 1:K, size = 1,prob = post_prob )
  }
  return(Z)
}




#' Finite_mixture_of_VONF
#'
#' @examples
#' data=Generate_mix_of_VMF(n=1000, K = 4)
#' Y=data$data
#' mix_vonf_fit<-Finite_mixture_of_VONF(Y, K=4)
#' @export
Finite_mixture_of_VONF<-function(Y,K, max_iter=100, starting_val=NULL){


  MCMC_samples= vector("list", max_iter)
  if(   is.null(starting_val)   ){ starting_val=initial_cluster_generation(Y=Y, cluster_size = K)  }
  ##Starting point
  MCMC_samples[[1]]$mu_list= starting_val$mu_list
  MCMC_samples[[1]]$kappa_list= starting_val$kappa_list
  MCMC_samples[[1]]$Z= starting_val$Z
  MCMC_samples[[1]]$pi_all_cluster= starting_val$pi_all_cluster
  MCMC_samples[[1]]$pi_vec=starting_val$Pi_Initial

  Mcmc_Iteration=1
  browser()

  for(Mcmc_Iteration in 1:(max_iter-1)){

  print(Mcmc_Iteration)
    #browser()
  ############## Sample Z ###################
Z=updateZ(Y, MCMC_samples[[Mcmc_Iteration]]$mu_list,MCMC_samples[[Mcmc_Iteration]]$kappa_list, MCMC_samples[[Mcmc_Iteration]]$pi_vec)

MCMC_samples[[Mcmc_Iteration]]$Z=Z

    #### Sample Mu Kappa ###################
       mu_list_updated=list(); kappa_list_updated=list()
       for( j in 1:K){
      current_mu_j= MCMC_samples[[Mcmc_Iteration]]$mu_list[[j]]
      current_kappa_j=MCMC_samples[[Mcmc_Iteration]]$kappa_list[[j]]
      Data_for_j_th_group=subset(Y, Z==j)
      Mu_Kappa_next= Update_parameter_for_single_cluster( Data_for_j_th_group , curr_param_j_th=list(mu=current_mu_j,kappa=current_kappa_j))
      mu_list_updated[[j]]=Mu_Kappa_next$mu
      kappa_list_updated[[j]]=Mu_Kappa_next$kappa
    }
    MCMC_samples[[Mcmc_Iteration+1]]$mu_list = mu_list_updated  #
    MCMC_samples[[Mcmc_Iteration+1]]$kappa_list=kappa_list_updated

        ####################### Sample pi ######################
    eta=1     # why eta is 1?  # We are initializing eta to 1
    cluster_id_cnt = sapply(1:K, function(x) sum(Z == x))
    pi_vec = rdirichlet(1,eta+cluster_id_cnt-1)
    MCMC_samples[[Mcmc_Iteration+1]]$pi_vec=pi_vec
    print(pi_vec)

  }
  ##############################################################
  return( MCMC_samples )
}

