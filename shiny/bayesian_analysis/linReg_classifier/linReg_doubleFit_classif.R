library(rjags)
library(MASS)
library(parallel)
source("helper_functions.R", local=TRUE)

folder = "linReg_doubleFit_classifier"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

modelstring = "
 model {
  # Likelihood of data given model parameters
  for (i in 1:N){
   Yobs[i] ~ dnorm(m_hat[i]*Xobs[i]+c_hat[i], tau_hat[i])
   
   tau_hat[i] = ifelse(class[i]==1, tau2, tau1)
   m_hat[i] = ifelse(class[i]==1, m2, m1)
   c_hat[i] = ifelse(class[i]==1, c2, c1)
   
   class[i] ~ dbern(probdiff)
  }
  
  for(j in 1:Nsyn){
   Ysyn_norm[j] ~ dnorm(m1*Xsyn[j]+c1, tau1)
   Ysyn_def[j] ~ dnorm(m2*Xsyn[j]+c2, tau2)
  }
  
  # Specify prior beliefs about parameters
  m1 ~ dnorm(mu_m1, tau_m1)
  m2 ~ dnorm(mu_m2, tau_m2)
  c1 ~ dnorm(mu_c1, tau_c1)
  c2 ~ dnorm(mu_c2, tau_c2)
  
  tau1 ~ dgamma(shape_tau1, rate_tau1)
  var2 ~ dgamma(shape_tau2, rate_tau2)
  tau2 = 1/var2
  
  p ~ dbeta(alpha, beta)
  probdiff = ifelse(ctrl_ind==1, 0, p)
 }
"

# inference function
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, pts=pat)
    
    ctrl_mat = data_mats$ctrl
    pat_mat = data_mats$pts
    Nctrl = nrow(ctrl_mat)
    Npat = nrow(pat_mat)
    
    # prior parameters for control data
    c1_est = 0
    tau_c1 = 1/2^2
    m1_est = 0
    tau_m1 = 1/2^2
    
    c2_est = 0
    tau_c2 = 1/2^2
    m2_est = 0
    tau_m2 = 1/2^2
    
    shape_tau1 = 10
    rate_tau1 = 1
    shape_tau2 = 20 
    rate_tau2 = 10 
    
    N_syn = 1e4
    X_syn = seq(0, max(c(ctrl_mat[,1], pat_mat[,1]))*1.5, length.out=N_syn)
    
    ## control inference
    data_ctrl = list( Xobs=ctrl_mat[,1], Yobs=ctrl_mat[,2], N=Nctrl, 
                      Nsyn=N_syn, Xsyn=X_syn,
                      mu_m1=m1_est, tau_m1=tau_m1, mu_c1=c1_est, tau_c1=tau_c2,
                      mu_m2=m2_est, tau_m2=tau_m2, mu_c2=c2_est, tau_c2=tau_c2,
                      shape_tau1=shape_tau1, rate_tau1=rate_tau1,
                      shape_tau2=shape_tau2, rate_tau2=rate_tau2,
                      alpha=1.0, beta=1.0,
                      ctrl_ind=1)
    
    data_ctrl_priorpred = data_ctrl
    data_ctrl_priorpred$Yobs = NULL
    data_ctrl_priorpred$N = 0
    
    model_ctrl = jags.model(textConnection(modelstring), data=data_ctrl, n.chains=n.chains)
    model_ctrl_priorpred = jags.model(textConnection(modelstring), data=data_ctrl_priorpred)
    
    update(model_ctrl,n.iter=MCMCBurnin)
    output_ctrl = coda.samples(model=model_ctrl, n.iter=MCMCOut*MCMCThin,thin=MCMCThin,
                               variable.names=c("m1", "c1", "tau1", "m2", "c2", "tau2" ,"Ysyn_norm", "Ysyn_def", "class", "probdiff"))
    output_ctrl_prior = coda.samples(model=model_ctrl_priorpred, n.iter=MCMCOut,thin=1,
                                     variable.names=c("m1", "c1", "tau1", "m2", "c2", "tau2", "Ysyn_norm", "Ysyn_def", "probdiff"))
    
    posterior_ctrl = as.data.frame(output_ctrl[[1]])
    prior_ctrl = as.data.frame(output_ctrl_prior[[1]])
    summ_ctrl = summary(output_ctrl)
    classifs_ctrl = summ_ctrl$statistics[grepl("class",rownames(summ_ctrl$statistics)),"Mean"]
    
    ### patient inference
    # pateint priors
    flex = 10
    c1_est = mean(posterior_ctrl$c1)
    tau_c1 = flex/(sd(posterior_ctrl$c1)^2)
    m1_est = mean(posterior_ctrl$m1)
    tau_m1 = flex/(sd(posterior_ctrl$m1)^2)
    
    c2_est = 0
    tau_c2 = 1/2^2
    m2_est = 0
    tau_m2 = 1/2^2
    
    # delta = 1.5*as.numeric(quantile(posterior_ctrl$tau_ctrl,0.5)) # Choose this value so that Tiago's replication dataset never predicts over-expression of CI or CIV
    tau1_mean = mean(posterior_ctrl$tau1)  
    tau1_sd = sd(posterior_ctrl$tau1) # Deviation from prior tau should require a lot of contradictory data
    tau1_shape = (tau1_mean^2)/(tau1_sd^2)
    tau1_rate = tau1_mean/(tau1_sd^2)
    
    tau2_mode = 20 
    tau2_mean = 20
    tau2_sd = 5
    # I have no idea where these have come from..?
    # rate_gamma = (gamma_mode+sqrt(gamma_mode^2+4*gamma_mode^2))/(2*gamma_sd^2)
    # shape_gamma = 1+gamma_mode*rate_gamma
    tau2_rate = tau2_mean^2 / tau2_sd^2
    tau2_shape = tau2_mean / tau2_sd^2
    
    data_pat = list(Xobs=pat_mat[,1], Yobs=pat_mat[,2], N=Npat, Nsyn=N_syn,
                    Xsyn=X_syn,
                    mu_m1=m1_est, tau_m1=tau_m1, mu_m2=m2_est, tau_m2=tau_m2, 
                    mu_c1=c1_est, tau_c1=tau_c1, mu_c2=c2_est, tau_c2=tau_c2,
                    shape_tau1=tau1_shape, rate_tau1=tau1_rate, 
                    shape_tau2=tau2_shape, rate_tau2=tau2_rate, 
                    alpha=1, beta=1,
                    ctrl_ind=0)
    
    data_pat_priorpred = data_pat
    data_pat_priorpred$Yobs = NULL
    data_pat_priorpred$N = 0
    
    model_pat = jags.model(textConnection(modelstring), data=data_pat, n.chains=n.chains)
    model_pat_priorpred = jags.model(textConnection(modelstring), data=data_pat_priorpred)
    update(model_pat, n.iter=MCMCBurnin)
    # converge_pat=coda.samples(model=model_pat,variable.names=c("m","c","tau_par","class","probdiff"),n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
    output_pat_post = coda.samples(model=model_pat, n.iter=MCMCOut*MCMCThin, thin=MCMCThin,
                                   variable.names=c("m1", "c1", "tau1", "m2", "c2", "tau2","Ysyn_norm", "Ysyn_def", "class", "probdiff"))
    output_pat_prior = coda.samples(model=model_pat_priorpred,n.iter=MCMCOut, thin=1,
                                    variable.names=c("m1", "c1", "tau1", "m2", "c2", "tau2","Ysyn_norm", "Ysyn_def", "probdiff"))
    
    posterior_pat = as.data.frame(output_pat_post[[1]])
    prior_pat = as.data.frame(output_pat_prior[[1]])
    
    summ_pat = summary(output_pat_post)
    classifs_pat = summ_pat$statistics[grepl("class",rownames(summ_pat$statistics)),"Mean"]
    
    posterior_ctrl_names = colnames(posterior_ctrl)
    post_ctrl = posterior_ctrl[,c("m1", "c1", "tau1", "m2", "c2", "tau2", "probdiff")]
    postpred_ctrl_norm = colQuantiles( posterior_ctrl[,grepl("Ysyn_norm", posterior_ctrl_names)], probs=c(0.025, 0.5, 0.975) )
    postpred_ctrl_def = colQuantiles( posterior_ctrl[,grepl("Ysyn_def", posterior_ctrl_names)], probs=c(0.025, 0.5, 0.975) )
    postpred_ctrl = cbind(X_syn, postpred_ctrl_norm, postpred_ctrl_def)
    colnames(postpred_ctrl) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    prior_ctrl_names = colnames(prior_ctrl)
    prior_control = prior_ctrl[,c("m1", "c1", "tau1", "m2", "c2", "tau2", "probdiff")]
    priorpred_ctrl_norm = colQuantiles(prior_ctrl[, grepl("Ysyn_norm", prior_ctrl_names)], probs=c(0.025,0.5,0.975))
    priorpred_ctrl_def = colQuantiles(prior_ctrl[, grepl("Ysyn_def", prior_ctrl_names)], probs=c(0.025,0.5,0.975))
    priorpred_ctrl = cbind(X_syn, priorpred_ctrl_norm, priorpred_ctrl_def)
    colnames(priorpred_ctrl) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    posterior_pat_names = colnames(posterior_pat)
    post_pat = posterior_pat[,c("m1", "c1", "tau1", "m2", "c2", "tau2", "probdiff")]
    postpred_pat_norm = colQuantiles(posterior_pat[,grepl("Ysyn_norm", posterior_pat_names)], probs=c(0.025,0.5,0.975))
    postpred_pat_def = colQuantiles(posterior_pat[,grepl("Ysyn_def", posterior_pat_names)], probs=c(0.025,0.5,0.975))
    postpred_pat = cbind(X_syn, postpred_pat_norm, postpred_pat_def)
    colnames(postpred_pat) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    
    prior_pat_names = colnames(prior_pat)
    prior_patient = prior_pat[,c("m1", "c1", "tau1", "m2", "c2", "tau2", "probdiff")]
    priorpred_pat_norm = colQuantiles(prior_pat[,grepl("Ysyn_norm", prior_pat_names)], probs=c(0.025,0.5,0.975))
    priorpred_pat_def = colQuantiles(prior_pat[,grepl("Ysyn_def", prior_pat_names)], probs=c(0.025,0.5,0.975))
    priorpred_pat = cbind(X_syn, priorpred_pat_norm, priorpred_pat_def)
    colnames(priorpred_pat) = c("mitochan", "lwr_norm", "mid_norm", "upr_norm", "lwr_def", "mid_def", "upr_def")
    
    ctrl_list = list(post=post_ctrl, postpred=postpred_ctrl, 
                     prior=prior_control, priorpred=priorpred_ctrl,
                     classif=classifs_ctrl)
    pat_list = list(post=post_pat, postpred=postpred_pat,
                    prior=prior_patient, priorpred=priorpred_pat,
                    classif=classifs_pat)
    
    return( list(ctrl=ctrl_list, pat=pat_list) )
  })
}

cord = c("NDUFB8", "NDUFA13", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1", "OSCP")
mitochan = "VDAC1"

dat = getData("../../dat.txt", cord)

ptsAll = unique(dat$patient_id)
pts = sort( ptsAll[grepl("P", ptsAll)] )

inputs = list()
{
  input0 = list()
  input0$MCMCOut = 1000
  input0$MCMCBurnin = 1000
  input0$MCMCThin = 50
  input0$n.chains = 1
  
  for(chan in cord){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      inputs[[outroot]] = input0
      inputs[[outroot]]$chan = chan
      inputs[[outroot]]$pat = pat
    } # pts
  } # chans
}

ncores = detectCores() - 1
cl  = makeCluster(ncores)
{
  clusterExport(cl, c("modelstring"))
  clusterEvalQ(cl, {
    library("rjags")
    source("helper_functions.R", local=TRUE)
  })
  linreg_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

for(outroot in names(linreg_output)){
  output_saver(linreg_output[[outroot]], folder, outroot)
}


