library("parallel")
library("rstan")
source("helper_functions_v0.R", local=TRUE)

folder = "naiveBayes_GMM"

dir.create("Output", showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

cord = c("NDUFB8", "NDUFA13", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1", "OSCP")
mitochan = "VDAC1"

dat = getData("../../dat.txt", cord)

ptsAll = unique(dat$patient_id)
pts = sort( ptsAll[grepl("P", ptsAll)] )

# inferenecce 
inference = function(input){
  with(c(input),{
    data_mats = getData_mats(chan=chan, pts=pts, data_transform=myData_transform, 
                             get_patindex=TRUE)
    
    Yctrl = data_mats$ctrl
    Ypat = data_mats$pts
    Nctrl = nrow(Yctrl)
    Npat = nrow(Ypat)
    data_list = list(Yctr=Yctrl, Ypat=Ypat, Nctrl=Nctrl, Npat=Npat)
    
    output = stan("naiveBayes_GMM.stan", data=c(prior_list, data_list), chains=1, 
                  iter=(MCMCOut+MCMCBurnin)*MCMCThin, warmup=MCMCBurnin, thin=MCMCThin )
    
    outmat = as.matrix(output)
    outcols = colnames(outmat)
    
    classifs_df = outmat[, grepl("classif", outcols)]
    classifs_avg = colMeans(classifs_df)
    post = outmat[, !(grepl("classif", outcols)|grepl("lp__", outcols)|
                        grepl("probDef_vec", outcols)|grepl("dens", outcols))]
    
    return( list(post=post, classifs=classifs_avg) )
  })
}

{
  p = 2
  # prior parameters
  mean_mu1 = double(p)
  mean_mu2 = double(p)
  var_mu1=  matrix(c(0.5^2, 0, 0, 0.5^2), ncol=p, nrow=p, byrow=TRUE) 
  var_mu2 = 5^2*diag(p)
  
  n_1 = 1000
  S_1 = matrix(c(3^2, 0, 0, 1^2), nrow=p, ncol=p, byrow=TRUE)*(n_1-p-1)
  n_2 = 50
  S_2 = matrix(c(10^2,0,0,10^2), nrow=p, ncol=p, byrow=TRUE)*(n_2-p-1)
  
  alpha = 1
  beta = 1
}
prior_list = list(mean_mu1=mean_mu1, var_mu1=var_mu1,
                  mean_mu2=mean_mu2, var_mu2=var_mu2,
                  n_1=n_1, n_2=n_2, S_1=S_1, S_2=S_2,
                  alpha=alpha, beta=beta,
                  D=2, K=2)

# inputs for the inference function - which data file name, chain length, thinning etc.
inputs = list()
{
  input0 = list()
  input0$MCMCOut = 2000
  input0$MCMCBurnin = 2000
  input0$MCMCThin = 1
  input0$n.chains = 1
  for(chan in cord){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      inputs[[outroot]] = input0
      inputs[[outroot]]$chan = chan
      inputs[[outroot]]$pts = pat
      inputs[[outroot]]$prior_list = prior_list
    } # pts
  } # chans
}

ncores = detectCores()-1
cl  = makeCluster(ncores) 
{
  clusterEvalQ(cl, {
    library("rstan")
    source("helper_functions.R")
  })
  gmm_output = parLapply(cl, inputs, inference)
}
stopCluster(cl)

# save output
for(outroot in names(gmm_output)){
  output = gmm_output[[outroot]]
  write.table(output$post, file.path("Output", folder, paste0(outroot,"_POST.txt")), row.names=FALSE, col.names=TRUE)
  write.table(output$classif, file.path("Output", folder, paste0(outroot,"_CLASS.txt")), row.names=FALSE, col.names=TRUE)
}

### prior beliefs
output_prior = stan("naiveBayes_GMM_prior.stan", data=prior_list, chains=1, 
                    iter=1e4, warmup=0, thin=1)

prior_all= as.matrix(output_prior)
prior = prior_all[, !grepl("_", colnames(prior_all))]

write.table(prior, file.path("Output", folder, "PRIOR.txt"), 
            row.names=FALSE, col.names=TRUE)











