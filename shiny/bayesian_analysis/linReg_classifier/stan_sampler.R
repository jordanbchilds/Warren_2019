
# library("devtools")
# install_github("jordanbchilds/analysis2Dmito")

library("analysis2Dmito")
library("data.table")
library("parallel")
library("rstan")
source("stan_sampler_function.R", local=TRUE)

folder = "stan_sampler_multiChain"

dir.create("Output", showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

raw_data = fread("../../dat.txt")
raw_data[raw_data$channel=="GRIM19", "channel"] = "NDUFB13"

mitochan = "VDAC1"
channels = c("NDUFB8", "NDUFB13", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1", "OSCP")
nChains = length(channels)

data = as.data.frame(raw_data[raw_data$channel %in% c(mitochan, channels), ])

data$value = log(data$value)

sbjs = unique(data$patient_id)
ctrls = unique(data$patient_id[data$subject_group=="Control"])
pts = unique(data$patient_id[data$subject_group!="Control"])

# remove patients with nuclear encoded mutation because what's the point
pts = pts[!(pts %in% c("P01", "P02"))]

gradients = matrix(NA, nrow=length(channels), ncol=length(ctrls))
rownames(gradients) = channels
colnames(gradients) = ctrls

intercepts = gradients
precisions = gradients

xlim = range(data$value) 
for(chan in channels){
  op = par(mfrow=c(1,3))
  for( crl in ctrls ){
    x = data[data$patient_id==crl & data$channel==mitochan, "value"]
    y = data[data$patient_id==crl & data$channel==chan, "value"]
    
    df = data.frame(mitochan=x, chan=y)
    
    plot(df, pch=20, xlim=range(data$value), ylim=range(data$value),
         main=paste(chan, crl))
    
    mod = lm(chan~mitochan, data=df)
    xx = seq(xlim[1], xlim[2], length.out=1000)
    pred = predict.lm(mod, newdata=data.frame(mitochan=xx), interval="prediction")
    
    abline(a=mod$coefficients[1], b=mod$coefficients[2], lwd=3, col="red")
    lines(xx, pred[,"lwr"], lty=3, lwd=2, col="red")
    lines(xx, pred[,"upr"], lty=3, lwd=2, col="red")
    
    gradients[chan, crl] = mod$coefficients["mitochan"]
    intercepts[chan, crl] = mod$coefficients["(Intercept)"]
    precisions[chan, crl] = 1 / summary(mod)$sigma^2
  }
  par(op)
}

grad_mean = apply(gradients, 1, mean)
inter_mean = apply(intercepts, 1, mean)
precs_mean = apply(precisions, 1, mean)

grad_var = apply(gradients, 1, var)
inter_var = apply(intercepts, 1, var)
precs_var = apply(precisions, 1, var)

tau_m_vars = rep(50, nChains)^2
tau_c_vars = rep(50, nChains)^2

tau_vars = rep(10, nChains)^2

names(tau_m_vars) = channels
names(tau_c_vars) = channels
names(tau_vars) = channels

data_list = list()
for(chan in channels) {
  mean_mu_m = grad_mean[chan]
  prec_mu_m = 1 / 0.025 ^ 2
  
  mean_mu_c = inter_mean[chan]
  prec_mu_c = 1 / 0.25 ^ 2
  
  tau_m_mode = 1 / grad_var[chan]
  tau_m_var =  tau_m_vars[chan]
  rate_tau_m = 0.5 * (tau_m_mode + sqrt(tau_m_mode ^ 2 + 4 * tau_m_var)) / tau_m_var
  shape_tau_m = 1 + tau_m_mode * rate_tau_m
  
  tau_c_mode = 1 / inter_var[chan]
  tau_c_var = tau_c_vars[chan]
  rate_tau_c = 0.5 * (tau_c_mode + sqrt(tau_c_mode ^ 2 + 4 * tau_c_var)) / tau_c_var
  shape_tau_c = 1 + tau_c_mode * rate_tau_c
  
  tau_mode = precs_mean[chan]
  tau_var = tau_vars[chan]
  rate_tau = 0.5 * (tau_mode + sqrt(tau_mode ^ 2 + 4 * tau_var)) / tau_var
  shape_tau = 1 + tau_mode * rate_tau
  
  paramVals = list(
    shape_tau = shape_tau,
    rate_tau = rate_tau,
    shape_tau_c = shape_tau_c,
    rate_tau_c = rate_tau_c,
    shape_tau_m = shape_tau_m,
    rate_tau_m = rate_tau_m,
    mean_mu_m = mean_mu_m,
    prec_mu_m = prec_mu_m,
    mean_mu_c = mean_mu_c,
    prec_mu_c = prec_mu_c,
    alpha_pi = 0.0, 
    beta_pi = 1.0,
    slope_lb = 0.1,
    tau_def = 0.001
  )
  
  for (pat in pts) {
    rt = paste(chan, pat, sep="__")
    dd = data[, c("value", "channel", "id", "patient_id")]
    colnames(dd) = c("value", "channel", "fibreID", "sampleID")
    
    data_list[[rt]] = analysis2Dmito::getData_mats(
      dd,
      channels = c(mitochan, chan),
      ctrlID = ctrls,
      pts = pat
    )
    
    data_list[[rt]]$parameterVals = paramVals
  }
}

ncores = detectCores() - 1 
cl  = makeCluster(ncores)
{
  clusterExport(cl, "stan")
  output = parLapply(cl,
                     data_list,
                     stan_inference,
                     MCMCout = 4000, 
                     MCMCburnin = 2000,
                     nChains = 5,
                     max_logLik = FALSE
                     )
}
stopCluster(cl)

for (root in names(output)) {
  list_saver(output[[root]], file.path("Output", folder, root))
}





