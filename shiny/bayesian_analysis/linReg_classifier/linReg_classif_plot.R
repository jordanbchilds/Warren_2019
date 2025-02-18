library(MASS)
source("helper_functions.R", local = TRUE)

folder = "linReg_classifier_flexTen"

dir.create("PDF", showWarnings=FALSE)
dir.create(file.path("PDF",folder), showWarnings=FALSE)

cord = c("NDUFB8", "NDUFA13", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1", "OSCP")
mitochan = "VDAC1"
dat = getData("../../dat.txt", cord)
ptsAll = unique(dat$patient_id)
pts = sort( ptsAll[grepl("P", ptsAll)] )

######################
### the plots
######################

pdf(file.path("PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  op = par(mfrow=c(2,3), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  for(chan in cord){
    outroot_ctrl = paste(chan, "CONTROL", sep="_")

    # MCMCplot(folder, chan, lag=100, 
    #            title=paste(chan, "CONTROL"))
    for(pat in pts){
      outroot_pat = paste(chan, pat, sep="_")
      
      prior = output_reader(folder, chan, pat, out_type="PRIOR")
      post = output_reader(folder, chan, pat, out_type="POST")
      ctrl_prior = output_reader(folder, chan, pat="CONTROL", out_type="POST")
      
      MCMCplot(post, prior, ctrl_prior, lag=100,
               title=paste(chan, pat))
    } # channel
  } # patient
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "model_post.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    data = getData_mats(chan=chan)
    ctrl_mat =  data$ctrl
    xlims = range(c(ctrl_mat[,1], data$pts[,1]))
    ylims = range(c(ctrl_mat[,2], data$pts[,2]))
    for(pat in pts){
      pat_mat = getData_mats(chan=chan, pts=pat)$pts
      
      priorpost(ctrl_data=ctrl_mat, pat_data=pat_mat,
                classif=output_reader(folder, chan, pat, out_type="CLASSIF")[[1]],
                priorpred=output_reader(folder, chan, pat, "PRIORPRED"), 
                postpred=output_reader(folder, chan, pat, "POSTPRED"),
                chan=chan, mitochan="VDAC1", title=paste("\n", chan, pat),
                xlims=xlims, ylims=ylims)
    } # patients
  } # channels
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "compare_preds.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    data = getData_mats(chan=chan)
    ctrl_mat =  data$ctrl
    xlims = range(c(ctrl_mat[,1], data$pts[,1]))
    ylims = range(c(ctrl_mat[,2], data$pts[,2]))
    for(pat in pts){
      pat_mat = getData_mats(chan=chan, pts=pat)$pts
      
      compare_preds(ctrl_data=ctrl_mat, pat_data=pat_mat,
                classif=output_reader(folder, chan, pat, out_type="CLASSIF")[[1]],
                post=output_reader(folder, chan, pat, "POSTPRED"),
                chan=chan, mitochan="VDAC1", title=paste("\n", chan, pat),
                xlims=xlims, ylims=ylims)
    } # patients
  } # channels
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "classifs.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    data = getData_mats(chan=chan)
    ctrl_mat =  data$ctrl
    xlims = range(c(ctrl_mat[,1], data$pts[,1]))
    ylims = range(c(ctrl_mat[,2], data$pts[,2]))
    for(pat in pts){
      pat_mat = getData_mats(chan=chan, pts=pat)$pts
      classifs_pat = output_reader(folder, chan, pat, "classif")[[1]]
      classif_plot(ctrl_data=ctrl_mat, pat_data=pat_mat, 
                   classifs_pat=classifs_pat, 
                   chan=chan, mitochan=mitochan, pat=pat)
    } # patients
  } # channels
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "pi_post.pdf"), width=13, height=8)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    pipost_plotter(chan, pts=pts, folder=folder) 
  }
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "pi_comp.pdf"), width=13, height=8)
{
  for(chan in cord){
    pipost_plotter_v2(chan, pts=pts, folder=folder) 
  }
}
dev.off()






