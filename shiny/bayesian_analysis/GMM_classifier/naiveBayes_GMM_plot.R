library("MASS")
source("helper_functions_v0.R", local = TRUE)

folder = "naiveBayes_GMM"

dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PDF", folder), showWarnings = FALSE)

cord = c("NDUFB8", "NDUFA13", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1", "OSCP")
mitochan = "VDAC1"

dat = getData(file.path("..", "..", "dat.txt"), cord)
ptsAll = unique(dat$patient_id)
pts = sort( ptsAll[grepl("P", ptsAll)] )

prior = read.table(file.path("Output", folder, "PRIOR.txt"), 
                   header=TRUE, stringsAsFactors=FALSE)


pdf(file.path("PDF", folder, "MCMC.pdf"), width=13, height=8)
{ 
  for(chan in cord){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      MCMCplot(folder, outroot,  lag=100,
               title=chan)
    }
  }
}
dev.off()

pdf(file.path("PDF", folder, "model_post_RAW.pdf"), width=13, height=11)
{
  for(chan in cord){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      
      ctrl_data =  getData_mats(chan=chan, ctrl_only=TRUE, 
                                data_transform=myData_transform)
      
      pat_data = getData_mats(chan=chan, pts=pat,
                              data_transform=myData_transform)$pts
      
      classif_allpats = output_reader(folder, outroot, out_type="CLASS")[[1]]
      
      post = output_reader(folder, outroot, out_type="POST")
      
      priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                classif=classif_allpats,
                prior=prior, post=post,
                chan=chan, mitochan="porin", 
                title=chan)
    } # patients
  } # channels
}
dev.off()

pdf(file.path("PDF", folder, "model_post.pdf"), width=13, height=11)
{
  for(chan in cord){
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      ctrl_data = getData_mats(chan=chan, ctrl_only=TRUE)
      pat_data = getData_mats(chan=chan, pts=pat)$pts
      
      classif = output_reader(folder, outroot, out_type="CLASS")[[1]]
      post = output_reader(folder, outroot, out_type="POST")
      
      postpred = back_transform(post, ctrl_mat=ctrl_data, parameters=list(c("comp.1.1.", "comp.2.1."), c("comp.1.2.", "comp.2.2.")))
      priorpred = back_transform(prior, ctrl_mat=ctrl_data, parameters=list(c("comp.1.1.", "comp.2.1."), c("comp.1.2." ,"comp.2.1.")))
      
      priorpost(ctrl_data=ctrl_data, pat_data=pat_data,
                classif=classif,
                prior=priorpred, post=postpred,
                chan=chan, mitochan="VDAC1", title=outroot )
    }
  }
}
dev.off()

pdf(file.path("PDF", folder, "pi_post.pdf"), width=13, height=8)
{   
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(chan in cord){
    pipost_plotter(folder, chan) 
  }
  par(op)
}
dev.off()

pdf(file.path("PDF", folder, "classifs.pdf"), width=11, height=7)
{
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.lab=2, cex.main=2, cex.axis=1.5)
  for(chan in cord){
    data =  getData_mats(chan=chan, pts=pat )
    xlims = range(c(data$ctrl[,1], data$pts[,1]))
    ylims = range(c(data$ctrl[,2], data$pts[,2]))
    
    for(pat in pts){
      outroot = paste(chan, pat, sep="_")
      ctrl_data = data$ctrl
      pat_data = data$pts
      classif_pat = output_reader(folder, outroot, out_type="CLASS")[[1]]
      
      classif_plot(ctrl_data=ctrl_data, pat_data=pat_data,
                   classifs_pat=classif_pat,
                   chan=chan, mitochan="VDAC1", 
                   title=paste(chan, pat), 
                   xlims=xlims, ylims=ylims)
    }
  }
  par(op)
}
dev.off()




