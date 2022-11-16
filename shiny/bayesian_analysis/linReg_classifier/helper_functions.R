library(data.table)
library(MASS)

##
## DATA COLLEECTION
##

getData = function(fname, cord, mitochan="VDAC1"){
  dat = fread(fname, sep="\t", stringsAsFactors=FALSE, header=TRUE)
  
  dat$channel = gsub("GRIM19","NDUFA13",dat$channel)
  
  dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
  dat = dat[dat$ch %in% c(cord, mitochan, "Area","AspectRatio","Circularity","Perimeter","xCoord","yCoord"),]
  
  dat$type = "Mean intensity"
  dat$type[grepl("LOG_",dat$channel)] = "Log mean intensity"
  dat$type[grepl("MED_",dat$channel)] = "Median intensity"
  dat$type[grepl("R_",dat$channel)] = "Ratio mean intensity (VDAC1)"
  dat$type[grepl("R_MED_",dat$channel)] = "Ratio median intensity (VDAC1)"
  dat$type[grepl("R_LOG_",dat$channel)] = "Ratio log mean intensity (VDAC1)"
  dat$type[grepl("Z_",dat$channel)] = "z-score"
  dat$type[dat$channel=="Area"] = "Area"
  dat$type[dat$channel=="AspectRatio"] = "AspectRatio"
  dat$type[dat$channel=="Circularity"] = "Circularity"
  dat$type[dat$channel=="Perimeter"] = "Perimeter"
  dat$type[dat$channel=="xCoord"] = "xCoord"
  dat$type[dat$channel=="yCoord"] = "yCoord"
  
  dat$outlier_diff = "NODIFF"
  dat$regression_diff = "NODIFF"
  dat$z_diff = "NODIFF"
  dat$z = 0

  dat$chstr = dat$ch
  transform = log
  dat_r = dat[dat$type=="Mean intensity",]
  dat_r$type = "r (VDAC1)"
  dat_theta = dat[dat$type=="Mean intensity",]
  dat_theta$type = "theta (VDAC1)"

  for(pid in unique(dat$patrep_id)){
    for(ch in unique(dat$ch)){
      dt = dat[(dat$patrep_id==pid)&(dat$type=="Mean intensity"),]

      isch = as.character(dt$ch)==ch
      ismito = as.character(dt$ch)==mitochan
      prot = dt[isch,]
      mito = dt[ismito,]

      x = mito$value
      y = prot$value
      dat_r$value[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = sqrt(x^2+y^2)
      dat_r$channel[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = paste("RADIUS",ch,sep="_")
      dat_theta$value[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = 360*atan(y/x)/(2*pi)
      dat_theta$channel[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = paste("THETA",ch,sep="_")
    }
  }
  dat=rbind(dat,dat_r,dat_theta)
  return(dat)
}

getData_mats = function(chan, mitochan="VDAC1", 
                        pts=NULL, ctrl_only=FALSE, 
                        data_transform=log_transform){

  dat = getData("../../dat.txt", chan)
  
  sbj = unique(dat$patient_id)
  crl = sbj[grepl("C", sbj)]
  Xctrl = dat[dat$channel==paste0("LOG", "_", mitochan) & dat$patient_type=="control", "value"]
  Yctrl = dat[dat$channel==paste0("LOG", "_", chan) & dat$patient_type=="control", "value"]
  
  ctrl_mat = as.matrix( cbind( Xctrl, Yctrl ) )

  if(!ctrl_only){
    if(is.null(pts)){
      pts = sort(sbj[grepl("P", sbj)])
      Ypts = list()
      for(pat in pts){
        Xpat = dat[dat$channel==paste("LOG", mitochan, sep="_") & dat$patient_id == pat, "value"]
        Ypat = dat[dat$channel==paste("LOG", chan, sep="_") & dat$patient_id == pat, "value"]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
      }
    } else {
      Ypts = list()
      for(pat in pts){
        Xpat = dat[dat$channel==paste("LOG", mitochan, sep="_") & dat$patient_id == pat, "value"]
        Ypat = dat[dat$channel==paste("LOG", chan, sep="_") & dat$patient_id == pat, "value"]
        XY_pat = cbind(Xpat, Ypat)
        Ypts[[pat]] = XY_pat
      }
    }
    pat_mat = list2matrix(Ypts)
    pat_mat = as.matrix(pat_mat)
  } else { pat_mat=NULL }
  
  if(!is.null(data_transform)) return( data_transform(ctrl_mat, pat_mat) )
  
  if(ctrl_only) return( ctrl_mat )
  list(ctrl=ctrl_mat, pts=pat_mat)
}

getData_Npts = function(chan, mitochan="VDAC1"){
  dat = getData("../../dat.txt", chan)
  
  sbj = unique(dat$patient_id)
  pts = sort(sbj[grepl("P", sbj)])
  Npts = double(length(pts))
  names(Npts) = pts
  for(pat in pts){
    Npts[pat] = nrow(dat[dat$channel==paste("LOG", mitochan, sep="_") & dat$patient_id == pat, ])
  }
  Nptts
}

index_creator = function(Npops, ind){
  index = rep(FALSE, sum(Npops))
  if( is.character(ind) ){
    ind = which(names(Npops)==pat)
    names(Npops) = NULL
  } else if( is.null(names(Npops)) & is.character(ind)){
    stop("if ind is a character Npops must be a named vector")
  }
  if(ind>1 & ind<length(Npops)){
    return((sum(Npops[1:(ind-1)])+1):sum(Npops[1:ind]))
  } else if(ind==1){
    return(1:Npops[1])
  } else {
    return( (sum(Npops[1:(ind-1)])+1):sum(Npops) ) 
  }
}

##
## COLOURS
##

pal = palette()
palette(c(pal, "darkblue", "darkorange", "deeppink"))

myBlack = function(alpha) rgb(0,0,0, alpha)
myDarkGrey = function(alpha) rgb(169/255,169/255,159/255, alpha)
myGrey = function(alpha) rgb(66/255,66/255,66/255, alpha)
myBlue = function(alpha) rgb(0,0,128/255, alpha)
myRed = function(alpha) rgb(1,0,0, alpha)

myGreen = function(alpha) rgb(0,100/255,0, alpha)
myYellow = function(alpha) rgb(225/255,200/255,50/255, alpha)
myPink = function(alpha) rgb(255/255,62/255,150/255, alpha)
myPurple = function(alpha) rgb(160/255, 32/255, 240/255, alpha)

cramp = colorRamp(c(myBlue(0.2), myRed(0.2)), alpha=TRUE)
# rgb(...) specifies a colour using standard RGB, where 1 is the maxColorValue
# 0.25 determines how transparent the colour is, 1 being opaque 
# cramp is a function which generates colours on a scale between two specifies colours

classcols = function(classif){
  # A function using the cramp specified colours to generate rgb colour names
  # input: a number in [0,1]
  # output: rgb colour name
  rgbvals = cramp(classif)/255.0
  return(rgb(rgbvals[,1],rgbvals[,2],rgbvals[,3], alpha=rgbvals[,4]))
}


##
## GEN FUNCTIONS
##

colQuantiles = function(x, probs=0.5){
  quants = matrix(NA, nrow=ncol(x), ncol=length(probs))
  for(i in 1:ncol(x)){
    quants[i,] = quantile(x[,i], probs)
  }
  colnames(quants) = probs
  return(quants)
}

vec_rep = function(x, n, byRow=TRUE){
  out = matrix(x, nrow=n, ncol=length(x), byrow=TRUE)
  if(byRow) return( out )
  else return( t(out) )
}

list2matrix = function(X, rowBind=TRUE){
  # X must list of matrices or a matrix
  if(is.matrix(X)) return(X)
  if(!is.list(X)) stop("X must be a list of matrices or a matrix")
  if(length(X)==1) return(X[[1]])
  
  n = length(X)
  out = X[[1]]
  if(rowBind){
    for(i in 2:n){ out = rbind(out, X[[i]]) }
    return(out)
  } else {
    for(i in 2:n){ out = cbind(out, X[[i]]) }
    return(out)
  }
}

log_transform = function(ctrl_mat, pat_mat=NULL){
  if( !is.null(pat_mat) ) return( list(ctrl=log(ctrl_mat), pts=log(pat_mat)) )
  log( ctrl_mat )
}

##
## READERS & SAVERS AOUTPUT
##

output_saver = function(output, folder, outroot){
  split = strsplit(outroot, split="_")[[1]]
  chan = split[1]
  pat = split[2]
  
  for(ctrl_pat in names(output)){
    out_ctrlpat = ifelse(ctrl_pat=="ctrl", "CONTROL", pat)
    for(out_type in names(output[[ctrl_pat]])){
      filename = paste(chan, out_ctrlpat, toupper(out_type), sep="_")
      write.table(output[[ctrl_pat]][[out_type]], paste0(file.path("Output", folder, filename), ".txt"),
                  row.names=FALSE, col.names=TRUE)
    }
  }
}

output_reader = function(folder, chan, pat, out_type){
  outroot = paste(chan, pat, sep="_")
  
  fp = file.path("Output", folder, paste0(outroot, "_", out_type, ".txt"))
  if(file.exists(fp)) return( read.table(fp, header=TRUE, stringsAsFactors=FALSE) )
  else stop(paste("file does not exist\n", fp))
}

##
## PLOTTERS
##

classif_plot = function(ctrl_data, pat_data, classifs_pat, chan, mitochan, pat){
    
    xrange = range(c(ctrl_data[,1], pat_data[,1]))
    yrange = range(c(ctrl_data[,2], pat_data[,2]))
    
    plot(ctrl_data, pch=20, col=myDarkGrey(0.2), xlim=xrange, ylim=yrange, 
         xlab=paste(mitochan), ylab=paste(chan), 
         main=paste(pat), cex.lab=1.2, cex.main=1.4)
    points(pat_data, pch=20, col=classcols(classifs_pat))
  }
  
MCMCplot = function(post, prior, ctrl_prior=NULL, title="", lag=20){
    col.names = colnames(post)
    n.chains = length(post)

    for(param in col.names){
      post_vec = post[[param]]
      plot(ts(post_vec), xlab="Iteration", ylab=paste(param), 
           main="", cex.lab=1.2)
      if(sum(post_vec==post_vec[1])!=length(post_vec)){
        acf(post[[param]], xlab="lag index", ylab="ACF", main="",
            cex.lab=1.2, lag=lag)
      } else {
        plot(NA, type='n', xlim=c(0,lag), ylim=c(0,1), 
             xlab="lag index", ylab="ACF", main="")
      }
      
      xlims = range(post[[param]])
      ylims = c(0,max( density(post[[param]])$y ))
      if(param %in% colnames(prior) & is.null(ctrl_prior)){
        xlims = range(c(post[[param]], prior[[param]]))
      #   ylims = c(0,max(c(density(post[[param]])$y, density(prior[[param]])$y)))
      } else if( !(param %in% colnames(prior)) & (!is.null(ctrl_prior) & param %in% colnames(ctrl_prior))){
        xlims = range(c(post[[param]], ctrl_prior[[param]]))
      #ylims = c(0,max(c(density(post[[param]])$y, density(ctrl_prior[[param]])$y)))
      } else if( param %in% colnames(prior) & (!is.null(ctrl_prior) & param %in% colnames(ctrl_prior))){
        xlims = range(c(post[[param]], prior[[param]], ctrl_prior[[param]]))
      # ylims = c(0,c(density(post[[param]])$y, density(prior[[param]])$y, density(ctrl_prior[[param]])$y))
      }
      
      plot(density(post[[param]]), lwd=2, col="blue", xlab=paste(param), ylab="Density",
           main="", xlim=xlims, ylim=ylims*2)
      if(param %in% colnames(prior)) lines(density(prior[[param]]), lwd=2, col="forestgreen")
      if(!is.null(ctrl_prior) & param %in% colnames(ctrl_prior)){
        lines(density(ctrl_prior[[param]]), lwd=2, col="darkorange")
      }
      legend("topright", col=c("blue", "forestgreen", "darkorange"), 
             lty=1, lwd=2, legend=c("pat post", "pat prior", "ctrl post"))
      
      title(main=title, line=-1, outer=TRUE, cex.main=1.6)
    }
}

priorpost = function(ctrl_data, pat_data=NULL, priorpred, postpred,
                     classif=NULL, 
                     chan, mitochan="VDAC1", title="", xlims=NULL, ylims=NULL){

  op = par(mfrow=c(1,2), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  plot(ctrl_data, pch=20, cex=0.7, col=myGrey(0.1),
        xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Prior Predictive", xlim=xlims, ylim=ylims)
  if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=myYellow(0.2))
  lines(priorpred[,"mitochan"], priorpred[,"lwr_norm"], lty=2, col=myGreen(0.6), lwd=3)
  lines(priorpred[,"mitochan"], priorpred[,"mid_norm"], lty=1, col=myGreen(0.6), lwd=4)
  lines(priorpred[,"mitochan"], priorpred[,"upr_norm"], lty=2, col=myGreen(0.6), lwd=3)
  
  plot(ctrl_data, pch=20, col=myGrey(0.1),
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Posterior Predictive", xlim=xlims, ylim=ylims)
  if(!is.null(pat_data)) points(pat_data, pch=20, cex=1.2, col=classcols(classif))
  lines(postpred[,"mitochan"], postpred[,"lwr_norm"], lty=2, col=myPink(0.6), lwd=3)
  lines(postpred[,"mitochan"], postpred[,"mid_norm"], lty=1, col=myPink(0.6), lwd=4)
  lines(postpred[,"mitochan"], postpred[,"upr_norm"], lty=2, col=myPink(0.6), lwd=3)
  title(main=title, line=-2, outer=TRUE)
  
  par(op)
}

compare_preds = function(ctrl_data, pat_data, post, classif, 
                     chan, mitochan="VDAC1", title="", xlims=NULL, ylims=NULL){
  
  op = par(mfrow=c(1,2), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  plot(ctrl_data, pch=20, cex=0.7, col=myGrey(0.1),
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Prior Predictive", xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, cex=1.2, col=classcols(classif))
  lines(post[,"mitochan"], post[,"lwr_def"], lty=2, col=myRed(0.6), lwd=3)
  lines(post[,"mitochan"], post[,"mid_def"], lty=1, col=myRed(0.6), lwd=4)
  lines(post[,"mitochan"], post[,"upr_def"], lty=2, col=myRed(0.6), lwd=3)
  
  plot(ctrl_data, pch=20, col=myGrey(0.1),
       xlab=paste0("log(",mitochan,")"), ylab=paste0("log(",chan,")"), 
       main="Posterior Predictive", xlim=xlims, ylim=ylims)
  points(pat_data, pch=20, cex=1.2, col=classcols(classif))
  lines(post[,"mitochan"], post[,"lwr_norm"], lty=2, col=myBlue(0.6), lwd=3)
  lines(post[,"mitochan"], post[,"mid_norm"], lty=1, col=myBlue(0.6), lwd=4)
  lines(post[,"mitochan"], post[,"upr_norm"], lty=2, col=myBlue(0.6), lwd=3)
  title(main=title, line=-2, outer=TRUE)
  
  par(op)
}


pipost_plotter = function(chan, folder, pts, alpha=0.05){
  npat = length(pts)
  pis = list()

  for(pat in pts){
    pis[[pat]] = output_reader(folder, chan, pat, out_type="POST")[,"probdiff"]
  }

  title = paste( chan )

  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=myBlack(alpha),
             group.names=pts,
             at=1:length(pts),
             main=title, ylim=c(0,1), ylab="Deficiency Proportion", 
             xlab="Patient Sample")
}

pipost_plotter_v2 = function(chan, folder, pts, alpha=0.01){
  
  Npts = length(pts)
  pis = list()
  boots = list()
  op = par(mfrow=c(1,1), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  for(pat in pts){
      boots_df = read.table(file.path("CharlotteWarrenBootstrap", "BootstrapParticles", paste0(pat, "_", chan, ".txt")),
                       header=TRUE, stringsAsFactors=FALSE)
      boots[[pat]] = as.vector(1 - as.vector(boots_df[,"med"]))
      pis[[pat]] = output_reader(folder, chan, pat, out_type="POST")[,"probdiff"]
  }
  
  stripchart(pis, pch=20, method="jitter", vertical=TRUE, 
             col=myGreen(alpha),
             at=1:Npts,
             ylim=c(0,1),
             main=chan, ylab="Deficiency Proportion", 
             xlab="Patient" )
  stripchart(boots, pch=20, method="jitter", vertical=TRUE, add=TRUE, 
             at=1:Npts,
             col=myPink(alpha))
  
  par(op)
}







































































