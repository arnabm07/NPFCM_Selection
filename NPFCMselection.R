pprocess<-function(df){
  covlist<-list()
  n=length(unique(df$id))
  T<-list()
  for(i in 1:n)
  {ind<-which(df$id==i)
  T[[i]]<-df$time[ind]
  }
  m<-length(unique(unlist(T)))
  c=ncol(df)-2
  for(j in 1:c)
  {
    covlist[[j]]<-list()
    for(i in 1:n)
    {ind<-which(df$id==i)
    covlist[[j]][[i]]<-df[ind,j+2]  #rnorm(length(T[[i]]),0,1)#
    }
    lenX<-c()
    for (i in 1:n)
    {lenX[i]<-length(T[[i]])}
    .id=rep(1:n,lenX)
    .index<-unlist(T)
    Wdata1 <- data.frame(.id ,.index)
    x1<-covlist[[j]]
    Wdata1$.value<-unlist(x1)
    library(refund)
    w.sm1= fpca.sc(ydata=Wdata1, pve=.98, var=TRUE,center = TRUE)
    xhat1un= w.sm1$Yhat
    xhat1 = xhat1un
    timepoints<-as.numeric(attr(w.sm1$Yhat,"dimnames")[[2]])
    Xhat1_resp<-list()
    for(i in 1:n)
    {Xhat1_resp[[i]]<-vector(mode = "logical",length = length(T[[i]]))
    for(l in 1: length(T[[i]])){
      ind<-which( round(timepoints,4)==round(T[[i]],4)[l])
      Xhat1_resp[[i]][l]=xhat1[i,ind]
    }
    }
    df[,j+2]<-unlist(Xhat1_resp)
  }
  return(df)
}
NPFCM.select<-function(y,mydata,basist=c(4:10),basisx=c(4:10))
{
  name<-names(mydata)
  name<-name[-c(1,2)]
  p<- ncol(mydata)-2
  n=length(unique(mydata$id))
  T<-list()
  for(i in 1:n)
  {ind<-which(mydata$id==i)
  T[[i]]<-mydata$time[ind]
  }
  UT<-c()
  for( i in 1:n)
  {UT<-union(UT,T[[i]])}
  UT<-sort(UT)
  tlenbas<-basist-2
  xlenbas<-basisx-2
  kl<-expand.grid(tlenbas,xlenbas)
  library(parallel)
  ncores <-detectCores()-1
  cl <- makeCluster(ncores)
  ebiclen<-function(klindex)
  {
    library(fda)
    knots = seq(min(UT),max(UT),l=kl[klindex,1])
    norder = 4
    nbasis = length(knots) + norder - 2
    dayrng = c(min(UT),max(UT))
    bbasis = create.bspline.basis(dayrng,nbasis,norder,knots)
    #bbasis$nbasis    # number of basis functions
    #bbasis$rangeval   # basis range
    #plot(bbasis)
    bbasisMat<-list()
    for (i in 1:n)
    {bbasisMat[[i]] = eval.basis(T[[i]],bbasis)}
    basismat<-bbasisMat
    for (i in 1:n){
      colnames(basismat[[i]])<-NULL}
    B<-basismat
    for(i in 1:n)
    {B[[i]]<-cbind(rep(1,nrow(B[[i]])), B[[i]])
    }
    GrandB<-NULL
    for(i in 1:n)
    {GrandB<-rbind(GrandB, B[[i]])
    }
    GrandZlist<-list()
    ### make xj basises
    for(j in 1:p){
      library(fda)
      #timepoints<-as.numeric(attr(X.sm[[j]]$xhat,"dimnames")[[2]])
      Xhat_resp<-list()
      for(i in 1:n)
      {Xhat_resp[[i]]<-vector(mode = "logical",length = length(T[[i]]))
      ind<-which(mydata$id==i)
      Xhat_resp[[i]]<-mydata[ind,(j+2)]
      }
      xvec<-unlist(Xhat_resp)
      knotsx = seq(min(xvec),max(xvec),l=kl[klindex,2])
      norderx = 4
      # this implies the number of basis functions
      nbasisx = length(knotsx) + norderx - 2
      dayrngx = c(min(xvec),max(xvec))
      bbasisx = create.bspline.basis(dayrngx,nbasisx,norderx,knotsx)
      bbasisMatx<-list()
      for (i in 1:n)
      {bbasisMatx[[i]] = eval.basis(Xhat_resp[[i]],bbasisx)}
      XMatstar<-list()
      for(i in 1:n)
      {
        XMatstar[[i]]<-matrix(0,length(T[[i]]),(nbasis+1)*(nbasisx))
        for(k in 1:length(T[[i]]))
        {XMatstar[[i]][k,]<-bbasisMatx[[i]][k,]%x%B[[i]][k,]}
      }
      
      Grandz1<-NULL
      for(i in 1:n)
      {Grandz1<-rbind(Grandz1, XMatstar[[i]])
      }
      GrandZlist[[j]]<-Grandz1
    } 
    Grandcov<-Reduce(cbind,GrandZlist)
    GrandZ<-cbind(GrandB,Grandcov)
    group<-c()
    tbasdim<-(nbasis+1)*nbasisx
    nbasisint<-nbasis+1
    d<-tbasdim-nbasisint
    for (j in 1:(p+1))
    {if(j==1)
    {group[((nbasisint*(j-1))+1):(nbasisint*j)]<-0}
      if(j>1) 
      {group[((tbasdim*(j-1))+1-d):((tbasdim*j)-d)]<-j-1}
    } 
    library(grpreg)
    Groupvar<-as.factor(group)
    fit2 <- grpreg(GrandZ,GrandY,Groupvar,penalty="grSCAD")
    fit3 <- grpreg(GrandZ,GrandY,Groupvar,penalty="grMCP")
    cvfit2<-select(fit2,crit="EBIC")
    cvfit3<-select(fit3,crit="EBIC")
    gamma2<-cvfit2$beta
    gamma3<-cvfit3$beta
    indsel2<-which(gamma2[-1]!=0)
    Groupsel2<-group[indsel2]
    grind2<-unique(Groupsel2)
    indsel3<-which(gamma3[-1]!=0)
    Groupsel3<-group[indsel3]
    grind3<-unique(Groupsel3)
    varselected2<-grind2
    varselected3<-grind3
    if(is.finite(min(cvfit2$IC)))
    {ebicscad<-min(cvfit2$IC)}
    if(!is.finite(min(cvfit2$IC)))
    {ebicscad<-9999999}
    if(is.finite(min(cvfit3$IC)))
    {ebicmcp<-min(cvfit3$IC)}
    if(!is.finite(min(cvfit3$IC)))
    {ebicmcp<-9999999}
    #varselected2
    #varselected3
    lambdascad<-cvfit2$lambda
    lambdamcp<-cvfit3$lambda
    fit22 <- grpreg(GrandZ,GrandY,Groupvar,penalty="grSCAD",lambda = lambdascad)
    fit23 <- grpreg(GrandZ,GrandY,Groupvar,penalty="grMCP",lambda = lambdamcp)
    Yhats<-predict(fit22,GrandZ)
    Yhatm<-predict(fit23,GrandZ)
    #plot(GrandY,Yhatm)
    #a<-mean((GrandY-Yhatm)^2)
    #b<-mean((GrandY-mean(GrandY))^2)
    #1-(a/b)
    result<-list(varselected2,varselected3,ebicscad,ebicmcp,Yhats,Yhatm)
    names(result)<-c("varselected2","varselected3","ebicscad","ebicmcp","Yhats","Yhatm")
    return(result)
  }
  clusterExport(cl,c())
  clusterEvalQ(cl, library("fda"))
  clusterEvalQ(cl, library("grpreg"))
  clusterExport(cl, list("T", "n","p","GrandY","kl","mydata","UT"), envir=environment())
  klindex<-c(1:nrow(kl))
  cvresultfinal<- parLapply(cl,klindex,ebiclen)
  stopCluster(cl)
  len<-length(klindex)
  cves<-c()
  cvem<-c()
  for (k in 1:len)  
  {
    cves<-c(cves,cvresultfinal[[k]]$ebicscad)
    cvem<-c(cvem,cvresultfinal[[k]]$ebicmcp) 
  }
  indscad<-which.min(unlist(cves))
  indmcp<-which.min(unlist(cvem))
  varselected2<-cvresultfinal[[indscad]]$varselected2
  varselected3<-cvresultfinal[[indmcp]]$varselected3
  #Yhatscad<-cvresultfinal[[indscad]]$Yhats  ##return prediction if needed
  #Yhatmcp<-cvresultfinal[[indmcp]]$Yhatm
  namescad<-name[varselected2[-1]]
  namemcp<-name[varselected3[-1]]
  result<-list(namescad,namemcp)
  names(result)<-c("NFSCAD","NFMCP")
  return(result)
}