## Décomposition en gaussiennes : algorithme EM

graph_gauss <- function(valeurs,Nclust,Col){# valeurs est le tableau contenant ici les niveaux de gris ("v") index?s par leurs coordonn?es, aui figurent dans les deux premi?res colonnes.
  names(valeurs)[1]<-"x"
  names(valeurs)[2]<-"y"
  names(valeurs)[3]<-"v"
  valeurs.mix<-normalmixEM(valeurs$v,k=Nclust)
  mu_mix=valeurs.mix$mu
  sigma_cmix=valeurs.mix$sigma
  lambda_mix=valeurs.mix$lambda
  xx<-seq(from=0, to=length(valeurs$x), by=1)
  plot(x=xx,y=valeurs$v,type="o",col="black",xlab='ngris',ylab='fr?quence',main='Niveaux de gris')
  for (i in 1:Nclust){
    yy<-lambda_mix[i]*dnorm(x=xx,mean=mu_mix[i],sd=sigma_cmix[i])
    lines(yy,type="l",color=Col[i])
  }
}
# Difficile ? partir de 4 clusters : trop de tentatives

# Série d'histogrammes : jours et slices

hist_slice_essai <- function(Nclust,vars,slice){
  dir=paste("essai",vars,slice,sep='_')
  J<-read.table(paste(dir,slash,"jours.txt",sep=''))
  #J_num<-read.table(paste(dir,slash,"jours.txt",sep=''))
  l<-dim(J)[1]-1
  #J_char<-read.table(paste(dir,slash,"jours.txt",sep=''),as.is=rbind(rep(TRUE,length(J_num))))
  slice_corr<-matrix(nrow=l,ncol=l)
  for (j in 1:l){#J_char_2){# boucle sur les jours
    ngr_file_name=paste(dir,slash,"26-J",J[j,],"-",vars,"-","seg","-",slice,".txt",sep='')
    data=read.table(ngr_file_name)
    names(data)[1]<-"x"
    names(data)[2]<-"y"
    names(data)[3]<-"v"
    ### ? ###
    # On trace et on sauvegarde
    sous_dir<-paste(slash,"stat_pos",sep='')
    png(file = paste(dir,sous_dir,slash,"26-J",J[j,],"-hist-",toString(Nclust),".png",sep=''), width = 800, height = 700)
    par(mfrow=c(2,1))#, bg='white')
    hist(data$v, xlim=c(0,3000), freq = FALSE, main = paste("Histogramme",vars,"J",toString(J[j,]),sep=' '), xlab=vars, col="red", nclass=30)
    boxplot(data$v,Range=0,ylab=vars)
    dev.off()
  }
}


## Clusters avec Mclust
# Représentation graphique de clusters sur un jour et une slice

clust_slice <- function(Nclust,vars,slice){
  dir=paste("essai",vars,slice,sep='_')
  J<-read.table(paste(dir,slash,"jours.txt",sep=''))
  for (j in J){# boucle sur les jours
    ngr_file_name=paste(dir,slash,"26-J",j,"-",vars,"-","seg","-",slice,".txt",sep='')
    data=read.table(ngr_file_name)
    names(data)[1]<-"x"
    names(data)[2]<-"y"
    names(data)[3]<-"v"
    data.clust<-Mclust(data=data$v,G=Nclust,modelNames='V')
    # On trace et on sauvegarde
    png(file = paste(dir,slash,"26-J",j,"-clust-",toString(Nclust),".png",sep=''), width = 800, height = 700)
    par(mfrow=c(1,2))#, bg='white')
    plot(data$x,data$y,col=color.vector[data.clust$classification],pch=20,cex=(1-data.clust$uncertainty)^3,main=paste("Brain map of",vars,"data classification",sep=' '),xlab="x",ylab="y")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    points(data$x,data$y,col=color.vector[data.clust$classification], pch=20, cex=(1-data.clust$uncertainty)^3)
    plot(data.clust, what="classification", colors=color.vector, xlab=vars)
    dev.off()
  }
}

