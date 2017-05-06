### Paquets et librairies

install.packages("mclust")
library('mclust')
install.packages('mixtools')
library('mixtools')
install.packages('ggplot2')

# Pour le boot linux de Grenoble
library('cluster')
library('ggplot2')
library('mixtools')

### Donn?es ? analyser

# au 25 avril : on a test? exemple_image_ngris
# 25 avril et ensuite : on va tester avec exxai_slices_P26

### Fonctions ? utiliser

# Clusters avec Mclust, 
# Repr?sentation graphique de clusters pour la variable vars, sur une slice

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
  ## Notes :
  #4 centroides car on sait que le meilleur mod?le est un m?lange de 4 gaussiennes de variance 'V'ariable ##
  #break.range <- 50
  #d<-data#.vars
  #data.n <- length(d$v)
}

# Il faut faire courir cette fonction sur les jours

## D?composition en gaussiennes : algorithme EM

graph_gauss <- function(valeurs,Nclust,Col){# valeurs est le tableau contenant ici les niveaux de gris ("v") index?s par leurs coordonn?es, aui figurent dans les deux premi?res colonnes.
  names(valeurs)[1]<-"x"
  names(valeurs)[2]<-"y"
  names(valeurs)[3]<-"v"
  valeurs.mix<-normalmixEM(valeurs$v,k=Nclust)
  mu_mix=valeurs.mix$mu
  sigma_cmix=valeurs.mix$sigma
  lambda_mix=valeurs.mix$lambda
  xx=xx<-seq(from=0, to=length(valeurs$x), by=1)
  plot(x=xx,y=valeurs$v,type="o",col="black",xlab='ngris',ylab='fr?quence',main='Niveaux de gris')
  for (i in 1:Nclust){
    yy<-lambda_mix[i]*dnorm(x=xx,mean=mu_mix[i],sd=sigma_cmix[i])
    lines(yy,type="l",color=Col[i])
  }
}
# Difficile ? partir de 4 clusters : trop de tentatives

# S?rie d'histogrammes

slash="/"
slash="//" #selon les machines, on exprime diff?remment le slash dans les chemins

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






### Instructions ? ?x?cuter

## Un exemple pour faire tourner les fonctions

d_exemple <- read.table("exemple_image_ngris.txt", header=F)
d_exemple_2 <- read.table("exemple_image_ngris_2.txt", header=F)

## Un exemple pour la coloration des slices : patient 26, vars= 'ADC'.
## Dossier essai_slice_P26

color.vector <- c('darkblue','red','green','purple')
Nclust=3
#vars='ADC'
#slice='slice1"

J=c("0","3","4","08","15","22")

cl_slice1<-clust_slice(Nclust,'ADC','slice1')
cl_slice2<-clust_slice(Nclust,'BVf','slice2') # d'abord appliquer la macro ImageJ pour obttenir les tableaux xyv !




hi_slice1_ess<-hist_slice_essai(Nclust,'ADC','slice1')



dir_2<-"essai_BVf_slice2"
J_num_2<-read.table(paste(dir_2,slash,"jours.txt",sep=''))
J_char_2<-read.table(paste(dir_2,slash,"jours.txt",sep=''),as.is=rbind(rep(TRUE,length(J_num_2))))

hi_slice2_ess<-hist_slice_essai(Nclust,'BVf','slice2')
