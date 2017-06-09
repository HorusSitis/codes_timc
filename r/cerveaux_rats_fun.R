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



#---------------------------- Une fonction : segmentation d'un cerveau, une fonctionnalité, un rat ----------------------------#

cerveau_jfr <- function(day,fonc,rat){
  d.slices.day <- read.csv(paste('liste_R',rat,'_',fonc,'_J',day,'.csv',sep=''),check.names=F,header=T)
  d <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(d) <- c("x","y","z",fonc)
  for (slice.select in 1:length(d.slices.day[,1])){
    d.filename <- sprintf("%s-J%s-%s-bg-slice%i.txt",rat,day,fonc,d.slices.day[slice.select,1])
    d.increment <- read.table(d.filename,header=F,sep='\t')
    d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*d.slices.day[slice.select,1], d.increment[,3]))
    colnames(d.increment) <- c("x","y","z",fonc)
    write.table(d.increment, sprintf("%s-J%s-%s-bg-slice%i.dat",rat,day,fonc,d.slices.day[slice.select,1]), row.names=F, quote=F, sep='\t')
    d<-as.data.frame(rbind(d,d.increment))
  }
  return(d)
}


# Fonction à utiliser après avoir sélectionné le bon répertoire de travail.
# Ne pas oublier les .csv pour les slices !
adc_3d_rat <- function(rat){
  jours <- liste_jours_ADC[[rat]]
  for (j in 1:length(jours)){
    day=jours[j]
    d <- cerveau_jfr(day, 'ADC', rat)
    write.table(d, sprintf("%s-J%s-%s-bg-all.dat",rat,day,fonc), row.names=F, quote=F, sep='\t')
  }
}

# Fonction pour clusteriser un cerveau, pour une fonctionnalité quelconque. Utilise la librairie 'mclust'.
# Il faut décider du nombre de clusters voulus.
# Suppose que l'on exclut les fonctionnalités trop basses -inférieures à 10.

cluster_jfr_f10 <- function(data,cl_min,cl_max){
  #d <- read.table(sprintf("%s-J%s-%s-bg-all.dat",rat,day,fonc),header=T) # A l'avenir : chargement depuis le dossier fonctionnel_gris complet
  data.fonc <- data[,4]
  data <- data[data.fonc>10,] # cutoff à 10 de diffusion
  d.clust <- Mclust(data.fonc, G=cl_min:cl_max, modelNames="V")
  #plot(d.clust, what="classification")
  return(d.clust)
}

### Pour un suivi temporel des données tridimensionnelles : une représentation graphique par jour. Mêmes précautions à prendre aue pour adc_3d_rat.

suivi_adc_3d <- function(rat,cl_min,cl_max){
  jours <- liste_jours_ADC[[rat]]
  for (j in 1:length(jours)){
    day=jours[j]
    d <- read.table(sprintf("%s-J%s-ADC-bg-all.dat",rat,day),header=T)
    d.clust <- cluster_jfr_f10(d,cl_min,cl_max)
    # On passe aux représentations graphiques
    
    #plot(d.clust, what="classification")
    
    par(mfrow = c(2,2))
    plot(d.clust, what="BIC")
    plot(d.clust, what="classification", col=color.vector)
    
    ADC.breaks <- seq(min(d$ADC)-0.1*min(d$ADC), max(d$ADC)+0.1*max(d$ADC), length.out=100)
    d.hist <- hist(d$ADC,breaks=ADC.breaks, col='grey50')
    plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, cex=2*(1-d.clust$uncertainty)^4, xlab='x', ylab='y')
  }
}




##########################################################################################
