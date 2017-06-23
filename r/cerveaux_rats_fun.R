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



# ---------------------------- Une fonction : segmentation d'un cerveau, une fonctionnalité, un rat ---------------------------- #

cerveau_jfr <- function(jour,fonc,rat){
  
  d.filename <- sprintf("liste_R%s_%s_J%s.csv",rat,fonc,jour) #paste('liste_R',rat,'_',fonc,'_J',day,'.csv',sep='')
  day.slices <- read.csv(d.filename,check.names=F,header=T)
  d <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(d) <- c("x","y","z",fonc)
  
  for (slice.select in 1:length(day.slices[,1])){
    d.filename <- sprintf("%s-J%s-%s-bg-slice%i.txt",rat,jour,fonc,day.slices[slice.select,1])
    d.increment <- read.table(d.filename,header=F,sep='\t')
    d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*day.slices[slice.select,1], d.increment[,3]))
    colnames(d.increment) <- c("x","y","z",fonc)
    #write.table(d.increment, sprintf("%s-J%s-%s-bg-slice%i.dat",rat,jour,fonc,day.slices[slice.select,1]), row.names=F, quote=F, sep='\t')
    d <- as.data.frame(rbind(d,d.increment))
  }
  
  return(d)
}

# Fonction à utiliser après avoir sélectionné le bon répertoire de travail : fonctionnel_gris/ADC
# Ne pas oublier les .csv pour les slices !

FONC_3d_rat <- function(fonc,rat){
  jours_fonc <- liste_jfr[[fonc]] # étape intermédiaire
  jours <- jours_fonc[[rat]] # liste_jours_ADC[[rat]]
  
  for (j in 1:length(jours)){
    day=jours[j]
    d <- cerveau_jfr(day, fonc, rat)
    liste.nan <- is.na(d[,4])
    d <- d[!liste.nan,] # On retire les valeurs NaN de la dataframe
    write.table(d, sprintf("%s-J%s-%s-bg-all.dat",rat,day,fonc), row.names=F, quote=F, sep='\t')
  }
}

# ------------------- Clusterise un cerveau, pour une fonctionnalité quelconque. Utilise la librairie 'mclust'. ------------------- #

# Il faut décider du nombre de clusters voulus.
# Exclut les fonctionnalités trop basses -inférieures à 10.
# data contient la fonctionnalité

cluster_jfr_f10 <- function(data,cl_min,cl_max){
  #d <- read.table(sprintf("%s-J%s-%s-bg-all.dat",rat,day,fonc),header=T) # A l'avenir : chargement depuis le dossier fonctionnel_gris complet
  data.fonc <- data[,4]
  data <- data[data.fonc>10,] # cutoff à 10 de diffusion
  d.clust <- Mclust(data.fonc, G=cl_min:cl_max, modelNames="V")
  #plot(d.clust, what="classification")
  return(d.clust)
}

# ------------------- Pour un affichage systématique des données tridimensionnelles : une représentation graphique par jour. ------------------- #

rg_FONC_3d <- function(fonc,rat,cl_min,cl_max){
  jours_fonc <- liste_jfr[[fonc]] # étape intermédiaire
  jours <- jours_fonc[[rat]] # liste_jours_ADC[[rat]]
  for (j in 1:length(jours)){
    day=jours[j]
    d <- read.table(sprintf("%s-J%s-%s-bg-all.dat",rat,day,fonc),header=T)
    d.clust <- cluster_jfr_f10(d,cl_min,cl_max)
    d.fonc <- d[,4]
    
    # On passe aux représentations graphiques
    #plot(d.clust, what="classification")
    
    par(mfrow = c(2,2))
    plot(d.clust, what="BIC")
    plot(d.clust, what="classification", col=color.vector)
    
    FONC.breaks <- seq(min(d.fonc)-0.1*min(d.fonc), max(d.fonc)+0.1*max(d.fonc), length.out=100)
    d.hist <- hist(d.fonc,breaks=FONC.breaks, col='grey50',main=paste("Histogram of",fonc))
    #plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, cex=2*(1-d.clust$uncertainty)^4, xlab='x', ylab='y',main=paste("Cerveau ",rat,", J",day))
    plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, xlab='x', ylab='y',main=paste("Cerveau ",rat,", J",day))
  }
}



# ------------------- Extraction de clusters pour une segmentation. Ouvrir dans le répertoire correspondant à la fonctionnalité concernée ------------------- #

seg_cl_FONC <- function(day,fonc,rat,clust,hemi){
  d <- read.table(sprintf("%s-J%s-%s-bg-all.dat",rat,day,fonc),header=T)
  d.clust <- cluster_jfr_f10(d,2,5)
  
  l <- length(d[,4])
  d.label <- rep(0,l)
  d.label <- ifelse(d.clust$classification==1,1,d.label)
  d.label <- ifelse(d$x>hemi,2,d.label)
  
  d.seg <- as.data.frame(cbind(d$x,d$y,d$z,d.label),header=T)#d[d.clust$classification==clust,]
  colnames(d.seg) <- c("x","y","z","Label")
  return(d.seg)
}

# ------------------- Représente les niveaux de gris pour chaque fonctionnalité, sur un cerveau segmenté ------------------- #

gr_ngris_seg <- function(rat,jour){
  cerveau_seg <- read.table(sprintf("isch3d-fonc-%s-J%s.dat",rat,jour),header=T)
  
  for (ind_fonc in 1:7){
    fonc <- liste_fonc[[ind_fonc]]
    cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
    
    cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]
    cerveau_hem <- cerveau_fonc[cerveau_seg$Label==2,]
    
    entier <- cerveau_fonc[,4]
    isch <- cerveau_isch[,4]
    hem_sain <- cerveau_hem[,4]
    n <- length(entier)
    ni <- length(isch)
    ns <- length(hem_sain)
    
    dst <- density(entier)
    dsti <- density(isch)
    dsts <- density(hem_sain)
    
    plot.new()
    par(lend="butt")
    plot(dst$x,dst$y,type="n",main=sprintf("Rat %s jour %s %s",rat,jour,fonc))
    lines(dsti$x, ni/n*dsti$y, lwd = 2, col = "darkred")
    lines(dsts$x, ns/n*dsts$y, lwd = 2, lty = 2, col = "darkblue")
    lines(dst$x, dst$y, lwd = 3, col="gray70")
    
    legend("topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain","Cerveau entier"),
           col = c("darkred","darkblue","gray70"),
           lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
  }
}




# ------------------- Pour un suivi temporel des données bidimensionnelles : une représentation graphique par jour. ------------------- #


  



##########################################################################################
