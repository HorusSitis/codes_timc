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
  
  d.filename <- sprintf("%s/liste_R%s_%s_J%s.csv",fonc,rat,fonc,jour)
  day.slices <- read.csv(d.filename,check.names=F,header=T)
  d <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(d) <- c("x","y","z",fonc,"Slice")
  
  for (slice.select in 1:length(day.slices[,1])){
    d.filename <- sprintf("%s/%s-J%s-%s-bg-slice%i.txt",fonc,rat,jour,fonc,day.slices[slice.select,1])
    d.increment <- read.table(d.filename,header=F,sep='\t')
    d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*day.slices[slice.select,1], d.increment[,3], day.slices[slice.select,1]))
    colnames(d.increment) <- c("x","y","z",fonc,"Slice")
    d <- as.data.frame(rbind(d,d.increment))
  }
  
  return(d)
}

# ---------------------------- Une fonction : segmentation d'un cerveau, une fonctionnalité, un rat ---------------------------- #

tranches_fr_ker <- function(jour,fonc,rat){

  d.filename <- sprintf("%s/suivi_temp_R%s_%s.csv",fonc,rat,fonc)
  day.slices <- read.csv(d.filename,check.names=F,header=T)
  d <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(d) <- c("x","y","z",fonc,"Slice")
  
  for (slice.select in 1:length(day.slices[,1])){
    d.filename <- sprintf("%s/%s-J%s-%s-ker-slice%i.txt",fonc,rat,jour,fonc,day.slices[slice.select,1])
    d.increment <- read.table(d.filename,header=F,sep='\t')
    d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*day.slices[slice.select,1], d.increment[,3], day.slices[slice.select,1]))
    colnames(d.increment) <- c("x","y","z",fonc,"Slice")
    d <- as.data.frame(rbind(d,d.increment))
  }
  
  return(d)
}

# ------------------- Boucle pour enregistrer les cerveaux : coordonnées tridimensionnelles et nom des slices pour les niveaux de gris  ------------------- #

# Fonction à utiliser après avoir sélectionné le bon répertoire de travail : fonctionnel_gris
# Ne pas oublier les .csv pour les slices !

FONC_3d_rat <- function(fonc,rat,val){
  jours_fonc <- liste_jfr[[fonc]]
  jours <- jours_fonc[[rat]]
  
  for (j in 1:length(jours)){
    day=jours[j]
    d <- cerveau_jfr(day, fonc, rat)
    if (val=='sans'){
      liste.nan <- is.na(d[,4])
      d <- d[!liste.nan,] # On retire les valeurs NaN de la dataframe si demandé
    }
    write.table(d, sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,day,fonc), row.names=F, quote=F, sep='\t')
  }
}

# ------------------- Même chose, avec les tranches dénoyautées suivables temporellement  ------------------- #

# Fonction à utiliser après avoir sélectionné le bon répertoire de travail : fonctionnel_gris
# Ne pas oublier les .csv pour les slices !

dker_3d_rat <- function(fonc,rat,val){# on garde les valeurs manquantes
  jours_fonc <- liste_jfr[[fonc]]
  jours <- jours_fonc[[rat]]
  
  for (j in 1:length(jours)){
    day=jours[j]
    d <- tranches_fr_ker(day, fonc, rat)

    write.table(d, sprintf("%s/%s-J%s-%s-kertemp-all.dat",fonc,rat,day,fonc), row.names=F, quote=F, sep='\t')
  }
}

# ------------------- Clusterise un cerveau, pour une fonctionnalité quelconque. Utilise la librairie 'mclust'. ------------------- #

# Il faut décider du nombre de clusters voulus.
# Exclut les fonctionnalités trop basses -inférieures à 10.
# data contient la fonctionnalité

cluster_jfr_fmin <- function(data,cl_min,cl_max,min){# on adapte le rabotage à la fonctionnalité.
  data.fonc <- data[,4]
  data <- data[data.fonc>min,] # cutoff à 10 de diffusion
  d.clust <- Mclust(data.fonc, G=cl_min:cl_max, modelNames="V")
  #plot(d.clust, what="classification")
  return(d.clust)
}

# ------------------- Pour un affichage systématique des données tridimensionnelles : une représentation graphique par jour. ------------------- #

rg_FONC_3d <- function(dim,fonc,rat,cl_bounds,opt){# option : clusterisation sur les tranches dénoyautées suivables temporellement
  
  jours_fonc <- liste_jfr[[fonc]]
  jours <- jours_fonc[[rat]]
  min_fonc <- liste_min_fonc[[fonc]]
  
  cl_min <- cl_bounds[1]
  cl_max <- cl_bounds[2]
  
  for (j in 1:length(jours)){
    day=jours[j]
    if (opt=='ker'){
      d <- read.table(sprintf("%s/%s-J%s-%s-kertemp-all.dat",fonc,rat,day,fonc),header=T)
    }
    else if (opt=='dark'){
      d <- read.table(sprintf("%s/%s-J%s-%s-dark-all.dat",fonc,rat,day,fonc),header=T)
    }
    else{
      d <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,day,fonc),header=T)
    }

    liste.nan <- is.na(d[,4])
    d <- d[!liste.nan,] # On retire les valeurs NaN de la dataframe, pas besoin de réassembler ensuite

    d.clust <- cluster_jfr_fmin(d,cl_min,cl_max,min_fonc)
    d.fonc <- d[,4]
    
    # On passe aux représentations graphiques

    #par(mfrow = c(2,2))
    m <- rbind(c(1,2,3),c(4,4,4))
    layout(m, heights = c(1,3), width=c(1,1,1), respect = TRUE)

    plot(d.clust, what="BIC")
    plot(d.clust, what="classification", col=color.vector)
    
    FONC.breaks <- seq(min(d.fonc)-0.1*min(d.fonc), max(d.fonc)+0.1*max(d.fonc), length.out=100)
    d.hist <- hist(d.fonc,breaks=FONC.breaks, col='grey50',main=paste("Histogram of",fonc))
    
    if (dim==2){ # Cerveaux en deux ou trois dimensions
      plot(d$x, d$y, col=color.vector[d.clust$classification],
           pch=20, cex=2*(1-d.clust$uncertainty)^4,
           xlab='x', ylab='y',
           main=paste("Cerveau ",rat,", J",day))
      #plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, xlab='x', ylab='y',main=paste("Cerveau ",rat,", J",day))
    }
    else{
      cerveau.clust <- cbind(d,d.clust$classification)
      colnames(cerveau.clust) <- c("x","y","z",'ADC',"Slice","clust")
      scatterplot3d(cerveau.clust$x,
                    cerveau.clust$y,
                    cerveau.clust$Slice,
                    #cerveau.clust$z, 
                    color= color.vector[cerveau.clust$clust],
                    pch=20,
                    #cex=2*(1-d.clust$uncertainty)^4,
                    xlab='x',
                    ylab='y',
                    zlab='Slice',#'z',
                    lab.z=1+cerveau.clust$Slice[length(cerveau.clust$Slice)]-cerveau.clust$Slice[1],
                    main=paste("Cerveau ",rat,", J",day)
      )
    }
    # Fin de la représentation graphique du jour
  }
}

# ------------------- Comparaison graphique entre clusterisation et segmentation ------------------- #

seg_clust_3d <- function(fonc,rat,cl_bounds,cl_seg,hemi){

  liste_jr <- liste_jfr[[fonc]]
  liste_j <- liste_jr[[rat]]
  min_fonc <- liste_min_fonc[[fonc]]
  
  cl_min <- cl_bounds[1]
  cl_max <- cl_bounds[2]

  for (jour in liste_j){

    # Importation des données : cerveau tridimensionnel et niveaux de gris de la fonctionnalité, tranches indexées.
    t <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)

    # Mise à part des valeurs manquantes.
    liste.nan <- is.na(t[,4])
    d <- t[!liste.nan,]
    e <- t[liste.nan,]
    
    # Clusterisation
    d.clust <- cluster_jfr_fmin(d,cl_min,cl_max,min_fonc)
    d.fonc <- d[,4]

    # Etiquetage des valeurs manquantes : pour utiliser color_seg.
    m <- length(e[,4])
    e.label <- rep(3,m)
    e.seg <- as.data.frame(cbind(e$x,e$y,e$z,e.label,e$Slice),header=T)
    colnames(e.seg) <- c("x","y","z","Label","Slice")
    
    # Segmentation à l'aide des clusters sur lesquels on effectue le test.
    clusters <- cl_seg[[jour]]
    l <- length(d[,4])
    d.label <- rep(0,l)
    for (cl in clusters){
      d.label <- ifelse(d.clust$classification==cl,1,d.label)
    }
    d.label <- ifelse(hemi[1]*d$x+hemi[2]>d$y,2,d.label)
    d.seg <- as.data.frame(cbind(d$x,d$y,d$z,d.label,d$Slice),header=T)
    colnames(d.seg) <- c("x","y","z","Label","Slice")
    
    # Réassemblage de la dataframe : étiquetage et inclusion des valeurs manquantes
    d.seg <- as.data.frame(rbind(d.seg,e.seg))
    # Seulement utilisé dans la dernière fenêtre graphique
    
    # On passe aux représentations graphiques

    #get( getOption( "device" ) )()
    plot.new()
    split.screen( figs = c( 2, 1 ) )
    split.screen( figs = c( 1, 3 ), screen = 1 )
    split.screen( figs = c( 1, 2 ), screen = 2 )

    screen( 3 )
    plot(d.clust, what="BIC")
    screen( 4 )
    plot(d.clust, what="classification", col=color.vector)
    
    FONC.breaks <- seq(min(d.fonc)-0.1*min(d.fonc), max(d.fonc)+0.1*max(d.fonc), length.out=100)
    screen( 5 )
    d.hist <- hist(d.fonc,breaks=FONC.breaks, col='grey50',main=paste("Histogram of",fonc))
    
    screen( 6 )
    cerveau.clust <- cbind(d,d.clust$classification)
    colnames(cerveau.clust) <- c("x","y","z",fonc,"Slice","clust")
    scatterplot3d(cerveau.clust$x,
                  cerveau.clust$y,
                  cerveau.clust$Slice, 
                  color= color.vector[cerveau.clust$clust],
                  pch=20,
                  #cex=2*(1-d.clust$uncertainty)^4,
                  xlab='x',
                  ylab='y',
                  zlab='Slice',
                  lab.z=1+cerveau.clust$Slice[length(cerveau.clust$Slice)]-cerveau.clust$Slice[1],
                  main=paste("Cerveau ",rat, ", J",jour)
    )
    screen( 7 )
    scatterplot3d(d.seg$x,
                  d.seg$y,
                  d.seg$z,
                  color=color_seg[1+d.seg$Label],
                  pch=20,
                  xlab='x',
                  ylab='y',
                  zlab='z',
                  main="Segmentation"
    )
    #close.screen( all = TRUE )
  }
}

# ------------------- Comparaison graphique : tranches d'intérêt clusterisées au sein, ou non, de cerveaux entiers. ------------------- #

comp_2vs3d_clust <- function(fonc,rat,cl_bounds,num_tranche){
  
  cl_min <- cl_bounds[1]
  cl_max <- cl_bounds[2]
  
  min_fonc <- liste_min_fonc[[fonc]]
  
  liste_jr <- liste_jfr[[fonc]]
  liste_j <- liste_jr[[rat]]

  for (jour in liste_j){
    
    cerveau <- read.table(sprintf('%s/%s-J%s-%s-bg-all.dat',fonc,rat,jour,fonc),header=T)
    c.nan <- is.na(cerveau[,4])
    cerveau <- cerveau[!c.nan,]
    
    #if (rat=="30"){# on retire les structures visibles à l'oeil, pour garder la substance blanche
    #  c.clust.v0 <- cluster_jfr_fmin(cerveau,cl_min,cl_max,min_fonc)
    #  cerveau <- cerveau[c.clust.v0$uncertainty>0.3,]
    #}

    tranche <- cerveau[cerveau$Slice==num_tranche,]
    #t.nan <- is.na(tranche[,4])
    #tranche <- tranche[!t.nan,]

    c.clust <- cluster_jfr_fmin(cerveau,cl_min,cl_max,min_fonc)
    
    c.lb.clust <- as.data.frame(cbind(cerveau,c.clust$classification,c.clust$uncertainty),header=T)
    colnames(c.lb.clust) <- c("x","y","z",fonc,"Slice","ClLabel","UncLabel")
    tc <- c.lb.clust[c.lb.clust$Slice==num_tranche,] # tranche clusterisée extraite
    colnames(tc) <- c("x","y","z",fonc,"Slice","ClLabel","UncLabel")

    t.clust <- cluster_jfr_fmin(tranche,cl_min,cl_max,min_fonc) # tranche clusterisée seule
    
    #plot.new()
    m <- rbind(c(1,4),c(2,0),c(3,5))
    layout(m, heights = c(1,1,1), width=c(1,1), respect = TRUE)
    
    # ------------------------------ Cerveau entier projeté, tranche d'intérêt, classification globale ------------------------------ #
    
    plot(tc$x, tc$y, col=color.vector[tc$ClLabel], # tranche extraite du cerveau entier clusterisé
         pch=20,
         cex=2*(1-tc$UncLabel)^4,
         xlab='x', ylab='y',
         main=sprintf("Tranche %i clusterisation induite",num_tranche),
         cex.main = 1.5,
         sub = sprintf("%s : Rat %s, J%s",fonc,rat,jour),
         col.sub = 'red',
         cex.sub=1.5
         )
    plot(cerveau$x, cerveau$y, col=color.vector[c.clust$classification], # cerveau entier clusterisé
         pch=20,
         cex=2*(1-c.clust$uncertainty)^4,
         xlab='x', ylab='y',
         main = sprintf("Cerveau entier clusterisé",num_tranche),
         cex.main = 1.5
         )
    plot(c.clust, what="classification", col=color.vector)
    title(sub = "Cerveau entier",
          col.sub = 'red',
          #cex.main = 1.5,
          cex.sub=1.5
          )
    
    # ------------------------------ Classification 2d pour la tranche d'intérêt ------------------------------ #
    
    plot(tranche$x, tranche$y, col=color.vector[t.clust$classification], # tranche seule
         pch=20,
         cex=2*(1-t.clust$uncertainty)^4,
         xlab='x', ylab='y',
         main=sprintf("Tranche %i clusterisée seule",num_tranche),
         cex.main = 1.5,
         sub = sprintf("%s : Rat %s, J%s",fonc,rat,jour),
         col.sub = 'red',
         cex.sub=1.5
         )
    plot(t.clust, what="classification", col=color.vector)
    title(sub = sprintf("Tranche %i clusterisée seule",num_tranche),
          col.sub = 'red',
          #cex.main = 1.5,
          cex.sub=1.5
          )
  }
}

# ------------------- Extraction de clusters pour une segmentation. Ouvrir dans le répertoire fonctionnel_gris ------------------- #

seg_cl_FONC <- function(jour,fonc,rat,clusters,hemi){# clust : vecteur des clusters utilisables pour la segmentation, estimés à l'oeil à partir de la fonction rg_FONC_3d.

  # Importation des données : cerveau tridimensionnel, tranches indexées
  t <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
    
  # Extraction des valeurs manquantes
  liste.nan <- is.na(t[,4])
  d <- t[!liste.nan,]
  e <- t[liste.nan,]
  
  # Clusterisation
  min_fonc <- liste_min_fonc[[fonc]]
  d.clust <- cluster_jfr_fmin(d,2,5,min_fonc)
  
  # Etiquetage pour la zone ischémiée : on utilise les clusters repérés à l'étape 2.
  l <- length(d[,4])
  d.label <- rep(0,l)
  for (cl in clusters){
    d.label <- ifelse(d.clust$classification==cl,1,d.label)
  }
  d.label <- ifelse(hemi[1]*d$x+hemi[2]>d$y,2,d.label)
  d.seg <- as.data.frame(cbind(d$x,d$y,d$z,d.label,d$Slice))#,header=T)#d[d.clust$classification==clust,]
  colnames(d.seg) <- c("x","y","z","Label","Slice")
  
  # Etiquetage des valeurs manquantes
  m <- length(e[,4])
  e.label <- rep(-1,m)
  e.seg <- as.data.frame(cbind(e$x,e$y,e$z,e.label,e$Slice))#,header=T)
  colnames(e.seg) <- c("x","y","z","Label","Slice")
  
  # Réassemblage de la dataframe : étiquetage et inclusion des valeurs manquantes
  d.seg <- as.data.frame(rbind(d.seg,e.seg))
  return(d.seg)
}

# ------------------- Enregistrement systématique des segmentations produites avec seg_cl_FONC. Utilise une liste de clusters établie à l'étape 2. ------------------- #

record_seg_cl <- function(rat, liste_clust, hemi){
  
  for (fonc in liste_fonc){
    liste_jr <- liste_jfr[[fonc]]
    liste_j <- liste_jr[[rat]]
    
    liste_clust_fonc <- liste_clust[[fonc]] # liste des familles de clusters permettant la segmentation, par jour.
    
    for (jour in liste_j){
      d_seg <- seg_cl_FONC(jour,fonc,rat,liste_clust_fonc[[jour]],hemi)
      write.table(d_seg, sprintf("%s/isch3d-%s-%s-J%s.dat",fonc,fonc,rat,jour), row.names=F, quote=F, sep='\t')
    }
    
  }
  
}

# ------------------- Enregistrement du volume, lésé en ADC au jour 00, examiné au jour et à la fonctionnalité en cours, avec les tranches du vecteur tr ------------------- #

vol_lesADC00 <- function(rat,tr){
  repertoires <- list('ADC'="fonctionnel_gris",# répertoire du rat courant
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  liste_jr <- liste_jfr[['ADC']]
  jours <- liste_jr[[rat]]
  for (jour in jours){
    #d.filename <- sprintf("%s/liste_R%s_%s_J00.csv",'ADC',rat,'ADC')
    #dr.filename <- paste(repertoires[[fonc]],"/",d.filename)
    day.slices <- tr
    d <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(d) <- c("x","y","z",'ADC',"Slice")
    
    for (slice in day.slices){
      d.filename <- sprintf("%s/%s-J%s-%s-dark-slice%i.txt",'ADC',rat,jour,'ADC',slice)
      dr.filename <- paste(repertoires[['ADC']],"/",d.filename,sep='') # on va chercher les niveaux de gris de la zone segmentée ADC dane le bon répertoire
      d.increment <- read.table(dr.filename,header=T,sep='\t')
      d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*slice, d.increment[,3],slice))
      colnames(d.increment) <- c("x","y","z",'ADC',"Slice")
      d <- as.data.frame(rbind(d,d.increment))
    }
    w.filename <- sprintf("%s/%s-J%s-%s-dark-all.dat",'ADC',rat,jour,'ADC')
    wr.filename <- paste(repertoires[['ADC']],"/",w.filename,sep='')# on met le volume segmenté dans le bon répertoire
    write.table(d, wr.filename, row.names=F, quote=F, sep='\t')
  }
}

# ------------------- Enregistrement du volume, lésé en CBF au jour 00, examiné au jour et à la fonctionnalité en cours, avec les tranches du vecteur tr. Répertoire benjamin_antoine_labo. ------------------- #

vol_lesCBF00 <- function(rat,tr){
  repertoires <- list('ADC'="fonctionnel_gris",# répertoire du rat courant
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  fonc <- 'CBF'
  liste_jr <- liste_jfr[['CBF']]
  jours <- liste_jr[[rat]]
  for (jour in c("00")){#jours){
    day.slices <- tr
    d <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(d) <- c("x","y","z",fonc,"Slice")
    
    for (slice in day.slices){
      # Lecture du fichier texte des valeurs de la modalité -> niveaux de gris pour ImageJ
      d.filename <- sprintf("%s/%s-J%s-%s-dark-slice%i.txt",'CBF',rat,jour,'CBF',slice)
      dr.filename <- paste("R",rat,"/",repertoires[[fonc]],"/",d.filename,sep='') # on va chercher les niveaux de gris de la zone segmentée ADC dane le bon répertoire
      d.increment <- read.table(dr.filename,header=T,sep='\t')
      
      # Remplissage de la dataframe à enregistrer
      d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*slice, d.increment[,3],slice))
      colnames(d.increment) <- c("x","y","z",fonc,"Slice")
      d <- as.data.frame(rbind(d,d.increment))
      print(sprintf("Rat %s, tranche%s, taille%i",rat,slice,length(d.increment[,3])))
    }
    w.filename <- sprintf("%s/%s-J%s-%s-dark-all.dat",fonc,rat,jour,fonc)
    wr.filename <- paste("R",rat,"/",repertoires[[fonc]],"/",w.filename,sep='')# on met le volume segmenté dans le bon répertoire
    write.table(d, wr.filename, row.names=F, quote=F, sep='\t')
  }
}

# ------------------- Représente graphiquement les volumes lésés, enregistrés avec ... , clusterisés. Répertoire benjamin_antoine_labo.  ------------------- #

# Option 1 : segmentation utilisée.
# Option 2 : sortie. pdf si 'pdf', affichage dans RStudio sinon.

# Attention : faire commencer liste_fonc par la modalité utilisée pour la segmentation.

comp_clust_vol00 <- function(rat,liste_s_slice,cl,lm_fonc,opt_1,opt_2){
  repertoires <- list('ADC'="fonctionnel_gris",# répertoire du rat courant
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  
  cl_min <- cl[1]
  cl_max <- cl[length(cl)]
  
  lsauv_seg_clust <- list("00"=list("vol"='',"val"='',"clust"=''),
                               "03"=list("vol"='',"val"='',"clust"=''),
                               "08"=list("vol"='',"val"='',"clust"=''),
                               "15"=list("vol"='',"val"='',"clust"=''),
                               "22"=list("vol"='',"val"='',"clust"='')
                               )#liste des adc_cerveau_les, adc_les et adc.l.clust
  
  for (fonc in liste_fonc){
    segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
    if (opt_1=='ADCdark00'){
      fonc_seg <- 'ADC'
    }
    else if (opt_1=='brightAnat00'){
      fonc_seg <- 'Anat'
    }
    else if (opt_1=='CBFdark00'){
      fonc_seg <- 'CBF'
    } # la segmentation avec fonc_seg sera celle utilisée
    
    tranches <- liste_s_slice[[fonc_seg]]
    
    subtitle <- sprintf('Rat %s, tranche(s) ',rat)
    for (tr in tranches){
      subtitle <- paste(subtitle,'-',tr)
    }
    
    liste_jr <- liste_jfr[[fonc]]
    jours <- liste_jr[[rat]]
    
    for (jour in jours){
      cerveau_fonc <- read.table(sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,jour,fonc,'bg'),header=T)
      cerveau_seg <- read.table(sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark'),header=T)
      
      cerveau_les <- cerveau_fonc # on va délimiter une enveloppe rectangulaire de la partie lésée sur l'image observée, fonctionnalité et jour courants
      
      if (fonc=='Anat'){# anatomique : trop de valeurs manquantes qui allongent la boucle
        list.nan <- is.na(cerveau_les[,4])
        cerveau_les <- cerveau_les[!list.nan,]
      }
      
      l <- length(cerveau_les[,4])
      liste_tr <- rep(FALSE,l)
      m <- length(cerveau_seg[,4])
      for (i in c(1:m)){
        ix <- cerveau_seg[i,1]
        iy <- cerveau_seg[i,2]
        sli <- cerveau_seg[i,5]
        liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
      }
      cerveau_les <- cerveau_les[liste_tr,]
      # on crée la liste de niveaux de gris exploitable par density() : zone lésée
      les <- cerveau_les[,4]
      liste.nan <- is.na(les)
      cerveau_les <- cerveau_les[!liste.nan,]
      les <- les[!liste.nan] # on retire les valeurs manquantes
      
      l.clust <- cluster_jfr_fmin(cerveau_les,cl_min,cl_max,lm_fonc[[fonc]])
      
      if (fonc==fonc_seg){
        #adc_cerveau_les <- cerveau_les
        #adc_les <- les
        #jouradc.l.clust <- l.clust
        
        lsauv_seg_clust[[jour]] <- list("vol"=cerveau_les,"val"=les,"clust"=l.clust)
      }
      
      if (opt_2=='pdf'){
        # on verra
      }# une lésion imprimée
      else{# On passe à la représentation graphique
        
        if (fonc==fonc_seg){
          m <- rbind(c(1,2,3),c(4,4,4))
          layout(m, heights = c(1,3), width=c(1,1,1), respect = TRUE)
          
          plot(l.clust, what="BIC")
          plot(l.clust, what="classification", col=color.vector)
          
          FONC.breaks <- seq(min(les)-0.1*min(les), max(les)+0.1*max(les), length.out=100)
          l.hist <- hist(les,breaks=FONC.breaks, col='grey50',main=paste("Histogram of",fonc))
          
          les.clust <- cbind(cerveau_les,l.clust$classification)
          colnames(les.clust) <- c("x","y","z",fonc,"Slice","clust")
          scatterplot3d(les.clust$x,
                        les.clust$y,
                        les.clust$Slice,
                        #cerveau.clust$z, 
                        color= color.vector[les.clust$clust],
                        pch=20,
                        cex.symbols=5*(1-l.clust$uncertainty)^4,
                        xlab='x',
                        ylab='y',
                        zlab='Slice',#'z',
                        lab.z=1+les.clust$Slice[length(les.clust$Slice)]-les.clust$Slice[1],
                        main=paste("Cerveau ",rat,", J",jour)
          )
          title(sprintf("Zone lésée clusterisée : rat %s, modalité %s, jour %s",rat,fonc,jour),outer=TRUE)
        }# modalité à traiter en premier : attention à liste_fonc
        else{
          liste_sauv <- lsauv_seg_clust[[jour]]
          
          seg_cerveau_les <-liste_sauv[["vol"]]
          seg_les <- liste_sauv[["val"]]
          seg.l.clust <- liste_sauv[["clust"]]
          
          #get( getOption( "device" ) )()
          #plot.new()
          #split.screen( figs = c( 2, 1 ) )
          #split.screen( figs = c( 1, 3 ), screen = 1 )
          #split.screen( figs = c( 1, 2 ), screen = 2 )
          layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2, 6,byrow = TRUE),
                 widths=c(1,1,1,1,1,1),
                 heights=c(1,2)
                 )
          
          #screen( 3 )
          plot(l.clust, what="BIC")
          #screen( 4 )
          plot(l.clust, what="classification", col=color.vector)
          
          FONC.breaks <- seq(min(les)-0.1*min(les), max(les)+0.1*max(les), length.out=100)
          #screen( 5 )
          l.hist <- hist(les,breaks=FONC.breaks, col='grey50',main=paste("Histogram of",fonc))
          
          #screen( 6 )
          les.clust <- cbind(cerveau_les,l.clust$classification)
          colnames(les.clust) <- c("x","y","z",fonc,"Slice","clust")
          scatterplot3d(les.clust$x,
                        les.clust$y,
                        les.clust$Slice,
                        #cerveau.clust$z, 
                        color= color.vector[les.clust$clust],
                        pch=20,
                        cex.symbols=2.5*(1-l.clust$uncertainty)^3,
                        xlab='x',
                        ylab='y',
                        zlab='Slice',#'z',
                        lab.z=1+les.clust$Slice[length(les.clust$Slice)]-les.clust$Slice[1],
                        main=sprintf("Cerveau %s, J%s, %s",rat,jour,fonc)
          )
          
          #screen( 7 )
          seg.les.clust <- cbind(seg_cerveau_les,seg.l.clust$classification)
          colnames(seg.les.clust) <- c("x","y","z",fonc_seg,"Slice","clust")
          scatterplot3d(seg.les.clust$x,
                        seg.les.clust$y,
                        seg.les.clust$Slice,
                        #cerveau.clust$z, 
                        color= color.vector[seg.les.clust$clust],
                        pch=20,
                        cex.symbols=2.5*(1-seg.l.clust$uncertainty)^3,
                        xlab='x',
                        ylab='y',
                        zlab='Slice',#'z',
                        lab.z=1+seg.les.clust$Slice[length(seg.les.clust$Slice)]-seg.les.clust$Slice[1],
                        main=paste("Comparaison ",rat,", J",jour,fonc_seg)
          )
          title(sprintf("Zone lésée clusterisée : rat %s, modalité %s, jour %s",rat,fonc,jour),outer=TRUE)
          #close.screen( all = TRUE )
        }# comparaison entre la fonctionnalité courante et la modalité de segmentation
      }# lésion affichée
    }# cerveaux terminés, fonctionnalité courante
  }# fonctionnalités toutes traitées
}

# ------------------- Représente les niveaux de gris pour chaque fonctionnalité, sur un cerveau segmenté avec FONC ------------------- #

gr_ngris_seg <- function(jour,FONC,rat){
  
  cerveau_seg <- read.table(sprintf("%s/isch3d-%s-%s-J%s.dat",FONC,FONC,rat,jour),header=T)
  
  liste_F <- list()

  for (fonc in liste_fonc){# On parcourt en largeur d'abord l'arborescence de la base de données pour savoir quelles fonctionnalités sont disponibles pour le jour courant.
    liste_jf <- liste_jfr[[fonc]]
    liste_j <- liste_jf[[rat]]
    if (any(liste_j==jour)){
      liste_F <- cbind(liste_F,list(fonc))# Liste des fonctionnalités disponibles
    }
  }

  plot.new()
  par(mfrow=c(4,2))
  
  for (fonc in liste_F){# Boucle pour représenter conjointement les densités pour toutes les fonctionnalités
    cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
    
    cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]
    cerveau_hem <- cerveau_fonc[cerveau_seg$Label==2,]
    
    entier <- cerveau_fonc[,4]
    liste.nan <- is.na(entier)
    entier <- entier[!liste.nan,] # on retire les valeurs manquantes de la fonctionnalité pour évaluer sa densité
    
    isch <- cerveau_isch[,4]
    hem_sain <- cerveau_hem[,4]
    
    n <- 1#length(entier)
    ni <- 1#length(isch)
    ns <- 1#length(hem_sain)
    
    dst <- density(entier)
    dsti <- density(isch)
    dsts <- density(hem_sain)
    
    #plot.new()
    #par(lend="butt")
    plot(dst$x,dst$y,type="n",main=sprintf("Rat %s jour %s %s",rat,jour,fonc))
    lines(dsti$x, ni/n*dsti$y, lwd = 2, col = "darkred")
    lines(dsts$x, ns/n*dsts$y, lwd = 2, lty = 2, col = "darkblue")
    lines(dst$x, dst$y, lwd = 3, col="gray70")
    
    legend("topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain","Cerveau entier"),
           col = c("darkred","darkblue","gray70"),
           lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
  }
  return(liste_F)# vérification
}

# ------------------- Suivi temporel de tranches clusterisées ou segmentées. Répertoire fonctionnel_gris. Une représentation graphique par jour. ------------------- #

suivi_temp_fonc <- function(rat, fonc){#, class){# représentations graphiques tridimensionnelles
  
  filename <- sprintf("%s/suivi_temp_R%s_%s.csv",fonc,rat,fonc) # tranches du suivi temporel, indices de la première boucle.
  slices_temp <- read.csv(filename,check.names=F,header=T)
  slices <- slices_temp[,1]
  
  jours_fonc <- liste_jfr[[fonc]]
  jours <- jours_fonc[[rat]]
  
  for (num.slice in slices){
    plot.new()
    par(mfrow=c(2,3))
    #if (class == "clust"){
      for (jour in jours){
        # On représente graphiquement, jour après jour, les images de la tranche courante.
        #jour <- jours[j]
        name <- sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc)
        d <- read.table(name, header=T, check.names = F)
        # On clusterise avant d'extraire la tranche qui nous intéresse
        d.clust <- cluster_jfr_f10(d,3,5)
        d <- as.data.frame(cbind(d,d.clust$classification))
        d <- d[d$Slice==num.slice,]
        
        plot(d$x, d$y,
             col=color.vector[d[,6]],
             pch=20,
             cex=2*(1-d.clust$uncertainty)^4, 
             xlab='x', ylab='y',
             main=sprintf("Rat %s, tranche %i J%s",rat,num.slice,jour)
        )
      }# images de la tranche courante toutes faites
      #plot(d.clust, what="classification", col=color.vector) # légende commune ... pertinence ?
    #}
  }# suivi des tranches fait
}

# ------------------- Suivi temporel des densités de niveaux de gris sur un cerveau segmenté, par fonctionnalité. Répertoire fonctionnel_gris du rat concerné. ------------------- #

# Option 1 : chaîne. Segmentation choisie.

# Option 2 : liste.
# Indique que les statistiques seront effectuées sur des cerveaux entiers, selon un suivi temporel de volume maximal ou sur une seule tranche.

# Option 3 : chaîne, sortie. pdf si 'pdf', affichage dans RStudio sinon.


dgris_temp_fonc <- function(rat,hemi,opt_1,opt_2,opt_3){#liste_s_slice,
  repertoires <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")

  if (any(opt_2=='cer')){# suivi sur le cerveau entier, boucle sur les jours existant pour chaque fonctionnalité.
    for (fonc in liste_fonc){
      #segtitle <- '' # indique si nécessaire la fonctionnalité utilisée pour la segmentation
      fonc_seg <- fonc # la segmentation avec fonc_seg sera celle utilisée
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]# plus tard --------#
      
      # Représentation graphique : fonctionnalité courante
      plot.new()
      par(mfrow=c(2,3),cex.main=1.7, cex.sub=1.2,col.main="black", col.sub="red")
      
      for (jour in jours){# une fenêtre pour la fonctionnalité courante
        
        cerveau_seg <- read.table(sprintf('%s/isch3d-%s-%s-J%s.dat',fonc_seg,fonc_seg,rat,jour),header=T)#,checknames=F)
        cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
        cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]
        
        if (jour=="00"){
          cerveau_hem <- cerveau_fonc[cerveau_seg$Label==2,]
          hem_sain <- cerveau_hem[,4]
          liste.nan <- is.na(hem_sain)
          hem_sain <- hem_sain[!liste.nan] # on retire les valeurs manquantes de la fonctionnalité pour évaluer sa densité
        }
        # on définit hem_sain au jour 00, on ne le modifie plus par la suite
        
        entier <- cerveau_fonc[,4]
        liste.nan <- is.na(entier)
        entier <- entier[!liste.nan] # on retire les valeurs manquantes
        
        isch <- cerveau_isch[,4]
        liste.nan <- is.na(isch)
        isch <- isch[!liste.nan] # on retire les valeurs manquantes

        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        n <- 1#length(entier)
        ni <- 1#length(isch)
        ns <- 1#length(hem_sain)
        
        dst <- density(entier)
        if (length(isch)!=0){
          dsti <- density(isch)
        }
        dsts <- density(hem_sain)
        
        #plot.new()
        #par(lend="butt")
        title <- sprintf("Rat %s jour %s %s",rat,jour,fonc)
        #title <- paste(title,segtitle)
        
        if (length(isch)!=0){
          plot(dst$x,dst$y,type="n",main=title,sub='Cerveau entier')
          lines(dsti$x, ni/n*dsti$y, lwd = 2, col = "darkred")
          lines(dsts$x, ns/n*dsts$y, lwd = 2, lty = 2, col = "darkblue")
          lines(dst$x, dst$y, lwd = 3, col="gray70")
          
          legend("topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain J00","Cerveau entier"),
                 col = c("darkred","darkblue","gray70"),
                 lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
        }
      }
      # sub-plot fait
    }
    # fonctionnalité vue
  }# option 2 codée : plus grand volume disponible pour chaque jour
  else if (length(opt_2)==1){# zones segmentées à l'oeil : sombres T2 ou ? lumineuses T1.
    for (fonc in liste_fonc){
      
      segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
      if (opt_1=='ADCdark00'){
        fonc_seg <- 'ADC'
      }
      else if (opt_1=='brightAnat00'){
        fonc_seg <- 'Anat'
      }
      else if (opt_1=='CBFdark00'){
        fonc_seg <- 'CBF'
      }
      
      tranches <- opt_2
      
      subtitle <- sprintf('Segmentation %s, Tr',fonc_seg)
      for (tr in tranches){
        subtitle <- paste(subtitle,'-',tr)
      }
      
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]# plus tard --------#
      
      # Représentation graphique : fonctionnalité courante
      plot.new()
      par(mfrow=c(2,3),cex.main=1.7, cex.sub=1.5,col.main="black", col.sub="red")
      
      for (jour in jours){
        name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
        name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
        
        cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
        cerveau_seg <- read.table(name_cerveau_seg,header=T)
        
        # on définit hem_sain au jour 00, on ne le modifie plus par la suite
        if (jour=="00"){
          # sélection de l'hémisphère sain au jour 00
          l <- length(cerveau_fonc[,4])
          liste_hem <- rep(FALSE,l)
          liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
          cerveau_hem <- cerveau_fonc[liste_hem,]
          # sélection des tranches pour le suivi temporel
          l <- length(cerveau_hem[,5])
          liste_tr <- rep(FALSE,l)
          for (tr in tranches){
            liste_tr <- ifelse(tr==cerveau_hem$Slice,TRUE,liste_tr)
          }
          tranches_hem <- cerveau_hem[liste_tr,]
          d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],'Hem sain J00'),stringsAsFactors=FALSE)
          colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
          # on crée la liste de niveaux de gris exploitable par density() : hémisphère sain au jour 0
          sain <- cerveau_hem[,4]
          liste.nan <- is.na(sain)
          sain <- sain[!liste.nan] # on retire les valeurs manquantes
        }
        
        cerveau_les <- cerveau_fonc # on va délimiter une enveloppe rectangulaire de la partie lésée sur l'image observée, fonctionnalité et jour courants
        
        if (fonc=='Anat'){# anatomique : trop de valeurs manquantes qui allongent la boucle
          list.nan <- is.na(cerveau_les[,4])
          cerveau_les <- cerveau_les[!list.nan,]
        }
        
        l <- length(cerveau_les[,4])
        liste_tr <- rep(FALSE,l)
        m <- length(cerveau_seg[,4])
        for (i in c(1:m)){#(ix in cerveau_seg$x){
          ix <- cerveau_seg[i,1]
          iy <- cerveau_seg[i,2]
          sli <- cerveau_seg[i,5]
          liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
        }
        cerveau_les <- cerveau_les[liste_tr,] # on garde les abscisses de la partie lésée
        
        l <- length(cerveau_fonc[,4])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(any(tr==cerveau_fonc$Slice),TRUE,liste_tr)
        }
        tranches_fonc <- cerveau_fonc[liste_tr,]
        
        # on crée la liste de niveaux de gris exploitable par density() : cerveau entier
        entieres <- tranches_fonc[,4]
        liste.nan <- is.na(entieres)
        entieres <- entieres[!liste.nan] # on retire les valeurs manquantes
        
        tranches_les <- cerveau_les
        # on crée la liste de niveaux de gris exploitable par density() : zone lésée
        les <- tranches_les[,4]
        liste.nan <- is.na(les)
        les <- les[!liste.nan] # on retire les valeurs manquantes
        
        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        n <- 1#length(entieres)
        nl <- 1#length(les)
        ns <- 1#length(sain)
        
        dst <- density(entieres)
        if (length(les)!=0){dstl <- density(les)}
        dsts <- density(sain)
        
        #plot.new()
        #par(lend="butt")
        title <- sprintf("Rat %s jour %s %s",rat,jour,fonc)
        #title <- paste(title,subtitle,segtitle)
        
        if (opt_3=='pdf'){# on imprime dans un fichier
          
          if (length(les)!=0){
            pdf(file = sprintf("%s/%s_suivi_dens_vol%s_%s-%s.pdf","segmentation_manuelle",num_rat,'ADC',fonc,jour))
            plot(dst$x,dst$y,type="n",main=title,sub=paste(subtitle,segtitle),ylim = c(-0.01*max(dsts$y),2*max(dsts$y)))
            lines(dstl$x, nl/n*dstl$y, lwd = 2, col = "darkred")
            lines(dsts$x, ns/n*dsts$y, lwd = 2, lty = 2, col = "darkblue")
            lines(dst$x, dst$y, lwd = 3, col="gray70")
            
            legend("topright", inset = 0.01, 
                   legend = c("Zone segmentée", "Hémisphère sain J00",
                              "Cerveau entier"),
                   col = c("darkred","darkblue",
                           "gray70"),
                   lty = c(1, 2),#, 1),
                   lwd = 2, pt.cex = 2)
            #print(p)
            dev.off()
          }
        }
        else{# on affiche systématiquement
          if (length(les)!=0){
            plot(dst$x,dst$y,type="n",main=title,sub=paste(subtitle,segtitle),ylim = c(-0.01*max(dsts$y),2*max(dsts$y)))
            lines(dstl$x, nl/n*dstl$y, lwd = 2, col = "darkred")
            lines(dsts$x, ns/n*dsts$y, lwd = 2, lty = 2, col = "darkblue")
            lines(dst$x, dst$y, lwd = 3, col="gray70")
            
            legend("topright", inset = 0.01, 
                   legend = c("Zone segmentée", "Hémisphère sain J00",
                              "Cerveau entier"),
                   col = c("darkred","darkblue",
                           "gray70"),
                   lty = c(1, 2),#, 1),
                   lwd = 2, pt.cex = 2)
          }
        }
      }# sub-plot fait
    }# fonctionnalité vue
  }# option 2 codée : suivi d'une tranche individuelle
  else{
    tranches <- opt_2
    
    if (opt_1=='ADCdark00'){
      fonc_seg <- 'ADC'
    }
    else if (opt_1=='brightAnat00'){
      fonc_seg <- 'Anat'
    }
    else if (opt_1=='CBFdark00'){
      fonc_seg <- 'CBF'
    }
    
    for (fonc in liste_fonc){
      
      segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
      #fonc_seg <- fonc # la segmentation avec fonc_seg sera celle utilisée
      
      #tranches <- liste_s_slice[[fonc_seg]]
      
      subtitle <- 'Tranches ' # indique si nécessaire la fonctionnalité utilisée pour la segmentation : refaire
      for (tr in tranches){
        subtitle <- paste(subtitle,'-',tr)
      }
      
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]# plus tard --------#
      
      # Représentation graphique : fonctionnalité courante
      plot.new()
      par(mfrow=c(2,3),cex.main=1.7, cex.sub=1.2,col.main="black", col.sub="red")
      
      for (jour in jours){
        cerveau_seg <- read.table(sprintf('%s/isch3d-%s-%s-J%s.dat',fonc_seg,fonc_seg,rat,jour),header=T)#,checknames=F)
        cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
        cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]
        
        if (jour=="00"){
          # sélection de l'hémisphère sain au jour 00
          cerveau_hem <- cerveau_fonc[cerveau_seg$Label==2,]
          # sélection des tranches pour le suivi temporel
          l <- length(cerveau_hem[,4])
          liste_tr <- rep(FALSE,l)
          for (tr in tranches){
            liste_tr <- ifelse(any(tr==cerveau_hem$Slice),TRUE,liste_tr)
          }
          tranches_hem <- cerveau_hem[liste_tr,]
          # on crée la liste de niveaux de gris exploitable par density()
          hem_sain <- tranches_hem[,4]
          liste.nan <- is.na(hem_sain)
          hem_sain <- hem_sain[!liste.nan] # on retire les valeurs manquantes de la fonctionnalité pour évaluer sa densité
        }
        # on définit hem_sain au jour 00, on ne le modifie plus par la suite
        
        l <- length(cerveau_fonc[,4])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(any(tr==cerveau_fonc$Slice),TRUE,liste_tr)
        }
        tranches_fonc <- cerveau_fonc[liste_tr,]
        # on crée la liste de niveaux de gris exploitable par density()
        entieres <- tranches_fonc[,4]
        liste.nan <- is.na(entieres)
        entieres <- entieres[!liste.nan] # on retire les valeurs manquantes
        
        l <- length(cerveau_isch[,4])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(any(tr==cerveau_isch$Slice),TRUE,liste_tr)
        }
        tranches_isch <- cerveau_isch[liste_tr,]
        # on crée la liste de niveaux de gris exploitable par density()
        isch <- tranches_isch[,4]
        liste.nan <- is.na(isch)
        isch <- isch[!liste.nan] # on retire les valeurs manquantes
        
        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        n <- 1#length(entieres)
        ni <- 1#length(isch)
        ns <- 1#length(hem_sain)
        
        dst <- density(entieres)
        if (length(isch)!=0){dsti <- density(isch)}
        dsts <- density(hem_sain)
        
        #plot.new()
        #par(lend="butt")
        title <- sprintf("Rat %s jour %s %s",rat,jour,fonc)
        #title <- paste(title,subtitle,segtitle)
        
        if (length(isch)!=0){
          plot(dst$x,dst$y,type="n",main=title,sub=paste(subtitle,segtitle))
          lines(dsti$x, ni/n*dsti$y, lwd = 2, col = "darkred")
          lines(dsts$x, ns/n*dsts$y, lwd = 2, lty = 2, col = "darkblue")
          lines(dst$x, dst$y, lwd = 3, col="gray70")
          
          legend("topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain J00","Cerveau entier"),
                 col = c("darkred","darkblue","gray70"),
                 lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
        }
      }# sub-plot fait
    }# fonctionnalité vue
  }# option 2 codée : suivi temporel d'un volume, cette option est utile seulement pour les rats 26 et 30 sinon cf cas précédent.
}

# ------------------- Suivi temporel des répartitions de niveaux de gris sur un cerveau segmenté, par fonctionnalité. Répertoire fonctionnel_gris du rat concerné. ------------------- #

# On utilise seulement les zones segmentées au jour 00.
# Option 1 : chaîne. Modalité utilisée pour effectuer la segmentation
# Option 2 : sortie. pdf si 'pdf', affichage dans RStudio sinon.

ngris_box_fonc <- function(rat, hemi, opt_1, liste_s_slice,opt_2){
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  repertoires <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  
  for (fonc in liste_fonc){
    
    # On crée la dataframe pour le suivi de la fonctionnalioté courante
    d <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
    if (opt_1=='ADCdark00'){
      fonc_seg <- 'ADC'
    }
    else if (opt_1=='brightAnat00'){
      fonc_seg <- 'Anat'
    }
    else if (opt_1=='CBFdark00'){
      fonc_seg <- 'CBF'
    }
    
    #print(fonc)
    #print(fonc_seg)
    
    tranches <- liste_s_slice[[fonc]]
    subtitle <- sprintf("Rat%s, tranche(s) ",rat)
    for (tr in tranches){
      subtitle <- paste(subtitle,'-',tr)
    }
    
    liste_jr <- liste_jfr[[fonc]]
    jours <- liste_jr[[rat]]
    
    for (jour in jours){
      name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
      name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
      
      cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
      cerveau_seg <- read.table(name_cerveau_seg,header=T)
      
      #print(name_cerveau_fonc)
      #print(name_cerveau_seg)
      
      # on définit hem_sain au jour 00, on ne le modifie plus par la suite
      if (jour=="00"){
        # sélection de l'hémisphère sain au jour 00
        l <- length(cerveau_fonc[,4])
        liste_hem <- rep(FALSE,l)
        liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
        cerveau_hem <- cerveau_fonc[liste_hem,]
        # sélection des tranches pour le suivi temporel
        l <- length(cerveau_hem[,5])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(tr==cerveau_hem$Slice,TRUE,liste_tr)
        }
        tranches_hem <- cerveau_hem[liste_tr,]
        d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],'Hem sain J00'),stringsAsFactors=FALSE)
        colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
      }
      
      # --------- On va délimiter la partie lésée sur l'image observée, fonctionnalité et jour courants --------- #
      cerveau_les <- cerveau_fonc
      
      l <- length(cerveau_les[,4])
      liste_tr <- rep(FALSE,l)
      
      m <- length(cerveau_seg[,4])
      for (i in c(1:m)){
        ix <- cerveau_seg[i,1]
        iy <- cerveau_seg[i,2]
        sli <- cerveau_seg[i,5]
        # On garde seulement, dans cerveau_les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt
        liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
        # Seules les tranches de liste_s_slice sont présentes dans cerveau_seg : pas besoin de trier les tranches du cerveau_les obtenu
      }
      cerveau_les <- cerveau_les[liste_tr,]
      
      # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
      d.les <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lesion'))
      colnames(d.les) <- c("x","y",fonc,"Slice","Zone")
      # --------- Dataframe créée --------- #
      
      l <- length(cerveau_fonc[,4])
      liste_tr <- rep(FALSE,l)
      for (tr in tranches){
        liste_tr <- ifelse(tr==cerveau_fonc$Slice,TRUE,liste_tr)
      }
      cerveau_fonc <- cerveau_fonc[liste_tr,]
      # Volume analysé délimité, avec les valeurs manquantes
      d.entier <- as.data.frame(cbind(cerveau_fonc[,1:2],cerveau_fonc[,4:5],'Cerveau entier'))#liste_lbl))
      colnames(d.entier) <- c("x","y",fonc,"Slice","Zone")
      # dataframe créée
      d.increment <- as.data.frame(rbind(d.entier,d.les,d.sain))#,header=T)
      # données pour le jour courant, fonctionnalité fonc
      num_jour <- jour#num_jours[[jour]]
      d.increment <- as.data.frame(cbind(d.increment,num_jour))
      colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
      
      d <- as.data.frame(rbind(d,d.increment))
    }# dataframe d remplie pour la fonctionnalité courante
    
    liste.nan <- is.na(d[,3])
    d <- d[!liste.nan,]
    
    gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
    
    if (opt_2=='pdf'){
      file_name <- sprintf("%s/%s_suivi_box_vol%s_%s.pdf","segmentation_manuelle",num_rat,opt_1,fonc)
      
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA)
      p <- p + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      #print("dd")
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA)
      p <- p + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      print(p)
      #dev.off()
    }
  }# fonctionnalité vue
  #return(d)
}# opt_1 : segmentation utilisée

# ------------------- Compare des tranches individuelles pour les quatre rats, pour une fonctionnalité. Répertoire benjamin_antoine_labo. ------------------- #

# Option 1 : segmentation, non complètement pris en compte pour le moment.
# Option 2 : sortie dans RStudio ou en pdf, dans le répertoire courant.

# Eventuellement : modifier pour avoir une liste de vecteurs à la place de liste_tr, et de nouveau Tranches(s) en sous titre.

comp_rats_fonc <- function(hemi,liste_r_tr,fonc,opt_1,opt_2){
  repertoires_rats <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  
  segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
  if (opt_1=='ADCdark00'){
    fonc_seg <- 'ADC'
  }
  else if (opt_1=='brightAnat00'){
    fonc_seg <- 'Anat'
  }
  else if (opt_1=='CBFdark00'){
    fonc_seg <- 'CBF'
  }
  
  for (rat in noms_rats){
    
    tranches <- liste_r_tr[[rat]]
    subtitle <- sprintf("Rat%s, tranche(s) ",rat)
    for (tr in tranches){
      subtitle <- paste(subtitle,'-',tr)
    }
    
    liste_jr <- liste_jfr[[fonc]]
    jours <- liste_jr[[rat]]
    
    # On crée la dataframe pour le suivi de la fonctionnalioté courante
    d <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    for (jour in jours){
      name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc]],fonc,rat,jour,fonc,'bg')
      name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
      
      cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
      cerveau_seg <- read.table(name_cerveau_seg,header=T) # on segmentera par rapport à ADC au jour 00.
      
      #hemi <- l_hemi[[rat]]
      # On définit hem_sain au jour 00, on ne le modifie plus par la suite
      if (jour=="00"){
        # sélection de l'hémisphère sain au jour 00
        l <- length(cerveau_fonc[,4])
        liste_hem <- rep(FALSE,l)
        liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
        cerveau_hem <- cerveau_fonc[liste_hem,]
        # sélection des tranches pour le suivi temporel
        l <- length(cerveau_hem[,5])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(tr==cerveau_hem$Slice,TRUE,liste_tr)
        }
        tranches_hem <- cerveau_hem[liste_tr,]
        d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],'Hem sain J00'),stringsAsFactors=FALSE)
        colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
      }
      
      # --------- On va délimiter la partie lésée sur l'image observée, fonctionnalité et jour courants --------- #
      cerveau_les <- cerveau_fonc
      
      l <- length(cerveau_les[,4])
      liste_tr <- rep(FALSE,l)
      
      m <- length(cerveau_seg[,4])
      for (i in c(1:m)){
        ix <- cerveau_seg[i,1]
        iy <- cerveau_seg[i,2]
        sli <- cerveau_seg[i,5]
        # On garde seulement, dans cerveau_les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt
        liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
        # Seules les tranches de liste_s_slice sont présentes dans cerveau_seg : pas besoin de trier les tranches du cerveau_les obtenu
      }
      cerveau_les <- cerveau_les[liste_tr,]
      
      # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
      d.les <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lesion'))
      colnames(d.les) <- c("x","y",fonc,"Slice","Zone")
      # --------- Dataframe créée --------- #
      
      l <- length(cerveau_fonc[,4])
      liste_tr <- rep(FALSE,l)
      for (tr in tranches){
        liste_tr <- ifelse(tr==cerveau_fonc$Slice,TRUE,liste_tr)
      }
      cerveau_fonc <- cerveau_fonc[liste_tr,]
      # Volume analysé délimité, avec les valeurs manquantes
      d.entier <- as.data.frame(cbind(cerveau_fonc[,1:2],cerveau_fonc[,4:5],'Cerveau entier'))#liste_lbl))
      colnames(d.entier) <- c("x","y",fonc,"Slice","Zone")
      # dataframe créée
      d.increment <- as.data.frame(rbind(d.entier,d.les,d.sain))#,header=T)
      # données pour le jour courant, fonctionnalité fonc
      num_jour <- jour#num_jours[[jour]]
      d.increment <- as.data.frame(cbind(d.increment,num_jour))
      colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
      
      d <- as.data.frame(rbind(d,d.increment))
    }# dataframe d remplie pour le rat courant
    
    liste.nan <- is.na(d[,3])
    d <- d[!liste.nan,]
    
    gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
    
    if (opt_2=='pdf'){
      file_name <- sprintf("R%s/%s/%s_suivi_box_vol%s_%s.pdf",num_rat,"segmentation_manuelle",num_rat,opt_1,fonc)
      
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA)
      p <- p + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      #print("dd")
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA)
      p <- p + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      print(p)
      #dev.off()
    }
    
  }# Un jeu de graphiques réalisé : les quatre rats pour la fonctionnalité passée en argument.
}

# ------------------- Même principe, mais on clusterise. ------------------- #

# Option 1 : segmentation, non complètement pris en compte pour le moment.
# Option 2 : sortie dans RStudio ou en pdf, dans le répertoire courant.

# Eventuellement : modifier pour avoir une liste de vecteurs à la place de liste_tr, et de nouveau Tranches(s) en sous titre.

comp_rats_clust_fonc <- function(hemi,cl,lm_fonc,liste_r_tr,fonc,opt_1,opt_2){
  repertoires_rats <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                           'BVf'="fonctionnel_gris",
                           'CBF'="fonctionnel_gris",
                           'CMRO2'="fonctionnel_gris",
                           'SO2map'="fonctionnel_gris",
                           'T1map'="fonctionnel_gris",
                           'VSI'="fonctionnel_gris",
                           'Anat'="anatomique_gris")
  cl_min <- cl[1]
  cl_max <- cl[length(cl)]
  
  segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
  if (opt_1=='ADCdark00'){
    fonc_seg <- 'ADC'
  }
  else if (opt_1=='brightAnat00'){
    fonc_seg <- 'Anat'
  }
  else if (opt_1=='CBFdark00'){
    fonc_seg <- 'CBF'
  }
  
  for (rat in noms_rats){
    
    # --------- Initialisation pour le graphique du rat courant --------- #
    block <- {
      tranches <- liste_r_tr[[rat]]
      subtitle <- sprintf("Rat%s, tranche(s) ",rat)
      for (tr in tranches){
        subtitle <- paste(subtitle,'-',tr)
      }
      
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]
      
      # On crée la dataframe pour le suivi de la fonctionnalioté courante
      d <- data.frame(matrix(ncol = 6, nrow = 0))
      colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
    }
    
    for (jour in jours){
      # --------- On charge les résultats de l'examen du jour pour fonc, et la zone segmentée en fonc_seg --------- #
      block <- {
        name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc]],fonc,rat,jour,fonc,'bg')
        name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
        
        cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
        cerveau_seg <- read.table(name_cerveau_seg,header=T) # on segmentera par rapport à ADC au jour 00. 
      }
      
      # --------- On définit hem_sain au jour 00, on ne le modifie plus par la suite --------- #
      if (jour=="00"){
        # sélection de l'hémisphère sain au jour 00
        l <- length(cerveau_fonc[,4])
        liste_hem <- rep(FALSE,l)
        liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
        cerveau_hem <- cerveau_fonc[liste_hem,]
        # sélection des tranches pour le suivi temporel
        l <- length(cerveau_hem[,5])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(tr==cerveau_hem$Slice,TRUE,liste_tr)
        }
        tranches_hem <- cerveau_hem[liste_tr,]
        d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],'Hem sain J00'),stringsAsFactors=FALSE)
        colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
      }
      
      # --------- On va délimiter la partie lésée sur l'image observée, fonctionnalité et jour courants --------- #
      block <- {
        cerveau_les <- cerveau_fonc
        
        l <- length(cerveau_les[,4])
        liste_tr <- rep(FALSE,l)
        
        m <- length(cerveau_seg[,4])
        for (i in c(1:m)){
          ix <- cerveau_seg[i,1]
          iy <- cerveau_seg[i,2]
          sli <- cerveau_seg[i,5]
          # On garde seulement, dans cerveau_les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt
          liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
          # Seules les tranches de liste_s_slice sont présentes dans cerveau_seg : pas besoin de trier les tranches du cerveau_les obtenu
        }
        cerveau_les <- cerveau_les[liste_tr,]
        
        # On retire les valeurs manquantes pour pouvoir ensuite clusteriser le volume lésé
        list.nan <- is.na(cerveau_les[,4])
        cerveau_les <- cerveau_les[!list.nan,]
        
        # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
        d.les <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lesion'))
        colnames(d.les) <- c("x","y",fonc,"Slice","Zone")
      }
      # --------- On clusterise puis on étiquette les pixels correspondant au cluster 1, et aux suivants respectivement. --------- #
      l.clust <- cluster_jfr_fmin(cerveau_les,cl_min,cl_max,lm_fonc[[fonc]])
      d.les$Zone <- ifelse(l.clust$classification >1,'Lesion 2','Lesion 1')
      
      # --------- On remplit la dataframe à représenter graphiquement, au jour courant. --------- #
      block <- {
        d.increment <- as.data.frame(rbind(d.les,d.sain))#,header=T)
        # données pour le jour courant, fonctionnalité fonc
        num_jour <- jour#num_jours[[jour]]
        d.increment <- as.data.frame(cbind(d.increment,num_jour))
        colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
        
        d <- as.data.frame(rbind(d,d.increment))
      }
    }# dataframe d remplie pour le rat courant
    
    liste.nan <- is.na(d[,3])
    d <- d[!liste.nan,]
    
    gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
    
    if (opt_2=='pdf'){
      file_name <- sprintf("%s_suivi_box_vol%s_%s.pdf",
                           #rat,"segmentation_manuelle",
                           rat,opt_1,fonc)
      
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA,varwidth = TRUE)
      p <- p + scale_fill_manual(values = alpha(c("cyan","red","orange"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA,varwidth = TRUE)
      p <- p + scale_fill_manual(values = alpha(c("cyan","red","orange"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      print(p)
      #dev.off()
    }
    
  }# Un jeu de graphiques réalisé : les quatre rats pour la fonctionnalité passée en argument.
}

# ------------------- Compare pour un rat, une tranche et une fonctionnalité, les statistiques sur les segmentations ADC et CBF. ------------------- #

# Option : sortie dans RStudio ou en pdf, dans le répertoire courant.

# Eventuellement : modifier pour avoir une liste de vecteurs à la place de liste_tr, et de nouveau Tranches(s) en sous titre.

comp_seg_ADCvsCBF <- function(hemi,liste_r_tr,fonc,opt){
  repertoires_rats <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                           'BVf'="fonctionnel_gris",
                           'CBF'="fonctionnel_gris",
                           'CMRO2'="fonctionnel_gris",
                           'SO2map'="fonctionnel_gris",
                           'T1map'="fonctionnel_gris",
                           'VSI'="fonctionnel_gris",
                           'Anat'="anatomique_gris")
  
  segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
  
  for (rat in noms_rats){
    
    tranches <- liste_r_tr[[rat]]
    subtitle <- sprintf("Rat%s, tranche(s) ",rat)
    for (tr in tranches){
      subtitle <- paste(subtitle,'-',tr)
    }
    
    liste_jr <- liste_jfr[[fonc]]
    jours <- liste_jr[[rat]]
    
    # On crée la dataframe pour le suivi de la fonctionnalité courante, segmentation ADC
    d.adc <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d.adc) <- c("x","y",fonc,"Slice","Zone","Jour")
    # On crée la dataframe pour le suivi de la fonctionnalité courante, segmentation CBF
    d.cbf <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d.cbf) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    for (jour in jours){
      name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc]],fonc,rat,jour,fonc,'bg')
      fonc_seg <- 'ADC'
      name_cerveau_segADC <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
      fonc_seg <- 'CBF'
      name_cerveau_segCBF <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
      
      cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
      cerveau_segADC <- read.table(name_cerveau_segADC,header=T) # on segmentera par rapport à ADC au jour 00.
      cerveau_segCBF <- read.table(name_cerveau_segCBF,header=T) # on segmentera par rapport à CBF au jour 00.
      
      # --------- On définit hem_sain au jour 00 ---------#
      if (jour=="00"){
        # sélection de l'hémisphère sain au jour 00
        l <- length(cerveau_fonc[,4])
        liste_hem <- rep(FALSE,l)
        liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
        cerveau_hem <- cerveau_fonc[liste_hem,]
        # sélection des tranches pour le suivi temporel
        l <- length(cerveau_hem[,5])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(tr==cerveau_hem$Slice,TRUE,liste_tr)
        }
        tranches_hem <- cerveau_hem[liste_tr,]
        d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],'Hem sain J00'),stringsAsFactors=FALSE)
        colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
      }
      # --------- On ne le modifie plus par la suite --------- #
      
      # --------- On va délimiter la partie lésée en ADC sur l'image observée, fonctionnalité et jour courants --------- #
      block <- {
        cerveau_les <- cerveau_fonc
        
        l <- length(cerveau_les[,4])
        liste_tr <- rep(FALSE,l)
        
        m <- length(cerveau_segADC[,4])
        for (i in c(1:m)){
          ix <- cerveau_segADC[i,1]
          iy <- cerveau_segADC[i,2]
          sli <- cerveau_segADC[i,5]
          # On garde seulement, dans cerveau_les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt
          liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
          # Seules les tranches de liste_s_slice sont présentes dans cerveau_seg : pas besoin de trier les tranches du cerveau_les obtenu
        }
        cerveau_les <- cerveau_les[liste_tr,]
        
        # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
        d.les.adc <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lesion'))
        colnames(d.les.adc) <- c("x","y",fonc,"Slice","Zone")
      }
      # --------- Dataframe créée --------- #
      
      # --------- On va délimiter la partie lésée en CBF sur l'image observée, fonctionnalité et jour courants --------- #
      block <- {
        cerveau_les <- cerveau_fonc
        
        l <- length(cerveau_les[,4])
        liste_tr <- rep(FALSE,l)
        
        m <- length(cerveau_segCBF[,4])
        for (i in c(1:m)){
          ix <- cerveau_segCBF[i,1]
          iy <- cerveau_segCBF[i,2]
          sli <- cerveau_segCBF[i,5]
          # On garde seulement, dans cerveau_les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt
          liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
          # Seules les tranches de liste_s_slice sont présentes dans cerveau_seg : pas besoin de trier les tranches du cerveau_les obtenu
        }
        cerveau_les <- cerveau_les[liste_tr,]
        
        # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
        d.les.cbf <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lesion'))
        colnames(d.les.cbf) <- c("x","y",fonc,"Slice","Zone")
      }
      # --------- Dataframe créée --------- #
      
      # --------- On construit une dataframe pour le cerveau entier, jour et fonctionnalité courants --------- #
      block <- {
        l <- length(cerveau_fonc[,4])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(tr==cerveau_fonc$Slice,TRUE,liste_tr)
        }
        cerveau_fonc <- cerveau_fonc[liste_tr,]
        # Volume analysé délimité, avec les valeurs manquantes
        d.entier <- as.data.frame(cbind(cerveau_fonc[,1:2],cerveau_fonc[,4:5],'Cerveau entier'))#liste_lbl))
        colnames(d.entier) <- c("x","y",fonc,"Slice","Zone")
      }
      # --------- Dataframe créée --------- #
      
      # On finit de remplir les dataframes utilisées pour l'affichage des statistiques sur la fonctionnalité courante
      # --- Segmentation ADC --- #
      block <- {
        d.increment <- as.data.frame(rbind(d.entier,d.les.adc,d.sain))#,header=T)
        # Fonnées pour le jour courant, fonctionnalité fonc
        num_jour <- jour
        d.increment <- as.data.frame(cbind(d.increment,num_jour))
        colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
        d.adc <- as.data.frame(rbind(d.adc,d.increment))
      }
      # --- Segmentation CBF --- #
      block <- {
        d.increment <- as.data.frame(rbind(d.entier,d.les.cbf,d.sain))#,header=T)
        # Fonnées pour le jour courant, fonctionnalité fonc
        num_jour <- jour
        d.increment <- as.data.frame(cbind(d.increment,num_jour))
        colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
        d.cbf <- as.data.frame(rbind(d.cbf,d.increment))
      }
    }# dataframes d.adc et d.cbf remplies pour le rat courant
    
    liste.nan <- is.na(d.adc[,3])
    d.adc <- d.adc[!liste.nan,]
    liste.nan <- is.na(d.cbf[,3])
    d.cbf <- d.cbf[!liste.nan,]
    
    if (opt=='pdf'){
      fonc_seg <- 'ADC'
      block <- {
        gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
        file_name <- sprintf("R%s/%s/%s_suivi_box_vol%s_%s.pdf",num_rat,"segmentation_manuelle",num_rat,'ADCvsCBF',fonc)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée ADC, et comparer avec ...
        p.adc <- ggplot(d.adc,
                        aes(x=d.adc$Jour,y=d.adc[,3],fill=d.adc$Zone
                        )
        )
        p.adc <- p.adc + geom_boxplot(outlier.shape = NA)
        p.adc <- p.adc + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
        p.adc <- p.adc + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      }
      
      fonc_seg <- 'CBF'
      block <- {
        gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
        file_name <- sprintf("R%s/%s/%s_suivi_box_vol%s_%s.pdf",num_rat,"segmentation_manuelle",num_rat,opt_1,fonc)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée CBF, et comparer avec ...
        p.cbf <- ggplot(d.cbf,
                        aes(x=d.cbf$Jour,y=d.cbf[,3],fill=d.cbf$Zone
                        )
        )
        p.cbf <- p.cbf + geom_boxplot(outlier.shape = NA)
        p.cbf <- p.cbf + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
        p.cbf <- p.cbf + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      }

      p <- grid.arrange(p.adc,p.cbf,ncol=1,nrow=2)
      
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      #print("dd")
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      fonc_seg <- 'ADC'
      block <- {
        gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée ADC, et comparer avec ...
        p.adc <- ggplot(d.adc,
                        aes(x=d.adc$Jour,y=d.adc[,3],fill=d.adc$Zone
                        )
        )
        p.adc <- p.adc + geom_boxplot(outlier.shape = NA)
        p.adc <- p.adc + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
        p.adc <- p.adc + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      }
      
      fonc_seg <- 'CBF'
      block <- {
        gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée CBF, et comparer avec ...
        p.cbf <- ggplot(d.cbf,
                        aes(x=d.cbf$Jour,y=d.cbf[,3],fill=d.cbf$Zone
                        )
        )
        p.cbf <- p.cbf + geom_boxplot(outlier.shape = NA)
        p.cbf <- p.cbf + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
        p.cbf <- p.cbf + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      }
      
      p <- grid.arrange(p.adc,p.cbf,ncol=1,nrow=2)
      print(p)
      #dev.off()
    }
    
  }# Un jeu de graphiques réalisé : les quatre rats pour la fonctionnalité passée en argument.
}

# ------------------- Suivi temporel des densités de niveaux de gris sur un volume lésé clusterisé, par fonctionnalité. Répertoire du rat concerné. ------------------- #

# Mêmes options que ngris_box_clust.

ngris_box_clust <- function(rat, hemi, cl, lm_fonc, opt_1, liste_s_slice, opt_2){
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  repertoires <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  cl_min <- cl[1]
  cl_max <- cl[length(cl)]
  
  for (fonc in liste_fonc){
    
    # On crée la dataframe pour le suivi de la fonctionnalioté courante
    d <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
    fonc_seg <- 'ADC' # la segmentation avec fonc_seg sera celle utilisée
    if (opt_1=='ADCdark00'){
      fonc_seg <- 'ADC'
    }
    else if (opt_1=='brightAnat00'){
      fonc_seg <- 'Anat'
    }
    else if (opt_1=='CBFdark00'){
      fonc_seg <- 'CBF'
    }
    
    tranches <- liste_s_slice[[fonc]]
    subtitle <- sprintf("Rat%s, tranche(s) ",rat)
    for (tr in tranches){
      subtitle <- paste(subtitle,'-',tr)
    }
    
    liste_jr <- liste_jfr[[fonc]]
    jours <- liste_jr[[rat]]
    
    for (jour in jours){
      name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
      name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
      
      cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
      cerveau_seg <- read.table(name_cerveau_seg,header=T)
      
      # on définit hem_sain au jour 00, on ne le modifie plus par la suite
      if (jour=="00"){
        # sélection de l'hémisphère sain au jour 00
        l <- length(cerveau_fonc[,4])
        liste_hem <- rep(FALSE,l)
        liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
        cerveau_hem <- cerveau_fonc[liste_hem,]
        # sélection des tranches pour le suivi temporel
        l <- length(cerveau_hem[,5])
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(tr==cerveau_hem$Slice,TRUE,liste_tr)
        }
        tranches_hem <- cerveau_hem[liste_tr,]
        d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],'Hem sain J00'),stringsAsFactors=FALSE)
        colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
      }
      # --------- Dataframe créée : hémisphère sain --------- #
      
      cerveau_les <- cerveau_fonc
      
      l <- length(cerveau_les[,4])
      liste_tr <- rep(FALSE,l)
      
      # Segmentation du volume lésé, à l'aide d'ADCdark00
      m <- length(cerveau_seg[,4])
      for (i in c(1:m)){
        ix <- cerveau_seg[i,1]
        iy <- cerveau_seg[i,2]
        sli <- cerveau_seg[i,5]
        # On garde seulement, dans cerveau_les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt
        liste_tr <- ifelse(ix==cerveau_les$x&iy==cerveau_les$y&sli==cerveau_les$Slice,TRUE,liste_tr)
        # Seules les tranches de liste_s_slice sont présentes dans cerveau_seg : pas besoin de trier les tranches du cerveau_les obtenu
      }
      cerveau_les <- cerveau_les[liste_tr,]
      
      # On retire les valeurs manquantes pour pouvoir ensuite clusteriser le volume lésé
      list.nan <- is.na(cerveau_les[,4])
      cerveau_les <- cerveau_les[!list.nan,]
      
      # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
      d.les <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lesion 1'))
      colnames(d.les) <- c("x","y",fonc,"Slice","Zone")
      # --------- Dataframe créée : volume lésé --------- #
      
      l.clust <- cluster_jfr_fmin(cerveau_les,cl_min,cl_max,lm_fonc[[fonc]])
      #les.clust <- cbind(cerveau_les,l.clust$classification)
      
      d.les$Zone <- ifelse(l.clust$classification >1,'Lesion 2','Lesion 1')
      # --------- Caractérisation sur la dataframe : volume lésé clusterrisé --------- #
      
      d.increment <- as.data.frame(rbind(d.les,d.sain))#,header=T)
      # données pour le jour courant, fonctionnalité fonc
      num_jour <- jour#num_jours[[jour]]
      d.increment <- as.data.frame(cbind(d.increment,num_jour))
      colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
      
      d <- as.data.frame(rbind(d,d.increment))
    }# dataframe d remplie pour la fonctionnalité courante
    
    liste.nan <- is.na(d[,3])
    d <- d[!liste.nan,]
    
    gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
    
    if (opt_2=='pdf'){
      file_name = sprintf("%s/%s_suivi_clust_box_vol%s_%s.pdf","segmentation_manuelle",num_rat,opt_1,fonc)
      #pdf(file_name)
      
      # On réalise le ggplot
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA,varwidth = TRUE)
      p <- p + scale_fill_manual(values = alpha(c("cyan","red","orange"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      #print("dd")
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA,varwidth = TRUE)
      p <- p + scale_fill_manual(values = alpha(c("cyan","red","orange"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      print(p)
      #dev.off()
    }
  }# fonctionnalité vue
}

# ------------------- Cartographie la distribution de valeurs d'une fonctionnalité, dans l'hémisphère sain, piour un rat. Option : pdf ou sortie RStudio . ------------------- #

# Une fonction auxiliaire : merci à ...

ggplotGrid <- function (l, path, ncol = 1, nrow = 1,
                        width = 8, height = 11, res = 300,
                        pdf.cairo = TRUE, onefile = TRUE, ...) {
  
  # test the classes
  lclass <- sapply(l, function(i) class(i)[1] == "gg")
  if (!all(lclass))
    stop("Provide list with only ggplots!")
  
  # presets
  type <- 'int'
  n <- nrow * ncol # per pdf page
  ggempty <- list(
    ggplot(data.frame()) +
      geom_point() + theme_bw() +
      theme(panel.border = element_rect(color = "white")) +
      scale_x_discrete(breaks=NULL) +
      scale_y_discrete(breaks=NULL))
  
  # get device
  if (!missing(path)) {
    type <- strpart(path, "\\.", 123, roll=T)
    #type <- tolower(substr(path, nchar(path)-2, nchar(path)))
    if (type == 'pdf' && pdf.cairo)
      cairo_pdf(path, width = width, height = height, onefile = onefile, ...)
    else if (type == "ps")
      cairo_ps(path, width = width, height = height, onefile = onefile, ...)
    else if (type == "svg")
      svg(path, width = width, height = height, onefile = onefile, ...)
    else if (type == 'pdf')
      pdf(path, width = width, height = height)
    else if (type == 'png')
      png(path, width = width, height = height, units = 'in', res = res)
    else if (type == "eps") {
      setEPS()
      postscript(file = path, width = width, height = height, ...)
    } else
      stop(paste0('.', type, ' is an invalid ending, use: .pdf .png .eps .svg or .ps.'))
  }
  
  if (type == 'pdf') {
    avail <- TRUE
    while (avail) {
      ln <- length(l)
      # fill with empty to match grid
      if (ln < n) {
        empty <- (ln+1):n
        l[empty] <- ggempty
      }
      do.call(grid.arrange, c(l[1:n], ncol=ncol))
      l <- l[-c(1:n)]
      if (length(l) < 1) avail <- FALSE
    }
  } else {
    do.call(grid.arrange, c(l, ncol=ncol))
  }
  
  if (!missing(path)) dev.off()
  
  invisible(NULL)
  
}

# Une fonction auxiliaire : merci à ...

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Option 1 : jour seul ou suivi
# Option 2 : sortie dans RStudio ou en pdf, dans le répertoire courant.

carte_fonc_rat <- function(rat,fonc,hemi,opt_1,opt_2){
  repertoires_rats <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                           'BVf'="fonctionnel_gris",
                           'CBF'="fonctionnel_gris",
                           'CMRO2'="fonctionnel_gris",
                           'SO2map'="fonctionnel_gris",
                           'T1map'="fonctionnel_gris",
                           'VSI'="fonctionnel_gris",
                           'Anat'="anatomique_gris")
  if (opt_1=="Suivi"){
    # --- Mise en page du ggplot : hémisphère contralatéral examiné sur tous les jours. Pas encore au point : grid arrange et liste de ggplots.
    
    liste_data_jours <- list("00"='',"03"='',"08"='',"15"='',"22"='')# certains jours ne seront pas toujours utilisés
    liste_sub <- list("00"='',"03"='',"08"='',"15"='',"22"='')
    
    liste_jr <- liste_jfr[[fonc]]
    jours <- liste_jr[[rat]]
    
    for (jour in jours){# --- On génère une dataframe par jour
      tr.filename <- sprintf("R%s/%s/%s/liste_R%s_%s_J%s.csv",rat,repertoires_rats[[fonc]],fonc,rat,fonc,jour)
      j.tranches <- read.csv(tr.filename,check.names=F,header=T)
      tranches <- j.tranches[[jour]]
      
      # On crée la dataframe pour le suivi de la fonctionnalité courante
      d <- data.frame(matrix(ncol = 4, nrow = 0))
      colnames(d) <- c("x","y",fonc,"Slice")
      
      subtitle <- sprintf("Rat%s, J%s, tranche(s) ",rat,jour)
      for (tr in tranches){
        subtitle <- paste(subtitle,'-',tr)
      }
      liste_sub[[jour]] <- subtitle
      
      name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc]],fonc,rat,jour,fonc,'bg')
      cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
      
      l <- length(cerveau_fonc[,4])
      liste_hem <- rep(FALSE,l)
      liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
      cerveau_hem <- cerveau_fonc[liste_hem,]
      
      d <- as.data.frame(cbind(cerveau_hem[,1:2],cerveau_hem[,4:5]),stringsAsFactors=FALSE)#,'Hem sain J00'
      colnames(d) <- c("x","y",fonc,"Slice")#,"Zone")
      # Jour courant
      # Seulement la région saine
      liste_data_jours[[jour]] <- d
    }
    # --- if else pour afficher ou enregistrer : comme dans comp_seg_ADCvsCBF
    if (opt_2=='pdf'){
      # nommage du fichier de sortie
      file_name <- sprintf("R%s/segmentation_manuelle/suivi_contralateral_%s.pdf",rat,fonc)
      
      # initialisation de la liste des subggplots
      d0 <- data.frame()
      plot0 <- ggplot(d0) + geom_point()
      plot0 <- plot0 + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL)# + xlim(0, 10) + ylim(0, 100)
      plot0 <- plot0 + ggtitle(bquote(atop(.(sprintf("%s indisponible",fonc)), atop(italic(.(liste_sub[[jour]])), ""))))
      
      plots <- list("00"=plot0,"03"=plot0,"08"=plot0,"15"=plot0,"22"=plot0)
      
      for (jour in jours){
        
        d <- liste_data_jours[[jour]]
        d$Slice <- as.character(d$Slice)
        
        gg_title <- sprintf("Distribution spatiale de %s",fonc)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée ADC, et comparer avec ...
        plots[[jour]] <- ggplot(d,
                                na.rm=TRUE,
                                aes(x=d$Slice,y=d[,3]#,fill=d$Zone
                                )
        )
        plots[[jour]] <- plots[[jour]] + scale_x_discrete(limits=as.character(tranches))
        plots[[jour]] <- plots[[jour]] + geom_boxplot(outlier.shape = NA, fill="lightblue")
        plots[[jour]] <- plots[[jour]] + ggtitle(bquote(atop(.(gg_title), atop(italic(.(liste_sub[[jour]])), "")))) + xlab("Jours") + ylab(fonc)
      }
      
      p <- grid.arrange(plots[["00"]],
                        plots[["03"]],plots[["08"]],plots[["15"]],plots[["22"]],
                        ncol=3,nrow=2)
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      # initialisation de la liste des subggplots
      d0 <- data.frame()
      plot0 <- ggplot(d0) + geom_point()
      plot0 <- plot0 + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL)# + xlim(0, 10) + ylim(0, 100)
      plot0 <- plot0 + ggtitle(bquote(atop(.(sprintf("%s indisponible",fonc)), atop(italic(.(liste_sub[[jour]])), ""))))
      
      plots <- list("00"=plot0,"03"=plot0,"08"=plot0,"15"=plot0,"22"=plot0)
      plots1 <- plots
      
      for (jour in jours){
        
        d <- liste_data_jours[[jour]]
        d$Slice <- as.character(d$Slice)
        
        gg_title <- sprintf("Distribution spatiale de %s",fonc)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée ADC, et comparer avec ...
        plots[[jour]] <- ggplot(d,
                                na.rm=TRUE,
                                aes(x=d$Slice,y=d[,3]#,fill=d$Zone
                                )
        )
        plots[[jour]] <- plots[[jour]] + scale_x_discrete(limits=as.character(tranches))
        plots[[jour]] <- plots[[jour]] + geom_boxplot(outlier.shape = NA, fill="lightblue")
        plots[[jour]] <- plots[[jour]] + ggtitle(bquote(atop(.(gg_title), atop(italic(.(liste_sub[[jour]])), "")))) + xlab("Jours") + ylab(fonc)
      }
      
      #p <- ggplotGrid(plots,#[["00"]],
      #plots[["03"]],plots[["08"]],plots[["15"]],plots[["22"]],
      #                  ncol=3,nrow=2)
      
      #p <- multiplot(plots1,rows=2)#plots[["00"]],plots[["03"]],rows=2)
      
      p <- plot_grid(plots[["00"]],plots[["03"]],plots[["08"]],plots[["15"]],plots[["22"]],nrow=2)
      print(p)
      #dev.off()
    }
  }
  else{# --- Option : jour.
    jour <- opt_1
    tr.filename <- sprintf("R%s/%s/%s/liste_R%s_%s_J%s.csv",rat,repertoires_rats[[fonc]],fonc,rat,fonc,jour)
    j.tranches <- read.csv(tr.filename,check.names=F,header=T)
    tranches <- j.tranches[[jour]]
    
    # On crée la dataframe pour le suivi de la fonctionnalité courante
    d <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(d) <- c("x","y",fonc,"Slice")
    
    subtitle <- sprintf("Rat%s, J %s, tranche(s) ",rat,jour)
    for (tr in tranches){
      subtitle <- paste(subtitle,'-',tr)
    }
    
    name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires_rats[[fonc]],fonc,rat,jour,fonc,'bg')
    cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
    
    l <- length(cerveau_fonc[,4])
    liste_hem <- rep(FALSE,l)
    liste_hem <- ifelse(hemi[1]*cerveau_fonc$x+hemi[2]>cerveau_fonc$y,TRUE,liste_hem)
    cerveau_hem <- cerveau_fonc[liste_hem,]
    
    d <- as.data.frame(cbind(cerveau_hem[,1:2],cerveau_hem[,4:5]),stringsAsFactors=FALSE)#,'Hem sain J00'
    colnames(d) <- c("x","y",fonc,"Slice")#,"Zone")
    # Un seul jour : opt_1
    # Seulement la région saine
    
    liste.nan <- is.na(d[,3])
    d <- d[!liste.nan,]
    
    gg_title <- sprintf("Distribution spatiale de %s",fonc)
    
    if (opt_2=='pdf'){
      d$Slice <- as.character(d$Slice)
      file_name <- #sprintf(".pdf","segmentation_manuelle",num_rat,opt,fonc)
        
        # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
        p <- ggplot(d,
                    aes(x=d$Slice,y=d[,3]#,fill=d$Zone
                    )
        )
      p <- p + geom_boxplot(outlier.shape = NA, fill="lightblue")
      #p <- p + scale_fill_manual(values = alpha(c("blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Tranches") + ylab(fonc)
      
      # On imprime pour que le graphique soit exporté
      print(p)
      dev.off()
      
      # On enregistre le dernier ggplot réalisé
      ggsave(filename=file_name,plot = p)
    }
    else{
      d$Slice <- as.character(d$Slice)
      #print(d$Slice[20:50])
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Slice,
                      y=d[,3]#,
                      #fill=d$Zone,
                      #group = cut_width(x, 0.1)
                  )
      )
      p <- p + scale_x_discrete(limits=as.character(tranches))
      p <- p + geom_boxplot(outlier.shape = NA, fill="lightblue")
      #p <- p + scale_fill_manual(values = alpha(c("blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Tranches") + ylab(fonc)
      print(p)
      #dev.off()
    }
  }
}

# ------------------- Suivi temporel de l'hémisphère contralatéral ------------------- #

#suivi_contra_rat

# ------------------- Suivi temporel de l'étendue d'une zone anormale, graphiques des diférentes fonctionnalités superposés ------------------- #

# Option : chaîne ou liste.
# 'cer' : suivi effectué sur le cerveau entier, ou plutôt sur les slices disponibles chaque jour pour le rat courant ;
# 'tranches' : suivi temporel sur les tranches pour lesquelles cela est possible, toutes fonctionnalités confondues ;
# Sinon : numéro de la tranche choisie pour le suivi temporel.

suivi_etendue_fonc <- function(rat,opt){# cerveau entier ou slices suivables temporellement, représentatives.
  
  vect.fonc <- c('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')
  vect.fonc.color <- c('red','green','blue','orange','purple','black','grey70')
  
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  color.fonc.list <- list('ADC'='red','BVf'='green','CBF'='blue','CMRO2'='orange','SO2map'='purple','T1map'='black','VSI'='grey70')
  
  abs_suivi <- c(0,3,8,15,22)
  # liste des temps disponibles par fonctionnalité
  liste_abs_seg <- list('ADC'=c(),'BVf'=c(),'CBF'=c(),'CMRO2'=c(),'SO2map'=c(),'T1map'=c(),'VSI'=c())
  # liste des aires ou volumes par fonctionnalité
  liste_etendues <- list('ADC'='','BVf'='','CBF'='','CMRO2'='','SO2map'='','T1map'='','VSI'='')
  
  if (opt=='cer'){# suivi sur le cerveau entier, boucle sur les jours existant pour chaque fonctionnalité.
    # Plus grand volume, toutes fonctionnalités et tous jours confondus
    vol_max <- 0
    for (fonc in liste_fonc){
      fonc_seg <- fonc

      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]
      
      # Abscisses : jours de segmentation
      abs_seg <- c()
      # Vecteur des volumes de la zone ischémiée pour la fonctionnalité courante, initialisé
      vols_fonc <- c()
 
      for (jour in jours){# une fenêtre pour la fonctionnalité courante
        
        cerveau_seg <- read.table(sprintf('%s/isch3d-%s-%s-J%s.dat',fonc_seg,fonc_seg,rat,jour),header=T)#,checknames=F)
        cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
        cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]
        
        #tranches_isch <- cerveau_isch[liste_tr,]
        ## on crée la liste de niveaux de gris exploitable par density()
        isch <- cerveau_isch[,4]
        liste.nan <- is.na(isch)
        isch <- isch[!liste.nan] # on retire les valeurs manquantes
        
        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        #n <- length(entieres)
        vol <- length(isch)
        njour <- num_jours[[jour]]
        
        if(vol!=0){
          abs_seg <- cbind(abs_seg,c(njour))
          vols_fonc <- cbind(vols_fonc,c(vol))
        }
      }
      vol_max <- max(max(vols_fonc,vol_max))
      # vecteur rempli pour la fonctionnalité
      liste_etendues[[fonc]] <- vols_fonc
      # liste des temps disponibles remplie pour fonctionnalité courante
      liste_abs_seg[[fonc]] <- abs_seg
    }
    # listes toutes remplies pour les nuages de points correspondant aux fonctionnalités
    # on passe à la représentation graphique
    plot.new()
    par(mfrow=c(1,1),new=T)
    sain <- rep(0,length(abs_suivi))
    plot(x=abs_suivi,
         y=sain,
         col='blue',ylim = range(c(-10, vol_max)),
         xlab = "",
         ylab = "",
         main = sprintf("Zones ischémiées, toutes les tranches observables, rat %s",rat)
         )
    for (fonc in liste_fonc){
      #print(liste_suivi[[fonc]])
      lines(x=liste_abs_seg[[fonc]],y=liste_etendues[[fonc]],col=color.fonc.list[[fonc]],cex=1.5,lwd=1.8)
    }
    legend(title="Suivi temporel sur les jours disponibles",
           "topright", inset = 0.01,
           legend = liste_fonc,#vect.fonc,
           col = vect.fonc.color,#color.fonc.list,
           #lty = c(1, 2, 1),
           lwd = 2, pt.cex = 2
    )
  }
  else if (opt=='tranches'){# suivi sur les slices sélectionnées, boucle sur les jours existant pour chaque fonctionnalité.
    liste_s_slice <- liste_suivi_slice
    # Plus grand volume, toutes fonctionnalités et tous jours confondus
    vol_max <- 0
    
    for (fonc in liste_fonc){
      fonc_seg <- fonc
      tranches <- liste_s_slice[[fonc_seg]]
      
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]
      
      # Abscisses : jours de segmentation
      abs_seg <- c()
      # Vecteur des volumes de la zone ischémiée pour la fonctionnalité courante, initialisé
      vols_fonc <- c()
      
      for (jour in jours){# une fenêtre pour la fonctionnalité courante
        
        cerveau_seg <- read.table(sprintf('%s/isch3d-%s-%s-J%s.dat',fonc_seg,fonc_seg,rat,jour),header=T)#,checknames=F)
        cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
        cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]
        
        l <- length(cerveau_isch$Slice)
        liste_tr <- rep(FALSE,l)
        for (tr in tranches){
          liste_tr <- ifelse(cerveau_isch$Slice==tr,TRUE,liste_tr)
        }
        
        tranches_isch <- cerveau_isch[liste_tr,]
        # on crée la liste de niveaux de gris exploitable par density()
        isch <- tranches_isch[,4]
        liste.nan <- is.na(isch)
        isch <- isch[!liste.nan] # on retire les valeurs manquantes
        
        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        #n <- length(entieres)
        vol <- length(isch)
        njour <- num_jours[[jour]]
        
        if(vol!=0){
          abs_seg <- cbind(abs_seg,c(njour))
          vols_fonc <- cbind(vols_fonc,c(vol))
        }
      }
      vol_max <- max(max(vols_fonc,vol_max))
      # vecteur rempli pour la fonctionnalité
      liste_etendues[[fonc]] <- vols_fonc
      # liste des temps disponibles remplie pour fonctionnalité courante
      liste_abs_seg[[fonc]] <- abs_seg
    }
    # listes toutes remplies pour les nuages de points correspondant aux fonctionalités
    # on passe à la représentation graphique
    plot.new()
    par(mfrow=c(1,1),new=T)
    sain <- rep(0,length(abs_suivi))
    plot(x=abs_suivi,y=sain,col='blue',ylim = range(c(-10, vol_max)),
         main = sprintf("Suivi temporel, tranches d'intérêt, rat %s.",rat),
         xlab = "",
         ylab = "")
    for (fonc in liste_fonc){
      lines(x=liste_abs_seg[[fonc]],y=liste_etendues[[fonc]],col=color.fonc.list[[fonc]],cex=1.5,lwd=1.8)
    }
    legend(title="Suivi temporel sur les jours disponibles",
           "topright", inset = 0.01,
           legend = liste_fonc,#vect.fonc,
           col = vect.fonc.color,#color.fonc.list,
           #lty = c(1, 2, 1),
           lwd = 2, pt.cex = 2
    )
  }
  else{
    # suivi sur la tranche sélectionnée, boucle sur les jours existant pour chaque fonctionnalité.
    num_tranche <- opt
    # Plus grand volume, toutes fonctionnalités et tous jours confondus
    vol_max <- 0
    
    for (fonc in liste_fonc){
      fonc_seg <- fonc

      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]
      
      # Abscisses : jours de segmentation
      abs_seg <- c()
      # Vecteur des volumes de la zone ischémiée pour la fonctionnalité courante, initialisé
      vols_fonc <- c()
      
      for (jour in jours){# une fenêtre pour la fonctionnalité courante
        
        cerveau_seg <- read.table(sprintf('%s/isch3d-%s-%s-J%s.dat',fonc_seg,fonc_seg,rat,jour),header=T)#,checknames=F)
        cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
        cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]

        tranches_isch <- cerveau_isch[cerveau_isch$Slice==num_tranche,]
        # on crée la liste de niveaux de gris exploitable par density()
        isch <- tranches_isch[,4]
        liste.nan <- is.na(isch)
        isch <- isch[!liste.nan] # on retire les valeurs manquantes
        
        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        #n <- length(entieres)
        vol <- length(isch)
        njour <- num_jours[[jour]]
        
        if(vol!=0){
          abs_seg <- cbind(abs_seg,c(njour))
          vols_fonc <- cbind(vols_fonc,c(vol))
        }
      }
      vol_max <- max(max(vols_fonc,vol_max))
      # vecteur rempli pour la fonctionnalité
      liste_etendues[[fonc]] <- vols_fonc
      # liste des temps disponibles remplie pour fonctionnalité courante
      liste_abs_seg[[fonc]] <- abs_seg
    }
    # listes toutes remplies pour les nuages de points correspondant aux fonctionalités
    # on passe à la représentation graphique
    plot.new()
    par(mfrow=c(1,1),new=T)
    sain <- rep(0,length(abs_suivi))
    plot(x=abs_suivi,y=sain,col='blue',ylim = range(c(-10, 800)),
         xlab = "",
         ylab = "",
         main=sprintf("Zones ischémiées, tranche %i, rat %s",num_tranche,rat)
    )
    for (fonc in liste_fonc){
      lines(x=liste_abs_seg[[fonc]],y=liste_etendues[[fonc]],col=color.fonc.list[[fonc]],lwd=1.8)
    }
    legend(title="Suivi temporel sur les jours disponibles",
           "topright", inset = 0.01,
           legend = liste_fonc,#vect.fonc,
           col = vect.fonc.color,#color.fonc.list,
           #lty = c(1, 2, 1),
           lwd = 2, pt.cex = 2
    )
  }# fin pour la troisième option : tranche unique, valeur numérique de l'option.
}




##########################################################################################







