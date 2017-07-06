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

# ------------------- Mâme chose, avec les tranches dénoyautées suivables temporellement  ------------------- #

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

cluster_jfr_fmin <- function(data,cl_min,cl_max,min){# on adapte le rabotage à la fonctionnalité
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
    else{
      d <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,day,fonc),header=T)
    }
    d <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,day,fonc),header=T)
    
    #print(d[1:50,4])
    liste.nan <- is.na(d[,4])
    #print(liste.nan[1:50])
    d <- d[!liste.nan,] # On retire les valeurs NaN de la dataframe, pas besoin de réassembler ensuite
    #print(d[1:50,4])
    
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

    tranche <- cerveau[cerveau$Slice==num_tranche,]
    t.nan <- is.na(tranche[,4])
    tranche <- tranche[!t.nan,]

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

# Option 1 : chaîne ou liste.
# 'cer' : suivi effectué sur le cerveau entier, ou plutôt sur les slices disponibles chaque jour pour le rat courant.
# Sinon : liste des numéros des slices utilisées pour le suivi temporel.

# Option 2 : chaîne de caractères, liste pour une version améliorée de la fonction.
# Désigne une éventuelle fonctionnalité, par exemple l'ADC, dont la segmentation choisie à l'étape 2 sera utilisée pour suivre les valeurs de toutes les fonctionnalités.
# Vide '' : le suivi de chaque fonctionnalité est effectué avec la segmentation qu'elle a elle-même permis de réaliser.


dgris_temp_fonc <- function(rat, opt_1, opt_2){

  if (any(opt_1=='cer')){# suivi sur le cerveau entier, boucle sur les jours existant pour chaque fonctionnalité.
    for (fonc in liste_fonc){
      #segtitle <- '' # indique si nécessaire la fonctionnalité utilisée pour la segmentation
      fonc_seg <- fonc # la segmentation avec fonc_seg sera celle utilisée
      #if (opt_2!=''){
      #  fonc_seg <- opt_2 # fonc prend la valeur de l'argument optionnel opt_1 si celui-ci est non vide
      #  segtitle <- paste("Segmentation ",opt_2)
      #}
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
          
          legend("topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain","Cerveau entier"),
                 col = c("darkred","darkblue","gray70"),
                 lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
        }
      }
      # sub-plot fait
    }
    # fonctionnalité vue
  }
  else{# suivi sur les slices sélectionnées, boucle sur les jours existant pour chaque fonctionnalité.
    liste_s_slice <- opt_1

    #subtitle <- paste("Tr ",tr) # indiquie les tranches du suivi
    
    for (fonc in liste_fonc){

      segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
      fonc_seg <- fonc # la segmentation avec fonc_seg sera celle utilisée
      #if (opt_2!=''){
      #  fonc_seg <- opt_2 # fonc prend la valeur de l'argument optionnel opt_1 si celui-ci est non vide
      #  segtitle <- paste("SEG : ",opt_2)
      #}
      tranches <- liste_s_slice[[fonc_seg]]
      
      subtitle <- 'Tranches ' # indique si nécessaire la fonctionnalité utilisée pour la segmentation
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
          liste_tr <- ifelse(any(tranches==cerveau_hem$Slice),TRUE,FALSE)
          tranches_hem <- cerveau_hem[liste_tr,]
          # on crée la liste de niveaux de gris exploitable par density()
          hem_sain <- tranches_hem[,4]
          liste.nan <- is.na(hem_sain)
          hem_sain <- hem_sain[!liste.nan] # on retire les valeurs manquantes de la fonctionnalité pour évaluer sa densité
        }
        # on définit hem_sain au jour 00, on ne le modifie plus par la suite
        
        liste_tr <- ifelse(any(tranches==cerveau_fonc$Slice),TRUE,FALSE)
        tranches_fonc <- cerveau_fonc[liste_tr,]
        # on crée la liste de niveaux de gris exploitable par density()
        entieres <- tranches_fonc[,4]
        liste.nan <- is.na(entieres)
        entieres <- entieres[!liste.nan] # on retire les valeurs manquantes
        
        liste_tr <- ifelse(any(tranches==cerveau_isch$Slice),TRUE,FALSE)
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
          
          legend("topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain","Cerveau entier"),
                 col = c("darkred","darkblue","gray70"),
                 lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
        }
      }# sub-plot fait

      
      }# fonctionnalité vue
    }
}

# ------------------- Suivi temporel de l'étendue d'une zone anormale, graphiques des diférentes fonctionnalités superposés ------------------- #

# Option : chaîne ou liste.
# 'cer' : suivi effectué sur le cerveau entier, ou plutôt sur les slices disponibles chaque jour pour le rat courant.
# Sinon : liste des numéros des slices utilisées pour le suivi temporel.

suivi_etendue_fonc <- function(rat,opt){# cerveau entier ou slices suivables temporellement, représentatives.
  
  vect.fonc <- c('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')
  vect.fonc.color <- c('red','green','blue','orange','purple','black','grey70')
  
  liste_suivi <- list('ADC'='','BVf'='','CBF'='','CMRO2'='','SO2map'='','T1map'='','VSI'='')
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  color.fonc.list <- list('ADC'='red','BVf'='green','CBF'='blue','CMRO2'='orange','SO2map'='purple','T1map'='black','VSI'='grey70')
  
  if (any(opt=='cer')){# suivi sur le cerveau entier, boucle sur les jours existant pour chaque fonctionnalité.
    for (fonc in liste_fonc){
      
      fonc_seg <- fonc # la segmentation avec fonc_seg sera celle utilisée
      #if (opt_2!=''){
      #  fonc_seg <- opt_2 # fonc prend la valeur de l'argument optionnel opt_1 si celui-ci est non vide
      #}
      
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]
      
      
      # Abscisses : jours de segmentation
      abs_seg <- c()
      # Vecteur des aires de la zone ischémiée pour la fonctionnalité courante, initialisé
      aires_fonc <- c()
      
      for (jour in jours){# une fenêtre pour la fonctionnalité courante
        cerveau_seg <- read.table(sprintf('%s/isch3d-%s-%s-J%s.dat',fonc_seg,fonc_seg,rat,jour),header=T)#,checknames=F)
        cerveau_fonc <- read.table(sprintf("%s/%s-J%s-%s-bg-all.dat",fonc,rat,jour,fonc),header=T)
        cerveau_isch <- cerveau_fonc[cerveau_seg$Label==1,]

        isch <- cerveau_isch[,4]
        liste.nan <- is.na(isch)
        isch <- isch[!liste.nan] # on retire les valeurs manquantes
        
        ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
        aire <- length(isch)
        njour <- num_jours[[jour]]
        #aires_fonc[njour] <- aire
        
        if(aire!=0){
          abs_seg <- cbind(abs_seg,c(num_jours[[jour]]))
          aires_fonc <- cbind(aires_fonc,c(aire))
        }

        #legend()#"topright", inset = 0.01, legend = c("Zone ischémiée", "Hémisphère sain","Cerveau entier"),
        #col = c("darkred","darkblue","gray70"),
        #lty = c(1, 2, 1), lwd = 2, pt.cex = 2)
      }
      # vecteur rempli pour la fonctionnalité
      liste_suivi[[fonc]] <- aires_fonc
    }
    # listes toutes remplies pour les nuages de points correspondant aux fonctionnalités
    # on passe à la représentation graphique
    plot.new()
    par(mfrow=c(1,1),new=T)
    abs <- c(1,2,3,4,5)
    sain <- rep(0,length(abs_seg))
    plot(x=abs,y=sain,col='blue',ylim = range(c(-10, 5000)), xlab = "", ylab = "")
    for (fonc in liste_fonc){
      #print(liste_suivi[[fonc]])
      lines(x=abs_seg,y=liste_suivi[[fonc]],col=color.fonc.list[[fonc]])
      #polygon(x=abs,y=liste_suivi[[fonc]],col='grey70',border=NA) # on superpose les nuages de points
    }
    legend(title="Suivi temporel : étendue des zones ischémiées",
           "topright", inset = 0.01,
           legend = liste_fonc,#vect.fonc,
           col = vect.fonc.color,#color.fonc.list,
           #lty = c(1, 2, 1),
           lwd = 2, pt.cex = 2
    )
  }
  else{# suivi sur les slices sélectionnées, boucle sur les jours existant pour chaque fonctionnalité.

    for (fonc in liste_fonc){
      fonc_seg <- fonc
      liste_suivi_slice <- opt
      print(opt)
      tranches <- liste_suivi_slice[[fonc_seg]]
      print(tranches)
      
      liste_jr <- liste_jfr[[fonc]]
      jours <- liste_jr[[rat]]
      
      # Abscisses : jours de segmentation
      abs_seg <- c()
      # Vecteur des aires de la zone ischémiée pour la fonctionnalité courante, initialisé
      aires_fonc <- c()
      
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
        aire <- length(isch)
        njour <- num_jours[[jour]]
        
        if(aire!=0){
          abs_seg <- cbind(abs_seg,c(njour))
          aires_fonc <- cbind(aires_fonc,c(aire))
        }
        
        #aires_fonc[njour] <- aire
        #print(aire)
        #ns <- length(hem_sain)

      }
      # vecteur rempli pour la fonctionnalité
      liste_suivi[[fonc]] <- aires_fonc
    }
    # listes toutes remplies pour les nuages de points correspondant aux fonctionalités
    # on passe à la représentation graphique
    plot.new()
    par(mfrow=c(1,1),new=T)
    sain <- rep(0,length(abs_seg))
    plot(x=abs_seg,y=sain,col='blue',ylim = range(c(-10, 5000)), xlab = "", ylab = "")
    for (fonc in liste_fonc){
      #print(liste_suivi[[fonc]])
      lines(x=abs_seg,y=liste_suivi[[fonc]],col=color.fonc.list[[fonc]])
    }
    legend(title="Suivi temporel : étendue des zones ischémiées",
           "topright", inset = 0.01,
           legend = liste_fonc,#vect.fonc,
           col = vect.fonc.color,#color.fonc.list,
           #lty = c(1, 2, 1),
           lwd = 2, pt.cex = 2
    )
  } # deuxième option prévue : avec les tranches représentatives
}




##########################################################################################
