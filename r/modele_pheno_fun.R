### Première volée de fonctions d'évolution : scalaires ; rat 11, une modalité, un jour. ###


# dep et arr sont les deux vecteurs qui contiennent respectivement les antécédents et images qui définissent la fonction affine par morceaux ...

trans_box <- function(x,dep,arr){
  if (x<=dep[1]){
    res <- arr[1]
  }
  else if (x>=dep[length(dep)]){
    res <- arr[length(arr)]
  }
  else{# on construit une fonction affine par morceaux
    c <- TRUE
    i <- 1
    while (c==TRUE){
      i <- i+1
      c <- (dep[i]<x)
    }
    a <- dep[i-1]
    b <- dep[i]
    theta <- (x-a)/(b-a)
    res <- (1-theta)*arr[i-1]+theta*arr[i]
  }
  return(res)
}

## Paramètres pour les fonctions de transition

# ------- Rat 11 : émanent des statistiques sur la segmentation ADC. ------- #

block <- {
  liste_succ_ADC_11 <- list("00"=list("dep"=c(600,700,725,750,850),"arr"=c(500,750,850,1000,1400)),
                            "03"=list("dep"=c(500,750,850,1000,1400),"arr"=c(450,1000,1200,1450,2100)),
                            "08"=list("dep"=c(450,1000,1200,1450,2100),"arr"=c(800,1200,2100,2300,3300)),
                            "15"=list("dep"=c(300,1200,2100,2300,3300),"arr"=c(800,1800,2200,2550,3300))
  )
  
  liste_succ_BVf_11 <- list("00"=list("dep"=c(0,0.5,1,1.5,2),"arr"=c(0,2.5,4,6.5,11.5)),
                            "03"=list("dep"=c(0,2.5,4,6.5,11.5),"arr"=c(0.5,5,7,12,23)),
                            "08"=list("dep"=c(0.5,5,7,12,23),"arr"=c(0,4,5.5,6.5,12)),
                            "15"=list("dep"=c(0,4,5.5,6.5,12),"arr"=c(0,0.5,2,3,8))
  )
  
  liste_succ_CBF_11 <- list("00"=list("dep"=c(0,15,25,40,80),"arr"=c(0,45,75,150,310)),
                            "03"=list("dep"=c(0,45,75,150,310),"arr"=c(0,10,20,30,60)),
                            "08"=list("dep"=c(0,10,20,30,60),"arr"=c(0,45,70,100,190)),
                            "15"=list("dep"=c(0,45,70,100,190),"arr"=c(0,20,40,80,170))
  )
  
  liste_succ_SO2map_11<- list("00"=list("dep"=c(0,10,25,40,72),"arr"=c(2,47,65,75,100)),
                              "03"=list("dep"=c(5,47,65,75,100),"arr"=c(18,53,73,80,100)),
                              "15"=list("dep"=c(18,53,73,80,100),"arr"=c(0,35,60,74,100))
  )
  
  liste_succ_T1map_11 <- list("00"=list("dep"=c(1150,1450,1550,1700,1950),"arr"=c(1155,1455,1555,1750,2050)),
                              "03"=list("dep"=c(1155,1455,1555,1750,2050),"arr"=c(1200,1600,1800,2100,2700)),
                              "08"=list("dep"=c(1200,1600,1800,2100,2700),"arr"=c(1600,2400,2600,2800,3500)),
                              "15"=list("dep"=c(1600,2400,2600,2800,3500),"arr"=c(850,2100,2700,3100,3500))
  )
  
  liste_succ_VSI_11 <- list("00"=list("dep"=c(1,7,12,20,40),"arr"=c(0,8,13,19,35)),
                            "03"=list("dep"=c(0,8,13,19,35),"arr"=c(2,11,15,21,34)),
                            "08"=list("dep"=c(2,11,15,21,34),"arr"=c(3,10,14,20,34)),
                            "15"=list("dep"=c(3,10,14,20,34),"arr"=c(0,6,10,15,29))
  )
}

liste_box_11 <- list("ADC"=liste_succ_ADC_11 ,
                     "BVf"=liste_succ_BVf_11 ,
                     "CBF"=liste_succ_CBF_11 ,
                     "SO2map"=liste_succ_SO2map_11 ,
                     "T1map"=liste_succ_T1map_11 ,
                     "VSI"=liste_succ_VSI_11
)

# ADC : pour le moment aucune relation avec les autres modalités. #

succ11_ADC_00 <- function(x){
  dep_arr <- liste_succ_ADC_11[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_ADC_03 <- function(x){
  dep_arr <- liste_succ_ADC_11[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_ADC_08 <- function(x){
  dep_arr <- liste_succ_ADC_11[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_ADC_15 <- function(x){
  dep_arr <- liste_succ_ADC_11[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_ADC_11 <- list("00"=succ11_ADC_00,
                           "03"=succ11_ADC_03,
                           "08"=succ11_ADC_08,
                           "15"=succ11_ADC_15
)

# BVf #

succ11_BVf_00 <- function(x){
  dep_arr <- liste_succ_BVf_11[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_BVf_03 <- function(x){
  dep_arr <- liste_succ_BVf_11[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_BVf_08 <- function(x){
  dep_arr <- liste_succ_BVf_11[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_BVf_15 <- function(x){
  dep_arr <- liste_succ_BVf_11[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_BVf_11 <- list("00"=succ11_BVf_00,
                           "03"=succ11_BVf_03,
                           "08"=succ11_BVf_08,
                           "15"=succ11_BVf_15
)

# CBF : même chose que précédemment. #

succ11_CBF_00 <- function(x){
  dep_arr <- liste_succ_CBF_11[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_CBF_03 <- function(x){
  dep_arr <- liste_succ_CBF_11[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_CBF_08 <- function(x){
  dep_arr <- liste_succ_CBF_11[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_CBF_15 <- function(x){
  dep_arr <- liste_succ_CBF_11[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_CBF_11 <- list("00"=succ11_CBF_00,
                           "03"=succ11_CBF_03,
                           "08"=succ11_CBF_08,
                           "15"=succ11_CBF_15
)

# SO2map #

succ11_SO2map_00 <- function(x){
  dep_arr <- liste_succ_SO2map_11[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_SO2map_03 <- function(x){
  dep_arr <- liste_succ_SO2map_11[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_SO2map_15 <- function(x){
  dep_arr <- liste_succ_SO2map_11[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_SO2map_11 <- list("00"=succ11_SO2map_00,
                           "03"=succ11_SO2map_03,
                           "15"=succ11_SO2map_15
)

# T1map #

succ11_T1map_00 <- function(x){
  dep_arr <- liste_succ_T1map_11[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_T1map_03 <- function(x){
  dep_arr <- liste_succ_T1map_11[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_T1map_08 <- function(x){
  dep_arr <- liste_succ_ADC_11[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_T1map_15 <- function(x){
  dep_arr <- liste_succ_T1map_11[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_T1map_11 <- list("00"=succ11_T1map_00,
                           "03"=succ11_T1map_03,
                           "08"=succ11_T1map_08,
                           "15"=succ11_T1map_15
)

# VSI : #

succ11_VSI_00 <- function(x){
  dep_arr <- liste_succ_VSI_11[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_VSI_03 <- function(x){
  dep_arr <- liste_succ_VSI_11[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_VSI_08 <- function(x){
  dep_arr <- liste_succ_VSI_11[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_VSI_15 <- function(x){
  dep_arr <- liste_succ_VSI_11[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_VSI_11 <- list("00"=succ11_VSI_00,
                           "03"=succ11_VSI_03,
                           "08"=succ11_VSI_08,
                           "15"=succ11_VSI_15
                           )
# Toutes fonctionnalités confondues

liste_fsucc_11 <- list("ADC"=liste_fsucc_ADC_11 ,
                       "BVf"=liste_fsucc_BVf_11 ,
                       "CBF"=liste_fsucc_CBF_11 ,
                       "SO2map"=liste_fsucc_SO2map_11 ,
                       "T1map"=liste_fsucc_T1map_11 ,
                       "VSI"=liste_fsucc_VSI_11
                       )


# ------- Rat 19 ------- #
# ------- Rat 26 ------- #
# ------- Rat 30 ------- #


## ------- Fonctions de transitikon : tous rats confondus ------- ##

liste_fsucc <- list("11"=liste_fsucc_11,
                    "19"='',#liste_fsucc_19,
                    "26"='',#liste_fsucc_26,
                    "30"=''#liste_fsucc_30
)

### Deuxième volée de fonctions : dans un répertoire de benjamin_antoine_labo, utilise des versions vectorisées des fonctions précédentes. ###
## Prennent en arguments des coordonnées pour le coin inférieur gauche, et la taille du carré de pixels. ##

# Etape préliminaire : fonction graphique pour tester des fonctions de transition
# options : 1- segmentation choisie ; 2- représentation choisie, densités ou boîtes ; 3- sortie choisie, plot R ou pdf dans un répertoire.
comp_succ_suivi <- function(rat,hemi,fonc,liste_s_slice,opt_1,opt_2,opt_3){
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  repertoires <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  
  segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
  fonc_seg <- 'ADC' # la segmentation avec fonc_seg sera celle utilisée
  if (opt_1=='ADCdark00'){
    fonc_seg <- 'ADC'
  }
  else if (opt_1=='brightAnat00'){
    fonc_seg <- 'Anat'
  }
  else if (opt_2=='CBFdark00'){
    fonc_seg <- 'CBF'
  }
  
  tranches <- liste_s_slice[[fonc]]
  subtitle <- sprintf("Rat %s, tranche(s) ",rat)
  for (tr in tranches){
    subtitle <- paste(subtitle,'-',tr)
  }
  
  liste_jr <- liste_jfr[[fonc]]
  jours <- liste_jr[[rat]]
  
  if (opt_2=='dens'){
    plot.new()
    par(mfrow=c(2,3),cex.main=1.7, cex.sub=1.5,col.main="black", col.sub="red")
    
    for (jour in jours){
      # On crée la dataframe auxiliaire : valeurs calculées au jour d'examen précédent
      if(jour!="00"){
        cerveau_aux <- cerveau_pr
      }
      name_cerveau_fonc <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
      name_cerveau_seg <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc_seg]],fonc_seg,rat,jour,fonc_seg,'dark')
      
      cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
      cerveau_seg <- read.table(name_cerveau_seg,header=T)
      
      # On définit hem_sain au jour 00, on ne le modifie plus par la suite.
      # Structure de données : tranches_hem
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
        
        # on crée la liste de niveauxs de gris exploitable par density() : hémisphère sain
        sain <- tranches_hem[,4]
        liste.nan <- is.na(sain)
        sain <- sain[!liste.nan]
      }
      
      # --------- On va délimiter la partie lésée sur l'image observée. --------- #
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
      # --------- Dataframe créée --------- #
      
      # --------- On initialise : cerveau_pr, auquel on appliquera successivement les fonctions de l_fsucc, et cerveau_aux, dont on représente la densité au juor courant. --------- #
      if (jour=="00"){
        cerveau_pr <- cerveau_les
        cerveau_aux <- cerveau_les
      }
      if (jour != "22"){#"15" pour le rat 30
        #cerveau_pr <- cerveau_les ----> on recommence avec le cerveau_pr calculé à l'étape précédente - boucle sur les jours d'examen.
        # Cerveau lésé pour l'examen suivant : à construire
        liste.nan <- is.na(cerveau_pr[,4])
        cerveau_pr <- cerveau_pr[!liste.nan,]
        
        ll_fsucc <- liste_fsucc[[rat]] 
        l_fsucc <- ll_fsucc[[fonc]]
        fsucc <- l_fsucc[[jour]]
        #print(fsucc(1))
        # Fonction utilisée pour prédire les résultats des prochains examens à l'aide de celui du jours courant
        cerveau_pr[[fonc]] <- sapply(cerveau_pr[[fonc]],fsucc)
        # Cet objet sera utilisé au prochain tour de boucle
      }
      # --------- --------- #
      
      # on crée la liste de niveaux de gris exploitable par density() : zone lésée prévue
      aux <- cerveau_aux[,4]
      liste.nan <- is.na(aux)
      aux <- aux[!liste.nan] # on retire les valeurs manquantes
      
      # on crée la liste de niveaux de gris exploitable par density() : zone lésée
      les <- cerveau_les[,4]
      liste.nan <- is.na(les)
      les <- les[!liste.nan] # on retire les valeurs manquantes
      
      ## Eventuellement, on oublie la normalisation pour représenter les courbes d'effectifs.
      na <- 1#length(aux)
      nl <- 1#length(les)
      ns <- 1#length(sain)
      
      dsta <- density(aux)
      if (length(les)!=0){dstl <- density(les)}
      dsts <- density(sain)
      
      #plot.new()
      #par(lend="butt")
      title <- sprintf("Rat %s jour %s %s",rat,jour,fonc)
      #title <- paste(title,subtitle,segtitle)
      
      if (opt_3=='pdf'){# on imprime dans un fichier
        
        if (length(les)!=0){
          pdf(file = sprintf("%s/%s_suivi_dens_vol%s_%s-%s.pdf","modele_ph",num_rat,'ADC',fonc,jour))
          plot(dsts$x,dsts$y,type="n",main=title,sub=paste(subtitle,segtitle),ylim = c(-0.01*max(dsts$y),2.5*max(dsts$y)))
          lines(dstl$x, nl/ns*dstl$y, lwd = 2, col = "darkred")
          lines(dsts$x, ns/ns*dsts$y, lwd = 2, lty = 2, col = "darkblue")
          lines(dsta$x, na/ns*dsta$y, lwd = 3, col="orange")
          
          legend("topright", inset = 0.01, 
                 legend = c("Zone segmentée : valeurs mesurées", "Hémisphère sain J00",
                            "Zone segmentée, valeurs calculées"),
                 col = c("darkred","darkblue",
                         "orange"),
                 lty = c(1, 2),#, 1),
                 lwd = 2, pt.cex = 2)
          #print(p)
          dev.off()
        }
      }
      else{# on affiche systématiquement
        if (length(les)!=0){
          plot(dsts$x,dsts$y,type="n",main=title,sub=paste(subtitle,segtitle),ylim = c(-0.01*max(dsts$y),2.5*max(dsts$y)))
          lines(dstl$x, nl/ns*dstl$y, lwd = 2, col = "darkred")
          lines(dsts$x, ns/ns*dsts$y, lwd = 2, lty = 2, col = "darkblue")
          lines(dsta$x, na/ns*dsta$y, lwd = 3, col="orange")
          
          #legend("topright", inset = 0.01, 
          #legend = c("Zone segmentée : valeurs mesurées", "Hémisphère sain J00","Zone segmentée, valeurs calculées"),
          #col = c("darkred","darkblue","gray70"),
          #lty = c(1, 2),#, 1),
          #lwd = 2, pt.cex = 2
          #)
        }
        
      }
    }
  }
  else if (opt_2=='box'){
    # On crée la dataframe pour le suivi de la fonctionnalité courante
    d <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    for (jour in jours){
      # On crée la dataframe auxiliaire : valeurs calculées au jour d'examen précédent
      if(jour!="00"){
        d.aux <- as.data.frame(cbind(cerveau_pr[,1:2],cerveau_pr[,4:5],'Lésion prévue'))
        colnames(d.aux) <- c("x","y",fonc,"Slice","Zone")
      }
      name_cerveau_fonc <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
      name_cerveau_seg <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc_seg]],fonc_seg,rat,jour,fonc_seg,'dark')
      
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
      
      # --------- On va délimiter la partie lésée sur l'image observée. --------- #
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
      
      # --------- On construit une data.frame contenant les valeurs prédites pour la fonctionnalité --------- #
      if (jour=="00"){
        cerveau_pr <- cerveau_les
        d.aux <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lésion prévue'))
        colnames(d.aux) <- c("x","y",fonc,"Slice","Zone")
      }
      if (jour != "22"){#"15" pour le rat 30
        #cerveau_pr <- cerveau_les ----> on recommence avec le cerveau_pr calculé à l'étape précédente - boucle sur les jours d'examen.
        # Cerveau lésé pour l'examen suivant : à construire
        liste.nan <- is.na(cerveau_pr[,4])
        cerveau_pr <- cerveau_pr[!liste.nan,]
        
        ll_fsucc <- liste_fsucc[[rat]] 
        l_fsucc <- ll_fsucc[[fonc]]
        fsucc <- l_fsucc[[jour]]
        #print(fsucc(1))
        # Fonction utilisée pour prédire les résultats des prochains examens à l'aide de celui du jours courant
        cerveau_pr[[fonc]] <- sapply(cerveau_pr[[fonc]],fsucc)
        # Cet objet sera utilisé au prochain tour de boucle
      }
      # --------- Dataframe construite : prédistions ou valeurs au jour 00. --------- #
      
      d.increment <- as.data.frame(rbind(d.les,d.aux,d.sain))#,header=T)
      # --------- Assemblage, jour courant --------- #
      num_jour <- jour
      d.increment <- as.data.frame(cbind(d.increment,num_jour))
      colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
      d <- as.data.frame(rbind(d,d.increment))
    }# dataframe d remplie
    
    liste.nan <- is.na(d[,3])
    d <- d[!liste.nan,]
    
    gg_title <- sprintf("Evolution de %s, segmentation %s",fonc,fonc_seg)
    
    if (opt_3=='pdf'){
      file_name <- sprintf("%s/%s_suivi_box_vol%s_%s.pdf","segmentation_manuelle",num_rat,opt_1,fonc)
      
      # On va représenter l'évolution des valeurs de fonc sur la zone lésée, et comparer avec ...
      p <- ggplot(d,
                  aes(x=d$Jour,y=d[,3],fill=d$Zone
                  )
      )
      p <- p + geom_boxplot(outlier.shape = NA)
      p <- p + scale_fill_manual(values = alpha(c("red","orange","blue"), .3))
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
      p <- p + scale_fill_manual(values = alpha(c("red","orange","blue"), .3))
      p <- p + ggtitle(bquote(atop(.(gg_title), atop(italic(.(subtitle)), "")))) + xlab("Jours") + ylab(fonc)
      print(p)
      #dev.off()
    }
  }
  #
}

# Dans un premier temps : fonctions qui génèrent et enregistrent éventuellement des dataframe pour les suivis temporels.
# Dans un deuxième temps : Affichage, comparaisons avec les données répertoriées pour les zones lésées.

suivi_voxels <- function(rat,fonc,slice,sommet,mesure){
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  repertoires <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  
  segtitle <- "" # indique si nécessaire la fonctionnalité utilisée pour la segmentation
  fonc_seg <- 'ADC' # la segmentation avec fonc_seg sera celle utilisée
  
  subtitle <- sprintf("Rat %s, tranche %i",rat,slice)
  
  liste_jr <- liste_jfr[[fonc]]
  jours <- liste_jr[[rat]]
  
  # On crée la dataframe pour le suivi de la fonctionnalité courante
  d <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(d) <- c("x","y",fonc,"Slice","Type","Jour")
  
  for (jour in jours){
    # On crée la dataframe auxiliaire : valeurs calculées au jour d'examen précédent
    if(jour!="00"){
      d.aux <- as.data.frame(cbind(vox_pr[,1:2],vox_pr[,4:5],'Prévisions'))
      colnames(d.aux) <- c("x","y",fonc,"Slice","Type")
    }
    name_cerveau_fonc <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
    #name_cerveau_seg <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc_seg]],fonc_seg,rat,jour,fonc_seg,'dark')
    
    cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
    #cerveau_seg <- read.table(name_cerveau_seg,header=T)
    
    
    # --------- On va délimiter la partie lésée sur l'ensemble des voxels exploités. --------- #
    vox_mes <- cerveau_fonc
    
    # On sélectionne la tranche qui contient les voxels : coordonnée z
    l <- length(vox_mes[,4])
    liste_tr <- rep(FALSE,l)
    liste_tr <- ifelse(vox_mes$Slice==slice,TRUE,liste_tr)
    vox_mes <- vox_mes[liste_tr,]
    
    # On sélectionne les coordonnées des pixels correspondant aux voxels étudiés : coordonnées x et y
    l <- length(vox_mes[,4])
    liste_tr <- rep(FALSE,l)
    xc <- sommet[1]
    yc <- sommet[2]
    for (i in c(1:mesure)){
      for (j in c(1:mesure)){
        ix <- xc + i-1
        iy <- yc + j-1
        liste_tr <- ifelse(ix==vox_mes$x&iy==vox_mes$y,TRUE,liste_tr)
      }
    }
    vox_mes <- vox_mes[liste_tr,]
    
    # On ajoute la colonne "Zone", qui qualifie les voxels lésés.
    d.mes <- as.data.frame(cbind(vox_mes[,1:2],vox_mes[,4:5],'Mesures'))
    colnames(d.mes) <- c("x","y",fonc,"Slice","Type")
    # --------- Dataframe créée --------- #
    
    # --------- On construit une data.frame contenant les valeurs prédites pour la fonctionnalité --------- #
    if (jour=="00"){
      vox_pr <- vox_mes
      d.aux <- as.data.frame(cbind(vox_mes[,1:2],vox_mes[,4:5],'Prévisions'))
      colnames(d.aux) <- c("x","y",fonc,"Slice","Type")
    }
    if (jour != "22"){#"15" pour le rat 30
      #cerveau_pr <- cerveau_les ----> on recommence avec le cerveau_pr calculé à l'étape précédente - boucle sur les jours d'examen.
      # Cerveau lésé pour l'examen suivant : à construire
      liste.nan <- is.na(vox_pr[,4])
      cerveau_pr <- vox_pr[!liste.nan,]
      
      ll_fsucc <- liste_fsucc[[rat]] 
      l_fsucc <- ll_fsucc[[fonc]]
      fsucc <- l_fsucc[[jour]]
      #print(fsucc(1))
      # Fonction utilisée pour prédire les résultats des prochains examens à l'aide de celui du jours courant
      vox_pr[[fonc]] <- sapply(vox_pr[[fonc]],fsucc)
      # Cet objet sera utilisé au prochain tour de boucle
    }
    # --------- Dataframe construite : prédictions ou valeurs au jour 00. --------- #
    
    d.increment <- as.data.frame(rbind(d.mes,d.aux))#,header=T)
    # --------- Assemblage, jour courant --------- #
    num_jour <- jour
    d.increment <- as.data.frame(cbind(d.increment,num_jour))
    colnames(d.increment) <- c("x","y",fonc,"Slice","Type","Jour")
    d <- as.data.frame(rbind(d,d.increment))
  }# dataframe d remplie
  #
  # --> Sortie ?
  #
  return(d)
}

# Pour une sortie en niveaux de gris : valeurs relatives représentées.

niveau_rel_gris <- function(g,min,max){
  p <- 100*(g-min)/(max-min)
  q <- floor(p)
  res <- sprintf("gray%i",q)
  return(res)
}

# Affichage : carré de pixels, prévisions vs mesures. Sortie pdf à faire.

# Option 1 : voxels ou hist, pour afficher les niveaux de gris en carré ou un histogramme des valeurs de la modalité sur le  carré de voxels.
# Option 2 : sortie dans RStudio ou en pdf, dans le répertoire courant.

aff_suivi_voxels <- function(rat,fonc,slice,sommet,mesure,opt_1,opt_2){
  num_jours <- list("00"=0,"03"=3,"08"=8,"15"=15,"22"=22)
  repertoires <- list('ADC'="fonctionnel_gris",# on peut ajouter ici les autres modalités
                      'BVf'="fonctionnel_gris",
                      'CBF'="fonctionnel_gris",
                      'CMRO2'="fonctionnel_gris",
                      'SO2map'="fonctionnel_gris",
                      'T1map'="fonctionnel_gris",
                      'VSI'="fonctionnel_gris",
                      'Anat'="anatomique_gris")
  
  d <- suivi_voxels(rat,fonc,slice,sommet,mesure)
  
  name_cerveau_fonc <- sprintf("%s/%s/%s-J%s-%s-%s-all.dat",repertoires[[fonc]],fonc,rat,"00",fonc,'bg')
  cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
  
  val_max <- max(cerveau_fonc[,4],na.rm=TRUE)
  val_min <- min(cerveau_fonc[,4],na.rm=TRUE)
  
  liste_jr <- liste_jfr[[fonc]]
  jours <- liste_jr[[rat]]
  njours <- length(jours)
  
  if (opt_1=='voxels'){
    if (opt_2=='pdf'){}
    else{
      #get( getOption( "device" ) )()
      plot.new()
      par(mfrow=c(2,njours))
      
      for (type in c('Mesures','Prévisions')){
        for (jour in jours){
          l <- length(d[,3])
          liste_sel <- rep(TRUE,l)
          liste_sel <- ifelse(d$Type==type,liste_sel,FALSE)
          liste_sel <- ifelse(d$Jour==jour,liste_sel,FALSE)
          e <- d[liste_sel,]
          
          # on crée un vecteur des niveaux de gris du carré de voxels
          vol <- mesure*mesure
          col_vec <- rep("",vol)
          
          for (m in c(1:vol)){
            col_vec[m] <- niveau_rel_gris(d[m,3],val_min,val_max)
          }
          
          plot(e$x, e$y,
               col=col_vec,
               pch=20,
               cex=3,#*(1-d.clust$uncertainty)^4, 
               xlab='x', ylab='y',
               main=sprintf("J%s, %s",jour,type)
          )
        }
      }
      title(sprintf("Suivi temporel : rat %s, modalité %s",rat,fonc),outer=TRUE)
    }
  }
  else if (opt_1=='hist'){
    if (opt_2=='pdf'){
      # histogrammes renvoyés en pdf
    }
    else{
      # histogrammes affichés : valeurs des modalités sur le carré.
      #get( getOption( "device" ) )()
      plot.new()
      par(mfrow=c(2,njours))
      
      for (type in c('Mesures','Prévisions')){
        for (jour in jours){
          # on extrait une sous-dataframe e correspondant au type et au jour courants
          l <- length(d[,3])
          liste_sel <- rep(TRUE,l)
          liste_sel <- ifelse(d$Type==type,liste_sel,FALSE)
          liste_sel <- ifelse(d$Jour==jour,liste_sel,FALSE)
          e <- d[liste_sel,]
          
          # on paramètre et on trace l'histogramme
          vol <- mesure*mesure
          d.fonc <- d[,3]
          liste.nan <- is.na(d.fonc)
          d.fonc <- d.fonc[!liste.nan]
          e.fonc <- e[,3]
          FONC.breaks <- seq(min(d.fonc)-0.1*min(d.fonc), max(d.fonc)+0.1*max(d.fonc), length.out=100)
          e.hist <- hist(e.fonc,
                         breaks=FONC.breaks,
                         ylim=c(0,vol),
                         xlab="",
                         ylab="Nombre de voxels",
                         col='grey50',
                         main=type,sub=sprintf("%s, jour %s",fonc,jour)
                         )
        }
      }
      title(sprintf("Suivi temporel : rat %s, modalité %s",rat,fonc),outer=TRUE)
    }
  }
}

