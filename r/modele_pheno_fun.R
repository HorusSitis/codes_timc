### Première volée de fonctions d'évolution : scalaires ; rat 11, une modalité, un jour. ###


# dep et arr sont les deux vecteurs qui contiennent respectivement les antécédents et images qui définissent le polynôme d'interpolation de Lagrange recherché.
# On utilisera des polynômes de degré 4, pour associer les médiane-quartiles-déciles.

trans_box <- function(x,dep,arr){
  res <- 0
  for (j in (1:length(dep))){
    b <- dep[j]
    a_j <- arr[j]
    res_1 <- b
    for (i in (1:length(arr))){
      a <- arr[i]
      if (i!=j){
        res_1 <- res_1*(x-a)/(a_j-a)
      }
    }
    res <- res + res_1
  }
  return(res)
}

# ADC : pour le moment aucune relation avec les autres modalités. #

liste_succ_ADC <- list("00"=list("dep"=c(600,700,725,750,850),"arr"=c(500,750,850,1000,1400)),
                       "03"=list("dep"=c(500,750,850,1000,1400),"arr"=c(450,1000,1200,1450,2100)),
                       "08"=list("dep"=c(450,1000,1200,1450,2100),"arr"=c(800,1200,2100,2300,3300)),
                       "15"=list("dep"=c(300,1200,2100,2300,3300),"arr"=c(800,1800,2200,2550,3300))
                       )

liste_succ_BVf <- list("00"=list("dep"=,"arr"=),
                       "03"=list("dep"=,"arr"=),
                       "08"=list("dep"=,"arr"=),
                       "15"=list("dep"=,"arr"=)
)

liste_succ_CBF <- list("00"=list("dep"=,"arr"=),
                       "03"=list("dep"=,"arr"=),
                       "08"=list("dep"=,"arr"=),
                       "15"=list("dep"=,"arr"=)
)

liste_succ_SO2map <- list("00"=list("dep"=,"arr"=),
                          "03"=list("dep"=,"arr"=),
                          "08"=list("dep"=,"arr"=),
                          "15"=list("dep"=,"arr"=)
                          )

liste_succ_T1map <- list("00"=list("dep"=c(),"arr"=),
                       "03"=list("dep"=,"arr"=),
                       "08"=list("dep"=,"arr"=),
                       "15"=list("dep"=,"arr"=)
)

liste_succ_VSI <- list("00"=list("dep"=c(1,7,12,20,40),"arr"=c(0,8,13,19,35)),
                       "03"=list("dep"=c(0,8,13,19,35),"arr"=c(2,11,15,21,34)),
                       "08"=list("dep"=c(2,11,15,21,34),"arr"=c(3,10,14,20,34)),
                       "15"=list("dep"=c(3,10,14,20,34),"arr"=c(0,6,10,15,29))
)

liste_box_11 <- list("ADC"=liste_succ_ADC ,
                     "BVf"=liste_succ_BVf ,
                     "CBF"=liste_succ_CBF ,
                     "SO2map"=liste_succ_SO2map ,
                     "T1map"=liste_succ_T1map ,
                     "VSI"=liste_succ_VSI
                     )

succ11_ADC_00 <- function(x){
  dep_arr <- liste_succ_ADC[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_ADC_03 <- function(x){
  dep_arr <- liste_succ_ADC[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_ADC_08 <- function(x){
  dep_arr <- liste_succ_ADC[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_ADC_15 <- function(x){
  dep_arr <- liste_succ_ADC[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_ADC_11 <- list("00"=succ11_ADC_00,
                           "03"=succ11_ADC_03,
                           "08"=succ11_ADC_08,
                           "15"=succ11_ADC_15
)


# CBF : même chose que précédemment. #


# VSI : #

succ11_VSI_00 <- function(x){
  dep_arr <- liste_succ_VSI[["00"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_VSI_03 <- function(x){
  dep_arr <- liste_succ_VSI[["03"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_VSI_08 <- function(x){
  dep_arr <- liste_succ_VSI[["08"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

succ11_VSI_15 <- function(x){
  dep_arr <- liste_succ_VSI[["15"]]
  return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
}

liste_fsucc_VSI_11 <- list("00"=succ11_VSI_00,
                           "03"=succ11_VSI_03,
                           "08"=succ11_VSI_08,
                           "15"=succ11_VSI_15
                           )

liste_fsucc_11 <- list("ADC"=liste_fsucc_ADC_11 ,
                       "BVf"='',#liste_fsucc_BVf_11 ,
                       "CBF"='',#liste_fsucc_CBF_11 ,
                       "SO2map"='',#liste_fsucc_SO2map_11 ,
                       "T1map"='',#liste_fsucc_T1map_11 ,
                       "VSI"=liste_fsucc_VSI_11
                       )

liste_fsucc <- list("11"=liste_fsucc_11,
                    "19"='',#liste_fsucc_19,
                    "26"='',#liste_fsucc_26,
                    "30"=''#liste_fsucc_30
                    )


### Deuxième volée de fonctions : dans un répertoire de benjamin_antoine_labo, utilise des versions vectorisées des fonctions précédentes. ###
## Prennent en arguments des coordonnées pour le coin inférieur gauche, et la taille du carré de pixels. ##

# Dans un premier temps : fonctions qui génèrent et enregistrent éventuellement des dataframe pour les suivis temporels.
# Dans un deuxième temps : Affichage, comparaisons avec les données répertoriées pour les zones lésées.

comp_succ_suivi <- function(rat,hemi,fonc,liste_s_slice,opt_1,opt_2){
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
  
  # On crée la dataframe pour le suivi de la fonctionnalité courante
  d <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
  
  tranches <- liste_s_slice[[fonc]]
  subtitle <- sprintf("Rat%s, tranche(s) ",rat)
  for (tr in tranches){
    subtitle <- paste(subtitle,'-',tr)
  }
  
  liste_jr <- liste_jfr[[fonc]]
  jours <- liste_jr[[rat]]
  
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
    
    if (jour != "22"){#"15" pour le rat 30
      cerveau_pr <- cerveau_les
      # Cerveau lésé pour l'examen suivant : à construire
      ll_fsucc <- liste_fsucc[[rat]]
      l_fsucc <- ll_fsucc[[fonc]]
      fsucc <- l_fsucc[[jour]]
      print(fsucc(1))
      # Fonction utilisée pour prédire les résultats des prochains examens à l'aide de celui du jours courant
      cerveau_pr[[fonc]] <- sapply(cerveau_pr[[fonc]],fsucc)
      # Cet objet sera utilisé au prochain tour de boucle
    }
    
    if (jour=="00"){
      d.aux <- as.data.frame(cbind(cerveau_les[,1:2],cerveau_les[,4:5],'Lésion prévue'))
      colnames(d.aux) <- c("x","y",fonc,"Slice","Zone")
    }
    d.increment <- as.data.frame(rbind(d.les,d.aux,d.sain))#,header=T)
    
    num_jour <- jour
    d.increment <- as.data.frame(cbind(d.increment,num_jour))
    colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    d <- as.data.frame(rbind(d,d.increment))
  }# dataframe d remplie
  
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












