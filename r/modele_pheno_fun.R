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

base_donnes_11 <- {
  coeff_aff_succ_11 <- {
    liste_succ_ADC_11 <- list("00"=list("dep"=c(500,700,725,800,1100),"arr"=c(500,750,850,1000,1400)),
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
  
  aff_fsucc_ADC_j <- {
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
  }
  
  liste_fsucc_ADC_11 <- list("00"=succ11_ADC_00,
                             "03"=succ11_ADC_03,
                             "08"=succ11_ADC_08,
                             "15"=succ11_ADC_15
  )
  
  # BVf #
  
  aff_fsucc_BVf_j <- {
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
  }
  
  liste_fsucc_BVf_11 <- list("00"=succ11_BVf_00,
                             "03"=succ11_BVf_03,
                             "08"=succ11_BVf_08,
                             "15"=succ11_BVf_15
  )
  
  # CBF : même chose que précédemment. #
  
  aff_fsucc_CBF_j <- {
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
  }
  
  liste_fsucc_CBF_11 <- list("00"=succ11_CBF_00,
                             "03"=succ11_CBF_03,
                             "08"=succ11_CBF_08,
                             "15"=succ11_CBF_15
  )
  
  # SO2map #
  
  aff_fsucc_SO2map_j <- {
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
  }
  
  liste_fsucc_SO2map_11 <- list("00"=succ11_SO2map_00,
                                "03"=succ11_SO2map_03,
                                "15"=succ11_SO2map_15
  )
  
  # T1map #
  
  aff_fsucc_T1map_j <- {
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
  }
  
  liste_fsucc_T1map_11 <- list("00"=succ11_T1map_00,
                               "03"=succ11_T1map_03,
                               "08"=succ11_T1map_08,
                               "15"=succ11_T1map_15
  )
  
  # VSI : #
  
  aff_fsucc_VSI_j <- {
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
  
  
}


# ------- Rat 19 ------- #

base_donnes_19 <- {
  coeff_aff_succ_19 <- {
    liste_succ_ADC_19 <- list("00"=list("dep"=c(500,700,750,800,1000),"arr"=c(550,750,800,950,1150)),
                              "03"=list("dep"=c(550,750,800,950,1150),"arr"=c(800,1050,1300,1550,2300)),
                              "08"=list("dep"=c(800,1050,1300,1550,2300),"arr"=c(750,1450,1950,2400,3700)),
                              "15"=list("dep"=c(750,1450,1950,2400,3700),"arr"=c(500,1800,2400,2800,4000))
    )
    
    liste_succ_BVf_19 <- list("00"=list("dep"=c(0,1,1.5,2.5,4.5),"arr"=c(0,2.5,4,5,9.5)),
                              "08"=list("dep"=c(0,2.5,4,5,9.5),"arr"=c(0,3,5,8.5,16)),
                              "15"=list("dep"=c(0,3,5,8.5,16),"arr"=c(0,2.5,4,6.5,13))
    )
    
    liste_succ_CBF_19 <- list("00"=list("dep"=c(0,10,15,20,40),"arr"=c(0,100,140,190,300)),
                              "03"=list("dep"=c(0,100,140,190,300),"arr"=c(0,90,145,220,420)),
                              "08"=list("dep"=c(0,90,145,220,420),"arr"=c(0,40,80,140,300)),
                              "15"=list("dep"=c(0,40,80,140,300),"arr"=c(0,40,105,280,460))
    )
    
    liste_succ_SO2map_19 <- list("00"=list("dep"=c(0,23,40,60,95),"arr"=c(12,48,60,70,100)),
                                 "08"=list("dep"=c(12,48,60,70,100),"arr"=c(25,52,77,88,100)),
                                 "15"=list("dep"=c(25,52,77,88,100),"arr"=c(13,55,70,80,100))
    )
    
    liste_succ_CMRO2_19 <- list("00"=list("dep"=c(0,1,2,4,6),"arr"=c(0,6,12,18,37)),
                                "08"=list("dep"=c(0,6,12,18,37),"arr"=c(0,1,2,5,8)),
                                "15"=list("dep"=c(0,1,2,5,8),"arr"=c(0,2,5,12,25))
    )
    
    liste_succ_VSI_19 <- list("00"=list("dep"=c(0,1.5,2,3,4),"arr"=c(0,3,4,5,9)),
                              "08"=list("dep"=c(0,3,4,5,9),"arr"=c(0,3,5,7,14)),
                              "15"=list("dep"=c(0,3,5,7,14),"arr"=c(0,2.5,4.5,7.5,15))
    )
  }
  
  liste_box_19 <- list("ADC"=liste_succ_ADC_19 ,
                       "BVf"=liste_succ_BVf_19 ,
                       "CBF"=liste_succ_CBF_19 ,
                       "SO2map"=liste_succ_SO2map_19 ,
                       "CMRO2"=liste_succ_CMRO2_19 ,
                       "VSI"=liste_succ_VSI_19
  )
  
  # ADC : pour le moment aucune relation avec les autres modalités. #
  
  aff_fsucc_ADC_j_19 <- {
    succ19_ADC_00 <- function(x){
      dep_arr <- liste_succ_ADC_19[["00"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_ADC_03 <- function(x){
      dep_arr <- liste_succ_ADC_19[["03"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_ADC_08 <- function(x){
      dep_arr <- liste_succ_ADC_19[["08"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_ADC_15 <- function(x){
      dep_arr <- liste_succ_ADC_19[["15"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
  }
  
  liste_fsucc_ADC_19 <- list("00"=succ19_ADC_00,
                             "03"=succ19_ADC_03,
                             "08"=succ19_ADC_08,
                             "15"=succ19_ADC_15
  )
  
  # BVf #
  
  aff_fsucc_BVf_j19 <- {
    succ19_BVf_00 <- function(x){
      dep_arr <- liste_succ_BVf_19[["00"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_BVf_08 <- function(x){
      dep_arr <- liste_succ_BVf_19[["08"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_BVf_15 <- function(x){
      dep_arr <- liste_succ_BVf_19[["15"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
  }
  
  liste_fsucc_BVf_19 <- list("00"=succ19_BVf_00,
                             #"03"=succ19_BVf_03,
                             "08"=succ19_BVf_08,
                             "15"=succ19_BVf_15
  )
  
  # CBF : même chose que précédemment. #
  
  aff_fsucc_CBF_j_19 <- {
    succ19_CBF_00 <- function(x){
      dep_arr <- liste_succ_CBF_19[["00"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_CBF_03 <- function(x){
      dep_arr <- liste_succ_CBF_19[["03"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_CBF_08 <- function(x){
      dep_arr <- liste_succ_CBF_19[["08"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_CBF_15 <- function(x){
      dep_arr <- liste_succ_CBF_19[["15"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
  }
  
  liste_fsucc_CBF_19 <- list("00"=succ19_CBF_00,
                             "03"=succ19_CBF_03,
                             "08"=succ19_CBF_08,
                             "15"=succ19_CBF_15
  )
  
  # SO2map #
  
  aff_fsucc_SO2map_j19 <- {
    succ19_SO2map_00 <- function(x){
      dep_arr <- liste_succ_SO2map_19[["00"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_SO2map_08 <- function(x){
      dep_arr <- liste_succ_SO2map_19[["08"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_SO2map_15 <- function(x){
      dep_arr <- liste_succ_SO2map_19[["15"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
  }
  
  liste_fsucc_SO2map_19 <- list("00"=succ19_SO2map_00,
                             #"03"=succ19_BVf_03,
                             "08"=succ19_SO2map_08,
                             "15"=succ19_SO2map_15
  )
  
  # CMRO2 #
  
  aff_fsucc_CMRO2_j_19 <- {
    succ19_CMRO2_00 <- function(x){
      dep_arr <- liste_succ_CMRO2_19[["00"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_CMRO2_08 <- function(x){
      dep_arr <- liste_succ_CMRO2_19[["08"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_CMRO2_15 <- function(x){
      dep_arr <- liste_succ_CMRO2_19[["15"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
  }
  
  liste_fsucc_CMRO2_19 <- list("00"=succ19_CMRO2_00,
                               #"03"=succ19_CMRO2_03,
                               "08"=succ19_CMRO2_08,
                               "15"=succ19_CMRO2_15
  )
  
  # VSI : #
  
  aff_fsucc_VSI_j_19 <- {
    succ19_VSI_00 <- function(x){
      dep_arr <- liste_succ_VSI_19[["00"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_VSI_08 <- function(x){
      dep_arr <- liste_succ_VSI_19[["08"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
    succ19_VSI_15 <- function(x){
      dep_arr <- liste_succ_VSI_19[["15"]]
      return(trans_box(x,dep_arr[["dep"]],dep_arr[["arr"]]))
    }
    
  }
  
  liste_fsucc_VSI_19 <- list("00"=succ19_VSI_00,
                             #"03"=succ19_VSI_03,
                             "08"=succ19_VSI_08,
                             "15"=succ19_VSI_15
  )
  # Toutes fonctionnalités confondues
  
  liste_fsucc_19 <- list("ADC"=liste_fsucc_ADC_19 ,
                         "BVf"=liste_fsucc_BVf_19 ,
                         "CBF"=liste_fsucc_CBF_19 ,
                         "SO2map"=liste_fsucc_SO2map_19 ,
                         "CMRO2"=liste_fsucc_CMRO2_19 ,
                         "VSI"=liste_fsucc_VSI_19
  )
  
  
}





# ------- Rat 26 ------- #
# ------- Rat 30 ------- #


## ------- Fonctions de transition : tous rats confondus ------- ##

liste_fsucc <- list("11"=liste_fsucc_11,
                    "19"=liste_fsucc_19,
                    "26"='',#liste_fsucc_26,
                    "30"=''#liste_fsucc_30
)

### Deuxième volée de fonctions : dans un répertoire de benjamin_antoine_labo, utilise des versions vectorisées des fonctions précédentes. ###
## Prennent en arguments des coordonnées pour le coin inférieur gauche, et la taille du carré de pixels. ##

# ------------ Etape préliminaire : fonction graphique pour tester des fonctions de transition. Echelle du cerveau entier : lésion, hémisphère contralatéral. ------------ #
# Effectue des calculs, sans renvoyer de dataframe de valeurs prévues.
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
  else if (opt_1=='CBFdark00'){
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
      name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
      name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
      
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
        
        # on crée la liste de niveaux de gris exploitable par density() : hémisphère sain
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
      
      # --------- Densités ; 'aux' : valeurs calculées pour la lésion. --------- #
      dsta <- density(aux)
      if (length(les)!=0){dstl <- density(les)}
      dsts <- density(sain)
      
      #plot.new()
      #par(lend="butt")
      title <- sprintf("%s, J%s segmentation %s",fonc,jour,fonc_seg)
      #title <- paste(title,subtitle,segtitle)
      
      if (opt_3=='pdf'){# on imprime dans un fichier
        
        if (length(les)!=0){
          pdf(file = sprintf("R%s/%s/%s_suivi_dens_vol%s_%s-%s.pdf",rat,"modele_ph",num_rat,fonc_seg,fonc,jour))
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

# ------------ Renvoie une dataframe constituée des valeurs prévues par le modèle, sur n carré de voxels. ------------ #
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
    name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,jour,fonc,'bg')
    #name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,epertoires[[fonc_seg]],fonc_seg,rat,"00",fonc_seg,'dark')
    
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

# ------------ Affichage : carré de pixels, prévisions vs mesures. Sortie pdf à faire. Utilise une dataframe qui contient les résultats des calculs. ------------ #

# Option 1 : voxels ou hist, pour afficher les niveaux de gris en carré ou un histogramme des valeurs de la modalité sur le carré de voxels.
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
  
  name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,repertoires[[fonc]],fonc,rat,"00",fonc,'bg')
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
      title(sprintf("Mesures vs prévisions : rat %s, modalité %s",rat,fonc),outer=TRUE)
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
      title(sprintf("Mesures vs prévisions : rat %s, modalité %s",rat,fonc),outer=TRUE)
    }
  }
}

#################################################################################
###----------------- Modèle phénomonologique : rat numéro 19 -----------------###
#################################################################################

# ------------ Construction de deux bases de données : Jour initial 00 ou 08, toutes fonctionnalités confondues. ------------ #
### Stockage dans .../R19/automate_2 3demi. ###

# Option : liste des modalités pour lesquelles on effectue le remplissage gaussien.

#liste_fonc_modele <- list('CBF','CMRO2','SO2map','BVf','VSI','ADC')
cerveau_multipar <- function(rat,automate,liste_tranches,liste_multipar,opt){
  # Chaque modèle a son jour de départ
  if (automate=="automate_2"){
    jour <- "00"
  }
  else{
    jour <- "08"
  }
  fonc <- liste_fonc[1]
  
  # On initialise avec la première table, à laquelle s'adjoindront les colonnes correspondant aux autres modalités.
  nom_table_ini <- sprintf('R%s/%s/%s/%s-J%s-%s-bg-all.dat',rat,"fonctionnel_gris",fonc,rat,jour,fonc)
  d <- read.table(nom_table_ini,header=T)
  # On ne garde que les coordonnées spatiales
  d <- d[,-c(4)]
  # On ne garde que les tranches d'intérêt
  l <- length(d[,4])
  liste_lignes <- rep(FALSE,l)
  for (tr in liste_tranches){
    liste_lignes <- ifelse(d$Slice==tr,TRUE,liste_lignes)
  }
  d <- d[liste_lignes,]
  
  # On intersecte les cerveaux à chaque itération et on ajoute la colonne correspondant à fonc.
  for (fonc in liste_multipar){
    nom_table_increment <- sprintf('R%s/%s/%s/%s-J%s-%s-bg-all.dat',rat,"fonctionnel_gris",fonc,rat,jour,fonc)
    e <- read.table(nom_table_increment,header=T)
    
    
    if (any(opt==fonc)){
      # On crée la base de données correspondant à la lésion dont on va compléter les valeurs.
      name_lesCBF <- sprintf('R%s/%s/%s/%s-J%s-%s-dark-all.dat',rat,"fonctionnel_gris",'CBF',rat,jour,'CBF')
      cerveau_lesCBF <- read.table(name_lesCBF,header=T)
      
      # On délimite la lésion pour la fonctionnalit\'e courante
      l <- length(e[,4])
      liste_lignes <- rep(FALSE,l)
      m <- length(cerveau_lesCBF[,4])
      for (i in 1:m){
        ix <- cerveau_lesCBF[i,1]
        iy <- cerveau_lesCBF[i,2]
        sli <- cerveau_lesCBF[i,5]
        # On garde seulement, dans d et pour chaque tranche, les coordonnées disponibles sur l'image fonc, enregistrées dans e.
        liste_lignes <- ifelse(ix==e$x&iy==e$y&sli==e$Slice,TRUE,liste_lignes)
        # Luxe en ce qui concerne les tranches, dans le cas du suivi de la tranche 9 du rat 19 : liste_tranches est de longueur 1 concrètement.
      }
      lesion3d <- e[liste_lignes,]
      
      # On crée la future base de données complétée pour la lésion : image contractile
      f <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(f) <- c('x','y','z','CBF','Slice')
      
      for (tranche in liste_tranches){
        lesion <- lesion3d[lesion3d$Slice==tranche,]
        
        liste.nan <- is.na(lesion[[fonc]])
        les_val <- lesion[[fonc]][!liste.nan]
        
        l <- length(lesion[[fonc]])
        norm_l <- rnorm(l)*sqrt(var(les_val))+rep(mean(les_val),l)
        norm_l_0 <- ifelse(norm_l>=0,norm_l,0)
        #print(norm_l_0[1:10])
        
        lesion[[fonc]] <- ifelse(is.na(lesion[[fonc]]),norm_l_0,lesion[[fonc]])
        # Les valeurs manquantes correspondant à la lésion sur la tranche courante sont remplacées par les composantes d'un vecteur gaussien de mêmes espérance et écart-type que l'ensemble des valeurs disponibles.
        
        #print(lesion[1:5,])
        f <- as.data.frame(rbind(f,lesion))
        colnames(f) <- c('x','y','z','CBF','Slice')
      }
      l <- length(e[,4])
      liste_lignes <- rep(FALSE,l)
      m <- length(f[,4])
      for (i in 1:m){
        ix <- f[i,1]
        iy <- f[i,2]
        sli <- f[i,5]
        # On garde seulement, dans d et pour chaque tranche, les coordonnées disponibles sur l'image fonc, enregistrées dans e.
        liste_lignes <- ifelse(ix==e$x&iy==e$y&sli==e$Slice,TRUE,liste_lignes)
        # Luxe en ce qui concerne les tranches, dans le cas du suivi de la tranche 9 du rat 19 : liste_tranches est de longueur 1 concrètement.
      }
      e[liste_lignes,] <- f
    }
    
    # On enlève les (dernières) valeurs manquantes.
    liste.nan <- is.na(e[,4])
    e <- e[!liste.nan,]
    
    # On garde seulement les lignes de la dataframe déjà créée...
    l <- length(d[,1])
    liste_lignes <- rep(FALSE,l)
    m <- length(e[,1])
    for (i in 1:m){
      ix <- e[i,1]
      iy <- e[i,2]
      sli <- e[i,5]
      # On garde seulement, dans d et pour chaque tranche, les coordonnées disponibles sur l'image fonc, enregistrées dans e.
      liste_lignes <- ifelse(ix==d$x&iy==d$y&sli==d$Slice,TRUE,liste_lignes)
      # Luxe en ce qui concerne les tranches, dans le cas du suivi de la tranche 9 du rat 19 : liste_tranches est de longueur 1 concrètement.
    }
    d <- d[liste_lignes,]
    # ... et seulement celles correspondantes à la modalité courante..
    l <- length(e[,1])
    liste_lignes <- rep(FALSE,l)
    m <- length(d[,1])
    for (i in 1:m){
      ix <- d[i,1]
      iy <- d[i,2]
      sli <- d[i,4]
      # On garde seulement, dans d et pour chaque tranche, les coordonnées disponibles sur l'image fonc, enregistrées dans e.
      liste_lignes <- ifelse(ix==e$x&iy==e$y&sli==e$Slice,TRUE,liste_lignes)
      # Luxe en ce qui concerne les tranches, dans le cas du suivi de la tranche 9 du rat 19 : liste_tranches est de longueur 1 concrètement.
    }
    e <- e[liste_lignes,]
    
    colnames.inc <- colnames(d)
    d <- as.data.frame(cbind(d,e[,4]))
    colnames(d) <- c(colnames.inc,fonc)
  }
  
  # Il reste à qualifier l'état initial : on utilise pour cela une fonction vectorisée.
  if (automate=="automate_2"){
    l <- length(d[,1])
    etat <- rep('per',l)
    # on segmente la lésion en CBFdark00
    cerveau_seg <- read.table(sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,"fonctionnel_gris",'CBF',rat,"00",'CBF','dark'),header=T)
    m <- length(cerveau_seg[,4])
    for (i in c(1:m)){
      ix <- cerveau_seg[i,1]
      iy <- cerveau_seg[i,2]
      sli <- cerveau_seg[i,5]
      etat <- ifelse(ix==d$x&iy==d$y&sli==d$Slice,'cyto',etat)
    }
  }
  else{
    # pour un stade vasogénique : par défaut, l'état est perturbé. Les fluctuations autour d'un équilibre arrivent, en principe, rapidement.
    l <- length(d[,1])
    etat <- rep('per',l)
    # on segmente la lésion en CBFdark00
    cerveau_seg <- read.table(sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,"fonctionnel_gris",'CBF',rat,"00",'CBF','dark'),header=T)
    m <- length(cerveau_seg[,4])
    for (i in c(1:m)){
      ix <- cerveau_seg[i,1]
      iy <- cerveau_seg[i,2]
      sli <- cerveau_seg[i,5]
      etat <- ifelse(ix==d$x&iy==d$y&sli==d$Slice,'les_1_deb',etat)
    }
    # on utilise le résultat de la clusterisation CBF, c(1,2), au jour 08
    etat <- ifelse(etat=='les_1_deb'&d$CBF>100,'les_2_deb',etat)
  }
  colnames.inc <- colnames(d)
  d <- as.data.frame(cbind(d,etat),stringAsFactors=TRUE)
  colnames(d) <- c(colnames.inc,'Etat')
  
  # Sorties : dataframe, qui est enregistrée dans le dossier automate_2_3demi. On mentionne les rtanches retenues.
  nom_tranches <- ""
  for (tr in liste_tranches){
    nom_tranches <- paste(nom_tranches,tr,sep="-")
  }
  if (length(opt)==0){
    nom_table_res <- sprintf('R%s/%s/cerveau%s_multi_J%s%s_Slices%s.dat',rat,automate,rat,jour,"",nom_tranches)
  }
  else{
    nom_table_res <- sprintf('R%s/%s/cerveau%s_multi_J%s%s_Slices%s.dat',rat,automate,rat,jour,"_fill_norm",nom_tranches)
  }
  write.table(d, nom_table_res, row.names=F, quote=F, sep='\t')
  #return(d)
}

### Affichage des états : tranche 9 du rat 19. ###

# Option 1 : on utilise ou non une base de données avec des valeurs manquantes complétées.
# Option 2 : sortie. pdf si 'pdf', affichage dans RStudio sinon.

#list.color=list('cyto'='darkred','equi'='cyan','les_1_deb=red','les_1_fin'='brown','les_2_deb'='gold','les_2_fin'='orange','per'='blue')
couleur_etat <- function(etat){
  if (etat=='cyto'){
    color <- 'darkred'
  }
  else if (etat=='per'){
    color <- 'blue'
  }
  else if (etat=='les_1_deb'){
    color <- 'red'
  }
  else if (etat=='les_1_fin'){
    color <- 'darkmagenta'
  }
  else if (etat=='les_2_deb'){
    color <- 'darkorange3'
  }
  else if (etat=='les_2_fin'){
    color <- 'gold'
  }
  else{
    color <- 'grey30'
  }
  return(color)
}

affichage_etats_cerveau <- function(rat,automate,liste_tranches,jour,opt_1,opt_2){
  if (length(opt_1)==0){
    nom_opt <- ""
  }
  else{
    nom_opt <- "_fill_norm"
  }
  # Nom du nombre de tranches, numéro du modèle
  nom_tranches <- ""
  for (tr in liste_tranches){
    nom_tranches <- paste(nom_tranches,tr,sep="-")
  }
  if (automate=="automate_2"){
    nom_mod <- "2"
  }
  else{
    nom_mod <- "3demi"
  }
  
  nom_table <- sprintf("R%s/%s/cerveau%s_multi_J%s%s_Slices%s.dat",rat,automate,rat,jour,nom_opt,nom_tranches)
  d <- read.table(nom_table,header=T)
  
  d$Etat <- sapply(d$Etat,couleur_etat)
  
  if (opt_2=='pdf'){
    if (length(opt_1)>0){
      suffixe <- "Com"
    }
    else{
      suffixe <- "Brut"
    }
    
    nom_fichier <- sprintf("R%s/modele_ph/%s-J%s-modele%s_sim%s.pdf",rat,rat,jour,nom_mod,suffixe)
    pdf(file=nom_fichier)
    
    par(mfrow=c(1,1))
    scatterplot3d(d$x,
                  d$y,
                  d$Slice,
                  #cerveau.clust$z, 
                  color= d$Etat,
                  pch=20,
                  xlab='x',
                  ylab='y',
                  zlab='Slice',#'z',
                  lab.z=1+d$Slice[length(d$Slice)]-d$Slice[1],
                  main=sprintf("Modèle %s : rat %s, jour %s",nom_mod,rat,jour),
                  sub = sprintf("Tranches %s",nom_tranches)
    )
    
    dev.off()
  }
  else{
    #plot.new()
    par(mfrow=c(1,1))
    scatterplot3d(d$x,
                  d$y,
                  d$Slice,
                  #cerveau.clust$z, 
                  color= d$Etat,
                  pch=20,
                  xlab='x',
                  ylab='y',
                  zlab='Slice',#'z',
                  lab.z=1+d$Slice[length(d$Slice)]-d$Slice[1],
                  main=sprintf("Modèle %s : rat %s, jour %s",nom_mod,rat,jour),
                  sub = sprintf("Tranches %s",nom_tranches)
    )
    #title(sprintf("Cerveau modélisé : rat %s, %s, jour %s, modèle %f",rat,automate,jour,mod),outer=TRUE)
  }
}


# ------------ Suivi pour ... un pixel, carré, cerveau entier ? ------------ #

## - Coefficients pour les équations - ##

tau_cont <- function(delta_ini,var_seuil,n_jours){
  return(1/n_jours*log(var_seuil/(delta_ini)))
}

q_dis <- function(delta_ini,var_seuil,n_jours){
  return((var_seuil/delta_ini)^(1/n_jours))
}




cerveau_dyn_2 <- function(cerveau){
  #
}


cerveau_dyn_3demi <- function(rat,liste_coupes,t){# t : jour final
  jour <- "08"
  num_jour <- 8
  suff <- ""
  for (cp in liste_coupes){
    suff <- paste(suff,"-",cp,sep='')
  }
  # Chargement du cerveau dont les données seront utilisées
  nom_cerveau <- sprintf("R%s/automate_3demi/cerveau%s_multi_J%s_Slices%s.dat",rat,rat,jour,suff)
  cerveau <- read.table(nom_cerveau,header=T)
  l <- length(cerveau$x)
  cerveau <- as.data.frame(cbind(cerveau,rep(num_jour,l)))
  colnames(cerveau) <- c(c('x','y','z','Slice'),vect_fonc,c('Etat','Jour'))
  
  # Constitution de la base de données : constantes du modèle
  ## Distributions asymptotiques
  cbf_asy1 <- rnorm(l,50,40)
  cmro2_asy1 <- rnorm(l,5,1.5)
  so2_asy1 <- rnorm(l,75,25)
  vsi_asy1 <- rnorm(l,4,2)
  cbf_trans3 <- rnorm(l,120,40)
  vsi_asy2 <- rnorm(l,5,2)
  so2_asy2 <- rnorm(l,75,15)
  so2_sain <- rnorm(l,75,10)
  
  cbf_sain <- rnorm(l,130,30)
  
  adc_asy1 <- rnorm(l,3000,1000)
  adc_asy2 <- rnorm(l,2500,700)
  
  ## Coefficients multiplicateurs
  q_cbf_asy1 <- q_dis(50,40,7)#rnorm(l,50,50)
  q_cmro2_asy1 <- q_dis(8,1.5,7)#rnorm(l,5,5.5)
  q_so2_asy1 <- q_dis(35,25,7)#rnorm(l,75,35)
  q_vsi_asy1 <- q_dis(3,2,7)#rnorm(l,4,2.8)
  q_cbf_trans2 <- q_dis(50,40,7)#rnorm(l,120,50)
  q_vsi_asy2 <- q_dis(3,2,7)#rnorm(l,5,4.3)
  q_so2_asy2 <- q_dis(20,15,7)#rnorm(l,75,30)
  q_so2_sain <- q_dis(30,10,7)#rnorm(l,75,25)
  
  q_cbf_sain <- q_dis(50,30,5)#rnorm(l,130,30)
  
  #adc_asy1 <- rnorm(l,3000,1000)
  #adc_asy2 <- rnorm(l,2500,700)
  
  
  const <- as.data.frame(cbind(cbf_asy1,cmro2_asy1,so2_asy1,cbf_trans2,vsi_asy2,so2_asy2,so2_sain,cbf_sain))
  colnames(const)=c("CBF_asy1","CMRO2_asy1","SO2_asy1","CBF_trans2","VSI_asy2","SO2_asy2","SO2_sain","CBF_sain")
  
  # Initialisation pour la simulation : ...
  
  # Evolution à partir du jour 8 : on part de "cerveau", dont les états sont labélisés
  # On concatène successivement les tables de données obtenues à chaque itération de la boucle
  
  for (num_jour in c(9:t)){
    # Increment
    c_inc <- cerveau[cerveau$Jour==num_jour-1,]
    # Accroissement des modalités, en un jour
    c_t <- data.frame(matrix(ncol = 6, nrow = l))
    colnames(c_t) <- c('ADC','BVf','CBF','CMRO2','SO2map','VSI')
    
    # Calcul du bruit pour chaque variable : équations de type 0
    bruit_mod <- data.frame(matrix(ncol = 6, nrow = l))
    colnames(bruit_mod) <- c('ADC','BVf','CBF','CMRO2','SO2map','VSI')
    for (fonc in liste_fonc){
      bruit_mod[[fonc]] <- rnorm(l,0,liste_mod_bruit[[fonc]])
    }
    
    # Calcul, vectorisé et conditionnel, des composantes des vecteurs de cerveau_temp.
    # On commence par l'état d'équilibre sain, perturbé non lésé, lésé 1 début, lésé 1 fin etc.
    # Il existe des niveaux de priorités entre les modalités
    # SO2
    c_t$SO2map <- bruit_mod$SO2map# Etat non lésé perturbé : même chose
    c_t$SO2map <- ifelse(c_inc$Etat=='les_1_deb',(1-q_so2_asy1)*(so2_asy1-c_inc$SO2map),c_t$SO2map)
    c_t$SO2map <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$SO2map,c_t$SO2map)
    c_t$SO2map <- ifelse(c_inc$Etat=='les_2_deb',(1-q_so2_sain)*(so2_sain-c_inc$SO2map),c_t$SO2map)
    c_t$SO2map <- ifelse(c_inc$Etat=='les_2_ifin',(1-q_so2_asy2)*(so2_asy2-c_inc$SO2map),c_t$SO2map)
    
    # CBF
    c_t$CBF <- bruit_mod$CBF
    c_t$CBF <- ifelse(c_inc$Etat=='per',(1-q_cbf_sain)*(cbf_sain-c_inc$CBF),c_t$CBF)
    c_t$CBF <- ifelse(c_inc$Etat=='les_1_deb',(1-q_cbf_asy1)*(cbf_asy1-c_inc$CBF),c_t$CBF)
    c_t$CBF <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$CBF,c_t$CBF)
    c_t$CBF <- ifelse(c_inc$Etat=='les_2_deb',(1-q_cbf_trans2)*(cbf_trans2-c_inc$CBF),c_t$CBF)
    c_t$CBF <- ifelse(c_inc$Etat=='les_2_fin',-10*c_t$SO2map,c_t$CBF)
    
    # ADC
    c_t$ADC <- bruit_mod$ADC
    c_t$ADC <- ifelse(c_inc$Etat=='les_1_deb'|c_inc$Etat=='les_1_fin',tau_1*(c_inc$ADC)^alpha_1*(adc_asy1-c_inc$ADC),c_t$ADC)
    c_t$ADC <- ifelse(c_inc$Etat=='les_2_deb'|c_inc$Etat=='les_2_fin',tau_2*(c_inc$ADC)^alpha_2*(adc_asy2-c_inc$ADC),c_t$ADC)
    
    # CMRO2
    c_t$CMRO2 <- bruit_mod$CMRO2
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_1_deb',(1-q_cmro2_asy1)*(cmro2_asy1-c_inc$CMRO2),c_t$CMRO2)
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$CMRO2,c_t$CMRO2)
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_2_deb',c_t$SO2map,c_t$CMRO2)
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_2_fin',-1*c_t$SO2map,c_t$CMRO2)
    
    # VSI
    c_t$VSI <- bruit_mod$VSI
    c_t$VSI <- ifelse(c_inc$Etat=='les_1_deb',(1-q_vsi_asy1)*(vsi_asy1-c_inc$VSI),c_t$VSI)
    c_t$VSI <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$VSI,c_t$VSI)
    c_t$VSI <- ifelse(c_inc$Etat=='les_2_deb',(1-q_vsi_asy2)*(vsi_asy2-c_inc$VSI),c_t$VSI)
    c_t$VSI <- ifelse(c_inc$Etat=='les_2_fin',bruit_mod$VSI,c_t$VSI)
    
    # BVf
    c_t$BVf <- bruit_mod$BVf
    c_t$BVf <- ifelse(c_inc$Etat=='les_1_deb',0.5*c_t$VSI,c_t$BVf)
    c_t$BVf <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$BVf,c_t$BVf)
    c_t$BVf <- ifelse(c_inc$Etat=='les_2_deb',2*c_t$VSI,c_t$BVf)
    c_t$BVf <- ifelse(c_inc$Etat=='les_2_fin',bruit_mod$BVf,c_t$BVf)
    
    # Calcul de la valeurs courante de cerveau_inc : modalités
    for (fonc in liste_fonc){
      c_inc[[fonc]] <- c_inc[[fonc]]+c_t[[fonc]]
    }
    
    # Les modalités ne prennent pas de valeur négative
    for (fonc in liste_fonc){
      c_inc[[fonc]] <- ifelse(c_inc[[fonc]]<0,0,c_inc[[fonc]])
    }
    
    #print(c_inc$Etat[1:10])
    # Calcul des états de cerveau_inc
    c_inc$Etat <- ifelse((c_inc$Etat=='per'&(abs(c_inc$CBF-cbf_sain)<0.5*abs(bruit_mod$CBF))),'sain',as.character(c_inc$Etat))
    #c_inc$Etat <- ifelse((c_inc$Etat=='les_1_deb'&(abs(c_inc$CBF-cbf_sain)<0.5*abs(bruit_mod$CBF))),'les_1_fin',as.character(c_inc$Etat))
    c_inc$Etat <- ifelse(c_inc$Etat=='les_2_deb'&(abs(c_inc$SO2map-so2_sain)<0.5*abs(bruit_mod$SO2map)),'les_2_fin',as.character(c_inc$Etat))
    #print(c_inc$Etat[1:10])
    
    # Jour courant pour cerveau_inc
    c_inc$Jour <- rep(num_jour,l)
    
    # Concaténation
    cerveau <- as.data.frame(rbind(cerveau,c_inc))
    colnames(cerveau) <- c(c('x','y','z','Slice'),vect_fonc,c('Etat','Jour'))
    
    # Fin de la boucle
  }
  #return(cerveau)
  nom_table_dyn <- sprintf("R%s/automate_3demi/cerveau%s_multi_dyn_Slices%s_fin%i.dat",rat,rat,suff,t)
  write.table(cerveau, nom_table_dyn, row.names=F, quote=F, sep='\t')
}

# Option 1 : automate choisi comme modèle
# Option 2 : remlissage gaussien, dans le cas de l'automate 2
# Option 3 : densités ou diagrammes en boîte pour la sortie
# Option 4 : sortie. pdf si 'pdf', affichage dans RStudio sinon.

suivi_mod <- function(rat,liste_mod,liste_coupes,t,opt_1,opt_2,opt_3,opt_4){
  
  # Nommage de la base de données, constituée à l'aide d'une fonction cerveau_dyn
  coupes <- ""
  for (cp in liste_coupes){coupes <- paste(coupes,"-",cp,sep='')}
  if (opt_1=="automate_2"){
    jour <- "00"
  }
  else{
    jour <- "08"
  }
  nom_cerveaux <- sprintf("R%s/%s/cerveau%s_multi_dyn%s_Slices%s_fin%s.dat",rat,opt_1,rat,opt_2,coupes,t)
  # Importation des données calculées à l'aide du modèle opt_1
  cerveaux <- read.table(nom_cerveaux,header=T)
  
  if (opt_3=='box'){
    # On change les composantes de cerveaux_Jour en chaînes de caractères, pour que l'ensemble des abscisses du futur ggplot soit discret.
    # On ajoute un "0" devant les numéros de jours inférieurs à 10, pour que les graduations temporelles soient dans le bon ordre.
    Jour10 <- cerveaux$Jour < 10
    cerveaux$Jour <- as.character(cerveaux$Jour)
    cerveaux$Jour <- ifelse(Jour10,paste("0",cerveaux$Jour,sep=''),cerveaux$Jour)
    
    for (fonc in liste_mod){
      if (opt_4=='pdf'){
        file_name <- sprintf("R%s/%s/cerveau%s_%s%s_dyn%s_Slices%s_fin%s.pdf",rat,opt_1,rat,fonc,opt_3,opt_2,coupes,t)
        
        gg_title <- sprintf("Partie lésée 1, jours %s à %i",jour,t)
        subtitle <- sprintf("Coupes %s, modalité %s",coupes,fonc)
        cerveaux_i1 <- cerveaux[cerveaux$Etat=='les_1_deb'|cerveaux$Etat=='les_1_fin',]
        pi1 <- ggplot(cerveaux_i1,
                      aes(x=cerveaux_i1$Jour,y=cerveaux_i1[[fonc]],fill=cerveaux_i1$Etat
                      )
        )
        pi1 <- pi1 + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        pi1 <- pi1 + scale_fill_manual(values = alpha(c("red","darkmagenta"), .3))
        pi1 <- pi1 + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("Partie lésée 2, jours %s à %i",jour,t)
        subtitle <- sprintf("Coupes %s, modalité %s",coupes,fonc)
        cerveaux_i2 <- cerveaux[cerveaux$Etat=='les_2_deb'|cerveaux$Etat=='les_2_fin',]
        pi2 <- ggplot(cerveaux_i2,
                      aes(x=cerveaux_i2$Jour,y=cerveaux_i2[[fonc]],fill=cerveaux_i2$Etat
                      )
        )
        pi2 <- pi2 + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        pi2 <- pi2 + scale_fill_manual(values = alpha(c("darkorange3","gold"), .3))
        pi2 <- pi2 + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("Partie saine, jours %s à %i",jour,t)
        subtitle <- sprintf("Coupes %s, modalité %s",coupes,fonc)
        cerveaux_s <- cerveaux[cerveaux$Etat=='per'|cerveaux$Etat=='sain',]
        ps <- ggplot(cerveaux_s,
                     aes(x=cerveaux_s$Jour,y=cerveaux_s[[fonc]],fill=cerveaux_s$Etat
                     )
        )
        ps <- ps + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        ps <- ps + scale_fill_manual(values = alpha(c("blue","grey30"), .3))
        ps <- ps + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        p <- grid.arrange(pi1,pi2,ps,ncol=1,nrow=3)
        
        # On imprime pour que le graphique soit exporté
        print(p)
        dev.off()
        
        # On enregistre le fichier .pdf
        ggsave(filename=file_name,plot = p)
      }
      else{
        gg_title <- sprintf("Evolution de %s, jours %s à %i",fonc,jour,t)
        subtitle <- sprintf("Coupes %s, modalité %s, partie %s %i",coupes,fonc,"ischémiée",1)
        cerveaux_i1 <- cerveaux[cerveaux$Etat=='les_1_deb'|cerveaux$Etat=='les_1_fin',]
        pi1 <- ggplot(cerveaux_i1,
                    aes(x=cerveaux_i1$Jour,y=cerveaux_i1[[fonc]],fill=cerveaux_i1$Etat
                    )
        )
        pi1 <- pi1 + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        pi1 <- pi1 + scale_fill_manual(values = alpha(c("red","darkmagenta"), .3))
        pi1 <- pi1 + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("Evolution de %s, jours %s à %i",fonc,jour,t)
        subtitle <- sprintf("Coupes %s, modalité %s, partie %s %i",coupes,fonc,"ischémiée",2)
        cerveaux_i2 <- cerveaux[cerveaux$Etat=='les_2_deb'|cerveaux$Etat=='les_2_fin',]
        pi2 <- ggplot(cerveaux_i2,
                      aes(x=cerveaux_i2$Jour,y=cerveaux_i2[[fonc]],fill=cerveaux_i2$Etat
                      )
        )
        pi2 <- pi2 + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        pi2 <- pi2 + scale_fill_manual(values = alpha(c("darkorange3","gold"), .3))
        pi2 <- pi2 + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("Evolution de %s, jours %s à %i",fonc,jour,t)
        subtitle <- sprintf("Coupes %s, modalité %s, partie %s",coupes,fonc,"saine")
        cerveaux_s <- cerveaux[cerveaux$Etat=='per'|cerveaux$Etat=='sain',]
        ps <- ggplot(cerveaux_s,
                     aes(x=cerveaux_s$Jour,y=cerveaux_s[[fonc]],fill=cerveaux_s$Etat
                     )
        )
        ps <- ps + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        ps <- ps + scale_fill_manual(values = alpha(c("blue","grey30"), .3))
        ps <- ps + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        p <- grid.arrange(pi1,pi2,ps,ncol=1,nrow=3)
        print(p)
      }
    }
  }
  else if (opt_3=='dens'){
    #
  }
}

# Option 1 : remlissage gaussien, dans le cas de l'automate 2
# Option 2 : densités ou diagrammes en boîte pour la sortie
# Option 3 : sortie. pdf si 'pdf', affichage dans RStudio sinon.

# liste_t choisit l'automate utilisé : celui dont la condition initiale est issue de l'expérieuce au premier jour de la liste.

comp_mod_clust <- function(rat,hemi,liste_mod,liste_coupes,liste_t,opt_1,opt_2,opt_3){
  
  # Nommage de la base de données, constituée à l'aide d'une fonction cerveau_dyn
  coupes <- ""
  for (cp in liste_coupes){coupes <- paste(coupes,"-",cp,sep='')}
  
  if (liste_t[1]==0){
    modele <- "automate_2"
  }
  else{
    modele <- "automate_3demi"
    # Dans ce cas la valeur opt_1 est triviale.
  }
  t <- liste_t[length(liste_t)]
  
  # Importation des données calculées à l'aide du modèle opt_1
  nom_cerveaux <- sprintf("R%s/%s/cerveau%s_multi_dyn%s_Slices%s_fin%s.dat",rat,modele,rat,opt_1,coupes,t)
  cerveaux <- read.table(nom_cerveaux,header=T)
  
  # On ne garde que les prédictions pour les jours de liste_t
  l <- length(cerveaux$Jour)
  Jours <- rep(FALSE,l)
  for (j in liste_t){
    Jours <- ifelse(j==cerveaux$Jour,TRUE,Jours)
  }
  cerveaux <- cerveaux[Jours,]
  
  # Construction de la base de données d : données de l'expérience, clusterisées, modalité 'CBF'. Pas d'hémisphère sain.
  cer_clust_cbf <- {
    fonc <- 'CBF'
    # Initialisation
    d <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    jour <- "08"
    
    cl_min <- 1
    cl_max <- 2
    
    name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,"fonctionnel_gris",fonc,rat,jour,fonc,'bg')
    name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,"fonctionnel_gris",'CBF',rat,"00",'CBF','dark')
    
    cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
    cerveau_seg <- read.table(name_cerveau_seg,header=T)
    
    cerveau_les <- cerveau_fonc
    
    l <- length(cerveau_les[,4])
    liste_tr <- rep(FALSE,l)
    
    # Segmentation du volume lésé, à l'aide de fonc_seg
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
    
    l.clust <- cluster_jfr_fmin(cerveau_les,cl_min,cl_max,5)
    #les.clust <- cbind(cerveau_les,l.clust$classification)
    
    d.les$Zone <- ifelse(l.clust$classification >1,'Lesion 2','Lesion 1')
    
    #print(length(d.les$Zone))
    # --------- Caractérisation sur la dataframe : volume lésé clusterisé --------- #
    
    d.increment <- d.les#as.data.frame(rbind(d.les,d.sain))#,header=T)
    # données pour le jour courant, fonctionnalité fonc
    #num_jour <- jour#num_jours[[jour]]
    d.increment <- as.data.frame(cbind(d.increment,jour))
    #print(length(d.increment$Zone))
    colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
    
    d <- as.data.frame(rbind(d,d.increment))
  }
  d_fcl <- d
  
  if (opt_2=='box'){
    # On change les composantes de cerveaux_Jour en chaînes de caractères, pour que l'ensemble des abscisses du futur ggplot soit discret.
    # On ajoute un "0" devant les numéros de jours inférieurs à 10, pour que les graduations temporelles soient dans le bon ordre.
    Jour10 <- cerveaux$Jour < 10
    cerveaux$Jour <- as.character(cerveaux$Jour)
    cerveaux$Jour <- ifelse(Jour10,paste("0",cerveaux$Jour,sep=''),cerveaux$Jour)
    
    for (fonc in liste_mod){
      
      # Dataframe contenant les résultats des expériences pour tous les jours de litse_t 
      cer_exp_fonc <- {
        # On crée la dataframe pour le suivi de la fonctionnalioté courante
        d <- data.frame(matrix(ncol = 6, nrow = 0))
        colnames(d) <- c("x","y",fonc,"Slice","Zone","Jour")
        
        tranches <- liste_coupes
        jours10 <- liste_t < 10
        jours <- as.character(liste_t)
        jours <- ifelse(jours10,paste("0",jours,sep=''),jours)
        
        for (jour in jours){
          name_cerveau_fonc <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,"fonctionnel_gris",fonc,rat,jour,fonc,'bg')
          name_cerveau_seg <- sprintf("R%s/%s/%s/%s-J%s-%s-%s-all.dat",rat,"fonctionnel_gris",'CBF',rat,"00",'CBF','dark')
          
          cerveau_fonc <- read.table(name_cerveau_fonc,header=T)
          cerveau_seg <- read.table(name_cerveau_seg,header=T)
          
          # Sélection de l'hémisphère sain au jour courant.
          block <- {
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
            d.sain <- as.data.frame(cbind(tranches_hem[,1:2],tranches_hem[,4:5],"Hem sain courant",stringsAsFactors=FALSE))#sprintf('Hem sain J%s',jour))
            colnames(d.sain) <- c("x","y",fonc,"Slice","Zone")
          }
          # --------- Dataframe créée : hémisphère sain --------- #
          block <- {
            cerveau_les <- cerveau_fonc
            
            l <- length(cerveau_les[,4])
            liste_tr <- rep(FALSE,l)
            
            # Segmentation du volume lésé, à l'aide de fonc_seg : ADC, CBF ...
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
          }
          # --------- Dataframe créée : volume lésé --------- #
          
          # --------- On traduit la clusterisation effectuée avec fonc_cl, en appelant d_fcl  --------- #
          block <- {
            l <- length(d_fcl$Jour)
            liste_fclj <- rep(FALSE,l)
            if (opt_2==''){
              jour_cl <- jour
            }
            else{
              jour_cl <- "08"
            }
            liste_fclj <- ifelse(d_fcl$Jour==jour_cl,TRUE,liste_fclj)
            d.fcl.jour <- d_fcl[liste_fclj,]
            
            # --- Intersection entre les volumes lésés au jour courant : fonc et fonc_cl --- #
            
            # - Première inclusion - #
            l <- length(d.les$Zone)
            liste_xytr <- rep(FALSE,l)
            m <- length(d.fcl.jour$Zone)
            for (i in c(1:m)){
              ix <- d.fcl.jour[i,1]
              iy <- d.fcl.jour[i,2]
              sli <- d.fcl.jour[i,4]
              # On garde seulement, dans d.les et pour chaque tranche segmentée, les coordonnées de la zone d'intérêt pour fonc_cl
              liste_xytr <- ifelse(ix==d.les$x&iy==d.les$y&sli==d.les$Slice,TRUE,liste_xytr)
              # Seules les tranches de ... sont présentes dans d_les : pas besoin de trier les tranches du cerveau_les obtenu
            }
            d.les <- d.les[liste_xytr,]
            
            # - Deuxième inclusion - #
            l <- length(d.fcl.jour$Zone)
            liste_xytr <- rep(FALSE,l)
            m <- length(d.les$Zone)
            for (i in c(1:m)){
              ix <- d.les[i,1]
              iy <- d.les[i,2]
              sli <- d.les[i,4]
              # On garde seulement, dans d.fcl.jour et pour chaque tranche segmentée, les coordonnées de la lésion, vue depuis la fonctionnalité courante.
              liste_xytr <- ifelse(ix==d.fcl.jour$x&iy==d.fcl.jour$y&sli==d.fcl.jour$Slice,TRUE,liste_xytr)
              # Seules les tranches de ... sont présentes dans d_fcl : pas besoin de trier les tranches du cerveau_les obtenu
            }
            d.fcl.jour <- d.fcl.jour[liste_xytr,]
            
            d.les$Zone <- ifelse(d.fcl.jour$Zone=='Lesion 2','Lesion 2','Lesion 1')
          }
          # --------- Caractérisation sur la dataframe : volume lésé clusterisé avec fonc_cl --------- #
          
          d.increment <- as.data.frame(rbind(d.les,d.sain))#,header=T)
          # données pour le jour courant, fonctionnalité fonc
          d.increment <- as.data.frame(cbind(d.increment,jour))
          colnames(d.increment) <- c("x","y",fonc,"Slice","Zone","Jour")
          
          d <- as.data.frame(rbind(d,d.increment))
        }# dataframe d remplie pour la fonctionnalité courante
      }
      # On donne à cette dataframe une forme adaptée à ggplot.
      liste.nan <- is.na(d[,3])
      d <- d[!liste.nan,]
      
      if (opt_3=='pdf'){
        # plus tard
      }
      else{
        gg_title <- sprintf("Prévisions en %s",fonc)#,liste_t[1],t)
        subtitle <- sprintf("Coupe(s) %s, partie %s %i",coupes,"ischémiée",1)
        cerveaux_i1 <- cerveaux[cerveaux$Etat=='les_1_deb'|cerveaux$Etat=='les_1_fin',]
        pi1 <- ggplot(cerveaux_i1,
                      aes(x=cerveaux_i1$Jour,y=cerveaux_i1[[fonc]],fill=cerveaux_i1$Etat
                      )
        )
        pi1 <- pi1 + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        pi1 <- pi1 + scale_fill_manual(values = alpha(c("red","darkmagenta"), .3))
        pi1 <- pi1 + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("Prévisions en %s",fonc)#,liste_t[1],t)
        subtitle <- sprintf("Coupe(s) %s, partie %s %i",coupes,"ischémiée",2)
        cerveaux_i2 <- cerveaux[cerveaux$Etat=='les_2_deb'|cerveaux$Etat=='les_2_fin',]
        pi2 <- ggplot(cerveaux_i2,
                      aes(x=cerveaux_i2$Jour,y=cerveaux_i2[[fonc]],fill=cerveaux_i2$Etat
                      )
        )
        pi2 <- pi2 + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        pi2 <- pi2 + scale_fill_manual(values = alpha(c("darkorange3","gold"), .3))
        pi2 <- pi2 + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("Prévisions en %s",fonc)#,liste_t[1],t)
        subtitle <- sprintf("Coupe(s) %s, partie %s",coupes,"saine")
        cerveaux_s <- cerveaux[cerveaux$Etat=='per'|cerveaux$Etat=='sain',]
        ps <- ggplot(cerveaux_s,
                     aes(x=cerveaux_s$Jour,y=cerveaux_s[[fonc]],fill=cerveaux_s$Etat
                     )
        )
        ps <- ps + geom_boxplot(outlier.shape = NA, varwidth = TRUE)
        ps <- ps + scale_fill_manual(values = alpha(c("blue","grey30"), .3))
        ps <- ps + ggtitle(bquote(atop(.(gg_title),atop(italic(.(subtitle)), ""))))+ xlab("Jours")+ ylab(fonc)
        
        gg_title <- sprintf("%s : résultats expérimentaux",fonc)#,liste_t[1],t)
        subtitle <- sprintf("Coupe(s) %s, lésion et hémisphère contralatéral",coupes)
        pe <- ggplot(d,
                    aes(x=d$Jour,y=d[,3],fill=d$Zone
                    )
        )
        pe <- pe + geom_boxplot(outlier.shape = NA,varwidth = TRUE)
        pe <- pe + scale_fill_manual(values = alpha(c("cyan","red","orange"), .3))
        pe <- pe + theme(plot.title = element_text(color="black", size=16),
                       plot.subtitle = element_text(color="red", size=13, face="bold.italic")
        )
        pe <- pe + ggtitle(bquote(atop(.(gg_title), atop(.(subtitle))))) + xlab("Jours") + ylab(fonc)
        
        p <- grid.arrange(pi1,pi2,ps,pe,ncol=2,nrow=2)
        print(p)
      }
      #
    }
    #
  }
  #
}

## Evolution d'un carr\'e de voxels : suivi sur 22 jours

carre_dyn <- function(rat,fonc,liste_coupes,sommet,mesure,zone,t,opt){
  #jour <- "08"
  #num_jour <- 8
  suff <- ""
  for (cp in liste_coupes){
    suff <- paste(suff,"-",cp,sep='')
  }
  
  zone <- paste(zone,"-","1",sep='')
  
  # Chargement du cerveau dont les données seront utilisées
  cer00 <- {
    nom_cerveau <- sprintf("R%s/automate_3demi/cerveau%s_multi_J%s%s_Slices%s.dat",rat,rat,jour,opt,suff)
    cerveau <- read.table(nom_cerveau,header=T)
    l <- length(cerveau$x)
    cerveau <- as.data.frame(cbind(cerveau,rep(num_jour,l)))
    colnames(cerveau) <- c(c('x','y','z','Slice'),vect_fonc,c('Etat','Jour'))
  }
  # Constitution du carré de pixels : état zone J00. Pas d'ADC.
  cpix <- {
    pixels <- cerveau[,-cerveau$ADC]
    pixels$Etat <- rep(zone,l)
    liste_l <- rep(FALSE,l)
    xc <- sommet[1]
    yc <- sommet[2]
    for (i in c(1:mesure)){
      for (j in c(1:mesure)){
        ix <- xc + i-1
        iy <- yc + j-1
        liste_l <- ifelse(ix==pixels$x&iy==pixels$y,TRUE,liste_l)
      }
    }
    pixels <- pixels[liste_tr,]
  }
  # Ajout du temps caractéristique de différentiation des zones
  lp <- length(pixels$Etat)
  T_diff <- rnorm(lp,5,2)
  pixels <- as.data.frame(cbind(pixels,T_diff))
  colnames(pixels) <- c(colnames(pixels),'Tchar')
  
  # Constitution de la base de données : constantes du modèle
  ## Distributions asymptotiques
  distr <- {
    cbf_asy1 <- rnorm(l,50,40)
    cmro2_asy1 <- rnorm(l,5,1.5)
    so2_asy1 <- rnorm(l,75,25)
    vsi_asy1 <- rnorm(l,4,2)
    cbf_trans3 <- rnorm(l,120,40)
    vsi_asy2 <- rnorm(l,5,2)
    so2_asy2 <- rnorm(l,75,15)
    so2_sain <- rnorm(l,75,10)
    
    cbf_sain <- rnorm(l,130,30)
    
    adc_asy1 <- rnorm(l,3000,1000)
    adc_asy2 <- rnorm(l,2500,700)
  }
  const <- as.data.frame(cbind(cbf_asy1,cmro2_asy1,so2_asy1,cbf_trans2,vsi_asy2,so2_asy2,so2_sain,cbf_sain))
  colnames(const)=c("CBF_asy1","CMRO2_asy1","SO2_asy1","CBF_trans2","VSI_asy2","SO2_asy2","SO2_sain","CBF_sain")
  
  ## Coefficients multiplicateurs
  qqdis <- {
    q_cbf_asy1 <- q_dis(50,40,7)#rnorm(l,50,50)
    q_cmro2_asy1 <- q_dis(8,1.5,7)#rnorm(l,5,5.5)
    q_so2_asy1 <- q_dis(35,25,7)#rnorm(l,75,35)
    q_vsi_asy1 <- q_dis(3,2,7)#rnorm(l,4,2.8)
    q_cbf_trans2 <- q_dis(50,40,7)#rnorm(l,120,50)
    q_vsi_asy2 <- q_dis(3,2,7)#rnorm(l,5,4.3)
    q_so2_asy2 <- q_dis(20,15,7)#rnorm(l,75,30)
    q_so2_sain <- q_dis(30,10,7)#rnorm(l,75,25)
    
    q_cbf_sain <- q_dis(50,30,5)#rnorm(l,130,30)
  }
  
  # Boucle de simulations
  for (num_jour in c(1:t)){
    # Increment
    c_inc <- pixels[pixels$Jour==num_jour-1,]
    # Conversion : jours
    c_inc$Jour <- rep(num_jour,lp)
    # Accroissement des modalités, en un jour
    c_t <- data.frame(matrix(ncol = 6, nrow = l))
    colnames(c_t) <- c('BVf','CBF','CMRO2','SO2map','VSI')
    
    # Calcul, vectorisé et conditionnel, des composantes des vecteurs de cerveau_temp.
    # On commence par l'état d'équilibre sain, perturbé non lésé, lésé 1 début, lésé 1 fin etc.
    # Il existe des niveaux de priorités entre les modalités
    ## SO2
    # Lésion 1
    c_t$SO2map <- ifelse(c_inc$Etat=='les_1_deb',(1-q_so2_asy1)*(so2_asy1-c_inc$SO2map),c_t$SO2map)
    c_t$SO2map <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$SO2map,c_t$SO2map)
    # Lésion 2
    c_t$SO2map <- ifelse(c_inc$Etat=='les_2_deb',(1-q_so2_sain)*(so2_sain-c_inc$SO2map),c_t$SO2map)
    c_t$SO2map <- ifelse(c_inc$Etat=='les_2_ifin',(1-q_so2_asy2)*(so2_asy2-c_inc$SO2map),c_t$SO2map)
    
    ## CBF
    # Lésion 1
    c_t$CBF <- ifelse(c_inc$Etat=='les_1_deb',(1-q_cbf_asy1)*(cbf_asy1-c_inc$CBF),c_t$CBF)
    c_t$CBF <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$CBF,c_t$CBF)
    # Lésion 2
    c_t$CBF <- ifelse(c_inc$Etat=='les_2_deb',(1-q_cbf_trans2)*(cbf_trans2-c_inc$CBF),c_t$CBF)
    c_t$CBF <- ifelse(c_inc$Etat=='les_2_fin',-10*c_t$SO2map,c_t$CBF)
    # Tissu sain (contra..)
    c_t$CBF <- ifelse(c_inc$Etat=='per',(1-q_cbf_sain)*(cbf_sain-c_inc$CBF),c_t$CBF)
    
    ## CMRO2
    # Lésion 1
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_1_deb',(1-q_cmro2_asy1)*(cmro2_asy1-c_inc$CMRO2),c_t$CMRO2)
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$CMRO2,c_t$CMRO2)
    # Lésion 2
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_2_deb',c_t$SO2map,c_t$CMRO2)
    c_t$CMRO2 <- ifelse(c_inc$Etat=='les_2_fin',-1*c_t$SO2map,c_t$CMRO2)
    
    ## VSI
    # Lésion 1
    c_t$VSI <- ifelse(c_inc$Etat=='les_1_deb',(1-q_vsi_asy1)*(vsi_asy1-c_inc$VSI),c_t$VSI)
    c_t$VSI <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$VSI,c_t$VSI)
    # Lésion 2
    c_t$VSI <- ifelse(c_inc$Etat=='les_2_deb',(1-q_vsi_asy2)*(vsi_asy2-c_inc$VSI),c_t$VSI)
    c_t$VSI <- ifelse(c_inc$Etat=='les_2_fin',bruit_mod$VSI,c_t$VSI)
    
    ## BVf
    # Lésion 1
    c_t$BVf <- ifelse(c_inc$Etat=='les_1_deb',0.5*c_t$VSI,c_t$BVf)
    c_t$BVf <- ifelse(c_inc$Etat=='les_1_fin',bruit_mod$BVf,c_t$BVf)
    # Lésion 2
    c_t$BVf <- ifelse(c_inc$Etat=='les_2_deb',2*c_t$VSI,c_t$BVf)
    c_t$BVf <- ifelse(c_inc$Etat=='les_2_fin',bruit_mod$BVf,c_t$BVf)
    
    # Calcul de la valeurs courante de cerveau_inc : modalités
    for (fonc in liste_fonc){
      c_inc[[fonc]] <- c_inc[[fonc]]+c_t[[fonc]]
    }
    
    # Les modalités ne prennent pas de valeur négative
    for (fonc in liste_fonc){
      c_inc[[fonc]] <- ifelse(c_inc[[fonc]]<0,0,c_inc[[fonc]])
    }
    
    # Calcul des états de cerveau_inc
    c_inc$Etat <- ifelse((c_inc$Etat=='per'&(abs(c_inc$CBF-cbf_sain)<0.5*abs(bruit_mod$CBF))),'sain',as.character(c_inc$Etat))
    #c_inc$Etat <- ifelse((c_inc$Etat=='les_1_deb'&(abs(c_inc$CBF-cbf_sain)<0.5*abs(bruit_mod$CBF))),'les_1_fin',as.character(c_inc$Etat))
    c_inc$Etat <- ifelse(c_inc$Etat=='les_2_deb'&(abs(c_inc$SO2map-so2_sain)<0.5*abs(bruit_mod$SO2map)),'les_2_fin',as.character(c_inc$Etat))
    
    # Changements d'états pour les pixels de la première lignée
    c_inc$Etat <- ifelse((c_inc$Etat=='les_1_1')&c_inc$Jour > c_inc$Tchar,'les_1_2',c_inc$Etat)
    c_inc$Etat <- ifelse((c_inc$Etat=='les_1_2')&c_inc$CMRO2<10,'les_1_3',c_inc$Etat)
    # Changements d'états pour les pixels de la deuxième lignée
    c_inc$Etat <- ifelse((c_inc$Etat=='les_2_1')&c_inc$Jour > c_inc$Tchar,'les_2_2',c_inc$Etat)
    c_inc$Etat <- ifelse((c_inc$Etat=='les_2_2')&(c_inc$CBF>100|c_inc$CMRO2>10),'les_2_3',c_inc$Etat)
    #c_inc$Etat <- ifelse((c_inc$Etat=='les_2_3')&(abs(c_inc$SO2map-so2_sain)<...),'les_2_4',c_inc$Etat)
    
    
    
    # Jour courant pour cerveau_inc
    c_inc$Jour <- rep(num_jour,l)
    
    # Concaténation
    pixels <- as.data.frame(rbind(pixels,c_inc))
    #colnames(cerveau) <- c(c('x','y','z','Slice'),vect_fonc,c('Etat','Jour'))
    
    # Fin de la boucle
  }
  
  
  
  
  
  #return(pixels)
  nom_carre_dyn <- sprintf("R%s/automate_3demi/carre%s_multi_dyn_Slices%s_fin%i.dat",rat,rat,suff,t)
  write.table(cerveau, nom_table_dyn, row.names=F, quote=F, sep='\t')
}











