### Paquets et librairies

rm(list=ls())

#install.packages("mclust", dependencies=T)
#install.packages('mixtools', dependencies=T)
#install.packages('ggplot2')
#install.packages("plot3D",dependencies=T)
#install.packages("scatterplot3d",dependencies=T)

# Pour le boot linux de Grenoble
library('mclust')
library('mixtools')
library('cluster')
library('ggplot2')
library('scatterplot3d')

#################################### En tête : base de données constantes ####################################

color.vector <- c('darkblue','red','green','purple', 'orange')
color_seg <- c('grey50','red','blue','black')
d.slice.size <- 1/14*10000
source(file="cerveaux_rats_fun.R", encoding ="UTF-8")

noms_rats <- c("11","19","26","30")

#---------------------------- Base de données : ADC ----------------------------#

jours_R11_ADC <- c("00","03","08","15","22")
jours_R19_ADC <- c("00","03","08","15","22")
jours_R26_ADC <- c("00","03","08","15","22")
jours_R30_ADC <- c("00","08","15")

liste_jours_ADC <- list("11"=jours_R11_ADC,"19"=jours_R19_ADC,"26"=jours_R26_ADC,"30"=jours_R30_ADC)

#suivi_temp_ADC <- list("11"=c(10,10),"19"=c(9,10),"26"=c(0,0),"30"=c(0,0)) --> mettre dans les répertoires, .csv

#---------------------------- Base de données : BVf ----------------------------#

jours_R11_BVf <- c("00","03","08","15","22") # .csv OK, attention au jour 00 où il faudra peut-être éliminer des slices.
jours_R19_BVf <- c("00","08","15","22")
jours_R26_BVf <- c("00","03","08","15","22")
jours_R30_BVf <- c("00","08","15")

liste_jours_BVf <- list("11"=jours_R11_BVf,"19"=jours_R19_BVf,"26"=jours_R26_BVf,"30"=jours_R30_BVf)

#---------------------------- Base de données : CBF ----------------------------#

jours_R11_CBF <- c("00","03","08","15","22") # attension à la slice 8, J00 qui est dans le fichier .csv
jours_R19_CBF <- c("00","03","08","15","22")
jours_R26_CBF <- c("00","03","08","15","22")
jours_R30_CBF <- c("00","08","15")

liste_jours_CBF <- list("11"=jours_R11_CBF,"19"=jours_R19_CBF,"26"=jours_R26_CBF,"30"=jours_R30_CBF)

#---------------------------- Base de données : CMRO2 ----------------------------#

jours_R11_CMRO2 <- c("00","03","15","22") # attension à la slice 8, J00 qui est dans le fichier .csv
jours_R19_CMRO2 <- c("00","08","15","22")
jours_R26_CMRO2 <- c("00","03","08","15","22")
jours_R30_CMRO2 <- c("00","08","15")

liste_jours_CMRO2 <- list("11"=jours_R11_CMRO2,"19"=jours_R19_CMRO2,"26"=jours_R26_CMRO2,"30"=jours_R30_CMRO2)

#---------------------------- Base de données : SO2map ----------------------------#

jours_R11_SO2map <- c("00","03","15","22") # attension à la slice 8, J00 qui est dans le fichier .csv
jours_R19_SO2map <- c("00","08","15","22")
jours_R26_SO2map <- c("00","03","08","15","22")
jours_R30_SO2map <- c("00","08","15")

liste_jours_SO2map <- list("11"=jours_R11_SO2map,"19"=jours_R19_SO2map,"26"=jours_R26_SO2map,"30"=jours_R30_SO2map)

#---------------------------- Base de données : T1map ----------------------------#

jours_R11_T1map <- c("00","03","08","15","22")
jours_R19_T1map <- c("00","03","08","15","22")
jours_R26_T1map <- c("00","03","08","15","22")
jours_R30_T1map <- c("00","08","15")

liste_jours_T1map <- list("11"=jours_R11_T1map,"19"=jours_R19_T1map,"26"=jours_R26_T1map,"30"=jours_R30_T1map)

#---------------------------- Base de données : VSI ----------------------------#

jours_R11_VSI <- c("00","03","08","15","22") # attension à la slice 8, J00 qui est dans le fichier .csv
jours_R19_VSI <- c("00","08","15","22")
jours_R26_VSI <- c("00","03","08","15","22")
jours_R30_VSI <- c("00","08","15")

liste_jours_VSI <- list("11"=jours_R11_VSI,"19"=jours_R19_VSI,"26"=jours_R26_VSI,"30"=jours_R30_VSI)

#############---------------------------- Anatomique ----------------------------#############

jours_R11_Anat <- c("00","03","08","15","22")
jours_R19_Anat <- c("00","03","08","15","22")
jours_R26_Anat <- c("00","03","08","15","22")
jours_R30_Anat <- c("00","08","15")

liste_jours_Anat <- list("11"=jours_R11_Anat,"19"=jours_R19_Anat,"26"=jours_R26_Anat,"30"=jours_R30_Anat)

suivi_temp_Anat <- list("11"=c(4,12),"19"=c(6,12),"26"=c(4,11),"30"=c(4,13))

#---------------------------- Toutes les fonctionnalités et anatomique : accès aux bases de données ----------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')

liste_jfr <- list("Anat"=liste_jours_Anat, "ADC"=liste_jours_ADC, "BVf"=liste_jours_BVf,"CBF"=liste_jours_CBF,"CMRO2"=liste_jours_CMRO2,"SO2map"=liste_jours_SO2map,"T1map"=liste_jours_T1map,"VSI"=liste_jours_VSI) # liste des jours par fonctionnalité et par rat

suivi_temp <- list("Anat"=suivi_temp_Anat,"ADC"='',"BVf"='',"CBF"='',"CMRO2"='',"SO2map"='',"T1map"='',"VSI"='')

# ++++++++++++++++++++++++++++++++++ Approche inverse : fonctionnalités utilisables par rat, pour la segmentation ++++++++++++++++++++++++++++++++++ #

liste_R11_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')
liste_R19_seg_FONC <-list("00"='ADC',"03"='T1map',"08"='ADC',"15"='ADC',"22"='ADC')
liste_R26_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')
liste_R30_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')

liste_sfr <- list("11"=liste_R11_seg_FONC, "19"=liste_R19_seg_FONC, "26"=liste_R26_seg_FONC, "30"=liste_R30_seg_FONC)

#################################### Instructions : traitement systématique des images, anatomique ####################################

## Histogrammes
## Angélique veut un mélange de Gaussiennes
## On sépare la partie ischémiée et ...

# Comme pour les autres modalités :

# Etape 1 :

#FONC_3d_rat('Anat',"11")
#FONC_3d_rat('Anat',"19")
#FONC_3d_rat('Anat',"26")
#FONC_3d_rat('Anat',"30")

# Etape 2 :

rg_FONC_3d(2,'Anat',"11",2,4)
rg_FONC_3d(2,'Anat',"19",3,5)
rg_FONC_3d(2,'Anat',"26",3,5)
rg_FONC_3d(2,'Anat',"30",2,5)

# Etape 3 : pas de sementation possible avec des clusters

####################################### Rat numéro 11 #######################################



for (fonc in liste_fonc){# On parcourt en largeur d'abord l'arborescence de la base de données pour savoir quelles fonctionnalités sont disponibles pour le jour courant.
  liste_jf <- liste_jfr[[fonc]]
  liste_j <- liste_jf[[rat]]
  if (any(liste_j==jour)){
    print("ok")
    liste_F <- cbind(liste_F,list(fonc))# Liste des fonctionnalités disponibles
  }
  else{print("non")}
}


liste_jf <- liste_jfr[["ADC"]]



# Etape 5 : histogrammes ou courbes de niveaux de gris
liste_R11_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')

gr_ngris_seg("00",'ADC',"11")

####################################### Rat numéro 19 #######################################


ll <- gr_ngris_seg("00",'ADC',"19")
ll <- gr_ngris_seg("03","T1map","19")
ll <- gr_ngris_seg("08","ADC","19")
gr_ngris_seg("15","ADC","19")
gr_ngris_seg("22","ADC","19")

for (j in 1:length(jours_R19_ADC)){
  jour <- jours_R19_ADC[j]
  fonc <- liste_R19_seg_FONC[[jour]]
  gr_ngris_seg(jour,fonc,"19")
}


dgris_temp_fonc("19",'cer','ADC')



####################################### Protocole pour réaliser la présentation. IRM fonctionnelle #######################################

# 1- Génération de fichiers .dat avec la fonction cerveau_jfr, fichiers .dat par tranche et pour le cerveau entier générées.
## --> FONC_3d, répertoires fonctionnel_gris/fonc. fichiers all.dat avec z et slices.

# 2- Représentation graphique systématique des résultats d'examen, Utilisation de clusters pour la segmentation.
## --> rg_FONC_3d, figure tridimentionnelle ou projeté sur le plan frontal.
## --> seg_clust_3d.

# 3- Extraction de données. On s'appuie sur l'étude graphique de l'étape précédente.
## --> seg_cl_FONC.

# 4- Suivi des densités de niveaux de gris, par fonctionnalité, en utilisant des segmentations : commune à toutes les fonctionnalités ou respectivement réalisées avec les foctionnalités étudiées.
# Pour chaque jour, on évalue les niveaux de gris sur la zone segmentée pour l'ischémie, l'hémisphère sain au jour 00 et le cerveau entier.
## --> dgris_temp_fonc. Alterner avec les cerveaux - tranches pour les représentations graphiques ?

# 5- On représente l'évolution de l'aire ou volume de la zone ischémiée sur des tranches bien choisies ou le cerveau entier.
# Les nuages de points représentant l'évolution des différentes fonctionnalités figurent sur un même graphique.

## --> suivi_etendue_fonc.

###################################################### Rat numéro 11 ######################################################

num_rat <- "11"

#------------ Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. ------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

#------------ Etape 2 : répertoires des fonctionnalités. ------------#

fr_hemi <- c(4,-170) # inter-hémisphère, supposé plan

cl <- c(3,5)# encadrement du nombre de clusters

rg_FONC_3d(2,'ADC',num_rat,cl)
rg_FONC_3d(3,'ADC',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')#1,"03"=c(1,3),"08"=c(3,4),"15"=c(4,5),"22"=4)
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('ADC',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'BVf',num_rat,cl)
rg_FONC_3d(3,'BVf',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('BVF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CBF',num_rat,cl)
rg_FONC_3d(3,'CBF',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CBF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CMRO2',num_rat,cl)
rg_FONC_3d(3,'CMRO2',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CMRO2',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'SO2map',num_rat,cl)
rg_FONC_3d(3,'SO2map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('SO2map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'T1map',num_rat,cl)
rg_FONC_3d(3,'T1map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('T1map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'VSI',num_rat,cl)
rg_FONC_3d(3,'VSI',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('VSI',num_rat,cl,cl_se,fr_hemi)

#------------ Etape 3 : ------------#

liste_clust <- list('ADC'=list("00"=c(1),"03"=c(),"08"=c(3,4,5),"15"=c(4,5),"22"=c(3)), # clusters utilisables, cf étape 2
                    'BVf'=list(),
                    'CBF'=list(),
                    'CMRO2'=list(),
                    'SO2map'=list(),
                    'T1map'=list(),
                    'VSI'=list()
)

record_seg_cl(num_rat,liste_clust,fr_hemi)

#------------ Etape 4 : ------------#
#------------ Etape 5 : ------------#

###################################################### Rat numéro 19 ######################################################

num_rat <- "19"

#------------ Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. ------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

#------------ Etape 2 : répertoires des fonctionnalités. ------------#

fr_hemi <- c(4,-170) # inter-hémisphère, supposé plan

cl <- c(3,5)# encadrement du nombre de clusters
cl <- c(4,4)
cl <- c(3,3)

rg_FONC_3d(2,'ADC',num_rat,cl)
rg_FONC_3d(3,'ADC',num_rat,cl)
##

# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['ADC']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('ADC',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'BVf',num_rat,cl)
rg_FONC_3d(3,'BVf',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['BVf']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('BVf',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CBF',num_rat,cl)
rg_FONC_3d(3,'CBF',num_rat,cl)
##
cl <- c(3,3)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['CBF']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CBF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CMRO2',num_rat,cl)
rg_FONC_3d(3,'CMRO2',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['CMRO2']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CMRO2',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'SO2map',num_rat,cl)
rg_FONC_3d(3,'SO2map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['SO2map']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('SO2map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'T1map',num_rat,cl)
rg_FONC_3d(3,'T1map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['T1map']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('T1map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'VSI',num_rat,cl)
rg_FONC_3d(3,'VSI',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- liste_clust[['VSI']]
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('VSI',num_rat,cl,cl_se,fr_hemi)

#------------ Etape 3 :

liste_clust <- list('ADC'=list("00"=c(1),"03"=c(1),"08"=c(3,4),"15"=c(4),"22"=c(4)), # clusters utilisables, cf étape 2
                    'BVf'=list("00"=c(1),"08"=c(3,4),"15"=c(),"22"=c()),
                    'CBF'=list("00"=c(1),"03"=c(),"08"=c(4),"15"=c(2),"22"=c(1,3)), # 3 clusters
                    'CMRO2'=list("00"=c(1),"08"=c(3),"15"=c(1),"22"=c(3)),
                    'SO2map'=list("00"=c(0),"08"=c(1),"15"=c(4),"22"=c()),
                    'T1map'=list("00"=c(3),"03"=c(3),"08"=c(3,4,5),"15"=c(4,5),"22"=c(4)),
                    'VSI'=list("00"=c(1),"08"=c(2,3),"15"=c(3),"22"=c(3))
                    )

record_seg_cl(num_rat,liste_clust,fr_hemi)

#------------ Etape 4 : ------------#

liste_suivi_clust <- list('ADC'=list(9),
                          'BVf'=list(8,9,10),
                          'CBF'=list(9,10),
                          'CMRO2'=list(9,10),
                          'SO2map'=list(10),
                          'T1map'=list(10),
                          'VSI'=list()
                          )

#liste_jf <- #--------> adaptée à l'éventuelle absence de segmentation
liste_fonc <- list('T1map')# on peut tester avec une liste des fonctionnalités restreinte

dgris_temp_fonc(num_rat,'cer',"")
dgris_temp_fonc(num_rat,'cer',"T1map")

dgris_temp_fonc(num_rat,liste_suivi_clust[['T1map']],"")
dgris_temp_fonc(num_rat,,"T1map")

#------------ Etape 5 : ------------#

###################################################### Rat numéro 26 ######################################################

num_rat <- "26"

#------------ Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. ------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

#------------ Etape 2 : répertoires des fonctionnalités. ------------#

fr_hemi <- c(4,-170) # inter-hémisphère, supposé plan

cl <- c(3,5)# encadrement du nombre de clusters

rg_FONC_3d(2,'ADC',num_rat,cl)
rg_FONC_3d(3,'ADC',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')#1,"03"=c(1,3),"08"=c(3,4),"15"=c(4,5),"22"=4)
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('ADC',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'BVf',num_rat,cl)
rg_FONC_3d(3,'BVf',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('BVF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CBF',num_rat,cl)
rg_FONC_3d(3,'CBF',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CBF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CMRO2',num_rat,cl)
rg_FONC_3d(3,'CMRO2',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CMRO2',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'SO2map',num_rat,cl)
rg_FONC_3d(3,'SO2map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('SO2map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'T1map',num_rat,cl)
rg_FONC_3d(3,'T1map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('T1map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'VSI',num_rat,cl)
rg_FONC_3d(3,'VSI',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('VSI',num_rat,cl,cl_se,fr_hemi)

#------------ Etape 3 : ------------#

liste_clust <- list('ADC'=list("00"=c(1),"03"=c(),"08"=c(3,4,5),"15"=c(4,5),"22"=c(3)), # clusters utilisables, cf étape 2
                    'BVf'=list(),
                    'CBF'=list(),
                    'CMRO2'=list(),
                    'SO2map'=list(),
                    'T1map'=list(),
                    'VSI'=list()
)

record_seg_cl(num_rat,liste_clust,fr_hemi)

#------------ Etape 4 : ------------#

#------------ Etape 5 : ------------#

###################################################### Rat numéro 30 ######################################################

num_rat <- "30"

#------------ Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. ------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

#------------ Etape 2 : répertoires des fonctionnalités. ------------#

fr_hemi <- c(4,-170) # inter-hémisphère, supposé plan

cl <- c(3,5)# encadrement du nombre de clusters

rg_FONC_3d(2,'ADC',num_rat,cl)
rg_FONC_3d(3,'ADC',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')#1,"03"=c(1,3),"08"=c(3,4),"15"=c(4,5),"22"=4)
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('ADC',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'BVf',num_rat,cl)
rg_FONC_3d(3,'BVf',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('BVF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CBF',num_rat,cl)
rg_FONC_3d(3,'CBF',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CBF',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'CMRO2',num_rat,cl)
rg_FONC_3d(3,'CMRO2',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('CMRO2',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'SO2map',num_rat,cl)
rg_FONC_3d(3,'SO2map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('SO2map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'T1map',num_rat,cl)
rg_FONC_3d(3,'T1map',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('T1map',num_rat,cl,cl_se,fr_hemi)

rg_FONC_3d(2,'VSI',num_rat,cl)
rg_FONC_3d(3,'VSI',num_rat,cl)
##
cl <- c(4,4)#c(3,5)
# ---> bornes pour le nombre de clusters
cl_se <- list("00"='')
# ---> choix possible de CLusters pour la SEgmentation
seg_clust_3d('VSI',num_rat,cl,cl_se,fr_hemi)

#------------ Etape 3 : ------------#

liste_clust <- list('ADC'=list("00"=c(1),"03"=c(),"08"=c(3,4,5),"15"=c(4,5),"22"=c(3)), # clusters utilisables, cf étape 2
                    'BVf'=list(),
                    'CBF'=list(),
                    'CMRO2'=list(),
                    'SO2map'=list(),
                    'T1map'=list(),
                    'VSI'=list()
)

record_seg_cl(num_rat,liste_clust,fr_hemi)

#------------ Etape 4 : ------------#
#------------ Etape 5 : ------------#