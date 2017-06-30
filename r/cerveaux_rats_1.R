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
color_seg <- c('grey50','red','blue')
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





#################################### Instructions : traitement systématique des images, fonctionnalités ####################################

# Un rat, une fonctionnalité. On boucle sur les jours d'examens. On se place dans le répertoire correspondant : rat/fonctionnel_gris_fonctionnalité.
# Il faut retirer les valeurs NaN présentes sur plusieurs fichiers .txt.

# 1- Génération de fichiers .dat avec la fonction cerveau_jfr, fichiers .dat par tranche et pour le cerveau entier générées.
# 2- Représentation graphique systématique des résultats d'examen. Utilise la boucle rg_FONC_3d.
# 3- Extraction de données. Non disponible pour toutes les fonctionnalités ou rats.
# 4- Effectuer le suivi temporel sur les slices dont les images sont disponibles sur toute la durée des examens. Intervient nécessairement après appel de la fonctions cluster_jrf_f10.
# 5- Comparer la distribution des niveaux de gris entre zone ischémiée, hémisphère sain et cerveau entier pour toutes les fonctionnalités, après segmentation de la zone ischémiée à l'aide d'une fonctionnalité.

####################################### Rat numéro 11 #######################################

# ----------------- La fonctionnalité ADC ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('ADC',"11")
# Etape 2 :
rg_FONC_3d(2,'ADC',"11",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
liste_R11_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.

# ------ Jour 00 ------ #
d_seg_00 <- seg_cl_FONC("00",'ADC',"11",1,c(4,-170))
write.table(d_seg_00, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"11","00"), row.names=F, quote=F, sep='\t')
# ------ Jour 15 ------ #
d_seg_15 <- seg_cl_FONC("15",'ADC',"11",3,c(4,-170))
write.table(d_seg_15, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"11","15"), row.names=F, quote=F, sep='\t')
# ------ Jour 22 ------ #
d_seg_22 <- seg_cl_FONC("22",'ADC',"11",4,c(4,-170))
write.table(d_seg_22, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"11","22"), row.names=F, quote=F, sep='\t')




# Etape 4 : suivi temporel.

# Slice 10, tous les jours.
# On utilise la clusterisation aux jours 00, 15 ou 22.

# Slices 9 et 10 : jours 00 et 03.
# On utilise la segmentation au jour 00.

# Slices 10 à 12 : jours 08 à 22.
# On utilise la segmentation aux jours 15 ou 22.  








# ----------------- La fonctionnalité BVf ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('BVf',"11")
# Etape 2 :
rg_FONC_3d(2,'BVf',"11",3,3)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.

#seg_cl_FONC("00",'BVf',"11",1) pas terrible


# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CBF ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CBF',"11")
# Etape 2 :
rg_FONC_3d(2,'CBF',"11",3,4)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.


# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CMRO2 ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat . Attention : pas de jour 08.
#FONC_3d_rat('CMRO2',"11")
# Etape 2 :
rg_FONC_3d(2,'CMRO2',"11",3,3)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité SO2map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat . Attention : pas de jour 08.
#FONC_3d_rat('SO2map',"11")
# Etape 2 :
rg_FONC_3d(2,'SO2map',"11",3,3)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité T1map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('T1map',"11")
# Etape 2 :
rg_FONC_3d(2,'T1map',"11",3,3)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité VSI ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('VSI',"11")
# Etape 2 :
rg_FONC_3d(2,'VSI',"11",3,3)

# Etape 4 : suivi temporel.

####################################### Rat numéro 19 #######################################

# ----------------- La fonctionnalité ADC ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('ADC',"19")
# Etape 2 :
rg_FONC_3d(3,'ADC',"19",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.

cl <- c(4,4)#list("00"=c(3,5),"03"=c(3,5),"08"=c(3,5),"15"=c(3,5),"22"=c(3,5))
cl_se <- list("00"=1,"03"=c(1,3),"08"=c(3,4),"15"=c(4,5),"22"=4)
seg_clust_3d('ADC',"19",cl,cl_se,c(4,-170))

jours_ADC_seg <- c("00","08","15","22")


#liste_R19_seg_FONC <-list("00"='ADC',"03"='T1map',"08"='ADC',"15"='ADC',"22"='ADC')

# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.
d_seg_00 <- seg_cl_FONC("00",'ADC',"19",1,c(4,-170))
d_seg_08 <- seg_cl_FONC("08",'ADC',"19",c(3,4,5),c(4,-170))
d_seg_15 <- seg_cl_FONC("15",'ADC',"19",c(4,5),c(4,-170))
d_seg_22 <- seg_cl_FONC("22",'ADC',"19",3,c(4,-170))
# Retour au répertoire fonctionnel_gris
write.table(d_seg_00, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"19","00"), row.names=F, quote=F, sep='\t')
write.table(d_seg_08, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"19","08"), row.names=F, quote=F, sep='\t')
write.table(d_seg_15, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"19","15"), row.names=F, quote=F, sep='\t')
write.table(d_seg_22, sprintf("isch3d-%s-%s-J%s.dat",'ADC',"19","22"), row.names=F, quote=F, sep='\t')

# Etape 4 : suivi temporel.

# Slices 9 et 10, tous les jours.
# On utilise la clusterisation aux jours ...
# On peut utiliser la segmentation aux jours 00, 08, 15, et 22.

#liste.fonc.19 <- list("00"='ADC',"03"='T1map',"08"='ADC',"15"='ADC',"22"='ADC')
#suivi_temp_fonc("19","ADC","seg") --> pas d'intérêt.
suivi_temp_fonc("19","ADC")#,"clust")

# ----------------- La fonctionnalité BVf ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('BVf',"19")
# Etape 2 :
rg_FONC_3d(2,'BVf',"19",3,5)
rg_FONC_3d(3,'BVf',"19",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.




# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CBF ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CBF',"19")
# Etape 2 :
rg_FONC_3d(3,'CBF',"19",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CMRO2 ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CMRO2',"19")
# Etape 2 :
rg_FONC_3d(2,'CMRO2',"19",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité SO2map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('SO2map',"19")
# Etape 2 :
rg_FONC_3d(2,'SO2map',"19",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité T1map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('T1map',"19")
# Etape 2 :
rg_FONC_3d(2,'T1map',"19",3,5)
rg_FONC_3d(3,'T1map',"19",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.

d_seg_03 <- seg_cl_FONC("03",'T1map',"19",3,c(4,-170))
write.table(d_seg_03, sprintf("isch3d-%s-%s-J%s.dat",'T1map',"19","03"), row.names=F, quote=F, sep='\t')





# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité VSI ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('VSI',"19")
# Etape 2 :
rg_FONC_3d(2,'VSI',"19",3,5)

# Etape 4 : suivi temporel.

####################################### Rat numéro 26 #######################################

# ----------------- La fonctionnalité ADC ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('ADC',"26")
# Etape 2 :
rg_FONC_3d(2,'ADC',"26",3,5)


# Etape 4 : suivi temporel.









# ----------------- La fonctionnalité BVf ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('BVf',"26")
# Etape 2 :
rg_FONC_3d(2,'BVf',"26",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CBF ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CBF',"26")
# Etape 2 :
rg_FONC_3d(2,'CBF',"26",3,3)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.

# ------ Jour 00 ------ #


# ------ Jour 22 ------ #


# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CMRO2 ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CMRO2',"26")
# Etape 2 :
rg_FONC_3d(2,'CMRO2',"26",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité SO2map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('SO2map',"26")
# Etape 2 :
rg_FONC_3d(2,'SO2map',"26",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité T1map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('T1map',"26")
# Etape 2 :
rg_FONC_3d(2,'T1map',"26",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité VSI ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('VSI',"26")
# Etape 2 :
rg_FONC_3d(2,'VSI',"26",3,5)

# Etape 4 : suivi temporel.

####################################### Rat numéro 30 : pas d'étape 3 #######################################

# ----------------- La fonctionnalité ADC ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('ADC',"30")
# Etape 2 :
rg_FONC_3d(2,'ADC',"30",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.
# On enregistre les cerveaux en étiquetant ses pixels à 0, ceux de la zone ischémiée à 1 et enfin ceux de l'hémisphère sain à 2.






# Etape 4 : suivi temporel.

# Slices 11 et 12 : tous les jours. 








# ----------------- La fonctionnalité BVf ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('BVf',"30")
# Etape 2 :
rg_FONC_3d(2,'BVf',"30",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CBF ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CBF',"30")
# Etape 2 :
rg_FONC_3d(2,'CBF',"30",3,3)


# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité CMRO2 ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('CMRO2',"30")
# Etape 2 :
rg_FONC_3d(2,'CMRO2',"30",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité SO2map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('SO2map',"30")
# Etape 2 :
rg_FONC_3d(2,'SO2map',"30",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité T1map ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('T1map',"30")
# Etape 2 :
rg_FONC_3d('T1map',"30",3,5)

# Etape 4 : suivi temporel.

# ----------------- La fonctionnalité VSI ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .
#FONC_3d_rat('VSI',"30")
# Etape 2 :
rg_FONC_3d(2,'VSI',"30",3,5)

# Etape 4 : suivi temporel.

#################################### Traitement systématique des images, fonctionnalités, suite ####################################

# Etape 7/2 : on représente graphiquement les cerveaux segmentés, afin d'évaluer la pertinence de leur utilisation en statistiques -étape 5.

# Etape 5 : densités de niveaux de gris pour la fonctionnalité courante.
# On utilise des segmentations créées avec un fonctionnalité fixée : ADC ou CBF par exemple.
# Jusqu'à nouvel ordre : rats 11, 19 et 26.

# read.table... rat et jour. Trois dimensions.
# boucle sur la liste des fonctionnalités.

####################################### Rat numéro 11 #######################################

jour <- "03"
rat <- "19"

liste_F <- list()

#liste_F 


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

# Etape 7/2 : segmentation automatique, méthode de Nicolas.

# Etape 5 : histogrammes ou courbes de niveaux de gris, cerveau entier à une date donnée puis suivi temporel.
#liste_R19_seg_FONC <-list("00"='ADC',"03"='T1map',"08"='ADC',"15"='ADC',"22"='ADC')

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

####################################### Rat numéro 26 #######################################

# Etape 5 : histogrammes ou courbes de niveaux de gris, cerveau entier à une date donnée puis suivi temporel.
#liste_R26_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')

####################################### Rat numéro 30 #######################################

#liste_R30_seg_FONC <- list("00"='',"03"='',"08"='',"15"='',"22"='')






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

####################################### Rat numéro 11 #######################################

# Etape 1 : répertoires des fonctionnalités.

FONC_3d_rat('ADC',"11")
FONC_3d_rat('BVf',"11")
FONC_3d_rat('CBF',"11")
FONC_3d_rat('CMRO2',"11")
FONC_3d_rat('SO2map',"11")
FONC_3d_rat('T1map',"11")
FONC_3d_rat('VSI',"11")

# Etape 2 : répertoires des fonctionnalités.

rg_FONC_3d(2,'ADC',"11",3,5)
rg_FONC_3d(3,'ADC',"11",3,5)

cl <- c()
cl_se <- list()
seg_clust_3d('ADC',"11",c(3,5),cl,cl_se,c(4,-170))

rg_FONC_3d(2,'BVf',"11",3,5)
rg_FONC_3d(3,'BVf',"11",3,5)

rg_FONC_3d(2,'CBF',"11",3,5)
rg_FONC_3d(3,'CBF',"11",3,5)

rg_FONC_3d(2,'CMRO2',"11",3,5)
rg_FONC_3d(3,'CMRO2',"11",3,5)

rg_FONC_3d(2,'SO2map',"11",3,5)
rg_FONC_3d(3,'So2map',"11",3,5)

rg_FONC_3d(2,'T1map',"11",3,5)
rg_FONC_3d(3,'T1map',"11",3,5)

rg_FONC_3d(2,'VSI',"11",3,5)
rg_FONC_3d(3,'VSI',"11",3,5)

# Etape 3 :

cl_se <- list("00"=c(1),"03"=c(3,5))# etc #c(1) 
seg_cl_FONC("00",'ADC',"19",cl,c(4,-170))

liste_fr <- liste_jfr[[19]]
for (fonc in liste_fonc){
  liste_jour_fonc <- liste_fr[[fonc]]
  for (jour in liste_jours_fonc){
    seg_cl_FONC("00",'ADC',"19",cl$jour,c(4,-170))
    
  }
  
  
  
}

# Etape 4 :
# Etape 5 :

####################################### Rat numéro 19 #######################################

# Etape 1 : répertoires des fonctionnalités.

FONC_3d_rat('ADC',"19")
FONC_3d_rat('BVf',"19")
FONC_3d_rat('CBF',"19")
FONC_3d_rat('CMRO2',"19")
FONC_3d_rat('SO2map',"19")
FONC_3d_rat('T1map',"19")
FONC_3d_rat('VSI',"19")

# Etape 2 :répertoires des fonctionnalités.

rg_FONC_3d(2,'ADC',"19",3,5)
rg_FONC_3d(3,'ADC',"19",3,5)

cl <- c(4,4)#list("00"=c(3,5),"03"=c(3,5),"08"=c(3,5),"15"=c(3,5),"22"=c(3,5))
cl_se <- list("00"=1,"03"=c(1,3),"08"=c(3,4),"15"=c(4,5),"22"=4)
seg_clust_3d('ADC',"19",c(3,5),cl,cl_se,c(4,-170))

rg_FONC_3d(2,'BVf',"19",3,5)
rg_FONC_3d(3,'BVf',"19",3,5)

rg_FONC_3d(2,'CBF',"19",3,5)
rg_FONC_3d(3,'CBF',"19",3,5)

rg_FONC_3d(2,'CMRO2',"19",3,5)
rg_FONC_3d(3,'CMRO2',"19",3,5)

rg_FONC_3d(2,'SO2map',"19",3,5)
rg_FONC_3d(3,'So2map',"19",3,5)

rg_FONC_3d(2,'T1map',"19",3,5)
rg_FONC_3d(3,'T1map',"19",3,5)

rg_FONC_3d(2,'VSI',"19",3,5)
rg_FONC_3d(3,'VSI',"19",3,5)

# Etape 3 :
# Etape 4 :
# Etape 5 :

####################################### Rat numéro 26 #######################################

# Etape 1 : répertoires des fonctionnalités.

FONC_3d_rat('ADC',"26")
FONC_3d_rat('BVf',"26")
FONC_3d_rat('CBF',"26")
FONC_3d_rat('CMRO2',"26")
FONC_3d_rat('SO2map',"26")
FONC_3d_rat('T1map',"26")
FONC_3d_rat('VSI',"26")

# Etape 2 :répertoires des fonctionnalités.

rg_FONC_3d(2,'ADC',"26",3,5)
rg_FONC_3d(3,'ADC',"26",3,5)

cl <- c()
cl_se <- list()
seg_clust_3d('ADC',"26",c(3,5),cl,cl_se,c(4,-170))

rg_FONC_3d(2,'BVf',"26",3,5)
rg_FONC_3d(3,'BVf',"26",3,5)

rg_FONC_3d(2,'CBF',"26",3,5)
rg_FONC_3d(3,'CBF',"26",3,5)

rg_FONC_3d(2,'CMRO2',"26",3,5)
rg_FONC_3d(3,'CMRO2',"26",3,5)

rg_FONC_3d(2,'SO2map',"26",3,5)
rg_FONC_3d(3,'So2map',"26",3,5)

rg_FONC_3d(2,'T1map',"26",3,5)
rg_FONC_3d(3,'T1map',"26",3,5)

rg_FONC_3d(2,'VSI',"26",3,5)
rg_FONC_3d(3,'VSI',"26",3,5)

# Etape 3 :
# Etape 4 :
# Etape 5 :

####################################### Rat numéro 30 #######################################

# Etape 1 : répertoires des fonctionnalités.

FONC_3d_rat('ADC',"30")
FONC_3d_rat('BVf',"30")
FONC_3d_rat('CBF',"30")
FONC_3d_rat('CMRO2',"30")
FONC_3d_rat('SO2map',"30")
FONC_3d_rat('T1map',"30")
FONC_3d_rat('VSI',"30")

# Etape 2 :répertoires des fonctionnalités.

rg_FONC_3d(2,'ADC',"30",3,5)
rg_FONC_3d(3,'ADC',"30",3,5)

cl <- c()
cl_se <- list()
seg_clust_3d('ADC',"30",c(3,5),cl,cl_se,c(4,-170))

rg_FONC_3d(2,'BVf',"30",3,5)
rg_FONC_3d(3,'BVf',"30",3,5)

rg_FONC_3d(2,'CBF',"30",3,5)
rg_FONC_3d(3,'CBF',"30",3,5)

rg_FONC_3d(2,'CMRO2',"30",3,5)
rg_FONC_3d(3,'CMRO2',"30",3,5)

rg_FONC_3d(2,'SO2map',"30",3,5)
rg_FONC_3d(3,'So2map',"30",3,5)

rg_FONC_3d(2,'T1map',"30",3,5)
rg_FONC_3d(3,'T1map',"30",3,5)

rg_FONC_3d(2,'VSI',"30",3,5)
rg_FONC_3d(3,'VSI',"30",3,5)

# Etape 3 :
# Etape 4 :
# Etape 5 :


