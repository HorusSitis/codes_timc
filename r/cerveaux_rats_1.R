### Paquets et librairies

rm(list=ls())

#install.packages("mclust", dependencies=T)
#install.packages('mixtools', dependencies=T)
#install.packages('ggplot2')

# Pour le boot linux de Grenoble
library('mclust')
library('mixtools')
library('cluster')
library('ggplot2')

#################################### En tête : base de données constantes ####################################

color.vector <- c('darkblue','red','green','purple', 'orange')
d.slice.size <- 1/14*10000
source(file="cerveaux_rats_fun.R", encoding ="UTF-8")

noms_rats <- c("11","19","26","30")

#---------------------------- Base de données : ADC ----------------------------#

jours_R11_ADC <- c("00","03","08","15","22")
jours_R19_ADC <- c("00","03","08","15","22")
jours_R26_ADC <- c("00","03","08","15","22")
jours_R30_ADC <- c("00","08","15")

liste_jours_ADC <- list("11"=jours_R11_ADC,"19"=jours_R19_ADC,"26"=jours_R26_ADC,"30"=jours_R30_ADC)

#---------------------------- Base de données : BVf ----------------------------#

#jours_R11_BVf <- 
#liste_jours_BVf <- 

#---------------------------- Base de données : CBF ----------------------------#
#---------------------------- Base de données : CMRO2 ----------------------------#
#---------------------------- Base de données : SO2map ----------------------------#
#---------------------------- Base de données : T1map ----------------------------#
#---------------------------- Base de données : VSI ----------------------------#

liste_jfr <- list("ADC"=liste_jours_ADC, "BVF"='',"CBF"='',"CMRO2"='',"T1map"='',"VSI"='') # liste des jours par fonctionnalité et par rat

#################################### Instructions : traitement systématique des images, anatomique ####################################

# Histogrammes
# Angélique veut un mélange de Gaussiennes
# On sépare la partie ischémiée et ...

# Etape 1 :
# Etape 2 :



#################################### Instructions : traitement systématique des images, fonctionnalités ####################################

# Une fonctionnalité, un rat. On boucle sur les jours d'examens. On se place dans le répertoire correspondant : rat/fonctionnel_gris_fonctionnalité.
# Il faut retirer les valeurs NaN présentes sur plusieurs fichiers .txt.

# 1- Génération de fichiers .dat avec la fonction cerveau_jfr, fichiers .dat par tranche et pour le cerveau entier générées.
# 2- Représentation graphique systématique des résultats d'examen. Utilise la boucle rg_FONC_3d.
# 3- Extraction de données ... INCOMPLET.
# 4- Effectuer le suivi temporel sur les slices dont les images sont disponible sur toute la durée des examens. Intervient nécessairement après appel de la fonctions cluster_jrf_f10.

# ----------------- La fonctionnalité ADC ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .

#FONC_3d_rat('ADC',"11")
#FONC_3d_rat('ADC',"19")
#FONC_3d_rat('ADC',"26")
#FONC_3d_rat('ADC',"30")

# Etape 2 :

rg_FONC_3d('ADC',"11",3,5)
rg_FONC_3d('ADC',"19",3,4)
rg_FONC_3d('ADC',"26",3,5)
rg_FONC_3d('ADC',"30",2,5)

# Etape 3 :





# Etape 4 :







