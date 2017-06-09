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

## Fonctions à utiliser : voir _fun.R
#source(file="monFichier.R", encoding ="UTF-8")



### Instructions


##########################################################################################

color.vector <- c('darkblue','red','green','purple', 'orange')
d.slice.size <- 1/14*10000
source(file="cerveaux_rats_fun.R", encoding ="UTF-8")

num_rat = '19'
fonc='ADC'

#---------------------------- base de données : rat 19 ; ADC ----------------------------#

# liste_R19_ADC <- list("00"=c(8,9,10),"03"=c(7,8,9,10,11),"08"=c(7,8,9,10,11),"15"=c(8,9,10,11),"22"=c(9,10))
jours_R19_ADC <- c("00","03","08","15","22")

#------------- Autres rats -------------#
#source(file="cerveaux_rats_fun.R", encoding ="UTF-8")

jours_R11_ADC <- c("00","03","08","15","22")

jours_R26_ADC <- c("00","03","08","15","22")
jours_R30_ADC <- c("00","08","15")

liste_jours_ADC <- list("11"=jours_R11_ADC,"19"=jours_R19_ADC,"26"=jours_R26_ADC,"30"=jours_R30_ADC)
noms_rats <- c("11","19","26","30")

# Génération des fichiers .dat, pour une indexation tridimensionnelle des valeurs de la fonctionnalité -ici d'ADC.
# Utilise cerveau_jfr, appelée pour la fonctionnalité ADC.
# Attention au répertoire de travail.

#adc_3d_rat("11")
#adc_3d_rat("19")

adc_3d_rat("30")

# Suivi temporel : études tridimensionnelles par jour.

#suivi_adc_3d("11",3,5) #attention :
#suivi_adc_3d("19",2,5) #n_clusters ?

suivi_adc_3d("30",2,4)

# Suivi temporel par tranche : la clusterisation spatiale est prioritaire.
# Il s'agit d'extraire correctement les résultats de clustersisation, préalablement enregistrés dans un fichier.

