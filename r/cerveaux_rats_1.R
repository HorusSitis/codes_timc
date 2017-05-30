### Paquets et librairies

install.packages("mclust")
library('mclust')
install.packages('mixtools')
library('mixtools')
install.packages('ggplot2')

# Pour le boot linux de Grenoble
library('cluster')
library('ggplot2')
library('mixtools')

## Fonctions à utiliser : voir _fun.R
#source(file="monFichier.R", encoding ="UTF-8")



### Instructions

color.vector <- c('darkblue','red','green','purple')

# Pour des répertoires etc

slash="/"
slash="//" #selon les machines, on exprime différemment le slash dans les chemins
slash="\"
#slash="\\"

## Tests : avec des fichiers individuels



## Boucles : appels des fonctions statistiques sur des séries spatiotemporelles.
# Il faut charger le répertoire courant adéquat

J<-c("J00","J03","J08","J15","J22") # rats 11, 19 et 26
J_30<-c("J00","J08","J15") # rat 30, pour plus tard

# Anatomique :



# Fonctionnel :

var<-'ADC' ######

# Spatial : S00, S03 etc


# Temporel : JS[1] et JS[2]



var<-'BVf' ######


