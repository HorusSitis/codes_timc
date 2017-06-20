### Paquets et librairies

rm(list=ls())

#install.packages("mclust", dependencies=T)
#install.packages('mixtools', dependencies=T)
#install.packages('ggplot2')
#install.packages("plot3D",dependencies=T)

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

suivi_temp_ADC <- list("11"=c(10,10),"19"=c(0,0),"26"=c(0,0),"30"=c(0,0))

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
#---------------------------- Base de données : SO2map ----------------------------#
#---------------------------- Base de données : T1map ----------------------------#
#---------------------------- Base de données : VSI ----------------------------#

liste_jfr <- list("Anat"=liste_jours_Anat, "ADC"=liste_jours_ADC, "BVf"=liste_jours_BVf,"CBF"=liste_jours_CBF,"CMRO2"='',"T1map"='',"VSI"='') # liste des jours par fonctionnalité et par rat
suivi_temp <- list("Anat"=suivi_temp_Anat,"ADC"='',"BVf"='',"CBF"='',"CMRO2"='',"SO2map"='',"T1map"='',"VSI"='')

#################################### Instructions : traitement systématique des images, anatomique ####################################

## Histogrammes
## Angélique veut un mélange de Gaussiennes
## On sépare la partie ischémiée et ...

# Comme pour les autres modalités :





jours_R11_Anat <- c("00","03","08","15","22")
jours_R19_Anat <- c("00","03","08","15","22")
jours_R26_Anat <- c("00","03","08","15","22")
jours_R30_Anat <- c("00","08","15")

liste_jours_Anat <- list("11"=jours_R11_Anat,"19"=jours_R19_Anat,"26"=jours_R26_Anat,"30"=jours_R30_Anat)

suivi_temp_Anat <- list("11"=c(4,12),"19"=c(6,12),"26"=c(4,11),"30"=c(4,13))

# Etape 1 :

#FONC_3d_rat('Anat',"11")
#FONC_3d_rat('Anat',"19")
#FONC_3d_rat('Anat',"26")
#FONC_3d_rat('Anat',"30")

# Etape 2 :

rg_FONC_3d('Anat',"11",2,4)
rg_FONC_3d('Anat',"19",3,4)
rg_FONC_3d('Anat',"26",3,4)
rg_FONC_3d('Anat',"30",2,5)

# Etape 3 :

ddd <- read.table("11-J00-Anat-bg-all.dat",header=T)

ddf <- cluster_jfr_10(ddd,3,4)
ddg <- as.data.frame(cbind(ddd$x,ddd$y,ddd$z,ddf$data,ddf$classification))




ddh <- ddg[ddd$z>3*d.slice.size & ddd$z<13*d.slice.size,]

ddh <- ddg[ddd$z>=3*d.slice.size-0.001 & ddd$z<=13*d.slice.size+0.001,]
View(ddh)



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
rg_FONC_3d('ADC',"19",3,5)
rg_FONC_3d('ADC',"26",3,5)
rg_FONC_3d('ADC',"30",3,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.

seg_cl_FONC("00",'ADC',"11",1)
seg_cl_FONC("00",'ADC',"19",1)
#26
#30



##day <- "00"
fonc <- "ADC"
rat <- "11"



d <- read.table('11-J00-ADC-bg-all.dat',header=T)
cl_test <- cluster_jfr_f10(d,3,5)
d_seg_isch <- d[cl_test$classification==1,]
d_hem_sain <- d[d$x>60,]

adc_entier <- d[,"ADC"]
adc_isch <- d_seg_isch[,"ADC"]
adc_sain <- d_hem_sain[,"ADC"]

taille <- length(d[,4])
taille_isch <- length(adc_isch)

breaks <- seq(min(d[,4])-0.1*min(d[,4]), max(d[,4])+0.1*max(d[,4]), length.out=100)

par(mfrow = c(2,2))
#plot(cl_test, what="classification", col=color.vector)
plot(d$x,d$y,col=color.vector[cl_test$classification], pch=20, cex=2*(1-cl_test$uncertainty)^4, xlab='x', ylab='y',main=paste("Cerveau ","11",", J","00"))
hist_entier <- hist(adc_entier,breaks=breaks, col='grey50',main="Cerveau entier")

data <- adc_isch*(taille_isch/taille)
breaks <- seq(min(data)-0.1*min(data), max(data)+0.1*max(data), length.out=100)
hist_isch <- hist(data,breaks=breaks, col='red',main="Zone ischémiée")

hist_sain <- hist(adc_sain,breaks=breaks, col='blue',main="Hémisphère sain")

write.table(d_seg_isch, sprintf("%s-J%s-%s-isch.dat",rat,"00",fonc), row.names=F, quote=F, sep='\t')

aa <- read.table(sprintf("%s-J%s-%s-isch.dat",rat,"00",fonc),header=T)

hist_aa <- hist(aa[,'ADC'],breaks=breaks, col='red',main="Zone ischémiée aa")

# Etape 4 :








# ----------------- La fonctionnalité BVf ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .

#FONC_3d_rat('BVf',"11")
#FONC_3d_rat('BVf',"19")
#FONC_3d_rat('BVf',"26")
#FONC_3d_rat('BVf',"30")

# Etape 2 :

rg_FONC_3d('BVf',"11",3,4)
rg_FONC_3d('BVf',"19",3,3)
rg_FONC_3d('BVf',"26",3,5)
rg_FONC_3d('BVf',"30",2,5)

# Etape 3 : on segmente la zone ischémiée aux jours où cela est possible, avec les clusters choisis manuellement cf étape 2.

#seg_cl_FONC('BVf',"11",1) refaire les représentations graphiques
#19
#26
#30




# Etape 4 :







# ----------------- La fonctionnalité CBF ----------------- #

# Etape 1 : .csv complète la base de données. Retire les NaN des fichiers all.dat .

#FONC_3d_rat('CBF',"11")
#FONC_3d_rat('CBF',"19")
#FONC_3d_rat('CBF',"26")
#FONC_3d_rat('CBF',"30")

# Etape 2 :

rg_FONC_3d('CBF',"11",3,4)
rg_FONC_3d('CBF',"19",3,4)
rg_FONC_3d('CBF',"26",3,3)
rg_FONC_3d('CBF',"30",2,5)

# Etape 3 :

d.seg <- seg_cl_FONC("22",'CBF',"26",3)



# Etape 4 :







