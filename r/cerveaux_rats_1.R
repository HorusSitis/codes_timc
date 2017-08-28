### Paquets et librairies

rm(list=ls())

#install.packages("mclust", dependencies=T)
#install.packages('mixtools', dependencies=T)
#install.packages('ggplot2')
#install.packages("plot3D",dependencies=T)
#install.packages("scatterplot3d",dependencies=T)
#install.packages("gridExtra",dependencies=T)

# ---> install.packages("cowplot",dependencies=T) ---> non disponible pour R 3.2.3.

# Pour le boot linux de Grenoble
library('mclust')
library('mixtools')
library('cluster')
library('ggplot2')
library('scatterplot3d')
library('gridExtra')

# Non disponible : library('cowplot')

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

#FONC_3d_rat('Anat',"11",'avec')
#FONC_3d_rat('Anat',"19",'avec')
#FONC_3d_rat('Anat',"26",'avec')
#FONC_3d_rat('Anat',"30",'avec')

# Etape 2 :

liste_fonc <- list('Anat')
liste_min_fonc <- list('Anat'=10)

cl <- c(2,5)
cl <- c(4,4)

rg_FONC_3d(3,'Anat',"11",cl,'')
rg_FONC_3d(3,'Anat',"19",cl,'')
rg_FONC_3d(3,'Anat',"26",cl,'')
rg_FONC_3d(3,'Anat',"30",cl,'')

for (tr in c(4:12)){
  comp_2vs3d_clust('Anat',"11",cl,tr)
}

for (tr in c(6:12)){
  comp_2vs3d_clust('Anat',"19",cl,tr)
}

#-------------------------------- Etape 3 : segmentation ADC, par rat ; liste_suivi_slice. --------------------------------#

fr_hemi <- c(5,-240) # tranches Anat et autres modalités recalées

# -- Répertoires : ceus des rats
num_rat <- "11"
vol_lesADC00(num_rat,c(10))
num_rat <- "19"
vol_lesADC00(num_rat,c(10))
num_rat <- "26"
vol_lesADC00(num_rat,c(10))
num_rat <- "30"
vol_lesADC00("11",c(10))
# -- Déjà fait au 21/07/0217

# -- listes de suivi temporel ADC
l_s_s_11 <- list('ADC'=list(10),
                 'BVf'=list(10),
                 'CBF'=list(10),
                 'CMRO2'=list(10),
                 'SO2map'=list(10),
                 'T1map'=list(10),
                 'VSI'=list(10)
)
l_s_s_19 <- list('ADC'=list(9,10),
                 'BVf'=list(9,10),
                 'CBF'=list(9),
                 'CMRO2'=list(9,10),
                 'SO2map'=list(9,10),
                 'T1map'=list(9),
                 'VSI'=list(9,10)
)
l_s_s_26 <- list('ADC'=c(6,7,8),
                          'BVf'=c(6,7,8),
                          'CBF'=c(6,7,8),
                          'CMRO2'=c(6,7,8),
                          'SO2map'=c(6,7,8),
                          'T1map'=c(6,7,8),
                          'VSI'=c(6,7,8)
)
l_s_s_30 <- list('ADC'=list(11,12,13),
                 'BVf'=list(11,12,13),
                 'CBF'=list(11,12,13),
                 'CMRO2'=list(11,12,13),
                 'SO2map'=list(11,12,13),
                 'T1map'=list(11,12,13),
                 'VSI'=list(11,12,13)
)

liste_ss_rat <- list("11"=l_s_s_11,"19"=l_s_s_19,"26"=l_s_s_26,"30"=l_s_s_30)

#-------------------------------- Etape 4 : statistiques sur les niveaux de gris. Répertoire du rat concerné --------------------------------#

num_rat <- "11"
num_rat <- "19"
num_rat <- "26"
num_rat <- "30"

liste_fonc <- list('ADC')
liste_suivi_slice <- liste_ss_rat[[num_rat]]

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'pdf')

dgris_temp_fonc(num_rat,fr_hemi,liste_suivi_slice,'dark','') # pas d'option prévue pour une segmentation CBF ou Anat
dgris_temp_fonc(num_rat,fr_hemi,liste_suivi_slice,'dark','pdf')

#-------------------------------- Etape 5 : --------------------------------#

####################################### Protocole pour réaliser la présentation. IRM fonctionnelle #######################################

# 1- Génération de fichiers .dat avec la fonction cerveau_jfr, fichiers .dat par tranche et pour le cerveau entier générées.
## --> FONC_3d, répertoires fonctionnel_gris/fonc. fichiers all.dat avec z et slices.

# 2- Représentation graphique systématique des résultats d'examen, Utilisation de clusters pour la segmentation.
## --> rg_FONC_3d, figure tridimentionnelle ou projeté sur le plan frontal.
## --> seg_clust_3d.
## --> comp_2vs3d_clust pour comparer, sur une tranche d'intérêt, la classification induite depuis la 3d d'une classification issue de la clusterisation 2d.

# 3- Extraction de données. On s'appuie sur l'étude graphique de l'étape précédente..
## --> seg_cl_FONC.

# 4- Suivi des densités de niveaux de gris, par fonctionnalité, en utilisant des segmentations : commune à toutes les fonctionnalités ou respectivement réalisées avec les foctionnalités étudiées.
# Pour chaque jour, on évalue les niveaux de gris sur la zone segmentée pour l'ischémie, l'hémisphère sain au jour 00 et le cerveau entier.
## --> dgris_temp_fonc. Alterner avec les cerveaux - tranches pour les représentations graphiques ?

# 5- On représente l'évolution de l'aire ou volume de la zone ischémiée sur des tranches bien choisies ou le cerveau entier.
# Les nuages de points représentant l'évolution des différentes fonctionnalités figurent sur un même graphique.

## --> suivi_etendue_fonc.

###################################################### Rat numéro 11 ######################################################

num_rat <- "11"

#-------------------------------- Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. --------------------------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

for (fonc in liste_fonc){
  dker_3d_rat(fonc,num_rat)
}

vol_lesADC00(num_rat,c(10))
#vol_lesADC00(num_rat,c(9:10))

# Un seul jour pour la segmentation CBF
vol_lesCBF00(num_rat,c(10))
#vol_lesADC00(num_rat,c(9:10))

#-------------------------------- Etape 2 : répertoires des fonctionnalités. Attention au BVf, J00, slice8. --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI') # ordre de priorité pour les fonctionnalités.
liste_min_fonc <- list('ADC'=10,'BVf'=1,'CBF'=10,'CMRO2'=10,'SO2map'=10,'T1map'=10,'VSI'=10)

cl <- c(4,4)

for (fonc in liste_fonc){
  #rg_FONC_3d(2,fonc,num_rat,cl)
  rg_FONC_3d(3,fonc,num_rat,cl,'')
}

# Remplir liste_clust et liste_suivi_slice en observant les représentations graphiques obtenues

# Vérification : aprè savoir rempli liste_clust

fr_hemi <- c(5,-240) # inter-hémisphère, supposé plan

fonc <- 'ADC'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'BVf'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CBF'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CMRO2'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'SO2map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'T1map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'VSI'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

# trop lente
#for (fonc in liste_fonc){
#  cl_se <- liste_clust[[fonc]]
#  seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)
#}

#-------------------------------- Etape 3 : --------------------------------#

liste_clust <- list('ADC'=list("00"=c(1),"03"=c(),"08"=c(3),"15"=c(4),"22"=c(4)), # clusters utilisables, cf étape 2
                    'BVf'=list("00"=c(1),"03"=c(1),"08"=c(3),"15"=c(3),"22"=c(1)),
                    'CBF'=list("00"=c(1),"03"=c(1),"08"=c(3),"15"=c(2),"22"=c(1)),
                    'CMRO2'=list("00"=c(1),"03"=c(3),"08"=c(3),"15"=c(1),"22"=c(1)),
                    'SO2map'=list("00"=c(0),"08"=c(0,1),"15"=c(1),"22"=c(1)),# ? pour J15
                    'T1map'=list("00"=c(3),"03"=c(3),"08"=c(4),"15"=c(4),"22"=c(4)),
                    'VSI'=list("00"=c(4),"03"=c(4),"08"=c(4),"15"=c(4),"22"=c())
)

record_seg_cl(num_rat,liste_clust,fr_hemi)

liste_fonc <- list('ADC','BVf','CBF','SO2map','T1map','VSI')
liste_min_fonc <- list('ADC'=50,'BVf'=1,'CBF'=5,'SO2map'=5,'T1map'=10,'VSI'=1)#,'CMRO2'=10)
comp_clust_volADC00(num_rat,liste_suivi_slice,c(1:3),liste_min_fonc,'')

comp_clust_volADC00(num_rat,liste_suivi_slice,c(2,2),liste_min_fonc,'pdf')

#-------------------------------- Etape 4 : --------------------------------#

liste_suivi_slice <- list('ADC'=list(10),
                          'BVf'=list(10),
                          'CBF'=list(10),
                          'CMRO2'=list(10),
                          'SO2map'=list(10),
                          'T1map'=list(10),
                          'VSI'=list(10)
)

#liste_jf <- #--------> adaptée à l'éventuelle absence de segmentation
liste_fonc <- list('ADC','BVf','CBF','CMRO2','VSI','T1map','SO2map')# on peut tester avec une liste des fonctionnalités restreinte

#dgris_temp_fonc(num_rat,'cer',"")
dgris_temp_fonc(num_rat,'tranches',"")
#dgris_temp_fonc(num_rat,10,"")


dgris_temp_fonc(num_rat,'dark','')

ngris_box_fonc(num_rat,fr_hemi,'ADC',liste_suivi_slice)


ddd <- ngris_box_fonc(num_rat,fr_hemi,'ADC',liste_suivi_slice)

ppp <- ggplot(ddd,
            aes(x=ddd$Jour,y=ddd[,3],fill=ddd$Zone
            )
)
ppp <- ppp + geom_boxplot(outlier.shape = NA)
ppp <- ppp + scale_fill_manual(values = alpha(c("grey70","red","blue"), .3))
ppp <- ppp + ggtitle(sprintf("Evolution de %s, segmentation %s",'SO2map','ADC')) + xlab("Jours") + ylab('VSI')
print(ppp)

#-------------------------------- Etape 5 : --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')

suivi_etendue_fonc(num_rat,'cer')

#suivi_etendue_fonc(num_rat,'tranches')

suivi_etendue_fonc(num_rat,10)

###################################################### Rat numéro 19 ######################################################

num_rat <- "19"

#-------------------------------- Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. --------------------------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

# on continue avec une segmentation manuelle

vol_lesADC00(num_rat,c(9:10))

# Un seul jour pour la segmentation CBF
#vol_lesCBF00(num_rat,c(10))
vol_lesCBF00(num_rat,c(9:10))

#-------------------------------- Etape 2 : répertoires des fonctionnalités. --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','T1map','SO2map','VSI') # ordre de priorité pour les fonctionnalités.
liste_min_fonc <- list('ADC'=10,'BVf'=1,'CBF'=10,'CMRO2'=10,'SO2map'=10,'T1map'=10,'VSI'=10)

cl <- c(3,5)# encadrement du nombre de clusters
#cl <- c(4,4)
for (fonc in liste_fonc){
  #rg_FONC_3d(2,fonc,num_rat,cl)
  rg_FONC_3d(3,fonc,num_rat,cl,'')
}

# Remplir liste_clust et liste_suivi_slice en observant les représentations graphiques obtenues

# Vérification : aprè savoir rempli liste_clust

fr_hemi <- c(5,-240) # inter-hémisphère, supposé plan


fonc <- 'ADC'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'BVf'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CBF'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CMRO2'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'SO2map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'T1map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'VSI'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

# trop lente
#for (fonc in liste_fonc){
#  cl_se <- liste_clust[[fonc]]
#  seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)
#}

#-------------------------------- Etape 3 : --------------------------------#

liste_clust <- list('ADC'=list("00"=c(1),"03"=c(),"08"=c(3),"15"=c(4),"22"=c(4)), # clusters utilisables, cf étape 2
                    'BVf'=list("00"=c(1),"08"=c(3),"15"=c(3),"22"=c(0,1)),
                    'CBF'=list("00"=c(1),"03"=c(),"08"=c(1),"15"=c(1),"22"=c(1)), # 3 clusters
                    'CMRO2'=list("00"=c(0),"08"=c(3),"15"=c(1),"22"=c(0)),
                    'SO2map'=list("00"=c(1),"08"=c(1),"15"=c(),"22"=c(1)),
                    'T1map'=list("00"=c(3),"03"=c(3),"08"=c(4),"15"=c(4),"22"=c(4)),
                    'VSI'=list("00"=c(1),"08"=c(1),"15"=c(1),"22"=c(1))
                    )

record_seg_cl(num_rat,liste_clust,fr_hemi)

liste_fonc <- list('ADC','BVf','CBF','SO2map','T1map','VSI')
liste_min_fonc <- list('ADC'=50,'BVf'=1,'CBF'=5,'SO2map'=5,'T1map'=10,'VSI'=1)#,'CMRO2'=10)
comp_clust_volADC00(num_rat,liste_suivi_slice,c(1:3),liste_min_fonc,'')

comp_clust_volADC00(num_rat,liste_suivi_slice,c(2,2),liste_min_fonc,'pdf')

#-------------------------------- Etape 4 : --------------------------------#

liste_suivi_slice <- list('ADC'=list(9,10),
                          'BVf'=list(9,10),
                          'CBF'=list(9),
                          'CMRO2'=list(9,10),
                          'SO2map'=list(9,10),
                          'T1map'=list(9),
                          'VSI'=list(9,10)
                          )

#liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')
# --------> on peut tester avec une liste des fonctionnalités restreinte

#dgris_temp_fonc(num_rat,list('cer'),"")

dgris_temp_fonc(num_rat,fr_hemi,'ADCdark00',list(9),"")
dgris_temp_fonc(num_rat,fr_hemi,'CBFdark00',list(9),"")
#dgris_temp_fonc(num_rat,fr_hemi,'brightAnat00',list(9),"")
#dgris_temp_fonc(num_rat,fr_hemi,'T1map...00',list(9),"")

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'pdf')

#-------------------------------- Etape 5 : --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')

suivi_etendue_fonc(num_rat,list('cer'))

suivi_etendue_fonc(num_rat,'tranches')

suivi_etendue_fonc(num_rat,9)

###################################################### Rat numéro 26 ######################################################

num_rat <- "26"

#-------------------------------- Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. --------------------------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

vol_lesADC00(num_rat,c(6,7,8))

# Un seul jour pour la segmentation CBF
vol_lesCBF00(num_rat,c(6,7,8))


#-------------------------------- Etape 2 : répertoires des fonctionnalités. --------------------------------#

liste_fonc <- list('BVf','CBF','CMRO2','T1map','ADC','SO2map','VSI') # ordre de priorité pour les fonctionnalités.
liste_min_fonc <- list('ADC'=10,'BVf'=1,'CBF'=10,'CMRO2'=10,'SO2map'=10,'T1map'=10,'VSI'=10)

cl <- c(3,5)# encadrement du nombre de clusters
cl <- c(4,4)
for (fonc in liste_fonc){
  #rg_FONC_3d(2,fonc,num_rat,cl)
  rg_FONC_3d(3,fonc,num_rat,cl,'')
}

# Remplir liste_clust et liste_suivi_slice en observant les représentations graphiques obtenues

# Vérification : aprè savoir rempli liste_clust

fr_hemi <- c(5,-240) # inter-hémisphère, supposé plan

fonc <- 'ADC'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'BVf'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

fonc <- 'CBF'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CMRO2'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'SO2map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

comp_2vs3d_clust(fonc,num_rat,cl,10)

fonc <- 'T1map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

comp_2vs3d_clust(fonc,num_rat,cl,10)

fonc <- 'VSI'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}





rr <- read.table(sprintf('%s/%s-J%s-%s-bg-all.dat','VSI','11','03','VSI'),header=T)
lrr <- is.na(rr[,'VSI'])
rr <- rr[!lrr,]
cc <- cluster_jfr_fmin(data = rr,cl_min = 4,cl_max = 4,min = 10)

# trop lente
#for (fonc in liste_fonc){
#  cl_se <- liste_clust[[fonc]]
#  seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)
#}

#-------------------------------- Etape 3 : --------------------------------#

liste_clust <- list('ADC'=list("00"=c(),"03"=c(),"08"=c(),"15"=c(),"22"=c()), # clusters utilisables, cf étape 2
                    'BVf'=list("00"=c(1),"03"=c(),"08"=c(),"15"=c(),"22"=c()),
                    'CBF'=list("00"=c(1),"03"=c(4),"08"=c(3),"15"=c(),"22"=c(4)),
                    'CMRO2'=list("00"=c(1),"03"=c(3),"08"=c(3),"15"=c(3),"22"=c(3)),
                    'SO2map'=list("00"=c(),"03"=c(),"08"=c(),"15"=c(),"22"=c()),
                    'T1map'=list("00"=c(),"03"=c(),"08"=c(),"15"=c(),"22"=c()),
                    'VSI'=list("00"=c(4),"03"=c(3),"08"=c(),"15"=c(),"22"=c(4))
)

record_seg_cl(num_rat,liste_clust,fr_hemi)

liste_suivi_slice <- list('ADC'=c(6,7,8),
                          'BVf'=c(6,7,8),
                          'CBF'=c(6,7,8),
                          'CMRO2'=c(6,7,8),
                          'SO2map'=c(6,7,8),
                          'T1map'=c(6,7,8),
                          'VSI'=c(6,7,8)
)

liste_fonc <- list('BVf','CBF','CMRO2','T1map','ADC','SO2map','VSI')#,'Anat')
comp_clust_volADC00(num_rat,liste_suivi_slice,c(1:3),liste_min_fonc,'')

comp_clust_volADC00(num_rat,liste_suivi_slice,c(2,2),liste_min_fonc,'pdf')

#-------------------------------- Etape 4 : --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')# on peut tester avec une liste des fonctionnalités restreinte

# Je choisis les fonctionnalités ADC et T1map comme références

dgris_temp_fonc(num_rat,'cer',"")
#dgris_temp_fonc(num_rat,list('cer'),"T1map")# ----> fonctionne avec toutes les modalités ! Les données isch3d ont le bon format.
#dgris_temp_fonc(num_rat,list('cer'),"ADC")

dgris_temp_fonc(num_rat,6,"")
dgris_temp_fonc(num_rat,7,"")
dgris_temp_fonc(num_rat,8,"")

dgris_temp_fonc(num_rat,'tranches',"")
#dgris_temp_fonc(num_rat,liste_suivi_slice,"T1map")
#dgris_temp_fonc(num_rat,liste_suivi_slice,"ADC")

dgris_temp_fonc(num_rat,'dark',"")

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice)

#-------------------------------- Etape 5 : --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')

suivi_etendue_fonc(num_rat,'cer')

suivi_etendue_fonc(num_rat,'tranches')

suivi_etendue_fonc(num_rat,6)
suivi_etendue_fonc(num_rat,7)
suivi_etendue_fonc(num_rat,8)

###################################################### Rat numéro 30 ######################################################

num_rat <- "30"

#-------------------------------- Etape 1 : répertoires fonctionnel_gris. Jusque là : 'sans' les valeurs manquantes NaN. --------------------------------#

for (fonc in liste_fonc){
  FONC_3d_rat(fonc,num_rat,'avec')
}

#-------------------------------- Etape 2 : répertoires des fonctionnalités. Attention au BVf, J00, slice8. --------------------------------#

liste_fonc <- list('ADC','BVf','CBF','T1map','SO2map','VSI')#,'CMRO2') # ordre de priorité pour les fonctionnalités.
liste_min_fonc <- list('ADC'=50,'BVf'=1,'CBF'=5,'SO2map'=5,'T1map'=10,'VSI'=1)#,'CMRO2'=10)

cl <- c(3,5)# encadrement du nombre de clusters
#cl <- c(4,4)
for (fonc in liste_fonc){
  #rg_FONC_3d(2,fonc,num_rat,cl)
  rg_FONC_3d(3,fonc,num_rat,cl,'')
}

fonc <- 'ADC'
for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

cl <- c(4,4)
for (fonc in liste_fonc){
  #rg_FONC_3d(2,fonc,num_rat,cl)
  rg_FONC_3d(3,fonc,num_rat,cl,'')
}

# Remplir liste_clust et liste_suivi_slice en observant les représentations graphiques obtenues

# Vérification : aprè savoir rempli liste_clust

fr_hemi <- c(5,-240) # inter-hémisphère, supposé plan

fonc <- 'ADC'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'BVf'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CBF'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'CMRO2'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'SO2map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'T1map'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

fonc <- 'VSI'
cl <- c(4,4)
cl_se <- liste_clust[[fonc]]
seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)

for (sl in liste_suivi_slice[[fonc]]){
  comp_2vs3d_clust(fonc,num_rat,cl,sl)
}

# trop lente
#for (fonc in liste_fonc){
#  cl_se <- liste_clust[[fonc]]
#  seg_clust_3d(fonc,num_rat,cl,cl_se,fr_hemi)
#}

#-------------------------------- Etape 3 : --------------------------------#

liste_clust <- list('ADC'=list("00"=c(4),"08"=c(4),"15"=c(4)), # clusters utilisables, cf étape 2
                    'BVf'=list(),
                    'CBF'=list(),
                    'CMRO2'=list(),
                    'SO2map'=list(),
                    'T1map'=list(),
                    'VSI'=list()
)

record_seg_cl(num_rat,liste_clust,fr_hemi)

liste_suivi_slice <- list('ADC'=list(9,10,11,12,13),
                          'BVf'=list(11,12,13),
                          'CBF'=list(11,12,13),
                          'CMRO2'=list(11,12,13),
                          'SO2map'=list(11,12,13),
                          'T1map'=list(11,12,13),
                          'VSI'=list(11,12,13)
)

#-------------------------------- Etape pi : suivi sur les taches d'ADC, observées au jour 00 --------------------------------#

tache_ADC00_R30 <- function(jour,fonc){
  
  d.filename <- sprintf("%s/liste_R%s_%s_J00.csv",fonc,"30",fonc)
  day.slices <- c(11:13)#read.csv(d.filename,check.names=F,header=T)
  d <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(d) <- c("x","y","z",fonc,"Slice")

  for (slice in day.slices){
    d.filename <- sprintf("%s/%s-J%s-%s-dark-slice%i.txt",'ADC',"30",jour,fonc,slice)
    d.increment <- read.table(d.filename,header=T,sep='\t')
    d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*slice, d.increment[,3],slice))
    colnames(d.increment) <- c("x","y","z",fonc,"Slice")
    d <- as.data.frame(rbind(d,d.increment))
  }
  
  return(d)
}

fonc <- 'ADC'
for (jour in jours_R30_ADC){
  t <- tache_ADC00_R30(jour,'ADC')
  write.table(t, sprintf("%s/%s-J%s-%s-dark-all.dat",fonc,"30",jour,fonc), row.names=F, quote=F, sep='\t')
}

cl <- c(4,4)
for (fonc in list('ADC')){
  #rg_FONC_3d(2,fonc,num_rat,cl)
  rg_FONC_3d(3,fonc,num_rat,cl,'dark')
}

#-------------------------------- Etape 1+e : suivi sur les taches d'ADC, observées au jour 00 --------------------------------#

# Un seul jour pour la segmentation CBF
vol_lesCBF00(num_rat,c(11:13))

liste_suivi_slice <- list('ADC'=list(11,12,13),
                          'BVf'=list(11,12,13),
                          'CBF'=list(11,12,13),
                          'CMRO2'=list(11,12,13),
                          'SO2map'=list(11,12,13),
                          'T1map'=list(11,12,13),
                          'VSI'=list(11,12,13)
)

liste_fonc <- list('ADC','BVf','CBF','SO2map','T1map','VSI')#,'Anat','CMRO2')
liste_min_fonc <- list('ADC'=50,'BVf'=1,'CBF'=5,'SO2map'=5,'T1map'=10,'VSI'=1)

comp_clust_volADC00(num_rat,liste_suivi_slice,c(1:3),liste_min_fonc,'')

comp_clust_volADC00(num_rat,liste_suivi_slice,c(2,2),liste_min_fonc,'pdf')

#-------------------------------- Etape 4 : --------------------------------#

ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc,'ADCdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc,'ADCdark00', liste_suivi_slice, 'pdf')

dgris_temp_fonc(num_rat,'dark','')

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice)

#-------------------------------- Etape 5 : --------------------------------#



##################################################################################################
# -------------------------------- Etape 4 : tous rats confondus-------------------------------- #
##################################################################################################

fr_hemi <- c(5,-240)

liste_fonc <- list('VSI','T1map','SO2map','ADC','BVf','CBF')#,'Anat','CMRO2')

# Pour des clusterisations sur des zones segmentées
liste_fonc <- list('ADC','VSI','T1map','SO2map','CMRO2','CBF','BVf')
liste_fonc <- list('CBF','BVf',
                   'VSI',
                   'SO2map','CMRO2'#,
                   #'T1map'
                   ,'ADC'
                   )
###################################################

liste_fonc <- list('VSI','T1map','SO2map','CMRO2','CBF','BVf','ADC')
liste_min_fonc <- list('ADC'=50,'BVf'=1,'CBF'=5,'SO2map'=5,'T1map'=10,'VSI'=1)

l_s_s_11 <- list('ADC'=list(10),
                 'BVf'=list(10),
                 'CBF'=list(10),
                 'CMRO2'=list(10),
                 'SO2map'=list(10),
                 'T1map'=list(10),
                 'VSI'=list(10)
)
l_s_s_19 <- list('ADC'=list(9,10),
                 'BVf'=list(9,10),
                 'CBF'=list(9),
                 'CMRO2'=list(9,10),
                 'SO2map'=list(9,10),
                 'T1map'=list(9),
                 'VSI'=list(9,10)
)
l_s_s_26 <- list('ADC'=c(6,7,8),
                 'BVf'=c(6,7,8),
                 'CBF'=c(6,7,8),
                 'CMRO2'=c(6,7,8),
                 'SO2map'=c(6,7,8),
                 'T1map'=c(6,7,8),
                 'VSI'=c(6,7,8)
)
l_s_s_30 <- list('ADC'=list(11,12,13),
                 'BVf'=list(11,12,13),
                 'CBF'=list(11,12,13),
                 'CMRO2'=list(11,12,13),
                 'SO2map'=list(11,12,13),
                 'T1map'=list(11,12,13),
                 'VSI'=list(11,12,13)
)

liste_ss_rat <- list("11"=l_s_s_11,"19"=l_s_s_19,"26"=l_s_s_26,"30"=l_s_s_30)

liste_les_jf11_cl <- list('CBF'=list("00"=c(1,1),"03"=c(1,1),"08"=c(1,1),"15"=c(1,2),"22"=c(1,2)),# J03 : clusters ?
                          'CMRO2'=list("00"=c(1,1),"03"=c(1,1),"08"=c(1,1),"15"=c(1,1),"22"=c(1,1)),
                          'SO2map'=list("00"=c(1,1),"03"=c(1,1),"08"=c(1,1),"15"=c(1,2),"22"=c(1,2)),# J03 compatible avec le CBF.
                          'BVf'=list("00"=c(1,1),"03"=c(1,2),"08"=c(1,2),"15"=c(1,2),"22"=c(1,1))
)
liste_les_jf19_cl <- list('CBF'=list("00"=c(1,1),"03"=c(1,1),"08"=c(1,2),"15"=c(1,2),"22"=c(1,2)),
                          'CMRO2'=list("00"=c(1,1),"03"=c(1,1),"08"=c(1,2),"15"=c(1,1),"22"=c(1,2)),
                          'SO2map'=list("00"=c(1,1),"03"=c(1,1),"08"=c(1,1),"15"=c(1,2),"22"=c(1,2))
)
liste_les_jfr_cl <- list("11"=liste_les_jf11_cl,"19"=liste_les_jf19_cl,"26"='',"30"='')

num_rat <- "11"
liste_suivi_slice <- liste_ss_rat[[num_rat]]

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, 'pdf')

ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,3), liste_min_fonc, 'CBFdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'CBFdark00', liste_suivi_slice, 'pdf')

comp_clust_vol00(num_rat,liste_suivi_slice,c(1,2),liste_min_fonc,'CBFdark00','')

ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00',"15",'CBF',liste_min_fonc,liste_suivi_slice,'')
ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00',"15",'SO2map',liste_min_fonc,liste_suivi_slice,'')
ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00',"08",'BVf',liste_min_fonc,liste_suivi_slice,'')

num_rat <- "19"
liste_suivi_slice <- liste_ss_rat[[num_rat]]

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, 'pdf')

dgris_temp_fonc(num_rat,fr_hemi,'CBFdark00',list(9),"00","")
dgris_temp_fonc(num_rat,fr_hemi,'CBFdark00',list(9),"22","")
ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'pdf')

ngris_box_clust(num_rat, fr_hemi, c(1,3), liste_min_fonc, 'CBFdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'CBFdark00', liste_suivi_slice, 'pdf')

comp_clust_vol00(num_rat,liste_suivi_slice,c(1,2),liste_min_fonc,'ADCdark00','')
# Attention avec liste_fonc : commencer avec la modalité de segmentation.
comp_clust_vol00(num_rat,liste_suivi_slice,c(1,2),liste_min_fonc,'CBFdark00','')

suivi_hemC_rat19 <- function(hemi,fonc){
  num_rat <- "19"
  liste_jr <- liste_jfr[[fonc]]
  jours <- liste_jr[[num_rat]]
  for (jour in jours){
    carte_fonc_rat(num_rat,fonc,hemi,jour,'')
  }
}

suivi_hemC_rat19(fr_hemi,'ADC')
suivi_hemC_rat19(fr_hemi,'CBF')

suivi_hemC_rat19(fr_hemi,'VSI')
suivi_hemC_rat19(fr_hemi,'BVf')
suivi_hemC_rat19(fr_hemi,'SO2map')
suivi_hemC_rat19(fr_hemi,'CMRO2')


liste_fonc <- list('ADC','BVf','CBF','CMRO2','SO2map','T1map','VSI')
ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00','','CBF',liste_min_fonc,liste_suivi_slice,'')
#ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00','','CBF,liste_min_fonc,liste_suivi_slice,'pdf')

ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00',"08",'CBF',liste_min_fonc,liste_suivi_slice,'')
ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00',"08",'CMRO2',liste_min_fonc,liste_suivi_slice,'')
ngris_box_fonc_cl_les(num_rat,fr_hemi,'CBFdark00',"15",'SO2map',liste_min_fonc,liste_suivi_slice,'')

num_rat <- "26"
liste_suivi_slice <- liste_ss_rat[[num_rat]]

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, 'pdf')

ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,3), liste_min_fonc, 'CBFdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'CBFdark00', liste_suivi_slice, 'pdf')

comp_clust_vol00(num_rat,liste_suivi_slice,c(1,2),liste_min_fonc,'CBFdark00','')

num_rat <- "30"
liste_suivi_slice <- liste_ss_rat[[num_rat]]

ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'ADCdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'ADCdark00', liste_suivi_slice, 'pdf')

ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'')
ngris_box_fonc(num_rat,fr_hemi,'CBFdark00',liste_suivi_slice,'pdf')
ngris_box_clust(num_rat, fr_hemi, c(1,3), liste_min_fonc, 'CBFdark00', liste_suivi_slice, '')
ngris_box_clust(num_rat, fr_hemi, c(1,2), liste_min_fonc, 'CBFdark00', liste_suivi_slice, 'pdf')

comp_clust_vol00(num_rat,liste_suivi_slice,c(1,2),liste_min_fonc,'CBFdark00','')

# Répertoire benjamin_antoine_labo
liste_rats_tr <- list("11"=c(10),"19"=c(9),"26"=c(8),"30"=c(11))

comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'ADC','CBFdark00','')
comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'BVf','CBFdark00','')
comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'CBF','CBFdark00','')
comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'CMRO2','CBFdark00','')
comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'SO2map','CBFdark00','')
comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'VSI','CBFdark00','')

#comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'ADC','CBFdark00','pdf')
#comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'BVf','CBFdark00','pdf')
#comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'CBF','CBFdark00','pdf')
comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'CMRO2','CBFdark00','pdf')
#comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'SO2map','CBFdark00','pdf')
#comp_rats_clust_fonc(fr_hemi,c(1,2),liste_min_fonc,liste_rats_tr,'VSI','CBFdark00','pdf')



## Une autre version
#carte_fonc_rat(fonc,fr_hemi,"00",'')
#carte_fonc_rat(fonc,fr_hemi,"15",'')


#carte_fonc_rat(fonc,fr_hemi,'Suivi','')# --- > pas au point.
#carte_fonc_rat(fonc,fr_hemi,'Suivi','pdf')

comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'ADC','')
comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'BVf','')
comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'CBF','')
comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'SO2map','')
comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'T1map','')
comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'VSI','')

comp_seg_ADCvsCBF(fr_hemi,liste_rats_tr,'CMRO2','')

##############################################################################################################################################
###----------------- Simulations : modèle phénoménologique pour décrire les suites de l'ischémie, sur les jours d'examen. -----------------###
##############################################################################################################################################


###----------------- Premières tentatives : on note les indicateurs de position pour les fonctionnalités. Rat 11. -----------------###


source(file="modele_pheno_fun.R", encoding ="UTF-8")

fr_hemi <- c(5,-240)

opt_seg <- 'ADCdark00'
opt_seg <- 'CBFdark00'

num_rat <- "11"
liste_suivi_slice <- liste_ss_rat[[num_rat]]

comp_succ_suivi(num_rat,fr_hemi,'ADC',liste_suivi_slice,opt_seg,'box','')
comp_succ_suivi(num_rat,fr_hemi,'BVf',liste_suivi_slice,opt_seg,'box','')
comp_succ_suivi(num_rat,fr_hemi,'CBF',liste_suivi_slice,opt_seg,'box','')
comp_succ_suivi(num_rat,fr_hemi,'SO2map',liste_suivi_slice,opt_seg,'box','')
comp_succ_suivi(num_rat,fr_hemi,'T1map',liste_suivi_slice,opt_seg,'box','')
comp_succ_suivi(num_rat,fr_hemi,'VSI',liste_suivi_slice,opt_seg,'box','')

# il faut définir des fonctions de transition pour --- > # CMRO2 # comp_succ_suivi(num_rat,fr_hemi,'CMRO2',liste_suivi_slice,opt_seg,'box','')


comp_succ_suivi(num_rat,fr_hemi,'ADC',liste_suivi_slice,'ADCdark00','dens','')
comp_succ_suivi(num_rat,fr_hemi,'ADC',liste_suivi_slice,'ADCdark00','dens','pdf')

comp_succ_suivi(num_rat,fr_hemi,'BVf',liste_suivi_slice,'ADCdark00','dens','')
comp_succ_suivi(num_rat,fr_hemi,'CBF',liste_suivi_slice,'ADCdark00','dens','')
comp_succ_suivi(num_rat,fr_hemi,'SO2map',liste_suivi_slice,'ADCdark00','dens','')
comp_succ_suivi(num_rat,fr_hemi,'T1map',liste_suivi_slice,'ADCdark00','dens','')
comp_succ_suivi(num_rat,fr_hemi,'VSI',liste_suivi_slice,'ADCdark00','dens','')

comp_succ_suivi(num_rat,fr_hemi,'BVf',liste_suivi_slice,'ADCdark00','dens','pdf')
comp_succ_suivi(num_rat,fr_hemi,'CBF',liste_suivi_slice,'ADCdark00','dens','pdf')
comp_succ_suivi(num_rat,fr_hemi,'SO2map',liste_suivi_slice,'ADCdark00','dens','pdf')
comp_succ_suivi(num_rat,fr_hemi,'T1map',liste_suivi_slice,'ADCdark00','dens','pdf')
comp_succ_suivi(num_rat,fr_hemi,'VSI',liste_suivi_slice,'ADCdark00','dens','pdf')



#suivi_voxels(num_rat,'ADC',10,c(50,30),7,'ADCdark00')

aff_suivi_voxels(num_rat,'ADC',10,c(50,30),7,'voxels','')
aff_suivi_voxels(num_rat,'ADC',10,c(50,30),7,'hist','')

aff_suivi_voxels(num_rat,'BVf',10,c(50,30),7,'voxels','')
aff_suivi_voxels(num_rat,'BVf',10,c(50,30),7,'hist','')

aff_suivi_voxels(num_rat,'CBF',10,c(60,40),7,'voxels','')
aff_suivi_voxels(num_rat,'CBF',10,c(60,40),7,'hist','')

aff_suivi_voxels(num_rat,'SO2map',10,c(50,30),7,'voxels','')
aff_suivi_voxels(num_rat,'SO2map',10,c(50,30),7,'hist','')

aff_suivi_voxels(num_rat,'T1map',10,c(60,40),7,'voxels','')
aff_suivi_voxels(num_rat,'T1map',10,c(60,40),7,'hist','')

aff_suivi_voxels(num_rat,'VSI',10,c(60,40),7,'voxels','')
aff_suivi_voxels(num_rat,'VSI',10,c(60,40),7,'hist','')