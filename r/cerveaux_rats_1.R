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

color.vector <- c('darkblue','red','green','purple', 'orange')

#-------------------------------------------
# Code nico


# combiner les fichiers
d.day <- "03"
d.slices.d03 <- c(8,9,10,11)
d.slice.size <- 1/14*10000 # nb max de slices + 1 * 10000 --> en µm
slice.select <- 1


#-*-*-
# commencer par initialiser une data.frame "d" avec 4 colonnes (x, y, z, adc)
# Boucle sur slice.select
# charge fichier à slice.select dans d.slice
# bind la d.slice avec z
# rbind d <- rbind(d,d.slice)
# reboucle
# enregistre d
# charge d // pour tester
# clust ...

d8.filename <- sprintf("11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d8 <- read.table(d8.filename, header=F, sep='\t')
d8 <- as.data.frame(cbind(d8[,1:2],z=d.slice.size*d.slices.d03[slice.select], d8[,3]))
colnames(d8) <- c("x","y","z","adc")
write.table(d8, sprintf("11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

slice.select <- 2
d9.filename <- sprintf("11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d9 <- read.table(d9.filename, header=F, sep='\t')
d9 <- as.data.frame(cbind(d9[,1:2],z=d.slice.size*d.slices.d03[slice.select], d9[,3]))
colnames(d9) <- c("x","y","z","adc")
write.table(d9, sprintf("11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

slice.select <- 3
d10.filename <- sprintf("11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d10 <- read.table(d10.filename, header=F, sep='\t')
d10 <- as.data.frame(cbind(d10[,1:2],z=d.slice.size*d.slices.d03[slice.select], d10[,3]))
colnames(d10) <- c("x","y","z","adc")
write.table(d10, sprintf("11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

slice.select <- 4
d11.filename <- sprintf("11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d11 <- read.table(d11.filename, header=F, sep='\t')
d11 <- as.data.frame(cbind(d11[,1:2],z=d.slice.size*d.slices.d03[slice.select], d11[,3]))
colnames(d11) <- c("x","y","z","adc")
write.table(d11, sprintf("11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

d <- as.data.frame(rbind(d10,d11))
write.table(d, sprintf("11-J%s-ADC-bg-all.dat",d.day), row.names=F, quote=F, sep='\t')
#-*-*-


#d <- read.table("ADC/11-J03-ADC-bg-all.txt", header=F, sep='\t')
d <- d[d$adc>10,] # cutoff à 10 de diffusion
d.clust <- Mclust(d$adc, G=2:5, modelNames="V")
plot(d.clust, what="classification")

par(mfrow = c(2,2))
plot(d.clust, what="BIC")
plot(d.clust, what="classification", col=color.vector)
adc.breaks <- seq(min(d$adc)-0.1*min(d$adc), max(d$adc)+0.1*max(d$adc), length.out=100)
d.hist <- hist(d$adc,breaks=adc.breaks, col='grey50')
plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, cex=2*(1-d.clust$uncertainty)^4, xlab='x', ylab='y')



#-------------------------------------------


# Jour 08
# combiner les fichiers

d.day <- "08"
d.slices.d08 <- c(10,11,12)
d.slice.size <- 1/14*10000 # nb max de slices + 1 * 10000 --> en µm
#slice.select <- 1

d <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(d) <- c("x","y","z","adc")


for (slice.select in 1:length(d.slices.d08)){
  #slice.select <- i
  d.filename <- sprintf("11-J%s-ADC-bg-slice%g.txt",d.day,d.slices.d08[slice.select])
  d.increment <- read.table(d.filename,header=F,sep='\t')
  d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*d.slices.d08[slice.select], d.increment[,3]))
  colnames(d.increment) <- c("x","y","z","adc")
  write.table(d.increment, sprintf("11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d08[slice.select]), row.names=F, quote=F, sep='\t')
  d<-as.data.frame(rbind(d,d.increment))
}

write.table(d, sprintf("11-J%s-ADC-bg-all.dat",d.day), row.names=F, quote=F, sep='\t')

d <- d[d$adc>10,] # cutoff à 10 de diffusion
d.clust <- Mclust(d$adc, G=2:4, modelNames="V")
plot(d.clust, what="classification")

par(mfrow = c(2,2))
plot(d.clust, what="BIC")
plot(d.clust, what="classification", col=color.vector)
adc.breaks <- seq(min(d$adc)-0.1*min(d$adc), max(d$adc)+0.1*max(d$adc), length.out=100)
d.hist <- hist(d$adc,breaks=adc.breaks, col='grey50')
plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, cex=2*(1-d.clust$uncertainty)^4, xlab='x', ylab='y')

#---------------------------- base de données : rat 11 ; ADC

liste_R11_ADC <- list("00"=c(8,9,10),"03"=c(8,9,10,11),"08"=c(10,11,12),"15"=c(10,11,12,13),"22"=c(10,11,12))


#---------------------------- une fonction à appliquer le jour 15


d.slice.size <- 1/14*10000


cerveau_jour <- function(day){
  d.slices.day <- liste_R11_ADC[[day]]
  d <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(d) <- c("x","y","z","adc")
  for (slice.select in 1:length(d.slices.day)){
    d.filename <- sprintf("11-J%s-ADC-bg-slice%g.txt",day,d.slices.day[slice.select])
    d.increment <- read.table(d.filename,header=F,sep='\t')
    d.increment <- as.data.frame(cbind(d.increment[,1:2],z=d.slice.size*d.slices.day[slice.select], d.increment[,3]))
    colnames(d.increment) <- c("x","y","z","adc")
    write.table(d.increment, sprintf("11-J%s-ADC-bg-slice%i.dat",day,d.slices.day[slice.select]), row.names=F, quote=F, sep='\t')
    d<-as.data.frame(rbind(d,d.increment))
  }
  #write.table(d, sprintf("11-J%s-ADC-bg-all.dat",day), row.names=F, quote=F, sep='\t')
  return(d)
}

d <- cerveau_jour("15")
write.table(d, sprintf("11-J%s-ADC-bg-all.dat","15"), row.names=F, quote=F, sep='\t')

d <- d[d$adc>10,] # cutoff à 10 de diffusion
d.clust <- Mclust(d$adc, G=2:4, modelNames="V")
plot(d.clust, what="classification")

par(mfrow = c(2,2))
plot(d.clust, what="BIC")
plot(d.clust, what="classification", col=color.vector)
adc.breaks <- seq(min(d$adc)-0.1*min(d$adc), max(d$adc)+0.1*max(d$adc), length.out=100)
d.hist <- hist(d$adc,breaks=adc.breaks, col='grey50')
plot(d$x, d$y, col=color.vector[d.clust$classification], pch=20, cex=2*(1-d.clust$uncertainty)^4, xlab='x', ylab='y')










# Fonctionnel :

var<-'ADC' ######

# Spatial : S00, S03 etc


# Temporel : JS[1] et JS[2]



var<-'BVf' ######


