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

d8.filename <- sprintf("ADC/ADC/11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d8 <- read.table(d8.filename, header=F, sep='\t')
d8 <- as.data.frame(cbind(d8[,1:2],z=d.slice.size*d.slices.d03[slice.select], d8[,3]))
colnames(d8) <- c("x","y","z","adc")
write.table(d8, sprintf("ADC/11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

slice.select <- 2
d9.filename <- sprintf("ADC/ADC/11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d9 <- read.table(d9.filename, header=F, sep='\t')
d9 <- as.data.frame(cbind(d9[,1:2],z=d.slice.size*d.slices.d03[slice.select], d9[,3]))
colnames(d9) <- c("x","y","z","adc")
write.table(d9, sprintf("ADC/11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

slice.select <- 3
d10.filename <- sprintf("ADC/ADC/11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d10 <- read.table(d10.filename, header=F, sep='\t')
d10 <- as.data.frame(cbind(d10[,1:2],z=d.slice.size*d.slices.d03[slice.select], d10[,3]))
colnames(d10) <- c("x","y","z","adc")
write.table(d10, sprintf("ADC/11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

slice.select <- 4
d11.filename <- sprintf("ADC/ADC/11-J%s-ADC-bg-slice%i.txt",d.day,d.slices.d03[slice.select])
d11 <- read.table(d11.filename, header=F, sep='\t')
d11 <- as.data.frame(cbind(d11[,1:2],z=d.slice.size*d.slices.d03[slice.select], d11[,3]))
colnames(d11) <- c("x","y","z","adc")
write.table(d11, sprintf("ADC/11-J%s-ADC-bg-slice%i.dat",d.day,d.slices.d03[slice.select]), row.names=F, quote=F, sep='\t')

d <- as.data.frame(rbind(d8,d9, d10,d11))
write.table(d, sprintf("ADC/11-J%s-ADC-bg-all.dat",d.day), row.names=F, quote=F, sep='\t')
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


















clust_slice <- function(Mclust,day,slice){
    ngr_file_name<-paste("11-J",day,"-ADC-bg-slice",slice,sep='')
    data=read.table(paste(ngr_file_name,".txt",sep=''))
    names(data)[1]<-"x"
    names(data)[2]<-"y"
    names(data)[3]<-"v"
    data.clust<-Mclust(data=data$v,G=Nclust,modelNames='V')
    # On trace et on sauvegarde
    png(file = paste(ngr_file_name,".png",sep=''), width = 800, height = 700)
    par(mfrow=c(1,2))#, bg='white')
    plot(data$x,data$y,col=color.vector[data.clust$classification],pch=20,cex=(1-data.clust$uncertainty)^3,main=paste("Brain map of",vars,"data classification",sep=' '),xlab="x",ylab="y")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    points(data$x,data$y,col=color.vector[data.clust$classification], pch=20, cex=(1-data.clust$uncertainty)^3)
    plot(data.clust, what="classification", colors=color.vector, xlab=vars)
    dev.off()
}




cc<-clust_slice(4,"03","10")






















# Fonctionnel :

var<-'ADC' ######

# Spatial : S00, S03 etc


# Temporel : JS[1] et JS[2]



var<-'BVf' ######


