##Análisis

#Librerías
library(readxl)
library(sp)
library(sf)
library(tidyverse)
library(rgdal)
library(RColorBrewer)
library(spdep)
library(tmap)
library(tmaptools)
library(spatialreg)
library(epitools)
library(DCluster)
library(plotrix)
library(MASS)
library(mgcv)


indicadores <- read_excel("Datos/Datos Cantonales.xlsx")
cantones_sp <- read_sf(dsn="Datos",layer = "Cantones")

#Merge
datos_sf <- merge(cantones_sp,indicadores)
datos_sf <- datos_sf[,-c(1:4)]
rm(cantones_sp,indicadores)

#Quitar la Isla del Coco
new_bb = c(286803.0, 889158.2, 658864.2,1241118.1)
names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
attr(new_bb, "class") = "bbox"
attr(st_geometry(datos_sf), "bbox") <- new_bb
rm(new_bb)

#Primer ploteo
#Estadísticas descriptivas
pdf("Trabajo Escrito/F11.pdf")

tm_shape(datos_sf) +
  tm_polygons("dengue", palette=c("lightgreen","tomato"), legend.hist=TRUE)

dev.off()

pdf("Trabajo Escrito/F12.pdf")

tm_shape(datos_sf) +
  tm_polygons("dengue", palette=c("lightgreen","tomato"),style="quantile")

dev.off()

pdf("Trabajo Escrito/FA1.pdf")

tm_shape(datos_sf) +
  tm_polygons("casos", palette=c("lightgreen","tomato"), legend.hist=TRUE)

dev.off()

pdf("Trabajo Escrito/FA2.pdf")

tm_shape(datos_sf) +
  tm_fill("dengue",style="sd",palette=c("lightgreen","tomato")) +
  tm_borders()

dev.off()

pdf("Trabajo Escrito/FA3.pdf")

tm_shape(datos_sf) +
  tm_polygons("tugurio",n=6, palette="-Spectral")

dev.off()

pdf("Trabajo Escrito/FA4.pdf")

tm_shape(datos_sf) +
  tm_polygons("densidad",n=6, palette="-Spectral")

dev.off()

pdf("Trabajo Escrito/FA5.pdf")

tm_shape(datos_sf) +
  tm_polygons("residuos",n=6, palette="Spectral")

dev.off()

pdf("Trabajo Escrito/FA6.pdf")

tm_shape(datos_sf) +
  tm_polygons("acueducto",n=6, palette="Spectral")

dev.off()

datos_sp <- as(datos_sf,"Spatial")
datos_sp@bbox <- matrix(c(286803.0, 889158.2, 658864.2,1241118.1),ncol = 2,byrow = F)

rm(datos_sf)

#Análisis Vecinos
coords <- coordinates(datos_sp)
id <-row.names(datos_sp) 

nb.1 <- poly2nb(datos_sp,queen = T)
nb.2 <- poly2nb(datos_sp,queen = F)
nb.3 <- knn2nb(knearneigh(coords, k=2), row.names=id)
nb.4 <- knn2nb(knearneigh(coords, k=4), row.names=id)

pdf("Trabajo Escrito/F21.pdf")

plot(datos_sp, axes=F, border="gray")
plot(nb.1,coords, pch = 20, cex = 0.6, add = T, col = "red")

dev.off()

pdf("Trabajo Escrito/F22.pdf")

plot(datos_sp, axes=F, border="gray")
plot(nb.4,coords, pch = 20, cex = 0.6, add = T, col = "red")

dev.off()

#Matrices de pesos
w.11 <- nb2listw(nb.1,style = "W")
w.12 <- nb2listw(nb.1,style = "B")
w.13 <- nb2listw(nb.1,style = "S")

w.21 <- nb2listw(nb.2,style = "W")
w.22 <- nb2listw(nb.2,style = "B")
w.23 <- nb2listw(nb.2,style = "S")

w.31 <- nb2listw(nb.3,style = "W")
w.32 <- nb2listw(nb.3,style = "B")
w.33 <- nb2listw(nb.3,style = "S")

w.41 <- nb2listw(nb.4,style = "W")
w.42 <- nb2listw(nb.4,style = "B")
w.43 <- nb2listw(nb.4,style = "S")

#Test de Moran
moran.test(datos_sp$dengue,listw=w.11)
moran.test(datos_sp$dengue,listw=w.12)
moran.test(datos_sp$dengue,listw=w.13)
moran.test(datos_sp$dengue,listw=w.21)
moran.test(datos_sp$dengue,listw=w.22)
moran.test(datos_sp$dengue,listw=w.23)
moran.test(datos_sp$dengue,listw=w.31)
moran.test(datos_sp$dengue,listw=w.32)
moran.test(datos_sp$dengue,listw=w.33)
moran.test(datos_sp$dengue,listw=w.41)
moran.test(datos_sp$dengue,listw=w.42)
moran.test(datos_sp$dengue,listw=w.43)

#Se elije Reina y matriz W
rm(nb.2,nb.3,nb.4,w.12,w.13,w.21,w.22,w.23,w.31,w.32,w.33,w.41,w.42,w.43,coords,id)

#Casos de influencia

pdf("Trabajo Escrito/F31.pdf")

msp <- moran.plot(datos_sp$dengue, listw=w.11, quiet=TRUE,xlab="Tasa de Dengue",ylab = "Tasa espacialmente rezagada")

dev.off()

infl <- apply(msp$is.inf, 1, any)
x <- datos_sp$dengue
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=TRUE)
wx <- stats::lag(w.11, datos_sp$dengue)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4

pdf("Trabajo Escrito/F32.pdf")

plot(datos_sp, col=brewer.pal(4, "Accent")[cols],axes=T,xaxt="n",yaxt="n")
legend("topright", legend=c("Ninguno", "HL", "LH", "HH"), fill=brewer.pal(4, "Accent"), bty="n", cex=0.8, y.intersp=0.8)

dev.off()

rm(msp,cols,infl,lhlh,lhwx,lhx,wx,x)

#Supuestos

m1 <- localmoran(datos_sp$dengue, listw = w.11)
m2 <- as.data.frame(localmoran.sad(lm(dengue ~ 1, datos_sp), nb = nb.1))
m3 <- as.data.frame(localmoran.exact(lm(dengue ~ 1, datos_sp), nb = nb.1))

datos_sp$Normal <- m2[,3]
datos_sp$Aleatorizado <- m1[,5]
datos_sp$Punto_de_silla <- m2[,5]
datos_sp$Exacto <- m3[,5]
gry <- c(rev(brewer.pal(6, "Reds")), brewer.pal(6, "Blues"))

pdf("Trabajo Escrito/FA7.pdf")

spplot(datos_sp, c("Punto_de_silla", "Exacto", "Normal", "Aleatorizado"), 
       at=c(0,0.01,0.05,0.1,0.9,0.95,0.99,1), 
       col.regions=colorRampPalette(gry)(7))

dev.off()

rm(m1,m2,m3,gry)

datos_sp@data <- datos_sp@data[,c(1:8)]

# Modelos

#Lineal
m1 <- lm(dengue~tugurio+densidad+residuos+acueducto,data = datos_sp)
summary(m1)
step(m1)
m1 <- lm(sqrt(dengue)~residuos+acueducto,data = datos_sp)
summary(m1)
#plot(m1)
lm.morantest(m1, listw = w.11)

#SAR
m2 <- spautolm(dengue~tugurio+densidad+residuos+acueducto,data = datos_sp,listw=w.11)
summary(m2)
m2 <- spautolm(dengue~residuos+acueducto,data = datos_sp,listw=w.11)
summary(m2)
moran.mc(residuals(m2),w.11, 999)

#CAR
m3 <- spautolm(dengue~tugurio+densidad+residuos+acueducto,data = datos_sp,listw=w.11,family = "CAR")
summary(m3)
m3 <- spautolm(dengue~residuos+acueducto,data = datos_sp,listw=w.11,family = "CAR")
summary(m3)
moran.mc(residuals(m3),w.11, 999)

#GLM
datos_sp$x<-coordinates(datos_sp)[,1]/1000
datos_sp$y<-coordinates(datos_sp)[,2]/1000
m4 <- gam(as.integer(casos)~+residuos+acueducto+offset(log(pob))+s(x,y), data=datos_sp, family= "quasipoisson")
summary(m4)
moran.mc(residuals(m4),w.11, 999)

#Residuales
datos_sp$Lineal <- residuals(m1)
datos_sp$SAR <- residuals(m2)
datos_sp$CAR <- residuals(m3)
datos_sp$GAM <- residuals(m4)

pdf("Trabajo Escrito/F4.pdf")

spplot(datos_sp,c("SAR", "CAR", "Lineal","GAM"))

dev.off()
#Epidem

datos_sp$observados <- datos_sp$casos
r <- sum(datos_sp$observados)/sum(datos_sp$pob)
datos_sp$esperados <- datos_sp$pob*r
datos_sp$SMR <- datos_sp$observados/datos_sp$esperados

pdf("Trabajo Escrito/FA8.pdf")

spplot(datos_sp,c("observados","esperados"), col.regions=rev(brewer.pal(7, "RdYlGn")), cuts=6)

dev.off()

pdf("Trabajo Escrito/F5.pdf")

spplot(datos_sp,"SMR",col.regions=rev(brewer.pal(7, "RdYlGn")), cuts=6)

dev.off()

int <- pois.exact(datos_sp$SMR)
int <- cbind(int,datos_sp$NOM_CANT_1)
col <- 1*(int$lower>1)
col <- ifelse(col==0,"grey","red")
linea <- ifelse(col=="grey",4,1) 

pdf("Trabajo Escrito/F6.pdf")

plotCI(x = 1:81, y = int$x, ui = int$upper,li = int$lower,pch=18,err="y",
       col=col,sfrac = 0,xlab="Cantones",ylab="Riesgo Relativo",xaxt="n")
abline(h=1,col="grey",lty=2,lwd=1.75)

dev.off()

datos_sp$ch <- choynowski(datos_sp$observados,datos_sp$esperados)$pmap

pdf("Trabajo Escrito/FA9.pdf")

spplot(datos_sp,"ch")

dev.off()

#Empirical Bayes Estimates
eb1 <- EBest(datos_sp$observados,datos_sp$esperados)
unlist(attr(eb1, "parameters"))
datos_sp$EB_mm <- eb1$estmm

res <- empbaysmooth(datos_sp$observados,datos_sp$esperados)
unlist(res[2:3])
datos_sp$EB_ml <- res$smthrr
eb2 <- EBlocal(datos_sp$observados,datos_sp$esperados, nb.1)
datos_sp$EB_mm_local <- eb2$est

pdf("Trabajo Escrito/FA10.pdf")

spplot(datos_sp, c("SMR", "EB_ml", "EB_mm", "EB_mm_local"))

dev.off()

#Fin del análisis