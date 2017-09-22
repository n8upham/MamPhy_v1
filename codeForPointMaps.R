----------------
title: "Untitled"
author: "Adam M. Wilson"
date: "October 31, 2016"
output: html_document
----------------
https://github.com/cmerow/RDataScience/blob/gh-pages/07_Reproducible.md


library(dplyr)
library(ggplot2)
library(maps)
library(spocc)



## define which species to query
sp='Turdus migratorius'

## run the query and convert to data.frame()
d = occ(query=sp, from='ebird',limit = 100) %>% occ2df()

# Load coastline
map=map_data(map="world",region="south america")

map('worldHires')
map('world2Hires')
map.scale(160,0,relwidth = 0.15, metric = TRUE, ratio = TRUE)
map('worldHires','Italy')
map('worldHires','Argentina')

map('worldHires',
    c('UK', 'Ireland', 'Isle of Man','Isle of Wight', 'Wales:Anglesey'))
map('worldHires',
		c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
		xlim=c(-11,3), ylim=c(49,60.9))	



d = occ(query=c("Tympanoctomys","Octomys"), from='gbif',limit = 100) %>% occ2df()

# Plot it
ggplot(d,aes(x=longitude,y=latitude))+
  geom_polygon(aes(x=long,y=lat,group=group,order=order),data=map)+
  geom_point(col="red")+
  coord_equal()

e<-as.data.frame(d)
f<-na.omit(e[,1:4])

Tympa<-na.omit(e[1:38,1:4])
TympaPT<-Tympa[,2:3]
Octomys<-na.omit(e[39:54,1:4])
OctomysPT<-Octomys[,2:3]

pdf(file="map_TympaOctomys_Argentina.pdf")
map('worldHires','Argentina')
points(TympaPT,col=2,pch=18)
points(OctomysPT,col=3,pch=18)
dev.off()


                     name longitude  latitude prov
1  Tympanoctomys barrerae -67.19972 -33.09361 gbif
2  Tympanoctomys barrerae -67.19972 -33.09361 gbif
3  Tympanoctomys barrerae -67.19972 -33.09361 gbif
7  Tympanoctomys barrerae -67.20364 -33.16068 gbif
8  Tympanoctomys barrerae -67.20364 -33.16068 gbif
18 Tympanoctomys barrerae -67.20364 -33.16068 gbif
20 Tympanoctomys barrerae -69.58250 -35.48217 gbif
21 Tympanoctomys barrerae -69.58250 -35.48217 gbif
23 Tympanoctomys barrerae -68.66119 -35.17008 gbif
25 Tympanoctomys barrerae -69.58250 -35.48217 gbif
27 Tympanoctomys barrerae -69.58250 -35.48217 gbif
32 Tympanoctomys barrerae -69.58250 -35.48217 gbif
35 Tympanoctomys barrerae -68.66119 -35.17008 gbif
36 Tympanoctomys barrerae -69.58250 -35.48217 gbif
40          Octomys mimax -68.73333 -31.98333 gbif
45          Octomys mimax -67.97423 -30.07373 gbif
46          Octomys mimax -67.00250 -32.49633 gbif
47          Octomys mimax -67.97423 -30.07373 gbif
48          Octomys mimax -67.00250 -32.49633 gbif
49          Octomys mimax -67.97423 -30.07373 gbif
52          Octomys mimax -68.73333 -31.98333 gbif
53          Octomys mimax -67.97423 -30.07373 gbif
54          Octomys mimax -67.97423 -30.07373 gbif

