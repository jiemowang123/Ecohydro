setwd("E:/旱情项目")
rm(list = ls(all = TRUE))
source("E:/bishe/multiplot.R")#排序的代码，用于多文件的编程
library(ncdf4)
library(sp)
library(raster) 
library(lubridate)
library(ggplot2)
require(rasterVis)
require(plotrix)
require(gtable)
library(dplyr)
library(Cairo)
require(grid)
require(extraDistr)
library(trend)
library(rgdal)
library(SPEI)
library(PRSim)
library(grDevices)
library(forecast)
library(MASS)
library(lmom)
library(nortest)
library(maptools)
cc<-readOGR("E:/旱情项目/湖北省真实.shp",encoding="utf-8",use_iconv = T)
hb_xian<-readOGR("E:/旱情项目/Export_Output.shp",encoding="utf-8",use_iconv = T)
cc_raster <- fortify(hb_xian)
#cj_1 <- rgdal::readOGR(dsn = cc, stringsAsFactors = FALSE)
#cj_3 <- fortify(cj_1)
# P_world<-stack("cru_ts4.03.1901.2018.pre.dat.nc",varname="pre")
# P_cj<-crop(P_world,cc)
# P_cj<-mask(P_cj,cc)
# PET_world<-stack("cru_ts4.03.1901.2018.pet.dat.nc",varname="pet")
# PET_cj<-crop(PET_world,cc)
# PET_cj<-mask(PET_cj,cc)
# #####只是取1949开始的数据(1948年加入计算)
# P_cj<-P_cj[[565:1416]]
# PET_cj<-PET_cj[[565:1416]]
P_hb<-stack("cru_ts4.03.1901.2018.pre.dat.nc",varname="pre")
PET_hb<-stack("cru_ts4.03.1901.2018.pet.dat.nc",varname="pet")
P_hb<-mask(crop(P_hb,cc),cc)
PET_hb<-mask(crop(PET_hb,cc),cc)
P_hb<-P_hb[[(1416-840+1):1416]]
PET_hb<-PET_hb[[(1416-840+1):1416]]
P_value<-rasterToPoints(P_hb)
PET_value<-rasterToPoints(PET_hb)
P_value<-P_value[,-(1:2)]
PET_value<-PET_value[,-(1:2)]
a_value<-getValues(P_hb[[1]])
number_nona<-which(!is.na(a_value))
SPEI_y_value<-getValues(P_hb)
PET_y_value<-PET_value[,-(1:12)]
#####风险危险性
###多年平均降雨量
Pm<-P_hb[[(840-228+1):840]]
P_mean<-mean(Pm)*12
rtz<-c(expression(paste("108"^"°","E",sep="")),expression(paste("110"^"°","E",sep="")),expression(paste("112"^"°","E",sep="")),expression(paste("114"^"°","E",sep="")),expression(paste("116"^"°","E",sep="")))
ori<-c(expression(paste("29"^"°","N",sep="")),expression(paste("30"^"°","N",sep="")),expression(paste("31"^"°","N",sep="")),expression(paste("32"^"°","N",sep="")),expression(paste("33"^"°","N",sep="")))
gg_P<-rasterToPoints(P_mean)
P_m<-getValues(P_mean)
gg_P<-as.data.frame(gg_P)
PRE<-ggplot(gg_P,aes(x=x,y=y,fill=gg_P[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="降水量(mm)",colours = c("#f7fbff", "#9ecae1","#2171b5"),limits=c(800,1450),breaks = seq(800, 1400,200))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
PRE
fig_name = "PRE_value.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(PRE)
dev.off()
#########
for(i in 1:73)
{
  BAL<-P_value[i,]-PET_value[i,]
  spei1<-spei(BAL,1)
  vxv<-as.numeric(spei1[[2]])
  SPEI_y_value[number_nona[i],]<-vxv
}
SPEI_y_raster<-setValues(P_hb[[1:840]],SPEI_y_value)
####2000~2018
SPEI_m_raster<-SPEI_y_raster[[(840-228+1):840]]
SPEI_m_raster[SPEI_m_raster>=-0.5]<-NA
SPEI_hb<-getValues(SPEI_m_raster)
Frequency<-apply(SPEI_hb,1,function(x)
  {
  a<-length(na.omit(x))/length(x)
  return(a)
})
Frequency[Frequency==0]<-NA
F_SPEI<-setValues(P_hb[[1]],Frequency)
Severity<-apply(SPEI_hb,1,function(x)
{
  a<--0.5-sum(na.omit(x))/length(na.omit(x))
  return(a)
})
Severity[Severity==0]<-NA
S_SPEI<-setValues(P_hb[[1]],Severity)
gg_f<-rasterToPoints(F_SPEI)
gg_f<-as.data.frame(gg_f)
gg_f[,3]<-gg_f[,3]*100
Draw_F<-ggplot(gg_f,aes(x=x,y=y,fill=gg_f[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="干旱频率(%)",colours = c("#4575b4", "#ffffbf","#d73027"),limits=c(25,40),breaks = seq(25, 40,5))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=14),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
Draw_F
fig_name = "frequency.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(Draw_F)
dev.off()
#S
gg_S<-rasterToPoints(S_SPEI)
gg_S<-as.data.frame(gg_S)
Draw_S<-ggplot(gg_S,aes(x=x,y=y,fill=gg_S[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="干旱强度",colours = c("#4575b4", "#ffffbf","#d73027"),limits=c(0.48,0.65),breaks = seq(0.5,0.65,0.05))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=14),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
Draw_S
fig_name = "Severity.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(Draw_S)
dev.off()
####危险性
Frequency<-0.5+0.5*(Frequency-min(Frequency,na.rm=T))/(max(Frequency,na.rm=T)-min(Frequency,na.rm=T))
Severity<-0.5+0.5*(Severity-min(Severity,na.rm=T))/(max(Severity,na.rm=T)-min(Severity,na.rm=T))
P_m<--P_m
P_m<-0.5+0.5*(P_m-min(P_m,na.rm=T))/(max(P_m,na.rm=T)-min(P_m,na.rm=T))
####土地利用
Land_use<-raster("E:/旱情项目/ld20151.tif")
Land_use[Land_use<20]<-1
Land_use[Land_use<30&Land_use>=20]<-2
Land_use[Land_use<40&Land_use>=30]<-3
Land_use[Land_use<50&Land_use>=40]<-4
Land_use[Land_use>=50&Land_use<=60]<-5
Land_use[Land_use>=60]<-6
plot(Land_use)
bb<-rasterToPolygons(S_SPEI)
vv<-extract(Land_use,bb)
water<-0##水域
arable<-0##耕地
tree<-0##林地
grass<-0##草地
city<-0##城镇
zong<-0##总网格数
for(i in 1:73)
{
  water[i]<-length(which(vv[[i]]==4))
  arable[i]<-length(which(vv[[i]]==1))
  tree[i]<-length(which(vv[[i]]==2))
  grass[i]<-length(which(vv[[i]]==3))
  city[i]<-length(which(vv[[i]]==5))
  zong[i]<-length(vv[[i]])
}
water_value<-0##水域
arable_value<-0##耕地
tree_value<-0##林地
grass_value<-0##草地
city_value<-0
zong_value<-0
j=1
for(i in 1:135)
{
  if(i%in%number_nona)
  {
    water_value[i]<-water[j]
    arable_value[i]<-arable[j]
    tree_value[i]<-tree[j]
    grass_value[i]<-grass[j]
    city_value[i]<-city[j]
    zong_value[i]<-zong[j]
    j<-j+1
  }
  else
  {
    water_value[i]<-NA
    arable_value[i]<-NA
    tree_value[i]<-NA
    grass_value[i]<-NA
    zong_value[i]<-NA
    city_value[i]<-NA
    
  }
}
water_value_tj<-water_value
arable_value_tj<-arable_value
tree_value_tj<-tree_value
grass_value_tj<-grass_value
city_value_tj<-city_value
zong_value_tj<-zong_value
###########################
###孕灾环境脆弱性################
#####################
##tu
water_value<--water_value###水域越多，说明受干旱的影响越小
water_value<-0.5+0.5*(water_value-min(water_value,na.rm=T))/(max(water_value,na.rm=T)-min(water_value,na.rm=T))
water_raster<-setValues(S_SPEI,water_value)
gg_water<-rasterToPoints(water_raster)
gg_water<-as.data.frame(gg_water)
mm<-ggplot(gg_water,aes(x=x,y=y,fill=gg_water[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="脆弱性",colours = c("#4575b4", "#ffffbf","#d73027"),limits=c(0.5,1),breaks = seq(0.6, 1,0.1))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=14),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
mm
fig_name = "脆弱性.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(mm)
dev.off()
##水体###

baifenbi_water<-water_value_tj/zong_value
water_raster_tj<-setValues(water_raster,baifenbi_water)
gg_water<-rasterToPoints(water_raster_tj)
gg_water<-as.data.frame(gg_water)
gg_water[,3]<-gg_water[,3]*100
ll<-ggplot(gg_water,aes(x=x,y=y,fill=gg_water[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="水体百分比(%)",colours = c("#f7fbff", "#9ecae1","#2171b5"),limits=c(0,35),breaks = seq(0, 35,7))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
ll
fig_name = "waster_value.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(ll)
dev.off()
##############暴露性
#耕地
baifenbi_arable<-arable_value_tj/zong_value
arable_raster_tj<-setValues(water_raster,baifenbi_arable)
gg_arable<-rasterToPoints(arable_raster_tj)
gg_arable<-as.data.frame(gg_arable)
gg_arable[,3]<-gg_arable[,3]*100
ara<-ggplot(gg_arable,aes(x=x,y=y,fill=gg_arable[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="耕地百分比(%)",colours = c("#ffffe5", "#78c679","#006837"),limits=c(0,100),breaks = seq(0, 100,20))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
ara
fig_name = "arable_value.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(ara)
dev.off()
####林地
baifenbi_tree<-tree_value_tj/zong_value
tree_raster_tj<-setValues(water_raster,baifenbi_tree)
gg_tree<-rasterToPoints(tree_raster_tj)
gg_tree<-as.data.frame(gg_tree)
gg_tree[,3]<-gg_tree[,3]*100
tree<-ggplot(gg_tree,aes(x=x,y=y,fill=gg_tree[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="林地百分比(%)",colours = c("#ffffe5", "#78c679","#006837"),limits=c(0,100),breaks = seq(0, 100,20))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
tree
fig_name = "tree_value.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(tree)
dev.off()
##草地
baifenbi_grass<-grass_value_tj/zong_value
grass_raster_tj<-setValues(water_raster,baifenbi_grass)
gg_grass<-rasterToPoints(grass_raster_tj)
gg_grass<-as.data.frame(gg_grass)
gg_grass[,3]<-gg_grass[,3]*100
grass<-ggplot(gg_grass,aes(x=x,y=y,fill=gg_grass[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="草地百分比(%)",colours = c("#ffffe5", "#78c679","#006837"),limits=c(0,30),breaks = seq(0,30,10))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
grass
fig_name = "grass_value.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(grass)
dev.off()
##城镇
baifenbi_city<-city_value_tj/zong_value
city_raster_tj<-setValues(water_raster,baifenbi_city)
gg_city<-rasterToPoints(city_raster_tj)
gg_city<-as.data.frame(gg_city)
gg_city[,3]<-gg_city[,3]*100
city<-ggplot(gg_city,aes(x=x,y=y,fill=gg_city[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="城镇百分比(%)",colours = c("#e0ecf4", "#9ebcda","#8856a7"),limits=c(0,20),breaks = seq(0,20,5))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=14),
    legend.text = element_text(size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
city
fig_name = "city_value.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 300 )
print(city)
dev.off()
###暴露性
arable_value<-0.5+0.5*(arable_value-min(arable_value,na.rm=T))/(max(arable_value,na.rm=T)-min(arable_value,na.rm=T))
tree_value<-0.5+0.5*(-tree_value+max(tree_value,na.rm=T))/(max(tree_value,na.rm=T)-min(tree_value,na.rm=T))
grass_value<-0.5+0.5*(grass_value-min(grass_value,na.rm=T))/(max(grass_value,na.rm=T)-min(grass_value,na.rm=T))
city_value<-0.5+0.5*(city_value-min(city_value,na.rm=T))/(max(city_value,na.rm=T)-min(city_value,na.rm=T))
####################干旱风险################
##########致灾因子危险性###########
##和积法求最大特征向量-致灾因子危险性
A<-read.table("land.txt")
A<-as.matrix(A)
B<-apply(A,2,FUN = sum)
C<-A
for( i in 1:3)
  
{
  C[,i]<-A[,i]/B[i]
}
D<-apply(C, 1, FUN=sum)
D<-D/sum(D)##权重系数``
E<-A%*%D
lamida<-mean(E/D)
CR<-(lamida-4)/3
drought_prone<-D[1]*P_m+D[2]*Frequency+D[3]*Severity
prone_raster<-setValues(S_SPEI,drought_prone)
gg_prone<-rasterToPoints(prone_raster)
gg_prone<-as.data.frame(gg_prone)
prone<-ggplot(gg_prone,aes(x=x,y=y,fill=gg_prone[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="危险性",colours = c("#4575b4", "#ffffbf","#d73027"),limits=c(0.5,0.9),breaks = seq(0.6, 0.9,0.1))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=14),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
prone
fig_name = "危险性.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 600 )
print(prone)
dev.off()
##########承载体暴露性##########
A<-read.table("MATRIX.txt")
A<-as.matrix(A)
B<-apply(A,2,FUN = sum)
C<-A
for( i in 1:4)
  
{
  C[,i]<-A[,i]/B[i]
}
D<-apply(C, 1, FUN=sum)
D<-D/sum(D)##权重系数``
E<-A%*%D
lamida<-mean(E/D)
CR<-(lamida-4)/3
extend<-D[1]*arable_value+D[2]*tree_value+D[3]*grass_value+D[4]*city_value
extend_raster<-setValues(S_SPEI,extend)
gg_explore<-rasterToPoints(extend_raster)
gg_explore<-as.data.frame(gg_explore)
explore<-ggplot(gg_explore,aes(x=x,y=y,fill=gg_explore[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="暴露性",colours = c("#4575b4", "#ffffbf","#d73027"),limits=c(0.5,0.91),breaks = seq(0.5, 0.9,0.1))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=14),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
explore
fig_name = "暴露性.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 600 )
print(explore)
dev.off()
###########总风险
A<-read.table("zong.txt")
A<-as.matrix(A)
B<-apply(A,2,FUN = sum)
C<-A
for( i in 1:3)
  
{
  C[,i]<-A[,i]/B[i]
}
D<-apply(C, 1, FUN=sum)
D<-D/sum(D)##权重系数``
E<-A%*%D
lamida<-mean(E/D)
CR<-(lamida-4)/3
risk<-D[1]*drought_prone+D[2]*extend+D[3]*water_value
risk_raster<-setValues(S_SPEI,risk)
gg_E<-rasterToPoints(risk_raster)
gg_E<-as.data.frame(gg_E)
ee1<-ggplot(gg_E,aes(x=x,y=y,fill=gg_E[,3]))+
  geom_raster()+
  #  geom_text(data = NULL, aes(108.8, 33.1), 
  #           label=paste("(b)"), size = 10, colour = "black")+
  # labs(x="E(°)",
  #     y="N(°)")+
  scale_x_continuous(labels=rtz)+
  scale_y_continuous(labels=ori)+
  scale_fill_gradientn(n="干旱风险",colours = c("#4575b4", "#ffffbf","#d73027"),limits=c(0.6,0.9),breaks = seq(0.6, 0.9,0.1))+
  geom_polygon(data=cc_raster,aes(x = long, y = lat, group = group), 
               color = "black", fill=NA,size=0.5)+
  theme_bw()+
  theme(
    #  axis.title.x =element_text(size=18,face = "bold"),
    #  axis.title.y=element_text(size=18,face = "bold"),
    axis.title.x =element_blank(),
    axis.title.y=element_blank(),
    axis.text= element_text( size = 16, lineheight=.1, colour="black", face = "bold"),
    plot.margin = unit(c(0.5, 0.8, 0.3, 0.5), "cm" ),
    legend.title = element_text(size=18),
    legend.text = element_text(size=14),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) 
ee1
fig_name = "draw_zong.jpg"
png(fig_name, width = 8, height = 4, units = "in", res = 600 )
print(ee1)
dev.off()
#############
a<-read.table("xian.txt")
xian_risk<-extract(risk_raster,a)
which(is.na(xian_risk))
a[58,2]<-a[58,2]+0.5
a[76,1]<-a[76,1]-0.5
xian_risk<-extract(risk_raster,a)
xian_risk<-(xian_risk-min(xian_risk))/(max(xian_risk)-min(xian_risk))
write.csv(xian_risk,file="risk.csv")
######xian calculate
a<-read.table("xianfengxian.txt")
b<-read.table("xname.txt")
c<-apply(b,1,function(x)
  {
  d<-which(a[,1]==x)
  if(length(d)==0){s<-NA}
  else{s<-a[d,2]}
  return(s)
})
c[1]<-0.94058
which(is.na(c))
c[25]<-0.709539848
c[45]<-0.575905361
c[44]<-0.755749093
c[88]<-0.78
c[89]<-0.81
v<-cbind(c,b)
write.csv(v,file="xianrisk.csv")
#################输出
risk_x<-lapply(1:103,function(x)
  {
  # x=3
  dd<-raster::extract(risk_raster,hb_xian[x,],df=T)
  dd<-mean(dd[,2],na.rm=T)
   cdd<-mask(crop(risk_raster,hb_xian[x,]),hb_xian[x,])
   cdd<-cellStats(cdd,mean,na.rm=T)
  c<-c(dd,cdd)
  return(c)
})
risk_xx<-as.data.frame(do.call(rbind,risk_x))
names(risk_xx)<-c("risk","risk_1")
hb_xian@data<-cbind(hb_xian@data,risk_xx$risk)
names(hb_xian)<-c(names(hb_xian)[-(length(names(hb_xian)))],"risk")
writeOGR(hb_xian, dsn="E:/旱情项目/Export_Output - 副本.shp",layer="risk", driver="ESRI Shapefile")
names<-cbind(hb_xian@data["NAME"],risk_xx[,1])%>%as.data.frame()
quhua<-rep(NA,times=103)
quhua[which(names[,2]>0.666810&names[,2]<0.710030)]<-"低风险区"
quhua[which(names[,2]>0.710030&names[,2]<0.748684)]<-"中低风险区"
quhua[which(names[,2]>0.748684&names[,2]<0.776147)]<-"中风险区"
quhua[which(names[,2]>0.776147&names[,2]<0.805429)]<-"中高风险区"
quhua[which(names[,2]>0.805429&names[,2]<0.857305)]<-"高风险区"
names<-cbind(names,quhua)
write.csv(names,"quhua.csv")
name_xian<-read.table("xianming.txt")[["V1"]]%>%as.character()
ss<-lapply(1:117,function(x)
  {
  # x=1
  a<-quhua[which(names[,1]==name_xian[x])]
  if(length(a)<1){a<-NA}
  return(a)
}) %>%unlist()
write.csv(ss,"quhua.csv")
##############
shuiku<-read.csv("测站基本属性表(ST_STBPRP_B)-水库水文站.csv")
hb<-shuiku[which(shuiku$交换管理单位=="湖北水文"),]
shuiku<-read.table("shuiku.txt",header=T)
aa<-shuiku[which(shuiku[,2]=="大型"),]
aa<-aa[order(-aa[,6]),]
aa<-aa[,-(2)]
write.csv(aa,"shuiku.csv")
c<-stack("index/hb_VSWI.nc")
d<-as.vector(rasterToPoints(c)[,-(1:2)])%>%na.omit()
quantile<-rank(d)/length(d)
ss<-d[which.min(abs(quantile-0.34))]
ss10<-d[which.min(abs(quantile-0.10))]
ss05<-d[which.min(abs(quantile-0.05))]
ss02<-d[which.min(abs(quantile-0.02))]
