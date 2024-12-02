###################################
####Load required packages#########
#################################
library(data.table);library(lubridate);library(dplyr);library(dlnm);library(splines)
library(tsModel);library(forestplot);library(reshape2);library(mvmeta);
library(metafor);library(ggthemes);library(ggsci);library(tidyr)

###########################
#######preparation stage####
############################
set.seed(20240107)
###########identify sunrise and sunset###########
ephemeris <- function(lat, lon, date,tz) {
  
  # convert to the format we need
  lon.lat <- matrix(c(lon,lat), nrow=1)
  
  # make our sequence - using noon gets us around daylight saving time issues
  day <- as.POSIXct(as.character(date),tz=tz)
  
  # get our data
  sunrise <- sunriset(lon.lat, day, direction="sunrise", POSIXct.out=TRUE)
  sunset <- sunriset(lon.lat, day, direction="sunset", POSIXct.out=TRUE)
  #solar_noon <- solarnoon(lon.lat, sequence, POSIXct.out=TRUE)
  
  # build a data frame from the vectors
  data.frame(
    sunrise=format(round_date(as.POSIXct(sunrise$time,tz=tz),unit="hour"),format="%H"),
    #solarnoon=as.numeric(format(solar_noon$time, "%H%M")),
    sunset=format(round_date(as.POSIXct(sunset$time,tz=tz),unit="hour"),format="%H")
    #day_length=as.numeric(sunset$time-sunrise$time)
  )
}
####identify sunrise and sunset for each city################
dlist1<-readRDS("D:/.../ERA5_28cities_HW.rds")
loca<-read.csv("D:/.../city_lat_lon.csv")
cityname1 <- gsub("Asia/", "", cityname)
loca$city <- factor(loca$city, levels = cityname1)
loca <- arrange(loca, city)

dlist2<-dlist1
for(z in seq(dlist1)){
  sunrise<-sunset<-vector()
  xx<-dlist1[[z]]
  for (j in 1:length(xx$date1)) {
    if(z %in% 1:6){ 
      sunrise[z]<-ephemeris(loca$lat[z],loca$lon[z],xx$date[j],tz="Asia/Seoul")[1]
      sunset[z]<-ephemeris(loca$lat[z],loca$lon[z],xx$date[j],tz="Asia/Seoul")[2]
    } else 
      if(z %in% 7:21){ 
        sunrise[z]<-ephemeris(loca$lat[z],loca$lon[z],xx$date[j],tz="Asia/Shanghai")[1]
        sunset[z]<-ephemeris(loca$lat[z],loca$lon[z],xx$date[j],tz="Asia/Shanghai")[2]
      } else 
        if(z %in% 22:28){ 
          sunrise[z]<-ephemeris(loca$lat[z],loca$lon[z],xx$date[j],tz="Asia/Tokyo")[1]
          sunset[z]<-ephemeris(loca$lat[z],loca$lon[z],xx$date[j],tz="Asia/Tokyo")[2]} 
  }
  sunrise1<-unlist(sunrise)
  sunrise1<-as.numeric(as.character(sunrise1))
  
  sunset1<-unlist(sunset)
  sunset1<-as.numeric(as.character(sunset1))
  
  xx$sunrise<-sunrise1
  xx$sunset<-sunset1
  
  dlist2[[z]]<-xx
  print(z)
}

###########Match outcome with exposure##########
dlist3<-readRDS("D:/.../GRL_28cities.rds")

result <- rbindlist(lapply(seq_along(dlist3), function(i) {
  xx <- dlist3[[i]][month %in% 6:8]
  xx1 <- dlist2[[i]][, .SD, .SDcols = c(2, 32:40)]
  xx <- merge(xx, xx1, by.x = "date1", by.y = "time", all.x = TRUE)
  return(xx)
}), fill = TRUE)

dlist41 <- split(result, by = "citycode")
names(dlist41) <- names(dlist41)
####################

#################################
###########statistic analysis###########
#################################
#Define the lagged parameters
logknots  <- logknots(0:6,df=3)
arglag=list(fun="ns",knots=logknots)

#Define the cross basis
cross_basis <- function(sub, variable, arglag) {
  matx <- quantile(sub[[variable]][sub[[variable]] > 0], 0.99, na.rm = TRUE)
  sub[[variable]][sub[[variable]] > matx] <- NA
  crossbasis(sub[[variable]], lag = c(0, 6),
             argvar = list(fun = "ns",
                           knots = quantile(sub[[variable]][sub[[variable]] > 0], c(0.5, 0.9), na.rm = TRUE)),
             arglag = arglag)
}

#Define the AF calculation funcition
burden_calculate <- function(indices, variable, hwindex, dlist, blup, total, matsim, arraysim, burden) {
  k <- 0
  # Loop through the specified city indices
  for (i in indices) {
    k <- k + 1
    sub <- dlist[[i]]  # Extract data for the current city
    
    # Calculate the 99th percentile and set values above it to NA
    matx <- quantile(sub[[variable]][sub[[variable]] > 0], 0.99, na.rm = TRUE)
    sub[[variable]][sub[[variable]] > matx] <- NA
    
    # Build the crossbasis object
    HW_basis <- crossbasis(
      sub[[variable]], lag = c(0, 6),
      argvar = list(
        fun = "ns",
        knots = quantile(sub[[variable]][sub[[variable]] > 0], c(0.5, 0.9), na.rm = TRUE)),
      arglag = arglag)
    
    # Calculate the AF using attrdl function
    matsim[k, 1] <- attrdl(sub[[variable]], HW_basis, sub$all_tot, coef = blup[[k]]$blup,
                           vcov = blup[[k]]$vcov, type = "an", dir = "forw", cen = 0)
    arraysim[k, ] <- attrdl(sub[[variable]], HW_basis, sub$all_tot, coef = blup[[k]]$blup,
                            vcov = blup[[k]]$vcov, type = "an", dir = "forw", cen = 0, sim = TRUE, nsim = 1000)
  }
  
  # Summarize results
  burden[hwindex, 1] <- paste0(
    round(colSums(matsim) / sum(total) * 100, 3),
    " (", round(quantile(apply(arraysim, 2, sum), 0.025) / sum(total) * 100, 3),
    ", ", round(quantile(apply(arraysim, 2, sum), 0.975) / sum(total) * 100, 3), ")")
}

variable_names <- c("HDW", "HNW", "COM")

#meta infos
infos<-read.csv("D:/.../metacities_infos.csv")
info_pred <- data.frame(meantemp=mean(infos1$meantemp),
                        gdp_city=mean(infos1$gdp_city),
                        kgclzone=names(which(table(infos1$kgclzone) == max(table(infos1$kgclzone)))),
                        Population=mean(infos1$Population),
                        Trange=mean(infos1$Trange))

##HW distribution in all cities
datax<-bind_rows(dlist41)
datax[datax$HDW>=quantile(datax[datax$HDW>0,]$HDW,0.99,na.rm=T),]$HDW<-NA
predvar1 <-seq(min(datax$HDW,na.rm=T),max(datax$HDW,na.rm=T),by=0.1)
xvar1<-datax$HDW

datax[datax$HNW>=quantile(datax[datax$HNW>0,]$HNW,0.99,na.rm=T),]$HNW<-NA
predvar2 <-seq(min(datax$HNW,na.rm=T),max(datax$HNW,na.rm=T),by=0.1)
xvar2<-datax$HNW

datax[datax$COM>=quantile(datax[datax$COM>0,]$COM,0.99,na.rm=T),]$COM<-NA
predvar3 <-seq(min(datax$COM,na.rm=T),max(datax$COM,na.rm=T),by=0.1)
xvar3<-datax$COM

#################
##############
#Non-accidental death
#East Asia
#############
xlag  <- 0:600/10
est1<-est2<-est3<-matrix(NA,12,1)
coef<-yall<-matrix(NA,28,3)
vcov<-Sall<-vector("list",28)
total<-vector()
burden<-matrix(NA,12,1)
for (t in 1:3) {
  j<-0;
  HW.basis1 <- HW.basis2 <- HW.basis3 <- NULL
  
  for(i in c(1:28)){
    j<-j+1
    sub <- dlist41[[i]]
    HW.basis1 <- cross_basis(sub, "HDW", arglag)
    HW.basis2 <- cross_basis(sub, "HNW", arglag)
    HW.basis3 <- cross_basis(sub, "COM", arglag)
    
    model <- glm(all_tot ~ HW.basis1 + HW.basis2 + HW.basis3 +
                   ns(day,4) + ns(meanhumi,3) +
                   as.factor(year) + as.factor(dow),
                 family = quasipoisson, data = sub, na.action = "na.exclude")
    
    eval(parse(text = paste0(" pred.hw<- crossreduce(HW.basis",t,", model,
                             cen=0,lag=6)")))
    
    coef[j,] <- coef(pred.hw)
    vcov[[j]] <- vcov(pred.hw)
    
    eval(parse(text = paste0("crall<- crossreduce(HW.basis",t,",model,type='var',cen=0,value=quantile(sub[sub$HDW>0,]$HDW,0.9))")))
    yall[j,]  <- coef(crall)
    Sall[[j]] <- vcov(crall)
    
    total[j]<-sum(sub$all_tot)
    
  }
  ##ER curves
  mv<-mvmeta(coef~meantemp+Trange+gdp_city+
               kgclzone+Population,vcov,
             data=infos,method="reml")
  
  mvpred <- predict(mv,info_pred,vcov=T,format="list")
  
  
  eval(parse(text = paste0("argvar<-list(x=predvar",t,",fun='ns',
               knots=quantile(xvar",t,"[xvar",t,">0],c(50,90)/100,na.rm=T))")))
  eval(parse(text = paste0("cb.tm<-crossbasis(predvar",t,"maxlag=6,argvar=argvar,arglag=arglag)")))
  eval(parse(text = paste0("bvar<-do.call('onebasis',c(list(x=predvar",t,"),attr(cb.tm,'argvar')))")))
  eval(parse(text = paste0("all.cp<-crosspred(bvar,coef=mvpred$fit,vcov=mvpred$vcov,model.link='log',at=(predvar",t,"))")))
  
  MMT<-round(as.numeric(names(which.min(all.cp$allfit))),1)
  print(MMT)
  
  q90<-quantile(datax[datax$HDW>0,]$HDW,0.90,na.rm=T)
  q95<-quantile(datax[datax$HDW>0,]$HDW,0.95,na.rm=T)
  q99<-quantile(datax[datax$HDW>0,]$HDW,0.99,na.rm=T)
  eval(parse(text = paste0("all.cp<-crosspred(bvar,coef=mvpred$fit,vcov=mvpred$vcov,model.link='log',at=(predvar",t,",q90,q95,q99),cen=0)")))
  
  plot(all.cp, "overall", col = "#AD002AE5", cex.axis = 2.0, 
       lwd = 5, xlab="CEHWI (\u00B0C)", ylab = "              Relative risk", main = "",
       xlim = c(0, 7.2), ylim = c(0.58, 1.8),xaxt = "n", yaxt = "n",
       cex.main = 2.0, cex.lab = 2.5, ci.arg = list(col = gray(0.85)))
  axis(2, at = seq(1.0, 1.8, by = 0.2), cex.axis = 2.0)
  
  #RR
  est1[t,1]<-paste0(round(with(all.cp,cbind(allRRfit)[as.character(q90),]),3), " (", round(with(all.cp,cbind(allRRlow)[as.character(q90),]),3), ", ",round(with(all.cp,cbind(allRRhigh)[as.character(q90),]),3), ")")
  est1[t,2]<-paste0(round(with(all.cp,cbind(allRRfit)[as.character(q95),]),3), " (", round(with(all.cp,cbind(allRRlow)[as.character(q95),]),3), ", ",round(with(all.cp,cbind(allRRhigh)[as.character(q95),]),3), ")")
  est1[t,3]<-paste0(round(with(all.cp,cbind(allRRfit)[as.character(q99),]),3), " (", round(with(all.cp,cbind(allRRlow)[as.character(q99),]),3), ", ",round(with(all.cp,cbind(allRRhigh)[as.character(q99),]),3), ")")
  
  ##lag-pattern
  mv1<-mvmeta(yall~meantemp+Trange+gdp_city+
                kgclzone+Population,Sall,
              data=infos,method="reml")
  mv1pred <- predict(mv1,info_pred,vcov=T,format="list")
  
  eval(parse(text = paste0("blag  <- do.call('onebasis',c(list(x=xlag),attr(HW.basis",t,",'arglag')))")))
  cplag <- crosspred(blag,coef=mv1pred$fit,vcov=mv1pred$vcov,
                     model.link="log",at=0:600/100)
  
  plot(cplag,xlim=c(0,6),lwd=3,yaxt="n",ylim=c(0.98,1.15),ylab="Relative Risk",xlab="Lag (day)",cex.main=2.0,cex.lab=2.5,cex.axis=2,ci.arg=list(col=gray(0.85)))
  axis(2,at=c(1.00,1.05,1.10,1.15),cex.axis=2.0)
  mtext("a Non-accidental death", side = 3, line = 2.0, adj = -0.25, font = 2, cex = 1.6)
  lines(cplag,col="#AD002AE5",lwd=5)
  
  ####AF#######
  blup <-blup.mvmeta(mv,vcov=T)#blup
  
  matsim <- matrix(NA,28,1)
  nsim <- 1000
  arraysim <- array(NA,dim=c(28,1000))
  
  burden_calculate(1:28, variable_names[t], t, dlist41, blup, total, matsim, arraysim, burden)
  
}
