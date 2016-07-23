setwd("/home/prm/Documents/Upwork/Treasurary_Rate")

###########################################################################

#                         Data Preparation Part

###########################################################################
data1<-read.csv("CPIAUCSL.csv")
names(data1)<-c('DATE','CPI')
data1$DATE<-as.Date(data1$DATE,"%Y-%m-%d")

data2<-read.csv("DEXUSEU.csv")
names(data2)<-c('DATE','Exchange_Rate')
data2$DATE<-as.Date(data2$DATE,"%Y-%m-%d")

data3<-read.csv("DGS10.csv")
names(data3)<-c('DATE','Treasury_Rate')
data3$DATE<-as.Date(data3$DATE,"%Y-%m-%d")

data4<-read.csv("UNRATE.csv")
names(data4)<-c('DATE','Unemployment_Rate')
data4$DATE<-as.Date(data4$DATE,"%Y-%m-%d")

final_data<-merge(merge(merge(data1,data2,by.x="DATE"),data3,by.x = "DATE"),data4,by.x="DATE")

final_data$Exchange_Rate<-as.numeric(as.character(final_data$Exchange_Rate))

baseline<-data.frame(apply(final_data[-c(1)],2,function(x) as.numeric(as.character(x))))


#########################################################################################

#			Cross-Correlation Checking

#########################################################################################

lagind=3
ordr=5
cross_report_final<-data.frame()
all<-data.frame()
for(var in names(baseline)){
  for(lag in 1:lagind){
    print(lag)
    ####################################################################
    
    #   Data Preparation for lag correlation
    # Assuming immediate previous dates are last observations
    # Base Variable observation count: lagvalue + 1 to number of observations
    # lagged variable observations count: (1st observation to total number of observation - lag)
    
    ####################################################################
    others<-baseline[1:(nrow(baseline)-lag),!names(baseline) %in% var]
    names(others)<-paste(names(others),lag,sep='_')
    dep<-baseline[(lag+1):nrow(baseline),var]
    names(dep)<-var
    data<-data.frame(cbind(dep,others))
    names(data)[1]<-paste(var)
    
    #####################################################################
    
    # Correlation MAtrix Calculation and Cross-report preparation
    
    ####################################################################
    corMatrix<-data.frame(cor(data))
    ordercor<-apply(corMatrix,2,function(x) order(abs(x),decreasing = TRUE))
    cross_report<-matrix(nrow = (nrow(ordercor)-1),ncol=ncol(ordercor))
    for( i in 2:nrow(ordercor)){
      for(j in 1:ncol(ordercor)){
        cross_report[i-1,j] = paste(names(data[as.numeric(ordercor[i,j])]))
      }
    }
    #assign(paste("Cross_Correlation_report",lag,var,sep="_"),cross_report)
    #assign(paste("Cross_Correlation.",lag,var,sep="_"),corMatrix)
    mat<-corMatrix
    names(corMatrix)<-c('id',paste('V',seq(1,(ncol(data)-1),sep='_')))
    corMatrix<-data.frame(cbind('lag'=lag,'BaseVar'=var,'Attribute'=names(mat),corMatrix))
    all<-rbind(all,corMatrix)
  }
  subs<-all[all[,which(names(all) == 'id')] <1,]
  ord<-order(abs(subs[,which(names(all) == 'id')]),decreasing = T)
  name<-as.character(subs$Attribute[ord[1:ordr]])
  value<-subs[,which(names(all) == 'id')][ord[1:ordr]]
  row<-as.character(paste(name,'(',value,')',sep=""))
  cross_report_final<-rbind(cross_report_final,cbind(paste(var),t(row)))
  
}
names(cross_report_final)<-c('Attribute',paste('Dependency',seq(1:ordr),sep='_'))

#write.csv(cross_report_final,"cross_report_final.csv",row.names=F)

########################################################################################

#			 Stationarity Checking

########################################################################################



library(tseries)
report<-data.frame()
baseline<-data.frame(baseline)
for(idx in 1:ncol(baseline)){
  
  z= baseline[,idx]  
  
  for(i in 1:10){
    #print(i)
    j=0
    if(adf.test(z)$p.value>0.1){
      
      repeat{
        
        z<-diff(z)
        print("Non-St")
        
        lag=j+1
        j=j+1
        print(j)
        if(adf.test(z)$p.value<0.05){
          report<-rbind(report,j)
          break
        }
      }
    }
  }
  
  
}
report<-data.frame(cbind(as.matrix(names(baseline)),report))
names(report)<-c('Variable','Difference_for_stationarity')

#write.csv(report,"lag_report.csv",row.names=F)

orders<-max(report$Difference_for_stationarity)
stationary_data<-data.frame(apply(baseline,2,function(x) arima(x,c(0,orders,0))$residual))

########################################################################################

#		Extraction of Trend(Long-Term Movement)

########################################################################################

Trend_Part<- baseline - stationary_data


library(lubridate)
at<-Trend_Part$DATE

jpeg(filename = "Trend.jpeg",
     width = 480, height = 480)
plot(as.matrix(scale(Trend_Part$CPI)),type='l',ylim=c(-2,2))
lines(scale(Trend_Part$Exchange_Rate),col='red')
lines(scale(Trend_Part$Treasury_Rate),col='blue')
lines(scale(Trend_Part$Unemployment_Rate),col='green')
dev.off()

jpeg(filename = "Error.jpeg",
     width = 480, height = 480)
plot(scale(stationary_data$CPI),type='l',ylim=c(-2,2))
lines(scale(stationary_data$Exchange_Rate),col='red')
lines(scale(stationary_data$Treasury_Rate),col='blue')
lines(scale(stationary_data$Unemployment_Rate),col='green')
dev.off()



library(MASS)

var_trend<-vars:::VAR(Trend_Part, p=2)

#######################################################################################################

#		Testing Granger Causality in Long-Term Movement

#######################################################################################################

Causality_Report<-data.frame()
for(var in names(baseline)){
  causality<-vars:::causality(var_trend,cause=setdiff(names(baseline),var))
  Cause = paste(setdiff(names(baseline),var),collapse = '&')
  Effect = var
  Gr_P_value = round(as.numeric(causality$Granger$p.value),3)
  Instant_P_value = round(as.numeric(causality$Instant$p.value),3)
  rep<-data.frame(cbind(Cause,Effect,Gr_P_value,Instant_P_value))
  Causality_Report<- rbind(Causality_Report,rep)
}
write.csv(Causality_Report,"Causality_Report.csv",row.names = F)

#######################################################################################################

#			Checking correlation for short term dpendency

#######################################################################################################



cor_matrix<-cor(baseline)
write.csv(cor_matrix,"cor_matrix.csv",row.names=F)

#######################################################################################################

#			Extracting Pattern

#######################################################################################################

coeff<-data.frame((var_trend$coefficients))
residual<-data.frame(var_trend$residuals)
serial_long_term<-vars:::serial.test(var_trend)$serial$p.value
arch_long_term<-vars:::arch.test(var_trend)$arch.mul$p.value
fevd<-data.frame(vars:::fevd(var_trend)$Treasury_Rate)

jpeg(filename = "FEVD.jpeg",
     width = 480, height = 480)
plot(scale(fevd[,1]),type='l',ylim=c(-2,2))
lines(scale(fevd[,2]),col='red')
lines(scale(fevd[,3]),col='blue')
lines(scale(fevd[,4]),col='green')
dev.off()

irf<-data.frame(vars:::irf(var_trend)$irf$Treasury_Rate)

jpeg(filename = "IRF_long_term.jpeg",
     width = 480, height = 480)
plot(scale(irf[,1]),type='l',ylim=c(-2,2))
lines(scale(irf[,2]),col='red')
lines(scale(irf[,3]),col='blue')
lines(scale(irf[,4]),col='green')
dev.off()

#################################################################################

#		Co-Integration Checking

#################################################################################

library(urca)
c<-ca.jo(baseline)
critical_value<-data.frame(cbind(as.matrix(c@teststat),as.matrix(c@cval)[,1]))
names(critical_value)<-c('Teststat','CV')

write.csv(critical_value,"critical_value.csv",row.names=F)

##################################################################################

#		Exploring short term fluctuation part

#################################################################################

var_err<-vars:::VAR(stationary_data,orders)
summary(var_err)

#################################################################################

#		Extracting Pattern

#################################################################################

serial_short_term<-vars:::serial.test(var_err)$serial$p.value
arch_short_term<-vars:::arch.test(var_err)$arch.mul$p.value
fevd_err<-data.frame(vars:::fevd(var_err)$Treasury_Rate)
jpeg(filename = "FEVD.jpeg",
     width = 480, height = 480)
plot(scale(fevd_err[,1]),type='l',ylim=c(-2,2))
lines(scale(fevd_err[,2]),col='red')
lines(scale(fevd_err[,3]),col='blue')
lines(scale(fevd_err[,4]),col='green')
dev.off()



irf_err<- vars:::irf(var_err)$irf$Treasury_Rate

jpeg(filename = "IRF_short_term.jpeg",
     width = 480, height = 480)
plot(scale(irf_err[,1]),type='l',ylim=c(-2,2))
lines(scale(irf_err[,2]),col='red')
lines(scale(irf_err[,3]),col='blue')
lines(scale(irf_err[,4]),col='green')
dev.off()


################################################################################

#		Rolling Prediction

################################################################################


trend_Predict<-matrix(nrow=18,ncol=ncol(baseline))
trend_CI<-matrix(nrow=18,ncol=ncol(baseline))
Fluctuation_Predict<-matrix(nrow=18,ncol=ncol(baseline))
Fluctuation_CI<-matrix(nrow=18,ncol=ncol(baseline))
colnames(trend_Predict)=colnames(trend_CI)=colnames(Fluctuation_Predict)=names(baseline)
for( f in seq(1,18,2))
{
  var_trend<-vars:::VAR(Trend_Part, p=2)
  prediction<-predict(var_trend,n.ahead=2)
  
  trend_Predict[f:(f+(2-1)),1]<-prediction$fcst$CPI[,1]
  trend_Predict[f:(f+(2-1)),2]<-prediction$fcst$Exchange_Rate[,1]
  trend_Predict[f:(f+(2-1)),3]<-prediction$fcst$Treasury_Rate[,1]
  trend_Predict[f:(f+(2-1)),4]<-prediction$fcst$Unemployment_Rate[,1]
  
  trend_CI[f:(f+(2-1)),1]<-prediction$fcst$CPI[,4]
  trend_CI[f:(f+(2-1)),2]<-prediction$fcst$Exchange_Rate[,4]
  trend_CI[f:(f+(2-1)),3]<-prediction$fcst$Treasury_Rate[,4]
  trend_CI[f:(f+(2-1)),4]<-prediction$fcst$Unemployment_Rate[,4]
  
  Trend_Part<-rbind(Trend_Part,trend_Predict[(f:(f+(2-1))),])
  
  var_err<-vars:::VAR(stationary_data,2)
  prediction1<-predict(var_err,n.ahead=2)
  
  Fluctuation_Predict[f:(f+(2-1)),1]<-prediction1$fcst$CPI[,1]
  Fluctuation_Predict[f:(f+(2-1)),2]<-prediction1$fcst$Exchange_Rate[,1]
  Fluctuation_Predict[f:(f+(2-1)),3]<-prediction1$fcst$Treasury_Rate[,1]
  Fluctuation_Predict[f:(f+(2-1)),4]<-prediction1$fcst$Unemployment_Rate[,1]
  
  Fluctuation_CI[f:(f+(2-1)),1]<-prediction1$fcst$CPI[,4]
  Fluctuation_CI[f:(f+(2-1)),2]<-prediction1$fcst$Exchange_Rate[,4]
  Fluctuation_CI[f:(f+(2-1)),3]<-prediction1$fcst$Treasury_Rate[,4]
  Fluctuation_CI[f:(f+(2-1)),4]<-prediction1$fcst$Unemployment_Rate[,4]
  
  stationary_data<-rbind(stationary_data,Fluctuation_Predict[(f:(f+(2-1))),])
}


FINAL_Prediction<-trend_Predict + Fluctuation_Predict

Date=as.character(seq(as.Date("2016/7/1"), as.Date("2017/12/1"), "months"),"%Y/%m/%d")

FINAL_Prediction<-data.frame(cbind(Date,(FINAL_Prediction)))


jpeg(filename = "Prediction.jpeg",
     width = 480, height = 480)

plot(as.numeric(as.character(FINAL_Prediction$Treasury_Rate)),type='l',xlab='Date',ylab='%_value')
axis(1,at=1:length(FINAL_Prediction[["Treasury_Rate"]]),labels=FINAL_Prediction[["Date"]])
title(main="Treasury Rate Prediction",sub="July 2016 to December 2017")
dev.off()
