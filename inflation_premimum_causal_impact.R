
setwd("/home/dbadmin/Documents/Upwork/Treasurary_Rate")

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

data5<-read.csv("EXTAUS.csv")
names(data5)<-c("DATE",'Taiwan_US_Exchange_Rate')
data5$DATE<-as.Date(data5$DATE,"%Y-%m-%d")


final_data<-merge(merge(merge(merge(data1,data2,by.x="DATE"),data3,by.x = "DATE"),data4,by.x="DATE"),data5,by.x="DATE")
baseline<-data.frame(apply(final_data[-c(1)],2,function(x) as.numeric(as.character(x))))

#####################################################################################

#           Visualization - Taiwan-US dollar exchage rate

####################################################################################


jpeg(filename = "Taiwan-US dollar exchange rate.jpeg",
     width = 480, height = 480)

plot(final_data$Taiwan_US_Exchange_Rate,type='l',xaxt='n')
axis(1,at=1:nrow(final_data),labels=final_data$DATE)
dev.off()

#########################################################################################

#			                  Cross-Correlation Checking

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

#########################################################################################

#                 Partial Correlation Checking based on lag

########################################################################################

Variable = 'Taiwan_US_Exchange_Rate'
library(ppcor)
Report_Partial_Correlation<-data.frame()
for(var in setdiff(names(baseline),Variable)){
for(lag in 1:5){
  Dependent = baseline[(lag + 1):nrow(baseline),Variable]
  Independent = baseline[(1):(nrow(baseline)-lag),var]
  Others = baseline[(1):(nrow(baseline)-lag),setdiff(setdiff(names(baseline),Variable),var)]
  partial_cor<-pcor.test(Dependent,Independent,Others)
  Test_Partial_Correlation<-cbind(lag,var,partial_cor[["p.value"]])
  names(Test_Partial_Correlation)<-c('lag','Dependent','p-value')
  Report_Partial_Correlation<-rbind(Report_Partial_Correlation,Test_Partial_Correlation)
  }
  
 }
names(Test_Partial_Correlation)<-c('lag','Dependent','p-value')
 
########################################################################################

#               Kalman Filtering State  Space Decomposition

#######################################################################################

 set.seed(20140504)
# length of the series
len <- nrow(Others)
treatment <- 121  #beginning of the treatment####******************
postperiod <- len - treatment
# parameters of the state
statedim <- 1 + 2
T <- diag(statedim)  #parameter matrix of state
Q <- 0.01 * diag(statedim)  #variance of state
sigmasq <- 0.01

KF <- function(y, A, H, R, Q, x10, V10) {
    len <- length(y)
    statedim <- length(x10)
    # initialization of the variance V_{t|t-1}
    Vttm1 <- V10
    # collect the state vectors x_{t|t-1} and x_{t|t}
    xttm1 <- matrix(nrow = len + 1, ncol = statedim)
    xtt <- xttm1[-1, , drop = F]
    xttm1[1, ] <- x10
    # start KF recursions
    for (i in 1:len) {
        # Var(y_k|y_{1:k-1})
        Vy <- H[i, , drop = F] %*% Vttm1 %*% t(H[i, , drop = F]) + R
        # forecast error (out of sample)
        epsilon <- y[i] - H[i, ] %*% xttm1[i, ]
        # update x_{t|t} and x_{t+1|t}
        xtt[i, ] <- xttm1[i, ] + Vttm1 %*% t(H[i, , drop = F]) %*% solve(Vy) %*% 
            epsilon
        xttm1[i + 1, ] <- xtt[i, ]
        # update V_{t|t}
        Vtt <- Vttm1 - Vttm1 %*% t(H[i, , drop = F]) %*% solve(Vy) %*% H[i, 
            , drop = F] %*% Vttm1
        # update V_{t+1|t}
        Vttm1 <- A %*% Vtt %*% t(A) + Q
    }
    return(xtt)
}

H<- as.matrix((cbind(rep(1,nrow(Others)),Score_data)))
V10 <- 10000 * diag(statedim)
x10 <- rep(0, statedim)



zsim <- KF(Dependent[1:treatment], A = T, H = H[1:treatment, , drop = F], R = sigmasq, 
    Q, x10, V10)


#######################################################################################

#                        State Representation

######################################################################################

states<- data.frame(zsim)
names(states)<-c('Trend','Seasonal_fluctuation','Irregular_fluctuation')
combined<-data.frame(cbind(Dependent[1:treatment],states))
names(combined) <- c('Response',names(states))


jpeg(filename = "State_Generation.jpeg",
     width = 720, height = 720)
plot(scale(states$Trend), type = "l", col = "black", ylim=c(-2,2),xaxt="n",ann=FALSE)
axis(1,at=1:length(final_data$DATE[1:treatment]),labels=final_data$DATE[1:treatment])
title(main="State_Generation")
lines(scale(states$Seasonal_fluctuation), col = "red")
lines(scale(states$Irregular_fluctuation), col = "green")
dev.off()

jpeg(filename = "State_Evolution.jpeg",
     width = 720, height = 720)

# plot the result of the filter (burn the first 20 observations)
burn <- 20
ztrue <- Dependent[(burn + 1):treatment]
zsim0 <- zsim[(burn + 1):treatment, ]
matplot(ztrue, type = "l", col = "black", ylim = c(min(ztrue, zsim0), max(ztrue, 
                                                                          zsim0)), xlab = "time", ylab = "", main = "State_Evolution",xaxt="n",ldy=2)
axis(1,at=1:length(final_data$DATE[1:treatment ]),labels=final_data$DATE[1:treatment])

matlines(zsim0, col = "red")

dev.off()


#######################################################################################

#                           Premium Calulation

#######################################################################################


Update_non_causal<-data.frame()
Update_causal<-data.frame()
Predicted_Premimum<-data.frame()
for(f in 1:18)
{

lm_loading<-as.matrix(round(summary(lmm<-lm(Response~Seasonal_fluctuation + Irregular_fluctuation,combined))$coefficients,3))
beta<-as.matrix(lm_loading[,"Estimate"])
  
  
state_level_var<-(vars:::VAR(states,2))
predict_state<-predict(state_level_var,n.ahead=1)$fcst

state_list<-data.frame(cbind(predict_state$Trend[1,1],predict_state$Seasonal_fluctuation[1,1],predict_state$Irregular_fluctuation[1,1]))
names(state_list)<-names(states)

non_causal_list<-data.frame(cbind(predict_state$Seasonal_fluctuation[1,1],predict_state$Irregular_fluctuation[1,1]))
names(non_causal_list)<-c('seasonal_fluctuation','Irregular_fluctuation')

states<- rbind(states,state_list)

Update_non_causal<-rbind(Update_non_causal,non_causal_list)

combined_var<-vars:::VAR(combined,2)
l=2
A_coef<-eval(parse(text = paste('cbind(',(paste(paste('vars:::Acoef(combined_var)[[',1:l,']]',sep=''),collapse = ',')),')',sep='')))
forecast<-predict(combined_var,n.ahead=1)$fcst

lists<-data.frame(cbind(forecast$Seasonal_fluctuation[1,1],forecast$Irregular_fluctuation[1,1]))
names(lists)<-c('Seasonal_fluctuation','Irregular_fluctuation')

combine_update<-data.frame(cbind(forecast$Response[1,1],forecast$Trend[1,1],forecast$Seasonal_fluctuation[1,1],forecast$Irregular_fluctuation[1,1]))
names(combine_update)<-names(combined)

Premimum_list<- lists - non_causal_list
Premium_predict<-  as.matrix(Premimum_list) %*% (beta[-1,1])

Predicted_Premimum<-rbind(Predicted_Premimum,Premium_predict)

combined<-rbind(combined,combine_update)

Update_causal<-rbind(Update_causal,lists)

}




    
pps <- function(i, zn, Z, T, sigmasq, Q, postperiod) {
    print(paste("prediction: ", i, sep = ""))
    statedim <- length(zn)
    zsim <- matrix(nrow = postperiod, ncol = statedim)
    # simulate the state
    for (i in 1:postperiod) {
        zsim[i, ] <- zn %*% T + mvrnorm(1, rep(0, statedim), Q)
        zn <- zsim[i, ]
    }
    # simulate y in the future
    ysim <- vector(length = postperiod)
    for (i in 1:postperiod) {
        ysim[i] <- Z[treatment + i, ] %*% zsim[i, ] + sqrt(sigmasq) * rnorm(1)
    }
    return(list(y = ysim, z = zsim))
}

zn <- zsim[treatment,]
sim_series <- lapply(1:1000,pps,zn,Z=H,T=T,sigmasq=sigmasq,Q=Q,postperiod=20)


y_sims <- do.call(rbind, lapply(sim_series, function(x) x$y))
yh <- yl <- ysim <- Dependent
ysim[(treatment + 1):len] <- apply(y_sims, 2, mean)
yh[(treatment + 1):len] <- apply(y_sims, 2, quantile, 0.975)
yl[(treatment + 1):len] <- apply(y_sims, 2, quantile, 0.025)
ye<-Dependent
plot(Dependent, type = "l", col = "gray", ylim = c(min(y_sims, ye, yl), max(y_sims, ye, yh)), 
    xlab = "time", ylab = "", main = "Measurement, true and predicted counterfactual")
lines(ysim, col = "red")
lines(yh, col = "red", lty = 2)
lines(yl, col = "red", lty = 2)
lines(ye, col = "blue", lwd = 2)
abline(v = treatment)
legend("topleft", c("observed", "true (unobserved) counterfactual", "predicted counterfactual (mean)", 
    "95% credible interval"), lty = c(rep(1, 3), 2, 2), lwd = c(2, rep(1, 3)), 
    col = c("blue", "gray", rep("red", 2)))

    
impact <- ye[(treatment + 1):len] - y_sims
cumimpact <- t(apply(impact, 1, cumsum))
cuml <- apply(cumimpact, 2, quantile, 0.025)
cuml[length(cuml)]


cumimp <- function(j, e, nsim) {
    print(paste("---simulation: ", j, sep = ""))
    # simulate data as before
    len <- 500
    treatment <- ceiling(2 * len/3)  #beginning of the treatment
    postperiod <- len - treatment
    # parameters of the state
    statedim <- 3
    T <- diag(statedim)  #parameter matrix of state
    Q <- 0.01 * diag(statedim)  #variance of state
    x1 <- sin(2 * pi * seq_len(len)/90)
    x2 <- sin(2 * pi * seq_len(len)/360)
    # simulation of the state
    z0 <- c(10, 0, 0)  #initial conditions: mu_0=10, beta_10=0, beta_20=0
    z <- matrix(ncol = statedim, nrow = len)
    z[1, ] <- z0 %*% T + mvrnorm(1, rep(0, statedim), Q)
    for (i in 2:len) {
        z[i, ] <- z[(i - 1), ] %*% T + mvrnorm(1, rep(0, statedim), Q)
    }
    # observations
    sigmasq <- 0.01
    Z <- matrix(ncol = statedim, nrow = len)
    Z[, 1] <- rep(1, len)
    Z[, 2] <- x1
    Z[, 3] <- x2

    y <- vector(length = len)
    for (i in 1:len) {
        y[i] <- Z[i, ] %*% z[i, ] + sqrt(sigmasq) * rnorm(1)
    }
    ye <- y
    ye[treatment:len] <- y[treatment:len] * (1 + e)

    # KF with diffuse priors
    V10 <- 10000 * diag(statedim)
    x10 <- rep(0, statedim)
    zsim <- KF(y = y[1:treatment], A = T, H = Z[1:treatment, , drop = F], R = sigmasq, 
        Q, x10, V10)

    # posterior predictive simulation
    zn <- zsim[treatment, ]
    lot_sim_series <- lapply(1:nsim, pps, zn, Z = Z, T = T, sigmasq = sigmasq, 
        Q = Q, postperiod = postperiod)
    y_sims <- do.call(rbind, lapply(lot_sim_series, function(x) x$y))

    # calculate cumulative impact
    impact <- ye[(treatment + 1):len] - y_sims
    cumimpact <- t(apply(impact, 1, cumsum))
    cuml <- apply(cumimpact, 2, quantile, 0.025)

    return(cuml[length(cuml)])
}


effect.sizes <- c(0.05, 0.1, 0.25, 0.5, 1)
sensitivity <- vector(length=length(effect.sizes)) 
for (i in 1:length(effect.sizes)) {
    impacts <- lapply(1:256,cumimp,e=effect.sizes[i],nsim=200)
    has.effect <- sapply(impacts,function(x) x>0)
    sensitivity[i] <- mean(has.effect)
}


