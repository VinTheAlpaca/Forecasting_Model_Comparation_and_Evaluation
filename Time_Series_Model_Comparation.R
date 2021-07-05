############################################################################
############################ Regression Model ##############################
############################################################################

list.of.packages <- c("Mcomp", 
                      "ggplot2", 
                      "forecast", 
                      "fpp", 
                      "tseries", 
                      "lmtest", 
                      "ggpubr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Mcomp)
library(ggplot2)
library(forecast)
library(fpp)
library(tseries)
library(lmtest)
library(ggpubr)


horizon = 18
y = M3[[1907]]$x
M3[[1907]]$description

plot(decompose(y))

fit.lin1 = tslm(y ~ trend + season) #tslm frequency is disregarded
summary(fit.lin1)
checkresiduals(fit.lin1)

t <- time(M3[[1907]]$x)
knot1 = 1988 + 1/12         # January of 1988 is the first knot
knot2 = 1991 + 1/12         # January of 1991 is the second knot

dumVar1 = ts(12*pmax(0,t - knot1),start = c(1982,1))
dumVar2 = ts(12*pmax(0,t - knot2),start = c(1982,1))

fit.lin2 = tslm(y ~ trend + season + t + dumVar1)
fit.lin2 = tslm(y ~ trend + season + t + dumVar1 + dumVar2)
fit.lin2 = tslm(y ~ trend + season + t + dumVar1 + dumVar2 + dumVar3)
fit.lin2 = tslm(y ~ trend + season + t + I(dumVar1^(2)) + I(dumVar2^2) + dumVar3)

fit.lin.final = tslm(y ~ trend + season + t + dumVar1 + dumVar2)
summary(fit.lin.final)
checkresiduals(fit.lin.final)

t.new = t[length(t)] + seq(horizon)/12 

dumVar1.new = dumVar1[length(dumVar1)] + seq(horizon)
dumVar2.new = dumVar1[length(dumVar1)] + seq(horizon)

newdata = as.data.frame(cbind(t = t.new, dumVar1 = dumVar1.new, dumVar2 = dumVar2.new))

autoplot(y) + 
  autolayer(fitted(fit.lin.final),series="Regression") + 
  autolayer(forecast(fit.lin.final, newdata = newdata),series="Regression (Forecast)") + 
  autolayer(M3[[1907]]$xx, series="Real Data", ) +
  xlab("Year") +
  ylab("Shipments") + 
  scale_color_manual(values = c("green","red", "blue")) +
  #ggtitle("Shipments of Cement in Portland") +
  theme(legend.position = "top", legend.title = element_blank())


### Residuals Diagonsis
shapiro.test(fit.lin.final$residuals) #50.84% > 5%, it does follow normality
Box.test(fit.lin.final$residuals, lag = 12, type = "Box-Pierce") # 1.181e-06 < 5% indicates a part of residuals can be distinguished from white noise
Box.test(fit.lin.final$residuals, lag = 12, type = "Ljung-Box") # 3.667e-07 < 5% indicates a part of residuals can be distinguished from white noise

 
install.packages("lmtest")
library(lmtest)

dwtest(fit.lin.final) # p-Value = 0.76 positive autocorrelation 

bptest(fit.lin.final) # p-value = 41.13% > 5% it does not show heteroscedasticity
bgtest(fit.lin.final) # p-value = 4.402e-14 small value indicates a significant autocorrelation in residuals


### Error Measures
accuracy(forecast(fit.lin.final, newdata = newdata), M3[[1907]]$xx)
accuracy(snaive(y), M3[[1907]]$xx)

############################################################################
######################## Exponential Smoothing Model #######################
############################################################################
decompose(y)

plot(y, xlim=c(1982,1995))

# Simple exponential smoothing
fit.exp1 = ets(y, model="ANN", damped=FALSE)
summary(fit.exp1)

# Holt exponential smoothing
fit.exp2 = ets(y, model="AAN", damped=FALSE)
summary(fit.exp2)

# Seasonal exponential smoothing
fit.exp3 = ets(y, model="ANA", damped=FALSE)
summary(fit.exp3)

# Damped exponential smoothing
fit.exp4 = ets(y, model="AAN", damped=TRUE)
summary(fit.exp4)

# Holt-winter exponential smoothing
fit.exp5 = ets(y, model="AAA", damped=FALSE)
summary(fit.exp5)

# Damped Holt-winter exponential smoothing
fit.exp6 = ets(y, model="AAA", damped=TRUE)
summary(fit.exp6)
plot(stl(y, t.window = 20, s.window = 20))

forecast(fit.exp6)
accuracy(snaive(y))

autoplot(y) + 
  autolayer(fitted(fit.exp6),series="EXP") + 
  autolayer(forecast(fit.exp6, h=horizon),series="forecast2") + 
  autolayer(M3[[1907]]$xx, series="Real Data", ) +
  xlab("Year") +
  ylab("Shipments") + 
  scale_color_manual(values = c("red","blue", "green")) +
  theme(legend.position = "top", legend.title = element_blank())


checkresiduals(fit.exp6)
mean(fit.exp6$residuals)

### Residuals Diagonsis
shapiro.test(fit.exp6$residuals) #61.69% > 5%, it does follow normality
Box.test(fit.exp6$residuals, lag = 12, type = "Box-Pierce") # 8.12e-05 < 5% indicates a part of residuals can be distinguished from white noise
Box.test(fit.exp6$residuals, lag = 12, type = "Ljung-Box") # 2.906e-05 < 5% indicates a part of residuals can be distinguished from white noise

#dwtest(fit.exp6) # p-Value = 0.76 positive autocorrelation 

#bptest(fit.exp6) # p-value = 41.13% > 5% it does not show heteroscedasticity

### Error Measure

accuracy(forecast(fit.exp6), M3[[1907]]$xx)
accuracy(forecast(fit.exp5), M3[[1907]]$xx)
############################################################################
############################## ARIMA Model #################################
############################################################################

adf.test(y, k=12)
kpss.test(y, lshort = FALSE)

kpss.test(diff(y), lshort = FALSE)

ndiffs(y)
nsdiffs(y)

ggtsdisplay(y, lag=1)
ggtsdisplay(diff(y, lag=1))

tsdisplay(diff(diff(y, lag=12))) # ARIMA(4,1,0)(0,1,0)12 or ARIMA(0,1,5)(0,1,1)12

fit.ar1 = Arima(y, order = c(4,1,0), seasonal = c(0,1,0))
fit.ma2 = Arima(y, order = c(0,1,5), seasonal = c(0,1,1))


EMs.ar = array(NA,c(12, 6))
counter.ar = 1
for (p in 4:1) {
  for (P in 0:2){
    fit.ar = try(Arima(y, order = c(p,1,0), seasonal = c(p,1,0)), silent = TRUE)
    
    if (typeof(fit.ar) == "list"){
      accu.ar = accuracy(fit.ar)
      
      EMs.ar[counter.ar, 1] = p
      EMs.ar[counter.ar, 2] = P
      EMs.ar[counter.ar, 3] = accu.ar[4]  # MPE
      EMs.ar[counter.ar, 4] = accu.ar[5]  # MAPE
      EMs.ar[counter.ar, 5] = accu.ar[6]  # MASE
      EMs.ar[counter.ar, 6] = fit.ar$aicc # AICc
    }else{
      EMs.ar[counter.ar, 1] = p
      EMs.ar[counter.ar, 2] = P
    }
    counter.ar = counter.ar + 1
  }
}

EMs.ma = array(NA,c(15, 6))
counter.ma = 1
for (q in 5:1) {
  for (Q in 1:3){
    fit.ma = try(Arima(y, order = c(0,1,q), seasonal = c(0,1,Q)), silent = TRUE)
    
    if (typeof(fit.ar) == "list"){
      accu.ma = accuracy(fit.ma)
      
      EMs.ma[counter.ma, 1] = q
      EMs.ma[counter.ma, 2] = Q
      EMs.ma[counter.ma, 3] = accu.ma[4]  # MPE
      EMs.ma[counter.ma, 4] = accu.ma[5]  # MAPE
      EMs.ma[counter.ma, 5] = accu.ma[6]  # MASE
      EMs.ma[counter.ma, 6] = fit.ma$aicc # AICc
    }else{
      EMs.ma[counter.ma, 1] = q
      EMs.ma[counter.ma, 2] = Q
    }
    counter.ma = counter.ma + 1
  }
}


fit.ar11 = Arima(y, order = c(4,1,0), seasonal = c(0,1,0))
fit.ar12 = Arima(y, order = c(2,1,0), seasonal = c(0,1,0))
fit.ar13 = Arima(y, order = c(1,1,0), seasonal = c(0,1,0))

fit.ma21 = Arima(y, order = c(0,1,5), seasonal = c(0,1,2))
fit.ma22 = Arima(y, order = c(0,1,4), seasonal = c(0,1,2))
fit.ma23 = Arima(y, order = c(0,1,3), seasonal = c(0,1,2))
fit.ma24 = Arima(y, order = c(0,1,2), seasonal = c(0,1,2))
fit.ma25 = Arima(y, order = c(0,1,1), seasonal = c(0,1,2))

summary(fit.ar11)
summary(fit.ar12)
summary(fit.ar13)
summary(fit.ma21)
summary(fit.ma22)
summary(fit.ma23)
summary(fit.ma24)
summary(fit.ma25)

checkresiduals(fit.ma21, lag=12)


autoplot(y) + 
  autolayer(fitted(fit.ma21),series="ARIMA") + 
  autolayer(forecast(fit.ma21, h=18),series="ARIMA(Forecast)") + 
  autolayer(M3[[1907]]$xx, series="Real Data") +
  xlab("Year") +
  ylab("Shipments") + 
  scale_color_manual(values = c("red","blue", "green")) +
  theme(legend.position = "top", legend.title = element_blank())

accuracy(forecast(fit.arima21), M3[[1907]]$xx)


mean(fit.arima21$residuals)
var(fit.arima21$residuals)

### Residuals Diagonsis 
Box.test(fit.ma21$residuals, lag = 12, type = "Box-Pierce") # p-value= 0.6664 > 0.005 Autocorrelation is fine
Box.test(fit.ma21$residuals, lag = 12, type = "Ljung-Box") # p-value= 0.5973 > 0.005 Autocorrelation is fine

### Normality Diagnosis 
hist(fit.ar21$residuals)
plot(fit.ar21$residuals)
qqnorm(fit.ar21$residuals)
qqline(fit.ar21$residuals, col="red") # shows it is not quite normally distributed

shapiro.test(fit.arima21$residuals) #0.6% < 5%, it does not follow normality

checkresiduals(fit.arima21)









############################part 2#######################################
############################part 2#######################################
############################part 2#######################################
############################part 2#######################################
############################part 2#######################################





library(Mcomp)
library(MAPA)
library(forecast)
library(Metrics)
library(gridExtra)
library(grid)

lastNumOfID = 6
firstNum = 1500 + lastNumOfID
seqOfTimeSeries = seq(firstNum, 2800, by = 10)
horizon = 18


#############################################################################
################## Automatic Exponential Smoothing Model ####################
#############################################################################
start_time1 = Sys.time()

insampleMPEs.ets = array(0, c(length(seqOfTimeSeries),2))

fcs.ets.list = array(NA, c(length(seqOfTimeSeries), horizon+1))

counter.ets = 1
for (ts_no in seqOfTimeSeries){
  
  cat(sprintf("ETS: %s\n", counter.ets))
  
  insampleMPEs.ets[counter.ets,1] = ts_no
  
  fit.autoets = ets(M3[[ts_no]]$x)
  
  fcs.ets.list[counter.ets, 1] = ts_no
  fcs.ets.list[counter.ets, 2:(horizon+1)] = forecast(fit.autoets,horizon)$mean
  
  insampleMPEs.ets[counter.ets,2] = mean(100*(M3[[ts_no]]$x - fit.autoets$fitted)/M3[[ts_no]]$x)   ### 0104
  
  counter.ets = counter.ets + 1
}
end_time1 = Sys.time()
end_time1 - start_time1

############################################################################
########################### Automatic ARIMA Model ##########################
############################################################################
start_time2 = Sys.time()

insampleMPEs.arima = array(0, c(length(seqOfTimeSeries),2))

fcs.arima.list = array(NA, c(length(seqOfTimeSeries), horizon+1))

counter.arima = 1
for (ts_no in seqOfTimeSeries){
  
  cat(sprintf("ARIMA: %s\n", counter.arima))
  
  insampleMPEs.arima[counter.arima,1] = ts_no
  
  fit.autoarima = auto.arima(M3[[ts_no]]$x)
  
  fcs.arima.list[counter.arima, 1] = ts_no
  fcs.arima.list[counter.arima, 2:(horizon+1)] = forecast(fit.autoarima,horizon)$mean
  #fcs.arima.list[counter.arima, (horizon+2)] = fit.autoarima$aicc
  
  insampleMPEs.arima[counter.arima,2] = mean(100*(M3[[ts_no]]$x - fit.autoarima$fitted)/M3[[ts_no]]$x)
  
  counter.arima = counter.arima + 1
}
end_time2 = Sys.time()
end_time2 - start_time2
############################################################################
########################### Automatic MAPA Model ###########################
############################################################################
start_time3 = Sys.time()

insampleMPEs.mapa = array(0, c(length(seqOfTimeSeries),2))

fcs.mapa.list = array(NA, c(length(seqOfTimeSeries), horizon+1))

counter.mapa = 1
for (ts_no in seqOfTimeSeries){
  
  cat(sprintf("MAPA: %s\n", counter.mapa))
  
  insampleMPEs.mapa[counter.mapa,1] = ts_no
  
  fit.automapa = mapa(M3[[ts_no]]$x, fh = horizon)
  
  fcs.mapa.list[counter.mapa, 1] = ts_no
  fcs.mapa.list[counter.mapa, 2:(horizon+1)] = fit.automapa$outfor
  
  insampleMAPATS = ts(c(fit.automapa$infor), start=time(M3[[ts_no]]$x)[1], frequency=12)
  
  insampleMPEs.mapa[counter.mapa,2] = mean(100*(M3[[ts_no]]$x - insampleMAPATS)/M3[[ts_no]]$x, na.rm = TRUE)
  
  counter.mapa = counter.mapa + 1
}
end_time3 = Sys.time()
end_time3 - start_time3
############################################################################
########################## Self-defined Single Stretagy ####################
############################################################################

no_method = 3
seqOfTimeSeries = seq(firstNum, 2800, by = 10)
best_model = array(NA, c(length(seqOfTimeSeries) , 2))

############################### Cross-Validation ###########################

start_time4 = Sys.time()

counter.sfd = 1
for (tsi in seqOfTimeSeries) {
  
  cat(sprintf("SFD CROSS-VALIDATION: %s\n", counter.sfd))
  
  curTS = M3[[tsi]]$x
  curTsTime = time(curTS)
  origins1 = (length(curTS) - (horizon*2) + 1) : (length(curTS) - horizon)
  origins = c(origins1[6],origins1[12], origins1[18])
  FCs = array(0, c(no_method, length(origins), horizon))
  actuals = array(0, c(length(origins), horizon))
  #fcs.sfd.list = array(NA, c(length(seqOfTimeSeries), horizon+3))
  
  for (origin in origins){
    insample = ts(curTS[1:origin], start=curTsTime[1], frequency=12)
    actuals[which(origin==origins), ] = curTS[(origin+1):(origin + horizon)]
    
    
    # calculate the fit and store the forecast of ets, arima and mapa
    fit.sfdets = ets(insample)
    FCs[1, which(origin==origins), ] = forecast(fit.sfdets, h=horizon)$mean
    fit.sfdarima = auto.arima(insample)
    FCs[2, which(origin==origins), ] = forecast(fit.sfdarima, h=horizon)$mean
    fit.sfdmapa = mapa(insample, fh=horizon)
    FCs[3, which(origin==origins), ] = fit.sfdmapa$outfor
  }
  
  # define a array for storing MASE
  MAEs = array(0, no_method)
  
  # calculate the MASE for all models
  for (m in 1:no_method){
    # calculate the MAE for model m and save MASE
    MAEs[m] = mean(abs(actuals - FCs[m,,]))                  # mase(actuals, FCs[m,,], step_size = 12)
  }
  
  #Identify which method has the lowest MAPE and update the frequencies
  #freq_best[counter.sfd,which.min(MAPEs)] = freq_best[which.min(MAPEs)] + 1
  best_model[counter.sfd,1] = tsi
  #best_model[counter.sfd,2] = which.min(MAPEs)
  best_model[counter.sfd,2] = which.min(MAEs)
  counter.sfd = counter.sfd + 1
}
end_time4 = Sys.time()
end_time4 - start_time4



########################## forecast out-sample data #######################
start_time5 = Sys.time()

insampleMPEs.sfd = array(0, c(length(seqOfTimeSeries),2))

fcs.sfd.list = array(NA, c(length(seqOfTimeSeries), horizon+1))

counter.sfd2 = 1
for (tsi in 1:length(seqOfTimeSeries)){
  
  cat(sprintf("SFD FORECAST: %s\n", tsi))
  
  series_no = best_model[tsi,1]
  model_chosen = best_model[tsi,2]
  current_ts = M3[[series_no]]$x
  
  insampleMPEs.sfd[tsi,1] = series_no
  
  if (model_chosen == 1){                               #ets
    fit.sfd = ets(current_ts)
    fcast = forecast(fit.sfd, h = horizon)$mean
    
    insampleMPEs.sfd[tsi,2] = mean(100*(current_ts - fit.sfd$fitted)/current_ts)
    
  } 
  
  else if (model_chosen == 2){                          #arima
    fit.sfd = auto.arima(current_ts)
    fcast = forecast(fit.sfd, h = horizon)$mean
    
    insampleMPEs.sfd[tsi,2] = mean(100*(current_ts - fit.sfd$fitted)/current_ts)
    
  }else {                                               #mapa
    mapaout = mapa(current_ts, fh = horizon)
    
    fcast = mapaout$outfor
    insampleMAPATS = ts(c(mapaout$infor), start=time(current_ts)[1], frequency=12)
    
    insampleMPEs.sfd[tsi,2] = mean(100*(current_ts - insampleMAPATS)/current_ts, na.rm = TRUE)
    
  }
  fcs.sfd.list[tsi, 1] = ts_no
  fcs.sfd.list[tsi, 2:(horizon+1)] = fcast
  
}
end_time5 = Sys.time()
end_time5 - start_time5

############################################################################
########################## Combination Strategy ############################ 
############################################################################



########################## Combination Strategy 1 ##########################

fcs.comb.list = array(NA, c(length(seqOfTimeSeries), horizon + 1))

models = array(1/4, 4)
for (tsi in 1:length(seqOfTimeSeries)){
  
  series_no = seqOfTimeSeries[tsi]
  
  fcs.comb.list[tsi,1] = series_no
  
  ets_MPEs = insampleMPEs.ets[tsi,2] 
  arima_MPEs = insampleMPEs.arima[tsi,2]
  mapa_MPEs = insampleMPEs.mapa[tsi,2]
  sfd_MPEs = insampleMPEs.sfd[tsi,2]
  
  
  MPEs.all = array(c(c(1,2,3,4), c(ets_MPEs, arima_MPEs, mapa_MPEs, sfd_MPEs)), c(4,2))
  
  MPEs.po = array(NA, c(4,2))
  MPEs.ne = array(NA, c(4,2))
  
  counter.po = 1
  counter.ne = 1
  for (mpe in 1:4){
    if (MPEs.all[mpe,2] >= 0 ){
      MPEs.po[counter.po,1] = MPEs.all[mpe,1]
      MPEs.po[counter.po,2] = MPEs.all[mpe,2]
      counter.po = counter.po + 1
    }else{
      MPEs.ne[counter.ne,1] = MPEs.all[mpe,1]
      MPEs.ne[counter.ne,2] = MPEs.all[mpe,2]
      counter.ne = counter.ne + 1
    }
  }
  
  if (is.na(MPEs.po[1,1]) | is.na(MPEs.ne[1,1])){   ### there is no balance, equal weights for all four forecasts
    models[1] = 1/4
    models[2] = 1/4
    models[3] = 1/4
    models[4] = 1/4
    
  }else{                                            ### there is a balance, combine four forcasts
    
    if (counter.po == 2){                           ### one positive
      mean.ne = abs((MPEs.ne[1,2] + MPEs.ne[2,2] + MPEs.ne[3,2]) / 3)
      
      models[MPEs.po[1,1]] = mean.ne / (MPEs.po[1,2] + mean.ne)
      
      models[MPEs.ne[1,1]] = (1 - models[MPEs.po[1,1]])/3
      models[MPEs.ne[2,1]] = (1 - models[MPEs.po[1,1]])/3
      models[MPEs.ne[3,1]] = (1 - models[MPEs.po[1,1]])/3
    } 
    
    if (counter.ne == 2){                           ### one negative
      mean.po = abs((MPEs.po[1,2] + MPEs.po[2,2] + MPEs.po[3,2]) / 3)
      
      models[MPEs.ne[1,1]] = mean.po / (abs(MPEs.ne[1,2]) + mean.po)
      
      models[MPEs.po[1,1]] = (1 - models[MPEs.ne[1,1]])/3
      models[MPEs.po[2,1]] = (1 - models[MPEs.ne[1,1]])/3
      models[MPEs.po[3,1]] = (1 - models[MPEs.ne[1,1]])/3
    }
    
    if (counter.po > 2 & counter.ne > 2){
      MPEs.sum = sum(abs(MPEs.all[,2]))
      
      models[MPEs.po[which.min(MPEs.po[,2]),1]] = abs(MPEs.ne[which.min(MPEs.ne[,2]),2]) / MPEs.sum
      models[MPEs.po[which.max(MPEs.po[,2]),1]] = abs(MPEs.ne[which.max(MPEs.ne[,2]),2]) / MPEs.sum
      models[MPEs.ne[which.min(MPEs.ne[,2]),1]] = MPEs.po[which.min(MPEs.po[,2]),2] / MPEs.sum
      models[MPEs.ne[which.max(MPEs.ne[,2]),1]] = MPEs.po[which.max(MPEs.po[,2]),2] / MPEs.sum
    }
    
  }                                
  weighted_ets = fcs.ets.list[tsi,2:(horizon+1)]*models[1]
  weighted_arima = fcs.arima.list[tsi,2:(horizon+1)]*models[2]
  weighted_mapa = fcs.mapa.list[tsi,2:(horizon+1)]*models[3]
  weighted_sfd = fcs.sfd.list[tsi,2:(horizon+1)]*models[4] 
  
  fcs.comb.list[tsi,2:(horizon+1)] = weighted_ets + weighted_arima + weighted_mapa + weighted_sfd
  
}


########################## Combination Strategy 2 ##########################

models.ew = array(1/4, 4)

#equal weights combination
fcs.comb.list.ew = array(NA, c(length(seqOfTimeSeries), horizon + 1))

for (tsi in 1:length(seqOfTimeSeries)){
  
  series_no = seqOfTimeSeries[tsi]
  
  fcs.comb.list.ew[tsi,1] = series_no
  
  weight_ets.ew = fcs.ets.list[tsi,2:(horizon+1)]*models.ew[1]
  weight_arima.ew = fcs.arima.list[tsi,2:(horizon+1)]*models.ew[2]
  weight_mapa.ew = fcs.mapa.list[tsi,2:(horizon+1)]*models.ew[3]
  weight_sfd.ew = fcs.sfd.list[tsi,2:(horizon+1)]*models.ew[4] 
  
  fcs.comb.list.ew[tsi,2:(horizon+1)] = weight_ets.ew + weight_arima.ew + weight_mapa.ew + weight_sfd.ew
  
}

############################################################################
############################ Benchmark Methods ############################# 
############################################################################

### Benchmark model1 naive method
start_time_b1 = Sys.time()
benmak.fcs.naive.list = array(NA, c(length(seqOfTimeSeries), horizon + 1))
for (tsi in 1:length(seqOfTimeSeries)){
  
  cat(sprintf("Benchmark methods: %s\n", tsi))
  
  series_no = seqOfTimeSeries[tsi]
  current_ts = M3[[series_no]]$x
  
  benmak.fcs.naive.list[tsi, 1] = series_no
  benmak.fcs.naive.list[tsi, 2:(horizon+1)] = naive(current_ts, h = horizon)$mean
  
}
end_time_b1 = Sys.time()
end_time_b1 - start_time_b1

### Benchmark model2 snaive method
start_time_b2 = Sys.time()
benmak.fcs.snaive.list = array(NA, c(length(seqOfTimeSeries), horizon + 1))
for (tsi in 1:length(seqOfTimeSeries)){
  series_no = seqOfTimeSeries[tsi]
  current_ts = M3[[series_no]]$x
  
  benmak.fcs.snaive.list[tsi, 1] = series_no
  benmak.fcs.snaive.list[tsi, 2:(horizon+1)] = snaive(current_ts, h = horizon)$mean
}
end_time_b2 = Sys.time()
end_time_b2 - start_time_b2

### Benchmark model3 damped ets method
start_time_b3 = Sys.time()
benmak.fcs.dampets.list = array(NA, c(length(seqOfTimeSeries), horizon + 1))
for (tsi in 1:length(seqOfTimeSeries)){
  series_no = seqOfTimeSeries[tsi]
  current_ts = M3[[series_no]]$x
  
  benmak.fcs.dampets.list[tsi, 1] = series_no
  fit.dampets = ets(current_ts, model="ZZZ", damped=TRUE)
  benmak.fcs.dampets.list[tsi, 2:(horizon+1)] = forecast(fit.dampets, h = horizon)$mean
}
end_time_b3 = Sys.time()
end_time_b3 - start_time_b3
############################################################################
################################# Analysis ################################# 
############################################################################


### Error measure function, which will return a vector that contains the serier no, MAE, MASE, and sMAPE
errorMeasures <- function(tsi, fcs.list, fh){
  series_no = seqOfTimeSeries[tsi] 
  fcastdata = fcs.list[tsi,2:(fh+1)]
  fcastTS= ts(c(fcastdata), start=time(tail(M3[[series_no]]$x, 1)) + 1/12, frequency=12)
  
  MAE = mean(abs(M3[[series_no]]$xx[1:fh] - fcastdata))
  MASE = mase(M3[[series_no]]$xx, fcastTS, step_size = 12)
  sMAPE = mean(200*(abs(M3[[series_no]]$xx[1:horizon] - fcastdata)/(M3[[series_no]]$xx[1:horizon] + fcastdata)))
  return(c(series_no, MAE, MASE, sMAPE))
}


############################## Horizon Analysis ############################ 
meanMAEs = array(NA, c(3, 9))
meanMASEs = array(NA, c(3, 9))
meansMAPEs = array(NA, c(3, 9))

fh = 12

errMeasures.ets = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.arima = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.mapa = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.sfd = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.comb = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.comb.ew = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.benmak.naive = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.benmak.snaive = array(NA, c(length(seqOfTimeSeries), 4))
errMeasures.benmak.dampets = array(NA, c(length(seqOfTimeSeries), 4))

for (tsi in 1:length(seqOfTimeSeries)){
  
  ### ETS Error Measures
  errMeasures.ets[tsi,] = errorMeasures(tsi, fcs.ets.list, fh)
  ### ARIMA Error Measures
  errMeasures.arima[tsi,] = errorMeasures(tsi, fcs.arima.list, fh)
  ### MAPA Error Measures
  errMeasures.mapa[tsi,] = errorMeasures(tsi, fcs.mapa.list, fh)
  ### SFD Error Measures
  errMeasures.sfd[tsi,] = errorMeasures(tsi, fcs.sfd.list, fh)
  ### Combination Error Measures
  errMeasures.comb[tsi,] = errorMeasures(tsi, fcs.comb.list, fh)
  ### Combination (equal weights) Error Measures
  errMeasures.comb.ew[tsi,] = errorMeasures(tsi, fcs.comb.list.ew, fh)
  ### Benchmark Method Naive Error Measures
  errMeasures.benmak.naive[tsi,] = errorMeasures(tsi, benmak.fcs.naive.list, fh)
  ### Benchmark Method sNaive Error Measures
  errMeasures.benmak.snaive[tsi,] = errorMeasures(tsi, benmak.fcs.snaive.list, fh)
  ### Benchmark Method Damped ETS Error Measures
  errMeasures.benmak.dampets[tsi,] = errorMeasures(tsi, benmak.fcs.dampets.list, fh)
  
}

MAEs.list = array(NA, c(length(seqOfTimeSeries), 10))
MAEs.list[,1] = errMeasures.ets[,1]
MAEs.list[,2] = errMeasures.ets[,2]
MAEs.list[,3] = errMeasures.arima[,2]
MAEs.list[,4] = errMeasures.mapa[,2]
MAEs.list[,5] = errMeasures.sfd[,2]
MAEs.list[,6] = errMeasures.comb.ew[,2]
MAEs.list[,7] = errMeasures.comb[,2]
MAEs.list[,8] = errMeasures.benmak.snaive[,2]
MAEs.list[,9] = errMeasures.benmak.naive[,2]
MAEs.list[,10] = errMeasures.benmak.dampets[,2]


MASEs.list = array(NA, c(length(seqOfTimeSeries), 10))
MASEs.list[,1] = errMeasures.ets[,1]
MASEs.list[,2] = errMeasures.ets[,3]
MASEs.list[,3] = errMeasures.arima[,3]
MASEs.list[,4] = errMeasures.mapa[,3]
MASEs.list[,5] = errMeasures.sfd[,3]
MASEs.list[,6] = errMeasures.comb.ew[,3]
MASEs.list[,7] = errMeasures.comb[,3]
MASEs.list[,8] = errMeasures.benmak.snaive[,3]
MASEs.list[,9] = errMeasures.benmak.naive[,3]
MASEs.list[,10] = errMeasures.benmak.dampets[,3]


sMAPEs.list = array(NA, c(length(seqOfTimeSeries), 10))
sMAPEs.list[,1] = errMeasures.ets[,1]
sMAPEs.list[,2] = errMeasures.ets[,4]
sMAPEs.list[,3] = errMeasures.arima[,4]
sMAPEs.list[,4] = errMeasures.mapa[,4]
sMAPEs.list[,5] = errMeasures.sfd[,4]
sMAPEs.list[,6] = errMeasures.comb.ew[,4]
sMAPEs.list[,7] = errMeasures.comb[,4]
sMAPEs.list[,8] = errMeasures.benmak.snaive[,4]
sMAPEs.list[,9] = errMeasures.benmak.naive[,4]
sMAPEs.list[,10] = errMeasures.benmak.dampets[,4]


for (x in 1:9){
  if (fh == 3){
    a = 1
  }else if (fh == 12){
    a = 2
  }else if (fh == 18){
    a = 3
  }
  meanMAEs[a,x] = mean(MAEs.list[,x+1])
  meanMASEs[a,x] = mean(MASEs.list[,x+1])
  meansMAPEs[a,x] = mean(sMAPEs.list[,x+1])
}

palette = rainbow(9)
plot(meanMASEs[,1], xlab = "Horizon", ylab = "MASE", xlim = c(0.4, 3.2), ylim = c(0, 2), xaxt="n", type = "n", main = " MASE Changes")
axis(1, 1:3, c("3 month", "12 Month", "18 Month"))
for (x in 1:9){
  lines(meanMASEs[,x], col = palette[x])
}

legendName = c("EXP", "Arima", "MAPA", "SFD", "Comb1", "Comb2", "sNaive", "Naive", "Damped")
legend("topleft", legend=legendName, col =palette[1:9] , lty=1, cex=0.75)


names = c("ETS", "ARIMA", "MAPA", "SFD", "Comb1","Comb2", "sNaive")

boxplot(MASEs.list[,2:8], names = names,ylim=c(0,1.5), outpch = NA,medcol = "red", main="MASE")
boxplot(sMAPEs.list[,2:8], names = names,ylim=c(0,45), outpch = NA,medcol = "red", main="sMAPE%")


############################# Type Analysis ########################

MASEs.list.micro = array(NA, c(length(seqOfTimeSeries), 9))
MASEs.list.macro = array(NA, c(length(seqOfTimeSeries), 9))
MASEs.list.industry = array(NA, c(length(seqOfTimeSeries), 9))
MASEs.list.demographic = array(NA, c(length(seqOfTimeSeries), 9))
MASEs.list.finance = array(NA, c(length(seqOfTimeSeries), 9))
MASEs.list.other = array(NA, c(length(seqOfTimeSeries), 9))

for (tsi in 1:length(seqOfTimeSeries)) {
  series_no = seqOfTimeSeries[tsi]
  current_ts = M3[[series_no]]$x
  
  if (M3[[series_no]]$type == "MICRO"){
    MASEs.list.micro[tsi,] = MASEs.list[tsi,2:10]
    
  }else if (M3[[series_no]]$type == "MACRO"){
    MASEs.list.macro[tsi,] = MASEs.list[tsi,2:10]
    
  }else if (M3[[series_no]]$type == "INDUSTRY"){
    MASEs.list.industry[tsi,] = MASEs.list[tsi,2:10]
    
  }else if (M3[[series_no]]$type == "DEMOGRAPHIC"){
    MASEs.list.demographic[tsi,] = MASEs.list[tsi,2:10]
    
  }else if (M3[[series_no]]$type == "FINANCE"){
    MASEs.list.finance[tsi,] = MASEs.list[tsi,2:10]
    
  }else{
    MASEs.list.other[tsi,] = MASEs.list[tsi,2:10]
  }
}

typeMASE.list = array(NA, c(6, 9))
for (x1 in 1:9){
  
  typeMASE.list[1,x1] = mean(MASEs.list.micro[,x1], na.rm = TRUE)
  typeMASE.list[2,x1] = mean(MASEs.list.macro[,x1], na.rm = TRUE)
  typeMASE.list[3,x1] = mean(MASEs.list.industry[,x1], na.rm = TRUE)
  typeMASE.list[4,x1] = mean(MASEs.list.demographic[,x1], na.rm = TRUE)
  typeMASE.list[5,x1] = mean(MASEs.list.finance[,x1], na.rm = TRUE)
  typeMASE.list[6,x1] = mean(MASEs.list.other[,x1], na.rm = TRUE)
}


plot(typeMASE.list[,1], xlab = "Type", ylab = "MASE", xlim = c(0, 6), ylim = c(0, 2), xaxt="n", type = "n", main = " MASE Changes")
axis(1, 1:6, c("MICRO", "MACRO", "INDUSTRY", "DEMO", "FINANCE", "OTHER"))
for (x in 1:9){
  lines(typeMASE.list[,x], col = palette1[x])
}

legendName = c("EXP", "Arima", "MAPA", "SFD", "Comb1", "Comb2", "sNaive", "Naive", "Damped")
legend("topleft", legend=legendName, col =palette[1:9] , lty=1, cex=0.75)




################## characteristics Analysis ################

plot(decompose(M3[[1516]]$x)) #complex trend up
plot(decompose(M3[[1526]]$x)) #near-liner trend up
plot(decompose(M3[[1556]]$x)) #linear trend up
plot(decompose(M3[[1566]]$x)) #linear trend up


plot(decompose(M3[[1696]]$x)) # near-linear down
plot(decompose(M3[[1686]]$x)) # near-linear down
plot(decompose(M3[[1706]]$x)) #complex trend down
plot(decompose(M3[[1716]]$x)) #complex trend down

sets = c("1516", "1526", "1556", "1566", "1696", "1686", "1706", "1716")
plot(MAEs.list[which(MAEs.list[,1] == 1516),])
plot(MAEs.list[which(MAEs.list[,1] == 1516),], xlab = "Type", ylab = "MASE", xlim = c(0, 10), ylim = c(300, 2000), xaxt="n", type = "n", main = " MASE Changes")
axis(1, 1:9, c("EXP", "ARIMA", "MAPA", "SFD", "COMB1", "COMB2", "sNaive", "Naive", "Damped"))
for (x2 in 1:length(sets)){
  
  lines(MAEs.list[which(MAEs.list[,1] == sets[x2]),2:10], col = palette1[x2])
  
}

#legendName = c("EXP", "Arima", "MAPA", "SFD", "Comb1", "Comb2", "sNaive", "Naive", "Damped")
legend("topleft", legend=sets, col =palette[1:9] , lty=1, cex=0.75)


