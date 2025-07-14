#Zwei mögliche Datasets:


###Reiseverkehrsdaten 1971-2018. Von Österreich, Schweiz und Spanien.

#path of data set removed due to data safety. 

###Produktion, Umsatz und Auftragseingang 1991-2018. Saisonal bereinigt und original Zeitreihe:

                                                                                          

#Mit dem ersten Dataset:
x=dataset$AusgabenSpanien[1:567]
y=dataset$AusgabenÖsterreich[1:567]


#Hauptfunktion:
scoint = function(x, y, option="default"){       #option  = "manual" and defining M or option = "default" can be used.
  library(smooth)                                 
  #path of distr1 and path of distr2 removed due to data safety.
  ##################################################################################
  
  Diff = function(x,lag=1, ...) {
    y <- c(rep(NA, lag) , diff(x,lag=lag, ...))
    class(y) <- "ts"
    y <- ts(y, start=stats::start(x), frequency=stats::frequency(x))
    return(y)
  }
  Lag = function(x,k) {
    y<- c(rep(NA,k), x)[1:length(x)] 
    class(y) <- "ts"
    y <- ts(y, start=stats::start(x), frequency=stats::frequency(x))
    return(y)
  }
  ftest = function(testmodel, reg1) {
    ((deviance(testmodel) - deviance(reg1)) / (df.residual(testmodel) - df.residual(reg1))) / (deviance(reg1) / df.residual(reg1))
  } 
  pval = function(testval, distribution, twosided=T, side="left") {
    if(twosided) {
      a <- rowSums(as.numeric(testval) > t(distribution))/length(distribution)
      a <- ifelse(a < 0.5, a, 1-a)*2
    } else {
      if (side=="left") {a <- rowSums(as.numeric(testval) > t(distribution))/length(distribution)}
      if (side=="right") {a <- rowSums(as.numeric(testval) < t(distribution))/length(distribution)}
    }
    
    return(a)
  }
  interpolation = function (tvalue, distrib, side , i){
    t=length(x)
    #Intervall finden.
    a = matrix(c(0,60,60,120,120,240,240,1200,1200,100000),ncol=2, byrow=TRUE)
    u = which(apply(a,1,findInterval, x=t) == 1)
    
    if (1 < u & u < 5){
      pval_u = pval(as.numeric(tvalue), distrib[[(u-1)]][,i], twosided=F, side) #Untergrenze.
      pval_o = pval(as.numeric(tvalue), distrib[[u]][,i], twosided=F, side)     #Obergrenze.
      pval_in = pval_u + (pval_o - pval_u) * (t - a[u,1]) / (a[u,2] - a[u,1])   #Interpolation
      
    } else {
      if(u == 1){
        pval_in = pval(as.numeric(tvalue), distrib[[(1)]][,i], twosided=F, side="left") #u=1 Untergrenze
      } else {
        pval_in = pval(as.numeric(tvalue), distrib[[4]][,i], twosided=F, side="left") #u=4 Obergrenze
        
      }
    }
    
    return(pval_in)
  }
  
  ##################################################################################
  
  #Function for linear filters for respective frequency (HEGY)
  Filter = function(series){ #Beaulieu & Miron, S.308
    
    y1 = series + Lag(series, 1) + Lag(series,2) + Lag(series,3) + Lag(series,4)+ Lag(series,5) + Lag(series,6) + Lag(series,7) + Lag(series,8) + Lag(series,9) + Lag(series,10) + Lag(series,11) 
    y2 = -(series - Lag(series, 1) + Lag(series,2) - Lag(series,3) + Lag(series,4) - Lag(series,5) + Lag(series,6) - Lag(series,7) + Lag(series,8) - Lag(series,9) + Lag(series,10) - Lag(series,11) )
    y3 = - (Lag(series, 1) - Lag(series,3) +  Lag(series,5) - Lag(series,7) +  Lag(series,9)  - Lag(series,11) )
    y4 = - (series - Lag(series,2) +  Lag(series,4) - Lag(series,6) +  Lag(series,8)  - Lag(series,10) )
    y5 = -0.5*(series + Lag(series, 1) - 2*Lag(series,2) + Lag(series,3) + Lag(series,4) - 2*Lag(series,5) + Lag(series,6) + Lag(series,7) - 2*Lag(series,8) + Lag(series,9) + Lag(series,10) - 2*Lag(series,11))
    y6 = sqrt(3)/2 * (series - Lag(series,1) + Lag(series,3) - Lag(series,4) + Lag(series,6) - Lag(series,7) + Lag(series,9) - Lag(series,10))
    y7 = 0.5*(series - Lag(series, 1) - 2*Lag(series,2) - Lag(series,3) + Lag(series,4) + 2*Lag(series,5) + Lag(series,6) - Lag(series,7) - 2*Lag(series,8) - Lag(series,9) + Lag(series,10) + 2*Lag(series,11))
    y8 = -sqrt(3)/2 * (series + Lag(series,1) - Lag(series,3) - Lag(series,4) + Lag(series,6) + Lag(series,7) - Lag(series,9) - Lag(series,10))
    y9 = -0.5*(sqrt(3)*series - Lag(series,1) + Lag(series,3) - sqrt(3)*Lag(series,4) + 2*Lag(series,5) - sqrt(3)*Lag(series,6) + Lag(series,7) - Lag(series,9) + sqrt(3)*Lag(series,10) - 2*Lag(series,11))
    y10 = 0.5*(series - sqrt(3)*Lag(series,1) + 2*Lag(series,2) - sqrt(3)*Lag(series,3) + Lag(series,4) - Lag(series,6) + sqrt(3)*Lag(series,7) - 2*Lag(series,8) + sqrt(3)*Lag(series,9) - Lag(series,10))
    y11 = 0.5*(sqrt(3)*series + Lag(series,1) - Lag(series,3) - sqrt(3)*Lag(series,4) - 2*Lag(series,5) - sqrt(3)*Lag(series,6) - Lag(series,7) + Lag(series,9) + sqrt(3)*Lag(series,10) + 2*Lag(series,11))
    y12 = -0.5*(series + sqrt(3)*Lag(series,1) + 2*Lag(series,2) + sqrt(3)*Lag(series,3) + Lag(series,4) - Lag(series,6) - sqrt(3)*Lag(series,7) - 2*Lag(series,8) - sqrt(3)*Lag(series,9) - Lag(series,10))
    
    y13 = (series - Lag(series,12))
    
    filter = list(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13)
    return(filter)
  }
  
  #Function for linear filters for respective frequency (EGHL)
  Filter_eghl = function(x){
    y1 = x + Lag(x,1) + Lag(x,2) + Lag(x,3) +Lag(x,4) + Lag(x,5) + Lag(x,6) + Lag(x,7) + Lag(x,8) + Lag(x,9) + Lag(x,10) + Lag(x,11)
    y2 = - (x-Lag(x,1)+Lag(x,2)-Lag(x,3)+Lag(x,4)-Lag(x,5)+Lag(x,6)-Lag(x,7)+Lag(x,8)-Lag(x,9)+Lag(x,10)-Lag(x,11))
    y3 = - (x-Lag(x,2)+Lag(x,4)-Lag(x,6)+Lag(x,8)-Lag(x,10))
    y4 = - (x-Lag(x,1)+Lag(x,3)-Lag(x,4)+Lag(x,6)-Lag(x,7)+Lag(x,9)-Lag(x,10))
    y5 = - (x+Lag(x,1)-Lag(x,3)-Lag(x,4)+Lag(x,6)+Lag(x,7)-Lag(x,9)-Lag(x,10))
    y6 = - (x +2*Lag(x,2)+Lag(x,4)-Lag(x,6)-2*Lag(x,8)-Lag(x,10)-sqrt(3)*Lag(x,1)-sqrt(3)*Lag(x,3)+sqrt(3)*Lag(x,7)+sqrt(3)*Lag(x,9))
    y7 = - (x+2*Lag(x,2)+Lag(x,4)-Lag(x,6)-2*Lag(x,8)-Lag(x,10)+sqrt(3)*Lag(x,1)+sqrt(3)*Lag(x,3)-sqrt(3)*Lag(x,7)-sqrt(3)*Lag(x,9))
    y8 = (x-Lag(x,12))
    
    filter = list(y1,y2,y3,y4,y5,y6,y7,y8)
    return(filter)
    
  }
  
  #HEGY function.
  int = function(z){
    filz = Filter(z)
    seriesz = unlist(filz[13])
    au = data.frame("Reg1"=Lag(unlist(filz[1]), 1),"Reg2"=Lag(unlist(filz[2]), 1),"Reg3"=Lag(unlist(filz[3]), 1),"Reg4"=Lag(unlist(filz[4]), 1),"Reg5"=Lag(unlist(filz[5]), 1),"Reg6"=Lag(unlist(filz[6]), 1),"Reg7"=Lag(unlist(filz[7]), 1),"Reg8"=Lag(unlist(filz[8]), 1),"Reg9"=Lag(unlist(filz[9]), 1), "Reg10"=Lag(unlist(filz[10]), 1),"Reg11"=Lag(unlist(filz[11]), 1),"Reg12"=Lag(unlist(filz[12]), 1))
    az = data.frame("Lag1" =Lag(seriesz,1),"Lag2" =Lag(seriesz,2),"Lag3" =Lag(seriesz,3),"Lag4" =Lag(seriesz,4),"Lag5" =Lag(seriesz,5),"Lag6" =Lag(seriesz,6),"Lag7" =Lag(seriesz,7),"Lag8" =Lag(seriesz,8),"Lag9" =Lag(seriesz,9),"Lag10" =Lag(seriesz,10),"Lag11" =Lag(seriesz,11), "Lag12" =Lag(seriesz,12)
                    ,"Lag13" =Lag(seriesz,13),"Lag14" =Lag(seriesz,14),"Lag15" =Lag(seriesz,15),"Lag16" =Lag(seriesz,16),"Lag17" =Lag(seriesz,17),"Lag18" =Lag(seriesz,18),"Lag19" =Lag(seriesz,19),"Lag20" =Lag(seriesz,20),"Lag21" =Lag(seriesz,21),"Lag22" =Lag(seriesz,22),"Lag23" =Lag(seriesz,23), "Lag24" =Lag(seriesz,24))
    
    name = names(az); p = 1 
    #Augmentation selection procedure.
    for (i in 1:ncol(az)){
      az = data.frame(az)
      if (ncol(az)==1){
        colnames(az)=name}
     # reg1z = lm(seriesz ~., data=data.frame(au,az[1:ncol(az)])) ## Alte version
      reg1z = lm(seriesz ~., data=data.frame(au,az)) 
      p = coef(summary(reg1z))[-(1:13),4]
      
      weg = which(p==max(p))
      if (max(p) > 0.05 & ncol(az)>1){ ## Vielleicht eher p-value=0.01 als cut-off
        az = data.frame(az[,-weg]); name = name[-weg]} else if (max(p) > 0.05 & ncol(az)==1){
          reg1z = lm(seriesz ~ as.ts(au)) 
        }
    }
    
    #partial F-test. 
    chi2 = function(model, k) {
      mat = vcov(model)[(k:(k+1)),(k:(k+1))] 
      coef = model$coefficients[(k:(k+1))]
      chi2 = t(coef) %*% solve(mat) %*% coef 
      nd = length(model$residuals)
      f_krit =  as.numeric(chi2/(nrow(mat)) * (nd  - nrow(mat))/ (nd)) 
      return(f_krit)
    }
    
    #test statistics.                             
    t1z = coef(summary(reg1z))[2,'t value']          
    t2z = coef(summary(reg1z))[3,'t value']          
    f34z = chi2(reg1z,4)
    f56z = chi2(reg1z,6)
    f78z = chi2(reg1z,8)
    f910z = chi2(reg1z,10)
    f1112z = chi2(reg1z,12)
    
    lags1_names = names(coef(summary(reg1z))[,4])
    lags1 = as.numeric(gsub("Lag", "", lags1_names[-(1:13)]))
    
    teststats1 = list(t1 = t1z, t2 = t2z, f34 = f34z, f56 = f56z, f78 = f78z, f910 = f910z, f1112 = f1112z)
    output1 = list(teststats1,lags1)                           #Output consists of test statistics and selected augmentations.
    return(output1)
  }
  
  #Simulation for "manual" option:
  aug = function(augm, regres2){
    a=NULL
    names = gsub("$", "Lag", augm)
    if (length(augm)==0){
      a=NULL
    } else{
      for (i in 1:length(augm)){
        
        a=cbind(a,Lag(regres2,augm[i]))
        
      }
      a = data.frame(a)
      names(a) = c(names)}
    return(a)
  } #generate a data.frame object with all augmentations for simulation.
  
  #Simulation of the zero distribution if taking "manual" option (EGHL):
  simulation = function(lags){   
    
    x5 = data.frame(matrix(NA, length(x), 0))  
    y5 = data.frame(matrix(NA, length(y), 0))  
    full_eghl_rslts5 = NULL
    
      EGHL_sim = function(x,y){
      
      filx = Filter_eghl(x)
      fily = Filter_eghl(y)
      
      #frequency = 0
      reg1 = lm(fily[[1]] ~ filx[[1]]) 
      res1 = resid(reg1)
      aug1 = aug(lags[[1]], Diff(res1)) #generate a data.frame object with all selected augmentations. Needed for regression. 
      if (length(aug1)==0){
        resreg1 = lm(Diff(res1) ~ 0 + Lag(res1,1))
      } else {
        resreg1 = lm(Diff(res1) ~. + 0 + Lag(res1,1), data = aug1)}     
      t1 = coef(summary(resreg1))[1,'t value']
      
     
      #frequency = 6/12
      reg2 = lm(fily[[2]] ~ filx[[2]]) 
      res2 = resid(reg2)
      v = Filter_eghl(res2); regres2 = -v[[8]]/v[[2]]         # =(1+L)
      aug2 = aug(lags[[2]], regres2)                                
      if (length(aug2)==0){
        resreg2 = lm(regres2 ~ 0 + Lag(-res2,1)) 
      } else {resreg2 = lm(regres2 ~. + 0 + Lag(-res2,1), data = aug2)} 
      t2 = coef(summary(resreg2))[1,'t value']
      
    

      
      pair_reg = function(k){
        reg3 = lm(fily[[k]] ~ filx[[k]] + Lag(filx[[k]],1))
        res3 = resid(reg3)
        w = Filter_eghl(res3)
        regres = -w[[8]]/w[[k]]
        aug3 = aug(lags[[k]], regres)        
        if (length(aug3)==0){
          resreg = lm(regres ~ 0 + Lag(-res3,2) + Lag(-res3,1))
        } else {
          resreg = lm(regres ~.  + 0 + Lag(-res3,2) + Lag(-res3,1), data = aug3)} 
        return(resreg)
      }
      
      #frequency = 3/12 = 9/12
      resreg3 = pair_reg(3)       
      fstat3 = summary(resreg3)$f[1]
      
      #Theta = 
      resreg5 = pair_reg(4) 
      fstat5 = summary(resreg5)$f[1]
      
      #Theta =
      resreg7 = pair_reg(5) 
      fstat7 = summary(resreg7)$f[1]
      
      #Theta =
      resreg9 = pair_reg(6) 
      fstat9 = summary(resreg9)$f[1]
      
      #Theta = 
      resreg11 = pair_reg(7) 
      fstat11 = summary(resreg11)$f[1]
      
      teststat1_sim = data.frame(t1 = t1, t2 = t2, f3 = fstat3, f5 = fstat5, f7 = fstat7, f9 = fstat9, f11 = fstat11)
      return(teststat1_sim)
    }
    
    
    
    for (i in 1:M){
      x <- smooth::sim.ssarima(orders = list(ar=c(0,0), i=c(0,1), ma=c(0,0)), lags=c(1,12), obs=length(x))$data
      x <- x/10^(log10(abs(mean(x)))) ;   x <- (x - min(x)) + 100
      x5 <- cbind(x5, as.numeric(x))
      y <- smooth:sim.ssarima(orders = list(ar=c(0,0), i=c(0,1), ma=c(0,0)), lags=c(1,12), obs=length(x))$data
      y <- y/10^(log10(abs(mean(y)))) ;   y <- (y - min(y)) + 100
      y5 <- cbind(y5, as.numeric(y))
      full_eghl_rslts5 = rbind(full_eghl_rslts5, EGHL_sim(ts(x5[,i]), ts(y5[,i])))
      print(i)
      
    }
    return(full_eghl_rslts5)
  } 
  #Simulation if the zero distribution if taking "manual" option (HEGY):
  simulation2 = function(lags, M=M){
    
    x5 = data.frame(matrix(NA, length(x), 0))  
    full_hegy_rslts5 = NULL
    
    int_sim = function(x){
      filz = Filter(x)
      seriesz = unlist(filz[13])
      au = data.frame("Reg1"=Lag(unlist(filz[1]), 1),"Reg2"=Lag(unlist(filz[2]), 1),"Reg3"=Lag(unlist(filz[3]), 1),"Reg4"=Lag(unlist(filz[4]), 1),"Reg5"=Lag(unlist(filz[5]), 1),"Reg6"=Lag(unlist(filz[6]), 1),"Reg7"=Lag(unlist(filz[7]), 1),"Reg8"=Lag(unlist(filz[8]), 1),"Reg9"=Lag(unlist(filz[9]), 1), "Reg10"=Lag(unlist(filz[10]), 1),"Reg11"=Lag(unlist(filz[11]), 1),"Reg12"=Lag(unlist(filz[12]), 1))
      aug2 = aug(lags,seriesz)
      
      #Full model.
      if (length(aug2)==0){
        regression = lm(seriesz ~., data = au)
      } else {
        regression = lm(seriesz ~., data = data.frame(au,aug2))
      }
      
      #Test 
      chi2 = function(model, k) {
        mat = vcov(model)[(k:(k+1)),(k:(k+1))] 
        coef = model$coefficients[(k:(k+1))]
        chi2 = t(coef) %*% solve(mat) %*% coef 
        
        nd = length(model$residuals)
        f_krit =  as.numeric(chi2/(nrow(mat)) * (nd  - nrow(mat))/ (nd)) 
        return(f_krit)
      }
      
      #test statistics.                                #Frequenzen benennen.
      t1z = coef(summary(regression))[2,'t value']          #With intercept.
      t2z = coef(summary(regression))[3,'t value']          #With intercept.
      f34z = chi2(regression,4)
      f56z = chi2(regression,6)
      f78z = chi2(regression,8)
      f910z = chi2(regression,10)
      f1112z = chi2(regression,12)
      
      
      teststat2_sim = list(t1 = t1z, t2 = t2z, f34 = f34z, f56 = f56z, f78 = f78z, f910 = f910z, f1112 = f1112z)
      
      return(teststat2_sim)
    }
    
    for (i in 1:M){
      x <- smooth::sim.ssarima(orders = list(ar=c(0,0), i=c(0,1), ma=c(0,0)), lags=c(1,12), obs=length(x))$data
      x <- x/10^(log10(abs(mean(x)))) ;   x <- (x - min(x)) + 100
      x5 <- cbind(x5, as.numeric(x))
      full_hegy_rslts5 = rbind(full_hegy_rslts5, int_sim(ts(x5[,i])))
      
    }
    return(full_hegy_rslts5)
  }  
  
  ####################################################################################  
  
  #Main function: EGHL test for cointegration
  EGHL = function(x,y){
    
    filx = Filter_eghl(x)
    fily = Filter_eghl(y)
  
    #Frequency = 0.
    reg1 = lm(fily[[1]] ~ filx[[1]]) 
    res1 = resid(reg1)
    
    #Augmentation selection procedure for frequency 0.
    lagsel1 = function(regres,reg1){                        ###reg1!?
      at = data.frame("Lag1" = Lag(regres,1) , "Lag2" =  Lag(regres,2) ,"Lag3" =  Lag(regres,3) , "Lag4" =  Lag(regres,4), "Lag5" =  Lag(regres,5) ,"Lag6" = Lag(regres,6) , "Lag7" = Lag(regres,7) ,"Lag8" = Lag(regres,8) ,"Lag9" = Lag(regres,9) ,"Lag10" = Lag(regres,10) ,"Lag11" = Lag(regres,11) , "Lag12" =Lag(regres,12) 
                      , "Lag13" = Lag(regres,13) ,  "Lag14" =Lag(regres,14) ,"Lag15" = Lag(regres,15) , "Lag16" =Lag(regres,16) , "Lag17" =Lag(regres,17) , "Lag18" =Lag(regres,18) , "Lag19" =Lag(regres,19) , "Lag20" =Lag(regres,20) , "Lag21" =Lag(regres,21) , "Lag22" =Lag(regres,22) , "Lag23" =Lag(regres,23) , "Lag24" =Lag(regres,24))
      name = names(at); p = 1
      for (i in 1:ncol(at)){ 
        at = data.frame(at)
        if (ncol(at)==1){
          colnames(at)=name}
        resreg1 = lm(regres ~. + 0 + Lag(res1,1), data = at)                          
        p = coef(summary(resreg1))[-(ncol(at)+1),4]                                  
        
        weg = which(p==max(p))
        if (max(p) > 0.05 & ncol(at)>1){
          at = at[,-weg]; name = name[-weg]} else if (max(p) > 0.05 & ncol(at)==1){
            resreg1 = lm(regres ~ 0 + Lag(res1,1))}
        
      }
      
      
      return(resreg1)
    }
    
    resreg1 = lagsel1(Diff(res1), reg1)                           
    t1 = coef(summary(resreg1))[1,'t value']
    #Extracting names.
    lags1_names = names(coef(summary(resreg1))[,4])
    lags1 = lags1_names[-(length(lags1_names))]   
    lags1 = as.numeric(gsub("Lag", "", lags1))
    
    
    
  
    #Frequency = 6/12
    reg2 = lm(fily[[2]] ~ filx[[2]]) 
    res2 = resid(reg2)
    #Augmentation selection procedure for frequency 6/12.
    lagsel2 = function(regres){
      at = data.frame("Lag1" = Lag(regres,1) , "Lag2" =  Lag(regres,2) ,"Lag3" =  Lag(regres,3) , "Lag4" =  Lag(regres,4), "Lag5" =  Lag(regres,5) ,"Lag6" = Lag(regres,6) , "Lag7" = Lag(regres,7) ,"Lag8" = Lag(regres,8) ,"Lag9" = Lag(regres,9) ,"Lag10" = Lag(regres,10) ,"Lag11" = Lag(regres,11) , "Lag12" =Lag(regres,12) 
                      , "Lag13" = Lag(regres,13) ,  "Lag14" =Lag(regres,14) ,"Lag15" = Lag(regres,15) , "Lag16" =Lag(regres,16) , "Lag17" =Lag(regres,17) , "Lag18" =Lag(regres,18) , "Lag19" =Lag(regres,19) , "Lag20" =Lag(regres,20) , "Lag21" =Lag(regres,21) , "Lag22" =Lag(regres,22) , "Lag23" =Lag(regres,23) , "Lag24" =Lag(regres,24))
      name = names(at); p = 1
      for (i in 1:ncol(at)){ #ineffizient mit 24 Wdh.
        at = data.frame(at)
        if (ncol(at)==1 & min(p)<=0.05){
          colnames(at)=name}
        resreg2 = lm(regres ~. + 0 + Lag(-res2,1), at)
        p = coef(summary(resreg2))[-(ncol(at)+1),4]
        
        weg = which(p==max(p))
        if (max(p) > 0.05 & ncol(at)>1){
          at = at[,-weg]; name = name[-weg]} else if (max(p) > 0.05 & ncol(at)==1){#Kürzen des Vectors zu einem Null vector wird übersprungen und direkt regressiert.
            resreg2 = lm(regres ~ 0 + Lag(-res2,1))}                        #d.h. alle Lags sind insignifikant
      }
      
      
      
      
      
      return(resreg2)
      
    }
    
    v = Filter_eghl(res2)                                                        
    resreg2 = lagsel2(-v[[8]]/v[[2]])                                        
    t2 = coef(summary(resreg2))[1,'t value']
    #Extracting names.
    lags2_names = names(coef(summary(resreg2))[,4])
    lags2 = lags2_names[-(length(lags2_names))]              
    lags2 = as.numeric(gsub("Lag", "", lags2))
    
    
    
    
    #Augmentation selection procedure for remaining frequencies.
    legsel3 = function(k){
      reg3 = lm(fily[[k]] ~ filx[[k]] + Lag(filx[[k]],1))
      res3 = resid(reg3)
      w = Filter_eghl(res3)
      
      regres = -w[[8]]/w[[k]]
      at = data.frame("Lag1" = Lag(regres,1) , "Lag2" =  Lag(regres,2) ,"Lag3" =  Lag(regres,3) , "Lag4" =  Lag(regres,4), "Lag5" =  Lag(regres,5) ,"Lag6" = Lag(regres,6) , "Lag7" = Lag(regres,7) ,"Lag8" = Lag(regres,8) ,"Lag9" = Lag(regres,9) ,"Lag10" = Lag(regres,10) ,"Lag11" = Lag(regres,11) , "Lag12" =Lag(regres,12) 
                      , "Lag13" = Lag(regres,13) ,  "Lag14" =Lag(regres,14) ,"Lag15" = Lag(regres,15) , "Lag16" =Lag(regres,16) , "Lag17" =Lag(regres,17) , "Lag18" =Lag(regres,18) , "Lag19" =Lag(regres,19) , "Lag20" =Lag(regres,20) , "Lag21" =Lag(regres,21) , "Lag22" =Lag(regres,22) , "Lag23" =Lag(regres,23) , "Lag24" =Lag(regres,24))
      name = names(at); p = 1
      for (i in 1:ncol(at)){ #ineffizient mit 24 Wdh.
        at = data.frame(at) #Grund: Letzter Durchlauf bei keinen Lags.
        if (ncol(at)==1 & min(p)<=0.05){
          colnames(at)=name}
        resreg3 = lm(regres ~. + 0 + Lag(-res3,2) + Lag(-res3,1), at)
        p = coef(summary(resreg3))[-((ncol(at)+1):(ncol(at)+2)),4]
        weg = which(p==max(p))
        
        if (max(p) > 0.05 & ncol(at)>1){
          at = at[,-weg]; name = name[-weg]} else if (max(p) > 0.05 & ncol(at)==1){#Kürzen des Vectors zu einem Null vector wird übersprungen und direkt regressiert.
            resreg3 = lm(regres ~ 0 + Lag(-res3,2) + Lag(-res3,1))                   #d.h. alle Lags sind insignifikant
          }
      }
      
      return(resreg3)
    }
    
    #Frequency = 3/12; 9/12 or pi/2; -pi/2.
    resreg3 = legsel3(3)                                    
    fstat3 = summary(resreg3)$f[1]
    #Extracting names.                                                              
    lags3_names = names(coef(summary(resreg3))[,4])
    lags3 = lags3_names[-c((length(lags3_names)),(length(lags3_names)-1))]  
    lags3 = as.numeric(gsub("Lag", "", lags3))
    
    
    
    #Frequency = 4/12; 8/12 or 2pi/3; -2pi/3.
    resreg5 = legsel3(4)                
    fstat5 = summary(resreg5)$f[1]
    lags5_names = names(coef(summary(resreg5))[,4])
    lags5 = lags5_names[-c((length(lags5_names)),(length(lags5_names)-1))]  
    lags5 = as.numeric(gsub("Lag", "", lags5))
    
    
    #Frequency = 2/12; 10/12 or pi/3; -pi/3
    resreg7 = legsel3(5)
    fstat7 = summary(resreg7)$f[1]
    lags7_names = names(coef(summary(resreg7))[,4])
    lags7 = lags7_names[-c((length(lags7_names)),(length(lags7_names)-1))]  
    lags7 = as.numeric(gsub("Lag", "", lags7))
    
    
    #Frequency = 5/12; 7/12 or 5pi/6; -5pi/6
    resreg9 = legsel3(6)
    fstat9 = summary(resreg9)$f[1]
    lags9_names = names(coef(summary(resreg9))[,4])
    lags9 = lags9_names[-c((length(lags9_names)),(length(lags9_names)-1))]  
    lags9 = as.numeric(gsub("Lag", "", lags9))
    
    
    #Frequency = 1/12; 11/12 or pi/6; -pi/6
    resreg11 = legsel3(7)
    fstat11 = summary(resreg11)$f[1]
    lags11_names = names(coef(summary(resreg11))[,4])
    lags11 = lags11_names[-c((length(lags11_names)),(length(lags11_names)-1))]  
    lags11 = as.numeric(gsub("Lag", "", lags11))
    
    
    teststats2 = list(t1 = t1, t2 = t2, f34 = fstat3, f56 = fstat5, f78 = fstat7, f910 = fstat9, f1112 = fstat11)
    
    
    #Liste mit allen lags.
    lags2=list(lags1,lags2,lags3,lags5,lags7,lags9,lags11)
    
    output2 = list(teststats2,lags2)
    return(output2)
  }
  
  teststat1x = int(x)[[1]]            
  teststat1y = int(y)[[1]]
  teststat2 = EGHL(x,y)[[1]]
  
  if (option == "manual"){
    
    #HEGY:manual
    distr_man = simulation(EGHL(x,y)[[2]])
    distr2_man = simulation2(int(x)[[2]])
    distr3_man = simulation2(int(y)[[2]])
    
    pval1x = pval(teststat1x$t1, distr2_man[,1], twosided = F, side="left")
    pval2x = pval(teststat1x$t2, distr2_man[,2], twosided = F, side="left")
    pval3x = pval(teststat1x$f34, distr2_man[,3], twosided = F, side="right")
    pval4x = pval(teststat1x$f56, distr2_man[,4], twosided = F, side="right")
    pval5x = pval(teststat1x$f78, distr2_man[,5], twosided = F, side="right")
    pval6x = pval(teststat1x$f910, distr2_man[,6], twosided = F, side="right")
    pval7x = pval(teststat1x$f1112, distr2_man[,7], twosided = F, side="right")
    pvalx = list(round(pval1x,3),round(pval2x,3),round(pval3x,3),round(pval4x,3),round(pval5x,3),round(pval6x,3),round(pval7x,3))     
    
    pval1y = pval(teststat1y$t1, distr3_man[,1], twosided = F, side="left")
    pval2y = pval(teststat1y$t2, distr3_man[,2], twosided = F, side="left")
    pval3y = pval(teststat1y$f34, distr3_man[,3], twosided = F, side="right")
    pval4y = pval(teststat1y$f56, distr3_man[,4], twosided = F, side="right")
    pval5y = pval(teststat1y$f78, distr3_man[,5], twosided = F, side="right")
    pval6y = pval(teststat1y$f910, distr3_man[,6], twosided = F, side="right")
    pval7y = pval(teststat1y$f1112, distr3_man[,7], twosided = F, side="right")
    pvaly = list(round(pval1y,3),round(pval2y,3),round(pval3y,3),round(pval4y,3),round(pval5y,3),round(pval6y,3),round(pval7y,3))     
    
    #EGHL:manual
    
    pval1 = pval(teststat2$t1, distr_man[,1], twosided=F, side="left")
    pval2 = pval(teststat2$t2, distr_man[,2], twosided=F, side="left")
    pval3 = pval(teststat2$f34, distr_man[,3], twosided=F, side="right")
    pval4 = pval(teststat2$f56, distr_man[,4], twosided=F, side="right")
    pval5 = pval(teststat2$f78, distr_man[,5], twosided=F, side="right")
    pval6 =  pval(teststat2$f910, distr_man[,6], twosided=F, side="right")
    pval7 =  pval(teststat2$f1112, distr_man[,7], twosided=F, side="right")
    
    pval = list(round(pval1,3),round(pval2,3),round(pval3,3),round(pval4,3),round(pval5,3),round(pval6,3),round(pval7,3))     
    
  } else if (option == "default") {
    
    #HEGY:default
    pval1x = interpolation(teststat1x$t1, distr, side="left", 1)
    pval2x = interpolation(teststat1x$t2, distr, side="left", 2)
    pval3x = interpolation(teststat1x$f34, distr, side="right", 3)
    pval4x = interpolation(teststat1x$f56, distr, side="right", 4)
    pval5x = interpolation(teststat1x$f78, distr, side="right", 5)
    pval6x = interpolation(teststat1x$f910, distr, side="right", 6)
    pval7x = interpolation(teststat1x$f1112, distr, side="right", 7)
    pvalx = list(round(pval1x,3),round(pval2x,3),round(pval3x,3),round(pval4x,3),round(pval5x,3),round(pval6x,3),round(pval7x,3))        
    
    
    pval1y = interpolation(teststat1y$t1, distr, side="left", 1)
    pval2y = interpolation(teststat1y$t2, distr, side="left", 2)
    pval3y = interpolation(teststat1y$f34, distr, side="right", 3)
    pval4y = interpolation(teststat1y$f56, distr, side="right", 4)
    pval5y = interpolation(teststat1y$f78, distr, side="right", 5)
    pval6y = interpolation(teststat1y$f910, distr, side="right", 6)
    pval7y = interpolation(teststat1y$f1112, distr, side="right", 7)
    pvaly = list(round(pval1y,3),round(pval2y,3),round(pval3y,3),round(pval4y,3),round(pval5y,3),round(pval6y,3),round(pval7y,3))       
    
    #EGHL:default
    pval1 = interpolation(teststat2$t1, distr2, side="left", 1)
    pval2 = interpolation(teststat2$t2, distr2, side="left", 2)
    pval3 = interpolation(teststat2$f34, distr2, side="right", 3)
    pval4 = interpolation(teststat2$f56, distr2, side="right", 4)
    pval5 = interpolation(teststat2$f78, distr2, side="right", 5)
    pval6 = interpolation(teststat2$f910, distr2, side="right", 6)
    pval7 = interpolation(teststat2$f1112, distr2, side="right", 7)        
    
    pval = list(round(pval1,3),round(pval2,3),round(pval3,3),round(pval4,3),round(pval5,3),round(pval6,3),round(pval7,3))          
    
    
  }
  
  
  #######################################################################
  
  #Return all p-values.
  pvalues = list("X" = pvalx, "Y" = pvaly, "XY" = pval)
  class(pvalues) = "coint"
  return(pvalues)
}


#Ausgaben:
a=scoint(x,y,option="default")
print.coint = function(x){
  fre = c("0","6/12","3/12 and 9/12","4/12 and 8/12","2/12 and 10/12","5/12 and 7/12","1/12 and 11/12")
  for (i in 1:7){
    if (x$X[[i]] > 0.05){
      cat("X is integrated at frequency", fre[i], "\n")
    } 
  }
  cat("\n")
  for (i in 1:7){
    if (x$Y[[i]] > 0.05){
      cat("Y is integrated at frequency",fre[i],"\n")
    } 
  }
  cat("\n")
  for (i in 1:7){
    if (x$XY[[i]] <= 0.05 & x$X[[i]] > 0.05 & x$Y[[i]] > 0.05){
      cat("X and Y are cointegrated at frequency", fre[i], "\n")
    }
  else if (x$XY[[i]] <= 0.05 & x$X[[i]] < 0.05 & x$Y[[i]] > 0.05 | x$XY[[i]] <= 0.05 & x$X[[i]] > 0.05 & x$Y[[i]] < 0.05 | x$XY[[i]] <= 0.05 & x$X[[i]] < 0.05 & x$Y[[i]] < 0.05){
    cat("( X and Y are cointegrated at frequency", fre[i],")", "\n")
  }
    
    }
}

summary.coint = function(x){
  result1 = NULL
  result2 = NULL
  result3 = NULL
  for (i in 1:7){
    if (x$X[[i]] > 0.05){
      result1[i] = "present"
    } else (result1[i] = "not present")
  }
  for (i in 1:7){
    if (x$Y[[i]] > 0.05){
      result2[i] = "present"
    } else (result2[i] = "not present")
  }
  for (i in 1:7){
    if (x$XY[[i]]<=0.05 && result1[i] == "present" && result2[i] == "present"){
      result3[[i]] = "present"
    }   else if (result1[i] == "present" && result2[i] == "not present" && x$XY[[i]]<=0.05 | result1[i] == "not present" && result2[i] == "present" && x$XY[[i]]<=0.05 | result1[i] == "not present" && result2[i] == "not present" && x$XY[[i]]<=0.05){
      result3[[i]] = "(present)"
    } else (result3[i] = "not present")
  }
  
  #Integration for X.
  cat("     X\t         p-value\t integration\n \nFrequency=0","\t", x$X[[1]],"\t\t", result1[1],"\n")
  cat("Frequency=6/12","\t", x$X[[2]],"\t\t", result1[2],"\n")
  cat("Frequency=3/12","\t", x$X[[3]],"\t\t", result1[3],"\n")
  
  cat("Frequency=4/12","\t", x$X[[4]],"\t\t", result1[4],"\n")
  cat("Frequency=2/12","\t", x$X[[5]],"\t\t", result1[5],"\n")
  
  cat("Frequency=5/12","\t", x$X[[6]],"\t\t", result1[6],"\n")
  cat("Frequency=1/12","\t", x$X[[7]],"\t\t", result1[7],"\n")
  
  
  #Integration for Y.
  cat("\n\n")
  
  cat("     Y\t         p-value\t integration\n \nFrequency=0","\t", x$Y[[1]],"\t\t", result2[1],"\n")
  cat("Frequency=6/12","\t", x$Y[[2]],"\t\t", result2[2],"\n")
  cat("Frequency=3/12","\t", x$Y[[3]],"\t\t", result2[3],"\n")
  
  cat("Frequency=4/12","\t", x$Y[[4]],"\t\t", result2[4],"\n")
  cat("Frequency=2/12","\t", x$Y[[5]],"\t\t", result2[5],"\n")
  
  cat("Frequency=5/12","\t", x$Y[[6]],"\t\t", result2[6],"\n")
  cat("Frequency=1/12","\t", x$Y[[7]],"\t\t", result2[7],"\n")
  
  
  
  #Integration for XY.
  cat("\n\n")
  
  cat("     XY\t         p-value\t cointegration\n \nFrequency=0","\t", x$XY[[1]],"\t\t", result3[1],"\n")
  cat("Frequency=6/12","\t", x$XY[[2]],"\t\t", result3[2],"\n")
  cat("Frequency=3/12","\t", x$XY[[3]],"\t\t", result3[3],"\n")
  
  cat("Frequency=4/12","\t", x$XY[[4]],"\t\t", result3[4],"\n")
  cat("Frequency=2/12","\t", x$XY[[5]],"\t\t", result3[5],"\n")
  
  cat("Frequency=5/12","\t", x$XY[[6]],"\t\t", result3[6],"\n")
  cat("Frequency=1/12","\t", x$XY[[7]],"\t\t", result3[7],"\n")
}
summary(a)


