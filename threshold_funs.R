compute_threshold <- function(scores, baseline, time_stamps, type = 'static', perc = 0.95){
  
  #Check that the scores, baseline and time_stamps are all the same length
  if(length(scores)!=length(baseline) || length(scores)!=length(time_stamps)){
    stop("The length of scores, baseline and time_stamps all need to match")
  }
  
  #Check `type` argument was correctly specified
  type <- match.arg(tolower(type),c("static","dynamic","both"))
  if(type=="both"){
    type <- c("static","dynamic")
  }
  
  if("dynamic"%in%type){
    if(!require(dlm)){
      warning("The 'dlm' package is missing and must be installed to compute dynamic thresholds, therefore only static threhsolds will be returned.")
      type <- "static"
    }
  }
  
  
  
  #Order the scores and baseline by time_stamps
  qc_data_frame <- data.frame(Scores=scores, Baseline=baseline, Time_Stamp=time_stamps)
  qc_data_frame <- qc_data_frame%>%arrange(Time_Stamp)%>%mutate(Time=difftime(Time_Stamp,min(Time_Stamp),units='hours'))
  qc_data_frame$Time <- as.numeric(qc_data_frame$Time)
  
  #filter down to baseline data
  qc_baseline <- filter(qc_data_frame, Baseline)
  
  if("static"%in%type){
    ##----------------------------------------##
    ##----- Static threshold computation  ----##
    ##----------------------------------------##
    
    #Try linear model in time
    lin_mod <- lm(log(Scores)~Time,data=qc_baseline)
    
    #If the linear model doesn't fit the data well, just use the intercept model
    if(anova(lin_mod)[1,5]>0.05){
      lin_mod <- lm(log(Scores)~1,data=qc_baseline)
      lnorm_mean <- unname(lin_mod$coefficients[1])
      lnorm_sd <- sqrt(anova(lin_mod)[1,3])
    }else{
      lnorm_mean <- lin_mod$coefficients[1]+lin_mod$coefficients[2]*qc_data_frame$Time
      lnorm_sd <- sqrt(anova(lin_mod)[2,3])*sqrt(1+1/nrow(qc_baseline))
    }
    
    #The threshold is the 
    qc_data_frame$Static_Threshold <- qlnorm(perc, meanlog = lnorm_mean, sdlog = lnorm_sd)
    static_model <- lin_mod
  }else{
    static_model <- NULL
  }
  
  if("dynamic"%in%type){
    
    ##----------------------------------------##
    ##---- Dynamic threshold computation  ----##
    ##----------------------------------------##
    #Define the DLM model
    buildCutpt <- function(theta,xvar){
      dlmModReg(X=xvar, dV=exp(theta[1]), dW=exp(theta[2:3]),addInt=TRUE) 
    }
    npars <- 3
    
    #Estimate parameters based on baseline data only
    fit <- dlmMLE(log(qc_baseline$Scores),parm=rep(0,npars),buildCutpt, xvar=matrix(qc_baseline$Time,ncol=1))
    
    #Now apply the model to the full dataset
    modCutpt <- buildCutpt(fit$par, xvar=qc_data_frame$Time)
    
    #filter estiamte for full data set
    filtQC <- dlmFilter(log(qc_data_frame$Scores),modCutpt)
    filtPred <- filtQC$f
    filtSD <- residuals(filtQC,sd=TRUE)$sd
    
    qc_data_frame$Dynamic_Threshold <- qlnorm(perc,meanlog=filtPred,sdlog=filtSD)

    dynamic_model <- filtQC
    
  }else{
    dynamic_model <- NULL
  }
  
  return(list(Results=qc_data_frame, Static_Model = static_model, Dynamic_Model = dynamic_model))
}

# ##------ Testing -------##
# library(dlm)
# scores <- vorb_04$QC_Art
# baseline <- vorb_04$Baseline
# time_stamps <- vorb_04$Acq_Time_Start
# 
# qc_data_frame <- compute_threshold(scores = scores, baseline = baseline, time_stamps = time_stamps, type='static')
# qc_baseline <- filter(qc_data_frame,Baseline)
# 
# buildCutpt <- function(theta,xvar){
#   dlmModReg(X=xvar, dV=exp(theta[1]), dW=exp(theta[2:3]),addInt=TRUE) #+ dlmModARMA(ar = c(theta[4]))
# }
# npars <- 3
# 
# #Estimate parameters based on baseline data only
# fit <- dlmMLE(log(qc_baseline$Scores),parm=rep(0,npars),buildCutpt, xvar=matrix(qc_baseline$Time,ncol=1))
# 
# #Now apply the model to the full dataset
# modCutpt <- buildCutpt(fit$par, xvar=qc_data_frame$Time)
# 
# #filter estiamte for full data set
# filtQC <- dlmFilter(log(qc_data_frame$Scores),modCutpt)
# filtPred <- filtQC$f
# filtSD <- residuals(filtQC,sd=TRUE)$sd
# 
# qc_data_frame$Dynamic_Threshold <- pmin(qlnorm(0.9,meanlog=filtPred,sdlog=filtSD),200)
# 
# 
# qplot(Time_Stamp,Scores,data=qc_data_frame,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
#   geom_line(aes(Time_Stamp,Dynamic_Threshold),colour=1)




