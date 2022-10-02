
##################################################################

# SECTION C: comparison of field evidence at resting sites and 
#non resting sites #

##################################################################

# VARIABLE DEFINITIONS: RS (RESTING SITE = 1, NON-RESTING SITE=0), 
#SPRAINT1M (spraint count 0.1m), SPILES.BIN (presence of spraint piles),
# RUN (presence of run), BED(presence of bedding), LAT (presence of latrine)


rm(list=ls())

# PACKAGES REQUIRED FOR THIS ANALYSIS

install.packages("MuMIn")

install.packages("logistf")

library(MuMIn)

library(logistf)


# LOAD "FIELDMEANSDATA.csv"

main<-read.csv("FIELDMEANSDATA.csv")


# ESTIMATE MEANS (OR MEDIAN) PER SITE (ACROSS ALL VISITS) 
#(SEE MAIN TEXT FOR EXPLANATION), GIVING SUFFIX "NEW"

# ASSIGN EACH SITE AS EITHER RESTING SITE OR NON RESTING SITE
RS.NEW<-tapply(main$RS,main$SITE.NO,mean) 

SPRAINT1M.NEW<-tapply(main$SPRAINT1M,main$SITE.NO,median)

SPRAINTALL.NEW<-tapply(main$SP.ALL,main$SITE.NO,mean)

BEDDING.NEW<-tapply(main$BEDDING,main$SITE.NO,mean)

LATRINE.NEW<-tapply(main$LATRINE,main$SITE.NO,mean)

RUN.NEW<-tapply(main$RUN,main$SITE.NO,mean)

SPILES.NEW<-tapply(main$SPILES,main$SITE.NO,mean)

SPILES.BIN.NEW<-tapply(main$SPILES.BIN,main$SITE.NO,mean)



#RUN LOGISTIC REGRESSION MODELS using Firth's penalized maximum likelihood approach

md1<-logistf(RS.NEW~1)

md2<-logistf(RS.NEW~LATRINE.NEW)

md3<-logistf(RS.NEW~BEDDING.NEW)

md4<-logistf(RS.NEW~RUN.NEW)

md5<-logistf(RS.NEW~SPRAINT1M.NEW)

md6<-logistf(RS.NEW~SPILES.BIN.NEW)



#CREATE MODEL SELECTION TABLE

mods.tab<-model.sel(md1,md2,md3,md4,md5,md6)

mods.tab



# To extract SEs for model parameters, use the summary() function

# Example with latrine model (md2)

summary(md2)



# For 95% CI of odds-ratio, take upper and lower 95% confidence limits of slope value
# from the model summary table and take exponential of these

# Example with latrine model (md2)

summary(md2)

exp(0.05608595)

exp(6.450600)

#######################################################################
#SECTION D. ESTIMATING THE OPTIMAL CAMERA-TRAP SAMPLING STRATEGY TO DETECT
#A RESTING EVENT
#######################################################################

#SIMULATIONS TO CALCULATE THE SAMPLING DURATION (NUMBER OF DAYS)
#FOR EACH SITE TO HAVE A 95% PROBABILITY OF DETECTING A REST.

#SIMULATION 1 SIMULATES A SINGLE SAMPLING DURATION ACROSS WINTER AND 
#THE FOLLOWING SPRING COMBINED, FOLLOWED BY SIMULATION 2 WHICH SIMULATES TWO
#EQUAL SAMPLING DURATIONS, ONE IN WINTER AND ONE IN SPRING.

##For each simulation, each site is plotted with sampling duration in days on the X axis and
#probability of detecting a rest on the Y axis

#Example data from one site that includes two site periods, 5a and 5b.

#####################################################################
#VARIABLE DEFINITIONS: 

#DATE.CYCLE: enables data selection for a single date cycle (i.e. consecutive
#winter and spring period)
#CAMS: camera traps functional =1, camera not functioning due to floods or fault = 0
#ALL.RESTS: rest detected =1, no rest detected = 0
#WINTER.DAY: numerical index starting with 1 on the first day of 
#winter and ascending, ceasing on last day of winter
#SPRING.DAY:numerical index starting with 1 on the first day of 
#spring and ascending, ceasing on last day of spring
#WINTER.SPRING.DAY: numerical index starting with 1 on the first day of 
#winter and ascending through winter and into spring, 
#ceasing on last day of spring

#########################################################################
#SIMULATION 1. FOR A SINGLE PERIOD OF CAMERA-TRAPPING DURING THE 
#WINTER-SPRING PERIOD

rm(list=ls())

#load data for site

loops<-read.csv("SITE.PERIOD5ab.csv")

#select date cycle
loops<-loops[loops$DATE.CYCLE==1,]

# SET THE THRESHOLD VALUE OF PROPORTION OF DAYS WHERE CTs SHOULD BE OPERATIONAL
#0/7 used here to include all data. IF SET AT 3/7 THIS WOULD BE THE MINIMUM PROPORTION
#OF DAYS THAT THE CAMS MUST BE FUNCTIONING I.E. 3 DAYS OUT OF 7. IF A SIMULATION
#DOES NOT MEET THIS THRESHOLD IT IS EXCLUDED. 
threshold<-(0/7)

# THIS IS THE MAXIMUM SAMPLE DURATION
max.CT.window<-nrow(loops)

ct.window<-0 # WILL BE POPULATED WITH ALL POSSIBLE CT WINDOWS (1 to longest)
true.pos.prop<-0
null.period.prop<-0

# THIS LOOP RUNS FOR EACH POSSIBLE CT WINDOW
for(j in 1:max.CT.window){
  
  # THIS JUST RECORDS THE CT WINDOW THAT j IS ON
  ct.window[j+1]<-j  
  
  # BLANK VECTOR RECORDING IF REST IS DETECTED ON GIVEN START DATE FOR WINTER AND SPRING 
  is.rest<-vector()
  
  CTdays.prop<-vector()
  
  # THIS LOOP RUNS FOR EACH POSSIBLE START DAY WITHIN SEASON
  for(i in 1:(nrow(loops)-(j-1))){
    is.rest[i]<-sum(loops$ALL.RESTS[i:(i+j-1)])>0;
    
    CTdays.prop[i]<-(sum(loops$CAMS[i:(i+j-1)])/j)>=threshold
    
  }
  
  is.rest[is.rest==0 & CTdays.prop!=1]<-2 
  
  # THIS ASSIGNS ONE OF THREE 3 codes: 0 = no rest detected, 1 = rest detected, 
  #2 = not enough CT dats to tell (i.e. below threshold days) 
  null.period<-length(is.rest[is.rest==2])
  true.positive<-length(is.rest[is.rest==1])
  false.negative<-length(is.rest[is.rest==0])
  
  true.pos.prop[j+1]<-(true.positive)/(true.positive+false.negative)
  null.period.prop[j+1]<-(null.period)/(true.positive+false.negative+null.period)
  
}

data1<-na.omit(data.frame(ct.window,round(true.pos.prop,3),round(null.period.prop,3)))
names(data1)<-c("CT.window","True.positives","Proportion.nulls")
data1
plot(data1$True.positives~data1$CT.window,type="l",lty=5,col="black",lwd=1.5,xlim=c(1,150), cex=2,xlab="Sampling duration(d) for a single camera-trapping period",ylab="Probability of detecting a rest")
abline(h=0.95,lty=3)
sum(data1$True.positives<0.95) # MIN SAMPLE DURATION (PER SEASON) FOR 95% PROPBABILITY OF DETECTING REST

#TO REPEAT FOR DATE CYCLE 2, 
#FIRST RUN
#rm(list=ls())
#THEN SUBSTITUTE FOLLOWING CODE TO REPLACE line 140
#loops<-loops[loops$DATE.CYCLE==2,] 

######################################################################
######################################################################

#SIMULATION 2. FOR TWO EQUAL PERIODS OF CAMERA-TRAPPING, ONE IN WINTER AND
#ONE IN THE FOLLOWING SPRING 

rm(list=ls())

#load data for site
#loops<-read.csv("SITE.PERIOD5ab.csv")

#CREATE SEASONAL VARIABLES, WINTER AND SPRING
winter<-loops[loops$SEASON=="WINTER",]
spring<-loops[loops$SEASON=="SPRING",]

#select date cycle
loops<-loops[loops$DATE.CYCLE==1,]

# SET THE THRESHOLD VALUE OF PROPORTION OF DAYS WHERE CTs SHOULD BE OPERATIONAL
#0/7 used here to include all data. IF SET AT 3/7 THIS WOULD BE THE MINIMUM PROPORTION
#OF DAYS THAT THE CAMS MUST BE FUNCTIONING I.E. 3 DAYS OUT OF 7. IF A SIMULATION
#DOES NOT MEET THIS THRESHOLD IT IS EXCLUDED. 
threshold<-(0/7)

# SET THE MAXIMUM SAMPLE DURATION PER SEASON IT IS WHICH EVER SEASON IS SHORTEST IN THE DATASET
max.CT.window<-min(c(length(winter$WINTER.DAY),length(spring$SPRING.DAY)))

# ct.window WILL BE POPULATED WITH ALL POSSIBLE CT WINDOWS (1 to longest)
ct.window<-0 
true.pos.prop<-0
null.period.prop<-0

# THIS LOOP RUNS FOR EACH POSSIBLE CT WINDOW
for(j in 1:max.CT.window){
  
  # THIS RECORDS THE CT WINDOW THAT j IS ON
  ct.window[j+1]<-j  
  
  # CREATE BLANK MATRICES THAT WILL BE POPULATED IF A REST IS DETECTED ON A GIVEN
  #START DATE FOR WINTER AND SPRING 
  #THESE WILL BE POPULATED AS THE NEXT LOOP OPERATES
  is.rest.spring<-matrix(ncol=nrow(winter),nrow=nrow(spring))
  is.rest.winter<-matrix(ncol=nrow(winter),nrow=nrow(spring))
  CTdays.spring.prop<-matrix(ncol=nrow(winter),nrow=nrow(spring))
  CTdays.winter.prop<-matrix(ncol=nrow(winter),nrow=nrow(spring))
  
  # THIS LOOP RUNS FOR EACH POSSIBLE START DAY WITHING WINTER
  for(k in 1:(length(winter$WINTER.DAY)-(j-1))){
    
    # THIS LOOP RUNS FOR EACH POSSIBLE START DAY WITHIN SPRING
    for(i in 1:(length(spring$SPRING.DAY)-(j-1))){
      is.rest.spring[i,k]<-sum(spring$ALL.RESTS[i:(i+j-1)])>0;
      is.rest.winter[i,k]<-sum(winter$ALL.RESTS[k:(k+j-1)])>0;
      
      CTdays.spring.prop[i,k]<-(sum(spring$CAMS[i:(i+j-1)])/j)>=threshold
      CTdays.winter.prop[i,k]<-(sum(winter$CAMS[k:(k+j-1)])/j)>=threshold
    }
    
  }
  
  is.rest.spring[CTdays.spring.prop!=1 & is.rest.spring<1]<-(3)
  is.rest.winter[CTdays.winter.prop!=1 & is.rest.winter<1]<-(3)
  
  false.negative<-sum((is.rest.spring+is.rest.winter)==0,na.rm=T)
  true.positive1<-sum((is.rest.spring+is.rest.winter)==1,na.rm=T)
  true.positive2<-sum((is.rest.spring+is.rest.winter)==2,na.rm=T)
  null.period<-sum((is.rest.spring+is.rest.winter)>2,na.rm=T)
  
  true.pos.prop[j+1]<-(true.positive1+true.positive2)/(true.positive1+true.positive2+false.negative)
  null.period.prop[j+1]<-(null.period)/(true.positive1+true.positive2+false.negative+null.period)
  
}

(data1<-na.omit(data.frame(ct.window,round(true.pos.prop,3),round(null.period.prop,3))))
names(data1)<-c("CT.window","True.positives","Proportion.nulls")

#AS THERE ARE TWO EQUAL PERIODS OF CAMERA TRAPPING, THE NUMBER OF CT DAYS MUST BE MULTIPLIED
#BY 2 TO MAKE THE PLOTS COMPARABLE IN TERMS OF THE TOTAL SURVEY EFFORT FOR EACH SAMPLING APPROACH
two.periods<-2*data1$CT.window
data2<-cbind(two.periods,data1)

plot(data2$True.positives~data2$two.periods,type="l", lty=5,col="black",lwd=1.5,xlim=c(1,150),cex=2, xlab="Total sampling duration(d) for two equal camera-trapping periods",ylab="Probability of detecting a rest")

abline(h=0.95,lty=3)
# MIN SAMPLE DURATION (PER SEASON) FOR 95% PROPBABILITY OF DETECTING REST
sum(data1$True.positives<0.95)
sum(data1$True.positives<1)

#TO REPEAT FOR DATE CYCLE 2, 
#FIRST RUN
#rm(list=ls())
#THEN SUBSTITUTE FOLLOWING CODE TO REPLACE line 216
#loops<-loops[loops$DATE.CYCLE==2,] 

######################################################################
#END

