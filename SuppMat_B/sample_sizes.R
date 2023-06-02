#######################################################################################################################################
#### Sample Size Calculation and Optimal Design for Regression-Based Norming of Tests and Questionnaires
#######################################################################################################################################
# Francesco Innocenti, Frans E. S. Tan, Math J. J. M. Candel, and Gerard J. P. van Breukelen
#############################################################################################################################################
############# Sample size Calculation procedure : R code to obtain Table 4 #################################################################
############################################################################################################################################
rm(list=ls(all=T))

# Sample size calculation formulae for hypothesis testing

N_star_HypTest_Z<-function(alpha,gamma,k,Zc,delta){ # equation (14)

if(Zc<0){

Nstar<-((qnorm(1-alpha)*sqrt(k+1+((Zc^2)/2))+qnorm(1-gamma)*sqrt(k+1+(((Zc-delta)^2)/2)))/(delta))^2 # Ha: Zt<Zc

}else{

Nstar<-((qnorm(1-alpha)*sqrt(k+1+((Zc^2)/2))+qnorm(1-gamma)*sqrt(k+1+(((Zc+delta)^2)/2)))/(delta))^2 # Ha: Zt>Zc
    
}
  return(Nstar)
}



N_star_HypTest_PR<-function(alpha,gamma,k,PRc,delta){ #equation (15)
  
if(PRc<50){
  # Ha: PRt<PRc
Nstar<-((qnorm(1-alpha)*100*dnorm(qnorm(p=(PRc)/100))*sqrt(k+1+((qnorm(p=(PRc)/100))^2/2))+qnorm(1-gamma)*100*dnorm(qnorm(p=(PRc-delta)/100))*sqrt(k+1+((qnorm(p=(PRc-delta)/100))^2/2)))/delta)^2 
    
}else{
  # Ha: PRt>PRc
Nstar<-((qnorm(1-alpha)*100*dnorm(qnorm(p=(PRc)/100))*sqrt(k+1+((qnorm(p=(PRc)/100))^2/2))+qnorm(1-gamma)*100*dnorm(qnorm(p=(PRc+delta)/100))*sqrt(k+1+((qnorm(p=(PRc+delta)/100))^2/2)))/delta)^2 
    
  }
  return(Nstar)
}

################################################### Z-scores #############################################################
delta<-c(0.4,0.3,0.2)
Zc<-c(1.5,2,2.5)
Zcneg<--Zc

k<-c(2,3,4,5)

alpha<-0.05 # 0.01
gamma<-0.2 # 0.10
Z_alpha<-qnorm(1-alpha)
z_power<-qnorm(1-gamma)

N_Z<-N_Zneg<-array(0,dim=c(length(k),length(delta),length(Zc)))
# rows = k, column= delta, matrix= Z
for(l in 1:length(Zc)){
for(j in 1:length(delta)){
for(i in 1:length(k)){
N_Z[i,j,l]<-N_star_HypTest_Z(alpha=alpha,gamma=gamma,k=k[i],Zc=Zc[l],delta=delta[j])
N_Zneg[i,j,l]<-N_star_HypTest_Z(alpha=alpha,gamma=gamma,k=k[i],Zc=Zcneg[l],delta=delta[j])
}}}
all(N_Z==N_Zneg)
round(N_Z)<338
N_Z<-round(N_Z)
colnames(N_Z)<-c("delta=0.4","delta=0.3","delta=0.2")
N_Z

# Z=+/-1.5
data.frame(N_Z[,,1], row.names = c("model (1)","models (2) and (3)","model (4)","model (5)"))
# Z=+/-2
data.frame(N_Z[,,2], row.names = c("model (1)","models (2) and (3)","model (4)","model (5)"))
# Z=+/-2.5
data.frame(N_Z[,,3], row.names = c("model (1)","models (2) and (3)","model (4)","model (5)"))

#################################################### PR-scores ########################################################## 
delta<-c(2,1.5,1) 
PRc<-c(2.5,5,10)
PRabove50<-100-PRc

k<-c(2,3,4,5)

alpha<-0.05 # 0.01
gamma<-0.2 # 0.10
Z_alpha<-qnorm(1-alpha)
z_power<-qnorm(1-gamma)

N_PR<-N_PRabove50<-array(0,dim=c(length(k),length(delta),length(PRc)))
# rows = k, column= delta, matrix= PR

for(l in 1:length(PRc)){
for(j in 1:length(delta)){
for(i in 1:length(k)){
N_PR[i,j,l]<-N_star_HypTest_PR(alpha=alpha,gamma=gamma,k=k[i],PRc=PRc[l],delta=delta[j])
N_PRabove50[i,j,l]<-N_star_HypTest_PR(alpha=alpha,gamma=gamma,k=k[i],PRc=PRabove50[l],delta=delta[j])
}}}
all(round(N_PR,3)==round(N_PRabove50,3))
round(N_PR)<1690
N_PR<-round(N_PR)
N_PR

colnames(N_PR)<-c("delta=2","delta=1.5","delta=1")
N_PR

# PR=2.5 or 97.5
data.frame(N_PR[,,1], row.names = c("model (1)","models (2) and (3)","model (4)","model (5)"))
# PR=5 or 95
data.frame(N_PR[,,2], row.names = c("model (1)","models (2) and (3)","model (4)","model (5)"))
# PR=10 or 90
data.frame(N_PR[,,3], row.names = c("model (1)","models (2) and (3)","model (4)","model (5)"))

