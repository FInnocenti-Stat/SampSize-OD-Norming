#######################################################################################################################################
#### Sample Size Calculation and Optimal Design for Regression-Based Norming of Tests and Questionnaires
#######################################################################################################################################
# Francesco Innocenti, Frans E. S. Tan, Math J. J. M. Candel, and Gerard J. P. van Breukelen
#######################################################################################################################################
rm(list=ls(all=T))
set.seed(05092016)
#######################################################################################################################################
########################################## Simulation study for the variance estimators ###############################################
##############################  step 0: set parameters and Z values ##################################################################### 
#      Intercept Age    Sex   Age^2 AgexSex Age^2xSex
Beta<-c(18.579,-0.054,0.753,-0.002, -0.03,-0.001); sigma_epsilon<-4.818  # PNVFT

# number of regression coefficients in each of the five models
k_mod1<-3
k_mod2<-k_mod3<-4
k_mod4<-5
k_mod5<-6

Z<-seq(from=-3,to=3,by=0.5)
Z
M<-length(Z)
M
Epsilon<-sigma_epsilon*Z
Epsilon

Phi<-pnorm(Z)*100
Phi

AGE<-seq(from=20,to=80,by=5) 
Age<-AGE-mean(AGE)  # age is centered 
Age

Age_2<-Age^2
Age_2

Sex<-c(0,1)
Sex

Age_Sex<-c(Age*Sex[1],Age*Sex[2])           # age sex interaction
Age_Sex

Age_2_Sex<-c(Age_2*Sex[1],Age_2*Sex[2])     # age^2 sex interaction
Age_2_Sex


X_comb<-expand.grid(Age,Sex) # combinations of all levels of age and sex
X_comb
colnames(X_comb)<-c("Age","Sex")
Q<-dim(X_comb)[1]
Q

# for each combination of age-sex levels: residual errors, Z scores, and Percentile rank scores
Xe_comb<-data.frame(rep(X_comb$Age, times=T),rep(X_comb$Sex, times=T),rep(Epsilon,each=Q),rep(Z,each=Q),rep(Phi,each=Q))
Xe_comb
colnames(Xe_comb)<-c("Age","Sex","Epsilon","Z","Phi")
N1<-dim(Xe_comb)[1]
N1
N1==Q*M

# design matrix for model (5)
# introduce L replicate in each cell Age x Sex
L<-1 # 2, 5, 6, 10
X<-data.frame(1,rep(Xe_comb$Age,times=L),rep(Xe_comb$Sex, times=L), rep(Age_2, times=M*L), rep(Age_Sex, times=M*L), rep(Age_2_Sex, times=M*L) )
colnames(X)<-c("Int","Age","Sex","Age2","AgeSex","Age2Sex")
X

dim(X)[1]==Q*M*L
dim(X[Xe_comb$Age==20 & Xe_comb$Sex==0.5,])[1]==M*L

N<-Q*M*L
N  
####################################################################################################################################################
############################################# model 1 #############################################################################################
####################################################################################################################################################

##############################  step 1: generate the ys ############################################################################################
y_mod1<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+rep(Xe_comb$Epsilon,times=L)
y_mod1

###############################   step 2: compute equations (7) and (8)  ##########################################################################
Var.Z_DM_mod1<-diag(data.matrix(X[,1:3]) %*% solve(t(data.matrix(X[,1:3]))%*%data.matrix(X[,1:3])) %*% t(data.matrix(X[,1:3]))+(1/(2*(N-k_mod1)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod1

Var.Phi_DM_mod1<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod1
Var.Phi_DM_mod1

############################# step 3: generate S normative samples ###################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod1<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Epsilon_s
#############################  step 4: estimate the model parameters  ##############################################################################
Y_hat_s_mod1<-Z_hat_s_mod1<-Phi_hat_s_mod1<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod1<-matrix(0,1,ncol = S)
Betas_hat_mod1<-matrix(0,nrow=k_mod1,ncol = S)
Var.Z_hat_L_mod1<-Var.Phi_hat_L_mod1<-matrix(0,nrow=N,ncol = S)
Var.Z_hat_mod1<-Var.Phi_hat_mod1<-matrix(0,Q*length(Z),S)

for(i in 1:S){
  
model <-lm(y_s_mod1[,i] ~ Age+Sex, data=data.frame(X,y_s_mod1[,i])) 
Betas_hat_mod1[,i]<-model$coefficients
Y_hat_s_mod1[,i]<-model$fitted.values
sigma_epsilon_hat_mod1[i]<-sigma(model)  

#############################   step 5: estimate Z-scores, compute equations (9) and (10)   ##########################################################
Z_hat_s_mod1[,i]<-(y_mod1-Y_hat_s_mod1[,i])/sigma_epsilon_hat_mod1[i]
Phi_hat_s_mod1[,i]<-pnorm(Z_hat_s_mod1[,i])*100
Var.Z_hat_L_mod1[,i]<-diag(data.matrix(X[,1:3]) %*% solve(t(data.matrix(X[,1:3]))%*%data.matrix(X[,1:3])) %*% t(data.matrix(X[,1:3]))+(1/(2*(N-k_mod1)))* Z_hat_s_mod1[,i] %*% t(Z_hat_s_mod1[,i]))
Var.Phi_hat_L_mod1[,i]<-(100*dnorm(Z_hat_s_mod1[,i]))^2*Var.Z_hat_L_mod1[,i]
  
  
}

##############################    step 6: compute Monte Carlo averages and true variances   ###################################################################
Var.Z_S_mod1<-apply(Z_hat_s_mod1,MARGIN=1,FUN=var)  
Var.Phi_S_mod1<-apply(Phi_hat_s_mod1,MARGIN=1,FUN=var)
Var.Z_mod1<-Var.Phi_mod1<-matrix(0,Q*length(Z),1)

# take average across L replications (Note that they are all equals)
for(i in 1:S){
  
t<-0
  
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(w in 1:length(Age)){
      
t<-t+1
Var.Z_hat_mod1[t,i]<-mean(Var.Z_hat_L_mod1[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
Var.Phi_hat_mod1[t,i]<-mean(Var.Phi_hat_L_mod1[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])

    }}}}

t<-0 # take average across L replications
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
Var.Z_mod1[t]<-mean(Var.Z_S_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod1[t]<-mean(Var.Phi_S_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
}}}

# MC average and SD for equation (9)
Mean_Var.Z_hat_mod1<-apply(Var.Z_hat_mod1,MARGIN=1,FUN=mean) 
SD_Var.Z_hat_mod1<-apply(Var.Z_hat_mod1,MARGIN=1,FUN=sd)  
summary(SD_Var.Z_hat_mod1)
# MC average and SD for equation (10)
Mean_Var.Phi_hat_mod1<-apply(Var.Phi_hat_mod1,MARGIN=1,FUN=mean)
SD_Var.Phi_hat_mod1<-apply(Var.Phi_hat_mod1,MARGIN=1,FUN=sd)
summary(SD_Var.Phi_hat_mod1)


###############################     step 7: compute Bias and coverage rate   #####################################################################
# R.B. for V(Z_hat)_hat
RB_Var.Z_hat_wrt_TV_mod1_L_1<-(Mean_Var.Z_hat_mod1-Var.Z_mod1)/Var.Z_mod1
summary(RB_Var.Z_hat_wrt_TV_mod1_L_1)
# R.B. for V(PR(Z_hat))_hat
RB_Var.PR_hat_wrt_TV_mod1_L_1<-(Mean_Var.Phi_hat_mod1-Var.Phi_mod1)/Var.Phi_mod1
summary(RB_Var.PR_hat_wrt_TV_mod1_L_1)  

## coverage rate
rm(Mean_Var.Z_hat_mod1, SD_Var.Z_hat_mod1, Mean_Var.Phi_hat_mod1, SD_Var.Phi_hat_mod1,Y_hat_s_mod1,sigma_epsilon_hat_mod1,y_s_mod1,
   Var.Z_DM_mod1_L, Var.Phi_DM_mod1_L, 
   MC.CV_Var.Phi_hat_mod1_L_1, MC.CV_Var.Z_hat_mod1_L_1, Epsilon_s, RB_Var.Z_hat_mod1_L_1, RB_Var.Phi_hat_mod1_L_1)

# Z score 
InfL_Z<-Z_hat_s_mod1-qnorm(0.975)*sqrt(Var.Z_hat_L_mod1) # lower-bound
SupL_Z<-Z_hat_s_mod1+qnorm(0.975)*sqrt(Var.Z_hat_L_mod1) # upper-bound
# PR score
InfL_PR<-Phi_hat_s_mod1-qnorm(0.975)*sqrt(Var.Phi_hat_L_mod1) # lower-bound
SupL_PR<-Phi_hat_s_mod1+qnorm(0.975)*sqrt(Var.Phi_hat_L_mod1) # upper-bound

cbind(rep(Xe_comb$Epsilon,times=L),rep(Xe_comb$Z,times=L))
Ztrue<-rep(Xe_comb$Z,times=L) # it should be like # rep(Xe_comb$Epsilon,times=L)
PRtrue<-pnorm(Ztrue)*100

# indicator variable: 1 if true parameter is within the confidence interval, 0 otherwilse
Indicator_L_Z<-Indicator_L_PR<-matrix(0,nrow=N,ncol = S)

for(j in 1:N){for(i in 1:S){
  ifelse(Ztrue[j]>=InfL_Z[j,i] & Ztrue[j]<=SupL_Z[j,i], Indicator_L_Z[j,i]<-1, Indicator_L_Z[j,i]<-0)
  ifelse(PRtrue[j]>=InfL_PR[j,i] & PRtrue[j]<=SupL_PR[j,i], Indicator_L_PR[j,i]<-1, Indicator_L_PR[j,i]<-0)
}}

# average of the indicator variable across S normative sample
coverage_rate_Z_S<-apply(Indicator_L_Z,MARGIN=1,FUN=mean)
coverage_rate_PR_S<-apply(Indicator_L_PR,MARGIN=1,FUN=mean)
# average of the indicator variable across L replications 
coverage_rate_Z_mod1_L_1<-coverage_rate_PR_mod1_L_1<-matrix(0,Q*length(Z),1)
t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      coverage_rate_Z_mod1_L_1[t]<-mean(coverage_rate_Z_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      coverage_rate_PR_mod1_L_1[t]<-mean(coverage_rate_PR_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
    }}}
#coverage rate
coverage_rate_Z_mod1_L_1;summary(coverage_rate_Z_mod1_L_1)
coverage_rate_PR_mod1_L_1;summary(coverage_rate_PR_mod1_L_1)


### Z-scores
## Figure S.A.11 Relative Bias of Equation (9)
# Sex=0
FigureS.A.11_Sex0<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                   RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.11_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.11_Sex0
# Sex=1
FigureS.A.11_Sex1<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                   RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.11_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.11_Sex1

## Figure S.A.13 coverage rate 
# Sex=0
FigureS.A.13_Sex0<-round(cbind(AGE,coverage_rate_Z_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                   coverage_rate_Z_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.13_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.13_Sex0
# Sex=1
FigureS.A.13_Sex1<-round(cbind(AGE,coverage_rate_Z_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                   coverage_rate_Z_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_Z_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.13_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.13_Sex1

### PR-scores  # recall that for PR-scores L=2,5,6,10 and N=676,1690,2028,3380
## Figure S.B.1.35  Relative Bias of Equation (10)
# Sex=0
FigureS.B.1.35_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                     RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.35_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.35_Sex0
# Sex=1
FigureS.B.1.35_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                     RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.35_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.35_Sex1

## Figure S.B.1.40 coverage rate
# Sex=0
FigureS.B.1.40_Sex0<-round(cbind(AGE,coverage_rate_PR_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.40_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.40_Sex0
# Sex=1
FigureS.B.1.40_Sex1<-round(cbind(AGE,coverage_rate_PR_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_PR_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.40_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.40_Sex1

####################################################################################################################################################
############################################# model 2 #############################################################################################
####################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)

##############################  step 1: generate the ys ###########################################################################################
y_mod2<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+rep(Xe_comb$Epsilon,times=L) 
y_mod2

###############################   step 2: compute equations (7) and (8)  ##########################################################################
Var.Z_DM_mod2<-diag(data.matrix(X[,1:4]) %*% solve(t(data.matrix(X[,1:4]))%*%data.matrix(X[,1:4])) %*% t(data.matrix(X[,1:4]))+(1/(2*(N-k_mod2)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod2

Var.Phi_DM_mod2<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod2
Var.Phi_DM_mod2


############################# step 3: generate S normative samples ################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod2<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Epsilon_s 
#############################  step 4: estimate the model parameters ##############################################################################
Y_hat_s_mod2<-Z_hat_s_mod2<-Phi_hat_s_mod2<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod2<-matrix(0,1,ncol = S)
Var.Z_hat_L_mod2<-Var.Phi_hat_L_mod2<-matrix(0,nrow=N,ncol = S)
Var.Z_hat_mod2<-Var.Phi_hat_mod2<-matrix(0,Q*length(Z),S)
Betas_hat_mod2<-matrix(0,nrow=k_mod2,ncol = S)


for(i in 1:S){
  
  model <-lm(y_s_mod2[,i] ~ Age+Sex+Age2, data=data.frame(X,y_s_mod2[,i])) 
  Betas_hat_mod2[,i]<-model$coefficients
    Y_hat_s_mod2[,i]<-model$fitted.values
  sigma_epsilon_hat_mod2[i]<-sigma(model)  
  
#############################   step 5: estimate  Z-scores, compute equations (9) and (10)   #####################################################
  Z_hat_s_mod2[,i]<-(y_mod2-Y_hat_s_mod2[,i])/sigma_epsilon_hat_mod2[i]
  Phi_hat_s_mod2[,i]<-pnorm(Z_hat_s_mod2[,i])*100
  Var.Z_hat_L_mod2[,i]<-diag(data.matrix(X[,1:4]) %*% solve(t(data.matrix(X[,1:4]))%*%data.matrix(X[,1:4])) %*% t(data.matrix(X[,1:4]))+(1/(2*(N-k_mod2)))* Z_hat_s_mod2[,i] %*% t(Z_hat_s_mod2[,i]))
  Var.Phi_hat_L_mod2[,i]<-(100*dnorm(Z_hat_s_mod2[,i]))^2*Var.Z_hat_L_mod2[,i]
  
}

##############################    step 6: compute Monte Carlo averages and true variances   ###################################################################
Var.Z_S_mod2<-apply(Z_hat_s_mod2,MARGIN=1,FUN=var)  
Var.Phi_S_mod2<-apply(Phi_hat_s_mod2,MARGIN=1,FUN=var)
Var.Z_mod2<-Var.Phi_mod2<-matrix(0,Q*length(Z),1)

# take average across L replications (Note that they are all equals)
for(i in 1:S){
  
  t<-0
  
  for(l in 1:length(Z)){
    for(j in 1:length(Sex)){
      for(w in 1:length(Age)){
        
        t<-t+1
        Var.Z_hat_mod2[t,i]<-mean(Var.Z_hat_L_mod2[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
        Var.Phi_hat_mod2[t,i]<-mean(Var.Phi_hat_L_mod2[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
      }}}}


t<-0 # take average across L replications
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
t<-t+1
Var.Z_mod2[t]<-mean(Var.Z_S_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod2[t]<-mean(Var.Phi_S_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
 }}}


# MC average and SD for equation (9)
Mean_Var.Z_hat_mod2<-apply(Var.Z_hat_mod2,MARGIN=1,FUN=mean) 
SD_Var.Z_hat_mod2<-apply(Var.Z_hat_mod2,MARGIN=1,FUN=sd)  
summary(SD_Var.Z_hat_mod2)
# MC average and SD for equation (10)
Mean_Var.Phi_hat_mod2<-apply(Var.Phi_hat_mod2,MARGIN=1,FUN=mean)
SD_Var.Phi_hat_mod2<-apply(Var.Phi_hat_mod2,MARGIN=1,FUN=sd)
summary(SD_Var.Phi_hat_mod2)


###############################     step 7: compute Bias and coverage rate   #####################################################################
# R.B. for V(Z_hat)_hat
RB_Var.Z_hat_wrt_TV_mod2_L_1<-(Mean_Var.Z_hat_mod2-Var.Z_mod2)/Var.Z_mod2
summary(RB_Var.Z_hat_wrt_TV_mod2_L_1)
# R.B. for V(PR(Z_hat))_hat
RB_Var.PR_hat_wrt_TV_mod2_L_1<-(Mean_Var.Phi_hat_mod2-Var.Phi_mod2)/Var.Phi_mod2
summary(RB_Var.PR_hat_wrt_TV_mod2_L_1)  

## coverage rate
rm(Mean_Var.Z_hat_mod2, SD_Var.Z_hat_mod2, Mean_Var.Phi_hat_mod2, SD_Var.Phi_hat_mod2,Y_hat_s_mod2,sigma_epsilon_hat_mod2,y_s_mod2,
   Var.Z_DM_mod2_L, Var.Phi_DM_mod2_L, 
   MC.CV_Var.Phi_hat_mod2_L_1, MC.CV_Var.Z_hat_mod2_L_1, Epsilon_s, RB_Var.Z_hat_mod2_L_1, RB_Var.Phi_hat_mod2_L_1)

# Z score 
InfL_Z<-Z_hat_s_mod2-qnorm(0.975)*sqrt(Var.Z_hat_L_mod2) # lower-bound
SupL_Z<-Z_hat_s_mod2+qnorm(0.975)*sqrt(Var.Z_hat_L_mod2) # upper-bound
# PR score
InfL_PR<-Phi_hat_s_mod2-qnorm(0.975)*sqrt(Var.Phi_hat_L_mod2) # lower-bound
SupL_PR<-Phi_hat_s_mod2+qnorm(0.975)*sqrt(Var.Phi_hat_L_mod2) # upper-bound

cbind(rep(Xe_comb$Epsilon,times=L),rep(Xe_comb$Z,times=L))
Ztrue<-rep(Xe_comb$Z,times=L) # it should be like # rep(Xe_comb$Epsilon,times=L)
PRtrue<-pnorm(Ztrue)*100

# indicator variable: 1 if true parameter is within the confidence interval, 0 otherwilse
Indicator_L_Z<-Indicator_L_PR<-matrix(0,nrow=N,ncol = S)

for(j in 1:N){for(i in 1:S){
  ifelse(Ztrue[j]>=InfL_Z[j,i] & Ztrue[j]<=SupL_Z[j,i], Indicator_L_Z[j,i]<-1, Indicator_L_Z[j,i]<-0)
  ifelse(PRtrue[j]>=InfL_PR[j,i] & PRtrue[j]<=SupL_PR[j,i], Indicator_L_PR[j,i]<-1, Indicator_L_PR[j,i]<-0)
}}

# average of the indicator variable across S normative sample
coverage_rate_Z_S<-apply(Indicator_L_Z,MARGIN=1,FUN=mean)
coverage_rate_PR_S<-apply(Indicator_L_PR,MARGIN=1,FUN=mean)
# average of the indicator variable across L replications 
coverage_rate_Z_mod2_L_1<-coverage_rate_PR_mod2_L_1<-matrix(0,Q*length(Z),1)
t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      coverage_rate_Z_mod2_L_1[t]<-mean(coverage_rate_Z_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      coverage_rate_PR_mod2_L_1[t]<-mean(coverage_rate_PR_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
    }}}
#coverage rate
coverage_rate_Z_mod2_L_1;summary(coverage_rate_Z_mod2_L_1)
coverage_rate_PR_mod2_L_1;summary(coverage_rate_PR_mod2_L_1)

### Z-scores
## Figure S.B.1.13 Relative Bias of Equation (9)
# Sex=0
FigureS.B.1.13_Sex0<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.13_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.13_Sex0
# Sex=1
FigureS.B.1.13_Sex1<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.13_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.13_Sex1

## Figure S.B.1.17 coverage rate 
# Sex=0
FigureS.B.1.16_Sex0<-round(cbind(AGE,coverage_rate_Z_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 coverage_rate_Z_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.16_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.16_Sex0
# Sex=1
FigureS.B.1.16_Sex1<-round(cbind(AGE,coverage_rate_Z_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 coverage_rate_Z_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_Z_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.16_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.16_Sex1

### PR-scores # recall that for PR-scores L=2,5,6,10 and N=676,1690,2028,3380
## Figure S.B.1.36  Relative Bias of Equation (10)
# Sex=0
FigureS.B.1.36_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.36_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.36_Sex0
# Sex=1
FigureS.B.1.36_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.36_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.36_Sex1

## Figure S.B.1.41 coverage rate
# Sex=0
FigureS.B.1.41_Sex0<-round(cbind(AGE,coverage_rate_PR_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.41_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.41_Sex0
# Sex=1
FigureS.B.1.41_Sex1<-round(cbind(AGE,coverage_rate_PR_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_PR_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.41_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.41_Sex1
####################################################################################################################################################
############################################# model 3 ##############################################################################################
####################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)

##############################  step 1: generate the ys ############################################################################################
y_mod3<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[5]*X$AgeSex+rep(Xe_comb$Epsilon,times=L) 
y_mod3

###############################   step 2: compute equations (7) and (8)  ###########################################################################
Var.Z_DM_mod3<-diag(data.matrix(X[,c(1:3,5)]) %*% solve(t(data.matrix(X[,c(1:3,5)]))%*%data.matrix(X[,c(1:3,5)])) %*% t(data.matrix(X[,c(1:3,5)]))+(1/(2*(N-k_mod3)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod3

Var.Phi_DM_mod3<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod3
Var.Phi_DM_mod3

############################# step 3: generate S normative samples ##################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod3<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[5]*X$AgeSex+Epsilon_s 
#############################  step 4: estimate the model parameters  ###############################################################################
Y_hat_s_mod3<-Z_hat_s_mod3<-Phi_hat_s_mod3<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod3<-matrix(0,1,ncol = S)
Var.Z_hat_L_mod3<-Var.Phi_hat_L_mod3<-matrix(0,nrow=N,ncol = S)
Var.Z_hat_mod3<-Var.Phi_hat_mod3<-matrix(0,Q*length(Z),S)
Betas_hat_mod3<-matrix(0,nrow=k_mod3,ncol = S)

for(i in 1:S){
  
  model <-lm(y_s_mod3[,i] ~ Age+Sex+AgeSex, data=data.frame(X,y_s_mod3[,i])) 
  Betas_hat_mod3[,i]<-model$coefficients
  Y_hat_s_mod3[,i]<-model$fitted.values
  sigma_epsilon_hat_mod3[i]<-sigma(model)
  
#############################   step 5: estimate Z-scores, compute equations (9) and (10)  #######################################################
  Z_hat_s_mod3[,i]<-(y_mod3-Y_hat_s_mod3[,i])/sigma_epsilon_hat_mod3[i]
  Phi_hat_s_mod3[,i]<-pnorm(Z_hat_s_mod3[,i])*100
  Var.Z_hat_L_mod3[,i]<-diag(data.matrix(X[,c(1:3,5)]) %*% solve(t(data.matrix(X[,c(1:3,5)]))%*%data.matrix(X[,c(1:3,5)])) %*% t(data.matrix(X[,c(1:3,5)]))+(1/(2*(N-k_mod3)))* Z_hat_s_mod3[,i] %*% t(Z_hat_s_mod3[,i]))
  Var.Phi_hat_L_mod3[,i]<-(100*dnorm(Z_hat_s_mod3[,i]))^2*Var.Z_hat_L_mod3[,i]
  
}


##############################    step 6: compute Monte Carlo averages and true variances   ###################################################################
Var.Z_S_mod3<-apply(Z_hat_s_mod3,MARGIN=1,FUN=var)  
Var.Phi_S_mod3<-apply(Phi_hat_s_mod3,MARGIN=1,FUN=var)
Var.Z_mod3<-Var.Phi_mod3<-matrix(0,Q*length(Z),1)

# take average across L replications (Note that they are all equals)
for(i in 1:S){
  
  t<-0
  
  for(l in 1:length(Z)){
    for(j in 1:length(Sex)){
      for(w in 1:length(Age)){
        
        t<-t+1
        Var.Z_hat_mod3[t,i]<-mean(Var.Z_hat_L_mod3[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
        Var.Phi_hat_mod3[t,i]<-mean(Var.Phi_hat_L_mod3[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
 }}}}

t<-0 # take average across L replications
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
Var.Z_mod3[t]<-mean(Var.Z_S_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod3[t]<-mean(Var.Phi_S_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
}}}

# MC average and SD for equation (9)
Mean_Var.Z_hat_mod3<-apply(Var.Z_hat_mod3,MARGIN=1,FUN=mean) 
SD_Var.Z_hat_mod3<-apply(Var.Z_hat_mod3,MARGIN=1,FUN=sd)  
summary(SD_Var.Z_hat_mod3)
# MC average and SD for equation (10)
Mean_Var.Phi_hat_mod3<-apply(Var.Phi_hat_mod3,MARGIN=1,FUN=mean)
SD_Var.Phi_hat_mod3<-apply(Var.Phi_hat_mod3,MARGIN=1,FUN=sd)
summary(SD_Var.Phi_hat_mod3)


###############################     step 7: compute Bias and coverage rate   #####################################################################
# R.B. for V(Z_hat)_hat
RB_Var.Z_hat_wrt_TV_mod3_L_1<-(Mean_Var.Z_hat_mod3-Var.Z_mod3)/Var.Z_mod3
summary(RB_Var.Z_hat_wrt_TV_mod3_L_1)
# R.B. for V(PR(Z_hat))_hat
RB_Var.PR_hat_wrt_TV_mod3_L_1<-(Mean_Var.Phi_hat_mod3-Var.Phi_mod3)/Var.Phi_mod3
summary(RB_Var.PR_hat_wrt_TV_mod3_L_1)  

## coverage rate
rm(Mean_Var.Z_hat_mod3, SD_Var.Z_hat_mod3, Mean_Var.Phi_hat_mod3, SD_Var.Phi_hat_mod3,Y_hat_s_mod3,sigma_epsilon_hat_mod3,y_s_mod3,
   Var.Z_DM_mod3_L, Var.Phi_DM_mod3_L, 
   MC.CV_Var.Phi_hat_mod3_L_1, MC.CV_Var.Z_hat_mod3_L_1, Epsilon_s, RB_Var.Z_hat_mod3_L_1, RB_Var.Phi_hat_mod3_L_1)

# Z score 
InfL_Z<-Z_hat_s_mod3-qnorm(0.975)*sqrt(Var.Z_hat_L_mod3) # lower-bound
SupL_Z<-Z_hat_s_mod3+qnorm(0.975)*sqrt(Var.Z_hat_L_mod3) # upper-bound
# PR score
InfL_PR<-Phi_hat_s_mod3-qnorm(0.975)*sqrt(Var.Phi_hat_L_mod3) # lower-bound
SupL_PR<-Phi_hat_s_mod3+qnorm(0.975)*sqrt(Var.Phi_hat_L_mod3) # upper-bound

cbind(rep(Xe_comb$Epsilon,times=L),rep(Xe_comb$Z,times=L))
Ztrue<-rep(Xe_comb$Z,times=L) # it should be like # rep(Xe_comb$Epsilon,times=L)
PRtrue<-pnorm(Ztrue)*100

# indicator variable: 1 if true parameter is within the confidence interval, 0 otherwilse
Indicator_L_Z<-Indicator_L_PR<-matrix(0,nrow=N,ncol = S)

for(j in 1:N){for(i in 1:S){
  ifelse(Ztrue[j]>=InfL_Z[j,i] & Ztrue[j]<=SupL_Z[j,i], Indicator_L_Z[j,i]<-1, Indicator_L_Z[j,i]<-0)
  ifelse(PRtrue[j]>=InfL_PR[j,i] & PRtrue[j]<=SupL_PR[j,i], Indicator_L_PR[j,i]<-1, Indicator_L_PR[j,i]<-0)
}}

# average of the indicator variable across S normative sample
coverage_rate_Z_S<-apply(Indicator_L_Z,MARGIN=1,FUN=mean)
coverage_rate_PR_S<-apply(Indicator_L_PR,MARGIN=1,FUN=mean)
# average of the indicator variable across L replications 
coverage_rate_Z_mod3_L_1<-coverage_rate_PR_mod3_L_1<-matrix(0,Q*length(Z),1)
t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      coverage_rate_Z_mod3_L_1[t]<-mean(coverage_rate_Z_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      coverage_rate_PR_mod3_L_1[t]<-mean(coverage_rate_PR_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
    }}}
#coverage rate
coverage_rate_Z_mod3_L_1;summary(coverage_rate_Z_mod3_L_1)
coverage_rate_PR_mod3_L_1;summary(coverage_rate_PR_mod3_L_1)

### Z-scores
## Figure S.B.1.14 Relative Bias of Equation (9)
# Sex=0
FigureS.B.1.14_Sex0<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                               RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.14_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.14_Sex0
# Sex=1
FigureS.B.1.14_Sex1<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                               RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.14_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.14_Sex1

## Figure S.B.1.17 coverage rate 
# Sex=0
FigureS.B.1.17_Sex0<-round(cbind(AGE,coverage_rate_Z_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                               coverage_rate_Z_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.17_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.17_Sex0
# Sex=1
FigureS.B.1.17_Sex1<-round(cbind(AGE,coverage_rate_Z_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                               coverage_rate_Z_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_Z_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.17_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.17_Sex1

### PR-scores # recall that for PR-scores L=2,5,6,10 and N=676,1690,2028,3380
## Figure S.B.1.37  Relative Bias of Equation (10)
# Sex=0
FigureS.B.1.37_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.37_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.37_Sex0
# Sex=1
FigureS.B.1.37_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.37_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.37_Sex1

## Figure S.B.1.42 coverage rate
# Sex=0
FigureS.B.1.42_Sex0<-round(cbind(AGE,coverage_rate_PR_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.42_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.42_Sex0
# Sex=1
FigureS.B.1.42_Sex1<-round(cbind(AGE,coverage_rate_PR_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_PR_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.42_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.42_Sex1
########################################################################################################################################################
############################################# model 4 #################################################################################################
#######################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)

##############################  step 1: generate the ys ###############################################################################################
y_mod4<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+rep(Xe_comb$Epsilon,times=L) 
y_mod4

###############################   step 2: compute equations (7) and (8)  ##############################################################################
Var.Z_DM_mod4<-diag(data.matrix(X[,-6]) %*% solve(t(data.matrix(X[,-6]))%*%data.matrix(X[,-6])) %*% t(data.matrix(X[,-6]))+(1/(2*(N-k_mod4)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod4

Var.Phi_DM_mod4<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod4
Var.Phi_DM_mod4

############################# step 3: generate S normative samples #####################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod4<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+Epsilon_s 
#############################  step 4: estimate the model parameters ####################################################################################
Y_hat_s_mod4<-Z_hat_s_mod4<-Phi_hat_s_mod4<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod4<-matrix(0,1,ncol = S)
Var.Z_hat_L_mod4<-Var.Phi_hat_L_mod4<-matrix(0,nrow=N,ncol = S)
Var.Z_hat_mod4<-Var.Phi_hat_mod4<-matrix(0,Q*length(Z),S)
Betas_hat_mod4<-matrix(0,nrow=k_mod4,ncol = S)


for(i in 1:S){
  
  model <-lm(y_s_mod4[,i] ~ Age+Sex+Age2+AgeSex, data=data.frame(X,y_s_mod4[,i])) 
  Betas_hat_mod4[,i]<-model$coefficients
  Y_hat_s_mod4[,i]<-model$fitted.values
  sigma_epsilon_hat_mod4[i]<-sigma(model)
  
#############################   step 5: estimate Z-scores, compute equations (9) and (10)   ##########################################################
  Z_hat_s_mod4[,i]<-(y_mod4-Y_hat_s_mod4[,i])/sigma_epsilon_hat_mod4[i]
  Phi_hat_s_mod4[,i]<-pnorm(Z_hat_s_mod4[,i])*100
  Var.Z_hat_L_mod4[,i]<-diag(data.matrix(X[,-6]) %*% solve(t(data.matrix(X[,-6]))%*%data.matrix(X[,-6])) %*% t(data.matrix(X[,-6]))+(1/(2*(N-k_mod4)))* Z_hat_s_mod4[,i] %*% t(Z_hat_s_mod4[,i]))
  Var.Phi_hat_L_mod4[,i]<-(100*dnorm(Z_hat_s_mod4[,i]))^2*Var.Z_hat_L_mod4[,i]
  
}



##############################    step 6: compute Monte Carlo averages and true variances   ###################################################################
Var.Z_S_mod4<-apply(Z_hat_s_mod4,MARGIN=1,FUN=var)  
Var.Phi_S_mod4<-apply(Phi_hat_s_mod4,MARGIN=1,FUN=var)
Var.Z_mod4<-Var.Phi_mod4<-matrix(0,Q*length(Z),1)

# take average across L replications (Note that they are all equals)
for(i in 1:S){
  
t<-0
  
for(l in 1:length(Z)){
for(j in 1:length(Sex)){
for(w in 1:length(Age)){
        
t<-t+1
Var.Z_hat_mod4[t,i]<-mean(Var.Z_hat_L_mod4[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
Var.Phi_hat_mod4[t,i]<-mean(Var.Phi_hat_L_mod4[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
}}}}

t<-0 # take average across L replications
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
Var.Z_mod4[t]<-mean(Var.Z_S_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod4[t]<-mean(Var.Phi_S_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
 }}}

# MC average and SD for equation (9)
Mean_Var.Z_hat_mod4<-apply(Var.Z_hat_mod4,MARGIN=1,FUN=mean) 
SD_Var.Z_hat_mod4<-apply(Var.Z_hat_mod4,MARGIN=1,FUN=sd)  
summary(SD_Var.Z_hat_mod4)
# MC average and SD for equation (10)
Mean_Var.Phi_hat_mod4<-apply(Var.Phi_hat_mod4,MARGIN=1,FUN=mean)
SD_Var.Phi_hat_mod4<-apply(Var.Phi_hat_mod4,MARGIN=1,FUN=sd)
summary(SD_Var.Phi_hat_mod4)


###############################     step 7: compute Bias and coverage rate   #####################################################################
# R.B. for V(Z_hat)_hat
RB_Var.Z_hat_wrt_TV_mod4_L_1<-(Mean_Var.Z_hat_mod4-Var.Z_mod4)/Var.Z_mod4
summary(RB_Var.Z_hat_wrt_TV_mod4_L_1)
# R.B. for V(PR(Z_hat))_hat
RB_Var.PR_hat_wrt_TV_mod4_L_1<-(Mean_Var.Phi_hat_mod4-Var.Phi_mod4)/Var.Phi_mod4
summary(RB_Var.PR_hat_wrt_TV_mod4_L_1)  

## coverage rate
rm(Mean_Var.Z_hat_mod4, SD_Var.Z_hat_mod4, Mean_Var.Phi_hat_mod4, SD_Var.Phi_hat_mod4,Y_hat_s_mod4,sigma_epsilon_hat_mod4,y_s_mod4,
   Var.Z_DM_mod4_L, Var.Phi_DM_mod4_L, 
   MC.CV_Var.Phi_hat_mod4_L_1, MC.CV_Var.Z_hat_mod4_L_1, Epsilon_s, RB_Var.Z_hat_mod4_L_1, RB_Var.Phi_hat_mod4_L_1)

# Z score 
InfL_Z<-Z_hat_s_mod4-qnorm(0.975)*sqrt(Var.Z_hat_L_mod4) # lower-bound
SupL_Z<-Z_hat_s_mod4+qnorm(0.975)*sqrt(Var.Z_hat_L_mod4) # upper-bound
# PR score
InfL_PR<-Phi_hat_s_mod4-qnorm(0.975)*sqrt(Var.Phi_hat_L_mod4) # lower-bound
SupL_PR<-Phi_hat_s_mod4+qnorm(0.975)*sqrt(Var.Phi_hat_L_mod4) # upper-bound

cbind(rep(Xe_comb$Epsilon,times=L),rep(Xe_comb$Z,times=L))
Ztrue<-rep(Xe_comb$Z,times=L) # it should be like # rep(Xe_comb$Epsilon,times=L)
PRtrue<-pnorm(Ztrue)*100

# indicator variable: 1 if true parameter is within the confidence interval, 0 otherwilse
Indicator_L_Z<-Indicator_L_PR<-matrix(0,nrow=N,ncol = S)

for(j in 1:N){for(i in 1:S){
  ifelse(Ztrue[j]>=InfL_Z[j,i] & Ztrue[j]<=SupL_Z[j,i], Indicator_L_Z[j,i]<-1, Indicator_L_Z[j,i]<-0)
  ifelse(PRtrue[j]>=InfL_PR[j,i] & PRtrue[j]<=SupL_PR[j,i], Indicator_L_PR[j,i]<-1, Indicator_L_PR[j,i]<-0)
}}

# average of the indicator variable across S normative sample
coverage_rate_Z_S<-apply(Indicator_L_Z,MARGIN=1,FUN=mean)
coverage_rate_PR_S<-apply(Indicator_L_PR,MARGIN=1,FUN=mean)
# average of the indicator variable across L replications 
coverage_rate_Z_mod4_L_1<-coverage_rate_PR_mod4_L_1<-matrix(0,Q*length(Z),1)
t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      coverage_rate_Z_mod4_L_1[t]<-mean(coverage_rate_Z_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      coverage_rate_PR_mod4_L_1[t]<-mean(coverage_rate_PR_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
    }}}
#coverage rate
coverage_rate_Z_mod4_L_1;summary(coverage_rate_Z_mod4_L_1)
coverage_rate_PR_mod4_L_1;summary(coverage_rate_PR_mod4_L_1)


### Z-scores
## Figure S.B.1.15 Relative Bias of Equation (9)
# Sex=0
FigureS.B.1.15_Sex0<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                               RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.15_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.15_Sex0
# Sex=1
FigureS.B.1.15_Sex1<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                               RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.15_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.15_Sex1

## Figure S.B.1.18 coverage rate 
# Sex=0
FigureS.B.1.18_Sex0<-round(cbind(AGE,coverage_rate_Z_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                               coverage_rate_Z_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.18_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.18_Sex0
# Sex=1
FigureS.B.1.18_Sex1<-round(cbind(AGE,coverage_rate_Z_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                               coverage_rate_Z_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_Z_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.18_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.18_Sex1

### PR-scores # recall that for PR-scores L=2,5,6,10 and N=676,1690,2028,3380
## Figure S.B.1.38  Relative Bias of Equation (10)
# Sex=0
FigureS.B.1.38_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.38_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.38_Sex0
# Sex=1
FigureS.B.1.38_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.38_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.38_Sex1

## Figure S.B.1.43 coverage rate
# Sex=0
FigureS.B.1.43_Sex0<-round(cbind(AGE,coverage_rate_PR_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.43_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.43_Sex0
# Sex=1
FigureS.B.1.43_Sex1<-round(cbind(AGE,coverage_rate_PR_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_PR_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.43_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.43_Sex1

#############################################################################################################################################################
############################################# model 5 #######################################################################################################
#############################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)

##############################  step 1: generate the ys #####################################################################################################
y_mod5<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+Beta[6]*X$Age2Sex+rep(Xe_comb$Epsilon,times=L)
y_mod5
###############################   step 2: compute equations (7) and (8)  ###################################################################################
Var.Z_DM_mod5<-diag(data.matrix(X) %*% solve(t(data.matrix(X))%*%data.matrix(X)) %*% t(data.matrix(X))+(1/(2*(N-k_mod5)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod5

Var.Phi_DM_mod5<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod5
Var.Phi_DM_mod5

############################# step 3: generate S normative samples ##########################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod5<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+Beta[6]*X$Age2Sex+Epsilon_s 

#############################  step 4: estimate the model parameters  ########################################################################################
Y_hat_s_mod5<-Z_hat_s_mod5<-Phi_hat_s_mod5<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod5<-matrix(0,1,ncol = S)
Var.Z_hat_L_mod5<-Var.Phi_hat_L_mod5<-matrix(0,nrow=N,ncol = S)
Var.Z_hat_mod5<-Var.Phi_hat_mod5<-matrix(0,Q*length(Z),S)
Betas_hat_mod5<-matrix(0,nrow=k_mod5,ncol = S)

for(i in 1:S){
  
  model <-lm(y_s_mod5[,i] ~ Age+Sex+Age2+AgeSex+Age2Sex, data=data.frame(X,y_s_mod5[,i])) 
  Betas_hat_mod5[,i]<-model$coefficients
    Y_hat_s_mod5[,i]<-model$fitted.values
  sigma_epsilon_hat_mod5[i]<-sigma(model)
  
#############################   step 5: estimate Z-scores, compute equations (9) and (10)   ##############################################################
  Z_hat_s_mod5[,i]<-(y_mod5-Y_hat_s_mod5[,i])/sigma_epsilon_hat_mod5[i]
  Phi_hat_s_mod5[,i]<-pnorm(Z_hat_s_mod5[,i])*100
  Var.Z_hat_L_mod5[,i]<-diag(data.matrix(X) %*% solve(t(data.matrix(X))%*%data.matrix(X)) %*% t(data.matrix(X))+(1/(2*(N-k_mod5)))* Z_hat_s_mod5[,i] %*% t(Z_hat_s_mod5[,i]))
  Var.Phi_hat_L_mod5[,i]<-(100*dnorm(Z_hat_s_mod5[,i]))^2*Var.Z_hat_L_mod5[,i]
  
}

##############################    step 6: compute Monte Carlo averages and true variances   ###################################################################
Var.Z_S_mod5<-apply(Z_hat_s_mod5,MARGIN=1,FUN=var)  
Var.Phi_S_mod5<-apply(Phi_hat_s_mod5,MARGIN=1,FUN=var)
Var.Z_mod5<-Var.Phi_mod5<-matrix(0,Q*length(Z),1)

# take average across L replications (Note that they are all equals)
for(i in 1:S){
  
  t<-0
  
  for(l in 1:length(Z)){
    for(j in 1:length(Sex)){
      for(w in 1:length(Age)){
        
        t<-t+1
        Var.Z_hat_mod5[t,i]<-mean(Var.Z_hat_L_mod5[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
        Var.Phi_hat_mod5[t,i]<-mean(Var.Phi_hat_L_mod5[X$Age==Age[w] & X$Sex==Sex[j] & Xe_comb$Z==Z[l],i])
  
 }}}}

t<-0 # take average across L replications
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
Var.Z_mod5[t]<-mean(Var.Z_S_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod5[t]<-mean(Var.Phi_S_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
 }}}

# MC average and SD for equation (9)
Mean_Var.Z_hat_mod5<-apply(Var.Z_hat_mod5,MARGIN=1,FUN=mean) 
SD_Var.Z_hat_mod5<-apply(Var.Z_hat_mod5,MARGIN=1,FUN=sd)  
summary(SD_Var.Z_hat_mod5)
# MC average and SD for equation (10)
Mean_Var.Phi_hat_mod5<-apply(Var.Phi_hat_mod5,MARGIN=1,FUN=mean)
SD_Var.Phi_hat_mod5<-apply(Var.Phi_hat_mod5,MARGIN=1,FUN=sd)
summary(SD_Var.Phi_hat_mod5)


###############################     step 7: compute Bias and coverage rate   #####################################################################
# R.B. for V(Z_hat)_hat
RB_Var.Z_hat_wrt_TV_mod5_L_1<-(Mean_Var.Z_hat_mod5-Var.Z_mod5)/Var.Z_mod5
summary(RB_Var.Z_hat_wrt_TV_mod5_L_1)
# R.B. for V(PR(Z_hat))_hat
RB_Var.PR_hat_wrt_TV_mod5_L_1<-(Mean_Var.Phi_hat_mod5-Var.Phi_mod5)/Var.Phi_mod5
summary(RB_Var.PR_hat_wrt_TV_mod5_L_1)  

## coverage rate
rm(Mean_Var.Z_hat_mod5, SD_Var.Z_hat_mod5, Mean_Var.Phi_hat_mod5, SD_Var.Phi_hat_mod5,Y_hat_s_mod5,sigma_epsilon_hat_mod5,y_s_mod5,
   Var.Z_DM_mod5_L, Var.Phi_DM_mod5_L,  Epsilon_s)

# Z score 
InfL_Z<-Z_hat_s_mod5-qnorm(0.975)*sqrt(Var.Z_hat_L_mod5) # lower-bound
SupL_Z<-Z_hat_s_mod5+qnorm(0.975)*sqrt(Var.Z_hat_L_mod5) # upper-bound
# PR score
InfL_PR<-Phi_hat_s_mod5-qnorm(0.975)*sqrt(Var.Phi_hat_L_mod5) # lower-bound
SupL_PR<-Phi_hat_s_mod5+qnorm(0.975)*sqrt(Var.Phi_hat_L_mod5) # upper-bound

cbind(rep(Xe_comb$Epsilon,times=L),rep(Xe_comb$Z,times=L))
Ztrue<-rep(Xe_comb$Z,times=L) # it should be like # rep(Xe_comb$Epsilon,times=L)
PRtrue<-pnorm(Ztrue)*100

# indicator variable: 1 if true parameter is within the confidence interval, 0 otherwilse
Indicator_L_Z<-Indicator_L_PR<-matrix(0,nrow=N,ncol = S)

for(j in 1:N){for(i in 1:S){
  ifelse(Ztrue[j]>=InfL_Z[j,i] & Ztrue[j]<=SupL_Z[j,i], Indicator_L_Z[j,i]<-1, Indicator_L_Z[j,i]<-0)
  ifelse(PRtrue[j]>=InfL_PR[j,i] & PRtrue[j]<=SupL_PR[j,i], Indicator_L_PR[j,i]<-1, Indicator_L_PR[j,i]<-0)
}}

# average of the indicator variable across S normative sample
coverage_rate_Z_S<-apply(Indicator_L_Z,MARGIN=1,FUN=mean)
coverage_rate_PR_S<-apply(Indicator_L_PR,MARGIN=1,FUN=mean)
# average of the indicator variable across L replications 
coverage_rate_Z_mod5_L_1<-coverage_rate_PR_mod5_L_1<-matrix(0,Q*length(Z),1)
t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      coverage_rate_Z_mod5_L_1[t]<-mean(coverage_rate_Z_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      coverage_rate_PR_mod5_L_1[t]<-mean(coverage_rate_PR_S[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
    }}}
#coverage rate
coverage_rate_Z_mod5_L_1;summary(coverage_rate_Z_mod5_L_1)
coverage_rate_PR_mod5_L_1;summary(coverage_rate_PR_mod5_L_1)


### Z-scores
## Figure S.A.12 Relative Bias of Equation (9)
# Sex=0
FigureS.A.12_Sex0<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                               RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.12_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.12_Sex0
# Sex=1
FigureS.A.12_Sex1<-round(cbind(AGE,RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                               RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.Z_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.12_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.12_Sex1

## Figure S.A.14 coverage rate 
# Sex=0
FigureS.A.14_Sex0<-round(cbind(AGE,coverage_rate_Z_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                               coverage_rate_Z_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.14_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.14_Sex0
# Sex=1
FigureS.A.14_Sex1<-round(cbind(AGE,coverage_rate_Z_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                               coverage_rate_Z_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_Z_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.14_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.14_Sex1

### PR-scores # recall that for PR-scores L=2,5,6,10 and N=676,1690,2028,3380
## Figure S.B.1.39  Relative Bias of Equation (10)
# Sex=0
FigureS.B.1.39_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.39_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.39_Sex0
# Sex=1
FigureS.B.1.39_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var.PR_hat_wrt_TV_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.39_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.39_Sex1

## Figure S.B.1.44 coverage rate
# Sex=0
FigureS.B.1.44_Sex0<-round(cbind(AGE,coverage_rate_PR_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.44_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.44_Sex0
# Sex=1
FigureS.B.1.44_Sex1<-round(cbind(AGE,coverage_rate_PR_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 coverage_rate_PR_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],coverage_rate_PR_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.44_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.44_Sex1


################################################################################################################################################################
################################## results for PR(Z)=1, 2.5, 5, 10, 90,95,97.5,99 #############################################################################
################################################################################################################################################################
rm(list=ls(all=T))
set.seed(05092016)

# Re-run steps 0-5 of the chosen model (this is needed to get the model parameters estimates)

Z_additional<-round(qnorm(p=c(0.01,0.025,0.05,0.1,0.9,0.95,0.975,0.99)),4)
Z_additional
### compute Epsilon_0
Epsilon_additional<-sigma_epsilon*Z_additional
Epsilon_additional
### compute Y_0
Xe_comb_additional<-data.frame(rep(X_comb$Age, times=T),rep(X_comb$Sex, times=T),rep(Epsilon_additional,each=Q),rep(Z_additional,each=Q),rep(c(0.01,0.025,0.05,0.1,0.9,0.95,0.975,0.99)*100,each=Q))
colnames(Xe_comb_additional)<-c("Age","Sex","Epsilon","Z","Phi")
Xe_comb_additional

X_additional<-data.frame(1,rep(Xe_comb_additional$Age,times=L),rep(Xe_comb_additional$Sex, times=L), rep(Age_2, times=length(Z_additional)*L), rep(Age_Sex, times=length(Z_additional)*L), rep(Age_2_Sex, times=length(Z_additional)*L) )
colnames(X_additional)<-c("Int","Age","Sex","Age2","AgeSex","Age2Sex")
X_additional

# model (1)
y_mod1_additional<-Beta[1]+Beta[2]*X_additional$Age+Beta[3]*X_additional$Sex+rep(Xe_comb_additional$Epsilon,times=L)
y_mod1_additional
# model (2)
#y_mod2_additional<-Beta[1]+Beta[2]*X_additional$Age+Beta[3]*X_additional$Sex+Beta[4]*X_additional$Age2+rep(Xe_comb_additional$Epsilon,times=L) 
#y_mod2_additional
# model (3)
#y_mod3_additional<-Beta[1]+Beta[2]*X_additional$Age+Beta[3]*X_additional$Sex+Beta[5]*X_additional$AgeSex+rep(Xe_comb_additional$Epsilon,times=L)
#y_mod3_additional
# model (4)
#y_mod4_additional<-Beta[1]+Beta[2]*X_additional$Age+Beta[3]*X_additional$Sex+Beta[4]*X_additional$Age2+Beta[5]*X_additional$AgeSex+rep(Xe_comb_additional$Epsilon,times=L) 
#y_mod4_additional
# model (5)
#y_mod5_additional<-Beta[1]+Beta[2]*X_additional$Age+Beta[3]*X_additional$Sex+Beta[4]*X_additional$Age2+Beta[5]*X_additional$AgeSex+Beta[6]*X_additional$Age2Sex+rep(Xe_comb_additional$Epsilon,times=L)
#y_mod5_additional

# compute Z, PR, V(Z)_hat, V(PR)_hat
Z_hat_s_mod1_additional<-Phi_hat_s_mod1_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)
Var.Z_DM_mod1_additional<-Var.Phi_DM_mod1_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)

for(i in 1:S){
  
  # model (1)
  Z_hat_s_mod1_additional[,i]<-(y_mod1_additional-(data.matrix(X_additional[,1:3])%*%Betas_hat_mod1[,i]))/sigma_epsilon_hat_mod1[i]
  Phi_hat_s_mod1_additional[,i]<-pnorm(Z_hat_s_mod1_additional[,i])*100
  # model (2)
  #Z_hat_s_mod2_additional[,i]<-(y_mod2_additional-(data.matrix(X_additional[,1:4])%*%Betas_hat_mod2[,i]))/sigma_epsilon_hat_mod2[i]
  #Phi_hat_s_mod2_additional[,i]<-pnorm(Z_hat_s_mod2_additional[,i])*100
  # model (3)
  #Z_hat_s_mod3_additional[,i]<-(y_mod3_additional-(data.matrix(X_additional[,c(1:3,5)])%*%Betas_hat_mod3[,i]))/sigma_epsilon_hat_mod3[i]
  #Phi_hat_s_mod3_additional[,i]<-pnorm(Z_hat_s_mod3_additional[,i])*100
  # model (4)
  #Z_hat_s_mod4_additional[,i]<-(y_mod4_additional-(data.matrix(X_additional[,-6])%*%Betas_hat_mod4[,i]))/sigma_epsilon_hat_mod4[i]
  #Phi_hat_s_mod4_additional[,i]<-pnorm(Z_hat_s_mod4_additional[,i])*100
  # model (5)
  #Z_hat_s_mod5_additional[,i]<-(y_mod5_additional-(data.matrix(X_additional)%*%Betas_hat_mod5[,i]))/sigma_epsilon_hat_mod5[i]
  #Phi_hat_s_mod5_additional[,i]<-pnorm(Z_hat_s_mod5_additional[,i])*100
  
  # model (1)
  Var.Z_DM_mod1_additional[,i]<-diag(data.matrix(X_additional[,1:3]) %*% solve(t(data.matrix(X[,1:3]))%*%data.matrix(X[,1:3])) %*% t(data.matrix(X_additional[,1:3]))+(1/(2*(N-k_mod1)))* Z_hat_s_mod1_additional[,i] %*% t(Z_hat_s_mod1_additional[,i]))
  Var.Phi_DM_mod1_additional[,i]<-(100*dnorm(Z_hat_s_mod1_additional[,i]))^2*Var.Z_DM_mod1_additional[,i]
  # model (2)
  #  Var.Z_DM_mod2_additional[,i]<-diag(data.matrix(X_additional[,1:4]) %*% solve(t(data.matrix(X[,1:4]))%*%data.matrix(X[,1:4])) %*% t(data.matrix(X_additional[,1:4]))+(1/(2*(N-k_mod2)))* Z_hat_s_mod2_additional[,i] %*% t(Z_hat_s_mod2_additional[,i]))
  #  Var.Phi_DM_mod2_additional[,i]<-(100*dnorm(Z_hat_s_mod2_additional[,i]))^2*Var.Z_DM_mod2_additional[,i]
  # model (3)
  #  Var.Z_DM_mod3_additional[,i]<-diag(data.matrix(X_additional[,c(1:3,5)]) %*% solve(t(data.matrix(X[,c(1:3,5)]))%*%data.matrix(X[,c(1:3,5)])) %*% t(data.matrix(X_additional[,c(1:3,5)]))+(1/(2*(N-k_mod3)))* Z_hat_s_mod3_additional[,i] %*% t(Z_hat_s_mod3_additional[,i]))
  #  Var.Phi_DM_mod3_additional[,i]<-(100*dnorm(Z_hat_s_mod3_additional[,i]))^2*Var.Z_DM_mod3_additional[,i]
  # model (4)
  #  Var.Z_DM_mod4_additional[,i]<-diag(data.matrix(X_additional[,-6]) %*% solve(t(data.matrix(X[,-6]))%*%data.matrix(X[,-6])) %*% t(data.matrix(X_additional[,-6]))+(1/(2*(N-k_mod4)))* Z_hat_s_mod4_additional[,i] %*% t(Z_hat_s_mod4_additional[,i]))
  #  Var.Phi_DM_mod4_additional[,i]<-(100*dnorm(Z_hat_s_mod4_additional[,i]))^2*Var.Z_DM_mod4_additional[,i]
  # model (5)
  #  Var.Z_DM_mod5_additional[,i]<-diag(data.matrix(X_additional) %*% solve(t(data.matrix(X))%*%data.matrix(X)) %*% t(data.matrix(X_additional))+(1/(2*(N-k_mod5)))* Z_hat_s_mod5_additional[,i] %*% t(Z_hat_s_mod5_additional[,i])
  #  Var.Phi_DM_mod5_additional[,i]<-(100*dnorm(Z_hat_s_mod5_additional[,i]))^2*Var.Z_DM_mod5_additional[,i]
  
}


# compute true variances across S normative samples
Var.Z_S_mod1_additional<-apply(Z_hat_s_mod1_additional,MARGIN=1,FUN=var)  
Var.Phi_S_mod1_additional<-apply(Phi_hat_s_mod1_additional,MARGIN=1,FUN=var)
Var.Z_mod1_additional<-Var.Phi_mod1_additional<-matrix(0,Q*length(Z_additional),1)

t<-0 # take average across L replications
for(l in 1:length(Z_additional)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
t<-t+1
  # true variance
Var.Z_mod1_additional[t]<-mean(Var.Z_S_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
Var.Phi_mod1_additional[t]<-mean(Var.Phi_S_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      
    }}}

Var.Z_hat_mod1_additional<-Var.Phi_hat_mod1_additional<-matrix(0,Q*length(Z_additional),S)

for(i in 1:S){
  t<-0
  for(l in 1:length(Z_additional)){    for(j in 1:length(Sex)){      for(w in 1:length(Age)){
    t<-t+1
    # equations (9) and (10)
    Var.Z_hat_mod1_additional[t,i]<-mean(Var.Z_DM_mod1_additional[X_additional$Age==Age[w] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l],i])
    Var.Phi_hat_mod1_additional[t,i]<-mean(Var.Phi_DM_mod1_additional[X_additional$Age==Age[w] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l],i])
  }}}}

# MC average for equation (9)
Mean_Var.Z_hat_mod1_additional<-apply(Var.Z_hat_mod1_additional,MARGIN=1,FUN=mean) 
# MC average for equation (10)
Mean_Var.Phi_hat_mod1_additional<-apply(Var.Phi_hat_mod1_additional,MARGIN=1,FUN=mean)

# R.B. for V(Z_hat)_hat
RB_Var.Z_hat_wrt_TV_mod1_L_1_additional<-(Mean_Var.Z_hat_mod1_additional-Var.Z_mod1_additional)/Var.Z_mod1_additional
summary(RB_Var.Z_hat_wrt_TV_mod1_L_1_additional)
# R.B. for V(PR(Z_hat))_hat
RB_Var.PR_hat_wrt_TV_mod1_L_1_additional<-(Mean_Var.Phi_hat_mod1_additional-Var.Phi_mod1_additional)/Var.Phi_mod1_additional
summary(RB_Var.PR_hat_wrt_TV_mod1_L_1_additional)  

## coverage rate
rm(SD_Var.Z_hat_mod1, SD_Var.Phi_hat_mod1,Y_hat_s_mod1,sigma_epsilon_hat_mod1,y_s_mod1,
   Var.Z_DM_mod1_L, Var.Phi_DM_mod1_L, 
   MC.CV_Var.Phi_hat_mod1_L_1, MC.CV_Var.Z_hat_mod1_L_1, Epsilon_s, RB_Var.Z_hat_mod1_L_1, RB_Var.Phi_hat_mod1_L_1)

# Z score
InfL_Z_additional<-Z_hat_s_mod1_additional-qnorm(0.975)*sqrt(Var.Z_DM_mod1_additional) # lower-bound
SupL_Z_additional<-Z_hat_s_mod1_additional+qnorm(0.975)*sqrt(Var.Z_DM_mod1_additional) # upper-bound
# PR score
InfL_PR_additional<-Phi_hat_s_mod1_additional-qnorm(0.975)*sqrt(Var.Phi_DM_mod1_additional) # lower-bound
SupL_PR_additional<-Phi_hat_s_mod1_additional+qnorm(0.975)*sqrt(Var.Phi_DM_mod1_additional) # upper-bound

cbind(rep(Xe_comb_additional$Epsilon,times=L),rep(Xe_comb_additional$Z,times=L))
Ztrue_additional<-rep(Xe_comb_additional$Z,times=L) # it should be like # rep(Xe_comb$Epsilon,times=L)
PRtrue_additional<-pnorm(Ztrue_additional)*100

# indicator variable: 1 if true parameter is within the confidence interval, 0 otherwilse
Indicator_L_Z_additional<-Indicator_L_PR_additional<-matrix(0,nrow=Q*length(Z_additional)*L,ncol = S)

for(j in 1:(Q*length(Z_additional)*L)){for(i in 1:S){
  ifelse(Ztrue_additional[j]>=InfL_Z_additional[j,i] & Ztrue_additional[j]<=SupL_Z_additional[j,i], Indicator_L_Z_additional[j,i]<-1, Indicator_L_Z_additional[j,i]<-0)
  ifelse(PRtrue_additional[j]>=InfL_PR_additional[j,i] & PRtrue_additional[j]<=SupL_PR_additional[j,i], Indicator_L_PR_additional[j,i]<-1, Indicator_L_PR_additional[j,i]<-0)
}}

# average of the indicator variable across S normative sample
coverage_rate_Z_S_additional<-apply(Indicator_L_Z_additional,MARGIN=1,FUN=mean)
coverage_rate_PR_S_additional<-apply(Indicator_L_PR_additional,MARGIN=1,FUN=mean)
# average of the indicator variable across L replications 
coverage_rate_Z_mod1_L_1_additional<-coverage_rate_PR_mod1_L_1_additional<-matrix(0,Q*length(Z_additional),1)
t<-0
for(l in 1:length(Z_additional)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      coverage_rate_Z_mod1_L_1_additional[t]<-mean(coverage_rate_Z_S_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      coverage_rate_PR_mod1_L_1_additional[t]<-mean(coverage_rate_PR_S_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      
    }}}
#coverage rate
coverage_rate_Z_mod1_L_1_additional;summary(coverage_rate_Z_mod1_L_1_additional)
coverage_rate_PR_mod1_L_1_additional;summary(coverage_rate_PR_mod1_L_1_additional)


### PR-scores
################### model (1) (for other models replace "mod1" with "mod#ofmodel#)
## Figure S.A.15  Relative Bias of Equation (10)
# Sex=0
FigureS.A.15_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                              RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                              RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                              RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.15_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.15_Sex0
# Sex=1
FigureS.A.15_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                              RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                              RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                              RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],RB_Var.PR_hat_wrt_TV_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.15_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.15_Sex1

## Figure S.A.17 coverage rate
# Sex=0
FigureS.A.17_Sex0<-round(cbind(AGE,coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                               coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                               coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                               coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.17_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.17_Sex0
# Sex=1
FigureS.A.17_Sex1<-round(cbind(AGE,coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                               coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                               coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                               coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],coverage_rate_PR_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.17_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.17_Sex1

################### model (5)
## Figure S.A.16  Relative Bias of Equation (10)
# Sex=0
FigureS.A.16_Sex0<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                              RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                              RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                              RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.16_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.16_Sex0
# Sex=1
FigureS.A.16_Sex1<-round(cbind(AGE,RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                              RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                              RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                              RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],RB_Var.PR_hat_wrt_TV_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.16_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.16_Sex1

## Figure S.A.18 coverage rate
# Sex=0
FigureS.A.18_Sex0<-round(cbind(AGE,coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                               coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                               coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                               coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.18_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.18_Sex0
# Sex=1
FigureS.A.18_Sex1<-round(cbind(AGE,coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                               coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                               coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                               coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],coverage_rate_PR_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.18_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.18_Sex1
