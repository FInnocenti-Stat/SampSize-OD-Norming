#######################################################################################################################################
#### Sample Size Calculation and Optimal Design for Regression-Based Norming of Tests and Questionnaires
#######################################################################################################################################
# Francesco Innocenti, Frans E. S. Tan, Math J. J. M. Candel, and Gerard J. P. van Breukelen
#################################################################################################################################
rm(list=ls(all=T))
set.seed(05092016)
#################################################################################################################################
########################################## Simulation study for the variances ###################################################
#################################################################################################################################
##############################  step 0: set parameters and Z values ############################# 
#      Intercept Age    Sex   Age^2 AgexSex Age^2xSex   SD(epsilon)
Beta<-c(18.579,-0.054,0.753,-0.002, -0.03,-0.001); sigma_epsilon<-4.818  # PNVFT
#Beta<-c(-3.316,0.231,0.794,-0.004, 0.05, 0.001); sigma_epsilon<-2.605   # DKEFS
#Beta<-c(10.854,0.040,-0.192,-0.002, 0.05,0.004); sigma_epsilon<-2.658  # SDMT


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
L<-1 # 2, 5, 10
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
which(y_mod1==y_mod1[1])
length(which(y_mod1==y_mod1[1]))==L

###############################   step 2: compute equations (7) and (8)  #############################################################################
Var.Z_DM_mod1<-diag(data.matrix(X[,1:3]) %*% solve(t(data.matrix(X[,1:3]))%*%data.matrix(X[,1:3])) %*% t(data.matrix(X[,1:3]))+(1/(2*(N-k_mod1)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod1

Var.Phi_DM_mod1<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod1
Var.Phi_DM_mod1

############################# step 3: generate S normative samples ####################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod1<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Epsilon_s
#############################  step 4: estimate the model parameters ##################################################################################
Y_hat_s_mod1<-Epsilon_hat_S_mod1<-Z_hat_s_mod1<-Phi_hat_s_mod1<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod1<-matrix(0,1,ncol = S)
Betas_hat_mod1<-matrix(0,nrow=k_mod1,ncol = S)

for(i in 1:S){

model1 <-lm(y_s_mod1[,i] ~ Age+Sex, data=data.frame(X,y_s_mod1[,i])) 
Betas_hat_mod1[,i]<-model1$coefficients
Y_hat_s_mod1[,i]<-model1$fitted.values
Epsilon_hat_S_mod1[,i]<-model1$residuals
sigma_epsilon_hat_mod1[i]<-sigma(model1)  

#############################   step 5: estimate Z-scores and PR-scores   ##############################################################################
Z_hat_s_mod1[,i]<-(y_mod1-Y_hat_s_mod1[,i])/sigma_epsilon_hat_mod1[i]
Phi_hat_s_mod1[,i]<-pnorm(Z_hat_s_mod1[,i])*100
}

##############################    step 6: compute true values    #######################################################################################
Z_bar_S_mod1<-apply(Z_hat_s_mod1,MARGIN=1,FUN=mean) 
Var.Z_S_mod1<-apply(Z_hat_s_mod1,MARGIN=1,FUN=var)  

Phi_bar_S_mod1<-apply(Phi_hat_s_mod1,MARGIN=1,FUN=mean)
Var.Phi_S_mod1<-apply(Phi_hat_s_mod1,MARGIN=1,FUN=var)

# take average across L replications (Note that they are all equals)
Z_bar_mod1<-Var.Z_mod1<-Phi_bar_mod1<-Var.Phi_mod1<-matrix(0,Q*length(Z),1)
Var.Z_DM_mod1_L<-matrix(0,Q*length(Z),1)
Var.Phi_DM_mod1_L<-matrix(0,Q*length(Z),1)

t<-0
for(l in 1:length(Z)){
for(j in 1:length(Sex)){
for(i in 1:length(Age)){
      
t<-t+1
Z_bar_mod1[t]<-mean(Z_bar_S_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_mod1[t]<-mean(Var.Z_S_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_DM_mod1_L[t]<-mean(Var.Z_DM_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
Phi_bar_mod1[t]<-mean(Phi_bar_S_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod1[t]<-mean(Var.Phi_S_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_DM_mod1_L[t]<-mean(Var.Phi_DM_mod1[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])

}}}

###############################     step 7: compute Bias    ###############################################################################################
# Absolute Bias for Z_hat
AB_Z_hat_mod1_L_1<-(Z_bar_mod1-Xe_comb$Z)
summary(AB_Z_hat_mod1_L_1)  
# Absolute Bias for PR(Z_hat)
AB_Phi_hat_mod1_L_1<-(Phi_bar_mod1-Xe_comb$Phi)
summary(AB_Phi_hat_mod1_L_1)  
# R.B. for V(Z_hat)
RB_Var_Z_DM_mod1_L_1<-(Var.Z_DM_mod1_L-Var.Z_mod1)/Var.Z_mod1
RB_Var_Z_DM_mod1_L_1
summary(RB_Var_Z_DM_mod1_L_1)  
# R.B. for V(PR(Z_hat))
RB_Var_Phi_DM_mod1_L_1<-(Var.Phi_DM_mod1_L-Var.Phi_mod1)/Var.Phi_mod1
RB_Var_Phi_DM_mod1_L_1
summary(RB_Var_Phi_DM_mod1_L_1)  

### Z-scores
## Figure S.A.3 Relative Bias of Equation (7)
# Sex=0
FigureS.A.3_Sex0<-round(cbind(AGE,RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
            RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.3_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.3_Sex0
# Sex=1
FigureS.A.3_Sex1<-round(cbind(AGE,RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                              RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Z_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.3_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.3_Sex1

## Figure S.A.5 Absolute Bias of Z_hat
# Sex=0
FigureS.A.5_Sex0<-round(cbind(AGE,AB_Z_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Z_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Z_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                              AB_Z_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Z_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Z_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.5_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.5_Sex0
# Sex=1
FigureS.A.5_Sex1<-round(cbind(AGE,AB_Z_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Z_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Z_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                              AB_Z_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Z_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Z_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.5_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.5_Sex1

### PR-scores
## Figure S.B.1.25  Relative Bias of Equation (8)
# Sex=0
FigureS.B.1.25_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.25_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.25_Sex0
# Sex=1
FigureS.B.1.25_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Phi_DM_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.25_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.25_Sex1

## Figure S.B.1.30 Absolute Bias of Equation (6)
# Sex=0
FigureS.B.1.30_Sex0<-round(cbind(AGE,AB_Phi_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                              AB_Phi_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.30_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.30_Sex0
# Sex=1
FigureS.B.1.30_Sex1<-round(cbind(AGE,AB_Phi_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                              AB_Phi_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Phi_hat_mod1_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.30_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.30_Sex1

##########################################################################################################################################################
############################################# model 2 ####################################################################################################
##########################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)
##############################  step 1: generate the ys ##################################################################################################
y_mod2<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+rep(Xe_comb$Epsilon,times=L) 
y_mod2

###############################   step 2: compute equations (7) and (8)  #################################################################################
Var.Z_DM_mod2<-diag(data.matrix(X[,1:4]) %*% solve(t(data.matrix(X[,1:4]))%*%data.matrix(X[,1:4])) %*% t(data.matrix(X[,1:4]))+(1/(2*(N-k_mod2)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod2

Var.Phi_DM_mod2<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod2
Var.Phi_DM_mod2

############################# step 3: generate S normative samples #######################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod2<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Epsilon_s 

#############################  step 4: estimate the model parameters  #####################################################################################
Y_hat_s_mod2<-Epsilon_hat_S_mod2<-Z_hat_s_mod2<-Phi_hat_s_mod2<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod2<-matrix(0,1,ncol = S)
Betas_hat_mod2<-matrix(0,nrow=k_mod2,ncol = S)

for(i in 1:S){

model2 <-lm(y_s_mod2[,i] ~ Age+Sex+Age2, data=data.frame(X,y_s_mod2[,i])) 
Betas_hat_mod2[,i]<-model2$coefficients
Y_hat_s_mod2[,i]<-model2$fitted.values
Epsilon_hat_S_mod2[,i]<-model2$residuals
sigma_epsilon_hat_mod2[i]<-sigma(model2)

#############################   step 5: estimate Z scores and  PR scores   ###############################################################################
Z_hat_s_mod2[,i]<-(y_mod2-Y_hat_s_mod2[,i])/sigma_epsilon_hat_mod2[i]
Phi_hat_s_mod2[,i]<-pnorm(Z_hat_s_mod2[,i])*100

}

##############################    step 6: compute true values    #########################################################################################
Z_bar_S_mod2<-apply(Z_hat_s_mod2,MARGIN=1,FUN=mean)
Var.Z_S_mod2<-apply(Z_hat_s_mod2,MARGIN=1,FUN=var)

Phi_bar_S_mod2<-apply(Phi_hat_s_mod2,MARGIN=1,FUN=mean)
Var.Phi_S_mod2<-apply(Phi_hat_s_mod2,MARGIN=1,FUN=var)

# take average across L replications (Note that they are all equals)

Z_bar_mod2<-Var.Z_mod2<-Phi_bar_mod2<-Var.Phi_mod2<-matrix(0,Q*length(Z),1)
Var.Z_DM_mod2_L<-matrix(0,Q*length(Z),1)
Var.Phi_DM_mod2_L<-matrix(0,Q*length(Z),1)

t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
  
Z_bar_mod2[t]<-mean(Z_bar_S_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_mod2[t]<-mean(Var.Z_S_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_DM_mod2_L[t]<-mean(Var.Z_DM_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
      
Phi_bar_mod2[t]<-mean(Phi_bar_S_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod2[t]<-mean(Var.Phi_S_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_DM_mod2_L[t]<-mean(Var.Phi_DM_mod2[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
  }}}
###############################     step 7: compute  Bias    ###############################################################################################
# Absolute bias for Z_hat
AB_Z_hat_mod2_L_1<-(Z_bar_mod2-Xe_comb$Z)
summary(AB_Z_hat_mod2_L_1)  
# Absolute bias for PR(Z_hat)
AB_Phi_hat_mod2_L_1<-(Phi_bar_mod2-Xe_comb$Phi)
summary(AB_Phi_hat_mod2_L_1)  
# R.B. for V(Z_hat)
RB_Var_Z_DM_mod2_L_1<-(Var.Z_DM_mod2_L-Var.Z_mod2)/Var.Z_mod2
RB_Var_Z_DM_mod2_L_1
summary(RB_Var_Z_DM_mod2_L_1)  
# R.B. for V(PR(Z_hat))
RB_Var_Phi_DM_mod2_L_1<-(Var.Phi_DM_mod2_L-Var.Phi_mod2)/Var.Phi_mod2
RB_Var_Phi_DM_mod2_L_1
summary(RB_Var_Phi_DM_mod2_L_1)  

### Z-scores
## Figure S.B.1.1   Relative Bias of Equation (7)
# Sex=0
FigureS.B.1.1_Sex0<-round(cbind(AGE,RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                              RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.1_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.1_Sex0
# Sex=1
FigureS.B.1.1_Sex1<-round(cbind(AGE,RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                              RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Z_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.1_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.1_Sex1

## Figure S.B.1.4   Absolute Bias of Z_hat
# Sex=0
FigureS.B.1.4_Sex0<-round(cbind(AGE,AB_Z_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Z_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Z_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                              AB_Z_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Z_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Z_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.4_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.4_Sex0
# Sex=1
FigureS.B.1.4_Sex1<-round(cbind(AGE,AB_Z_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Z_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Z_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                              AB_Z_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Z_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Z_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.4_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.4_Sex1

### PR-scores
## Figure S.B.1.26   Relative Bias of Equation (8)
# Sex=0
FigureS.B.1.26_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.26_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.26_Sex0
# Sex=1
FigureS.B.1.26_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Phi_DM_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.26_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.26_Sex1

## Figure S.B.1.31   Absolute Bias of Equation (6)
# Sex=0
FigureS.B.1.31_Sex0<-round(cbind(AGE,AB_Phi_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.31_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.31_Sex0
# Sex=1
FigureS.B.1.31_Sex1<-round(cbind(AGE,AB_Phi_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Phi_hat_mod2_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.31_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.31_Sex1

############################################################################################################################################################
############################################# model 3 ######################################################################################################
############################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)
##############################  step 1: generate the ys ####################################################################################################
y_mod3<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[5]*X$AgeSex+rep(Xe_comb$Epsilon,times=L) 
y_mod3

###############################   step 2: compute equations (7) and (8)  ##################################################################################
Var.Z_DM_mod3<-diag(data.matrix(X[,c(1:3,5)]) %*% solve(t(data.matrix(X[,c(1:3,5)]))%*%data.matrix(X[,c(1:3,5)])) %*% t(data.matrix(X[,c(1:3,5)]))+(1/(2*(N-k_mod3)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod3

Var.Phi_DM_mod3<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod3
Var.Phi_DM_mod3

############################# step 3: generate S normative samples #########################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod3<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[5]*X$AgeSex+Epsilon_s 
#############################  step 4: estimate the model parameters ######################################################################################
Y_hat_s_mod3<-Epsilon_hat_S_mod3<-Z_hat_s_mod3<-Phi_hat_s_mod3<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod3<-matrix(0,1,ncol = S)
Betas_hat_mod3<-matrix(0,nrow=k_mod3,ncol = S)

for(i in 1:S){

model3 <-lm(y_s_mod3[,i] ~ Age+Sex+AgeSex, data=data.frame(X,y_s_mod3[,i])) 
Betas_hat_mod3[,i]<-model3$coefficients
Y_hat_s_mod3[,i]<-model3$fitted.values
Epsilon_hat_S_mod3[,i]<-model3$residuals
sigma_epsilon_hat_mod3[i]<-sigma(model3)

#############################   step 5: estimate Z scores and PR scores   ###################################################################################
Z_hat_s_mod3[,i]<-(y_mod3-Y_hat_s_mod3[,i])/sigma_epsilon_hat_mod3[i]
Phi_hat_s_mod3[,i]<-pnorm(Z_hat_s_mod3[,i])*100

}


##############################    step 6: compute true values    #############################################################################################
Z_bar_S_mod3<-apply(Z_hat_s_mod3,MARGIN=1,FUN=mean)
Var.Z_S_mod3<-apply(Z_hat_s_mod3,MARGIN=1,FUN=var)

Phi_bar_S_mod3<-apply(Phi_hat_s_mod3,MARGIN=1,FUN=mean)
Var.Phi_S_mod3<-apply(Phi_hat_s_mod3,MARGIN=1,FUN=var)

# take average across L replications (Note that they are all equals)

Z_bar_mod3<-Var.Z_mod3<-Phi_bar_mod3<-Var.Phi_mod3<-matrix(0,Q*length(Z),1)
Var.Z_DM_mod3_L<-matrix(0,Q*length(Z),1)
Var.Phi_DM_mod3_L<-matrix(0,Q*length(Z),1)


t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
t<-t+1
Z_bar_mod3[t]<-mean(Z_bar_S_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_mod3[t]<-mean(Var.Z_S_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_DM_mod3_L[t]<-mean(Var.Z_DM_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
Phi_bar_mod3[t]<-mean(Phi_bar_S_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod3[t]<-mean(Var.Phi_S_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_DM_mod3_L[t]<-mean(Var.Phi_DM_mod3[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
 }}}
###############################     step 7: compute  Bias    #################################################################################################
# Absolute bias for Z_hat
AB_Z_hat_mod3_L_1<-(Z_bar_mod3-Xe_comb$Z)
summary(AB_Z_hat_mod3_L_1)  
# Absolute bias for PR(Z_hat)
AB_Phi_hat_mod3_L_1<-(Phi_bar_mod3-Xe_comb$Phi)
summary(AB_Phi_hat_mod3_L_1)  
# R.B. for V(Z_hat)
RB_Var_Z_DM_mod3_L_1<-(Var.Z_DM_mod3_L-Var.Z_mod3)/Var.Z_mod3
RB_Var_Z_DM_mod3_L_1
summary(RB_Var_Z_DM_mod3_L_1)  
# R.B. for V(PR(Z_hat))
RB_Var_Phi_DM_mod3_L_1<-(Var.Phi_DM_mod3_L-Var.Phi_mod3)/Var.Phi_mod3
RB_Var_Phi_DM_mod3_L_1
summary(RB_Var_Phi_DM_mod3_L_1)  

### Z-scores
## Figure S.B.1.2   Relative Bias of Equation (7)
# Sex=0
FigureS.B.1.2_Sex0<-round(cbind(AGE,RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.2_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.2_Sex0
# Sex=1
FigureS.B.1.2_Sex1<-round(cbind(AGE,RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Z_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.2_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.2_Sex1

## Figure S.B.1.5   Absolute Bias of Z_hat
# Sex=0
FigureS.B.1.5_Sex0<-round(cbind(AGE,AB_Z_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Z_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Z_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                AB_Z_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Z_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Z_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.5_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.5_Sex0
# Sex=1
FigureS.B.1.5_Sex1<-round(cbind(AGE,AB_Z_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Z_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Z_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                AB_Z_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Z_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Z_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.5_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.5_Sex1

### PR-scores
## Figure S.B.1.27   Relative Bias of Equation (8)
# Sex=0
FigureS.B.1.27_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.27_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.27_Sex0
# Sex=1
FigureS.B.1.27_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Phi_DM_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.27_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.27_Sex1

## Figure S.B.1.32   Absolute Bias of Equation (6)
# Sex=0
FigureS.B.1.32_Sex0<-round(cbind(AGE,AB_Phi_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.32_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.32_Sex0
# Sex=1
FigureS.B.1.32_Sex1<-round(cbind(AGE,AB_Phi_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Phi_hat_mod3_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.32_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.32_Sex1

################################################################################################################################################################
############################################# model 4 #########################################################################################################
##############################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)
##############################  step 1: generate the ys #######################################################################################################
y_mod4<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+rep(Xe_comb$Epsilon,times=L) 
y_mod4

###############################   step 2: compute equations (7) and (8)  #####################################################################################
Var.Z_DM_mod4<-diag(data.matrix(X[,-6]) %*% solve(t(data.matrix(X[,-6]))%*%data.matrix(X[,-6])) %*% t(data.matrix(X[,-6]))+(1/(2*(N-k_mod4)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod4

Var.Phi_DM_mod4<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod4
Var.Phi_DM_mod4

############################# step 3: generate S normative samples ##########################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod4<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+Epsilon_s 
#############################  step 4: estimate the model parameters #########################################################################################
Y_hat_s_mod4<-Epsilon_hat_S_mod4<-Z_hat_s_mod4<-Phi_hat_s_mod4<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod4<-matrix(0,1,ncol = S)
Betas_hat_mod4<-matrix(0,nrow=k_mod4,ncol = S)

for(i in 1:S){

model4 <-lm(y_s_mod4[,i] ~ Age+Sex+Age2+AgeSex, data=data.frame(X,y_s_mod4[,i])) 
Betas_hat_mod4[,i]<-model4$coefficients
Y_hat_s_mod4[,i]<-model4$fitted.values
Epsilon_hat_S_mod4[,i]<-model4$residuals
sigma_epsilon_hat_mod4[i]<-sigma(model4)
 
#############################   step 5: estimate Z-scores and PR-scores   #####################################################################################
Z_hat_s_mod4[,i]<-(y_mod4-Y_hat_s_mod4[,i])/sigma_epsilon_hat_mod4[i]
Phi_hat_s_mod4[,i]<-pnorm(Z_hat_s_mod4[,i])*100
 
}

##############################    step 6: compute true values    #############################################################################################
Z_bar_S_mod4<-apply(Z_hat_s_mod4,MARGIN=1,FUN=mean)
Var.Z_S_mod4<-apply(Z_hat_s_mod4,MARGIN=1,FUN=var)

Phi_bar_S_mod4<-apply(Phi_hat_s_mod4,MARGIN=1,FUN=mean)
Var.Phi_S_mod4<-apply(Phi_hat_s_mod4,MARGIN=1,FUN=var)

# take average across L replications (Note that they are all equals)

Z_bar_mod4<-Var.Z_mod4<-Phi_bar_mod4<-Var.Phi_mod4<-matrix(0,Q*length(Z),1)
Var.Z_DM_mod4_L<-matrix(0,Q*length(Z),1)
Var.Phi_DM_mod4_L<-matrix(0,Q*length(Z),1)


t<-0
for(l in 1:length(Z)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
t<-t+1

Z_bar_mod4[t]<-mean(Z_bar_S_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_mod4[t]<-mean(Var.Z_S_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_DM_mod4_L[t]<-mean(Var.Z_DM_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
Phi_bar_mod4[t]<-mean(Phi_bar_S_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod4[t]<-mean(Var.Phi_S_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_DM_mod4_L[t]<-mean(Var.Phi_DM_mod4[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
    }}}
###############################     step 7: compute  Bias    ###################################################################################################
# Absolute Bias for Z_hat
AB_Z_hat_mod4_L_1<-(Z_bar_mod4-Xe_comb$Z)
summary(AB_Z_hat_mod4_L_1)  
# Absolute Bias for PR(Z_hat)
AB_Phi_hat_mod4_L_1<-(Phi_bar_mod4-Xe_comb$Phi)
summary(AB_Phi_hat_mod4_L_1)  
# R.B. for V(Z_hat)
RB_Var_Z_DM_mod4_L_1<-(Var.Z_DM_mod4_L-Var.Z_mod4)/Var.Z_mod4
RB_Var_Z_DM_mod4_L_1
summary(RB_Var_Z_DM_mod4_L_1)  
# R.B. for V(PR(Z_hat))
RB_Var_Phi_DM_mod4_L_1<-(Var.Phi_DM_mod4_L-Var.Phi_mod4)/Var.Phi_mod4
RB_Var_Phi_DM_mod4_L_1
summary(RB_Var_Phi_DM_mod4_L_1)  

### Z-scores
## Figure S.B.1.3   Relative Bias of Equation (7)
# Sex=0
FigureS.B.1.3_Sex0<-round(cbind(AGE,RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.3_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.3_Sex0
# Sex=1
FigureS.B.1.3_Sex1<-round(cbind(AGE,RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Z_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.3_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.3_Sex1

## Figure S.B.1.6   Absolute Bias of Z_hat
# Sex=0
FigureS.B.1.6_Sex0<-round(cbind(AGE,AB_Z_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Z_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Z_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                AB_Z_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Z_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Z_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.6_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.6_Sex0
# Sex=1
FigureS.B.1.6_Sex1<-round(cbind(AGE,AB_Z_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Z_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Z_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                AB_Z_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Z_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Z_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.6_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.B.1.6_Sex1

### PR-scores
## Figure S.B.1.28   Relative Bias of Equation (8)
# Sex=0
FigureS.B.1.28_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.28_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.28_Sex0
# Sex=1
FigureS.B.1.28_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Phi_DM_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.28_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.28_Sex1

## Figure S.B.1.33   Absolute Bias of Equation (6)
# Sex=0
FigureS.B.1.33_Sex0<-round(cbind(AGE,AB_Phi_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.33_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.33_Sex0
# Sex=1
FigureS.B.1.33_Sex1<-round(cbind(AGE,AB_Phi_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Phi_hat_mod4_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.33_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.33_Sex1

##############################################################################################################################################################
############################################# model 5 #######################################################################################################
##############################################################################################################################################################
rm(list= ls()[!(ls() %in% c('N','X','L','N1','M', 'Xe_comb','Q','X_comb','AGE','Age_2_Sex','Age_Sex','Age','Age_2','Sex','Phi','Epsilon','Z','Beta','sigma_epsilon','k_mod1','k_mod2','k_mod3','k_mod4', 'k_mod5'))])
set.seed(05092016)
##############################  step 1: generate the ys ######################################################################################################
y_mod5<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+Beta[6]*X$Age2Sex+rep(Xe_comb$Epsilon,times=L)
y_mod5
###############################   step 2: compute equations (7) and (8)  #####################################################################################
Var.Z_DM_mod5<-diag(data.matrix(X) %*% solve(t(data.matrix(X))%*%data.matrix(X)) %*% t(data.matrix(X))+(1/(2*(N-k_mod5)))* rep(Xe_comb$Z,times=L) %*% t(rep(Xe_comb$Z,times=L)))
Var.Z_DM_mod5

Var.Phi_DM_mod5<-(100*dnorm(Xe_comb$Z))^2*Var.Z_DM_mod5
Var.Phi_DM_mod5

############################# step 3: generate S normative samples #########################################################################################
S<-20000
Epsilon_s<-replicate(S,rnorm(n=N,0,sigma_epsilon),simplify =T)
dim(Epsilon_s)
y_s_mod5<-Beta[1]+Beta[2]*X$Age+Beta[3]*X$Sex+Beta[4]*X$Age2+Beta[5]*X$AgeSex+Beta[6]*X$Age2Sex+Epsilon_s 

#############################  step 4: estimate the model parameters ##########################################################################################

Y_hat_s_mod5<-Epsilon_hat_S_mod5<-Z_hat_s_mod5<-Phi_hat_s_mod5<-matrix(0,nrow=N,ncol = S)
sigma_epsilon_hat_mod5<-matrix(0,1,ncol = S)
Betas_hat_mod5<-matrix(0,nrow=k_mod5,ncol = S)

for(i in 1:S){

model5 <-lm(y_s_mod5[,i] ~ Age+Sex+Age2+AgeSex+Age2Sex, data=data.frame(X,y_s_mod5[,i])) 
Betas_hat_mod5[,i]<-model5$coefficients
Y_hat_s_mod5[,i]<-model5$fitted.values
Epsilon_hat_S_mod5[,i]<-model5$residuals
sigma_epsilon_hat_mod5[i]<-sigma(model5)

#############################   step 5: estimate  Z-scores and PR-scores   ###################################################################################

Z_hat_s_mod5[,i]<-(y_mod5-Y_hat_s_mod5[,i])/sigma_epsilon_hat_mod5[i]
Phi_hat_s_mod5[,i]<-pnorm(Z_hat_s_mod5[,i])*100

}

##############################    step 6: compute true values    ##############################################################################################
Z_bar_S_mod5<-apply(Z_hat_s_mod5,MARGIN=1,FUN=mean)
Var.Z_S_mod5<-apply(Z_hat_s_mod5,MARGIN=1,FUN=var)

Phi_bar_S_mod5<-apply(Phi_hat_s_mod5,MARGIN=1,FUN=mean)
Var.Phi_S_mod5<-apply(Phi_hat_s_mod5,MARGIN=1,FUN=var)

# take average across L replications (Note that they are all equals)

Z_bar_mod5<-Var.Z_mod5<-Phi_bar_mod5<-Var.Phi_mod5<-matrix(0,Q*length(Z),1)
Var.Z_DM_mod5_L<-matrix(0,Q*length(Z),1)
Var.Phi_DM_mod5_L<-matrix(0,Q*length(Z),1)


t<-0
for(l in 1:length(Z)){
for(j in 1:length(Sex)){
for(i in 1:length(Age)){
      
t<-t+1
Z_bar_mod5[t]<-mean(Z_bar_S_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_mod5[t]<-mean(Var.Z_S_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Z_DM_mod5_L[t]<-mean(Var.Z_DM_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
Phi_bar_mod5[t]<-mean(Phi_bar_S_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_mod5[t]<-mean(Var.Phi_S_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
Var.Phi_DM_mod5_L[t]<-mean(Var.Phi_DM_mod5[X$Age==Age[i] & X$Sex==Sex[j] & Xe_comb$Z==Z[l]])
      
 }}}

###############################     step 7: compute  Bias   ###################################################################################################

# Absolute bias for Z_hat
AB_Z_hat_mod5_L_1<-(Z_bar_mod5-Xe_comb$Z)
summary(AB_Z_hat_mod5_L_1)  
# Absolute bias for PR(Z_hat)
AB_Phi_hat_mod5_L_1<-(Phi_bar_mod5-Xe_comb$Phi)
summary(AB_Phi_hat_mod5_L_1)  
# R.B. for V(Z_hat)
RB_Var_Z_DM_mod5_L_1<-(Var.Z_DM_mod5_L-Var.Z_mod5)/Var.Z_mod5
RB_Var_Z_DM_mod5_L_1
cbind(Var.Z_DM_mod5_L,Var.Z_mod5)
# R.B. for V(PR(Z_hat))
RB_Var_Phi_DM_mod5_L_1<-(Var.Phi_DM_mod5_L-Var.Phi_mod5)/Var.Phi_mod5
RB_Var_Phi_DM_mod5_L_1
summary(RB_Var_Phi_DM_mod5_L_1)  

### Z-scores
## Figure S.A.4   Relative Bias of Equation (7)
# Sex=0
FigureS.A.4_Sex0<-round(cbind(AGE,RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.4_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.4_Sex0
# Sex=1
FigureS.A.4_Sex1<-round(cbind(AGE,RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Z_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.4_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.4_Sex1

## Figure S.A.6   Absolute Bias of Z_hat
# Sex=0
FigureS.A.6_Sex0<-round(cbind(AGE,AB_Z_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Z_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Z_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                AB_Z_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Z_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Z_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.6_Sex0)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.6_Sex0
# Sex=1
FigureS.A.6_Sex1<-round(cbind(AGE,AB_Z_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Z_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Z_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                AB_Z_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Z_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Z_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.A.6_Sex1)<-c("Age","Z=-2.5","Z=-2","Z=-1.5","Z=2.5","Z=2","Z=1.5")
FigureS.A.6_Sex1

### PR-scores
## Figure S.B.1.29   Relative Bias of Equation (8)
# Sex=0
FigureS.B.1.29_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.29_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.29_Sex0
# Sex=1
FigureS.B.1.29_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],RB_Var_Phi_DM_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.29_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.29_Sex1

## Figure S.B.1.34   Absolute Bias of Equation (6)
# Sex=0
FigureS.B.1.34_Sex0<-round(cbind(AGE,AB_Phi_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2.5],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-2],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2.5],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==2],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==0 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.34_Sex0)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.34_Sex0
# Sex=1
FigureS.B.1.34_Sex1<-round(cbind(AGE,AB_Phi_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2.5],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-2],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==-1.5],
                                 AB_Phi_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2.5],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==2],AB_Phi_hat_mod5_L_1[Xe_comb$Sex==1 & Xe_comb$Z==1.5]),3)
colnames(FigureS.B.1.34_Sex1)<-c("Age","PR(-2.5)","PR(-2)","PR(-1.5)","PR(2.5)","PR(2)","PR(1.5)")
FigureS.B.1.34_Sex1

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
### compute equations (7) and (8)
# model (1)
Var.Z_DM_mod1_additional<-diag(data.matrix(X_additional[,1:3]) %*% solve(t(data.matrix(X[,1:3]))%*%data.matrix(X[,1:3])) %*% t(data.matrix(X_additional[,1:3]))+(1/(2*(N-k_mod1)))* rep(Xe_comb_additional$Z,times=L) %*% t(rep(Xe_comb_additional$Z,times=L)))
Var.Z_DM_mod1_additional
Var.Phi_DM_mod1_additional<-(100*dnorm(Xe_comb_additional$Z))^2*Var.Z_DM_mod1_additional
Var.Phi_DM_mod1_additional
# model (2)
#Var.Z_DM_mod2_additional<-diag(data.matrix(X_additional[,1:4]) %*% solve(t(data.matrix(X[,1:4]))%*%data.matrix(X[,1:4])) %*% t(data.matrix(X_additional[,1:4]))+(1/(2*(N-k_mod2)))* rep(Xe_comb_additional$Z,times=L) %*% t(rep(Xe_comb_additional$Z,times=L)))
#Var.Z_DM_mod2_additional
#Var.Phi_DM_mod2_additional<-(100*dnorm(Xe_comb_additional$Z))^2*Var.Z_DM_mod2_additional
#Var.Phi_DM_mod2_additional
# model (3)
#Var.Z_DM_mod3_additional<-diag(data.matrix(X_additional[,c(1:3,5)]) %*% solve(t(data.matrix(X[,c(1:3,5)]))%*%data.matrix(X[,c(1:3,5)])) %*% t(data.matrix(X_additional[,c(1:3,5)]))+(1/(2*(N-k_mod3)))* rep(Xe_comb_additional$Z,times=L) %*% t(rep(Xe_comb_additional$Z,times=L)))
#Var.Z_DM_mod3_additional
#Var.Phi_DM_mod3_additional<-(100*dnorm(Xe_comb_additional$Z))^2*Var.Z_DM_mod3_additional
#Var.Phi_DM_mod3_additional
# model (4)
#Var.Z_DM_mod4_additional<-diag(data.matrix(X_additional[,-6]) %*% solve(t(data.matrix(X[,-6]))%*%data.matrix(X[,-6])) %*% t(data.matrix(X_additional[,-6]))+(1/(2*(N-k_mod4)))* rep(Xe_comb_additional$Z,times=L) %*% t(rep(Xe_comb_additional$Z,times=L)))
#Var.Z_DM_mod4_additional
#Var.Phi_DM_mod4_additional<-(100*dnorm(Xe_comb_additional$Z))^2*Var.Z_DM_mod4_additional
#Var.Phi_DM_mod4_additional
# model (5)
#Var.Z_DM_mod5_additional<-diag(data.matrix(X_additional) %*% solve(t(data.matrix(X))%*%data.matrix(X)) %*% t(data.matrix(X_additional))+(1/(2*(N-k_mod5)))* rep(Xe_comb_additional$Z,times=L) %*% t(rep(Xe_comb_additional$Z,times=L)))
#Var.Z_DM_mod5_additional
#Var.Phi_DM_mod5_additional<-(100*dnorm(Xe_comb_additional$Z))^2*Var.Z_DM_mod5_additional
#Var.Phi_DM_mod5_additional
### calculate Z_hat and PR_hat using beta_hat sigma_e_hat from the S normative samples
Z_hat_s_mod1_additional<-Phi_hat_s_mod1_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)
#Z_hat_s_mod2_additional<-Phi_hat_s_mod2_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)
#Z_hat_s_mod3_additional<-Phi_hat_s_mod3_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)
#Z_hat_s_mod4_additional<-Phi_hat_s_mod4_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)
#Z_hat_s_mod5_additional<-Phi_hat_s_mod5_additional<-matrix(0,nrow=length(Z_additional)*Q*L,ncol = S)

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
}

###### model 1 (for other models replace "mod1" with "mod#ofmodel#)
### compute true values   
Z_bar_S_mod1_additional<-apply(Z_hat_s_mod1_additional,MARGIN=1,FUN=mean) 
Var.Z_S_mod1_additional<-apply(Z_hat_s_mod1_additional,MARGIN=1,FUN=var)  
Phi_bar_S_mod1_additional<-apply(Phi_hat_s_mod1_additional,MARGIN=1,FUN=mean)
Var.Phi_S_mod1_additional<-apply(Phi_hat_s_mod1_additional,MARGIN=1,FUN=var)

Z_bar_mod1_additional<-Var.Z_mod1_additional<-Phi_bar_mod1_additional<-Var.Phi_mod1_additional<-matrix(0,Q*length(Z_additional),1)
Var.Z_DM_mod1_L_additional<-matrix(0,Q*length(Z_additional),1)
Var.Phi_DM_mod1_L_additional<-matrix(0,Q*length(Z_additional),1)

t<-0
for(l in 1:length(Z_additional)){
  for(j in 1:length(Sex)){
    for(i in 1:length(Age)){
      
      t<-t+1
      
      Z_bar_mod1_additional[t]<-mean(Z_bar_S_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      Var.Z_mod1_additional[t]<-mean(Var.Z_S_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      Var.Z_DM_mod1_L_additional[t]<-mean(Var.Z_DM_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      
      Phi_bar_mod1_additional[t]<-mean(Phi_bar_S_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      Var.Phi_mod1_additional[t]<-mean(Var.Phi_S_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      Var.Phi_DM_mod1_L_additional[t]<-mean(Var.Phi_DM_mod1_additional[X_additional$Age==Age[i] & X_additional$Sex==Sex[j] & Xe_comb_additional$Z==Z_additional[l]])
      
    }}}

### compute Relative Biases    
### Absolute Bias for Z_hat
AB_Z_hat_mod1_L_1_additional<-(Z_bar_mod1_additional-Xe_comb_additional$Z)
summary(AB_Z_hat_mod1_L_1_additional)  
### Absolute Bias for PR(Z_hat)
AB_Phi_hat_mod1_L_1_additional<-(Phi_bar_mod1_additional-Xe_comb_additional$Phi)
summary(AB_Phi_hat_mod1_L_1_additional)  
### R.B. for V(Z_hat)
RB_Var_Z_DM_mod1_L_1_additional<-(Var.Z_DM_mod1_L_additional-Var.Z_mod1_additional)/Var.Z_mod1_additional
RB_Var_Z_DM_mod1_L_1_additional
summary(RB_Var_Z_DM_mod1_L_1_additional)  
### R.B. for V(PR(Z_hat))
RB_Var_Phi_DM_mod1_L_1_additional<-(Var.Phi_DM_mod1_L_additional-Var.Phi_mod1_additional)/Var.Phi_mod1_additional
RB_Var_Phi_DM_mod1_L_1_additional
summary(RB_Var_Phi_DM_mod1_L_1_additional)  

### PR-scores
################### model (1)
## Figure S.A.7  Relative Bias of Equation (8)
# Sex=0
FigureS.A.7_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                                  RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                                  RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                                  RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.7_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.7_Sex0
# Sex=1
FigureS.A.7_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                                  RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                                  RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                                  RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],RB_Var_Phi_DM_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.7_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.7_Sex1

## Figure S.A.9 Absolute Bias of Equation (6)
# Sex=0
FigureS.A.9_Sex0<-round(cbind(AGE,AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                                  AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                                  AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                                  AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.9_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.9_Sex0
# Sex=1
FigureS.A.9_Sex1<-round(cbind(AGE,AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                                  AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                                  AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                                  AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],AB_Phi_hat_mod1_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.9_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.9_Sex1

################### model (5)
## Figure S.A.8  Relative Bias of Equation (8)
# Sex=0
FigureS.A.8_Sex0<-round(cbind(AGE,RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                              RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                              RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                              RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.8_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.8_Sex0
# Sex=1
FigureS.A.8_Sex1<-round(cbind(AGE,RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                              RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                              RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                              RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],RB_Var_Phi_DM_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.8_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.8_Sex1

## Figure S.A.10 Absolute Bias of Equation (6)
# Sex=0
FigureS.A.10_Sex0<-round(cbind(AGE,AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-2.3263],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.9600],
                              AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.6449],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==-1.2816],
                              AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.2816],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.6449],
                              AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==1.9600],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==0 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.10_Sex0)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.10_Sex0
# Sex=1
FigureS.A.10_Sex1<-round(cbind(AGE,AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-2.3263],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.9600],
                              AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.6449],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==-1.2816],
                              AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.2816],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.6449],
                              AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==1.9600],AB_Phi_hat_mod5_L_1_additional[Xe_comb_additional$Sex==1 & Xe_comb_additional$Z==2.3263]),3)
colnames(FigureS.A.10_Sex1)<-c("Age","PR=1","PR=2.5","PR=5","PR=10","PR=90","PR=95", "PR=97.5", "PR=99")
FigureS.A.10_Sex1
