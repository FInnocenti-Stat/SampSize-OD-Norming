#######################################################################################################################################
#### Sample Size Calculation and Optimal Design for Regression-Based Norming of Tests and Questionnaires
#######################################################################################################################################
# Francesco Innocenti, Frans E. S. Tan, Math J. J. M. Candel, and Gerard J. P. van Breukelen
#######################################################################################################################################
##################################### Robustness of the optimal design: Maximin Design ################################################
########################################################################################################################################
rm(list=ls(all=T))
N<-16*6*5*13

################################# Generate the design matrix for each design under each model ##########################################

########### Optimal designs in Table 2
######## true model = model 1 = E(Y)= b0+b1X1+b2X2
p<-3 # number of parameters in the model
# d1= (1,0), (1,1), (-1,0), (-1,1), with equal weight=1/4
X_d1_M1<-matrix(cbind(rep(1,times=N),rep(c(1,-1),each=2),rep(c(0,1),times=2)),N,3)
X_d1_M1
# d2= (1,0), (0,0), (-1,0), (1,1), (0,1), (-1,1) with equal weight=1/6
X_d2_M1<-matrix(cbind(rep(1,times=N),rep(c(1,0,-1),each=2),rep(c(0,1),times=2)),N,3)
X_d2_M1
# d3= (-1,0), (1,0), (-1,1), (1,1) with equal weight=3/16 and (0,1), (0,0) with equal weight=2/16
X_d3_M1<-matrix(cbind(rep(1,times=N),c(rep(c(1,-1), times=36),rep(0,times=24)),c(rep(c(0,1),each=36),rep(c(0,1),each=12))),N,3)
X_d3_M1

l=0;g=0;d=0;e=0
for(i in 1:N){
 
  ifelse(all(X_d1_M1[i,]==rbind(c(1,1,0))),l<-l+1,l<-l)
  ifelse(all(X_d2_M1[i,]==rbind(c(1,1,0))),g<-g+1,g<-g)
  ifelse(all(X_d3_M1[i,]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d3_M1[i,]==rbind(c(1,0,0))),e<-e+1,e<-e)
  
}
l==N/4 # N times weigth
g==N/6
d==N*(3/16)
e==N*(2/16)



######## true model = model 2 = E(Y)= b0+b1X1+b2X2+b3x1x1
p<-4 # number of parameters in the model
# d1= (1,0), (1,1), (-1,0), (-1,1), with equal weight=1/4
X_d1_M2<-matrix(cbind(rep(1,times=N),rep(c(1,-1),each=2),rep(c(0,1),times=2),rep(c(1,-1),each=2)^2),N,p)
X_d1_M2
# d2= (1,0), (0,0), (-1,0), (1,1), (0,1), (-1,1) with equal weight=1/6
X_d2_M2<-matrix(cbind(rep(1,times=N),rep(c(1,0,-1),each=2),rep(c(0,1),times=2),X_d2_M1[,2]^2),N,p)
X_d2_M2  
# d3= (-1,0), (1,0), (-1,1), (1,1) with equal weight=3/16 and (0,1), (0,0) with equal weight=2/16
X_d3_M2<-matrix(cbind(rep(1,times=N),c(rep(c(1,-1), times=36),rep(0,times=24)),c(rep(c(0,1),each=36),rep(c(0,1),each=12)),
                                     c(rep(c(1,-1), times=36),rep(0,times=24))^2),N,p)
X_d3_M2

l=0;g=0;d=0;e=0
for(i in 1:N){
  
  ifelse(all(X_d2_M2[i,]==rbind(c(1,1,0,1))),l<-l+1,l<-l)
  ifelse(all(X_d2_M2[i,]==rbind(c(1,1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d3_M2[i,]==rbind(c(1,1,0,1))),d<-d+1,d<-d)
  ifelse(all(X_d3_M2[i,]==rbind(c(1,0,0,0))),e<-e+1,e<-e)
  
}

l==N/6  # N times weigth
g==N/6
d==N*(3/16)
e==N*(2/16)


######## true model = model 3 = E(Y)= b0+b1X1+b2X2+b4x1x2
p<-4 # number of parameters in the model
# d1= (1,0), (1,1), (-1,0), (-1,1), with equal weight=1/4
X_d1_M3<-matrix(cbind(rep(1,times=N),rep(c(1,-1),each=2),rep(c(0,1),times=2),rep(c(1,-1),each=2)*rep(c(0,1),times=2)),N,p)
X_d1_M3
# d2= (1,0), (0,0), (-1,0), (1,1), (0,1), (-1,1) with equal weight=1/6
X_d2_M3<-matrix(cbind(rep(1,times=N),rep(c(1,0,-1),each=2),rep(c(0,1),times=2),X_d2_M1[,2]*X_d2_M1[,3]),N,p)
X_d2_M3

# d3= (-1,0), (1,0), (-1,1), (1,1) with equal weight=3/16 and (0,1), (0,0) with equal weight=2/16
X_d3_M3<-matrix(cbind(rep(1,times=N),c(rep(c(1,-1), times=36),rep(0,times=24)),c(rep(c(0,1),each=36),rep(c(0,1),each=12)),
                      c(rep(c(1,-1), times=36),rep(0,times=24))*c(rep(c(0,1),each=36),rep(c(0,1),each=12))),N,p)
X_d3_M3
         
l=0;g=0;d=0;e=0
for(i in 1:N){
  
  ifelse(all(X_d1_M3[i,1:3]==rbind(c(1,1,0))),l<-l+1,l<-l)
  ifelse(all(X_d2_M3[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d3_M3[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d3_M3[i,1:3]==rbind(c(1,0,0))),e<-e+1,e<-e)
  
}

l==N/4  # N times weigth
g==N/6
d==N*(3/16)
e==N*(2/16)


######## true model = model 4 = E(Y)= b0+b1X1+b2X2+b3x1x1+b4x1x2
p<-5 # number of parameters in the model
# d1= (1,0), (1,1), (-1,0), (-1,1), with equal weight=1/4
X_d1_M4<-matrix(cbind(rep(1,times=N),rep(c(1,-1),each=2),rep(c(0,1),times=2),rep(c(1,-1),each=2)^2,X_d1_M2[,2]*X_d1_M2[,3]),N,p)
X_d1_M4
# d2= (1,0), (0,0), (-1,0), (1,1), (0,1), (-1,1) with equal weight=1/6
X_d2_M4<-matrix(cbind(rep(1,times=N),rep(c(1,0,-1),each=2),rep(c(0,1),times=2),X_d2_M1[,2]^2,X_d2_M2[,2]*X_d2_M2[,3]),N,p)
X_d2_M4

# d3= (-1,0), (1,0), (-1,1), (1,1) with equal weight=3/16 and (0,1), (0,0) with equal weight=2/16
X_d3_M4<-matrix(cbind(rep(1,times=N),c(rep(c(1,-1), times=36),rep(0,times=24)),c(rep(c(0,1),each=36),rep(c(0,1),each=12)),
                      c(rep(c(1,-1), times=36),rep(0,times=24))^2,X_d3_M2[,2]*X_d3_M2[,3]),N,p)
X_d3_M4

l=0;g=0;d=0;e=0
for(i in 1:N){
  
  ifelse(all(X_d2_M4[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d3_M4[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d3_M4[i,1:3]==rbind(c(1,0,0))),e<-e+1,e<-e)
  
}

# N times weigth
g==N/6
d==N*(3/16)
e==N*(2/16)



######## true model = model 5 = E(Y)= b0+b1X1+b2X2+b3x1x1+b4x1x2+b5x1x1x2
p<-6 # number of parameters in the model
# d1= (1,0), (1,1), (-1,0), (-1,1), with equal weight=1/4
X_d1_M5<-matrix(cbind(rep(1,times=N),rep(c(1,-1),each=2),rep(c(0,1),times=2),rep(c(1,-1),each=2)^2,X_d1_M2[,2]*X_d1_M2[,3],X_d1_M2[,4]*X_d1_M2[,3]),N,p)
X_d1_M5
# d2= (1,0), (0,0), (-1,0), (1,1), (0,1), (-1,1) with equal weight=1/6
X_d2_M5<-matrix(cbind(rep(1,times=N),rep(c(1,0,-1),each=2),rep(c(0,1),times=2),X_d2_M1[,2]^2,X_d2_M2[,2]*X_d2_M2[,3],X_d2_M2[,4]*X_d2_M2[,3]),N,p)
X_d2_M5

# d3= (-1,0), (1,0), (-1,1), (1,1) with equal weight=3/16 and (0,1), (0,0) with equal weight=2/16
X_d3_M5<-matrix(cbind(rep(1,times=N),c(rep(c(1,-1), times=36),rep(0,times=24)),c(rep(c(0,1),each=36),rep(c(0,1),each=12)),
                      c(rep(c(1,-1), times=36),rep(0,times=24))^2,X_d3_M2[,2]*X_d3_M2[,3],X_d3_M2[,4]*X_d3_M2[,3]),N,p)
X_d3_M5

l=0;g=0;d=0;e=0
for(i in 1:N){
  
  ifelse(all(X_d2_M5[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d3_M5[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d3_M5[i,1:3]==rbind(c(1,0,0))),e<-e+1,e<-e)
  
}

#l==N/4  # N times weigth
g==N/6
d==N*(3/16)
e==N*(2/16)



########### Equidistant age levels designs (with more than 3 age levels) 
AGE<-seq(from=20,to=80,by=1)
d_tilde<-min(AGE)+(max(AGE)-min(AGE))/2
Age_scaled<-(AGE-d_tilde)/(max(AGE)-d_tilde)
cbind(AGE,Age_scaled)

# 4 age levels (equal weight 1/8 per sex's level)
seq(20,80,by=20)
FourEquiDistpts<-Age_scaled[c(which(AGE==seq(20,80,by=20)[1]),which(AGE==seq(20,80,by=20)[2]),which(AGE==seq(20,80,by=20)[3]),which(AGE==seq(20,80,by=20)[4]))]
FourEquiDistpts
X_4equip<-matrix(cbind(rep(1,times=N),rep(FourEquiDistpts,times=(N/4)),rep(c(0,1),each=4)),N,3)
X_4equip

X_d4_M1<-X_4equip # model 1 =true model
X_d4_M1
X_d4_M2<-cbind(X_4equip,X_4equip[,2]^2) # model 2 =true model
X_d4_M2
X_d4_M3<-cbind(X_4equip,X_4equip[,2]*X_4equip[,3]) # model 3 =true model
X_d4_M3
X_d4_M4<-cbind(X_4equip,X_4equip[,2]^2,X_4equip[,2]*X_4equip[,3]) # model 4 =true model
X_d4_M4
X_d4_M5<-cbind(X_4equip,X_4equip[,2]^2,X_4equip[,2]*X_4equip[,3],(X_4equip[,2]^2)*X_4equip[,3]) # model 5 =true model
X_d4_M5

l=0;g=0;d=0;e=0
for(i in 1:N){
  ifelse(all(X_d4_M5[i,1:3]==rbind(c(1,1,0))),l<-l+1,l<-l)
  ifelse(all(X_d4_M1[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d4_M2[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d4_M3[i,1:3]==rbind(c(1,1,0))),e<-e+1,e<-e)
  
}
l==N/8  # N times weigth
g==N/8
d==N/8
e==N/8

# 5 age levels (equal weight 1/10 per sex's level)
seq(20,80,by=15)
FiveEquiDistpts<-Age_scaled[c(which(AGE==seq(20,80,by=15)[1]),which(AGE==seq(20,80,by=15)[2]),which(AGE==seq(20,80,by=15)[3]),which(AGE==seq(20,80,by=15)[4]),which(AGE==seq(20,80,by=15)[5]))]
FiveEquiDistpts
X_5equip<-matrix(cbind(rep(1,times=N),rep(FiveEquiDistpts,times=(N/5)),rep(c(0,1),each=5)),N,3)
X_5equip

X_d5_M1<-X_5equip # model 1 =true model
X_d5_M1
X_d5_M2<-cbind(X_5equip,X_5equip[,2]^2) # model 2 =true model
X_d5_M2
X_d5_M3<-cbind(X_5equip,X_5equip[,2]*X_5equip[,3]) # model 3 =true model
X_d5_M3
X_d5_M4<-cbind(X_5equip,X_5equip[,2]^2,X_5equip[,2]*X_5equip[,3]) # model 4 =true model
X_d5_M4
X_d5_M5<-cbind(X_5equip,X_5equip[,2]^2,X_5equip[,2]*X_5equip[,3],(X_5equip[,2]^2)*X_5equip[,3]) # model 5 =true model
X_d5_M5

l=0;g=0;d=0;e=0
for(i in 1:N){
  ifelse(all(X_d5_M5[i,1:3]==rbind(c(1,1,0))),l<-l+1,l<-l)
  ifelse(all(X_d5_M1[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d5_M2[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d5_M3[i,1:3]==rbind(c(1,1,0))),e<-e+1,e<-e)
  
}
l==N/10 # N times weigth
g==N/10
d==N/10
e==N/10

# 6 age levels (equal weight 1/12 per sex's level)
seq(20,80,by=12)
SixEquiDistpts<-Age_scaled[c(which(AGE==seq(20,80,by=12)[1]),which(AGE==seq(20,80,by=12)[2]),which(AGE==seq(20,80,by=12)[3]),which(AGE==seq(20,80,by=12)[4]),which(AGE==seq(20,80,by=12)[5]),which(AGE==seq(20,80,by=12)[6]))]
SixEquiDistpts
X_6equip<-matrix(cbind(rep(1,times=N),rep(SixEquiDistpts,times=(N/6)),rep(c(0,1),each=6)),N,3)
X_6equip

X_d6_M1<-X_6equip # model 1 =true model
X_d6_M1
X_d6_M2<-cbind(X_6equip,X_6equip[,2]^2) # model 2 =true model
X_d6_M2
X_d6_M3<-cbind(X_6equip,X_6equip[,2]*X_6equip[,3]) # model 3 =true model
X_d6_M3
X_d6_M4<-cbind(X_6equip,X_6equip[,2]^2,X_6equip[,2]*X_6equip[,3]) # model 4 =true model
X_d6_M4
X_d6_M5<-cbind(X_6equip,X_6equip[,2]^2,X_6equip[,2]*X_6equip[,3],(X_6equip[,2]^2)*X_6equip[,3]) # model 5 =true model
X_d6_M5

l=0;g=0;d=0;e=0
for(i in 1:N){
  ifelse(all(X_d6_M5[i,1:3]==rbind(c(1,1,0))),l<-l+1,l<-l)
  ifelse(all(X_d6_M1[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d6_M2[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d6_M3[i,1:3]==rbind(c(1,1,0))),e<-e+1,e<-e)
  
}
l==N/12 # 1/12 is the weight
g==N/12
d==N/12
e==N/12

# 13 age levels (equal weight 1/26 per sex's level)
seq(20,80,by=5)
ThirteenEquiDistpts<-Age_scaled[c(which(AGE==seq(20,80,by=5)[1]),which(AGE==seq(20,80,by=5)[2]),which(AGE==seq(20,80,by=5)[3]),which(AGE==seq(20,80,by=5)[4]),which(AGE==seq(20,80,by=5)[5]),which(AGE==seq(20,80,by=5)[6]),
                                  which(AGE==seq(20,80,by=5)[7]),which(AGE==seq(20,80,by=5)[8]),which(AGE==seq(20,80,by=5)[9]),which(AGE==seq(20,80,by=5)[10]),which(AGE==seq(20,80,by=5)[11]),which(AGE==seq(20,80,by=5)[12]),
                                  which(AGE==seq(20,80,by=5)[13]))]
ThirteenEquiDistpts
X_13equip<-matrix(cbind(rep(1,times=N),rep(ThirteenEquiDistpts,times=(N/13)),rep(c(0,1),each=13)),N,3)
X_13equip

X_d13_M1<-X_13equip # model 1 =true model
X_d13_M1
X_d13_M2<-cbind(X_13equip,X_13equip[,2]^2) # model 2 =true model
X_d13_M2
X_d13_M3<-cbind(X_13equip,X_13equip[,2]*X_13equip[,3]) # model 3 =true model
X_d13_M3
X_d13_M4<-cbind(X_13equip,X_13equip[,2]^2,X_13equip[,2]*X_13equip[,3]) # model 4 =true model
X_d13_M4
X_d13_M5<-cbind(X_13equip,X_13equip[,2]^2,X_13equip[,2]*X_13equip[,3],(X_13equip[,2]^2)*X_13equip[,3]) # model 5 =true model
X_d13_M5

l=0;g=0;d=0;e=0
for(i in 1:N){
  ifelse(all(X_d13_M5[i,1:3]==rbind(c(1,1,0))),l<-l+1,l<-l)
  ifelse(all(X_d13_M1[i,1:3]==rbind(c(1,1,1))),g<-g+1,g<-g)
  ifelse(all(X_d13_M2[i,1:3]==rbind(c(1,1,0))),d<-d+1,d<-d)
  ifelse(all(X_d13_M3[i,1:3]==rbind(c(1,1,0))),e<-e+1,e<-e)
  
}
l==N/26 # 1/26 is the weight
g==N/26
d==N/26
e==N/26


########################################## Compute (X'X)^-1 for each design under each model ###########################################################

# under model (1)
V_d1_M1<-solve(t(X_d1_M1)%*%X_d1_M1)   # design with 2 age levels
V_d2_M1<-solve(t(X_d2_M1)%*%X_d2_M1)   # balanced design with 3 age levels
V_d3_M1<-solve(t(X_d3_M1)%*%X_d3_M1)   # unbalanced design with 3 age levels
V_d4_M1<-solve(t(X_d4_M1)%*%X_d4_M1)   # design with 4 age levels
V_d5_M1<-solve(t(X_d5_M1)%*%X_d5_M1)   # design with 5 age levels
V_d6_M1<-solve(t(X_d6_M1)%*%X_d6_M1)   # design with 6 age levels
V_d7_M1<-solve(t(X_d13_M1)%*%X_d13_M1) # design with 13 age levels

# under model (2)
#V_d1_M2<-solve(t(X_d1_M2)%*%X_d1_M2)  # the effect of Age^2 is not identifiable for a design with only 2 age levels
V_d2_M2<-solve(t(X_d2_M2)%*%X_d2_M2)
V_d3_M2<-solve(t(X_d3_M2)%*%X_d3_M2)
V_d4_M2<-solve(t(X_d4_M2)%*%X_d4_M2)
V_d5_M2<-solve(t(X_d5_M2)%*%X_d5_M2)
V_d6_M2<-solve(t(X_d6_M2)%*%X_d6_M2)
V_d7_M2<-solve(t(X_d13_M2)%*%X_d13_M2)

# under model (3)
V_d1_M3<-solve(t(X_d1_M3)%*%X_d1_M3)
V_d2_M3<-solve(t(X_d2_M3)%*%X_d2_M3)
V_d3_M3<-solve(t(X_d3_M3)%*%X_d3_M3)
V_d4_M3<-solve(t(X_d4_M3)%*%X_d4_M3)
V_d5_M3<-solve(t(X_d5_M3)%*%X_d5_M3)
V_d6_M3<-solve(t(X_d6_M3)%*%X_d6_M3)
V_d7_M3<-solve(t(X_d13_M3)%*%X_d13_M3)

# under model (4)
#V_d1_M4<-solve(t(X_d1_M4)%*%X_d1_M4) # the effect of Age^2 is not identifiable for a design with only 2 age levels
V_d2_M4<-solve(t(X_d2_M4)%*%X_d2_M4)
V_d3_M4<-solve(t(X_d3_M4)%*%X_d3_M4)
V_d4_M4<-solve(t(X_d4_M4)%*%X_d4_M4)
V_d5_M4<-solve(t(X_d5_M4)%*%X_d5_M4)
V_d6_M4<-solve(t(X_d6_M4)%*%X_d6_M4)
V_d7_M4<-solve(t(X_d13_M4)%*%X_d13_M4)

# under model (5)
#V_d1_M5<-solve(t(X_d1_M5)%*%X_d1_M5) # the effect of Age^2 is not identifiable for a design with only 2 age levels
V_d2_M5<-solve(t(X_d2_M5)%*%X_d2_M5)
V_d3_M5<-solve(t(X_d3_M5)%*%X_d3_M5)
V_d4_M5<-solve(t(X_d4_M5)%*%X_d4_M5)
V_d5_M5<-solve(t(X_d5_M5)%*%X_d5_M5)
V_d6_M5<-solve(t(X_d6_M5)%*%X_d6_M5)
V_d7_M5<-solve(t(X_d13_M5)%*%X_d13_M5)


############################################ Maximin Design based on the RE criterion #######################################################################

# generate x0
AGE<-seq(from=20,to=80,by=1)
d_tilde<-min(AGE)+(max(AGE)-min(AGE))/2
Age_scaled<-(AGE-d_tilde)/(max(AGE)-d_tilde)
cbind(AGE,Age_scaled)

Sex<-c(0,1)

X_comb<-expand.grid(Age_scaled,Sex)
colnames(X_comb)<-c("Age","Sex")
X_comb

x0<-data.frame(1,X_comb$Age,X_comb$Sex,X_comb$Age^2,X_comb$Age*X_comb$Sex,(X_comb$Age^2)*X_comb$Sex)
colnames(x0)<-c("Int", "Age","Sex", "Age2", "AgeSex", "Age2Sex")
x0

### Find the minimum RE (i.e. equation (11) for Z=0) across x0 values --> results in Table 3

# note: d1 = 2 age levels design, d2 = 3 age levels balanced design, d3 = 3 age levels unbalanced design
# d4 = 4 age levels design, d5 = 5 age levels design, d6 = 6 age levels design, d7 = 13 age levels design

# under model (1)
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3])))))
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d2_M1%*%t(data.matrix(x0[,1:3])))))
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d3_M1%*%t(data.matrix(x0[,1:3])))))
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d4_M1%*%t(data.matrix(x0[,1:3])))))
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d5_M1%*%t(data.matrix(x0[,1:3])))))
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d6_M1%*%t(data.matrix(x0[,1:3])))))
min((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))/(N*diag(data.matrix(x0[,1:3])%*%V_d7_M1%*%t(data.matrix(x0[,1:3])))))


# under model (2): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
#min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d1_M2%*%t(data.matrix(x0[,1:4])))))
min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4])))))
min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d3_M2%*%t(data.matrix(x0[,1:4])))))
min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d4_M2%*%t(data.matrix(x0[,1:4])))))
min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d5_M2%*%t(data.matrix(x0[,1:4])))))
min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d6_M2%*%t(data.matrix(x0[,1:4])))))
min((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))/(N*diag(data.matrix(x0[,1:4])%*%V_d7_M2%*%t(data.matrix(x0[,1:4])))))

# under model (3)
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)])))))
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d2_M3%*%t(data.matrix(x0[,c(1:3,5)])))))
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d3_M3%*%t(data.matrix(x0[,c(1:3,5)])))))
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d4_M3%*%t(data.matrix(x0[,c(1:3,5)])))))
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d5_M3%*%t(data.matrix(x0[,c(1:3,5)])))))
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d6_M3%*%t(data.matrix(x0[,c(1:3,5)])))))
min((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d7_M3%*%t(data.matrix(x0[,c(1:3,5)])))))

# under model (4): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
#min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d1_M4%*%t(data.matrix(x0[,-6])))))
min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d2_M4%*%t(data.matrix(x0[,-6])))))
min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6])))))
min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d4_M4%*%t(data.matrix(x0[,-6])))))
min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d5_M4%*%t(data.matrix(x0[,-6])))))
min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d6_M4%*%t(data.matrix(x0[,-6])))))
min((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))/(N*diag(data.matrix(x0[,-6])%*%V_d7_M4%*%t(data.matrix(x0[,-6])))))

# under model (5): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
#min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d1_M5%*%t(data.matrix(x0)))))
min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0)))))
min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d3_M5%*%t(data.matrix(x0)))))
min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d4_M5%*%t(data.matrix(x0)))))
min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d5_M5%*%t(data.matrix(x0)))))
min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d6_M5%*%t(data.matrix(x0)))))
min((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))/(N*diag(data.matrix(x0)%*%V_d7_M5%*%t(data.matrix(x0)))))



### Find the minimum RE (i.e. equation (11)) for Z=+/-2 across x0 values --> results in Table 3

# note: d1 = 2 age levels design, d2 = 3 age levels balanced design, d3 = 3 age levels unbalanced design
# d4 = 4 age levels design, d5 = 5 age levels design, d6 = 6 age levels design, d7 = 13 age levels design

Z=2 # 0.5, 1, 1.5, 2.5, 3
  
# under model (1)
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d2_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d3_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d4_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d5_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d6_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:3])%*%V_d7_M1%*%t(data.matrix(x0[,1:3])))+((Z^2)/2)))


# under model (2): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
min(((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:4])%*%V_d3_M2%*%t(data.matrix(x0[,1:4])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:4])%*%V_d4_M2%*%t(data.matrix(x0[,1:4])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:4])%*%V_d5_M2%*%t(data.matrix(x0[,1:4])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:4])%*%V_d6_M2%*%t(data.matrix(x0[,1:4])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,1:4])%*%V_d7_M2%*%t(data.matrix(x0[,1:4])))+((Z^2)/2)))

# under model (3)
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d2_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d3_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d4_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d5_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d6_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d7_M3%*%t(data.matrix(x0[,c(1:3,5)])))+((Z^2)/2)))

# under model (4): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
min(((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,-6])%*%V_d2_M4%*%t(data.matrix(x0[,-6])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,-6])%*%V_d4_M4%*%t(data.matrix(x0[,-6])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,-6])%*%V_d5_M4%*%t(data.matrix(x0[,-6])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,-6])%*%V_d6_M4%*%t(data.matrix(x0[,-6])))+((Z^2)/2)))
min(((N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))+((Z^2)/2))/(N*diag(data.matrix(x0[,-6])%*%V_d7_M4%*%t(data.matrix(x0[,-6])))+((Z^2)/2)))

# under model (5): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
min(((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))+((Z^2)/2))/(N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0)))+((Z^2)/2)))
min(((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))+((Z^2)/2))/(N*diag(data.matrix(x0)%*%V_d3_M5%*%t(data.matrix(x0)))+((Z^2)/2)))
min(((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))+((Z^2)/2))/(N*diag(data.matrix(x0)%*%V_d4_M5%*%t(data.matrix(x0)))+((Z^2)/2)))
min(((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))+((Z^2)/2))/(N*diag(data.matrix(x0)%*%V_d5_M5%*%t(data.matrix(x0)))+((Z^2)/2)))
min(((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))+((Z^2)/2))/(N*diag(data.matrix(x0)%*%V_d6_M5%*%t(data.matrix(x0)))+((Z^2)/2)))
min(((N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))+((Z^2)/2))/(N*diag(data.matrix(x0)%*%V_d7_M5%*%t(data.matrix(x0)))+((Z^2)/2)))



################################################### Maximin design based on the efficiency criterion #########################################################
# To obtain the maximin design based on the efficiency criterion, 
# find the maximum standardized prediction variance over x0 for each design under each model --> results of Table S.A.3

# note: d1 = 2 age levels design, d2 = 3 age levels balanced design, d3 = 3 age levels unbalanced design
# d4 = 4 age levels design, d5 = 5 age levels design, d6 = 6 age levels design, d7 = 13 age levels design


# under model (1)
max(N*diag(data.matrix(x0[,1:3])%*%V_d1_M1%*%t(data.matrix(x0[,1:3]))))
max(N*diag(data.matrix(x0[,1:3])%*%V_d2_M1%*%t(data.matrix(x0[,1:3]))))
max(N*diag(data.matrix(x0[,1:3])%*%V_d3_M1%*%t(data.matrix(x0[,1:3]))))
max(N*diag(data.matrix(x0[,1:3])%*%V_d4_M1%*%t(data.matrix(x0[,1:3]))))
max(N*diag(data.matrix(x0[,1:3])%*%V_d5_M1%*%t(data.matrix(x0[,1:3]))))
max(N*diag(data.matrix(x0[,1:3])%*%V_d6_M1%*%t(data.matrix(x0[,1:3]))))
max(N*diag(data.matrix(x0[,1:3])%*%V_d7_M1%*%t(data.matrix(x0[,1:3]))))

# under model (2): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
#max(N*diag(data.matrix(x0[,1:4])%*%V_d1_M2%*%t(data.matrix(x0[,1:4]))))
max(N*diag(data.matrix(x0[,1:4])%*%V_d2_M2%*%t(data.matrix(x0[,1:4]))))
max(N*diag(data.matrix(x0[,1:4])%*%V_d3_M2%*%t(data.matrix(x0[,1:4]))))
max(N*diag(data.matrix(x0[,1:4])%*%V_d4_M2%*%t(data.matrix(x0[,1:4]))))
max(N*diag(data.matrix(x0[,1:4])%*%V_d5_M2%*%t(data.matrix(x0[,1:4]))))
max(N*diag(data.matrix(x0[,1:4])%*%V_d6_M2%*%t(data.matrix(x0[,1:4]))))
max(N*diag(data.matrix(x0[,1:4])%*%V_d7_M2%*%t(data.matrix(x0[,1:4]))))

# under model (3)
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d1_M3%*%t(data.matrix(x0[,c(1:3,5)]))))
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d2_M3%*%t(data.matrix(x0[,c(1:3,5)]))))
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d3_M3%*%t(data.matrix(x0[,c(1:3,5)]))))
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d4_M3%*%t(data.matrix(x0[,c(1:3,5)]))))
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d5_M3%*%t(data.matrix(x0[,c(1:3,5)]))))
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d6_M3%*%t(data.matrix(x0[,c(1:3,5)]))))
max(N*diag(data.matrix(x0[,c(1:3,5)])%*%V_d7_M3%*%t(data.matrix(x0[,c(1:3,5)]))))

# under model (4): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
#max(N*diag(data.matrix(x0[,-6])%*%V_d1_M4%*%t(data.matrix(x0[,-6]))))
max(N*diag(data.matrix(x0[,-6])%*%V_d2_M4%*%t(data.matrix(x0[,-6]))))
max(N*diag(data.matrix(x0[,-6])%*%V_d3_M4%*%t(data.matrix(x0[,-6]))))
max(N*diag(data.matrix(x0[,-6])%*%V_d4_M4%*%t(data.matrix(x0[,-6]))))
max(N*diag(data.matrix(x0[,-6])%*%V_d5_M4%*%t(data.matrix(x0[,-6]))))
max(N*diag(data.matrix(x0[,-6])%*%V_d6_M4%*%t(data.matrix(x0[,-6]))))
max(N*diag(data.matrix(x0[,-6])%*%V_d7_M4%*%t(data.matrix(x0[,-6]))))

# under model (5): recall that the effect of Age^2 is not identifiable for d1 (i.e. design with 2 age levels)
#max(N*diag(data.matrix(x0)%*%V_d1_M5%*%t(data.matrix(x0))))
max(N*diag(data.matrix(x0)%*%V_d2_M5%*%t(data.matrix(x0))))
max(N*diag(data.matrix(x0)%*%V_d3_M5%*%t(data.matrix(x0))))
max(N*diag(data.matrix(x0)%*%V_d4_M5%*%t(data.matrix(x0))))
max(N*diag(data.matrix(x0)%*%V_d5_M5%*%t(data.matrix(x0))))
max(N*diag(data.matrix(x0)%*%V_d6_M5%*%t(data.matrix(x0))))
max(N*diag(data.matrix(x0)%*%V_d7_M5%*%t(data.matrix(x0))))
