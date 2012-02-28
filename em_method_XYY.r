/*
 *  em.h
 *  
 *
 *  Created by xyy on 12/16/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

######################################################################
####02-09-2012 update for glasso with par penalize.diagonal= FALSE####
######################################################################



setwd("/Users/xyy/Research/Yufeng_Will_project/Tissue Data/")
library(glasso)
library(Matrix)
library(MASS) 
source("correg.R")



#y1 = read.csv("adipose.matrix_ml.csv", header=T)
#y2 = read.csv("hypo.matrix_ml.csv", header=T)
#y3 = read.csv("liver.matrix_ml.csv", header=T)
#y4 = read.csv("muscle.matrix_ml.csv", header=T)


y1 <- read.csv("y1.csv", header=T)
y2 <- read.csv("y2.csv", header=T)
y3 <- read.csv("y3.csv", header=T)
y4 <- read.csv("y4.csv", header=T)

y1 = as.matrix(y1)
y2 = as.matrix(y2)
y3 = as.matrix(y3)
y4 = as.matrix(y4)


## very nice trick
for(ii in 1:4){
	EvalStatement = paste("rm(Y", ii, ")", sep = "");
	cat(EvalStatement)
	eval(parse(text = EvalStatement));
}




NormMatrix = function(X){
	p = dim(X)[1]
	A = matrix(NA, nrow = p, ncol = dim(X)[2])
	for (i in 1:p)
	{u = mean(X[i, ])
		A[i, ] = (X[i, ]-u)/(sqrt(var(X[i, ])))
	}
	A 
}



for(ii in 1:4){
	EvalStatement = paste("y", ii, "=","NormMatrix(y",ii,")", sep="");
	cat(EvalStatement)
	eval(parse(text=EvalStatement));
}



Mscale = function(X) {
	X/max(X)
} 





GetTheta <- function(x1, x2, x3, x4, Z, lambda.i=0.4, lambda.z=0.1){
	s1 = var(t(x1))	
	s2 = var(t(x2))	
	s3 = var(t(x3))	
	s4 = var(t(x4))	
	sz = var(t(Z))
	
	T1 = glasso(s1, rho = lambda.i, penalize.diagonal = FALSE, maxit = 100)$wi
	T2 = glasso(s2, rho = lambda.i, penalize.diagonal = FALSE, maxit = 100)$wi
	T3 = glasso(s3, rho = lambda.i, penalize.diagonal = FALSE, maxit = 100)$wi
	T4 = glasso(s4, rho = lambda.i, penalize.diagonal = FALSE, maxit = 100)$wi
	Tz = glasso(sz, rho = lambda.z, penalize.diagonal = FALSE, maxit = 100)$wi
	
	theta = list(theta1 = T1, theta2 = T2, theta3 = T3, theta4 = T4, thetaz = Tz)	
	
	theta
	
}


OffdiagMean=function(X){
	n = dim(X)[1]
	lID = lower.tri(X)
	uID = upper.tri(X)
	a = mean(abs(X[lID]))
	b = mean(abs(X[uID]))
	m = (a+b)/2
	m
}




########################################
###Alan way to find the index for 0#####
########################################
Get.index.g = function(X, a){
	jMat <- matrix(rep(1:NROW(X), NROW(X)), NROW(X), NROW(X));
    iMat <- t(jMat);
    index.result <- cbind(jMat[X > a], iMat[X > a])
	index.result
}

Get.index.l = function(X, a){
	jMat <- matrix(rep(1:NROW(X), NROW(X)), NROW(X), NROW(X));
	iMat <- t(jMat);
	index.result <- cbind(jMat[X <= a], iMat[X <= a])
	index.result
}











## test = EMg(y1, y2, y3, y4, lambda.i = 1.1, lambda.z = 1.1);
EMg <- function(y1, y2, y3, y4, lambda.i = 0.4, lambda.z = 0.1){
	
	Z = (y1+y2+y3+y4)/4
	theta = GetTheta(y1-Z, y2-Z, y3-Z, y4-Z,Z, lambda.i, lambda.z)

	
	count = 1
    max_iter = 100
	tol_value = 1e-2
	diff_value = 1
	diff_value1 = NA
	diff.ETz = NA
	diff.sum = NA
	n = dim(x1)[2]
	ET1 = theta[[1]]
	ET2 = theta[[2]]
	ET3 = theta[[3]]
	ET4 = theta[[4]]
	ETz = theta[[5]]
	
	
	
	while((diff_value > tol_value) & (count < max_iter))
	{
		A = ET1 + ET2 + ET3 + ET4 + ETz
		AA = solve(A)
		CC = ET1%*%y1 + ET2%*%y2 + ET3%*%y3 + ET4%*%y4
		Ez = AA%*%CC
		Ezz = AA + Ez%*%t(Ez);
			
		if(max(ET1) <= 50){
			SS1 = (y1 %*% t(y1) - y1 %*% t(Ez) - Ez %*% t(y1) + Ezz)/n}else{
				SS1 = (y1 %*% t(y1) - y1 %*% t(Ez) - Ez %*% t(y1) + Ezz)/n + diag(rep(.001,NROW(ET1)))}	
        if(max(ET2) <= 50){
			SS2 = (y2 %*% t(y2) - y2 %*% t(Ez) - Ez %*% t(y2) + Ezz)/n}else{
				SS2 = (y2 %*% t(y2) - y2 %*% t(Ez) - Ez %*% t(y2) + Ezz)/n + diag(rep(.001,NROW(ET2)))}
		if(max(ET3) <= 50){
		   SS3 = (y3 %*% t(y3) - y3 %*% t(Ez) - Ez %*% t(y3) + Ezz)/n}else{
		   SS3 = (y3 %*% t(y3) - y3 %*% t(Ez) - Ez %*% t(y3) + Ezz)/n + diag(rep(.001,NROW(ET3)))}
		   if(max(ET4) <= 50){
			  SS4 = (y4 %*% t(y4) - y4 %*% t(Ez) - Ez %*% t(y4) + Ezz)/n}else{
			  SS4 = (y4 %*% t(y4) - y4 %*% t(Ez) - Ez %*% t(y4) + Ezz)/n + diag(rep(.001,NROW(ET4)))}
		if(max(ETz) <= 50){
		   SSz = Ezz/n}else{SSz = Ezz/n + diag(rep(.001,NROW(ETz)))} 
		if (kappa(SS1) > 1e+15)
        {  
            SS1 <- SS1 + 0.001*diag(dim(SS1)[1]) 
        }
		if (kappa(SS2) > 1e+15)
        {  
            SS2 <- SS2 + 0.001*diag(dim(SS2)[1]) 
        }
		if (kappa(SS3) > 1e+15)
        {  
            SS3 <- SS3 + 0.001*diag(dim(SS3)[1]) 
        }
		if (kappa(SS4) > 1e+15)
        {  
            SS4 <- SS4 + 0.001*diag(dim(SS4)[1]) 
        }
		if (kappa(SSz) > 1e+15)
        {  
            SSz <- SSz + 0.001*diag(dim(SSz)[1]) 
        }
		
		M1 = glasso(SS1, rho=lambda.i, thr=0.001, maxit=400, penalize.diagonal= FALSE)
		M2 = glasso(SS2, rho=lambda.i, thr=0.001, maxit=400, penalize.diagonal= FALSE)
		M3 = glasso(SS3, rho=lambda.i, thr=0.001, maxit=400, penalize.diagonal= FALSE)
		M4 = glasso(SS4, rho=lambda.i, thr=0.001, maxit=400, penalize.diagonal= FALSE)
		Mz = glasso(SSz, rho=lambda.z, thr=0.001, maxit=400, penalize.diagonal= FALSE)
		
		ET1.new = (M1$wi + t(M1$wi)) / 2 
		ET2.new = (M2$wi + t(M2$wi)) / 2
		ET3.new = (M3$wi + t(M3$wi)) / 2
		ET4.new = (M4$wi + t(M4$wi)) / 2
		ETz.new = (Mz$wi + t(Mz$wi)) / 2
		
		
#diff_value1[count] = max(OffdiagMean(ET1 - ET1.new), OffdiagMean(ET2 - ET2.new),
#								 OffdiagMean(ET3 - ET3.new), OffdiagMean(ET4 - ET4.new),
#								 OffdiagMean(ETz - ETz.new))

		
		diff_value1[count] = (sum(abs(ET1 - ET1.new) + abs(ET2 - ET2.new) + abs(ET3 - ET3.new) + abs(ET4 - ET4.new) +
							  abs(ETz - ETz.new))) / sum(abs(ET1) + abs(ET2) + abs(ET3) + abs(ET4) + abs(ETz))
		
		diff.ETz[count] = sum(abs(ETz - ETz.new))
		diff.sum[count] = (sum(abs(ET1 - ET1.new) + abs(ET2 - ET2.new) + abs(ET3 - ET3.new) + abs(ET4 - ET4.new) +
						abs(ETz - ETz.new)))
		
		diff_value = diff_value1[count]
		
		ET1 = ET1.new 
		ET2 = ET2.new
		ET3 = ET3.new
		ET4 = ET4.new
		ETz = ETz.new
		
		cat(sep = "", "[",count,"]","dif=", round(diff_value, 3) ,"\n")
		
		count = count + 1
		
	}
	
	EMT = list(T1 = ET1, T2 = ET2, T3 = ET3, T4 = ET4, Tz = ETz, diff.ETz = diff.ETz, diff.sum = diff.sum)
	
	
####### Filter the noise learned from Ji Zhu's code
   #for ( k in seq(1,4))
  #{ ome = EMT[[1]]
    
#}
	EMT
}


