library(glmnet)
library(MASS)
library(fastDummies)
library(lsa)
library(pROC)
library(tictoc)
library(rsvd)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(foreach)
library(doSNOW)
library(pryr)
library(grplasso)
library(dplyr)

################################ Preparation for the main MUGS alg ################################
############ Parameters w.r.t tuning & parallel computing 
TUNE = F
Lambda <- c(2500)
Lambda.delta <- c(7000)
K1 <- length(Lambda)
K2 <- length(Lambda.delta)
n.core <- 10
tol <- 1
############  Create arrays to store results
CV_overall_BCH <- array(0, dim = c(K1, K2, 2))
CV_overall_MGB <- array(0, dim = c(K1, K2, 2))
delta.spst.ovl <- matrix(0,K1, K2)
############ Read in PPMI matrices and initial embeddings from MGB and BCH
load('Emb.BCH.rollup.500.Rdata')
load('Emb.MGB.rollup.scale.500.Rdata')
load('PPMI.MGB.rollup.Rdata')
load('PPMI.BCH.rollup.Rdata')
p <- dim(V.BCH)[2]
n1 <- dim(V.MGB)[1]
n2 <- dim(V.BCH)[1]
S.1 <- PPMI.MGB[sort(rownames(PPMI.MGB)), sort(colnames(PPMI.MGB))]
S.2 <- PPMI.BCH[sort(rownames(PPMI.BCH)), sort(colnames(PPMI.BCH))]
S.1 <- S.1/1.8
############ Align two sets of initial embeddings 
common_codes <- intersect(rownames(V.BCH), rownames(V.MGB))
common_codes <- sort(common_codes)
n.common<-length(common_codes) 
idx_MGB <- rownames(V.MGB)%in%common_codes
idx_BCH <- rownames(V.BCH)%in%common_codes
comm.embed.MGB <- V.MGB[idx_MGB,]
comm.embed.MGB <- comm.embed.MGB[sort(rownames(comm.embed.MGB), index.return=T)$ix,]
comm.embed.BCH <- V.BCH[idx_BCH,]
comm.embed.BCH <- comm.embed.BCH[sort(rownames(comm.embed.BCH), index.return=T)$ix,]
#### orthogonal Procrustes problem based on overlapped parts
A <- t(comm.embed.MGB)%*%comm.embed.BCH
eigenA<-eigen(t(A)%*%A + 0.0000000001*diag(p), symmetric = T)
W <- A%*%eigenA$vectors%*%diag(eigenA$values^(-1/2), p, p)%*%t(eigenA$vectors)
embed.MGB.ot<- V.MGB%*%W
temp <- sort(rownames(embed.MGB.ot), index.return=T)
embed.MGB.init <- embed.MGB.ot[temp$ix,]
temp <- sort(rownames(V.BCH), index.return=T)
embed.BCH.init <- V.BCH[temp$ix,]
############ Group code based on hierarchical ontology
#### Group PheCodes 
## MGB
num.MGB <-sapply(strsplit(rownames(embed.MGB.init),':',fixed=TRUE), "[",2)
group.MGB <-sapply(strsplit(rownames(embed.MGB.init),':',fixed=TRUE), "[",1)
embed.phecode.MGB <- embed.MGB.init[group.MGB == "PheCode",]
int.MGB<-sapply(strsplit(num.MGB[group.MGB == "PheCode"],'.',fixed=TRUE), "[",1)
## BCH
num.BCH <-sapply(strsplit(rownames(embed.BCH.init),':',fixed=TRUE), "[",2)
group.BCH <-sapply(strsplit(rownames(embed.BCH.init),':',fixed=TRUE), "[",1)
embed.phecode.BCH <- embed.BCH.init[group.BCH == "PheCode",]
int.BCH<-sapply(strsplit(num.BCH[group.BCH == "PheCode"],'.',fixed=TRUE), "[",1)
Phecode.MGB.group.1 <- names(table(as.factor(int.MGB))[table(as.factor(int.MGB))==1])
Phecode.MGB.group.2 <- names(table(as.factor(int.MGB))[table(as.factor(int.MGB))>1])
Phecode.BCH.group.1 <- names(table(as.factor(int.BCH))[table(as.factor(int.BCH))==1])
Phecode.BCH.group.2 <- names(table(as.factor(int.BCH))[table(as.factor(int.BCH))>1])
group.Phecode.2 <- union(Phecode.MGB.group.2, Phecode.BCH.group.2)
n.group.Phecode.2 <- length(group.Phecode.2) 
group.Phecode.no.1 <- intersect(Phecode.MGB.group.1, Phecode.BCH.group.1)
n.group.Phecode <-n.group.Phecode.2
group.Phecode  <- group.Phecode.2 
MGB.group.Phecode<- union(Phecode.MGB.group.2, Phecode.MGB.group.1[Phecode.MGB.group.1 %in% Phecode.BCH.group.2])
BCH.group.Phecode <- union(Phecode.BCH.group.2, Phecode.BCH.group.1[Phecode.BCH.group.1 %in% Phecode.MGB.group.2])
#### Group LOINC
load("LOINC.Hierarchy.Rdata")
## MGB
embed.LOINC.MGB <- embed.MGB.init[!(group.MGB%in% c("CCS", "PheCode", 'RXNORM')),]
LOINC.MGB.names <- as.data.frame(rownames(embed.LOINC.MGB))
colnames(LOINC.MGB.names) <- 'CODE'
LOINC.new$CODE <- LOINC.new$LOINC
LOINC.new$CODE[LOINC.new$rollup==T] <- LOINC.new$PARENT_LOINC[LOINC.new$rollup==T]
LOINC.new.1 <-  LOINC.new[,c(8, 4)]
LOINC.new.2 <- LOINC.new.1[!duplicated(LOINC.new.1$CODE),]
LOINC.MGB.group <-  merge(LOINC.MGB.names, LOINC.new.2,  by = 'CODE', all.x = T)
MGB.unmapped.LOINC.new <- LOINC.MGB.group[is.na(LOINC.MGB.group$group),1] # 241 MGB local LOINCs
LOINC.MGB.group$group[is.na(LOINC.MGB.group$group)] <- paste('MGB', LOINC.MGB.group$CODE[is.na(LOINC.MGB.group$group)])
## BCH
embed.LOINC.BCH <- embed.BCH.init[!(group.BCH%in% c("CCS", "PheCode", 'RXNORM')),]
LOINC.BCH.names <- as.data.frame(rownames(embed.LOINC.BCH))
colnames(LOINC.BCH.names) <- 'CODE'
LOINC.BCH.group <-  merge(LOINC.BCH.names, LOINC.new.2,  by = 'CODE', all.x = T)
LOINC.MGB.group.1 <- names(table(as.factor(LOINC.MGB.group$group))[table(as.factor(LOINC.MGB.group$group))==1])
LOINC.MGB.group.2 <- names(table(as.factor(LOINC.MGB.group$group))[table(as.factor(LOINC.MGB.group$group))>1])
LOINC.BCH.group.1 <- names(table(as.factor(LOINC.BCH.group$group))[table(as.factor(LOINC.BCH.group$group))==1])
LOINC.BCH.group.2 <- names(table(as.factor(LOINC.BCH.group$group))[table(as.factor(LOINC.BCH.group$group))>1])
group.LOINC.2 <- union(LOINC.MGB.group.2, LOINC.BCH.group.2)
group.LOINC.no.1 <- intersect(LOINC.MGB.group.1, LOINC.BCH.group.1)
group.LOINC.1 <- group.LOINC.no.1[LOINC.BCH.group$CODE[LOINC.BCH.group$group%in%group.LOINC.no.1]!=LOINC.MGB.group$CODE[LOINC.MGB.group$group%in%group.LOINC.no.1 ]]
group.LOINC <- union(group.LOINC.2, group.LOINC.1)
n.group.LOINC <- length(group.LOINC) 
MGB.group.LOINC <- union(union(LOINC.MGB.group.2, LOINC.MGB.group.1[LOINC.MGB.group.1 %in% LOINC.BCH.group.2]), LOINC.MGB.group.1[LOINC.MGB.group.1 %in%group.LOINC.1])
BCH.group.LOINC <- union(union(LOINC.BCH.group.2, LOINC.BCH.group.1[LOINC.BCH.group.1 %in% LOINC.MGB.group.2]), LOINC.BCH.group.1[LOINC.BCH.group.1 %in%group.LOINC.1])
#### Group RXNORM 
RXNORM.Hierarchy <-  read.csv("RXNORM_ING_ATC_version2_28MAR2023vid.csv")
RXNORM.Hierarchy.l4 <- RXNORM.Hierarchy[,c(1,11)]
## MGB
embed.RXNORM.MGB <- embed.MGB.init[group.MGB == "RXNORM",]
num.RXNORM.MGB <- num.MGB[group.MGB == "RXNORM"]
RXNORM.MGB.names <- as.data.frame(num.RXNORM.MGB)
colnames(RXNORM.MGB.names) <- 'RXCUI'
RXNORM.MGB.group <-  merge(RXNORM.MGB.names, RXNORM.Hierarchy.l4,  by = 'RXCUI', all.x = T)
UnmappedRxNorm_MGB <- RXNORM.MGB.group$RXCUI[RXNORM.MGB.group$LEVEL4=='']
RXNORM.MGB.group$LEVEL4[RXNORM.MGB.group$LEVEL4 == ''] <- paste('MGB', RXNORM.MGB.group$RXCUI[RXNORM.MGB.group$LEVEL4==''])
## BCH
embed.RXNORM.BCH <- embed.BCH.init[group.BCH == "RXNORM",]
num.RXNORM.BCH <- num.BCH[group.BCH == "RXNORM"]
RXNORM.BCH.names <- as.data.frame(num.RXNORM.BCH)
colnames(RXNORM.BCH.names) <- 'RXCUI'
RXNORM.BCH.group <-  merge(RXNORM.BCH.names, RXNORM.Hierarchy.l4,  by = 'RXCUI', all.x = T)
UnmappedRxNorm_BCH <- RXNORM.BCH.group$RXCUI[RXNORM.BCH.group$LEVEL4=='']
RXNORM.BCH.group$LEVEL4[RXNORM.BCH.group$LEVEL4 == ''] <- paste('BCH', RXNORM.BCH.group$RXCUI[RXNORM.BCH.group$LEVEL4==''])
RXNORM.MGB.group.1 <- names(table(as.factor(RXNORM.MGB.group$LEVEL4))[table(as.factor(RXNORM.MGB.group$LEVEL4))==1])
RXNORM.MGB.group.2 <- names(table(as.factor(RXNORM.MGB.group$LEVEL4))[table(as.factor(RXNORM.MGB.group$LEVEL4))>1])
RXNORM.BCH.group.1 <- names(table(as.factor(RXNORM.BCH.group$LEVEL4))[table(as.factor(RXNORM.BCH.group$LEVEL4))==1])
RXNORM.BCH.group.2 <- names(table(as.factor(RXNORM.BCH.group$LEVEL4))[table(as.factor(RXNORM.BCH.group$LEVEL4))>1])
group.RXNORM.2 <- union(RXNORM.MGB.group.2, RXNORM.BCH.group.2)
group.RXNORM.no.1 <- intersect(RXNORM.MGB.group.1, RXNORM.BCH.group.1)
group.RXNORM.1 <- group.RXNORM.no.1[RXNORM.BCH.group$RXCUI[RXNORM.BCH.group$LEVEL4%in%group.RXNORM.no.1 ]!=RXNORM.MGB.group$RXCUI[RXNORM.MGB.group$LEVEL4%in%group.RXNORM.no.1 ]]
group.RXNORM <- union(group.RXNORM.2, group.RXNORM.1)
n.group.RXNORM <- length(group.RXNORM) 
MGB.group.RXNORM <- union(union(RXNORM.MGB.group.2, RXNORM.MGB.group.1[RXNORM.MGB.group.1 %in% RXNORM.BCH.group.2]),  RXNORM.MGB.group.1[RXNORM.MGB.group.1 %in% group.RXNORM.1 ])
BCH.group.RXNORM <- union(union(RXNORM.BCH.group.2, RXNORM.BCH.group.1[RXNORM.BCH.group.1 %in% RXNORM.MGB.group.2]),  RXNORM.BCH.group.1[RXNORM.BCH.group.1 %in% group.RXNORM.1 ])


################################ Initialization for the main MUGS alg ################################
############  Calculate initial group effects 
n.group = n.group.Phecode + n.group.LOINC + n.group.RXNORM # 1341
embed.CCS.MGB <- embed.MGB.init[group.MGB == "CCS",]
embed.CCS.BCH <- embed.BCH.init[group.BCH == "CCS",]
group.names.MGB <- c(rownames(embed.CCS.MGB), LOINC.MGB.group$group, int.MGB, RXNORM.MGB.group$LEVEL4)
group.names.BCH <- c(rownames(embed.CCS.BCH), LOINC.BCH.group$group, int.BCH, RXNORM.BCH.group$LEVEL4)
group.names <- c(group.names.MGB, group.names.BCH) 
X.group<- fastDummies::dummy_cols(group.names)
X.group1 <- X.group[,-1]
colnames.X.group1 <- sapply(strsplit(colnames(X.group1),'_',fixed=TRUE), "[",2)
beta.names <- c(group.LOINC, group.Phecode, group.RXNORM)
X.group2 <- X.group1[,colnames.X.group1%in%beta.names]
beta.names.1 <- colnames(X.group2)
beta.names.2 <- sapply(strsplit(beta.names.1 ,'_',fixed=TRUE), "[",2)
X.group.MGB.full <- X.group2[1: n1,]
X.group.BCH.full <- X.group2[(n1 + 1):(n1 + n2),]

embed.MGB <- rbind(embed.CCS.MGB, embed.LOINC.MGB, embed.phecode.MGB, embed.RXNORM.MGB)
embed.BCH <- rbind(embed.CCS.BCH, embed.LOINC.BCH, embed.phecode.BCH, embed.RXNORM.BCH)
embed <- rbind(embed.MGB, embed.BCH)
code.names.MGB <- rownames(embed.MGB)
code.names.BCH <- rownames(embed.BCH)
code.names <- rownames(embed)

beta.int <- matrix(0, n.group, p)
rownames(beta.int) <- beta.names.1
names.in.group.MGB <- matrix(NA, n.group, max(colSums(X.group.MGB.full)))
names.in.group.BCH <- matrix(NA, n.group, max(colSums(X.group.BCH.full)))
for (i in 1:n.group){
  if (sum(X.group.MGB.full[,i]==T)>0){
    names.in.group.MGB[i, 1:sum(X.group.MGB.full[,i]==T)] <- code.names.MGB[X.group.MGB.full[,i]==T]
  }
  if (sum(X.group.BCH.full[,i]==T)>0){
    names.in.group.BCH[i, 1:sum(X.group.BCH.full[,i]==T)] <- code.names.BCH[X.group.BCH.full[,i]==T]
  }
  beta.int[i, ] <- colMeans(embed[X.group2[,i]==T ,]) 
}

############  Calculate initial code effect 
n.code <- n1 + n2 - n.common
## MGB
n.nogroup.MGB <-length(setdiff(code.names.MGB, names.in.group.MGB)) 
names.nogroup.MGB <- setdiff(code.names.MGB, names.in.group.MGB)
n.group.MGB <- sum(!(is.na(names.in.group.MGB))) 
names.group.MGB <- names.in.group.MGB[!(is.na(names.in.group.MGB))]
names.nogroup.ovl.MGB <- intersect(names.nogroup.MGB, common_codes)
n.nogroup.ovl.MGB <- length(names.nogroup.ovl.MGB ) 
names.nogroup.novl.MGB <- setdiff(names.nogroup.MGB, names.nogroup.ovl.MGB )
n.nogroup.novl.MGB <- length(names.nogroup.novl.MGB ) 
names.group.ovl.MGB <- intersect(names.group.MGB, common_codes)
n.group.ovl.MGB <-length(names.group.ovl.MGB) 
names.group.novl.MGB <- setdiff(names.group.MGB, names.group.ovl.MGB)
n.group.novl.MGB <-length(names.group.novl.MGB) 
## BCH
n.nogroup.BCH <-length(setdiff(code.names.BCH, names.in.group.BCH))
names.nogroup.BCH <- setdiff(code.names.BCH, names.in.group.BCH)
n.group.BCH <- sum(!(is.na(names.in.group.BCH)))
names.group.BCH <- names.in.group.BCH[!(is.na(names.in.group.BCH))]
names.nogroup.ovl.BCH <- intersect(names.nogroup.BCH, common_codes)
n.nogroup.ovl.BCH <- length(names.nogroup.ovl.BCH )
names.nogroup.novl.BCH <- setdiff(names.nogroup.BCH, names.nogroup.ovl.BCH)
n.nogroup.novl.BCH <- length(names.nogroup.novl.BCH) 
names.group.ovl.BCH <- intersect(names.group.BCH, common_codes)
n.group.ovl.BCH <-length(names.group.ovl.BCH) 
names.group.novl.BCH <- setdiff(names.group.BCH, names.group.ovl.BCH)
n.group.novl.BCH <-length(names.group.novl.BCH)
#### Codes belonging to a group and present at MGB and BCH
zeta.group.ovl <- matrix(0,n.group.ovl.BCH,p)
zeta.group.ovl.temp <- (embed.MGB[rownames(embed.MGB)%in%names.group.ovl.MGB,] + embed.BCH[rownames(embed.BCH)%in%names.group.ovl.BCH,])/2
rownames(zeta.group.ovl) <- rownames(zeta.group.ovl.temp)
for (i in 1:n.group.ovl.BCH){
  group.name <- group.names.BCH[code.names.BCH==rownames(zeta.group.ovl.temp)[i]] 
  zeta.group.ovl[i, ] <- zeta.group.ovl.temp[i,] - beta.int[beta.names.2==group.name,]
}
#### Codes not belonging to any groups and present at MGB and BCH
zeta.nogroup.ovl <- (embed.MGB[rownames(embed.MGB)%in%names.nogroup.ovl.MGB,] + embed.BCH[rownames(embed.BCH)%in%names.nogroup.ovl.BCH,])/2
#### Codes belonging to a group but only present at MGB or BCH
zeta.group.novl.MGB <- matrix(0,n.group.novl.MGB,p)
zeta.group.novl.BCH <- matrix(0,n.group.novl.BCH,p)
zeta.group.novl.MGB.temp <- embed.MGB[rownames(embed.MGB)%in%names.group.novl.MGB,]
zeta.group.novl.BCH.temp <- embed.BCH[rownames(embed.BCH)%in%names.group.novl.BCH,]
rownames(zeta.group.novl.MGB) <- rownames(zeta.group.novl.MGB.temp)
rownames(zeta.group.novl.BCH) <- rownames(zeta.group.novl.BCH.temp)
for (i in 1:n.group.novl.MGB){
  group.name <- group.names.MGB[code.names.MGB==rownames(zeta.group.novl.MGB.temp)[i]] 
  zeta.group.novl.MGB[i, ] <- zeta.group.novl.MGB.temp [i,] - beta.int[beta.names.2==group.name,]
}
for (i in 1:n.group.novl.BCH){
  group.name <- group.names.BCH[code.names.BCH==rownames(zeta.group.novl.BCH.temp)[i]] 
  zeta.group.novl.BCH[i, ] <- zeta.group.novl.BCH.temp [i,] - beta.int[beta.names.2==group.name,]
}
#### Codes not belonging to any groups and only present at MGB or BCH
zeta.nogroup.novl.MGB <- embed.MGB[rownames(embed.MGB)%in%names.nogroup.novl.MGB,]
zeta.nogroup.novl.BCH <- embed.BCH[rownames(embed.BCH)%in%names.nogroup.novl.BCH,]
zeta.MGB.temp <- rbind(zeta.group.ovl, zeta.nogroup.ovl, zeta.group.novl.MGB, zeta.nogroup.novl.MGB)
zeta.MGB <- zeta.MGB.temp[match(code.names.MGB, rownames(zeta.MGB.temp)),]
zeta.BCH.temp <- rbind(zeta.group.ovl, zeta.nogroup.ovl, zeta.group.novl.BCH, zeta.nogroup.novl.BCH)
zeta.BCH <- zeta.BCH.temp[match(code.names.BCH, rownames(zeta.BCH.temp)),]
zeta.int <- rbind(zeta.MGB, zeta.BCH)
zeta.ovl.temp <- rbind(zeta.group.ovl, zeta.nogroup.ovl)
zeta.int.ovl <- zeta.ovl.temp[match(common_codes, rownames(zeta.ovl.temp)), ]
zeta.novl.MGB.temp <- rbind(zeta.group.novl.MGB, zeta.nogroup.novl.MGB)
zeta.int.novl.MGB <- zeta.novl.MGB.temp[sort(rownames(zeta.novl.MGB.temp)),]
zeta.novl.BCH.temp <- rbind(zeta.group.novl.BCH, zeta.nogroup.novl.BCH)
zeta.int.novl.BCH <- zeta.novl.BCH.temp[sort(rownames(zeta.novl.BCH.temp)),]

############ Calculate initial code-site effect 
delta.int <- embed - as.matrix(X.group2)%*%beta.int - zeta.int
## By default that delta is zero for codes only present at MGB or BCH


################################ Main MUGS Alg ################################
###### Load in key functions of MUGS alg
source('MUGSFun.R')
###### Load in functions for validation the embeddings 
source('Embed_Eval_Pediatric.R')
###### Silver-standard labels for 
AllRelationPairs <- readRDS("new.CV.rds")


tic('main loop')
for (k1 in 1:K1) {
  cat('\n k1=', k1)
  lambda = Lambda[k1]
  for (k2 in 1:K2) {
    cat('\n k2=', k2)
    lambda.delta = Lambda.delta[k2]
    #### initialization
    U.1 <- embed.MGB
    U.2 <- embed.BCH
    V.1 <- U.1 
    V.2 <- U.2 
    dif_loss <- 10
    delta.int.V <- delta.int
    zeta.int.V <- zeta.int 
    beta.int.V <-beta.int
    delta.int.U <- delta.int
    zeta.int.U <-zeta.int 
    beta.int.U <-beta.int
    while (abs(dif_loss) > tol){
      #### updating V_j
      loss_main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2) 
      
      loss_zeta_V_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.V[1:n1, ])[rownames(zeta.int.V[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
      loss_zeta_V_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.V[1:n1, ])[!(rownames(zeta.int.V[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
      loss_zeta_V_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.V[(1+n1):(n1+n2), ])[!(rownames(zeta.int.V[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
      loss_zeta_V <- loss_zeta_V_ovl + loss_zeta_V_1 + loss_zeta_V_2
      loss_zeta_U_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.U[1:n1, ])[rownames(zeta.int.U[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
      loss_zeta_U_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.U[1:n1, ])[!(rownames(zeta.int.U[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
      loss_zeta_U_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.U[(1+n1):(n1+n2), ])[!(rownames(zeta.int.U[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
      loss_zeta_U <- loss_zeta_U_ovl + loss_zeta_U_1 + loss_zeta_U_2
      
      loss_delta_V.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.V[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) 
      loss_delta_V.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.V[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2) 
      loss_delta_V <- loss_delta_V.1 + loss_delta_V.2 
      loss_delta_U.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.U[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) 
      loss_delta_U.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.U[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2) 
      loss_delta_U <- loss_delta_U.1 + loss_delta_U.2 
    
      loss_p_V <- loss_zeta_V +  loss_delta_V
      loss_p_U <- loss_zeta_U +  loss_delta_U
      loss <- loss_main +  loss_p_V + loss_p_U
      cat('\n loss:',loss)
      loss.V  <- loss 
      dif_loss_V <-10
      while (abs(dif_loss_V) > tol)
      {
        ## Update group effects 
        tic('group')
        GroupEff.out <- GroupEff_par(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, X.group.MGB.full, X.group.BCH.full, n.group,code.names, beta.int.V, 0, p, n.core)
        toc()
        beta.int.V <- as.matrix(GroupEff.out$beta)
        V.1 <- as.matrix(GroupEff.out$V.MGB.new)
        V.2 <- as.matrix(GroupEff.out$V.BCH.new)
        temp.group <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2) 
        cat('\n diff_groupEff=', temp.group - loss_main)
        ## Update code effects
        CodeEff.out <-CodeEff_Matrix(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int.V, lambda, p)
        zeta.int.V <- as.matrix(CodeEff.out$zeta)
        V.1 <- as.matrix(CodeEff.out$V.1.new)
        V.2 <-as.matrix(CodeEff.out$V.2.new)
        temp.code.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2) 
        temp.code.p.ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.V[1:n1, ])[rownames(zeta.int.V[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2) 
        temp.code.p.1 <-lambda*sqrt(log(p)/n1)*norm((zeta.int.V[1:n1, ])[!(rownames(zeta.int.V[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        temp.code.p.2 <-lambda*sqrt(log(p)/n2)*norm((zeta.int.V[(1+n1):(n1+n2), ])[!(rownames(zeta.int.V[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)  
        cat('\n diff_CodeEff=', temp.code.main-temp.group+temp.code.p.ovl+temp.code.p.1+temp.code.p.2-loss_zeta_V)
        ## Update code-site effects
        tic('CodeSiteEff_l2')
        CodeSiteEff.out <- CodeSiteEff_l2_par(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, delta.int.V, lambda.delta, p, common_codes, n.common, n.core)
        toc()
        delta.int.V <- as.matrix(CodeSiteEff.out$delta) 
        V.1 <- as.matrix(CodeSiteEff.out$V.1.new)
        V.2 <- as.matrix(CodeSiteEff.out$V.2.new)
        temp.codesite.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
        temp.codesite.p <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.V[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) +lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.V[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2) 
        cat('\n diff_CodeSiteEff=', temp.codesite.main - temp.code.main + temp.codesite.p -loss_delta_V)

        loss_main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
        loss_zeta_V_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.V[1:n1, ])[rownames(zeta.int.V[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_V_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.V[1:n1, ])[!(rownames(zeta.int.V[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_V_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.V[(1+n1):(n1+n2), ])[!(rownames(zeta.int.V[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_V <- loss_zeta_V_ovl + loss_zeta_V_1 + loss_zeta_V_2
        
        loss_delta_V.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.V[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) 
        loss_delta_V.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.V[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2) 

        loss_delta_V <- loss_delta_V.1 + loss_delta_V.2 
        loss_p_V <- loss_zeta_V +  loss_delta_V
        loss.new.V <- loss_main + loss_p_V + loss_p_U 
        dif_loss_V <- (loss.new.V - loss.V)/loss.V 
        cat('\n dif_loss_V',dif_loss_V)
        loss.V <- loss.new.V
      }
      ## Update U_i
      loss.U <- loss.new.V
      dif_loss_U <- 10
      while (abs(dif_loss_U) > 1)
      {
        ## Update group effect
        GroupEff.out <- GroupEff_par(t(S.1), t(S.2), n1, n2, V.1, V.2, U.1, U.2, X.group.MGB.full, X.group.BCH.full, n.group, code.names, beta.int.U, 0, p, n.core)
        beta.int.U <- as.matrix(GroupEff.out$beta)
        U.1 <- as.matrix(GroupEff.out$V.MGB.new)
        U.2 <- as.matrix(GroupEff.out$V.BCH.new)
        temp.group <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2) 
        cat('\n diff_groupEff=', temp.group - loss_main)
        
        ## Update code effect
        CodeEff.out <-CodeEff_Matrix(t(S.1), t(S.2), n1, n2, V.1, V.2, U.1, U.2, common_codes, zeta.int.U, lambda, p)
        zeta.int.U <- as.matrix(CodeEff.out$zeta)
        U.1 <- as.matrix(CodeEff.out$V.1.new)
        U.2 <-as.matrix(CodeEff.out$V.2.new)
        temp.code.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2) 
        temp.code.p.ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.U[1:n1, ])[rownames(zeta.int.U[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2) 
        temp.code.p.1 <-lambda*sqrt(log(p)/n1)*norm((zeta.int.U[1:n1, ])[!(rownames(zeta.int.U[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        temp.code.p.2 <-lambda*sqrt(log(p)/n2)*norm((zeta.int.U[(1+n1):(n1+n2), ])[!(rownames(zeta.int.U[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)  
        cat('\n diff_codeEff=', temp.code.main-temp.group+temp.code.p.ovl+temp.code.p.1+temp.code.p.2-loss_zeta_U)
        
        ## Update code-site effect
        CodeSiteEff.out <- CodeSiteEff_l2_par(t(S.1), t(S.2), n1, n2, V.1, V.2, U.1, U.2, delta.int.U, lambda.delta, p, common_codes, n.common, n.core)
        delta.int.U <- as.matrix(CodeSiteEff.out$delta) 
        U.1 <- as.matrix(CodeSiteEff.out$V.1.new)
        U.2 <- as.matrix(CodeSiteEff.out$V.2.new)
        temp.codesite.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
        temp.codesite.p <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.U[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) +lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.U[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2) 
        cat('\n diff_CodeSiteEff=', temp.codesite.main - temp.code.main + temp.codesite.p - loss_delta_U)
        
        loss_main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
        loss_zeta_U_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.U[1:n1, ])[rownames(zeta.int.U[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_U_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.U[1:n1, ])[!(rownames(zeta.int.U[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_U_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.U[(1+n1):(n1+n2), ])[!(rownames(zeta.int.U[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_U <- loss_zeta_U_ovl + loss_zeta_U_1 + loss_zeta_U_2

        loss_delta_U.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.U[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) 
        loss_delta_U.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.U[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2) 

        loss_delta_U <- loss_delta_U.1 + loss_delta_U.2 
        
        loss_p_U <- loss_zeta_U + loss_delta_U
        loss.new.U <- loss_main + loss_p_V + loss_p_U
        dif_loss_U <- (loss.new.U - loss.U)/loss.U 
        cat('\n dif_loss_U',dif_loss_U)
        loss.U <- loss.new.U
      }
      dif_loss <-  (loss.new.U - loss)/loss
      cat('\n dif_loss=', dif_loss) 
    } 
    delta.spst.ovl[k1, k2]  <- sum(rowSums(abs(delta.int.U[rownames(delta.int.U)%in%common_codes,]))!=0)
    #### validation 
    ans = Evaluate(U.2[rowSums(U.2)!=0 ,]/apply(U.2[rowSums(U.2)!=0 ,],1,norm,'2')) 
    n.label.pp <- (ans$fulltable$`codi-codi`[c(1:6,8,41,42),])[4, 1] 
    n.label.pr <- (ans$fulltable$`codi-codi`[c(1:6,8,41,42),])[5, 1] 
    CV_overall_BCH[k1,k2,] <- round(ans$fulltable$`codi-codi`[c(1:6,8,41,42),],2)[4:5, 2] 
    ans = Evaluate(U.1[rowSums(U.1)!=0 ,]/apply(U.1[rowSums(U.1)!=0 ,],1,norm,'2')) 
    CV_overall_MGB[k1,k2,] <- round(ans$fulltable$`codi-codi`[c(1:6,8,41,42),],2)[4:5, 2] 
  }
}
toc()
#### Find the results with optimal lambdas ####
if (TUNE==T){
  ## Weighted average 
  CV_overall_BCH.wavr <- (n.label.pp*CV_overall_BCH[, , 1] + n.label.pr*CV_overall_BCH[, , 2])/(n.label.pp + n.label.pr)
  idx <- which(CV_overall_BCH.wavr == max(CV_overall_BCH.wavr), arr.ind = TRUE)
  lambda.opt <- Lambda[idx[1]]
  lambda.delta.opt <- Lambda.delta[idx[2]]
  #### Output
  cat('\n lambda1=', lambda.opt)
  cat('\n lambda2=', lambda.delta.opt)
  ## Sparsity of code-site effect
  cat('\n delta.spst.ovl=', delta.spst.ovl)
}


if (TUNE == F){
  ######## Output ########
  #### Embedding matrices and the cosine similarity matrices 
  save(U.1, file = 'U.MGB.Rdata')
  save(U.2, file = 'U.BCH.Rdata')
  CS.1 <- (U.1/apply(U.1,1,norm,'2'))%*%t(U.1/apply(U.1,1,norm,'2'))
  CS.2 <- (U.2/apply(U.2,1,norm,'2'))%*%t(U.2/apply(U.2,1,norm,'2'))
  save(CS.1, file = 'CS.MGB.Rdata')
  save(CS.2, file = 'CS.BCH.Rdata')
  
  #### Group effects 
  rownames(beta.int.U) <- beta.names.2
  save(beta.int.U, file = 'beta.Rdata')
  
  #### Similar codes between sites vs dissimilar codes between sites
  similar.codes <- rownames(delta.int.U)[rowSums(abs(delta.int.U))==0]
  save(similar.codes, file='similar.codes.Rdata')
  dissimilar.codes <- rownames(delta.int.U)[rowSums(abs(delta.int.U))!=0]
  save(dissimilar.codes, file='dissimilar.codes.Rdata')
  #### Sparsity of code-site effect
  cat('\n delta.spst.ovl=', delta.spst.ovl)
}



