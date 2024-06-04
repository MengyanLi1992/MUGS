GroupEff_par <- function(S.MGB, S.BCH, n.MGB, n.BCH, U.MGB, U.BCH, V.MGB, V.BCH, X.MGB.group, X.BCH.group, n.group, name.list, beta.int, lambda,p, n.core){
  ####  function used to estimate group effects parallelly
  n.j.MGB <- apply(X.MGB.group>0, 2, sum)
  n.j.BCH <- apply(X.BCH.group>0, 2, sum)
  BCH.only <- names(which(n.j.MGB == 0))
  MGB.only <- names(which(n.j.BCH == 0))
  n.BCH.only <- length(BCH.only) 
  n.MGB.only <- length(MGB.only)
  BETA.BCHonly <- matrix(0, p, n.BCH.only)
  colnames(BETA.BCHonly) <- BCH.only
  BETA.MGBonly <- matrix(0, p, n.MGB.only )
  colnames(BETA.MGBonly) <- MGB.only
  if (n.BCH.only >0){
    for (j in 1:n.BCH.only){
      nn <- n.j.BCH[names(n.j.BCH)==BCH.only[j]]
      temp.BCH <- matrix(0,nn,p)
      for (jj in 1: nn){
        temp.BCH[jj,] <- as.vector((V.BCH[X.BCH.group[, colnames(X.BCH.group)==BCH.only[j] ]>0,])[jj,] - beta.int[rownames(beta.int)==BCH.only[j],] )
      }
      name.BCH <- (name.list[(1+n.MGB):(n.MGB + n.BCH)])[X.BCH.group[,colnames(X.BCH.group)==BCH.only[j]]>0]
      Y.BCH.1 <-as.vector(S.BCH[,rownames(S.BCH)%in%name.BCH])
      Y.BCH <- Y.BCH.1 - as.vector(as.matrix(U.BCH)%*%t(temp.BCH)) 
      V.j.BCH <- do.call(rbind, replicate(nn, U.BCH, simplify=FALSE))
      m <- glmnet(V.j.BCH, Y.BCH, alpha = 0, family = 'gaussian', lambda = lambda*2/length(Y.BCH), intercept=F)
      BETA.BCHonly[,j] <- as.vector(m$beta)
    }
  }
  if (n.MGB.only >0){
    for (j in 1:n.MGB.only){
      nn <- n.j.MGB[names(n.j.MGB)==MGB.only[j]]
      temp.MGB <- matrix(0,nn,p)
      for (jj in 1: nn){
        temp.MGB[jj,] <- as.vector((V.MGB[X.MGB.group[, colnames(X.MGB.group)==MGB.only[j] ]>0,])[jj,] - beta.int[rownames(beta.int)==MGB.only[j],] )
      }
      name.MGB <- (name.list[1:n.MGB])[X.MGB.group[,colnames(X.MGB.group)==MGB.only[j]]>0]
      Y.MGB.1 <- as.vector(S.MGB[,rownames(S.MGB)%in%name.MGB])
      Y.MGB <- Y.MGB.1 - as.vector(as.matrix(U.MGB)%*%t(temp.MGB)) 
      V.j.MGB <- do.call(rbind, replicate(nn, U.MGB, simplify=FALSE))
      m <- glmnet(V.j.MGB, Y.MGB, alpha = 0, family = 'gaussian', lambda = lambda*2/length(Y.MGB), intercept=F)
      BETA.MGBonly[,j] <- as.vector(m$beta)
    }
  }
  cl <- makeCluster(n.core, type="SOCK") 
  registerDoSNOW (cl)
  n.both <- n.group - n.BCH.only - n.MGB.only
  name.both <- setdiff(colnames(X.MGB.group),c(BCH.only, MGB.only))
  BETA = foreach(j = 1:n.both, .packages = c("glmnet"), .combine = 'cbind') %dopar%{
    name.MGB <- (name.list[1:n.MGB])[X.MGB.group[, colnames(X.MGB.group)==name.both[j]]>0]
    name.BCH <- (name.list[(1+n.MGB):(n.MGB + n.BCH)])[X.BCH.group[,colnames(X.BCH.group)==name.both[j]]>0]
    temp.MGB <- V.MGB[X.MGB.group[,colnames(X.MGB.group)==name.both[j]]>0,] - do.call(rbind, replicate(length(name.MGB),  beta.int[rownames(beta.int)==name.both[j],], simplify=FALSE))
    Y.MGB <- as.vector(S.MGB[,rownames(S.MGB)%in%name.MGB]) - as.vector(as.matrix(U.MGB)%*%t(temp.MGB))
    temp.BCH<- V.BCH[X.BCH.group[,colnames(X.BCH.group)==name.both[j]]>0,] - do.call(rbind, replicate(length(name.BCH),  beta.int[rownames(beta.int)==name.both[j],], simplify=FALSE))
    Y.BCH <- as.vector(S.BCH[,rownames(S.BCH)%in%name.BCH]) - as.vector(as.matrix(U.BCH)%*%t(temp.BCH)) 
    Y <- c(Y.MGB, Y.BCH)
    V.j <- rbind(do.call(rbind, replicate(length(name.MGB), U.MGB, simplify=FALSE)), do.call(rbind, replicate(length(name.BCH), U.BCH, simplify=FALSE)))
    m <- glmnet(V.j, Y, alpha = 0, family = 'gaussian', lambda = lambda*2/length(Y), intercept=F)
    return(m$beta)
  }  
  stopCluster(cl)
  colnames(BETA) <- name.both
  BETA.1 <- t(cbind(BETA.BCHonly, BETA.MGBonly, BETA))
  beta <- BETA.1[match(rownames(beta.int),rownames(BETA.1)),]
  dif_F <- norm(beta- beta.int, type = "F")^2/(p*n.group) 
  X.group <- rbind(X.MGB.group,X.BCH.group)
  V.new <- rbind(V.MGB, V.BCH) - as.matrix(X.group)%*%(beta.int- beta)
  V.MGB.new  <- V.new[1:n.MGB,]
  V.BCH.new  <- V.new[(n.MGB+1):(n.MGB + n.BCH),]
  output <- list("beta" = beta, "dif_F" = dif_F,'V.MGB.new' = V.MGB.new, 'V.BCH.new' = V.BCH.new)
  return(output)
}

CodeEff_Matrix<- function(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda, p){
  #### function used to estimate code effects 
  zeta.int.1 <- zeta.int[1:n1,]
  zeta.int.2 <- zeta.int[(1+n1):(n1+n2),]
  Y.1 <- S.1 - U.1%*%t(V.1 -  zeta.int.1)
  Y.2 <- S.2 - U.2%*%t(V.2 -  zeta.int.2)
  zeta.1.only <- t(solve(t(U.1)%*%U.1 + lambda*sqrt(log(p)/n1)*diag(p))%*%t(U.1)%*%Y.1[,!(colnames(Y.1)%in%common_codes)])
  zeta.2.only <- t(solve(t(U.2)%*%U.2 + lambda*sqrt(log(p)/n2)*diag(p))%*%t(U.2)%*%Y.2[,!(colnames(Y.2)%in%common_codes)])
  X <- rbind(U.1, U.2)
  Y <- rbind(Y.1[,colnames(Y.1)%in%common_codes], Y.2[,colnames(Y.2)%in%common_codes])
  zeta.common <- t(solve(t(X)%*%X + lambda*sqrt(log(p)/(n1+n2))*diag(p))%*%t(X)%*%Y)
  zeta.1.temp <- rbind(zeta.1.only, zeta.common)
  zeta.1 <- zeta.1.temp[match(rownames(U.1),rownames(zeta.1.temp)),]
  zeta.2.temp <- rbind(zeta.2.only, zeta.common)
  zeta.2 <- zeta.2.temp[match(rownames(U.2),rownames(zeta.2.temp)),]
  zeta <- rbind(zeta.1, zeta.2)
  
  dif_F <- norm(zeta- zeta.int, type = "F")^2/(p*(n1+n2))
  V.new <- rbind(V.1, V.2) - zeta.int + zeta
  V.1.new  <- V.new[1:n1,]
  V.2.new  <- V.new[(n1+1):(n1 + n2),]
  output <- list("zeta" = zeta, "dif_F" = dif_F,'V.1.new' = V.1.new, 'V.2.new' = V.2.new)
  return(output)
}


CodeSiteEff_l2_par <- function(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, delta.int, lambda.delta, p, common_codes, n.common, n.core){
  #### function used to estimate code-site effects parallelly
  delta.int.1 <- delta.int[1:n1,]
  Y.1 <- S.1 - U.1 %*% t(V.1 - delta.int.1)
  delta.1.common <- matrix(0, n.common, p)
  Y.1.common <- Y.1[,colnames(Y.1)%in%common_codes]
  cl <- makeCluster(n.core, type="SOCK") 
  registerDoSNOW (cl)
  DELTA.1 = foreach(i = 1:n.common, .packages = c("grplasso"), .combine = 'cbind') %dopar%{
    fit <- grplasso(U.1, y =  Y.1.common[,i], index = rep(1, p), lambda = lambda.delta*sqrt(log(p)/n1)/sqrt(p), model = LinReg(),
                    penscale = sqrt, center = FALSE, standardize = FALSE)
    return(fit$coefficients)
    }
  stopCluster(cl)
  delta.1.common <- t(as.matrix(DELTA.1))
  rownames(delta.1.common) <- rownames(delta.int.1)[rownames(delta.int.1)%in%common_codes]
  delta.1.temp <- rbind(delta.1.common, delta.int.1[!(rownames(delta.int.1)%in%common_codes),])
  delta.1 <- delta.1.temp[match(rownames(delta.int.1), rownames(delta.1.temp)),]
  
  ## BCH
  delta.int.2 <- delta.int[(1+n1):(n1 + n2),]
  Y.2 <- S.2 - U.2%*% t(V.2 -delta.int.2)
  delta.2.common <- matrix(0, n.common, p)
  Y.2.common <- Y.2[,colnames(Y.2)%in%common_codes]
  cl <- makeCluster(n.core, type="SOCK") 
  registerDoSNOW (cl)
  DELTA.2= foreach(i = 1:n.common, .packages = c("grplasso"), .combine = 'cbind') %dopar%{
    fit <- grplasso(U.2, y =  Y.2.common[,i], index = rep(1, p), lambda = lambda.delta*sqrt(log(p)/n2)/sqrt(p), model = LinReg(),
                    penscale = sqrt, center = FALSE, standardize = FALSE)
    return(fit$coefficients)
    }
  stopCluster(cl)
  delta.2.common <- t(as.matrix(DELTA.2))
  rownames(delta.2.common) <- rownames(delta.int.2)[rownames(delta.int.2)%in%common_codes]
  delta.2.temp <- rbind(delta.2.common, delta.int.2[!(rownames(delta.int.2)%in%common_codes),])
  delta.2 <- delta.2.temp[match(rownames(delta.int.2), rownames(delta.2.temp)),]
  
  delta <- rbind(delta.1, delta.2)
  V.1.new <- V.1 - delta.int.1 + delta.1
  V.2.new <- V.2 - delta.int.2 + delta.2
  output <- list("delta" = delta, 'V.1.new' = V.1.new, 'V.2.new' = V.2.new)
  return(output)
}







