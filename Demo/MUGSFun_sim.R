get_embed = function(mysvd, d=2000, normalize=TRUE){
  # Function for getting embedding from svd
  # d: dim of the final embedding
  # mysvd: the (managed) svd result (adding an element with 'names')
  id = which(sign(mysvd$u[1,])==sign(mysvd$v[1,]))
  id = id[1:min(d,length(id))]
  embed = mysvd$u[,id]%*%diag(sqrt(mysvd$d[id]))
  if(normalize){
    embed = embed/apply(embed,1,norm,'2')
  }
  rownames(embed) = mysvd$names
  return(embed)
}


DataGen_rare_group <- function(seed, p, n1, n2, n.common, n.group, sigma.eps.1, sigma.eps.2, ratio.delta, network.k, rho.beta, rho.U0, rho.delta, sigma.rare){
  #### Function used to generate SPPMIs, dummy matrices based on prior group structures, and code-code pairs for tuning and evaluation
  set.seed(seed)
  N = n1 + n2 - n.common
  n1.no = n1 - n.common
  n2.no = n2 - n.common
  #### Group effect
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  M.beta.1 <- ar1_cor(3,rho.beta)
  M.beta.2 <- bdiag(M.beta.1, diag(7))
  Sigma.beta <- bdiag(replicate(40, M.beta.2 , simplify=FALSE))
  beta <- mvrnorm(p, rep(0,n.group), Sigma.beta)
  beta <- t(beta)
  #### Code effect
  M.U0.1 <- ar1_cor(6,rho.U0)
  Sigma.U0.1 <- bdiag(replicate(N/6, M.U0.1 , simplify=FALSE))
  M.U0.2 <- matrix(0.05, 6, 6)
  diag(M.U0.2) <- 0
  Sigma.U0.2.temp <-  bdiag(replicate(300/6, M.U0.2 , simplify=FALSE))
  Sigma.U0.2 <- matrix(0, N, N)
  Sigma.U0.2[501:800, 501:800] <- as.matrix(Sigma.U0.2.temp)
  Sigma.U0 <- Sigma.U0.1 + Sigma.U0.2
  u0 <- rmvnorm.rcpp(p, rep(0,N), as.matrix(Sigma.U0))
  u0 <- t(u0)
  temp <- seq(501, 900, 2)
  group.only.codes <- rep(0, length(temp)*group.size)
  for (i in 1: length(temp)){
    group.only.codes[(group.size*(i-1)+1):(group.size*i) ] <- seq(temp[i], 2500, n.group)
  }
  u0[group.only.codes,] =0
  #### code-site effect
  Sigma.delta.0 <- ar1_cor(network.k, rho.delta)
  Sigma.delta.1 <- bdiag(replicate(n1*ratio.delta/network.k, Sigma.delta.0, simplify=FALSE))
  Sigma.delta.2 <- bdiag(replicate(n2*ratio.delta/network.k, Sigma.delta.0, simplify=FALSE))
  delta1.temp <- mvrnorm(p, rep(0,n1*ratio.delta), Sigma.delta.1)
  delta2.temp <- mvrnorm(p, rep(0,n2*ratio.delta), Sigma.delta.2)
  set.seed(1) 
  idx1 <- sample(n.common, ncol(delta1.temp)) + 1000
  idx2 <- sample(n.common, ncol(delta2.temp)) 
  delta1 <- matrix(0, n1, p)
  delta2 <- matrix(0, n2, p)
  delta1[as.vector(idx1),] <-  t(delta1.temp) 
  delta2[as.vector(idx2),] <-  t(delta2.temp) 
  #### Code-Code pairs for tuning and evaluation
  ## similar pairs 
  name.beta.full.1 <- matrix(501:2500, 400, 5)
  name.beta.full <-c(t(name.beta.full.1 ))
  supp.cov.beta <- bdiag(replicate(40, bdiag(matrix(1, 3*5, 3*5), bdiag(replicate(7, matrix(1, 1*5, 1*5), simplify=FALSE))), simplify=FALSE))
  supp.cov.beta.full <- matrix(0, N, N)
  supp.cov.beta.full[name.beta.full, name.beta.full] <- as.matrix(supp.cov.beta)
  supp.cov.beta.full.upper <- supp.cov.beta.full
  supp.cov.beta.full.upper[lower.tri(supp.cov.beta.full.upper, diag = TRUE)] <- 0
  pairs.sim <- which(supp.cov.beta.full.upper!=0, arr.ind = T)
  pairs.sim <- as.data.frame(pairs.sim)
  pairs.sim$type <- rep('similarity', length(pairs.sim$row))
  pairs.sim$row <- as.character(pairs.sim$row)
  pairs.sim$col <- as.character(pairs.sim$col)
  ## related pairs
  supp.cov.u0 <- as.matrix(Sigma.U0)
  supp.cov.u0[as.matrix(Sigma.U0)!=0] <- 1
  supp.cov.u0.upper <- supp.cov.u0
  supp.cov.u0.upper[lower.tri(supp.cov.u0.upper, diag = TRUE)] <- 0
  pairs.rel.new <- matrix(0,n2*ratio.delta/network.k*(network.k*(network.k-1)/2),2)
  for (j in 1: (n2*ratio.delta/network.k)){
    pairs.rel.new[((network.k*(network.k-1)/2)*(j-1) + 1) :(j*network.k*(network.k-1)/2),] <- t(combn(idx2[((j-1)*network.k+1):(j*network.k)],2))
  }
  pairs.rel.new <- as.data.frame(pairs.rel.new+ (n1 - (n1-n.common)))
  supp.cov.delta.2 <- matrix(0, n2, n2)
  supp.cov.delta.2[cbind(as.vector((pairs.rel.new-(n1 - (n1-n.common)))[,1]), as.vector((pairs.rel.new-(n1 - (n1-n.common)))[,2]))] <- 1
  supp.cov.delta.2 <- supp.cov.delta.2 + t(supp.cov.delta.2 )
  supp.cov.delta.2.upper <- supp.cov.delta.2 
  supp.cov.delta.2.upper[lower.tri(supp.cov.delta.2.upper , diag = TRUE)] <- 0
  n.rel <- sum(supp.cov.u0.upper!=0) + n2*ratio.delta/network.k*(network.k*(network.k-1)/2)  
  pairs.rel.shared <- which(supp.cov.u0.upper!=0, arr.ind = T)
  pairs.rel.2 <- which(supp.cov.delta.2.upper!=0, arr.ind = T) + 1000
  ## Merge similar and related pairs
  pairs.rel <- rbind(pairs.rel.shared, pairs.rel.2) 
  pairs.rel <- as.data.frame(pairs.rel)
  pairs.rel$type <- rep('related', dim(pairs.rel)[1])
  pairs.rel$row <- as.character(pairs.rel$row)
  pairs.rel$col <- as.character(pairs.rel$col)
  pairs.rel <- pairs.rel[!(duplicated(pairs.rel)),]
  pairs.rel.full <- rbind(pairs.sim, pairs.rel)
  pairs.rel.full$type <- 'related'
  set.seed(1)
  idx.rel <- sample(1:dim( pairs.rel.full)[1])
  pairs.rel.CV <-  pairs.rel.full[idx.rel[1:floor(dim( pairs.rel.full)[1]/2) ], ]
  pairs.rel.EV <-  pairs.rel.full[idx.rel[(floor(dim(pairs.rel)[1]/2)+ 1):dim( pairs.rel.full)[1]], ]
  idx.rare <- as.numeric(names(sort(table(c(pairs.rel.EV[,1], pairs.rel.EV[,2])), decreasing = T)))
  idx.rare <- (intersect(idx.rare, 1001:3000))[(1+100):(n.rare+100)]
  #### Embeddings
  beta.full.1 = do.call(rbind, replicate(5, beta, simplify=FALSE))
  dim(beta.full.1)
  beta.full <- rbind(rbind(matrix(0, 500,p), beta.full.1), matrix(0,500,p))
  # site 1 
  u.1 = beta.full[1:n1,] + u0[1:n1,] + delta1
  ## site 2
  u.2 = beta.full[(n1 - n1.no +1):N,] + u0[(n1 - n1.no +1):N,] +  delta2 
  ## True S matrices
  S.1.0 <- u.1%*%t(u.1)
  S.2.0 <- u.2%*%t(u.2)
  S.miss.0 <- u.1[1:n1.no, ]%*%t(u.2[(n.common+1):n2, ])
  ## Add noise
  set.seed(seed)
  err.1 <- matrix(rnorm(n1^2, 0, sigma.eps.1), n1, n1)
  S.1 <- S.1.0 + err.1
  #### Add noises for rare codes in target
  idx.freq <- setdiff(1001:3000, idx.rare)
  n.rare <- length(idx.rare)
  n.freq <- length(idx.freq)
  err.rare.1 <- matrix(rnorm(n.rare^2, 0, sigma.rare), n.rare, n.rare)
  err.rare.2 <- matrix(rnorm(n.rare*n.freq, 0, sqrt(sigma.rare*sigma.eps.2)), n.rare, n.freq)
  err.freq <- matrix(rnorm(n.freq^2, 0, sigma.eps.2), n.freq, n.freq)
  err.2.temp <- matrix(0, n2, n2)
  err.2.temp[1:n.rare, 1:n.rare] <- err.rare.1
  err.2.temp[1:n.rare, (1+n.rare):n2] <- err.rare.2
  err.2.temp[(1+n.rare):n2, 1:n.rare] <- t(err.rare.2)
  err.2.temp[(1+n.rare):n2, (1+n.rare):n2] <- err.freq
  err.2 <- err.2.temp[match(1001:3000, c(idx.rare, idx.freq)), match(1001:3000, c(idx.rare, idx.freq))]
  S.2 <- S.2.0 + err.2
  S.1[lower.tri(S.1)] = t(S.1)[lower.tri(S.1)]
  S.2[lower.tri(S.2)] = t(S.2)[lower.tri(S.2)]
  rownames(S.1) <- as.character(1:n1)
  colnames(S.1) <- as.character(1:n1)
  rownames(S.2) <- as.character((n1.no +1):N)
  colnames(S.2) <- as.character((n1.no +1):N)
  ######## Dummy matrices based on prior group structures ########
  group.names <- as.character(1:n.group)
  group.names.full <- rep(group.names, group.size)
  X.group.0 <- dummy_cols(group.names.full)
  X.group.1 <- X.group.0[,-1]
  X.group <- rbind(rbind(matrix(0,500,n.group), as.matrix(X.group.1)),matrix(0,500,n.group))
  X.group.source <- X.group[1:n1, ]
  X.group.target <- X.group[1001:N, ]
  rownames(X.group.source) <-  as.character(1:n1)
  rownames(X.group.target) <- as.character(1001:N)
  output <- list('delta1' = delta1, 'delta2' = delta2, 'u.1' = u.1, 'u.2' = u.2, 'S.1' = S.1, 'S.2' = S.2, 
                 'S.1.0' = S.1.0, 'S.2.0' = S.2.0, 'X.group.source' = X.group.source, 'X.group.target' = X.group.target, 
                 'pairs.rel.CV' = pairs.rel.CV, 'pairs.rel.EV' = pairs.rel.EV)
  return(output)
}


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

evaluation.sim <- function(pairs.rel, U){
  #### Function used for tuning and evaluation 
  names <- rownames(U)
  pairs <- data.frame(matrix(0, length(pairs.rel$row), 2) )
  colnames(pairs) <- c('id1', 'id2')
  pairs$id1 <- match(pairs.rel$row,names)
  pairs$id2 <- match(pairs.rel$col, names)
  pairs <- na.omit(pairs)
  n.rel.neg <- nrow(pairs) 
  pairs.temp <- pairs
  pairs$id1 <- apply(pairs.temp, 1, min)
  pairs$id2 <- apply(pairs.temp, 1, max)
  
  names.rel.1 <- (pairs.rel[,1])[(pairs.rel[,1]%in%names)&(pairs.rel[,2]%in%names)]
  names.rel.2 <-(pairs.rel[,2])[(pairs.rel[,1]%in%names)&(pairs.rel[,2]%in%names)]
  names.rel <- union(names.rel.1, names.rel.2) 
  set.seed(1)
  pairs.rel.neg = data.frame(id1 = sample(which(names%in%names.rel), 3*n.rel.neg, replace = TRUE),
                             id2 = sample(which(names%in%names.rel), 3*n.rel.neg, replace = TRUE))
  pairs.rel.neg = subset(pairs.rel.neg, id1!=id2)
  pairs.rel.neg_temp = pairs.rel.neg
  pairs.rel.neg$id1 = apply(pairs.rel.neg_temp, 1, min)
  pairs.rel.neg$id2 = apply(pairs.rel.neg_temp, 1, max)
  pairs.rel.neg = pairs.rel.neg[!(duplicated(pairs.rel.neg)), ]
  pairs.rel.neg = anti_join(pairs.rel.neg, pairs, by = c("id1", "id2")) #19273     2
  
  if(nrow(pairs.rel.neg)<n.rel.neg){
    cat("Warning! Something Wrong Here! Too little Null sample!\n")
    return(0)
  }
  
  pairs.rel.pos <- pairs
  pairs.rel.neg <- pairs.rel.neg[sample(nrow(pairs.rel.neg),n.rel.neg,replace = FALSE), c('id1','id2')]
  y = c(rep(1, nrow(pairs.rel.pos)), rep(0, nrow(pairs.rel.neg)))
  id1 = pairs.rel.pos$id1; id2 = pairs.rel.pos$id2; num = length(id1)
  p1 = unlist(lapply(1:num, function(i){
    return(sum(U[id1[i],]*U[id2[i],]))
  }))
  id1 = pairs.rel.neg$id1; id2 = pairs.rel.neg$id2; num = length(id1)
  p2 = unlist(lapply(1:num, function(i){
    return(sum(U[id1[i],]*U[id2[i],]))
  }))
  p = c(p1,p2)
  AUC.rel = auc(y,p)
  output <- list("n.rel" = n.rel.neg, "AUC.rel" = AUC.rel)
  return(output)
}







