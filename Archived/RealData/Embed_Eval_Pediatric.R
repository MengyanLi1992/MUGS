Cos_id_boot = function(pairs, dict, g1, g2, type = 1){
  # type: 1 or 2: 
  ## type 1: sample from pairs; type 2: sample from dict
  set.seed(1)
  pairs$id1 = match(pairs$feature_id1, dict$feature_id)
  pairs$id2 = match(pairs$feature_id2, dict$feature_id)
  pairs = na.omit(pairs)
  num = nrow(pairs)
  if(num==0) return(0)
  if(type==1){
    if(length(g1)==length(g2) & sum(g1==g2)==length(g1)){
      idg = unique(c(unique(pairs$id1), unique(pairs$id2)))
      pairs_null = data.frame(id1 = sample(idg, 3*num, replace = TRUE),
                              id2 = sample(idg, 3*num, replace = TRUE))
    }else{
      pairs_null = data.frame(id1 = sample(unique(pairs$id1), 3*num, replace = TRUE),
                              id2 = sample(unique(pairs$id2), 3*num, replace = TRUE))
    }
  }else if(type==2){
    pairs_null = data.frame(id1 = sample(which(dict$group==g1), 3*num, replace = TRUE),
                            id2 = sample(which(dict$group==g2), 3*num, replace = TRUE))
  }
  pairs_full = data.frame(id1 = c(pairs$id1, pairs$id2), 
                          id2 = c(pairs$id2, pairs$id1))
  pairs_null = setdiff(pairs_null, pairs_full)
  rm(pairs_full)
  pairs_null = subset(pairs_null, id1!=id2)
  if(nrow(pairs_null)<num){
    cat("Warning! Something Wrong Here! Too little Null sample!\n")
    return(0)
  }
  return(list(pairs_h = pairs[,c('id1','id2')], 
              pairs_null = pairs_null[1:num,c('id1','id2')]))
}
ROC_id = function(embed, pairs, dict, g1, g2, type = 1, savestore = TRUE){
  Boot = Cos_id_boot(pairs = pairs, dict=dict, g1 = g1, g2 = g2, 
                     type = type)
  if(length(Boot)==1){
    return(c('#pairs'=0,'auc'=NA,'cut/0.01'=NA,'cut/0.05'=NA,'cut/0.1'=NA,
             'TPR/0.01'=NA,'TPR/0.05'=NA,'TPR/0.1'=NA))
  }
  y = c(rep(1, nrow(Boot$pairs_h)), rep(0, nrow(Boot$pairs_null)))
  if(savestore){
    id1 = Boot$pairs_h$id1; id2 = Boot$pairs_h$id2; num = length(id1)
    p1 = unlist(lapply(1:num, function(i){
      return(sum(embed[id1[i],]*embed[id2[i],]))
    }))
    id1 = Boot$pairs_null$id1; id2 = Boot$pairs_null$id2; num = length(id1)
    p2 = unlist(lapply(1:num, function(i){
      return(sum(embed[id1[i],]*embed[id2[i],]))
    }))
    p = c(p1,p2)
  }else{
    Cos = embed%*%t(embed)
    p = c(Cos[cbind(Boot$pairs_h$id1,Boot$pairs_h$id2)],
          Cos[cbind(Boot$pairs_null$id1,Boot$pairs_null$id2)])
  }
  
  roc0 = roc(y, p)
  cut0 = as.numeric(quantile(p2,0.99))
  cut1 = as.numeric(quantile(p2,0.95))
  cut2 = as.numeric(quantile(p2,0.9))
  
  tpr0 = roc0$sensitivities[which(roc0$thresholds>cut0)[1]]
  tpr1 = roc0$sensitivities[which(roc0$thresholds>cut1)[1]]
  tpr2 = roc0$sensitivities[which(roc0$thresholds>cut2)[1]]
  return(c('#pairs'=num,'auc'=roc0$auc,'cut/0.01'=cut0[1],'cut/0.05'=cut1[1],'cut/0.1'=cut2[1],
           'TPR/0.01'=tpr0,'TPR/0.05'=tpr1,'TPR/0.1'=tpr2))
}
Get_subset = function(embed, dict, pairs, g, clust){
  dict = dict[match(rownames(embed),dict$feature_id),]
  if("codi"%in%g){
    idx = which(dict$group%in%g|dict$group!="CUI")
    embed = embed[idx,]
    dict = dict[idx,]
    dict$group[which(dict$group!="CUI")] = "codi"
  }else{
    idx = which(dict$group%in%g)
    embed = embed[idx,]
    dict = dict[idx,]
  }
  if(clust == 1){
    pairs = data.frame(feature_id1 = pairs$CUI1,
                       feature_id2 = pairs$CUI2)
  }else if(clust == 2){
    pairs = data.frame(feature_id1 = c(pairs$CUI1,pairs$CUI2),
                       feature_id2 = c(pairs$code2,pairs$code1))
  }else if(clust == 3){
    pairs = data.frame(feature_id1 = pairs$code1,
                       feature_id2 = pairs$code2)
  }
  pairs = pairs[!duplicated(pairs),]
  pairs = pairs[which(pairs$feature_id1!=pairs$feature_id2),]
  return(list(embed = embed, pairs = pairs, dict = dict))
}
evalu_part = function(embed, dict, AllRelationPairs, k, clust){
  # clust: 1: CUI-CUI; 2: CUI-codi; 3: codi-codi
  cuinamelist = c("may cause","inverse_isa","belong(s) to the category of",
                  "isa","dose_form_of","has_dose_form","ddx",
                  "has_active_ingredient","may be caused by",
                  "belongs to the drug family of","active_ingredient_of",
                  "procedure_diagnoses","is a subtype of","lab_diagnoses",
                  "is a risk factor for","interacts with",
                  "may contraindicate","see also",
                  "has_basis_of_strength_substance",
                  "basis_of_strength_substance_of","is associated with",
                  "possibly_equivalent_to","differential diagnosis","symptoms",
                  "complications","may be allelic with",
                  "precise_active_ingredient_of","has_precise_active_ingredient",
                  "risk factors","same_as","is a vector for","causes",
                  "belongs to drug super-family","is_modification_of",
                  "has_modification","is a category subset of","replaces",
                  "replaced_by")
  if(k==1){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="PheCode-PheCode" &
                                     AllRelationPairs$source == "PheCode-PheCode Hierachy"),]
    g1 = ifelse(clust<=2,"CUI","PheCode")
    g2 = ifelse(clust<=1,"CUI","PheCode")
    name = "PheCode-PheCode(sim)"
  }else if(k==2){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="RXNORM-RXNORM" &
                                     AllRelationPairs$source == "SNOMED-CT"),]
    g1 = ifelse(clust<=2,"CUI","RXNORM")
    g2 = ifelse(clust<=1,"CUI","RXNORM")
    name = "RXNORM-RXNORM(sim)"
  }else if(k==3){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="LAB-LAB" &
                                     AllRelationPairs$source == "Manual annotated"),]
    g1 = ifelse(clust<=2,"CUI","LOINC")
    g2 = ifelse(clust<=1,"CUI","LOINC")
    name = "LAB-LAB(sim)"
  }else if(k==4){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="PheCode-PheCode" &
                                     AllRelationPairs$source%in%c("Wikipedia", "BCH_Children_Hospital", "Cinci_Children_Hospital", "up_to_date", "Survey")),]
    g1 = ifelse(clust<=2,"CUI","PheCode")
    g2 = ifelse(clust<=1,"CUI","PheCode")
    name = "PheCode-PheCode"
  }else if(k==5){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="PheCode-RXNORM" &
                                     AllRelationPairs$source%in%c("MEDRT","SNOMED-CT","drug.com", "BCH_Children_Hospital", "Cinci_Children_Hospital", "up_to_date", "Survey")),]
    if(clust==1){
      g1 = g2 = "CUI"
    }else if(clust==2){
      g1 = "CUI"; g2 = c("PheCode","RXNORM")
    }else if(clust==3){
      g1 = "PheCode"; g2 = "RXNORM"
    }
    name = "PheCode-RXNORM"
  }else if(k==6){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="PheCode-CCS"),]
    if(clust==1){
      g1 = g2 = "CUI"
    }else if(clust==2){
      g1 = "CUI"; g2 = c("PheCode","CCS")
    }else if(clust==3){
      g1 = "PheCode"; g2 = "CCS"
    }
    name = "PheCode-CCS"
  }else if(k==7){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="LOINC-LOINC"),]
    g1 = ifelse(clust<=2,"CUI","LOINC")
    g2 = ifelse(clust<=1,"CUI","LOINC")
    name = "LOINC-LOINC"
  }else if(k==8){
    pairs = AllRelationPairs[which(AllRelationPairs$pair=="PheCode-LAB"),]
    if(clust==1){
      g1 = g2 = "CUI"
    }else if(clust==2){
      g1 = "CUI"; g2 = c("PheCode","LOINC")
    }else if(clust==3){
      g1 = "PheCode"; g2 = "LOINC"
    }
    name = "PheCode-LAB"
  }else if(k>=9 & k<=46){
    pairs = AllRelationPairs[which(AllRelationPairs$RELA==cuinamelist[k-8]),]
    g1 = ifelse(clust<=2,"CUI","codi")
    g2 = ifelse(clust<=1,"CUI","codi")
    name = cuinamelist[k-8]
  }
  edp = Get_subset(embed = embed, dict = dict, pairs = pairs, 
                   g = c(g1,g2), clust = clust)
  myauc = ROC_id(embed = edp$embed, pairs = edp$pairs, dict = edp$dict, 
                 g1 = g1, g2 = g2, type = ifelse(clust<=2,1,2))
  return(list(ans = myauc, name = name, 
              type = ifelse(k<=3,"similar","related")))
}

Evaluate = function(embed){
  dict = data.frame(feature_id = rownames(embed))
  dict$group = "CUI"
  idx = grep(":",dict$feature_id)
  dict$group[idx] = sapply(dict$feature_id[idx], function(s) strsplit(s,":")[[1]][1])
  fulltable = lapply(1:3, function(clust){
    anslist = lapply(1:40, function(k){
      evalu_part(embed, dict, AllRelationPairs, k, clust)
    })
    ans = do.call("rbind",lapply(anslist, function(x) x$ans))
    ans = as.data.frame(ans)
    rownames(ans) = unlist(lapply(anslist, function(x) x$name))
    pairtype = unlist(lapply(anslist, function(x) x$type))
    id = which(pairtype=="similar")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,1]))
    ans[nrow(ans),1] = sum(ans[id,1])
    rownames(ans)[nrow(ans)] = "weighted.sim"
    id = which(pairtype=="related")
    ans = rbind(ans, apply(ans[id,],2,weighted.mean,ans[id,1]))
    ans[nrow(ans),1] = sum(ans[id,1])
    rownames(ans)[nrow(ans)] = "weighted.rela"
    return(ans)
  })
  anstable = do.call(cbind, lapply(fulltable,function(x) x[,1:2]))
  colnames(anstable) = c("CUI-CUI","auc","CUI-codi","auc","codi-codi","auc")
  names(fulltable) = c("CUI-CUI","CUI-codi","codi-codi")
  anstable = na.omit(anstable)
  return(list(anstable = anstable, fulltable = fulltable))
}



