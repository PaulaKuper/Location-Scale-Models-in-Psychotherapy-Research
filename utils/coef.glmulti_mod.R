# model averaging: done through coef
# modified code from: https://github.com/cran/glmulti/blob/master/R/glmulti.R

#These lines define a function called coef.glmulti, which takes several arguments: 
# - res (results object from glmulti package), 
# - select (specifies which models to use), 
# - varweighting (method for weighting variance), 
# - icmethod (method for calculating confidence intervals), 
# - alphaIC (significance level for confidence intervals), and additional arguments
# - models (list of fitted location-scale models)
# - formulas_as_vector (list of model formulas)

coef.glmulti <- function(res, select="all", varweighting="Buckland", icmethod="Lukacs", alphaIC=0.05, models, formulas_as_vector, ...) 
{
  
  # calculate the weights for each model based on their respective AIC values
  # assign higher weights to models with lower AIC values
  ww = exp(-(res@crits - res@crits[1])/2)
  ww = ww/sum(ww)
  
  # sort weights according to original model formulas (not in descending order)
  res_weightable= weightable(res)[, c("model", "weights")]
  
  formulas_orig = map_chr(formulas_vec, ~paste("~", .x))
  formulas_orig = sapply(formulas_orig, function(x) formula(x))

  ww = res_weightable[match(formulas_orig, res_weightable$model),]$weights
  
  # define helper function cucu to calculate cumulative weights
  cucu = function(i) sum(ww[1:i])
  wwc = lapply(1:length(ww),cucu)
  
  # determine which models to include in the model averaging procedure based on 
  # the select argument
  if (length(select)>1)
    whom=select
  else if (select=="all")
    whom=1:length(res@crits)
  else if (is.integer(select))
    whom = 1:select
  else if (is.numeric(select)) {
    if (select <=1) 
      whom = which(wwc<=select)
    else 
      whom = which(res@crits <= (res@crits[1]+select))
  }
  
  # extract the AIC values (mods) and formulas (formo) of the selected models
  mods = res@crits[whom]
  formo = res@formulas[whom]
  
  # check if the results object (res) contains the refitted models (ress).
  hasobj = try(res@ress,TRUE)
  if (!inherits(hasobj,"try-error") && length(res@ress) > 0) {
    # model ress included: no need to refit
    coffee = res@ress[whom]
  } else {
    # store fitted models in coffee.
    coffee = list()
    coffee = models
    }
  
  # construct list of coefficients
  if (length(coffee)==1) {
    # only one model ! Do conditional inference for continuity
    warning("Only one candidate: standard conditional inference was performed.")
    return(coef(coffee[[1]]))
  }
  
  # define new getfit function to extract coefficients of an object of class 'rma.ls'
  getfit2 <- function(model) {
    if (model$test == "z") {
                cbind(estimate = c(coef(model)$beta,
                                   coef(model)$alpha), 
                      se = c(sqrt(diag(vcov(model)$beta)),
                             sqrt(diag(vcov(model)$alpha))), 
                      df = Inf)
              }
              else {
                cbind(estimate = c(coef(model)$beta,
                                   coef(model)$alpha), 
                      se = c(sqrt(diag(vcov(model)$beta)),
                             sqrt(diag(vcov(model)$alpha))),
                      df = model$k -model$p)
                }
  }
  
  # apply the getfit function to each model in coffee and store the results in coke
  coke=lapply(coffee, getfit2)          
              
  # separate location and scale intercept
  coke <- lapply(coke, function(x) {
    dimnames(x)[[1]][1] <- c("intercept_l")
    return(x)
  })
  
  # extracts the names of the coefficients (namou) by concatenating the 
  # coefficient names of each model
  namou=unique(unlist(lapply(coke,function(x) dimnames(x)[[1]])))
  
  # this equates synonymous notations (e.g. x:y and y:x)
  unique(sapply(namou, function(x) {
    sort(strsplit(x, ":")[[1]])-> pieces
    if (length(pieces)>1) paste(pieces[1],":",pieces[2],"<>",pieces[2],":",pieces[1],sep="")
    else x
  }))-> codamou 
  namou = sapply(codamou, function(x)  {
    sort(strsplit(x, "<>")[[1]])-> pieces
    if (length(pieces)>1) (pieces[1])
    else x
  })
  namou2 = sapply(codamou, function(x)  {
    sort(strsplit(x, "<>")[[1]])-> pieces
    if (length(pieces)>1) pieces[2]
    else x
  })
  
  # initialize matrices coconutM, coconutSE, and coconutN with zeros. 
  # these matrices will store the coefficient values, standard errors, and 
  # sample sizes for each model.
  coconutM=matrix(0,length(formo),length(namou))
  coconutSE=matrix(0,length(formo),length(namou))
  coconutN = numeric(length(namou))
  
  # define a helper function to match coefficient names (quiqui) with the indices 
  # in namou or namou2 based on their position to get values, deviations, 
  # presence/absence and sample sizes for all models
  matchou=function(quiqui) {
    match(quiqui,namou)-> w1
    if (is.na(w1)) {
      match(quiqui,namou2)
    } else w1
  }
  
  gettou=function(i) {
    ele=coke[[i]] # estimate, se, df
    nana = dimnames(ele)			
    if (length(nana)==2) nana=nana[[1]]
    mimi=numeric(4*length(namou)) 
    if (length(nana) > 0) {
      for (k in 1:(length(nana))) {
        mimi[matchou(nana[k])]=ele[k,1]
        mimi[matchou(nana[k])+length(namou)]=ele[k,2]
        mimi[matchou(nana[k])+2*length(namou)]=1
        mimi[matchou(nana[k])+3*length(namou)]=ele[k,3]
      }
    } 
    return(mimi)
  }
  
  # calculate the model-specific coefficient values (coconutM), standard errors 
  # (coconutSE), and sample sizes (coconutN) for each selected model. 
  # use the results obtained from the gettou function and assigns them to the 
  # corresponding matrices.
  lol=sapply(lapply(1:length(coke),gettou),rbind)
  
  coconutM = matrix(unlist(t(lol[1:length(namou),])),nrow=length(whom)) # rows = models, cols= coefs
  coconutSE = matrix(unlist(t(lol[(1:length(namou))+length(namou),])),nrow=length(whom))
  
  # NA are set to zero
  coconutM[is.na(coconutM)]=0
  coconutN =  matrix(unlist(t(lol[(1:length(namou))+2*length(namou),])),nrow=length(whom))
  
  # handle degrees of freedom by taking the maximum value of each row in coconutDF
  coconutDF = matrix(unlist(t(lol[(1:length(namou))+3*length(namou),])),nrow=length(whom))
  modelsdf= unlist(apply(coconutDF,1,max))
  
  #  calculates the total number of models (nene) by multiplying the number of 
  # rows with the total sample sizes
  nene = matrix(rep(1,length(whom))%*%coconutN, ncol=1, dimnames=list( namou, c("Nb models")))
 
  # construct weight vectors by dividing the weights of the selected models by 
  # the sum of the weights of the selected models to obtain normalized weights (waou). 
  waou=ww[whom]/sum(ww[whom])
  
  # expand waou to have the same dimensions as coconutN by repeating the rows. 
  # calculate the total weights for each coefficient by multiplying the 
  # normalized weights (waou) with the total sample sizes (coconutN).
  waouv=waou
  for (i in 2:length(namou)) waouv=rbind(waouv,waou) 
  waouv= t(waouv)*coconutN
  totwaou = waou%*%coconutN
 
  # weight estimates by multiplying the normalized weights (waouv) with the 
  # model-specific coefficient values (coconutM)
  averest = matrix((rep(1,length(whom))%*%(waouv*coconutM)), ncol=1, dimnames=list( namou, c("Estimate")))
  
  # create matrix to store the total weights for each coefficient.
  weighty =  matrix(totwaou, ncol=1, dimnames=list( namou, c("Importance")))
  
  # weight variances; use different methods depending on the 'varweighting' argument
  if (varweighting=="Johnson") {
    squaredevs = waou%*%(((coconutM-t(matrix(rep(averest,length(whom)), length(namou), length(whom))))^2))
    condivars =  waou%*%((coconutSE^2))
    avervar = matrix(condivars+squaredevs, ncol=1, dimnames=list( namou, c("Uncond. variance")))
  } else if (varweighting=="Buckland") {
    squaredevs = ((coconutM-t(matrix(rep(averest,length(whom)), length(namou), length(whom)))))^2
    condivars = coconutSE^2
    avervar = matrix(((waou)%*%(sqrt(squaredevs+condivars)))^2, ncol=1, dimnames=list( namou, c("Uncond. variance")))
  } else 
    avervar = matrix(rep(NA,length(namou)), ncol=1, dimnames=list( namou, c("Uncond. variance")))
  
  # now move on to confidence intervals
  if (icmethod=="Burnham") {
    # uses Burnham & Anderson (2002) suggestion
    stuvals = (as.numeric(lapply(modelsdf,  function(x) qt(1-alphaIC/2, x)))/qnorm(1-alphaIC/2))^2
    adjsem= (matrix(rep(stuvals, length(namou)), nrow=length(whom)))*coconutSE^2
    adjsem = adjsem + (coconutM-t(matrix(rep(averest, length(whom)), ncol=length(whom))))^2
    adjse = qnorm(1-alphaIC/2)*(waou%*%sqrt(adjsem))
    uncondIC = matrix(adjse, ncol=1,  dimnames=list(namou, c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
  } else if (icmethod=="Lukacs") {
    # uses Lukacs et al. (2008) student-like method
    # get degrees of freedom for each model
    averddf = sum(waou*modelsdf)
    uncondIC = matrix(sqrt(avervar)*qt(1-alphaIC/2,averddf), ncol=1,  dimnames=list(namou, c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
  } else {
    # uses standard gaussian interval 
    uncondIC = matrix(sqrt(avervar)*qnorm(1-alphaIC/2), ncol=1,  dimnames=list(namou, c(paste("+/- (alpha=", alphaIC, ")",sep=""))))
  }
  
  averaging = cbind(averest, avervar, nene, weighty, uncondIC)
  ordonat = order(weighty[,1])
  
  # remove null component if necessary
  which(namou=="NULLOS")-> nooz
  if (length(nooz)>0) ordonat=setdiff(ordonat,nooz)
  averaging[ordonat,]
}

