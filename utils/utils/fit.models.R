fitModels <- function(formula, att_scale, data){
  M <- rma(yi=.g, sei=.g_se, scale= formula, 
           data=data, method= "REML", att= att_scale,
           test="knha", control= list(optimizer="BFGS")) 
  return(M)
}