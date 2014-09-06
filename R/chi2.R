#' Computes asymptotic and/or parametric bootstrap \eqn{\chi^2} test of a general function of direct and indirect effects
#' 
#' @param model a \link{lavaan} \link{model.syntax}. The \link{model.syntax} may not include constraints as defined in \code{constraints}
#' @param data a data frame  containing the observed data from the sample.
#' @param constraints a general function of the model parameters stated in terms of (in)equality constraint. 
#' @param alpha significance level. Default is .05.
#' @param boot if \code{boot="none"} (default), asymptotic \eqn{\chi^2} test is performed. If \code{boot="parametric"}, parametric \eqn{\chi^2} test is produced. 
#' @param R number of bootstrap samples.
#' @param parallel The type of parallel operation to be used (if any). If missing, the default is "no".
#' @param ncpus Integer: number of processes to be used in parallel operation. Default is 2L.
#' @param ... other options. See \link{sem}.
#' @return a \link{list} of the following objects:
#' \item{LRT}{Asymptotic likelihood ratio test. One can also extract an asymptotic p value here. See \link{anova} for more detail.}
#' \item{p}{bootstrap p value for the \eqn{\chi^2}. This occurs only when \code{boot="parametric"}}
#' \item{fitH0}{Model estimate for the more restricted model. This is a \link{lavaan} object. One can extract all the information from this using \link{lavaan} functions.}
#' @examples
#' data(PoliticalDemocracy)
#' m1 <- '
#' ind60 =~ x1 + x2 + x3
#' dem60 =~ y1 + y2 + y3 + y4
#' dem65 =~ y5 + y6 + y7 + y8
#' dem60 ~ a*ind60
#' dem65 ~ cp*ind60 + b* dem60
#' y1 ~~ y5
#' y2 ~~ y4 + y6
#' y3 ~~ y7
#' y4 ~~ y8
#' y6 ~~ y8
#' '
#' # Fit the mode without estimation
#' ( chi2Test <- chi2(m1, data = PoliticalDemocracy ,"a*b==0", boot="parametric", R=20 , parallel = "multicore", ncpus=2) )


chi2 <- function(model, data, constraints = "a*b==0", alpha=.05, boot="none", R=200, parallel="none", ncpus=2L, ...){
  
  fm0 <- sem(model, data, constraints = constraints, ...) #more restricted model
  fm1 <- sem(model, data,  ... ) #less restricted model
  
  Saturated <- fitMeasures(fm1,"df")==0 #Is model saturated?
  
  if(Saturated){
    LRT <- anova(fm0) #as.list(fitMeasures(fm0,c("chisq","df","pvalue") ) )
  }
  else LRT <- lavTestLRT(fm1,fm0)
  
  if(boot=="parametric"){ 
    if(!Saturated)  { # if not saturated
      bootParamRes <- bootstrapLRT(fm0,fm1,R=R, type='parametric', return.LRT=TRUE, ncpus=ncpus, parallel=parallel)
      pboot <- bootParamRes[1] 
    }
    else{ #if saturated
      Tboot <- bootstrapLavaan(fm0, R=R, type='parametric', FUN=fitMeasures, fit.measures="chisq")
      T0 <- fitMeasures(fm0,"chisq") # chi2 for the less restricted model
      # compute a bootstrap based p-value
      pboot <- length(which(Tboot > T0))/length(Tboot)
    }
    
  }
  
  if(boot=="parametric"){
    return(list(LRT=LRT,p=pboot,Nullfit=fm0))
  }
  else return(list(LRT=LRT,fitH0=fm0))
}
