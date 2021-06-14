# this function returns a sample of num values from a 
# GG distribution with the hardcoded parameters

get_sample <- function(num) {
    # n is the number of samples that will be returned
    suppressMessages(library(gamlss))
    suppressMessages(library(gamlss.dist))
    suppressMessages(library(gamlss.add))
    
    # the fitted model parameters
    param_mu = 5.628154
    param_sigma = -3.645491
    param_nu = -576.9269
    
    # sampling 2x because some values will be Inf and
    # we must filter them out
    sample <- rGG(n = num*2, 
                  mu = exp(param_mu), 
                  sigma = exp(param_sigma), 
                  nu = param_nu
    )
    
    # removing the spurious Inf values
    sample <- sample[sample != Inf]
    
    return(as.integer(sample[1:num]))
}
