#' @title EM Clustering
#'
#' @description Preform EM clustering on a dataset.
#'
#' @usage EM(x, k)
#'
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix
#' (such as a numeric vector or a data frame with all numeric columns).
#' @param k the number of initial means
#'
#' @return A list containing cluster means, clustering vector, and covariance.
#'
#' @import tidyverse
#' @import mvtnorm
#'
#' @export

em_clust = function(x, k, show_cluster_probs = F){
    priors = rep(1/k, k)
    x = as.matrix(x)
    if(ncol(x) > 1){
        x = data.frame(x)
        means = sample_n(x, k)
    }else {
        means = sample(x, k)}
    x = as.matrix(x)
    cov.var = vector('list', length = k)
    for(i in 1:k){
        if(
            ncol(x) > 1){cov.var[[i]] = cov(x)
        }else {
            cov.var[[i]] = sd(x)}
    }
    probs2 = 0
    probs1 = 1
    while(sum((probs2 - probs1)^2) > 0.00001){
        probs1 = probs2
        probs2 = NULL
        if(ncol(x) > 1){
            for(i in 1:k){
                initprob = dmvnorm(x, as.numeric(means[i,]), cov.var[[i]])
                probs2 = cbind(probs2, initprob)}
        } else {
            for(i in 1:k){
                initprob = dnorm(x, means[i], cov.var[[i]])
                probs2 = cbind(probs2, initprob)
            }
        }
        totalpost = 0
        for(i in 1:k){
            tempP = priors[i]*probs2[,i]
            totalpost = totalpost + tempP}
        posts = NULL
        for(i in 1:k){
            initpost = (priors[i]*probs2[,i])/totalpost
            posts = cbind(posts, initpost)}
        priors = apply(posts, 2, mean)
        means = NULL
        if(ncol(x) > 1){
            for(i in 1:k){
                tempM = colSums(posts[,i]*x)/sum(posts[,i])
                means = rbind(means, tempM)}
        } else {
            for(i in 1:k){
                tempM = sum(posts[,i]*x)/sum(posts[,i])
                means = c(means, tempM)
            }
        }
        lambda = NULL
        for(i in 1:k){
            templambda = posts[,i]/(length(posts[,i])*priors[i])
            lambda = cbind(lambda, templambda)}
        cov.var = vector('list', length = k)
        if(ncol(x) > 1){
            for(i in 1:k){
                tempcov.var = 0
                for(l in 1:nrow(posts)){
                    initcov.var = lambda[l,i] * ((x[l,] - means[i,]) %*% t(x[l,] - means[i,]))
                    tempcov.var = tempcov.var + initcov.var}
                cov.var[[i]] = tempcov.var}
        } else {
            cov.var = vector('list', length = k)
            for(i in 1:k){
                tempcov.var = 0
                for(l in 1:nrow(posts)){
                    initcov.var = lambda[l,i] * ((x[l,] - means[i]) %*% t(x[l,] - means[i]))
                    tempcov.var =tempcov.var + initcov.var}
                cov.var[[i]] = sqrt(tempcov.var)
            }
        }

        if(show_cluster_probs == T){
        print(matrix(priors))}
    }
    clust = rep(NA, nrow(posts))
    for(i in 1:nrow(posts)){clust[i] = which(posts[i,] == max(posts[i,]), arr.ind = T)}
    rownames(means) = NULL
    return(list('Clusters' = clust, 'Means' = means, 'Covariance' = cov.var))
}


