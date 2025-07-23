#Dist. de probabilidades

#NORMAL

normal.bayes <- function(p, mu, cov, low, up){
  c_i <- solve(cov)
  options(warn = -1)
  d_prob <- d_prob_normal(p, cov, c_i, mu)

  int.d_prob <- as.character(as_r(d_prob))
  integrando.d_prob <- function(p){
    eval(parse(text=int.d_prob))
  }
  integral.d_prob <- hcubature(integrando.d_prob, low, up, tol=1e-4)
  d_prob <- d_prob/integral.d_prob[[1]]
  options(warn = 0)
  return(list(d_prob, low, up, mu, cov))
}

#UNIFORME

uniff.bayes <- function(low, up){
  d_prob <- prod(1/(up-low))
  return(list(d_prob, low, up))
}

#POLINOMIO FRACIONARIO

alpha.fp <- function(alpha, p.alpha){
  if(0 %in% alpha == TRUE){
    alpha.2 <- alpha[!alpha %in% 0]
    pos.0 <- which(alpha == 0)
    p.alpha.2 <- p.alpha[-pos.0]
    alpha <- c(alpha.2, alpha[pos.0])
    p.alpha <- c(p.alpha.2, p.alpha[pos.0])
    aux <- 0
  }
  else{aux <- 1}
  return(list(alpha, p.alpha, aux))
}

#distribuicao de probabilidade
d_prob_normal <- function(p, covariancia, c_i, medias){
  y_p <- ysym(p)
  l.p <- length(p)
  
  if(l.p != 1){
    y_p_mu <- ysym(p)-medias
    
    p_mu <- numeric()
    for(i in 1:l.p){
      p_mu[i] <- as.matrix(yac_str(y_p_mu[i]))
    }
    
    y_p_mu <- yac_symbol(as_y(as.matrix(p_mu)))
    y_p_mu_t <- yac_symbol(as_y(t(p_mu)))
    
    arg <- 0.5*y_p_mu_t%*%yac_symbol(c_i)%*%y_p_mu
    
    cte <- (2*pi)^(-l.p/2)*(det(covariancia))^(-1/2)
    norm_prob <- cte*exp(-arg)
    
    aux <- rep("p[i]", l.p)
    aux  <- ysym(aux)
    list.p <- list()
    for(j in 1:l.p){
      aux[j] <- y_eval(aux[j], i=j)
      list.p[[j]] <- aux[j]
      names(list.p)[j] <- p[j]
    }
    
    norm_prob_2 <- do.call(y_eval, c(list(norm_prob), list.p))
  }
  
  else{
    y_p_mu <- ysym(p)-medias
    
    cte <- 1/(sqrt(2*pi)*covariancia)
    arg <- y_p_mu^2/(2*covariancia)
    norm_prob <- cte*exp(-arg)
    
    aux <- "p[1]"
    aux  <- ysym(aux)
    list.p <- list()
    list.p[[1]] <- aux
    names(list.p)[1] <- p
    norm_prob_2 <- do.call(y_eval, c(list(norm_prob), list.p))
  }
  return(norm_prob_2)
}
