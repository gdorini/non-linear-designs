#FUNÇÕES BAYES CONTINUO

#Determina os delineamentos inciais
Del_0 <- function(I, m, N_TENT, a, b){
  V <- matrix(0, nr=m, nc=N_TENT) 
  Z <- matrix(0, I, N_TENT)
  z <- rep(0, I)
  aux2 <- numeric()
  v <- rep(0, m)
  for(i in 1:(N_TENT)){
    # R <- sample(0:100, size = I, replace=F)
    R <- sample(1:100, size = I, replace=F)
    
    #Sorteio dos pesos (todos sobre uma esfera I-dim de raio 10)
    aux1 <- 0
    for(h in 1:I){
      if(aux1 <= 100){
        aux2[h] <- sample(1:(100-aux1), size=1, replace=T)
        z[h] <- sqrt(aux2[h])
        aux1 <- aux1 + aux2[h]
      }}
    
    #sorteio dos níveis
    for(j in 1:m){
      if(j <= I){
        v[j] <- a+0.01*R[j]*(b-a)
      }
      else{
        v[j] <- z[(j-I)]^2/sum(z^2)
      }
    }
    V[,i] <- v
    Z[,i] <- z
  }
  return(list(V, Z))
}

#formata a expressao simbolica do modelo
ajusta.expr <- function(p, l.p, l.f, fat, d.mod){
  if(l.p != 1){
    auxiliar <- rep("p[i]", l.p)
    auxiliar  <- ysym(auxiliar)
    list.p <- list()
    list.f <- list()
    options(warn = -1)
    for(j in 1:l.p){
      auxiliar[j] <- y_eval(auxiliar[j], i=j)
      list.p[[j]] <- auxiliar[j]
      names(list.p)[j] <- p[j]
    }
    options(warn = 0)
    d.mod_2 <- do.call(y_eval, c(list(d.mod), list.p))
    for(i in 1:l.f){
      list.f[[i]] <- ysym("x")
      names(list.f)[i] <- fat[i]
    }
    d.mod_f <- do.call(y_eval, c(list(d.mod_2), list.f))
  }
  else{
    auxiliar <- "p[1]"
    auxiliar  <- ysym(auxiliar)
    list.p <- list()
    list.f <- list()
    list.p[[1]] <- auxiliar
    names(list.p)[1] <- p
    d.mod_2 <- do.call(y_eval, c(list(d.mod), list.p))
    for(i in 1:l.f){
      list.f[[i]] <- ysym("x")
      names(list.f)[i] <- fat[i]
    }
    d.mod_f <- do.call(y_eval, c(list(d.mod_2), list.f))
  }
  return(d.mod_f)
}

#controi a matriz de informacao

IM_SIMB <- function(d.mod_f, l.p){
  f <- numeric()
  
  if(l.p != 1){
    for(i in 1:l.p){
      f[i] <- yac_str(d.mod_f[i])
    }}
  else{
    f <- yac_str(d.mod_f)
  }
  
  y_f <- as_y(as.matrix(f))
  y_ft <- as_y(t(f))
  
  y_F <- yac_symbol(y_f)%*%yac_symbol(y_ft)
  w <- ysym("w")
  IM <- w*y_F
  return(IM)
}

IM_SIMB.2 <- function(d.mod_f, l.p){
  f <- numeric()
  
  if(l.p != 1){
  for(i in 1:l.p){
    f[i] <- yac_str(d.mod_f[i])
  }}
  else{
    f <- yac_str(d.mod_f)
  }
  
  y_f <- as_y(as.matrix(f))
  y_ft <- as_y(t(f))
  
  y_F <- yac_symbol(y_f)%*%yac_symbol(y_ft)
  w <- ysym("w")
  IM <- w*y_F
  return(list(y_f, y_ft))
}

#soma a matriz de informacao para a dimensao do delineamento
IM.sum <- function(I, IM, l.p){
  x_aux <- ysym(rep("x[i]", I))
  w_aux <- ysym(rep("w[i]", I))
  options(warn = -1)
  for(j in 1:I){
    x_aux[j] <- y_eval(x_aux[j], i=j)
    w_aux[j] <- y_eval(w_aux[j], i=j)
  }
  options(warn = 0)
  IM_del_aux <- as_y(matrix(0, l.p, l.p))
  
  for(i in 1:I){
    IM_del_aux <- IM_del_aux + y_eval(IM, x=x_aux[i], w=w_aux[i])
  }
  return(IM_del_aux)
}

#integrando 
criterio.int <- function(metodo, OTM, Ws=NULL, IM_del_aux, d_prob, l.p){
  if(metodo == "cubature"){
    #para delineamentos D-otimos
    if(OTM == "D"){
      IM_det <- det(IM_del_aux) #determinante da matriz de inf.
      
      integrando <- as.character(as_r(-log(IM_det)*d_prob)) #integrando
    }
    
    #para delineamentos A-otimos
    if(OTM == "A"){
      IM_del_aux.inv <- solve(IM_del_aux)
      sum.A <- sum(diag(IM_del_aux.inv))
      integrando <- as.character(as_r(sum.A*d_prob)) #integrando
      print(integrando)
    }
    
    #para delineamentos A-otimos ponderados
    if(OTM == "As"){
      IM_del_aux.inv <- solve(IM_del_aux) #matriz de inf. inversa
      IM_del_aux.inv <- Ws%*%IM_del_aux.inv
      sum.A <- sum(diag(IM_del_aux.inv))
      integrando <- as.character(as_r(sum.A*d_prob)) #integrando
    }
  }
  
  if(metodo == "MC"){
    #para delineamentos D-otimos
    if(OTM == "D"){
      IM_det <- det(IM_del_aux) #determinante da matriz de inf
      
      integrando <- as.character(as_r(-log(IM_det))) #integrando
    }
    
    #para delineamentos A-otimos
    if(OTM == "A"){
      IM_del_aux.inv <- solve(IM_del_aux)
      sum.A <- sum(diag(IM_del_aux.inv))
      integrando <- as.character(as_r(sum.A)) #integrando
    }
    
    #para delineamentos A-otimos aprox.
    if(OTM == "As"){
      IM_del_aux.inv <- solve(IM_del_aux) #matriz de inf. inversa
      IM_del_aux.inv <- Ws%*%IM_del_aux.inv
      sum.A <- sum(diag(IM_del_aux.inv))
      integrando <- as.character(as_r(sum.A)) #integrando
    }
  }
  return(integrando)
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

#pontos para integracao por monte carlo
amostra_normal_MC <- function(N=10^3, low, up, mu, sigma){
  dim <- length(mu)
  
  if(dim > 1){
  #Amostra normal multivariada
  amostra <- mvrnorm(N, mu, sigma)
  
  #selecionando os pontos dentro do dominio de integracao
  c <- list()
  nr.aux <- numeric()
  for(i in 1:dim){
    c.aux <- amostra[,i]
    c[[i]] <- c.aux[c.aux>=low[i] & c.aux<=up[i]]
    nr.aux[i] <- length(c[[i]])
  }
  nr <- min(nr.aux)
  N <- nr
  dominio <- matrix(0, nrow = nr, ncol = 1)
  for (i in 1:dim) {
    c[[i]] <- c[[i]][1:nr]
    dominio <- cbind(dominio, c[[i]])
  }
  dominio <- dominio[,-1]}
  else{
    dp <- sqrt(sigma)
    dominio <- rtruncnorm(N, low, up, mu, dp)
  }
  
  return(dominio)
}

amostra_uniff_MC <- function(N, low, up){
  aux <- list()
  aux.2 <- rep(0, N)
  for(i in 1:length(low)){
    aux[[i]] <- runif(N, min=low[i], max=up[i])
    aux.2 <- cbind(aux.2, aux[[i]])
  }
  aux.2 <- aux.2[,-1]
  return(aux.2)
}

#realiza a integraacao por quadratura(hcubature)
integra_crit_v <- function(integrando, low, up, X, W){
  
  func_int_v <- function(par, x, w){
    matrix(apply(par, 2, function(p){return(eval(parse(text=integrando)))}), ncol = ncol(par))
  }
  
  h_int_v <- hcubature(func_int_v, low, up, x = X, w = W, tol=1e-4, maxEval = 1000, vectorInterface = TRUE)
  return(h_int_v)
}

#realiza a integral por metodo de monte carlo
integral_MC_IS <- function(integrando, dominio, x, w){
  f_int.mc <- function(p, x, w){
    return(eval(parse(text=integrando)))
  }
  
  #numero de pontos
  N <- nrow(dominio)
  
  #avaliar a funcao em cada ponto
  f_dominio <- numeric()
  if(is.matrix(dominio)==TRUE){
  N <- nrow(dominio)
  for(j in 1:N){
    f_dominio[j] <- f_int.mc(dominio[j,], x, w)
  }}
  else{
    N <- length(dominio)
    for(j in 1:N){
      f_dominio[j] <- f_int.mc(dominio[j], x, w)
    }
  }
  #soma
  integral.mc <- sum(f_dominio)/N
  return(integral.mc)
}

#parametros do PSO

P_BEST_AJ <- function(cont, aux, crit, OPT, V, P, z, I, l, s, N_TENT){
  for(i in 1:N_TENT){
    if(crit[i] < OPT){OPT <- crit[i]
    l <- i
    s <- cont}
    if(crit[i] < aux[i]){aux[i] <- crit[i]
    for(j in 1:I){
      P[j,i] <- V[j, i]
      P[j+I,i] <- z[j,i]
    }}
  }
  return(list(P, l, s, OPT, aux))
}

cte.3 <- function(cont, Ni, p3, p4, V, G_BEST, a, b, I){
  ncol.V <- ncol(V) 
  aux <- matrix(0, I, ncol.V)
  c1 <- matrix(0, I, ncol.V)
  c2 <- matrix(0, I, ncol.V)
  for(j in 1:I){
    for(i in 1:ncol.V){
      #aux[j,i] <- abs((b-a)/(25*(G_BEST[j,i]-V[j,i])))
      aux[j,i] <- abs((b-a)/((G_BEST[j,i]-V[j,i])))
      if(G_BEST[j,i]-V[j,i]<(b-a)/50){
        if(cont < Ni/3){
          c1[j,i] <- 0.4
          c2[j,i] <- 0.1
          c3 <- p3
          c4 <- p4
        }
        else{
          c1[j,i] <- 0.1
          c2[j,i] <- 0.4
          c3 <- p3
          c4 <- p4
        }
      }
      else{
        if(cont < Ni/3){
          c1[j,i] <- 4*aux[j,i]/5
          c2[j,i] <- aux[j,i]/5
          c3 <- p3
          c4 <- p4
        }
        else{
          c1[j,i] <- aux[j,i]/5
          c2[j,i] <- 4*aux[j,i]/5
          c3 <- p4
          c4 <- p3
        }
      }}}
  C <- list(c1, c2, c3, c4)
  return(C)
}

DESLOCAMENTO.2 <- function(N_TENT, P, V, G_BEST, c1, c2, c3, c4, Z, I){
  m = 2*I
  DX <- matrix(0, I, N_TENT)
  DZ <- matrix(0, I, N_TENT)
  for(j in 1:(N_TENT)){
    R2 <- sample(0:10, size = 2, replace=T)
    r1 <- 0.1*R2[1]
    r2 <- 0.1*R2[2]
    for(i in 1:m){
      if(i <= I){
        #DX[i, j] <- c1*(P[i, j]-V[i, j]) + c2*(G_BEST[i,j]-V[i,j])
        DX[i, j] <- c1[i,j]*(P[i, j]-V[i, j]) + c2[i,j]*(G_BEST[i,j]-V[i,j])
      }
      else{
        DZ[i-I,j] <- c3*r1*(P[i, j]-Z[i-I,j]) + c4*r2*(G_BEST[i,j]-Z[i-I,j])
      }
    }}
  return(list(DX, DZ))
}

AJUSTE <- function(I, m, N_TENT, V, a, b){
  for(j in 1:(N_TENT)){
    for(i in 1:m){
      if(i < I+1){
        if(V[i,j] < a){
          V[i,j] <- a
        }
        if(V[i,j] > b){
          V[i,j] <- b
        }
      }
    }
  }
  return(V)
}

#realiza a mutacao do delineamento otimo na iteracao em questao
mutacao.otm <- function(CRIT, V, Z, cont, Ni, I, P, m, a, b, aux){
  CRIT.aux <- CRIT
  CRIT.aux.2 <- round(CRIT, 6)
  c.ord <- sort(CRIT)
  erro <- c.ord[1]+0.05*abs(c.ord[1])
  CRIT <- unique(round(c.ord[c.ord<erro], 6))
  #print(CRIT[1])
  l<-1
  pos.vec <- numeric()
  for(h in 1:length(CRIT)){
    pos.vec[h] <- which(CRIT.aux.2 == CRIT[h])[1]
  }
  V <- as.matrix(V[,pos.vec])
  Z <- as.matrix(Z[,pos.vec, drop=FALSE])
  N_TENT <- ncol(V)
  
  if(cont >= 3 && Ni-0.3*Ni >= cont){
    
    new.c1 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c2 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c3 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c4 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c5 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c6 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c7 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c8 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c9 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.c10 <- matrix(V[,1], nrow=2*I, ncol=2*I)
    new.z1 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z2 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z3 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z4 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z5 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z6 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z7 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z8 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z9 <- matrix(Z[,1], nrow=I, ncol=I)
    new.z10 <- matrix(Z[,1], nrow=I, ncol=I)
    ncol.aj <- 10*2*I
    rdm <- sample(5:10, 1)/10
    for(i in 1:I){
      new.c1[i,i] <- new.c1[i,i]+rdm*0.1*new.c1[i,i]
      new.c2[i,i] <- new.c2[i,i]-rdm*0.1*new.c2[i,i]
      new.c3[i,i] <- new.c3[i,i]+rdm*0.05*new.c3[i,i]
      new.c4[i,i] <- new.c4[i,i]-rdm*0.05*new.c4[i,i]
      new.c5[i,i] <- new.c5[i,i]+rdm*0.01*new.c5[i,i]
      new.c6[i,i] <- new.c6[i,i]-rdm*0.01*new.c6[i,i]
      new.c7[i,i] <- new.c7[i,i]+rdm*0.5*new.c7[i,i]
      new.c8[i,i] <- new.c8[i,i]-rdm*0.5*new.c8[i,i]
      new.c9[i,i] <- new.c9[i,i]+rdm*0.001*new.c9[i,i]
      new.c10[i,i] <- new.c10[i,i]-rdm*0.001*new.c10[i,i]
      new.z1[i,i] <- new.z1[i,i]+rdm*0.1*new.z1[i,i]
      new.z2[i,i] <- new.z2[i,i]-rdm*0.1*new.z2[i,i]
      new.z3[i,i] <- new.z3[i,i]+rdm*0.05*new.z3[i,i]
      new.z4[i,i] <- new.z4[i,i]-rdm*0.05*new.z4[i,i]
      new.z5[i,i] <- new.z5[i,i]+rdm*0.01*new.z5[i,i]
      new.z6[i,i] <- new.z6[i,i]-rdm*0.01*new.z6[i,i]
      new.z7[i,i] <- new.z7[i,i]+rdm*0.5*new.z7[i,i]
      new.z8[i,i] <- new.z8[i,i]-rdm*0.5*new.z8[i,i]
      new.z9[i,i] <- new.z9[i,i]+rdm*0.001*new.z9[i,i]
      new.z10[i,i] <- new.z10[i,i]-rdm*0.001*new.z10[i,i]
    }
    
    for(j in 1:I){
      for(i in 1:I){
        new.c1[i+I,j+I] <- new.z1[i,j]^2/sum(new.z1[,j]^2)
        new.c2[i+I,j+I] <- new.z2[i,j]^2/sum(new.z2[,j]^2)
        new.c3[i+I,j+I] <- new.z3[i,j]^2/sum(new.z3[,j]^2)
        new.c4[i+I,j+I] <- new.z4[i,j]^2/sum(new.z4[,j]^2)
        new.c5[i+I,j+I] <- new.z5[i,j]^2/sum(new.z5[,j]^2)
        new.c6[i+I,j+I] <- new.z6[i,j]^2/sum(new.z6[,j]^2)
        new.c7[i+I,j+I] <- new.z7[i,j]^2/sum(new.z7[,j]^2)
        new.c8[i+I,j+I] <- new.z8[i,j]^2/sum(new.z8[,j]^2)
        new.c9[i+I,j+I] <- new.z9[i,j]^2/sum(new.z9[,j]^2)
        new.c10[i+I,j+I] <- new.z10[i,j]^2/sum(new.z10[,j]^2)
      }}
    
    new.c <- cbind(new.c1, new.c2, new.c3, new.c4, new.c5, new.c6, new.c7, new.c8, new.c9, new.c10)
    z.aux <- matrix(0, I, I)
    for(s in 1:I){
      z.aux[,s] <- Z[,1]
    }
    new.z<- cbind(z.aux, new.z1, z.aux, new.z2, z.aux, new.z3, z.aux, new.z4, z.aux, new.z5, z.aux, new.z6, z.aux, new.z7, z.aux, new.z8, z.aux, new.z9, z.aux, new.z10)
    if(ncol(V)>1){
      ncv.aux <- ncol(V)
      V <- cbind(V[,1], new.c, V[,2:ncv.aux])
      ncz.aux <- ncol(Z)
      Z <- cbind(Z[,1],new.z, Z[,2:ncz.aux, drop=FALSE])}
    else{
      V <- cbind(V[,1], new.c)
      Z <- cbind(Z[,1], new.z)
    }
    #aux.z <- matrix(Z[,1], nrow=I, ncol=ncol.aj)
    N_TENT <- ncol(V)
    V <- AJUSTE(I, m, N_TENT, V, a, b)
    
    aux.aux <- rep(sort(aux)[1], ncol.aj)
    aux <- c(aux, aux.aux)
    p.aux <- matrix(P[,1], nrow=nrow(P), ncol=ncol.aj)
    P <- cbind(P, p.aux)
  }
  corte <- 10*2*I+5*I
  if(ncol(V)>corte){
    V <- V[,1:corte]
    Z <- Z[,1:corte]
    N_TENT <- ncol(V)
  }
  return(list(CRIT, V, Z, P, aux, N_TENT, l))
}
