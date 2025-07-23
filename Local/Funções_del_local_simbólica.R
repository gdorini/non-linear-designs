#FUNCOES DEL LOCAL CONTINUO

Del_0 <- function(I, m, N_TENT, a, b){
  V <- matrix(0, nr=m, nc=N_TENT) 
  Z <- matrix(0, I, N_TENT)
  z <- rep(0, I)
  aux2 <- numeric()
  v <- rep(0, m)
  for(i in 1:(N_TENT)){
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

ajusta.expr.local <- function(p, p_loc, l.p, l.f, fat, d.mod){
  if(l.p != 1){
    list.p <- list()
    list.f <- list()
    options(warn = -1)
    for(j in 1:l.p){
      list.p[[j]] <- p_loc[j]
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
    list.p <- list()
    list.f <- list()
    list.p[[1]] <- p_loc[1]
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
  return(list(IM, y_f, y_ft))
}

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
  # corte <- 10*2*I+5*I
  # if(ncol(V)>corte){
  #   V <- V[,1:corte]
  #   Z <- Z[,1:corte]
  #   N_TENT <- ncol(V)
  # }
  return(list(CRIT, V, Z, P, aux, N_TENT, l))
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

criterio <- function(OTM, N_TENT, V, l.p, IM_R, Ws, I){
  #D-otimo
  if(OTM == "D"){
    for(k in 1:(N_TENT)){
      #Obtenção da matriz de informação
      X <- V[1:I, k]
      W <- V[(I+1):(2*I),k]
      l.X <- length(X)
      IM_del <- matrix(0, l.p, l.p)
      for(i in 1:l.X){
        IM_del <- IM_del+eval(IM_R, list(x=X[i], w=W[i]))
      }
      options(warn = -1)
      IM_det <- det(IM_del)
      CRIT[k] <- -log(IM_det)
      options(warn = 0)
    }}
  #A-otimo
  if(OTM == "A"){
    for(k in 1:(N_TENT)){
      
      #Obtenção da matriz de informação
      X <- V[1:I, k]
      W <- V[(I+1):(2*I),k]
      l.X <- length(X)
      #IM_del <- as_y(matrix(0, l.X, l.X))
      IM_del <- matrix(0, l.p, l.p)
      for(i in 1:l.X){
        IM_del <- IM_del+eval(IM_R, list(x=X[i], w=W[i]))
      }
      IM_inv <- tryCatch(
        {
          # Tentar resolver o sistema linear
          solve(IM_del) 
        },
        error = function(e) {
          if (grepl("system is computationally singular", e$message) || grepl("system is exactly singular", e$message)) {
            # Matriz singular detectada
            #print("Matriz singular detectada. Substituída por matriz 10^25.")
            IM_sub <- diag(nrow(IM_del))
            for(n in 1:nrow(IM_del)){
              IM_sub[n, n] <- 10^10*IM_sub[n, n]} 
            return(IM_sub)
          } else {
            # Relançar outros erros
            stop(e) 
          }
        }
      )
      CRIT[k] <- sum(diag(IM_inv))
    }}
  
  #As-otimo
  if(OTM == "As"){
    if(is.null(Ws)==F){
      Ws <- diag(Ws)
    }
    else{
      Ws <- diag(1, l.p)
    }
    for(k in 1:(N_TENT)){
      
      #Obtenção da matriz de informação
      X <- V[1:I, k]
      W <- V[(I+1):(2*I),k]
      l.X <- length(X)
      #IM_del <- as_y(matrix(0, l.X, l.X))
      IM_del <- matrix(0, l.p, l.p)
      for(i in 1:l.X){
        IM_del <- IM_del+eval(IM_R, list(x=X[i], w=W[i]))
      }
      IM_inv <- tryCatch(
        {
          # Tentar resolver o sistema linear
          solve(IM_del) 
        },
        error = function(e) {
          if (grepl("system is computationally singular", e$message) || grepl("system is exactly singular", e$message)) {
            # Matriz singular detectada
            #print("Matriz singular detectada. Substituída por matriz 10^25.")
            IM_sub <- diag(nrow(IM_del))
            for(n in 1:nrow(IM_del)){
              IM_sub[n, n] <- 10^10*IM_sub[n, n]} 
            return(IM_sub)
          } else {
            # Relançar outros erros
            stop(e) 
          }
        }
      )
      CRIT[k] <- sum(diag(Ws%*%IM_inv))
    }}
  return(list(CRIT, l.X))
}