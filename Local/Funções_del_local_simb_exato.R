#FUNCOES DEL. LOCAL EXATO

Del_0.e <- function(N,N_TENT, a, b){
  V <- matrix(0, N, N_TENT) 
  v <- rep(0, N)
  for(i in 1:(N_TENT)){
    R <- sample(0:100, size = N, replace=T)
    
    #sorteio dos níveis
    for(j in 1:N){
      v[j] <- a+0.01*R[j]*(b-a)
      V[j, i] <- v[j]
    }
  }
  return(V)
}

P_BEST_AJ.e <- function(cont, aux, crit, OPT, V, P, I, l, s, N_TENT){
  for(i in 1:N_TENT){
    if(crit[i] < OPT){OPT <- crit[i]
    l <- i
    s <- cont}
    if(crit[i] < aux[i]){aux[i] <- crit[i]
    for(j in 1:I){
      P[j,i] <- V[j, i]
    }}
  }
  return(list(P, l, s, OPT, aux))
}


DESLOCAMENTO <- function(N_TENT, P, V, G_BEST, c1, c2, c3, c4, Z){
  DX <- matrix(0, I, N_TENT)
  DZ <- matrix(0, I, N_TENT)
  for(j in 1:(N_TENT)){
    R2 <- sample(0:10, size = 2, replace=T)
    r1 <- 0.1*R2[1]
    r2 <- 0.1*R2[2]
    for(i in 1:m){
      if(i <= I){
        DX[i, j] <- c1*r1*(P[i, j]-V[i, j]) + c2*r2*(G_BEST[i,j]-V[i,j])
      }
      else{
        DZ[i-I,j] <- c3*r1*(P[i, j]-Z[i-I,j]) + c4*r2*(G_BEST[i,j]-Z[i-I,j])
      }
    }}
  return(list(DX, DZ))
}

cte.3.e <- function(cont, Ni, V, G_BEST, a, b, I){
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
        }
        else{
          c1[j,i] <- 0.1
          c2[j,i] <- 0.4
        }
      }
      else{
        if(cont < Ni/3){
          c1[j,i] <- 4*aux[j,i]/5
          c2[j,i] <- aux[j,i]/5
        }
        else{
          c1[j,i] <- aux[j,i]/5
          c2[j,i] <- 4*aux[j,i]/5
        }
      }}}
  C <- list(c1, c2)
  return(C)
}

DESLOCAMENTO.2.e <- function(N_TENT, P, V, G_BEST, c1, c2, I){
  DX <- matrix(0, I, N_TENT)
  for(j in 1:(N_TENT)){
    R2 <- sample(0:10, size = 2, replace=T)
    r1 <- 0.1*R2[1]
    r2 <- 0.1*R2[2]
    for(i in 1:I){
      #DX[i, j] <- c1*(P[i, j]-V[i, j]) + c2*(G_BEST[i,j]-V[i,j])
      DX[i, j] <- c1[i,j]*(P[i, j]-V[i, j]) + c2[i,j]*(G_BEST[i,j]-V[i,j])
    }}
  return(DX)
}

AJUSTE.e <- function(I, N_TENT, V, a, b){
  for(j in 1:(N_TENT)){
    for(i in 1:I){
      if(V[i,j] < a){
        V[i,j] <- a
      }
      if(V[i,j] > b){
        V[i,j] <- b
      }
    }
  }
  return(V)
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

mutacao.otm.e <- function(CRIT, V, cont, Ni, N, P, a, b, aux){
  I=N
  CRIT.aux <- CRIT
  CRIT.aux.2 <- round(CRIT, 6)
  c.ord <- sort(CRIT)
  erro <- c.ord[1]+0.05*abs(c.ord[1])
  CRIT <- unique(round(c.ord[c.ord<erro], 6))
  l<-1
  pos.vec <- numeric()
  for(h in 1:length(CRIT)){
    pos.vec[h] <- which(CRIT.aux.2 == CRIT[h])[1]
  }
  V <- as.matrix(V[,pos.vec, drop=FALSE])
  N_TENT <- ncol(V)
  
  if(cont >= 3 && Ni-0.3*Ni >= cont){
    
    new.c1 <- matrix(V[,1], nrow=I, ncol=I)
    new.c2 <- matrix(V[,1], nrow=I, ncol=I)
    new.c3 <- matrix(V[,1], nrow=I, ncol=I)
    new.c4 <- matrix(V[,1], nrow=I, ncol=I)
    new.c5 <- matrix(V[,1], nrow=I, ncol=I)
    new.c6 <- matrix(V[,1], nrow=I, ncol=I)
    new.c7 <- matrix(V[,1], nrow=I, ncol=I)
    new.c8 <- matrix(V[,1], nrow=I, ncol=I)
    new.c9 <- matrix(V[,1], nrow=I, ncol=I)
    new.c10 <- matrix(V[,1], nrow=I, ncol=I)
    ncol.aj <- 10*I
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
    }
    
    new.c <- cbind(new.c1, new.c2, new.c3, new.c4, new.c5, new.c6, new.c7, new.c8, new.c9, new.c10)
    if(ncol(V)>1){
      ncv.aux <- ncol(V)
      V <- cbind(V[,1], new.c, V[,2:ncv.aux])}
    else{
      V <- cbind(V[,1], new.c)
    }
    N_TENT <- ncol(V)
    V <- AJUSTE.e(I=N, N_TENT, V, a, b)
    
    aux.aux <- rep(sort(aux)[1], ncol.aj)
    aux <- c(aux, aux.aux)
    p.aux <- matrix(P[,1], nrow=nrow(P), ncol=ncol.aj)
    P <- cbind(P, p.aux)
  }
  corte <- 10*I+5*I
  if(ncol(V)>corte){
    V <- V[,1:corte]
    N_TENT <- ncol(V)
  }
  return(list(CRIT, V, P, aux, N_TENT, l))
}