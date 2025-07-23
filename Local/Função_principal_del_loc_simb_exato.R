source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Local/Funções_del_local_simb_exato.R", encoding = 'UTF-8')
library(Ryacas)
library(cubature)

del.loc.e <- function(expr, p, fat, a, b, N, Ni, p_loc, OTM, Ws=NULL){
  
  l.p <- length(p) 
  l.f <- length(fat)
  mod <- ysym(expr) #modelo em yacas
  d.mod <- deriv(mod, p) #derivadas do modelo
  
  #matriz de pesos para o delineamento As-otimo
  if(OTM == "As"){
    if(is.null(Ws) == FALSE){
      Ws <- diag(Ws)
    }
    else{
      Ws <- diag(1, l.p)
    }
  }
  
  #altera p para p[i] e os fatores para x na expressão
  d.mod_f <- ajusta.expr.local(p, p_loc, l.p, l.f, fat, d.mod)
  
  #construcao da matriz de informacao simbolica (parametros locais substituidos)
  IM_S <- IM_SIMB(d.mod_f, l.p)
  IM <- IM_S[[1]]
  IM <- y_eval(IM, w=1)
  IM_R <- as_r(IM)
  
  #Delineamentos iniciais (população inicial do espaço)
  N_TENT = 20*N
  V <- Del_0.e(N, N_TENT, a, b)
  
  #definição de algumas variáveis
  
  CRIT <- numeric()
  aux <- numeric()
  DV <- matrix(0, N, N_TENT)
  P <- matrix(0, N, N_TENT)
  P <- V
  X <- numeric()
  W <- numeric()
  der1 <- numeric()
  der2 <- numeric()
  g_best <- numeric()
  AUX3 <- numeric()
  AUX4 <- numeric()
  
  #início das iterações
  cont<-1
  
  while (cont<=Ni) {
    TEMPO.1.MM <- Sys.time()
    
    if(cont > 3){
    AUX.mut <- mutacao.otm.e(CRIT, V, cont, Ni, N, P, a, b, aux)
    CRIT <- AUX.mut[[1]]
    V <- AUX.mut[[2]]
    P <- AUX.mut[[3]]
    aux <- AUX.mut[[4]]
    N_TENT <- AUX.mut[[5]]
    l <- AUX.mut[[6]]
    }
    
    if(OTM == "D"){
      for(k in 1:(N_TENT)){
        
        #Obtenção da matriz de informação
        X <- V[1:N, k]
        l.X <- length(X)
        IM_del <- matrix(0, l.p, l.p)
        for(i in 1:l.X){
          IM_del <- IM_del+eval(IM_R, list(x=X[i]))
        }
        IM_del <- IM_del/N
        options(warn = -1)
        IM_det <- det(IM_del)
        CRIT[k] <- -log(IM_det)
        options(warn = 0)
      }
    }
    
    if(OTM == "A"){
      for(k in 1:(N_TENT)){
        
        #Obtenção da matriz de informação
        X <- V[1:N, k]
        l.X <- length(X)
        IM_del <- matrix(0, l.p, l.p)
        for(i in 1:l.X){
          IM_del <- IM_del+eval(IM_R, list(x=X[i]))
        }
        IM_del.inverse <- tryCatch(
          {
            # Tentar resolver o sistema linear
            solve(IM_del) 
          },
          error = function(e) {
            if (grepl("system is computationally singular", e$message)) {
              # Matriz singular detectada
              print("Matriz singular detectada. Substituída por matriz 10^25.")
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
        #IM_del.inverse <- solve(IM_del)
        sum.A <- 0
        for(q in 1:l.p){
          sum.A <- sum.A+IM_del.inverse[q,q]
        }
        CRIT[k] <- sum.A
      }
    }
    
    if(OTM == "As"){
      #Ws <- diag(Ws)
      for(k in 1:(N_TENT)){
        
        #Obtenção da matriz de informação
        X <- V[1:N, k]
        l.X <- length(X)
        IM_del <- matrix(0, l.p, l.p)
        for(i in 1:l.X){
          IM_del <- IM_del+eval(IM_R, list(x=X[i]))
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
      }
    }
    
    #elimina os crit NAN
    T.F <- is.na(CRIT)
    CRIT[T.F] <- max(na.omit(CRIT))
    
    #obtenção do ótimo
    if(cont==1){aux = CRIT
    l <- 1
    s <- 1
    OPT <- CRIT[1]} #aux guarda os p-best
    
    aux.s <- s
    LISTA <- P_BEST_AJ.e(cont, aux=aux, crit=CRIT, OPT, V=V, P=P, I=N, l=l, s=s, N_TENT=N_TENT)
    P <- LISTA[[1]]
    l <- LISTA[[2]]
    s <- LISTA[[3]]
    OPT <- LISTA[[4]]
    aux <- LISTA[[5]]
    
    if(cont==1){
      for(e in 1:N){
        g_best[e] <- V[e, l]
      }}
    if(aux.s != s){
      for(e in 1:N){
        g_best[e] <- V[e, l]
      }
    }
    
    design.aux <- g_best[1:N]
    G_BEST <- matrix(g_best, N, N_TENT)
    
    #Cálculo das constantes do algoritmo
    C <- cte.3.e(cont, Ni, V, G_BEST, a, b, I=N)
    c1 <- C[[1]]
    c2 <- C[[2]]
    
    #Cálculo do deslocamento de cada indivíduo
    AUX2 <- DESLOCAMENTO.2.e(N_TENT, P, V, G_BEST, c1, c2, I=N)
    DX <- AUX2
    
    for(h in 1:N_TENT){
      for(i in 1:N){
        V[i,h] <- V[i,h]+DX[i,h]
      }}
    #Ajuste dos indivíduos deslocados para fora da região de delineamento
    V <- AJUSTE.e(I=N, N_TENT, V, a, b)
    
    design <- as.matrix(sort(g_best[1:N]))
    cat(paste0("\n", "Iteração = ", cont, "\n"))
    colnames(design) <- c("X")
    print(design)
    cat(paste0("Criterio = ", OPT, "\n"))
    FIM.MM <- Sys.time() - TEMPO.1.MM
    print(FIM.MM)
    
    cont <- cont+1
  }
}