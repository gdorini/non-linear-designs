source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Local/Funções_del_local_simbólica.R", encoding = 'UTF-8')
library(Ryacas)
library(cubature)

del.loc <- function(expr, p, fat, a, b, NP, Ni, p_loc, OTM, Ws=NULL){
  
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
  IM_R <- as_r(IM)
  
  I <- NP
    #Delineamentos iniciais (população inicial do espaço)
    m = 2*I
    N_TENT = 200*m
    AUX <- Del_0(I, m, N_TENT, a, b)
    V <- AUX[[1]]
    Z <- AUX[[2]]
    
    #definição de algumas variáveis
    
    CRIT <- numeric()
    aux <- numeric()
    DV <- matrix(0, m, N_TENT)
    P <- matrix(0, m, N_TENT)
    P <- rbind(V,Z)
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
      if(cont >= 3){
        AUX.mut <- mutacao.otm(CRIT, V, Z, cont, Ni, I, P, m, a, b, aux) #funcao geradora de mutacoes ao redor do otimo
        CRIT <- AUX.mut[[1]]
        V <- AUX.mut[[2]]
        Z <- AUX.mut[[3]]
        P <- AUX.mut[[4]]
        aux <- AUX.mut[[5]]
        N_TENT <- AUX.mut[[6]]
        l <- AUX.mut[[7]]
      }
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
          #print(IM_inv)
        }}

      #As-otimo
      if(OTM == "As"){
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
          #print(Ws%*%IM_inv)
        }}
      
      # AUX.CRIT <- criterio(OTM, N_TENT, V, l.p, IM_R, Ws, I)
      # CRIT <- AUX.CRIT[[1]]
      # l.X <- AUX.CRIT[[2]]
      
      #elimina os crit NAN
      T.F <- is.na(CRIT)
      CRIT[T.F] <- max(na.omit(CRIT))
      
      #obtenção do ótimo
      if(cont==1){aux = CRIT
      l <- 1
      s <- cont
      OPT <- CRIT[1]} #aux guarda os p-best
      
      LISTA <- P_BEST_AJ(cont, aux=aux, crit=CRIT, OPT, V=V, P=P, z=Z, I, l, s, N_TENT=N_TENT)
      P <- LISTA[[1]]
      l <- LISTA[[2]]
      s <- LISTA[[3]]
      OPT <- LISTA[[4]]
      aux <- LISTA[[5]]
      
      if(cont==1){
        for(e in 1:m){
          if(e <= I){
            g_best[e] <- V[e, l]}
          else{
            g_best[e] <- Z[e-I, l]
          }
        }}
      if(s != 1){
        for(e in 1:m){
          if(e <= I){
            g_best[e] <- V[e, l]}
          else{
            g_best[e] <- Z[e-I, l]
          }}
      }
      
      G_BEST <- matrix(g_best, m, N_TENT)
      
      #Cálculo das constantes do algoritmo
      C <- cte.3(cont, Ni, 8, 2, V, G_BEST, a, b, I)
      c1 <- C[[1]]
      c2 <- C[[2]]
      c3 <- C[[3]]
      c4 <- C[[4]]
      
      #Cálculo do deslocamento de cada indivíduo
      AUX2 <- DESLOCAMENTO.2(N_TENT, P, V, G_BEST, c1, c2, c3, c4, Z, I)
      DX <- AUX2[[1]]
      DZ <- AUX2[[2]]
      Z <- Z + DZ
      for(q in 1:N_TENT){
        if(sum(Z[,q]^2) == 0){
          Z[,q] <- sample(1:100, size=I, replace=T)
        }
      }
      
      for(h in 1:N_TENT){
        for(i in 1:I){
          V[i,h] <- V[i,h]+DX[i,h]
          V[i+I,h] <- Z[i,h]^2/sum(Z[,h]^2)
        }}
      
      #Ajuste dos indivíduos deslocados para fora da região de delineamento
      V <- AJUSTE(I, m, N_TENT, V, a, b)
      
      
      FIM.MM <- Sys.time() - TEMPO.1.MM
      
      AUX4 <- g_best[(I+1):(2*I)]
      design <- cbind(g_best[1:I], g_best[(I+1):(2*I)]^2/(sum(AUX4^2)))
      cat(paste0("\n", "Iteração = ", cont, "\n"))
      colnames(design) <- c("X", "W")
      print(design)
      cat(paste0("Criterio = ", OPT, "\n"))
      print(FIM.MM)
      cont <- cont+1
    }
    AUX4 <- g_best[(I+1):(2*I)]
    design <- cbind(g_best[1:I], g_best[(I+1):(2*I)]^2/(sum(AUX4^2)))
    colnames(design) <- c("X", "W")
    cat(paste0("\n", "for ", I, " points", "\n"))
    print(design)
    cat(paste0("Criterio = ", OPT, "\n"))
    
    IM_otm <- matrix(0, l.p, l.p)
    y_f <- IM_S[[2]]
    y_ft <- IM_S[[3]]
    for(i in 1:l.X){
      IM_otm <- IM_otm + eval(IM_R, list(x=design[i, 1], w=design[i, 2]))
    }
    IM_otm_inv <- solve(IM_otm)
    d <- as_r(yac_symbol(y_ft)%*%yac_symbol(as_y(IM_otm_inv))%*%yac_symbol(y_f))
    
    del.opt <- list(design, d, c(a,b), l.p, fat)
    return(del.opt)
}

var.pad <- function(del.opt){
  a <- del.opt[[3]][1]
  b <- del.opt[[3]][2]
  d <- del.opt[[2]]
  l.p <- del.opt[[4]]
  fat <- del.opt[[5]]
  
  ab <- seq(a, b, by=(b-a)/20000)
  ord <- eval(d, list(x = ab))
  ord.round <- round(ord, digits=4)
  max.x <- which(ord.round==l.p)
  plot(x=ab, y=ord, xlab = fat, ylab= "std. variance", type ="n", xlim = c(a, b))
  lines(x=ab, y=ord)
  abline(a=l.p,b=0,col=2)
  }