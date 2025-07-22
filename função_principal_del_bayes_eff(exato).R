source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/funções_del_bayes_exato.R", encoding = 'UTF-8')
source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Dist_prob/distribuições_de_probabilidade.R", encoding = 'UTF-8')
library(Ryacas)
library(cubature)
library(doFuture)
library(mvtnorm)
library(MASS)
library(truncnorm)
library(parallel)

del.bayes.exact <- function(expr, p, fat, a, b, N, OTM, Ni, param, p.loc = NULL, p.bayes = NULL, parallel = FALSE, metodo, Np = 10^3){
  
  l.p <- length(p) 
  l.f <- length(fat)
  mod <- ysym(expr) #modelo em yacas
  d.mod <- deriv(mod, p) #derivadas do modelo
  
  #distribuicao de probabilidade para hcubature
  if(metodo == "cubature"){
    #Dist. prob.
    d_prob <- param[[1]] #distribuicao de probabilidade
    low <- param[[2]] #limite inferior da dist.
    up <- param[[3]] #limite superior da dist.
  }
  
  if(metodo == "MC"){
    if(length(param) > 3){
      d_prob <- param[[1]] #distribuicao de probabilidade
      low <- param[[2]] #limite inferior da dist.
      up <- param[[3]] #limite superior da dist.
      medias <- param[[4]] #media dos parametros
      cov <- param[[5]] #matriz de covariancia dos parametros
      dominio <- amostra_normal_MC(N = Np, low = low, up = up, mu = medias, sigma = cov)
    }
    else{
      d_prob <- param[[1]] #distribuicao de probabilidade
      low <- param[[2]] #limite inferior da dist.
      up <- param[[3]] #limite superior da dist.
      dominio <- amostra_uniff_MC(Np, low, up)
    }
  }
  
  #Atribuicao dos valores aos parametros locais
  if(is.null(p.loc) == FALSE){
    d.mod <- do.call(y_eval, c(list(d.mod), p.loc))
    l.p.loc <- l.p-length(p.loc)
    
    #altera p para p[i] e os fatores para x na expressão
    d.mod_f <- ajusta.expr(p.bayes, l.p.loc, l.f, fat, d.mod) #derivadas do modelo ajustadas
  }
  else{
    #altera p para p[i] e os fatores para x na expressão
    d.mod_f <- ajusta.expr(p, l.p, l.f, fat, d.mod) #derivadas do modelo ajustadas
  }
  
  #construcao da matriz de informacao simbolica
  IM <- IM_SIMB(d.mod_f, l.p) #prototipo da matriz de informacao simbolica
  IM <- y_eval(IM, w=1)
  
  #Delineamentos iniciais (população inicial do espaço)
  N_TENT = 20*N
  V <- Del_0.e(N, N_TENT, a, b)
  
  #soma da matriz de inf. para a dimensao do delineamento
  IM_del_aux <- IM.sum.e(N, IM, l.p) #inf. matrix simbolica completa
  
  if(metodo == "cubature"){
  #para delineamentos D-otimos
  if(OTM == "D"){
    IM_det <- det(IM_del_aux) #determinante da matriz de inf.
    
    integrando <- as.character(as_r(-log(IM_det)*d_prob)) #integrando
  }
  return(integrando)
  #para delineamentos A-otimos
  if(OTM == "A"){
    IM_del_aux.inv <- solve(IM_del_aux) #matriz de inf. inversa
    sum.A <- 0
    for(q in 1:l.p){
      sum.A <- sum.A+IM_del_aux.inv[q,q]
    }
    integrando <- as.character(as_r(sum.A*d_prob)) #integrando
  }
  }
  
  if(metodo == "MC"){
    #para delineamentos D-otimos
    if(OTM == "D"){
      IM_det <- det(IM_del_aux) #determinante da matriz de inf.
      
      options(warn = -1)
      integrando <- as.character(as_r(-log(IM_det))) #integrando
      options(warn = 0)
    }
    
    #para delineamentos A-otimos
    if(OTM == "A"){
      IM_del_aux.inv <- solve(IM_del_aux) #matriz de inf. inversa
      sum.A <- 0
      for(q in 1:l.p){
        sum.A <- sum.A+IM_del_aux.inv[q,q]
      }
      integrando <- as.character(as_r(sum.A)) #integrando
    }
  }
  
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
  
  if(parallel ==  TRUE){
    num_cores <- detectCores() - 1
    my_cluster <- makeCluster(num_cores, type = "PSOCK")
    clusterEvalQ(my_cluster, library(cubature))
    
    #exportando variaveis para o cluster
    variables_to_export <- c("integrando", "low", "up", "integra_crit_v.e", "integral_MC_IS.e")
    clusterExport(my_cluster, varlist = variables_to_export, envir = environment())
    
  while (cont<=Ni) {
    TEMPO.1.MM <- Sys.time()
    
    if(cont >= 3){
      AUX.mut <- mutacao.otm.e(CRIT, V, cont, Ni, N, P, a, b, aux)
      CRIT <- AUX.mut[[1]]
      V <- AUX.mut[[2]]
      P <- AUX.mut[[3]]
      aux <- AUX.mut[[4]]
      N_TENT <- AUX.mut[[5]]
      l <- AUX.mut[[6]]
    }
    
    if(metodo == "cubature"){
      worker_function <- function(k){
      #Obtenção da matriz de informação
      X <- V[1:N, k]
      
      # CRIT[k] <- integra_crit_2(integrando, low, up, X, W)
      options(warn = -1)
      CRIT[k] <- integra_crit_v.e(integrando, low, up, X)[[1]]
      options(warn = 0)
      return(CRIT[k])
      }
      iterations <- 1:N_TENT
      aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
      
      CRIT <- unlist(aux_parallel_list)
    }
    
    if(metodo == "MC"){
      worker_function <- function(k){
        #Obtenção da matriz de informação
        X <- V[1:N, k]
        
        # CRIT[k] <- integra_crit_2(integrando, low, up, X, W)
        options(warn = -1)
        CRIT[k] <- integral_MC_IS.e(integrando, dominio, X)
        options(warn = 0)
        return(CRIT[k])
      }
      iterations <- 1:N_TENT
      aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
      
      CRIT <- unlist(aux_parallel_list)
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
    stopCluster(my_cluster)}
  else{
    while (cont<=Ni) {
      TEMPO.1.MM <- Sys.time()
      
      if(cont >= 3){
        AUX.mut <- mutacao.otm.e(CRIT, V, cont, Ni, N, P, a, b, aux)
        CRIT <- AUX.mut[[1]]
        V <- AUX.mut[[2]]
        P <- AUX.mut[[3]]
        aux <- AUX.mut[[4]]
        N_TENT <- AUX.mut[[5]]
        l <- AUX.mut[[6]]
      }
      
      if(metodo == "cubature"){
      for(k in 1:(N_TENT)){
        
        #Obtenção da matriz de informação
        X <- V[1:N, k]
        
        # CRIT[k] <- integra_crit_2(integrando, low, up, X, W)
        options(warn = -1)
        CRIT[k] <- integra_crit_v.e(integrando, low, up, X)[[1]]
        options(warn = 0)
      }
      }
      
      if(metodo == "MC"){
        for(k in 1:(N_TENT)){
          
          #Obtenção da matriz de informação
          X <- V[1:N, k]
          
          # CRIT[k] <- integra_crit_2(integrando, low, up, X, W)
          options(warn = -1)
          CRIT[k] <- integral_MC_IS.e(integrando, dominio, X)
          options(warn = 0)
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
}