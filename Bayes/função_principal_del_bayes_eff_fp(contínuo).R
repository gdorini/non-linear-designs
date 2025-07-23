source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/funções_del_bayes_cont_fp.R", encoding = 'UTF-8')
source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Dist_prob/distribuições_de_probabilidade.R", encoding = 'UTF-8')
library(Ryacas)
library(cubature)
library(parallel)
library(mvtnorm)
library(MASS)
library(truncnorm)

#funcao principal polinomio frac
del.bayes.cont.fp <- function(expr, p.lin, alpha, fat, a, b, N, OTM, Ws=NULL, Ni, param, param.alpha, p.loc = NULL, p.bayes = NULL, parallel = FALSE, metodo, Np = 10^3){
  comeco <- Sys.time()
  
  #derivadas do modelo
  p <- c(p.lin, alpha) #vetor de parametros com alpha em ultimo
  l.p <- length(p)
  l.f <- length(fat)
  l.p.bayes <- length(p.bayes)
  mod <- ysym(expr) #modelo em yacas
  d.mod <- deriv(mod, p) #derivadas 
  
  if(OTM=="As"){
    if(is.null(Ws) == FALSE){
      Ws <- ysym(as_y(diag(Ws)))}
    else{
      Ws <- ysym(as_y(diag(1, l.p)))
    }
  }
  
  #distribuicao de probabilidade para hcubature
  if(metodo == "cubature"){
    #Dist. prob.
    d_prob <- param[[1]] #distribuicao de probabilidade
    low <- param[[2]] #limite inferior da dist.
    up <- param[[3]] #limite superior da dist.
  }
  
  #distribuicao para Monte Carlo
  if(metodo == "MC"){
    if(length(param) > 3){
      d_prob <- param[[1]] #distribuicao de probabilidade
      low <- param[[2]] #limite inferior da dist.
      up <- param[[3]] #limite superior da dist.
      medias <- param[[4]] #media dos parametros
      cov <- param[[5]] #matriz de covariancia dos parametros
      dominio <- amostra_normal_MC(N = Np, low = low, up = up, mu = medias, sigma = cov) #pontos para a integracao
      }
    else{
      d_prob <- param[[1]] #distribuicao de probabilidade
      low <- param[[2]] #limite inferior da dist.
      up <- param[[3]] #limite superior da dist.
    }
  }
  
  #distribuicao alpha
  alph <- param.alpha[[1]] #valores de alpha
  l.alph <- length(alph)
  p.alph <- param.alpha[[2]] #probabilidade de cada alpha
  lgc <- param.alpha[[3]] #operador logico que indica se há alpha = 0
  
  #alpha vira a nas derivadas do modelo
  aux.alpha <- ajusta.expr.fp(p.lin, p, d.mod, a, b) #substituindo alpha por a
  d.mod <- aux.alpha[[1]] 
  d.mod.2 <- aux.alpha[[2]] 
  
  #ajusta expressoes
  d.mod.ajust <- ajusta.expr(p.bayes, l.p.bayes, l.f, fat, d.mod) #derivadas do modelo ajustadas
  d.mod.ajust.2 <- ajusta.expr(p.bayes, l.p.bayes, l.f, fat, d.mod.2) #para a=0
  
  #construcao da matriz de informacao simbolica
  IM <- IM_SIMB(d.mod.ajust, l.p) #prototipo da matriz de informacao simbolica
  
  IM.2 <- IM_SIMB(d.mod.ajust.2, l.p) #se a=0
  
  I <- N
  #Delineamentos iniciais (população inicial do espaço)
  m = 2*I 
  N_TENT = 10*m #numero inicial de delineamentos
  AUX <- Del_0(I, m, N_TENT, a, b)
  V <- AUX[[1]] #matriz dos delineamentos iniciais
  Z <- AUX[[2]] #matriz dos pesos
  
  #soma da matriz de inf. para a dimensao do delineamento
  IM_del_aux <- IM.sum(I, IM, l.p) #inf. matrix simbolica completa
  IM_del_aux.2 <- IM.sum(I, IM.2, l.p) #a=0
  
  integrando <- criterio.int(metodo, OTM, Ws, IM_del_aux, d_prob, l.p) #funcao que sera integrada
  integrando.2 <- criterio.int(metodo, OTM, Ws, IM_del_aux.2, d_prob, l.p) #para a=0
  
  
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
  #paralelismo
  if(parallel ==  TRUE){
    num_cores <- detectCores() - 1 #numero de nucleos utilizados no paralelismo
    my_cluster <- makeCluster(num_cores, type = "PSOCK") #deinicao do cluster
    clusterEvalQ(my_cluster, library(cubature)) #rodando o pacote cubature no cluster
    
    if(metodo == "cubature"){
      #exportando variaveis para o cluster
      variables_to_export <- c("integrando", "integrando.2", "low", "up", "alph", "l.alph", "p.alph", "lgc", "integra_crit_v.fp", "integra_crit_v")
      clusterExport(my_cluster, varlist = variables_to_export, envir = environment())
      
    }
    if(metodo == "MC"){
      #exportando variaveis para o cluster
      variables_to_export <- c("integrando", "integrando.2", "low", "up", "alph", "l.alph", "p.alph", "lgc", "integral_MC_IS.fp", "integral_MC_IS", "dominio")
      clusterExport(my_cluster, varlist = variables_to_export, envir = environment())
    }
    
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
      #exportando variaveis para o cluster
      variables_to_export.2 <- c("V", "I", "CRIT")
      clusterExport(my_cluster, varlist = variables_to_export.2, envir = environment())
      
      if(lgc == 0){
        if(metodo == "cubature"){
          worker_function <- function(k){
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I), k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            
            for(k.a in 1:(l.alph-1)){
              AUX.CRIT <- AUX.CRIT + integra_crit_v.fp(integrando, low, up, X, W, alpha = alph[k.a])[[1]]*p.alph[k.a]
            }
            CRIT[k] <- AUX.CRIT + integra_crit_v(integrando.2, low, up, X, W)[[1]]*p.alph[l.alph]
            options(warn = 0)
            return(CRIT[k])
          }
          iterations <- 1:N_TENT
          aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
          
          CRIT <- unlist(aux_parallel_list)
        }
        
        if(metodo == "MC"){
            worker_function <- function(k){
              print("t1")
              #Obtenção da matriz de informação
              X <- V[1:I, k]
              W <- V[(I+1):(2*I), k]
              
              options(warn = -1)
              AUX.CRIT <- 0
              for(k.a in 1:(l.alph-1)){
                AUX.CRIT <- AUX.CRIT + integral_MC_IS.fp(integrando, dominio, X, W, alpha = alph[k.a])*p.alph[k.a]
              }
              CRIT[k] <- AUX.CRIT + integral_MC_IS(integrando.2, dominio, X, W)*p.alph[l.alph]
              options(warn = 0)
              return(CRIT[k])
            }
            iterations <- 1:N_TENT
            aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
            
            CRIT <- unlist(aux_parallel_list)
          }
      }
      else{
        if(metodo == "cubature"){
          worker_function <- function(k){
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I), k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            
            for(k.a in 1:(l.alph)){
              AUX.CRIT <- AUX.CRIT + integra_crit_v.fp(integrando, low, up, X, W, alpha = alph[k.a])[[1]]*p.alph[k.a]
            }
            CRIT[k] <- AUX.CRIT
            options(warn = 0)
            return(CRIT[k])
          }
          iterations <- 1:N_TENT
          aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
          
          CRIT <- unlist(aux_parallel_list)
        }
        
        if(metodo == "MC"){
          worker_function <- function(k){
            print("t1")
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I), k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            for(k.a in 1:(l.alph)){
              AUX.CRIT <- AUX.CRIT + integral_MC_IS.fp(integrando, dominio, X, W, alpha = alph[k.a])*p.alph[k.a]
            }
            CRIT[k] <- AUX.CRIT 
            options(warn = 0)
            return(CRIT[k])
          }
          iterations <- 1:N_TENT
          aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
          
          CRIT <- unlist(aux_parallel_list)
        }
      }
      #elimina os crit NAN
      T.F <- is.na(CRIT)
      CRIT[T.F] <- max(na.omit(CRIT)) 
      
      #obtenção do ótimo
      if(cont==1){aux = CRIT #aux guarda os p-best
      l <- 1 #indica a coluna do delineamento otimo
      s <- cont #indica a iteracao na qual o otimo foi obtido
      OPT <- CRIT[1]} 
      
      LISTA <- P_BEST_AJ(cont, aux=aux, crit=CRIT, OPT, V=V, P=P, z=Z, I, l, s, N_TENT=N_TENT)
      P <- LISTA[[1]] #p-best
      l <- LISTA[[2]] #indica a coluna do delineamento otimo
      s <- LISTA[[3]] #indica a iteracao na qual o otimo foi obtido
      OPT <- LISTA[[4]] #indica o valor do criterio otimo
      aux <- LISTA[[5]] #aux guarda os p-best
      
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
    stopCluster(my_cluster)}
  
  #sem paralelismo
  else{
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
      if(lgc == 0){
        if(metodo == "cubature"){
          for(k in 1:(N_TENT)){
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I),k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            for(k.a in 1:(l.alph-1)){
              AUX.CRIT <- AUX.CRIT + integra_crit_v.fp(integrando, low, up, X, W, alpha = alph[k.a])[[1]]*p.alph[k.a]
            }
            CRIT[k] <- AUX.CRIT + integra_crit_v(integrando.2, low, up, X, W)[[1]]*p.alph[l.alph]
            options(warn = 0)
          }
        }
        
        if(metodo == "MC"){
          for(k in 1:(N_TENT)){
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I),k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            for(k.a in 1:(l.alph-1)){
              AUX.CRIT <- AUX.CRIT + integral_MC_IS.fp(integrando, dominio, X, W, alpha = alph[k.a])*p.alph[k.a]
            }
            CRIT[k] <- AUX.CRIT + integral_MC_IS(integrando.2, dominio, X, W)*p.alph[l.alph]
            options(warn = 0)
          }
        }
      }
      else{
        if(metodo == "cubature"){
          for(k in 1:(N_TENT)){
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I),k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            for(k.a in 1:l.alph){
              AUX.CRIT <- AUX.CRIT + integra_crit_v.fp(integrando, low, up, X, W, alpha = alph.2[k.a])[[1]]*p.alph.2[k.a]
            }
            CRIT[k] <- AUX.CRIT 
            options(warn = 0)
          }
        }
        
        if(metodo == "MC"){
          for(k in 1:(N_TENT)){
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I),k]
            
            options(warn = -1)
            AUX.CRIT <- 0
            for(k.a in 1:l.alph){
              AUX.CRIT <- AUX.CRIT + integral_MC_IS.fp(integrando, dominio, X, W, alpha = alph[k.a])*p.alph[k.a]
            }
            CRIT[k] <- AUX.CRIT 
            options(warn = 0)
          }
        }
      }
      #elimina os crit NAN
      T.F <- is.na(CRIT)
      CRIT[T.F] <- max(na.omit(CRIT)) 
      
      #obtenção do ótimo
      if(cont==1){aux = CRIT #aux guarda os p-best
      l <- 1 #indica a coluna do delineamento otimo
      s <- cont #indica a iteracao na qual o otimo foi obtido
      OPT <- CRIT[1]} 
      
      LISTA <- P_BEST_AJ(cont, aux=aux, crit=CRIT, OPT, V=V, P=P, z=Z, I, l, s, N_TENT=N_TENT)
      P <- LISTA[[1]] #p-best
      l <- LISTA[[2]] #indica a coluna do delineamento otimo
      s <- LISTA[[3]] #indica a iteracao na qual o otimo foi obtido
      OPT <- LISTA[[4]] #indica o valor do criterio otimo
      aux <- LISTA[[5]] #aux guarda os p-best
      
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
  }
  AUX4 <- g_best[(I+1):(2*I)]
  design <- cbind(g_best[1:I], g_best[(I+1):(2*I)]^2/(sum(AUX4^2)))
  colnames(design) <- c("X", "W")
  cat(paste0("\n", "for ", I, " points", "\n"))
  print(design)
  cat(paste0("Criterio = ", OPT, "\n"))
  fim <- Sys.time()-comeco
  print(fim)
}