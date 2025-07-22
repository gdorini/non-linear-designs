source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/funções_del_bayes_cont.R", encoding = 'UTF-8')
source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Dist_prob/distribuições_de_probabilidade.R", encoding = 'UTF-8')
library(Ryacas)
library(cubature)
library(mvtnorm)
library(MASS)
library(truncnorm)
library(parallel)

del.bayes.cont <- function(expr, p, fat, a, b, N, OTM, Ws=NULL, Ni, param, p.loc = NULL, p.bayes = NULL, parallel = FALSE, metodo, Np = 10^3){
  l.p <- length(p) #numero de parametros do modelo
  l.f <- length(fat) #numero de fatores do modelo
  mod <- ysym(expr) #modelo em yacas
  d.mod <- deriv(mod, p) #derivadas do modelo
  
  if(OTM == "As"){
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
  
  I <- N
    #Delineamentos iniciais (população inicial do espaço)
    m = 2*I 
    N_TENT = 10*m #numero inicial de delineamentos
    AUX <- Del_0(I, m, N_TENT, a, b)
    V <- AUX[[1]] #matriz dos delineamentos iniciais
    Z <- AUX[[2]] #matriz dos pesos
    
    #soma da matriz de inf. para a dimensao do delineamento
    IM_del_aux <- IM.sum(I, IM, l.p) #inf. matrix simbolica completa
    
    integrando <- criterio.int(metodo, OTM, Ws, IM_del_aux, d_prob, l.p)
    
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
    if(parallel ==  TRUE){
      num_cores <- detectCores() - 1
      my_cluster <- makeCluster(num_cores, type = "PSOCK")
      clusterEvalQ(my_cluster, library(cubature))
      
      #exportando variaveis para o cluster
      variables_to_export <- c("integrando", "low", "up", "integra_crit_v", "integral_MC_IS")
      clusterExport(my_cluster, varlist = variables_to_export, envir = environment())
      
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
      
      if(metodo ==  "cubature"){
        worker_function <- function(k){
          #Obtenção da matriz de informação
          X <- V[1:I, k]
          W <- V[(I+1):(2*I), k]
          
          options(warn = -1)
          CRIT[k] <- integra_crit_v(integrando, low, up, X, W)[[1]]
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
          X <- V[1:I, k]
          W <- V[(I+1):(2*I), k]
          
          options(warn = -1)
          CRIT[k] <- integral_MC_IS(integrando, dominio, x=X, w=W)
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
        
        if(metodo == "cubature"){
        for(k in 1:(N_TENT)){

          #Obtenção da matriz de informação
          X <- V[1:I, k]
          W <- V[(I+1):(2*I),k]

          options(warn = -1)
          CRIT[k] <- integra_crit_v(integrando, low, up, X, W)[[1]] #integral
          options(warn = 0)
        }
        }
        
        if(metodo == "MC"){
          for(k in 1:(N_TENT)){
            
            #Obtenção da matriz de informação
            X <- V[1:I, k]
            W <- V[(I+1):(2*I),k]
            
            options(warn = -1)
            CRIT[k] <- integral_MC_IS(integrando, dominio, x=X, w=W) #integral
            options(warn = 0)
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

  return(list(design, IM, l.p, d.mod_f, d_prob, a, b, low, up, fat))
}

std.var <- function(del.opt, n){
  ti_d <- Sys.time()
  
  #Definacao de variaveis
  design <- del.opt[[1]]
  IM <- del.opt[[2]]
  l.p <- del.opt[[3]]
  d.mod_f <- del.opt[[4]]
  d_prob <- del.opt[[5]]
  a <- del.opt[[6]]
  b <- del.opt[[7]]
  low <- del.opt[[8]]
  up <- del.opt[[9]]
  fat <- del.opt[[10]]
  if(is.numeric(d_prob)==TRUE){
    d_prob <- ysym(d_prob)
  }
  
  #Construcao da matriz de informacao do delineamento otimo
  IM_otm <- matrix(0, l.p, l.p)
  for(i in 1:nrow(design)){
    IM_otm <- IM_otm + y_eval(IM, x=design[i, 1], w=design[i, 2])
  }
  IM_otm_inv <- solve(IM_otm)
  #matriz das derivadas do modelo
  aux <- IM_SIMB.2(d.mod_f, l.p)
  y_f <- aux[[1]]
  y_ft <- aux[[2]]
  
  #variancia padrao dependente de theta
  #d_theta <- as_r(yac_symbol(y_ft)%*%yac_symbol(as_y(IM_otm_inv))%*%yac_symbol(y_f))
  parc1 <- yac_symbol(y_ft)%*%IM_otm_inv
  d_theta <- parc1%*%yac_symbol(y_f)
  
  #integrando
  integrando.1 <- as.character(as_r(d_theta*d_prob))
  
  #integral
  func.int.1 <- function(p, x){
    return(eval(parse(text=integrando.1)))
  }
  
  x_d <- seq(a, b, by=(b-a)/n)
  d <- numeric()
  
  num_cores <- detectCores() - 1
  my_cluster <- makeCluster(num_cores, type = "PSOCK")
  clusterEvalQ(my_cluster, library(cubature))
  
  #exportando variaveis para o cluster
  d.aux <- numeric()
  variables_to_export <- c("func.int.1", "low", "up", "x_d", "d.aux")
  clusterExport(my_cluster, varlist = variables_to_export, envir = environment())
  
  worker_function <- function(i){
    #Obtenção da matriz de informação
    d.aux[i] <- hcubature(func.int.1, low, up, x = x_d[i], tol=1e-4, maxEval = 1000)[[1]]
    return(d.aux[i])
  }
  iterations <- 1:length(x_d)
  aux_parallel_list <- parLapply(my_cluster, iterations, worker_function)
  
  d <- unlist(aux_parallel_list)
  stopCluster(my_cluster)
  
  plot(x = x_d, y = d, xlab = fat, ylab= "std. variance", type="n")
  lines(x=x_d, y=d)
  abline(a=l.p,b=0,col=2)
  tf_d <- Sys.time()-ti_d
  print(tf_d)
}
