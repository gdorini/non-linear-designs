#Função Principal
opt.design <- function(design.type, plot.var=c(FALSE, 0), expr, p, fat, a, b, N, OTM, Ws=NULL, Ni, param=NULL, p.loc = NULL, p.bayes = NULL, parallel = FALSE, metodo, Np = 10^3){
  if(design.type[1] == "bayes"){
    if(design.type[2] == "exact"){
      source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/função_principal_del_bayes_eff(exato).R", encoding = 'UTF-8')
      return(del.bayes.exact(expr = expr, p=p, fat=fat, a=a, b=b, N=N, OTM=OTM, Ni=Ni, param = param, p.loc = p.loc, p.bayes=p.bayes, parallel = parallel, metodo = metodo, Np = Np))
    }
    if(design.type[2] == "cont"){
      source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/função_principal_del_bayes_eff(contínuo).R", encoding = 'UTF-8')
      
      otimo <- del.bayes.cont(expr = expr, p=p, fat=fat, a=a, b=b, N=N, OTM=OTM, Ws=Ws, Ni=Ni, param = param, p.loc = p.loc, p.bayes = p.bayes, parallel = parallel, metodo = metodo, Np = Np)
      if(plot.var[1] == TRUE && OTM == "D"){
        std.var(del.opt = otimo, n = plot.var[2])
      }
    }
  }
  if(design.type[1] == "local"){
    if(design.type[2] == "exact"){
      source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Local/Função_principal_del_loc_simb_exato.R", encoding = 'UTF-8')
      del.loc.e(expr=expr, p=p, fat=fat, a=a, b=b, N=N, Ni=Ni, p_loc=p.loc, OTM=OTM, Ws=Ws)
    }
    if(design.type[2] == "cont"){
      source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Local/Função_principal_del_local_simb.R", encoding = 'UTF-8')
      otimo <- del.loc(expr=expr, p=p, fat=fat, a=a, b=b, NP=N, Ni=Ni, p_loc=p.loc, OTM=OTM, Ws=Ws)
      if(plot.var[1] == TRUE && OTM == "D"){
        var.pad(otimo)
      }
    }
  }
}