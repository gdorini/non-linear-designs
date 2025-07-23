source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/função_principal_del_bayes_eff_fp(contínuo).R", encoding = 'UTF-8')

#OBSERVAÇÃO: quando escrever a expressão do modelo, sempre definir o x_máximo como xm e o x_mínimo como x0.

#explicação para as variáveis

#expr: expressão para o modelo do polinômio fracionário já reparametrizado
#p.lin: vetor com os parametros lineares do modelo
#alpha: vetor com o parâmetro não linear do modelo (alpha)
#fat: vetor com o fator/covariavel do modelo
#a: limite inferior do delineamento
#b: limite superior do delineamento
#N: dimensão do delineamento
#OTM: criterio de otimização 
#Ws: caso OTM = As, Ws é um vetor de pesos

#Exemplo 1

#modelo linear pf
medias <- 2.5
dpadrao <- 1.5
variancia <- dpadrao^2
param.p1loc <- normal.bayes(p = c("p2"), mu = medias, cov = variancia, low=medias-3*dpadrao, up=medias+3*dpadrao)

#distribuicao alpha
alpha <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
p.alpha <- c(0.15, 0.25, 0.25, 0.15, 0.10, 0.07, 0.03)
d.alpha <- alpha.fp(alpha, p.alpha)

#polinomio fracionario linear
del.fp <- del.bayes.cont.fp(expr = "b0 + p1*x^p2/(xm^p2-x0^p2)", p.lin = c("b0", "p1"), alpha = c("p2"), fat = c("x"), a = 0.1, b = 1, N = 4, OTM = "D", Ni = 15,
                            param = param.p1loc, param.alpha =  d.alpha, p.loc = NULL, p.bayes = c("p1"), parallel = TRUE, metodo = "cubature", Np = 5)


#Exemplo 2

#modelo linear pf
medias <- 2.5
dpadrao <- 1.5
variancia <- dpadrao^2
param.p1loc <- normal.bayes(p = c("p2"), mu = medias, cov = variancia, low=medias-3*dpadrao, up=medias+3*dpadrao)

#distribuicao alpha
alpha <- c(-2, 0, 0.5, 1)
p.alpha <- rep(1/4, 4)
d.alpha <- alpha.fp(alpha, p.alpha)

#polinomio fracionario linear 2
del.fp <- del.bayes.cont.fp(expr = "b0 + p1*x^p2/(xm^p2-x0^p2)", p.lin = c("b0", "p1"), alpha = c("p2"), fat = c("x"), a = 0.1, b = 1, N = 4, OTM = "As", Ws=c(0, 0.137931, 0.862069), Ni = 15,
                            param = param.p1loc, param.alpha =  d.alpha, p.loc = NULL, p.bayes = c("p1"), parallel = TRUE, metodo = "MC", Np = 5)

#Exemplo 3

#modelo linear pf
medias <- 2.5
dpadrao <- 1
variancia <- dpadrao^2
param.p1loc <- normal.bayes(p = c("p2"), mu = medias, cov = variancia, low=medias-3*dpadrao, up=medias+3*dpadrao)

#distribuicao alpha
alpha <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
p.alpha <- c(0.15, 0.25, 0.25, 0.15, 0.10, 0.07, 0.03)
d.alpha <- alpha.fp(alpha, p.alpha)

Ws=c(0, 0.137931, 0.862069)
#polinomio fracionario linear 2
del.fp <- del.bayes.cont.fp(expr = "b0 + p1*x^p2/(xm^p2-x0^p2)", p.lin = c("b0", "p1"), alpha = c("p2"), fat = c("x"), a = 0.1, b = 1, N = 4, OTM = "As", Ws=c(0, 0.137931, 0.862069), Ni = 30,
                            param = param.p1loc, param.alpha =  d.alpha, p.loc = NULL, p.bayes = c("p1"), parallel = TRUE, metodo = "cubature", Np = 5)


#Exemplo 4

#modelo quadratico pf
medias <- c(1, -2.5)
dpadrao <- c(0.2, 0.6)
covariancia <- matrix(c(dpadrao[1]^2, 0, 0, dpadrao[2]^2), nrow = 2)
param <- normal.bayes(p = c("p2", "p4"), mu = medias, cov=covariancia, low = c(medias[1]-3*dpadrao[1], medias[2]-3*dpadrao[2]), up = c(medias[1]+3*dpadrao[1],medias[2]+3*dpadrao[2]))

#distribuicao alpha 2
alpha <- c( -1, -0.5, 0, 2)
p.alpha <- rep(0.25, 4)
d.alpha <- alpha.fp(alpha, p.alpha)

#polinomio fracionario quadratico
del.fp <- del.bayes.cont.fp(expr = "p1+p2*x^p3/(xm^p3-x0^p3)+4*p4*((x^(2*p3))-(xm^p3+x0^p3)*x^p3)/(xm^p3-x0^p3)^2", p.lin = c("p1", "p2", "p4"), alpha = c("p3"), fat = c("x"), a = 0.1, b = 1, N = 4, OTM = "A", Ni = 15,
                            param = param, param.alpha =  d.alpha, p.loc = NULL, p.bayes = c("p2", "p4"), parallel = TRUE, metodo = "cubature", Np = 5)


