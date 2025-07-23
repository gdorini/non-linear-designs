source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Função_principal/Função_principal_delineamentos.R", encoding = 'UTF-8')
source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Dist_prob/distribuições_de_probabilidade.R", encoding = 'UTF-8')

#Exemplos de uso da Função Principal

#Delineamentos Locais



#Contínuos

#Reacoes consec. (D-ótimo)
time.i <- Sys.time()
del.opt <- opt.design(design.type = c("local", "cont"), plot.var = c(TRUE, 100), expr = "p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p = c("p1", "p2"), fat = c("t"), a=0, b=10, N=2, Ni = 20, p.loc = c(1, 0.5), OTM = "D")
time.e <- Sys.time() - time.i
print(time.e)

#Reacoes consec. (As-ótimo)
time.i <- Sys.time()
del.opt <- opt.design(design.type = c("local", "cont"), plot.var = c(FALSE, 0), expr = "p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p = c("p1", "p2"), fat = c("t"), a=0, b=10, N=2, Ni = 20, p.loc = c(1, 0.5), OTM = "As", W=c(1, 0))
time.e <- Sys.time() - time.i
print(time.e)

#Modelo Compartimental (D-ótimo)
time.i <- Sys.time()
del.opt <- opt.design(design.type = c("local", "cont"), plot.var = c(TRUE, 100), expr = "p3*(Exp(-p2*t)-Exp(-p1*t))", p = c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=3, Ni = 50, p.loc = c(4.29, 0.0589, 21.80), OTM = "D")
time.e <- Sys.time() - time.i
print(time.e)

#Exemplo 1 parâmetro
time.i <- Sys.time()
del.opt <- opt.design(design.type = c("local", "cont"), plot.var = c(TRUE, 100), expr = "Exp(-p*t)", p = c("p"), fat = c("t"), a=0, b=5, N = 1, Ni = 50, p.loc = c(1), OTM = "D")
time.e <- Sys.time() - time.i
print(time.e)


#Exatos

#Michaelis-Menten (D-ótimo)
time.i <- Sys.time()
del.opt.e <- opt.design(design.type = c("local", "exact"), expr = "p1*t/(p2+t)", p=c("p1", "p2"), fat=c("t"), a=0, b=1.1, N=8, Ni=50, p.loc=c(216, 0.072), OTM="D")
time.e <- Sys.time() - time.i
print(time.e)

#Modelo Compartimental (As-ótimo)
time.i <- Sys.time()
del.opt <- opt.design(design.type = c("local", "exact"), expr = "p3*(Exp(-p2*t)-Exp(-p1*t))", p = c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=10, Ni = 100, p.loc = c(4.29, 0.0589, 21.80), OTM = "As", Ws=c(1, 10^4, 1/2))
time.e <- Sys.time() - time.i
print(time.e)




#Delineamentos pseudo-Bayesianos



#Contínuos

#Exemplo Michaelis-Menten (D-ótimo/cubature)

#Distribuicao normal 
medias <- c(216, 0.072)
dpadrao <- c(0.01, 0.05)
covariancia <- matrix(c(dpadrao[1]^2, 0, 0, dpadrao[2]^2), nrow = 2)
param <- normal.bayes(p = c("p1", "p2"), mu = medias, cov=covariancia, low = c(medias[1]-4*dpadrao[1], 0.01), up = c(medias[1]+4*dpadrao[1], 0.13))

time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "cont"), plot.var = c(TRUE, 100), expr="p1*t/(p2+t)", p=c("p1", "p2"), fat = c("t"), a=0, b=1.1, N=2, OTM="D", Ni=25, param=param, parallel = TRUE, metodo = "cubature")
time.e <- Sys.time() - time.i
print(time.e)

#Exemplo Michaelis-Menten com p1 linearizado (D-ótimo/cubature)

#Distribuicao normal p1 local 
medias <- 0.072
dpadrao <- 0.05
variancia <- dpadrao^2
param.p1loc <- normal.bayes(p = c("p2"), mu = medias, cov = variancia, low=0.01, up=0.13)

#Exemplo Michaelis-Menten (p1 loc)
time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "cont"), plot.var = c(TRUE, 100), expr="p1*t/(p2+t)", p=c("p1", "p2"), fat = c("t"), a=0, b=1.1, N=2, OTM="D", Ni=25, param=param.p1loc, p.loc = list(p1=216), p.bayes=c("p2"), parallel = TRUE, metodo = "MC", Np=100)
time.e <- Sys.time() - time.i
print(time.e)

#Exemplo - Modelo Compartimental (D-ótimo/cubature)

#Distribuicao uniforme
mu <- c(0.05884, 4.298, 21.80)
low <- c(mu[1]-0.01, mu[2]-1, mu[3]-0.001)
up <- c(mu[1]+0.01, mu[2]+1, mu[3]+0.001)
param.uniff <- uniff.bayes(low, up)

#Exemplo compartimental
time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "cont"), plot.var = c(FALSE, 100), expr="p3*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=3, OTM="D", Ni=30, param=param.uniff, parallel = TRUE, metodo = "cubature")
time.e <- Sys.time() - time.i
print(time.e)


#Exemplo - Modelo Compartimental (D-ótimo/Monte Carlo)

#Distribuicao uniforme
mu <- c(0.05884, 4.298, 21.80)
low <- c(mu[1]-0.01, mu[2]-1, mu[3]-0.001)
up <- c(mu[1]+0.01, mu[2]+1, mu[3]+0.001)
param.uniff <- uniff.bayes(low, up)

#Exemplo compartimental
time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "cont"), plot.var = c(FALSE, 100), expr="p3*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=3, OTM="D", Ni=30, param=param.uniff, parallel = TRUE, metodo = "MC", Np = 100)
time.e <- Sys.time() - time.i
print(time.e)


#Exemplo - Modelo Compartimental (D-ótimo/cubature) (p3 linearizado)

#Distribuicao uniforme
mu <- c(0.05884, 4.298)
low <- c(mu[1]-0.01, mu[2]-1)
up <- c(mu[1]+0.01, mu[2]+1)
param.uniff <- uniff.bayes(low, up)

#Exemplo compartimental
time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "cont"), plot.var = c(FALSE, 100), expr="p3*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=3, OTM="D", Ni=30, param=param.uniff, p.loc = list(p3=21.80), p.bayes = c("p1", "p2"), parallel = TRUE, metodo = "cubature")
time.e <- Sys.time() - time.i
print(time.e)


#Exemplo - Modelo Compartimental (D-ótimo/Monte Carlo) (p3 linearizado)

#Distribuicao uniforme
mu <- c(0.05884, 4.298)
low <- c(mu[1]-0.01, mu[2]-1)
up <- c(mu[1]+0.01, mu[2]+1)
param.uniff <- uniff.bayes(low, up)

#Exemplo compartimental
time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "cont"), plot.var = c(FALSE, 100), expr="p3*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=3, OTM="D", Ni=30, param=param.uniff, p.loc = list(p3=21.80), p.bayes = c("p1", "p2"), parallel = TRUE, metodo = "MC", Np=100)
time.e <- Sys.time() - time.i
print(time.e)



#Exatos 


# Exemplo Reacoes consecutivas (D-ótimo/cubature)

#Distribuicao normal 
medias.2 <- c(1, 0.5)
#dpadrao <- c(0.01, 0.05)
covariancia.2 <- matrix(c(1, 0.25, 0.25, 1), nrow = 2)
param.2 <- normal.bayes(p = c("p1", "p2"), mu = medias.2, cov=covariancia.2, low = c(0.8, 0.3), up = c(1.2, 0.7))

#Exemplo reacoes-consec
time.i <- Sys.time()
del.opt.bayes <- opt.design(design.type = c("bayes", "exact"), expr="p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p=c("p1", "p2"), fat = c("t"), a=0, b=10, N=5, OTM="D", Ni=25, param=param.2, parallel = TRUE, metodo = "cubature")
time.e <- Sys.time() - time.i
print(time.e)