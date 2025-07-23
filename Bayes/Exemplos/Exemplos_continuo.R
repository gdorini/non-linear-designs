source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Bayes/função_principal_del_bayes_eff(contínuo).R", encoding = 'UTF-8')

#Exemplo 1 - Michaelis Menten

#Distribuicao normal 1
medias <- c(216, 0.072)
dpadrao <- c(0.01, 0.05)
covariancia <- matrix(c(dpadrao[1]^2, 0, 0, dpadrao[2]^2), nrow = 2)
param <- normal.bayes(p = c("p1", "p2"), mu = medias, cov=covariancia, low = c(medias[1]-4*dpadrao[1], 0.01), up = c(medias[1]+4*dpadrao[1], 0.13))

#Exemplo Michaelis-Menten
del.opt.bayes <- del.bayes.cont(expr="p1*t/(p2+t)", p=c("p1", "p2"), fat = c("t"), a=0, b=1.1, N=2, OTM="D", Ni=25, param=param, parallel = TRUE, metodo = "cubature", Np=50)

#Variancia padrao exemplo anterior
std.var(del.opt = del.opt.bayes, n=100)

#Exemplo 2 - Michaelis-Menten com p1 linearizado

#Distribuicao normal p1 local (MM)
medias <- 0.072
dpadrao <- 0.05
variancia <- dpadrao^2
param.p1loc <- normal.bayes(p = c("p2"), mu = medias, cov = variancia, low=0.01, up=0.13)

#Exemplo Michaelis-Menten (p1 loc)
del.opt.bayes <- del.bayes.cont(expr="p1*t/(p2+t)", p=c("p1", "p2"), fat = c("t"), a=0, b=1.1, N=2, OTM="D", Ni=25, param=param.p1loc, p.loc = list(p1=216), p.bayes=c("p2"), parallel = TRUE, metodo = "cubature", Np=100)

#Variancia padrao exemplo anterior
std.var(del.opt = del.opt.bayes, n=100)

# Exemplo 3 - Reacoes consecutivas

#Distribuicao normal 2
medias.2 <- c(1, 0.5)
#dpadrao <- c(0.01, 0.05)
covariancia.2 <- matrix(c(1, 0.25, 0.25, 1), nrow = 2)
param.2 <- normal.bayes(p = c("p1", "p2"), mu = medias.2, cov=covariancia.2, low = c(0.8, 0.3), up = c(1.2, 0.7))

#Exemplo reacoes-consec
del.opt.bayes <- del.bayes.cont(expr="p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p=c("p1", "p2"), fat = c("t"), a=0, b=10, N=2, OTM="D", Ni=25, param=param.2, parallel = TRUE, metodo = "cubature", Np=100)

#Variancia padrao exemplo anterior
std.var(del.opt = del.opt.bayes, n=100)

#Exemplo 4 - Modelo Compartimental

#Distribuicao uniforme
mu <- c(0.05884, 4.298)
low <- c(mu[1]-0.01, mu[2]-1)
up <- c(mu[1]+0.01, mu[2]+1)
param.uniff <- uniff.bayes(low, up)

#Exemplo compartimental
del.opt.bayes <- del.bayes.cont(expr="p3*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=3, OTM="D", Ni=30, param=param.uniff, p.loc = list(p3=21.80), p.bayes=c("p1", "p2"), parallel = TRUE, metodo = "cubature")

#variancia padrao exemplo anterior
std.var(del.opt = del.opt.bayes, n=500)

#Exemplo 5 - Modelo Compartimental (Dist. 2)

#Distribuicao uniforme
mu <- c(0.05884, 4.298)
low <- c(mu[1]-0.04, mu[2]-4)
up <- c(mu[1]+0.04, mu[2]+4)
param.uniff <- uniff.bayes(low, up)

#Exemplo compartimental
del.opt.bayes <- del.bayes.cont(expr="p3*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, N=6, OTM="D", Ni=100, param=param.uniff, p.loc = list(p3=21.80), p.bayes=c("p1", "p2"), parallel = TRUE, metodo = "cubature")

#variancia padrao exemplo anterior
std.var(del.opt = del.opt.bayes, n=100)



#Exemplo 6 - Michaelis-Menten com p1 linearizado (As-otimalidade)

#Distribuicao normal p1 local (MM)
medias <- 0.072
dpadrao <- 0.05
variancia <- dpadrao^2
param.p1loc <- normal.bayes(p = c("p2"), mu = medias, cov = variancia, low=0.01, up=0.13)

#Exemplo Michaelis-Menten (p1 loc)
del.opt.bayes <- del.bayes.cont(expr="p1*t/(p2+t)", p=c("p1", "p2"), fat = c("t"), a=0, b=1.1, N=2, OTM="As", Ws=c(2,10^6), Ni=25, param=param.p1loc, p.loc = list(p1=216), p.bayes=c("p2"), parallel = TRUE, metodo = "cubature", Np=1000)



# Exemplo 7 - Reacoes consecutivas (As-otimalidade)

#Distribuicao normal 2
medias.2 <- c(1, 0.5)
#dpadrao <- c(0.01, 0.05)
covariancia.2 <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)
param.2 <- normal.bayes(p = c("p1", "p2"), mu = medias.2, cov=covariancia.2, low = c(0.8, 0.3), up = c(1.2, 0.7))

#Exemplo reacoes-consec
del.opt.bayes <- del.bayes.cont(expr="p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p=c("p1", "p2"), fat = c("t"), a=0, b=10, N=2, OTM="As", Ws=c(1/17,1/3), Ni=25, param=param.2, parallel = TRUE, metodo = "MC", Np=5)
