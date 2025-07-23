source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Local/Função_principal_del_local_simb.R", encoding = 'UTF-8')

inicio <- Sys.time()

#Reacoes consec.
#del.opt <- del.loc(expr = "p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p = c("p1", "p2"), fat = c("t"), a=0, b=10, NP=2, Ni = 20, p_loc = c(1, 0.5), OTM = "D")

#Michaelis-Menten
#del.opt <- del.loc(expr = "p1*t/(p2+t)", p = c("p1", "p2"), fat = c("t"), a=0, b=1.1, NP = 2, Ni = 50, p_loc = c(216, 0.072), OTM = "D")

#Michaelis-Menten (As)
#del.opt <- del.loc(expr = "p1*t/(p2+t)", p = c("p1", "p2"), fat = c("t"), a=0, b=1.1, NP = 2, Ni = 50, p_loc = c(216, 0.072), OTM = "A", Ws = c(1,1))

#Exemplo 1 parametro
#del.opt <- del.loc(expr = "Exp(-p*t)", p = c("p"), fat = c("t"), a=0, b=5, NP = 1, Ni = 50, p_loc = c(1), OTM = "D")

#Modelo Compartimental
#del.opt <- del.loc(expr = "p3*(Exp(-p2*t)-Exp(-p1*t))", p = c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, NP=3, Ni = 50, p_loc = c(4.29, 0.0589, 21.80), OTM = "D")

#Modelo Compartimental (As)
del.opt <- del.loc(expr = "p3*(Exp(-p2*t)-Exp(-p1*t))", p = c("p1", "p2", "p3"), fat = c("t"), a=0, b=30, NP=3, Ni = 50, p_loc = c(4.29, 0.0589, 21.80), OTM = "As", Ws=c(1, 10^4, 1/2))

#Reacoes consec.
del.opt <- del.loc(expr = "p1*(Exp(-p2*t)-Exp(-p1*t))/(p1-p2)", p = c("p1", "p2"), fat = c("t"), a=0, b=10, NP=2, Ni = 20, p_loc = c(1, 0.5), OTM = "As", Ws=c(1/16,1/3))

#exemplo reuniao
#del.opt <- del.loc(expr = "(p1*p2*4/(p3*(p1-p2)))*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat=c("t"), a=0, b=25, N=3, Ni=50, p_loc=c(exp(-0.3631), exp(-2.6042), exp(-3.4283)), OTM="D", Ws=NULL)

fim <- Sys.time()-inicio
print(fim)

#Gráfico da variância padrão da resposta esperada

plot.var <- var.pad(del.opt)
