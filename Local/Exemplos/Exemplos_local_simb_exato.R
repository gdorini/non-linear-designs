source("C:/Users/GustavoJoseDoriniJun/Desktop/Códigos em R 2025/Funções/Local/Função_principal_del_loc_simb_exato.R", encoding = 'UTF-8')

comeco <- Sys.time()

#Michaelis-Menten
#del.opt.e <- del.loc.e(expr = "p1*t/(p2+t)", p=c("p1", "p2"), fat=c("t"), a=0, b=1.1, N=4, Ni=25, p_loc=c(216, 0.072), OTM="D", Ws=NULL)

#modelo com 1 parametro
#Exemplo 1 parametro
#del.opt.e <- del.loc.e(expr = "Exp(-p*t)", p = c("p"), fat = c("t"), a=0, b=5, N = 5, Ni = 50, p_loc = c(1), OTM = "D")

#exemplo reuniao
del.opt.e <- del.loc.e(expr = "(p1*p2*4/(p3*(p1-p2)))*(Exp(-p2*t)-Exp(-p1*t))", p=c("p1", "p2", "p3"), fat=c("t"), a=0, b=25, N=11, Ni=100, p_loc=c(exp(-0.3631), exp(-2.6042), exp(-3.4283)), OTM="D", Ws=NULL)

fim <- Sys.time() - comeco
print(fim)

