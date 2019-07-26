price=read.csv("C:\\Users\\76783\\Desktop\\Data of Raphael\\data_cc_daily_prices.csv",header=T)
Date<-na.omit(price$X)
SP500<-na.omit(price$SP500)
GLD<-na.omit(price$GLD)

length(GLD)
length(SP500)

par(mfrow=c(2,1))
plot(log(GLD))
plot(log(SP500))

#stationary?
adfTest(diff(log(SP500)),lags=12)
adfTest(diff(log(GLD)),lags=12)
acf(diff(log(GLD)))
acf(log(SP500))

##########
# VECM(0)#
##########

end<-ncol(logprice)
dlogprice<-logprice[,2:end]-logprice[,1:(end-1)]

S_11<-logprice[,-end]%*%t(logprice[,-end])
S_00<-dlogprice%*%t(dlogprice)
S_01<-dlogprice%*%t(logprice[,-end])
#E<-eigen(S_01%*%solve(S_11)%*%t(S_01)%*%solve(S_00))

E<-eigen(solve(S_11)%*%t(S_01)%*%solve(S_00)%*%S_01)
E$values

beta=matrix(E$vectors[,1],2,1)
alpha<-S_01%*%beta%*%solve(t(beta)%*%S_01%*%beta)
alpha%*%t(beta)
longrun<-t(beta)%*%as.matrix(logprice[,-end])
adfTest(longrun,lags=10)     #p-value<0.01
kpss.test(as.vector(longrun))

##########
# VECM(3)#
##########

end<-ncol(logprice)
dlogprice<-logprice[,2:end]-logprice[,1:(end-1)]#Delta_X_t
End<-ncol(dlogprice)
#constant
VARselect(as.data.frame(logprice),lag.max=10,type=c("none"))  #AIC:p=2
?VAR

Z_0t<-dlogprice
Z_1t<-logprice[,-end]
constant<-rep(1,end-1)
 Z_2t_1<-cbind(logprice[,1],dlogprice[,-End])#delta_x_t-1
 Z_2t_2<-cbind(0,logprice[,1],dlogprice[,-c(ncol(dlogprice)-1,ncol(dlogprice))])#delta_x_t_2
 Z_2t_3<-cbind(0,0,logprice[,1],dlogprice[,-c(ncol(dlogprice)-2,ncol(dlogprice)-1,ncol(dlogprice))])#delta_x_t_3
 Z_2t_4<-cbind(0,0,0,logprice[,1],dlogprice[,-c(ncol(dlogprice)-3,ncol(dlogprice)-2,ncol(dlogprice)-1,ncol(dlogprice))])#delta_x_t_2
 
#Z_2t<-rbind(Z_2t_1)
Z_2t<-rbind(Z_2t_1,Z_2t_2,Z_2t_3,Z_2t_4)
            #,constant)
#Z_2t<-rbind(Z_2t_1,Z_2t_2,Z_2t_3,Z_2t_4,constant)
#Z_2t<-rbind(cbind(logprice[,1],dlogprice[,-ncol(dlogprice)]),constant)#Delta_X_t-1

M_02<-Z_0t%*%t(Z_2t)/(end-1)
M_12<-Z_1t%*%t(Z_2t)/(end-1)
M_22<-Z_2t%*%t(Z_2t)/(end-1)

R_0t<-Z_0t-M_02%*%solve(M_22)%*%Z_2t
R_1t<-Z_1t-M_12%*%solve(M_22)%*%Z_2t

S_00<-R_0t%*%t(R_0t)/(end-1)
S_01<-R_0t%*%t(R_1t)/(end-1)
S_11<-R_1t%*%t(R_1t)/(end-1)

E<-eigen(solve(S_11)%*%t(S_01)%*%solve(S_00)%*%S_01)
E$values
#0.0026533898 0.0001813212
#beta=cbind(matrix(E$vectors[,1],5,1),matrix(E$vectors[,2],5,1))
beta=matrix(E$vectors[,1],2,1)
# 0.99853896 0.05403653
##2019.0717  ECMvar()
logprice<-as.data.frame(t(logprice)) 

n3<-ECMvar(logprice,4,ibeta = beta,include.const = FALSE)
n3a<-refECMvar(n3,thres=0.5)
MTSdiag(n3a)
alpha<-n3a$alpha
Gamma<-n3a$Phip
#long-run relationship
longrun<-as.vector(t(beta)%*%as.matrix(logprice[,-end]))
adfTest(longrun,lags=12) #<0.01 GLD v.s. VIX
#p-value<0.01, one cointegration
#
png("acf_longrun_R2.png", width=20, height=15, units="cm", res=800, bg="transparent")
acf(longrun)
dev.off()


#################
# Johansen test3#
#################

logprice<-log(rbind(GLD,SP500))
model2<-ca.jo(logprice, ecdet = "none", type="eigen", K=4,spec="longrun")
summary(model2)

#Values of teststatistic and critical values of test:(K=4)
#         test  10pct 5pct  1pct
#r <= 1 | 0.82  6.50  8.18 11.65
#r = 0  | 9.40 12.91 14.90 19.19

#Values of teststatistic and critical values of test:(K=9)
#         test 10pct  5pct  1pct
#r <= 1 | 0.52  6.50  8.18 11.65
#r = 0  | 8.12 12.91 14.90 19.19

