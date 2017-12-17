setwd("D:/ChangePoint/source")
require(Rcpp)
sourceCpp("./slasso.cpp")

n=200
sigma=0.05
p=0
K=12
##Initial Setup
break_point = c(10, 13, 15, 23, 25, 40, 44, 65, 76, 78, 81)*n/100
break_point_ex = c(0,break_point,n)
diff_bp = diff(break_point_ex)
Beta = matrix(c(0  ,40 ,-10 , 20 ,-20 , 30, -12 ,  9  ,52  ,21 , 42  , 0)/20 -0.8,nrow=1,byrow = TRUE)
X = matrix(rnorm(n*p),n,p)
X = cbind(rep(1,n),X)
Y = matrix(0,n,1)

for (bp in 1:K){
  Y = Y + diag(c(rep(0,break_point_ex[bp]),rep(1,diff_bp[bp]),rep(0,n-break_point_ex[bp+1])))%*% X %*% Beta[1:(p+1),bp]
}


##Add Noise
E = matrix(rnorm(n,0,sigma),n,1)
Y = Y + E

Z = matrix(0,n,(n-1)*(p+1))
for (t in 1:(n-1)){
  Z[(t+1):n,(t-1)*(p+1)+1:(p+1)] = X[(t+1):n,]
}
sl = SLasso(Y,Z)


sifs=SIFS(X,Y,n,p,TRUE)

