require(lars)
require(MASS)
require(cpm)

## Function
Summary =function(S,Y,X,Z){
    if (is.null(S)){
        coef=lm(Y~X-1)$coefficients
        names(coef)=NULL
        Estimation=list(Break_Point=NULL,Beta=coef[1:(p+1)],diff=NULL,raw_beta=coef[1:(p+1)])
    } else {
        p=dim(X)[2]-1
        n=length(Y)
        coef=lm(Y~cbind(X,Z[,S])-1)$coefficients
        names(coef)=NULL
        
        ls0=p+1
        est_beta=coef[1:ls0]
        #### shrink some coefficient
        idx_shrinkage = abs(est_beta)<0.001
        est_beta[idx_shrinkage]=0
        ####
        est_kesi=rep(0,(n-1)*(p+1))
        est_kesi[S]=coef[-(1:ls0)]
        est_kesi =matrix(est_kesi,p+1,n-1)
        est_idx=rep(FALSE,(n-1)*(p+1))
        est_idx[S]=TRUE
        est_idx =matrix(est_idx,p+1,n-1)
        est_bp=which(apply(est_idx,2,any))
        
        temp_beta=est_beta
        for ( i in 1:length(est_bp)){
            temp_beta = temp_beta + est_kesi[,est_bp[i]]
            est_beta = cbind(est_beta,temp_beta)
        }
        colnames(est_beta)=NULL
        Estimation=list(Break_Point=est_bp,Beta=est_beta,diff=t(diff(t(est_beta))),raw_beta=cbind(coef[1:(p+1)],est_kesi))
    }
    return(Estimation)
}
dist_err=function(BP1,BP2){
    max(sapply(BP1, function(x) min(abs(x-BP2))))    
}
PDR_FDR=function(M1,M2){
    nrow=dim(M1)[1]
    ncol=dim(M1)[2]
    FP=0
    FN=0
    r1=0
    r2=0
    
    for ( i in 1:nrow){
        for ( j in 1:ncol){
            if (M2[i,j] == 0){
                if (M1[i,j]!=0) {
                    FP=FP+1
                    r1=r1+1
                }
            }
            else {
                r2=r2+1
                if (M1[i,j]==0) {
                    FN=FN+1
                } else {
                    r1=r1+1
                }
            }
        }
    }
    output=c(FP/r1,1-FN/r2)
    names(output)=c("FDR","PDR")
    return(output)
}
PDR_bias=function(M1,M2){
    nrow=dim(M1)[1]
    ncol=dim(M1)[2]
    FP=0
    FN=0
    r1=0
    r2=0
    
    for ( i in 1:nrow){
        for ( j in 1:ncol){
            if (M2[i,j] != 0){
                r2=r2+1
                if (M1[i,j]==0 & M1[i,max(j-1,1)]==0 & M1[i,min(j+1,ncol)]==0){
                    FN=FN+1
                } 
            }
        }
    }
    for ( i in 1:nrow){
        for ( j in 1:ncol){
            if (M1[i,j] != 0){
                r1=r1+1
                if (M2[i,j]==0 & M2[i,max(j-1,1)]==0 & M2[i,min(j+1,ncol)]==0){
                    FP=FP+1
                } 
            }
        }
    }
    output=c(FP/r1,1-FN/r2)
    names(output)=c("FDR","PDR")
    return(output)
}
makeTable = function(sim_list){
  mean_table=t(sapply(sim_list,colMeans))
  row.names(mean_table)=names(sim_list)
  print(round(mean_table,4))
  sd_table=t(sapply(sim_list,function(x) apply(x,2,sd)))
  row.names(sd_table)=names(sim_list)
  print(round(sd_table,4))
}

## SLasso Detection
SIFS = function(X,Y,n,p,normalize=TRUE){
    EBIC= function(S,projH){
        full=(n-1)*(p+1)
        gamma = (1-log(n)/(2*log(full)))
        return( n*log( sum((projH%*%Y)^2) /n)+length(S)*log(n)+2*gamma*log(choose(full,length(S))) )
    }
    
    projection_H = function(HH,update_idx){
        newHH=HH
        for (ii in update_idx){
            z = as.matrix(Z[,ii],n,1)
            newHH = newHH %*% (diag(1,n)- (z %*% t(z) %*% newHH) / as.numeric (t(z)%*% newHH %*% z) )
        }
        return(newHH)
    }
    
    Z = matrix(0,n,(n-1)*(p+1))
    for (t in 1:(n-1)){
        Z[(t+1):n,(t-1)*(p+1)+1:(p+1)] = X[(t+1):n,]
    }
    
    
    ## Standarized##
    if(normalize){
        X[,-1]=scale(X[,-1])
        Z=scale(Z)
    }
    ################
    S1 = NULL
    SA=1:(n-1)*(p+1)
    ebic=Inf
    H = diag(1,n)- X%*%solve(t(X)%*%X)%*%t(X)
    
    repeat{
        ybar = H%*%Y
        R=0
        for ( i in SA){
            tempr = abs(t(Z[,i])%*%ybar)
            if(tempr > R){
                R=tempr
                temp_s=i
            } 
        }
        new_S=c(S1,temp_s)
        new_S=new_S[order(new_S)]
        new_H= projection_H(H, temp_s)
        
        rm(tempr)
        
        new_ebic=EBIC(new_S,new_H)
        
        if(new_ebic>ebic | length(S1)==(n-1)*(p+1)){
            break
        } else {
            ebic=new_ebic
            S1=new_S
            H=new_H
            SA=SA[SA!=temp_s]
        }
    }
    rm(i,R)
    rm(temp_s,new_S,new_ebic,ebic,SA,ybar)
    estimation=Summary(S1,Y,X,Z)
}
SBFS = function(X,Y,n,p){
    EBIC= function(S,projH,v){
        projection  = sum((projH%*%Y)^2)
        gamma = (1-log(n)/(2*log((n-1)*(p+1))))
        if (length(v) == 0){
            return( n*log(projection/n)+(length(S))*log(n) )
        } else { 
            return( n*log(projection/n)+(length(S)+2*gamma*length(v))*log(n)+2*gamma*sum(sapply( v,function(x) log(choose(p+1,x)) )) )
        }
    }
    
    projection_H = function(HH,update_idx){
        newHH=HH
        for (ii in update_idx){
            z = as.matrix(Z[,ii],n,1)
            newHH = newHH %*% (diag(1,n)- (z %*% t(z) %*% newHH) / as.numeric (t(z)%*% newHH %*% z) )
        }
        return(newHH)
    }
    
    
    Z = matrix(0,n,(n-1)*(p+1))
    for (t in 1:(n-1)){
        Z[(t+1):n,(t-1)*(p+1)+1:(p+1)] = X[(t+1):n,]
    }
    
    ## Standarized##
    #X[,-1]=scale(X[,-1])
    Z=scale(Z)
    ################
    S = NULL
    TA=1:(n-p-1)
    v = NULL
    ebic=Inf
    
    H = diag(1,n)- X%*%solve(t(X)%*%X)%*%t(X)
    repeat{
        ybar = H%*%Y

        R2=0
        for ( i in TA){
            Zt = Z[,((i-1)*(p+1)+1):(i*(p+1))]
            tempR2 = t(ybar)%*%Zt%*%solve(t(Zt)%*%Zt)%*%t(Zt)%*%ybar
            if(tempR2 > R2){
                R2=tempR2
                temp_t=i
            } 
        }
        ST = ((temp_t-1)*(p+1)+1):(temp_t*(p+1))
        rm(tempR2)
        new_S=c(S,ST)
        new_S=new_S[order(new_S)]
        new_H=projection_H(H,ST)
        v = c(v,p+1)
        inter_ebic=EBIC(new_S,new_H,v)
        
        repeat{## deduction current block
            tebic=Inf
            tempv=v
            tempv[length(tempv)]=tempv[length(tempv)]-1
            for(i in ST){
                temp_new_S=new_S[new_S!=i]
                temp_ebic = EBIC(temp_new_S,projection_H(H,ST[ST!=i]),tempv)
                if(temp_ebic<tebic){
                    tebic=temp_ebic
                    temp_rms=i
                }
            }
            
            if(tebic<inter_ebic){
                ST = ST[ST!=temp_rms]
                new_S = new_S[new_S!=temp_rms]
                new_H = projection_H(H,ST)
                v=tempv
                inter_ebic=EBIC(new_S,new_H,v)
            } else {
                break    
            }
        }
        

        new_ebic=inter_ebic
        
        if(new_ebic < ebic){
            ebic=new_ebic
            S=new_S
            TA=TA[TA!=temp_t]
            H=new_H
        } else {
            break
        }
    }
    rm(i,v,R2)      
    rm(temp_rms,temp_t,inter_ebic,temp_ebic,new_ebic,ebic,tebic,temp_new_S,new_S,ST,TA,ybar)
    estimation=Summary(S,Y,X,Z)
}
## LS_TV
LS_TV= function(Y,n,K_max,K_chs){
    Z = matrix(0,n,n)
    for (t in 1:n){
        Z[t:n,t] = 1 
    }
    if (K_chs>0){
        ### LAR
        lar = lars(Z,Y,type='lar',intercept=FALSE,max.steps=K_max+1)
        cf_max = coef(lar)[K_max+2,]
        cp_max = which(cf_max!=0) [-1]
        
        ### DP
        path = matrix(NA,(K_max-K_chs+1),K_chs)
        ####step 0
        tmp = sapply(cp_max[1:(K_max-K_chs+1)], function(x) var(Y[1:x])*(x-1))
        path[,1] = 1: (K_max-K_chs+1)
        ####step 1 to K_chs-1
        if ( K_chs > 1 ){
            for ( k in 2:K_chs-1){
                new_path= matrix(NA,(K_max-K_chs+1), k+1)
                new_tmp=rep(0,(K_max-K_chs+1))
                for(x in 1:(K_max-K_chs+1) ) {
                    L = cp_max[x+k]
                    min_tmp=tmp[1:x] +sapply(cp_max[k:(x+k-1)], function(y) ifelse(y<L-2,var(Y[y:(L-1)])*(L-y-1),0) )
                    new_tmp[x]=min( min_tmp )
                    
                    new_path[x,]=c(path[which(new_tmp[x]==min_tmp)[1],1:k],x+k)
                }
                tmp=new_tmp
                path[,1:(k+1)]=new_path
            }
        }
        ####step K_chs
        min_tmp=tmp +sapply(cp_max[K_chs:K_max], function(y) ifelse(y<n,var(Y[y:n])*(n-y),0) )
        final_path=path[which(min(min_tmp)==min_tmp),]
        bp = cp_max[final_path]-1
        #LS = min(min_tmp)
    } else if (K_chs == 0) {
        bp=NULL
        #LS = var(Y)*(n-1)
    } else if (K_chs < 0) {
        ### LAR
        lar = lars(Z,Y,type='lar',intercept=FALSE,max.steps=K_max+1)
        cf_max = coef(lar)[K_max+2,]
        cp_max = which(cf_max!=0) [-1]
        
        
        library(snow)
        inputs <- 1:K_max
        rDP <- function(K_pt,K_max,cp_max) {
            path = matrix(NA,(K_max-K_pt+1),K_pt)
            ####step 0
            tmp = sapply(cp_max[1:(K_max-K_pt+1)], function(x) var(Y[1:x])*(x-1))
            path[,1] = 1: (K_max-K_pt+1)
            ####step 1 to K_chs-1
            if ( K_pt > 1 ){
                for ( k in 2:K_pt-1){
                    new_path= matrix(NA,(K_max-K_pt+1), k+1)
                    new_tmp=rep(0,(K_max-K_pt+1))
                    for(x in 1:(K_max-K_pt+1) ) {
                        L = cp_max[x+k]
                        min_tmp=tmp[1:x] +sapply(cp_max[k:(x+k-1)], function(y) ifelse(y<L-2,var(Y[y:(L-1)])*(L-y-1),0) )
                        new_tmp[x]=min( min_tmp )
                        
                        new_path[x,]=c(path[which(new_tmp[x]==min_tmp)[1],1:k],x+k)
                    }
                    tmp=new_tmp
                    path[,1:(k+1)]=new_path
                }
            }
            ####step K_chs
            min_tmp=tmp +sapply(cp_max[K_pt:K_max], function(y) ifelse(y<n,var(Y[y:n])*(n-y),0) )
            LS = min(min_tmp)
            final_path=path[which(LS==min_tmp)[1],]
            return(list(LS=LS,PATH=final_path))
        }
        numCores = 7 
        cl <- makeCluster(numCores)
        results = clusterApplyLB(cl, inputs, rDP,K_max=K_max,cp_max=cp_max)
        stopCluster(cl)
        
        ratio=1-0.05
        min_ratio = Inf
        for ( i in 1:(K_max-1)){
            tmp_ratio = results[[i+1]]$LS/results[[i]]$LS
            if (tmp_ratio>ratio & tmp_ratio<min_ratio) {
                min_ratio=tmp_ratio
                idx=i
            }    
        }
        final_path=results[[idx]]$PATH
        bp = cp_max[final_path]-1
    }
    ### final 
    estimation=Summary(bp,Y,matrix(1,n,1),Z[,-1])
    return(estimation)
}
## https://www.r-bloggers.com/change-point-detection-in-time-series-with-r-and-tableau/ other methods
CPM  = function(Y){
  Z = matrix(0,n,n)
  for (t in 1:n){
    Z[t:n,t] = 1 
  }
  
  bp = NULL
  y = Y
  last_bp = 0
  repeat{
    cpm_tmp = detectChangePoint(y,"Student")
    if (cpm_tmp$changeDetected ){
      last_bp = last_bp + cpm_tmp$detectionTime-1
      bp = c(bp, last_bp)
      y = Y[ (last_bp+1):n ]
    } else {
      break
    }	
  }
  
  estimation=Summary(bp,Y,matrix(1,n,1),Z[,-1])
  return(estimation)
}

#Y = matrix(0,n,1)
# Beta = NULL
# for (bp in 1:K){
#     beta = runif(p+1,bp,1+bp)  # random beta
#     Y = Y + diag(c(rep(0,break_point[bp]),rep(1,diff_bp[bp]),rep(0,n-break_point[bp+1])))%*% X %*% beta
#     Beta = cbind(Beta,beta)
# }


#############################################   No Covariates   ###################################################
##Settings
n=200
sigma=0.05
p=0
K=12
##Initial Setup
break_point = c(10, 13, 15, 23, 25, 40, 44, 65, 76, 78, 81)*n/100
break_point_ex = c(0,break_point,n)
diff_bp = diff(break_point_ex)
Beta = matrix(c(0  ,40 ,-10 , 20 ,-20 , 30, -12 ,  9  ,52  ,21 , 42  , 0)/20 -0.8,nrow=1,byrow = TRUE)
raw_Beta=matrix(0,p+1,n)
raw_Beta[,c(0,break_point)+1]=Beta[1:(p+1),]
result_SIFS_nor=NULL
result_SIFS=NULL
result_LSTV_1=NULL
result_LSTV_2=NULL
for (i in 1:100){
    X = matrix(rnorm(n*p),n,p)
    X = cbind(rep(1,n),X)
    Y = matrix(0,n,1)
    ##Setting 1
    for (bp in 1:K){
        Y = Y + diag(c(rep(0,break_point_ex[bp]),rep(1,diff_bp[bp]),rep(0,n-break_point_ex[bp+1])))%*% X %*% Beta[1:(p+1),bp]
    }
    ##Setting 2
    #idx=sample(1:12,n,replace=TRUE,prob=c(0.1, 0.03, 0.02, 0.08, 0.02, 0.15, 0.04, 0.21, 0.11, 0.02, 0.03,0.19))
    #idx=idx[order(idx)]
    #break_point=which(diff(idx)==1)
    #raw_Beta=matrix(0,p+1,n)
    #raw_Beta[,c(0,break_point)+1]=Beta[1:(p+1),]
    #    for ( j in 1:n){
    #    Y[j]=X[j,]%*% Beta[idx[j]]
    #}
    
    ##Add Noise
    E = matrix(rnorm(n,0,sigma),n,1)
    Y = Y + E
    
    ##Evaluation
    estimation_SIFS=SIFS(X,Y,n,p,TRUE)
    estimation_LSTV_1=LS_TV(Y,n,3*(K-1),K-1)
    estimation_LSTV_2=LS_TV(Y,n,3*(K-1),-1)
    result_SIFS=rbind(result_SIFS,c(
        dist_err(break_point,estimation_SIFS$Break_Point)/n,
        dist_err(estimation_SIFS$Break_Point,break_point)/n,
        PDR_FDR(estimation_SIFS$raw_beta,raw_Beta),
        PDR_bias(estimation_SIFS$raw_beta,raw_Beta))
    )
    result_LSTV_1=rbind(result_LSTV_1,c(
        dist_err(break_point,estimation_LSTV_1$Break_Point)/n,
        dist_err(estimation_LSTV_1$Break_Point,break_point)/n,
        PDR_FDR(estimation_LSTV_1$raw_beta,raw_Beta),
        PDR_bias(estimation_LSTV_1$raw_beta,raw_Beta))
    )
    result_LSTV_2=rbind(result_LSTV_2,c(
        dist_err(break_point,estimation_LSTV_2$Break_Point)/n,
        dist_err(estimation_LSTV_2$Break_Point,break_point)/n,
        PDR_FDR(estimation_LSTV_2$raw_beta,raw_Beta),
        PDR_bias(estimation_LSTV_2$raw_beta,raw_Beta))
    )
    print(i)
}
sim_list2=list(SIFS_normalized=result_SIFS_nor,LS_TV1=result_LSTV_1,LS_TV2=result_LSTV_2)
makeTable(sim_list2)



####################################  Covariates Involved   ##############################################
n=1000
sigma=0.1
p=1
K=3
break_point = c(25,80)*n/100
break_point_ex = c(0,break_point,n)
diff_bp = diff(break_point_ex)
Beta = cbind(
    c(  1,   1, 1, 0, 0,  1,   2, -1, 0, 0, 1.2, 0),    
    c(1.2, 1.1, 1, 0, 0, -1, 1.8,  1, 0, 0,   1, 0),
    c(1.8, 1.2, 1, 0, 0,  0, 1.6,  0, 0, 0, 0.8, 0)
)
raw_Beta=matrix(0,p+1,n)
raw_Beta[,1]=Beta[1:(p+1),1]
raw_Beta[,break_point+1]=t(diff(t(Beta[1:(p+1),])))

result_SIFS=NULL
result_SBFS=NULL
for (i in 1:100){
    X = matrix(rnorm(n*p),n,p)
    X = cbind(rep(1,n),X)
    Y = matrix(0,n,1)
    ##Setting 1
    for (bp in 1:K){
        Y = Y + diag(c(rep(0,break_point_ex[bp]),rep(1,diff_bp[bp]),rep(0,n-break_point_ex[bp+1])))%*% X %*% Beta[1:(p+1),bp]
    }
    E = matrix(rnorm(n,0,sigma),n,1)
    Y = Y + E
    
    ##Evaluation
    estimation_SIFS=SIFS(X,Y,n,p)
    estimation_SBFS=SBFS(X,Y,n,p)
    
    result_SIFS=rbind(result_SIFS,c(
        dist_err(break_point,estimation_SIFS$Break_Point)/n,
        dist_err(estimation_SIFS$Break_Point,break_point)/n,
        PDR_FDR(estimation_SIFS$raw_beta,raw_Beta))
    )
    result_SBFS=rbind(result_SBFS,c(
        dist_err(break_point,estimation_SBFS$Break_Point)/n,
        dist_err(estimation_SBFS$Break_Point,break_point)/n,
        PDR_FDR(estimation_SBFS$raw_beta,raw_Beta))
    )
    print(i)
}

result=result_SIFS
colMeans(result)
apply(result,2,sd)
result=result_SBFS
colMeans(result)
apply(result,2,sd)