
## This is the main workhorse for the tensor decompositions, and the wrapper 
## function RunSDA4D is provided for ease of use.  In particular providing 
## hyperparameters for the priors specified.

SparsePARAFAC<-function(params,
                        datatensor,
                        maxiter,
                        stopping=TRUE,
                        track=1,
                        debugging=FALSE){
    
    # port of matlab code to R, in 2016.
    
    #### readindata ####
    #datatensor<-readRDS(paste0(datapath,'/datatensor.RDS'))
    #system(paste0('mkdir -p ',outDIR))
    # init_function ####
    initialise_4D <- function(params){
        list_of_vars <- list()
        list_of_vars$Error=c(0)
        list_of_vars$Neg_FE=c(0)
        list_of_vars$A <- list(mu = matrix(rnorm(params$N * params$C),params$N,params$C),
                               precision = matrix(100,params$N,params$C))
        list_of_vars$A$mom2 <- list_of_vars$A$mu^2 + 0.01
        if(params$T>1){
            list_of_vars$B <- list(nu = matrix(rnorm(params$T * params$C),params$T,params$C),
                                   precision = matrix(100,params$T,params$C))
            list_of_vars$B$mom2 <- list_of_vars$B$nu^2 + 0.01
        }else{
            list_of_vars$B <- list(nu = matrix(c(1),params$T,params$C),
                                   precision = matrix(1,params$T,params$C))
            list_of_vars$B$mom2 <- list_of_vars$B$nu^2
        }
        if(params$M>1){
            list_of_vars$D <- list(alpha = matrix(rnorm(params$M * params$C),params$M,params$C),
                                   precision = matrix(100,params$M,params$C))
            list_of_vars$D$mom2 <- list_of_vars$D$alpha^2 + 0.01
        }else{
            list_of_vars$D <- list(alpha = matrix(c(1),params$M,params$C),
                                   precision = matrix(1,params$M,params$C))
            list_of_vars$D$mom2 <- list_of_vars$D$alpha^2
        }

        list_of_vars$Lam <- list(u = matrix(params$u,params$L,params$T),
                                 v = matrix(params$v,params$L,params$T),
                                 mom1 = matrix(params$u * params$v,params$L,params$T))

        list_of_vars$Beta <- list(e = matrix(params$e,params$C,1),
                                  f = matrix(params$f,params$C,1),
                                  mom1 = matrix(params$e * params$f,params$C,1))

        list_of_vars$Ph <- matrix(0.5,params$C,params$L)

        list_of_vars$Ps <- matrix(0.5,params$C,params$L)

        list_of_vars$Rho <- matrix(rbeta(params$C,params$r,params$z),1,params$C)

        #list_of_vars$S <- list(gamma = matrix(0.5,params$C,params$L))

        list_of_vars$WS <- list(gamma = matrix(0.5,params$C,params$L),
                                mom1 = matrix(0,params$C,params$L),
                                mom2 = matrix(0,params$C,params$L),
                                sigma = matrix(c(100),params$C,params$L),
                                m = matrix(list_of_vars$Beta$e *list_of_vars$Beta$f,params$C,params$L))

        list_of_vars$WS$mom1 = list_of_vars$WS$m * list_of_vars$WS$gamma

        list_of_vars$WS$mom2 = (list_of_vars$WS$m**2 + 1/list_of_vars$WS$sigma) * list_of_vars$WS$gamma

        list_of_vars$Lam = list(u=matrix(c(params$u),params$L,params$T),
                                v=matrix(c(params$v),params$L,params$T),
                                mom1=list_of_vars$Lam$u*list_of_vars$Lam$v)
        return(list_of_vars)
    }
    #### initialise variables ####
    vars<-initialise_4D(params)
    continue=TRUE
    iteration=1

    ## local functions ####
    ## FE ####
    Free_Energy<-function(params){
        #returns negative free energy up to the constant terms
        # requires vector of parameters...
        FE=0
        A2=apply(vars$A$mom2,2,sum) #  a 1 by C vector
        for(t in 1:params$T){
            summation=rep(0,params$L)
            for(m in 1:params$M){
                datamt<-datatensor[,,m,t]
                #   First component = sum_n (y_{nlmt} - sum_c (<a_{nc}><b_{tc}><d_{mc}><x_{cl}>))^2 a vector of length L
                #  Second component = sum_c <a_{nc}^2><b_{tc}^2><d_{mc}^2><x_{cl}^2> a vector of length L
                #  Third component = sum_c(sum_n <a_{nc}>^2)<b_{tc}>^2<d_{mc}>^2 <x_{cl}>^2
                FirstComponent=apply((datamt - vars$A$mu %*% diag(vars$B$nu[t,]*vars$D$alpha[m,]) %*% vars$WS$mom1)**2,2,sum)
                SecondComponent= (A2 * vars$B$mom2[t,]*vars$D$mom2[m,])%*% vars$WS$mom2
                ThirdComponent=(apply(vars$A$mu**2,2,sum) * (vars$B$nu[t,]**2) * (vars$D$alpha[m,]**2)) %*% (vars$WS$mom1**2)
                summation = summation + (FirstComponent +SecondComponent - ThirdComponent) # vector of length L
            }
            FE = FE + 0.5*sum(params$M * params$N * (digamma(vars$Lam$u[,t])+log(vars$Lam$v[,t]))) - 0.5*summation %*% vars$Lam$mom1[,t]
        }

        ## now the terms from the prior and approx posteriors
        #A,B,D
        FE = FE - 0.5*sum( vars$A$mom2 + log(abs(vars$A$precision)) ) + (params$N * params$C)/2
        FE = FE - 0.5*sum( vars$B$mom2 + log(abs(vars$B$precision)) ) + (params$T * params$C)/2
        FE = FE - 0.5*sum( vars$D$mom2 + log(abs(vars$D$precision)) ) + (params$M * params$C)/2
        #5
        # Wmom2 = <S_cl>(1/sigma_{cl} + m_{cl}^2) + (1-<S_{cl}>)<Beta_c> # C by L matrix
        Wmom2 = vars$WS$gamma * ((1/vars$WS$sigma) + (vars$WS$m)**2) + (1-vars$WS$gamma) * matrix(1/vars$Beta$mom1,params$C,params$L)
        FE = FE + 0.5*sum(
            matrix(digamma(vars$Beta$e) + log(vars$Beta$f),params$C,params$L) - matrix(vars$Beta$mom1,params$C,params$L)*Wmom2
        )

        #6 #checked
        FE = FE + sum(
            -0.5*(
                vars$WS$gamma * log(vars$WS$sigma) + (1-vars$WS$gamma)*log(matrix(vars$Beta$mom1,params$C,params$L))
            )
        )

        #7 Rho # checked
        FE = FE + sum(
            (params$r - 1)*log(vars$Rho) + (params$z - 1)*log(1-vars$Rho)
        )

        #8 psi # checked
        FE = FE + sum(
            (params$g - 1) * log(vars$Ps) + (params$h - 1)*log(1 - vars$Ps)
        )

        #9 # checked
        FE = FE + sum(
            vars$Ph * matrix(log(vars$Rho),params$C,params$L) + (1-vars$Ph)*matrix(log(1-vars$Rho),params$C,params$L)
        )

        #10 - s # checked
        FE = FE + sum(
            vars$WS$gamma*log(vars$Ph*vars$Ps) + (1-vars$WS$gamma)*log(1-vars$Ph*vars$Ps)
        )
        # checked
        Xtmp<- (-(1-vars$WS$gamma)*log(1-vars$WS$gamma) - vars$WS$gamma*log(vars$WS$gamma))
        if(any(vars$WS$gamma==0 | vars$WS$gamma==1)){
            Xtmp[vars$WS$gamma==0 | vars$WS$gamma==1]=0
        }
        FE = FE + sum(Xtmp)


        #11 beta_c # checked
        FE = FE + sum(
            params$e * log(abs(vars$Beta$f)) + (params$e - vars$Beta$e)*digamma(vars$Beta$e) +
                vars$Beta$e - vars$Beta$mom1/params$f + lgamma(vars$Beta$e)
        )

        #12 - lambda # checked
        FE = FE + sum(
            params$u*log(vars$Lam$v) + (params$u - vars$Lam$u)*digamma(vars$Lam$u) +
                vars$Lam$u - vars$Lam$u*vars$Lam$v/params$v + lgamma(vars$Lam$u)
        )
        return(FE)
    }
    # UPDATES ####

    # B is T by C
    # D is M by C
    # WS is C by L
    # Lam is L by T

    updateA=function(params){
        A=vars$A


        ## complexity NC + CLT +TC +TC +MC + C ## so linear in N,C,L,T,M

        A$precision = matrix(1,params$N,params$C) +
            t(matrix(
                colSums(t(vars$WS$mom2 %*% vars$Lam$mom1)*vars$B$mom2) *
                    colSums(vars$D$mom2),params$C,params$N))


        #term 1

        ## complexity NM*(LT + C + C + CLT + CT + CT) - linear in all N,L,M,T,C


        tmpmat1 = matrix(0,params$N,params$C) # N by C matrix

        for (m in 1:params$M){
            for (n in 1:params$N){
                datatmp = datatensor[n,,m,]*vars$Lam$mom1
                tmpmat1[n,] = tmpmat1[n,] + vars$D$alpha[m,]*rowSums(t(vars$B$nu)*(vars$WS$mom1%*%datatmp))
            }
        }



        #term 2

        ## initial complexity C*( (C-1)L^2 + (c-1)LT +(c-1)T +(c-1)T +M+ NC + N(c-1) +(c-1) + (N))
        ## so currently quadratic in L and C, linear in N,T,M

        # now its C(KL + KLT + KT +KT +KT ) so linear in LT, quadratic in C

        for (c in 1:params$C){
            Xc = vars$WS$mom1[c,] # length L
            Bc = vars$B$nu[,c] # length T

            ## this next part is the term quadratic in L because of the diag(Xc) term - lets remove it.
            XcXncLamB = colSums(
                t( t(Xc*t(vars$WS$mom1[-c,]))%*%vars$Lam$mom1)*
                    (matrix(Bc,params$T,params$C-1)*vars$B$nu[,-c])
            ) # length C-1


            Dc=vars$D$alpha[,c] # length M
            Dc=t(Dc)%*%vars$D$alpha[,-c] # length C-1

            A$mu[,c] = (tmpmat1[,c] - A$mu[,-c]%*%t(Dc*XcXncLamB))/A$precision[,c]
        }
        A$mom2 = 1/A$precision + A$mu^2
        return(A)
    }
    ## currently quadratic in L and C, linear in N,T,M
    ## now linear in N,L,M,T and quadratic in C


    updateB=function(params){
        B=vars$B
        ##precision:

        ## complexity:
        ## TC + NC + MC +C + CLT + C^2T
        ### now its linear in everything

        B$precision = matrix(1,params$T,params$C) +
            t((colSums(vars$A$mom2)*colSums(vars$D$mom2))*(vars$WS$mom2 %*% vars$Lam$mom1))
        ## term1

        ## complexity: MT*(C +C + CNL + CL + C) so linear in all terms

        tmpmat1 = matrix(0,params$T,params$C)
        for(t in 1:params$T){
            for(m in 1:params$M){
                datamt=datatensor[,,m,t]
                tmpmat1[t,] = tmpmat1[t,] + vars$D$alpha[m,]*t(((t(vars$A$mu)%*%datamt)*vars$WS$mom1)%*%vars$Lam$mom1[,t])
            }
        }
        ##term 2

        ## complexity: C*( L(C-1) + N(C-1) + N(C-1) + C + M(C-1) + M(C-1)  + (C-1)L + (c-1)LT +T(C-1))
        ## linear in L,N,M,T, but quadratic in C.


        tmpmat2 = matrix(0,params$T,params$C)
        for(c in 1:params$C){
            Xc = vars$WS$mom1[c,] # a vector length params$L
            Xc = t(matrix(Xc,params$L,params$C-1))*vars$WS$mom1[-c,] ## C-1 by L matrix
            Ac = vars$A$mu[,c]  # length N
            Dc = vars$D$alpha[,c] # length M
            DAc = colSums(matrix(Ac,params$N,params$C-1)*vars$A$mu[,-c])*
                colSums(matrix(Dc,params$M,params$C-1)*vars$D$alpha[,-c])  ## length C-1
            tmpmat2[,c] = colSums(
                ((matrix(DAc,params$C-1,params$L)*Xc)%*%vars$Lam$mom1)*
                    t(B$nu[,-c])
            )  ## length T
            B$nu[,c] = (tmpmat1[,c] - tmpmat2[,c])/B$precision[,c]
        }
        B$mom2 = 1/B$precision + B$nu^2
        return(B)
    }## linear in L,N,M,T, but quadratic in C.


    updateD=function(params){
        D=vars$D
        ##precision:

        ## complexity: MC + NC + C + CLT +CT + CT  - linear in MNCLT

        D$precision = matrix(1,params$M,params$C) +
            t(matrix(
                colSums(vars$A$mom2)*
                    colSums(vars$B$mom2*t(vars$WS$mom2%*%vars$Lam$mom1)),params$C,params$M))
        ## term1:

        ## complexity: MT( C + CNL + CL  +CL + C ) - linear in C, N,L,M,T

        tmpmat1 = matrix(0,params$M,params$C)
        for(m in 1:params$M){
            for (t in 1:params$T){
                #datamt = Y[,,m,t]
                tmpmat1[m,] = tmpmat1[m,] + vars$B$nu[t,] * #T by C
                    t(((t(vars$A$mu) %*% datatensor[,,m,t])*vars$WS$mom1) %*% vars$Lam$mom1[,t])
            }
        }
        ## term 2
        ## complexity: C( T(C-1) + (C-1)L + N(C-1) + M(C-1) + LT + (C-1)T +M(C-1))  so quadratic in C, linear in L,M,N,T

        for (c in 1:params$C){
            Bc = vars$B$nu[,c]  # length T
            Bc = matrix(Bc,params$T,params$C-1)*vars$B$nu[,-c] #T by C-1
            Xc = vars$WS$mom1[c,] # length L
            Xc = t(matrix(Xc,params$L,params$C-1))*vars$WS$mom1[-c,] #C-1 by L
            AcAnc = vars$A$mu[,c] # length N
            AcAnc = apply(
                matrix(AcAnc,params$N,params$C-1)*vars$A$mu[,-c]
                ,2,sum) ## length C-1
            tmpmat2 = D$alpha[,-c] %*% (AcAnc*(colSums( t(Xc%*%vars$Lam$mom1)*Bc)))
            ##
            D$alpha[,c] = (tmpmat1[,c] - tmpmat2)/D$precision[,c]

        }
        D$mom2 = 1/D$precision + D$alpha^2
        return(D)
    }## so quadratic in C, linear in L,M,N,T


    updateBeta=function(params){
        Beta=vars$Beta
        Wmom2 = vars$WS$gamma * (1/vars$WS$sigma  + vars$WS$m^2) +
            (1-vars$WS$gamma)*(matrix(1/Beta$mom1,params$C,params$L))

        Beta$e = (params$e + params$L/2)*matrix(1,params$C,1)

        for (c in 1:params$C){
            Beta$f[c] = 1/(1/params$f + 0.5*sum(Wmom2[c,]))
        }
        Beta$mom1 = Beta$e*Beta$f

        return(Beta)
    } ## linear in C and L


    updateLam=function(params){
        Lam=vars$Lam
        Lam$u = (params$u + 0.5*params$N*params$M)*matrix(1,params$L,params$T)

        A2 = colSums(vars$A$mom2)
        for(t in 1:params$T){
            summation=rep(0,params$L)
            for(m in 1:params$M){
                FirstComponent = colSums(
                    (datatensor[,,m,t] - vars$A$mu %*% diag(vars$B$nu[t,]*vars$D$alpha[m,])%*%vars$WS$mom1)^2)
                SecondComponent = (A2 * vars$B$mom2[t,]*vars$D$mom2[m,])%*%vars$WS$mom2
                ThirdComponent = (colSums(vars$A$mu^2)*(vars$B$nu[t,]^2)* (vars$D$alpha[m,]^2))%*%(vars$WS$mom1^2)
                summation = summation + (FirstComponent + SecondComponent - ThirdComponent)
            }
            Lam$v[,t] = 1.0/(1.0/params$v + 0.5*summation)
        }
        Lam$mom1 = Lam$u * Lam$v
        return(Lam)
    } ## complexity: LT + NC + TM( NL + C + C^2L + NCL  +NL  + C + C + CL + NC + NC +C + C +C +C + CL +CL + L + )
    ## so linear in N,L,M,T and quadratic in C


    updateRho=function(params){
        Rho=vars$Rho
        for(c in 1:params$C){
            Rho[c] = (params$r - 1 + sum(vars$Ph[c,]))/(params$L + params$r + params$z -2)
        }
        return(Rho)
    } # linear in C and L



    updateWS=function(params){
        WS=vars$WS

        ## precision

        ## complexity: CL + NC + MC + C + LTC - linear in C,N,M,T,L

        WS$sigma = matrix(vars$Beta$mom1,params$C,params$L) +
            matrix(colSums(vars$A$mom2)*colSums(vars$D$mom2),params$C,params$L)*
            t(vars$Lam$mom1%*%vars$B$mom2)
        ## term1
        ## complexity: MT(C + CL +CNL) - linear in M,T,N,L,C

        tmpmat1=matrix(0,params$C,params$L)
        for (m in 1:params$M){
            for(t in 1:params$T){
                tmpmat1 = tmpmat1 +
                    matrix(vars$B$nu[t,]*vars$D$alpha[m,],params$C,params$L)*
                    (t(vars$A$mu)%*%datatensor[,,m,t])*
                    t(matrix(vars$Lam$mom1[,t],params$L,params$C))
            }
        }
        ## term2
        ## complexity: C( T(C-1) + LT(C-1) +  N(C-1) + N + M(C-1) + (C-1) + L(C-1))
        tmpmat2 = matrix(0,params$C,params$L)
        for(c in 1:params$C){
            Bc = vars$B$nu[,c] # length T
            Bc = matrix(Bc,params$T,params$C-1)*vars$B$nu[,-c]  # T by C-1
            LamB = t(vars$Lam$mom1 %*% Bc) # C-1 by L
            Ac = vars$A$mu[,c] # length N
            Ac = colSums(matrix(Ac,params$N,params$C-1)*vars$A$mu[,-c]) # length C-1
            Dc = vars$D$alpha[,c] # length M
            Dc = colSums(matrix(Dc,params$M,params$C-1)*vars$D$alpha[,-c]) #length C-1
            tmpmat2[c,] = colSums(WS$mom1[-c,]*LamB*matrix(Ac*Dc,params$C-1,params$L))

            ##

            WS$m[c,] = (tmpmat1[c,] - tmpmat2[c,])/WS$sigma[c,]

            u = -0.5*log(WS$sigma[c,]) +
                0.5*(WS$m[c,]^2)*
                WS$sigma[c,] +
                log(vars$Ps[c,]*vars$Ph[c,]) +
                0.5*log(vars$Beta$mom1[c]) - log(1-vars$Ps[c,]*vars$Ph[c,])

            WS$gamma[c,] = 1/(1+exp(-u))
            WS$mom1[c,] = WS$m[c,]*WS$gamma[c,]
            WS$mom2[c,] = (1/WS$sigma[c,] + WS$m[c,]^2)*WS$gamma[c,]
        }

        return(WS)
    } ## linear in L,M,N,T, and quadratic in C



    updatePhiPsi=function(params){

        Grad=function(x,y,c,l){
            v1 = vars$WS$gamma[c,l]/x - (1-vars$WS$gamma[c,l])/(1/y -x) + log(vars$Rho[c]) - log(1-vars$Rho[c])
            v2 = vars$WS$gamma[c,l]/y - (1-vars$WS$gamma[c,l])/(1/x - y) + (params$g - 1)/y - (params$h - 1)/(1-y)
            vec=c(v1,v2)
            return(vec)
        }
        Hess = function(X,c,l){
            mat=matrix(0,2,2)
            mat[1,1] = -vars$WS$gamma[c,l]/(X[1]^2) -  (1-vars$WS$gamma[c,l])/((1/X[2] - X[1])^2)
            mat[1,2] = -(1-vars$WS$gamma[c,l])/(1-X[1]*X[2])^2
            mat[2,1] = mat[1,2]
            mat[2,2] = -vars$WS$gamma[c,l]/(X[2]^2) - (1 - vars$WS$gamma[c,l])/((1/X[1] - X[2])^2) - (params$g - 1)/(X[2]^2) - (params$h - 1)/((1-X[2])^2)
            return(mat)
        }
        tildeF = function(X,c,l){
            FE =0
            FE = FE + (vars$WS$gamma[c,l])*log(X[1]*X[2]) + (1-vars$WS$gamma[c,l])*log(1-X[1]*X[2]) +
                (params$g - 1)*log(X[2]) + (params$h -1)*log(1-X[2]) + X[1]*log(vars$Rho[c]) +(1-X[1])*log(1-vars$Rho[c])
            return(FE)
        }
        Phi=vars$Ph
        Psi=vars$Ps
        #use Newton's method for finding Ph, Ps (c,l)
        Xtol = 1e-6
        ftol = 1e-17
        for(c in 1:params$C){
            for (l in 1:params$L){
                tmpY = c(Phi[c,l],Psi[c,l])
                X=tmpY
                g = Grad(tmpY[1],tmpY[2],c,l)
                H=Hess(tmpY,c,l)
                dH = H[1,1]*H[2,2] - H[1,2]*H[2,1]
                Direction = (matrix(c(H[2,2],-H[2,1],-H[1,2],H[1,1]),2,2))%*%g * (1/dH)
                alpha = 0.1
                i=0
                current=tildeF(X,c,l)
                while(alpha^i>1e-10){
                    tmpY=c(X - (alpha^i)*Direction)
                    if(all(tmpY>1e-10) && all(tmpY<1-1e-10)){
                        if(tildeF(tmpY,c,l)>current){
                            Phi[c,l]=tmpY[1]
                            Psi[c,l]=tmpY[2]
                            break
                        }
                    }
                    i=i+1
                }
            }

        }
        return(list(Phi=Phi,Psi=Psi))
    } # linear in L and C

    ##  Iterate till convergence:
    FEcur=Free_Energy(params)
    FEold=FEcur-1
    vars$Neg_FE=rep(0,(maxiter+1)*8)
    if(params$M>1){stepsize=8}else{stepsize=7}
    #NegFEfile = paste0(outDIR,'/Neg_FE.txt')
    #write(FEcur,NegFEfile)
    vars$Neg_FE[1]=FEcur

    trackingvec=rep(10*track,track)


    while(iteration<=maxiter & continue){
        print(paste0('Iteration ',iteration))
        vars$Beta=updateBeta(params)
        print('Beta')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            vars$Neg_FE[(iteration-1)*stepsize + 2]=FEcur
            #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
            if(FEcur<FEold){
                vars$Error=c(1)
                #write('Error',NegFEfile,append=TRUE)
                break
            }
        }
        vars$WS=updateWS(params)
        print('WS')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
            vars$Neg_FE[(iteration-1)*stepsize + 3]=FEcur
            if(FEcur<FEold){
                vars$Error=c(1)
                #write('Error',NegFEfile,append=TRUE)
                break
            }
        }
        vars$Rho=updateRho(params)
        print('Rho')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
            vars$Neg_FE[(iteration-1)*stepsize + 4]=FEcur
            if(FEcur<FEold){
                vars$Error=c(1)
                #write('Error',NegFEfile,append=TRUE)
                break
            }
        }
        # vars$PhPs=updatePhiPsi(params) # NEED TO FIX THIS... tildeF throws NaNs when too many Inf values returned by log - not sure why this is happening...
        tmp = updatePhiPsi(params)
        vars$Ph = tmp$Phi
        vars$Ps = tmp$Psi
        print('PhiPsi')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
            vars$Neg_FE[(iteration-1)*stepsize + 5]=FEcur
            if(FEcur<FEold){
                vars$Error=c(1)
                #write('Error',NegFEfile,append=TRUE)
                break
            }
        }
        vars$A=updateA(params)
        print('A')
        if(debugging){
            FEold=FEcur
            FEcur=Free_Energy(params)
            #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
            vars$Neg_FE[(iteration-1)*stepsize + 6]=FEcur
            if(FEcur<FEold){
                vars$Error=c(1)
                #write('Error',NegFEfile,append=TRUE)
                break
            }
        }
        if(params$T>1){
            vars$B=updateB(params)
            print('B')
            if(debugging){
                FEold=FEcur
                FEcur=Free_Energy(params)
                #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
                vars$Neg_FE[(iteration-1)*stepsize + 7]=FEcur
                if(FEcur<FEold){
                    vars$Error=c(1)
                    #write('Error',NegFEfile,append=TRUE)
                    break
                }
            }
        }
        if(params$M>1){
            vars$D=updateD(params)
            print('D')
            if(debugging){
                FEold=FEcur
                FEcur=Free_Energy(params)
                #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
                vars$Neg_FE[iteration*stepsize]=FEcur
                if(FEcur<FEold){
                    vars$Error=c(1)
                    #write('Error',NegFEfile,append=TRUE)
                    break
                }
            }
        }
        vars$Lam=updateLam(params)
        print('Lam')
        FEold=FEcur
        FEcur=Free_Energy(params)
        #write(format(FEcur,nsmall=10),NegFEfile,append=TRUE)
        vars$Neg_FE[iteration*stepsize + 1]=FEcur
        if(FEcur<FEold){
            vars$Error=c(1)
            #write('Error',NegFEfile,append=TRUE)
            break
        }
        ##evaluate whether to stop...
        if(stopping){
            if(iteration>1){
                indexingvar=iteration%%track
                if(indexingvar==0){indexingvar=track}
                trackingvec[indexingvar]=sum(abs(round(vars$WS$gamma) - Sold))
                if(mean(trackingvec)<1){continue=FALSE}
                #     if(sum(abs(round(vars$WS$gamma) - Sold))<1){
                #        continue=FALSE
                #    }
            }
        }
        Sold=round(vars$WS$gamma)
        iteration=iteration+1
    }
    vars$maximumiteration=iteration-1
    return(vars)
    # for(i in 1:length(names(vars))){
    #     saveRDS(vars[[i]],paste0(outDIR,'/',names(vars)[i],'.RDS'))
    # }
    # saveRDS(params,paste0(outDIR,'/params.RDS'))
}### end of main function
