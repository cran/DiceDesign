
############################################################
# MAIN FUNCTIONS OF THIS FILE:
# discrepSA_LHS,discrepESE_LHS: low-discrepancy LHS
# maximinSA_LHS,maximinESE_LHS: Maximin LHS

# Optimization algorithms:
# SA  = Simulated Annealing 
# ESE = Enhanced Stochastic Evolutionary 

# EP  = Elementary Permutation

#####REFERENCES#####

#Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs:
#optimization algorithms and subprojection properties, Journal of Simulation, submitted.

#D.Morris and J.Mitchell. Exploratory designs for computationnal experiments. Journal of 
#Statistical Planning and Inference, 43:381-402, 1995.

#R.JIN, W.CHEN, A.SUDJIANTO. An efficient algorithm for constructing optimal design
#of computer experiments. Journal of Statistical Planning and Inference, 134:268-287, 2005

#A.MARREL Mise en oeuvre et utilisation du méta-modèle processus gaussien pour
#l'analyse de sensibilité de modèles numériques. PhD thesis, INSA Toulouse, France, 2008.

###############################################################

#####INCLUDING IN THE DICE_DESIGN PACKAGE#####

#require(DiceDesign)

#####discrepancyL2_STAR#####


#---------------------------------------------------------------------------|
#args :  i,j : two points x_i and x_j from the design m                     |      
#        m   : the design                                                   |
#depends :  "discrepancyTERM_L2_STAR"                                       |
#---------------------------------------------------------------------------|
   
discrepancyL2_STAR=function(m) 
{
  n<-nrow(m)
  d<-ncol(m)
  dL2<-0
  for (j in 1:n)
    {for (i in 1:n)
      {dL2<-dL2+discrepancyTERM_L2_STAR(i,j,m)}
    }  
dL2<-sqrt(3^(-d)+dL2)
return(dL2)
}

   

#####discrepancyTERM_L2_STAR#####
#---------------------------------------------------------------------------|
#depends :  both "term_oneDD" and "term_twoDD"                              |
#---------------------------------------------------------------------------|
  

discrepancyTERM_L2_STAR <- function(i,j,m)
{
  n<-nrow(m)
  d<-ncol(m)
  if(i!=j)
      { t<-c()
        for (l in 1:d)
          {t<-c(t,1-max(m[i,l],m[j,l]))}
        t<-(prod(t))/(n^2)
                      

      }
   else 
      {t<-term_oneDD(i,m)/(n^2)-((2^(1-d))/n)*term_twoDD(i,m)}
return(t)
}

#####term_twoDD#####
   
term_twoDD=function(i,m)
    
{
 t2<-1-m[i,]^2
 t2<-prod(t2)
 return(t2)
}

#####term_oneDD#####

term_oneDD=function(i,m)
      
{
 t1<-1-m[i,]
 t1<-prod(t1)
 return(t1)
}


#####FUNCTION PERFORMING ELEMENTARY PERMUTATION (EP) IN LHD##
#####USED IN SA ALGORITHMS#####


#####lhs_EP#####
#---------------------------------------------------------------------------|
#args :  m = the design                                                     |
#out  :  l = list including design after EP, ligns and columns defining EP  |          
#---------------------------------------------------------------------------|


lhs_EP<-function(m) 
{  
 G<-m
 if (!is.matrix(G)) G<-m[[1]]
 d<-ncol(G)
 n<-nrow(G)
 ligns<-trunc(runif(2,0,n))+1      
 column<-trunc(runif(1,0,d))+1
 x<-G[ligns[1],column]
 G[ligns[1],column]<-G[ligns[2],column]
 G[ligns[2],column]<-x 
 l<-list(G,ligns[1],ligns[2],column)   
 return(l)
} 

#####phiP_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the power in phiP criterion                                |
#out     l     : a list with new design and its mindist value               |
#---------------------------------------------------------------------------|

phiP_EP_ESE<-function(m,k,p) 
{  
 d<-ncol(m)
 n<-nrow(m)
 G<-m
 ligns<-trunc(runif(2,0,n))+1      
 x<-G[ligns[1],k]
 G[ligns[1],k]<-G[ligns[2],k]
 G[ligns[2],k]<-x 
 l<-list(phiP(G,p),G)   
 return(l)
} 



#####INTERMEDIATE FUNCTIONS#####

alpha_oneDD<-function(i1,i2,k,m)
{

alpha<-(1-m[i2,k])/(1-m[i1,k])
alpha
}

beta_oneDD<-function(i1,i2,k,m)
{

beta<-(1-m[i2,k]^2)/(1-m[i1,k]^2)
beta
}


gamma_oneDD<-function(i1,i2,k,j,m)
{

gamma<-(1-max(m[i2,k],m[j,k]))/(1-max(m[i1,k],m[j,k]))
gamma
}

#####DESIGN'S L2 STAR DISCREPANCY VALUE AFTER AN EP#####

#####discrepancyL2_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the L2_star_discrepancy of m                 |
#out     l : list with dL2_2 (the square of the L2_star_discrepancy         |
#             of the EP design and the matrix corresponding to the EP design|                    
#depends :  term_oneDD, term_twoDD, alpha_oneDD, beta_oneDD,                |
#           gamma_oneDD, discrepancyTERM_L2_STAR                            |
#---------------------------------------------------------------------------|
   

discrepancyL2_EP_ESE<-function(m,k,p)
{
   G<-m
   i<-trunc(runif(2,0,nrow(m)))+1   # On genère un autre LHS à partir de M, en le perturbant de façon minimale
   x<-G[i[1],k]
   G[i[1],k]<-G[i[2],k]
   G[i[2],k]<-x 
   i1<-i[1]
   i2<-i[2]


   n<-nrow(m)
   d<-ncol(m)
   dL2_2<-p+(term_oneDD(i1,m)*alpha_oneDD(i1,i2,k,m))/(n^2)-beta_oneDD(i1,i2,k,m)*term_twoDD(i1,m)*((2^(1-d))/n)-discrepancyTERM_L2_STAR(i1,i1,m)
   dL2_2<-dL2_2+(term_oneDD(i2,m)/(alpha_oneDD(i1,i2,k,m)*(n^2)))-((term_twoDD(i2,m)*((2^(1-d))/n))/(beta_oneDD(i1,i2,k,m)))-discrepancyTERM_L2_STAR(i2,i2,m)
   y<-rep(0,n-2)
   w<-0
   for (j in 1:n)
    {if (j!=i1 & j!=i2)
      {
     w<-w+1
     y[w]<-gamma_oneDD(i1,i2,k,j,m)*discrepancyTERM_L2_STAR(i1,j,m)-discrepancyTERM_L2_STAR(i1,j,m)+(discrepancyTERM_L2_STAR(i2,j,m)/(gamma_oneDD(i1,i2,k,j,m)))-discrepancyTERM_L2_STAR(i2,j,m)
      }

   }

     y<-sum(y)
     dL2_2<-dL2_2+2*y
     l<-list(dL2_2,G)
     return(l)
}


#####DESIGN'S L2 DISCREPANCY VALUE AFTER AN EP#####


#####discrepancyL2_EP#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        i1,i2 : the two ligns of the lhs elementary permutation            |              
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the L2_star_discrepancy of m                 |
#out     dL2_2 : the square of the L2_star_discrepancy of the EP design     |
#depends :  term_oneDD, term_twoDD, alpha_oneDD, beta_oneDD,                |
#           gamma_oneDD, discrepancyTERM_L2_STAR                            |
#---------------------------------------------------------------------------|
   
discrepancyL2_EP<-function(m,i1,i2,k,p)
{
   n<-nrow(m)
   d<-ncol(m)
   dL2_2<-p+(term_oneDD(i1,m)*alpha_oneDD(i1,i2,k,m))/(n^2)-beta_oneDD(i1,i2,k,m)*term_twoDD(i1,m)*((2^(1-d))/n)-discrepancyTERM_L2_STAR(i1,i1,m)
   dL2_2<-dL2_2+(term_oneDD(i2,m)/(alpha_oneDD(i1,i2,k,m)*(n^2)))-((term_twoDD(i2,m)*((2^(1-d))/n))/(beta_oneDD(i1,i2,k,m)))-discrepancyTERM_L2_STAR(i2,i2,m)
   y<-rep(0,n-2)
   w<-0
   for (j in 1:n)
    {if (j!=i1 & j!=i2)
       {
         w<-w+1
         y[w]<-gamma_oneDD(i1,i2,k,j,m)*discrepancyTERM_L2_STAR(i1,j,m)-discrepancyTERM_L2_STAR(i1,j,m)+(discrepancyTERM_L2_STAR(i2,j,m)/(gamma_oneDD(i1,i2,k,j,m)))-discrepancyTERM_L2_STAR(i2,j,m)
      }

   }

     y<-sum(y)
     dL2_2<-dL2_2+2*y

return(dL2_2)
}

#####gammaWDD#####

gammaWDD=function(i1,i2,k,j,m)
{
gamma<-((3/2)-(abs(m[i2,k]-m[j,k])*(1-abs(m[i2,k]-m[j,k]))))/((3/2)-(abs(m[i1,k]-m[j,k])*(1-abs(m[i1,k]-m[j,k]))))
gamma
}

#####ccWDD#####

ccWDD=function(i,j,m)
{
 n<-nrow(m)

 if(i!=j)
 {c<-(3/2)-((abs(m[i,]-m[j,]))*(1-abs(m[i,]-m[j,])))
 c<-(prod(c))/(n^2)}
c
}


#####discrepancyW2_EP#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        i1,i2 : the two ligns of the lhs elementary permutation            |              
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the W2_star_discrepancy of m                 |
#out     dW2_2 : the square of the W2_star_discrepancy of the new design    |
#depends       : ccWDD, gammaWDD                                            |
#---------------------------------------------------------------------------|
   


discrepancyW2_EP=function(m,i1,i2,k,p)        
#m=le LHS initial, i1 et i2=les deux lignes de la permutation, k la colonne,  DW=discrépance Wrap-around de m au carré
{

 n<-nrow(m)
 dW2_2<-rep(0,n-2)
 w<-0
 for (j in 1:n)
    {if (j!=i1 & j!=i2)
         
          {w<-w+1
           dW2_2[w]<-gammaWDD(i1,i2,k,j,m)*ccWDD(i1,j,m)-ccWDD(i1,j,m)+ccWDD(i2,j,m)/(gammaWDD(i1,i2,k,j,m))-ccWDD(i2,j,m)}

    }

dW2_2<-sum(dW2_2)
dW2_2<-p+2*dW2_2
return(dW2_2)

}

#####discrepancyW2_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the W2_star_discrepancy of m                 |
#out     l     : a list with new design and it W2 discrepancy value         |
#depends       : ccWDD, gammaWDD                                            |
#---------------------------------------------------------------------------|
   


discrepancyW2_EP_ESE=function(m,k,p)    
{
 n<-nrow(m)
 G<-m
 i<-trunc(runif(2,0,n))+1  
 i1=i[1]
 i2=i[2]
 x<-G[i1,k]
 G[i1,k]<-G[i2,k]
 G[i2,k]<-x 


 
 dW2<-rep(0,n-2)
 w<-0
 for (j in 1:n)
    {if (j!=i1 & j!=i2)
         
          {w<-w+1
           dW2[w]<-gammaWDD(i1,i2,k,j,m)*ccWDD(i1,j,m)-ccWDD(i1,j,m)+ccWDD(i2,j,m)/(gammaWDD(i1,i2,k,j,m))-ccWDD(i2,j,m)}

    }

 dW2<-sum(dW2)
 dW2<-sqrt(p+2*dW2)
 l<-list(dW2,G)
 return(l)

}



#####alphaDD#####

alphaDD<-function(i1,i2,k,m)
{

alpha<-(1+abs(m[i2,k]))/(1+abs(m[i1,k]))
alpha
}

#####betaDD#####

betaDD<-function(i1,i2,k,m)
{

beta<-(2-abs(m[i2,k]))/(2-abs(m[i1,k]))
beta
}

#####gammaDD#####

gammaDD<-function(i1,i2,k,j,m)
{

gamma<-(2+abs(m[i2,k])+abs(m[j,k])-abs(m[i2,k]-m[j,k]))/(2+abs(m[i1,k])+abs(m[j,k])-abs(m[i1,k]-m[j,k]))
gamma
}

#####gDD#####

gDD<-function(i,m)
{
 g<-1+abs(m[i,])
 g<-prod(g)
 g
}

#####hDD#####

hDD<-function(i,m)
{
 h<-1+(abs(m[i,])/2)-0.5*(m[i,]^2)
 h<-prod(h)
 h
}

#####ccDD#####

ccDD<-function(i,j,m)
{
 n<-nrow(m)

 if(i!=j)
 {c<-0.5*(2+abs(m[i,])+abs(m[j,])-abs(m[i,]-m[j,]))
 c<-(prod(c))/(n^2)}
else {
c<-gDD(i,m)/(n^2)-2*hDD(i,m)/n}
c
}


#####discrepancyC2_EP#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        i1,i2 : the two ligns of the lhs elementary permutation            |              
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the C2_star_discrepancy of m                 |
#out     dC2_2 : the square of the C2_star_discrepancy of the new design    |
#depends       : alphaDD, betaDD, gammaDD, gDD, hDD, ccDD                   |
#---------------------------------------------------------------------------|
   
discrepancyC2_EP=function(m,i1,i2,k,p)        #m=le LHS initial,     i1 et i2=les deux lignes de la permutation, k la colonne,  DW=discrépance Wrap-around de m au carré
{

 n<-nrow(m)
 m<-m-0.5
 dC2_2<-p+(gDD(i1,m)*alphaDD(i1,i2,k,m))/(n^2)-2*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)*hDD(i1,m)/n
 dC2_2<-dC2_2-ccDD(i1,i1,m)-ccDD(i2,i2,m) 
 dC2_2<-dC2_2+(gDD(i2,m)/((n^2)*alphaDD(i1,i2,k,m)))-(2*hDD(i2,m)/(n*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)))


 x<-rep(0,n-2)
 w<-0
  for (j in 1:n)
    {if (j!=i1 & j!=i2)
     {
      w<-w+1
      x[w]<-gammaDD(i1,i2,k,j,m)*ccDD(i1,j,m)-ccDD(i1,j,m)+ccDD(i2,j,m)/(gammaDD(i1,i2,k,j,m))-ccDD(i2,j,m)}
     } 

x<-sum(x)
dC2_2<-dC2_2+2*x
return(dC2_2)


}


#####discrepancyC2_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the C2_star_discrepancy of m                 |
#out     l     : a list with new design and it C2 discrepancy value         |
#depends       : alphaDD, betaDD, gammaDD, gDD, hDD, ccDD                   |
#---------------------------------------------------------------------------|

discrepancyC2_EP_ESE=function(m,k,p)        #m=le LHS initial,     i1 et i2=les deux lignes de la permutation, k la colonne,  DW=discrépance Wrap-around de m au carré
{
 n<-nrow(m)
 G<-m
 i<-trunc(runif(2,0,n))+1 
 i1=i[1]
 i2=i[2] 
 x<-G[i1,k]
 G[i1,k]<-G[i2,k]
 G[i2,k]<-x 
 
 m<-m-0.5
 dC2<-p+(gDD(i1,m)*alphaDD(i1,i2,k,m))/(n^2)-2*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)*hDD(i1,m)/n
 dC2<-dC2-ccDD(i1,i1,m)-ccDD(i2,i2,m) 
 dC2<-dC2+(gDD(i2,m)/((n^2)*alphaDD(i1,i2,k,m)))-(2*hDD(i2,m)/(n*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)))


 x<-rep(0,n-2)
 w<-0
  for (j in 1:n)
    {if (j!=i1 & j!=i2)
     {
      w<-w+1
      x[w]<-gammaDD(i1,i2,k,j,m)*ccDD(i1,j,m)-ccDD(i1,j,m)+ccDD(i2,j,m)/(gammaDD(i1,i2,k,j,m))-ccDD(i2,j,m)}
     } 

x<-sum(x)
dC2<-sqrt(dC2+2*x)
l=list(dC2,G)
return(l)


}


#####discrepSA_LHS#####
#####L2 DISCREPANCY LHS VIA SIMULATED ANNEALING OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  m     : the design                                                 |
#        T0    : the initial temperature                                    |
#        c     : paramater regulating the temperature                       |
#        it    : the number of iterations                                   |
#    criterion : criterion to be optimized ("DL2","DW2","DC2")              |
#       profile: "GEOM" or "LINEAR_MORRIS". By default : "GEOM"             |
#       Imax   :paramater adjusting number of iterations without improvement|
#               (if you choose "LINEAR_MORRIS" profile)                     |
#out           : low L2_STAR design                                         |
#depends :  discrepancyL2_STAR, discrepancyW2_EP, lhs_EP, discrepancyL2_EP  |
#           discrepancyC2_EP, lhs_EP                                        |
#---------------------------------------------------------------------------|
   

discrepSA_LHS<-function(design,T0=10,c=0.95,it=2000,criterion="DC2",profile="GEOM",Imax=100) 

{ 
  m<-design
  if(criterion=="DC2")
  
    {if(profile=="GEOM"){
         
           i<-0
           T<-T0
           v<-discrepancyCriteria(m)[[1]]	
           while(T>0 & i<it){
           
                   G<-lhs_EP(m)
                   g<-discrepancyC2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
                   g<-sqrt(g)
                 
                   diff<-exp((v-g)/T)
                         if (diff>1){m<-G[[1]]
                                     v<-g    
                                    }
                       else         {
                                    Bernoulli<-rbinom(1,1,diff)
                                    if (Bernoulli==1){m=G[[1]]
                                      v<-g    }}
              i<-i+1
              T<-(c^i)*(T0)
                       }
             return(m)
           }

    if(profile=="LINEAR_MORRIS"){

       
       T<-T0
       Dbest<-m
       v<-discrepancyCriteria(m)[[1]]
       ref<-v
       for (i in 1:it)
         {
          flag<-0
          I<-1
          while(I<Imax){
                     
             G<-lhs_EP(m)
             g<-discrepancyC2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
             g=sqrt(g)
             diff<-exp((v-g)/T)
             if (diff>1){m<-G[[1]]
                         v<-g   
                         flag=1
                       
                       if(v<ref){ref=v
                                  Dbest=m
                                  I=1}
                          else {I=I+1}}
                        
             else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G[[1]]
                                    v=g 
                                    flag=1  
                         if (v<ref){ref=v
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
         } 
          
       if (flag==1) {T=T*c}
               else {break}
        } 
        return(Dbest)}
   } # end if(design=="DC2")


   
   if(criterion=="DL2")
    {if(profile=="GEOM")
          {i<-0
           T<-T0
           v<-discrepancyL2_STAR(m)	
           while(T>0 & i<it){
           
                    G<-lhs_EP(m)
                    g<-discrepancyL2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
                    g<-sqrt(g)
                 
                    diff<-exp((v-g)/T)
                        if (diff>1){m<-G[[1]]
                                  v<-g    
                                 }
                    else         {
                                 Bernoulli<-rbinom(1,1,diff)
                                 if (Bernoulli==1){m=G[[1]]
                                    v<-g    }}
                  i<-i+1
                  T<-(c^i)*(T0)
                    }
               return(m)
             }
      if(profile=="LINEAR_MORRIS")

       {
        T<-T0
        Dbest<-m
        v<-discrepancyL2_STAR(m)
        ref<-v
        for (i in 1:it)
         {
          flag<-0
          I<-1
          while(I<Imax){
                     
             G<-lhs_EP(m)
             g<-discrepancyL2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
             g=sqrt(g)
             diff<-exp((v-g)/T)
             if (diff>1){m<-G[[1]]
                         v<-g   
                         flag=1
                       
                       if(v<ref){ref=v
                                  Dbest=m
                                  I=1}
                          else {I=I+1}}
                        
             else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G[[1]]
                                    v=g 
                                    flag=1  
                         if (v<ref){ref=v
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
         } 
          
       if (flag==1) {T=T*c}
               else {break}
        } 
        return(Dbest)}
   } # end if(design=="DL2")
 
       
  
   if(criterion=="DW2")
    {if(profile=="GEOM")
          {i<-0
           T<-T0
           v<-discrepancyCriteria(m)[[5]]	
           while(T>0 & i<it){
           
                G<-lhs_EP(m)
                g<-discrepancyW2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
                g<-sqrt(g)
                 
                diff<-exp((v-g)/T)
                      if (diff>1){m<-G[[1]]
                                  v<-g    
                                 }
                    else         {
                                 Bernoulli<-rbinom(1,1,diff)
                                 if (Bernoulli==1){m=G[[1]]
                                    v<-g    }}
                     i<-i+1
                     T<-(c^i)*(T0)
                    }
            return(m)
           }
     
       if(profile=="LINEAR_MORRIS")

       {
        T<-T0
        Dbest<-m
        v<-discrepancyCriteria(m)[[5]]
        ref<-v
        for (i in 1:it)
         {
          flag<-0
          I<-1
          while(I<Imax){
                     
             G<-lhs_EP(m)
             g<-discrepancyW2_EP(m,G[[2]],G[[3]],G[[4]],v^2)
             g=sqrt(g)
             diff<-exp((v-g)/T)
             if (diff>1){m<-G[[1]]
                         v<-g   
                         flag=1
                       
                       if(v<ref){ref=v
                                  Dbest=m
                                  I=1}
                          else {I=I+1}}
                        
             else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G[[1]]
                                    v=g 
                                    flag=1  
                         if (v<ref){ref=v
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
         } 
          
       if (flag==1) {T=T*c}
               else {break}
        } 
        return(Dbest)}
     } # end if(design=="DW2")
 
} # end function

#####Function to extract elements from lists composed with lists#####
extract_list<-function(l){return(l[[1]])}



#####discrepESE_LHS#####
#####L2 DISCREPANCY LHS VIA ESE OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  design     : the design                                            |
#        T0    : the initial temperature                                    |
#        inner_it  : number of iterations for inner loop                    |
#        J     : number of new proposed LHS in inner loop                   |
#        it    : number of iterations for outer loop                        |
#        design: "DC2", "DW2" or "DL2"                                      |
#out           : B=new design with low discrepancy                          |
#depends :  discrepancyL2_EP_ESE, discrepancyW2_EP_ESE, discrepancyL2_STAR  |
#           discrepancyC2_EP_ESE, discrepancyCriteria                       |
#---------------------------------------------------------------------------|


discrepESE_LHS<-function(design,T0=0.005*discrepancyCriteria(design)[[1]],inner_it=100,J=50,it=2,criterion="DC2")
{
  m<-design
  if(criterion=="DC2")
  {
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  dC2<-discrepancyCriteria(m)[[1]]                   
  best<-dC2              

  for (q in 1:it)
   {
     
     BOLD<-Best
     bold<-best                                   # Best=new LHS built at every step        
                                                # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) #liste de liste
     
     g<-lapply(l,discrepancyC2_EP_ESE,k=modulo+1,p=discrepancyCriteria(m)[[1]]^2)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]  
                                             
       Delta<-a-dC2
                                                
     if((Delta)<=(Temperature*runif(1)))                 # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                                # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           dC2<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio


  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step

    }}


  if(criterion=="DL2")
  {
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  dL2<-discrepancyL2_STAR(m)                   
  best<-dL2              

  for (q in 1:it)
   {

     BOLD<-Best 
     bold<-best                                   # Best=new LHS built at every step        
                                                # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) #liste de liste
     
     g<-lapply(l,discrepancyL2_EP_ESE,k=modulo+1,p=dL2^2)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]
                          
                                               
       Delta<-a-dL2
                                                
     if((Delta)<=(Temperature*runif(1)))                 # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                                # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           dL2<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio


  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step



    }}

 if(criterion=="DW2")
  {
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  dW2<-discrepancyCriteria(m)[[5]]                  
  best<-dW2              

  for (q in 1:it)
   {

     BOLD<-Best
     bold<-best                                    # Best=new LHS built at every step        
                                                # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) #liste de liste
     
     g<-lapply(l,discrepancyW2_EP_ESE,k=modulo+1,p=discrepancyCriteria(m)[[5]]^2)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]
                          
                                               
       Delta<-a-dW2
                                                
     if((Delta)<=(Temperature*runif(1)))                 # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                               # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           dW2<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio


  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step

    }}

return(Best)

}


#####maximinSA_LHS#####
#####Maximin LHS VIA SIMULATED ANNEALING OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  m     : the design                                                 |
#        T0    : the initial temperature                                    |
#        c     : paramater regulating the temperature                       |
#        it    : the number of iterations                                   |
#        p     : power required in phiP criterion                           |
#      profile : temperature down profile, "GEOM" by default.               |
#out           : maximin design                                             |
#depends :  phiP,lhs_EP                                                     |
#---------------------------------------------------------------------------|
   

maximinSA_LHS<-function(design,T0=10,c=0.95,it=2000,p=50,profile="GEOM",Imax=100) 

{ if(profile=="GEOM"){
  m<-design
  i<-0
  T<-T0
  fi_p<-phiP(m,p)
  
  
  while(T>0 & i<it)

         {           
               G<-lhs_EP(m)[[1]]
               fi_p_ep<-phiP(G,p)
               diff<-exp((fi_p-fi_p_ep)/T)
               if (diff>1){m<-G
                       fi_p<-fi_p_ep    }
               else { Bernoulli<-rbinom(1,1,diff)
                      if (Bernoulli==1){m<-G
                                    fi_p<-fi_p_ep    }}

           i<-i+1
           T<-(c^i)*(T0)
           }
  return(m)    
  }
  if(profile=="LINEAR_MORRIS") 
     
  {m<-design
   T<-T0
   Dbest<-m
   fi_p=phiP(m,p)
   ref<-fi_p
   for (i in 1:it)
   {
     flag=0
     I=1
     while(I<Imax){
                     
           G=lhs_EP(m)[[1]]
           fi_p_ep=phiP(G,p)
           diff=exp((fi_p-fi_p_ep)/T)
           if (diff>1){m=G
                       fi_p=fi_p_ep    
                       flag=1
                       
                       if(fi_p<ref){ref=fi_p
                                  Dbest=m
                                  I=1}
                          else {I=I+1}}
                        
           else {Bernoulli=rbinom(1,1,diff)
                   if (Bernoulli==1){m=G
                                    fi_p=fi_p_ep 
                                    flag=1  
                       if (fi_p<ref){ref=fi_p
                                    Dbest=m  
                                    I=1}
                            else {I=I+1}}
                                    
                   else {I=I+1}     
                 }
         } 
          
       if (flag==1) {T=T*c}
               else {break}
   } 
return(Dbest)
   }
   
}

#####maximinESE_LHS#####
#####Maximin LHS VIA ESE OPTIMIZATION#####

#---------------------------------------------------------------------------|
#args :  design: the design                                                 |
#        T0    : the initial temperature                                    |
#        inner_it  : number of iterations for inner loop                    |
#        J     : number of new proposed LHS in inner loop                   |
#        it    : number of iterations for outer loop                        |
#        p     : the power in phiP criterion                                |
#out           : a design with an improved mindist value                    |
#depends :     phiP_EP_ESE, phiP                                            |
#---------------------------------------------------------------------------|

maximinESE_LHS<-function(design,T0=0.005*phiP(design,p=50),inner_it=100,J=50,it=1,p=50)
{
  m<-design
  d<-ncol(m)
  Temperature<-T0
  Best<-m
  fi_p<-phiP(m,p)                   
  best<-fi_p             

  for (q in 1:it)
   {

     BOLD<-Best
     bold<-best                              # B=new LHS built at every step        
                                             # BOLD= new LHS built at each iteration q
     ni<-0
     count<-0
     na<-0
     while(count<=inner_it)                      # inner loop
     {
       count<-count+1
   

     modulo<-count%%d                         # d : number of columns of m
     l<-list(m)
     l<-rep(l,J) 
     
     g<-lapply(l,phiP_EP_ESE,k=modulo+1,p=p)
     values<-lapply(g,extract_list) 
     k<-which.min(values)
     a<-values[[k]]
                                               
       Delta<-a-fi_p
                                                
     if((Delta)<=(Temperature*runif(1)))        # higher is the temperature, higher is the probability of accepting a bad design.
                                                # if Delta is low, the probability is high of accepting a bad design.   
                                                # if Delta>Temperature, m is never accept.
       
     
          {m<-g[[k]][[2]]
           
           fi_p<-a
                    
          na<-na+1
             if(a<=best)  
                 {Best<-m
                  best<-a
                 ni<-ni+1}                       #if optimization is ok, ni=ni+1
           }
    }
          
  v1<-na/inner_it    # v1<-acceptance ratio
  v2<-ni/inner_it    # v2<-optimization ratio


  if (best-bold<0){f<-1
              if(v1>=0.1 & v2<=v1)
                  {Temperature<-0.8*Temperature}
                       else {if (v1>=0.1 & v2==v1){} 
                                      else {Temperature<-Temperature/0.8}
                            }  
                               }

                                               # if the criteria is optimized after the inner loop, then f equals 1
      else {f<-0
      if (v1<=0.1){Temperature<-Temperature/0.7}
              else {Temperature<-Temperature*0.9}
     }

                                               # else, it is the exploratory step
}
return(Best)
}
