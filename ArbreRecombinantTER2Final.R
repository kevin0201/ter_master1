#Arbre recombinant pour la méthode de Hull and White 

####Valeurs du modèle

a <- 0.05
deltat <- 1
volatilité2 <- 0.05
T <- 5 
P <-c(0.03824,0.04512,0.05086,0.05,0.6)
M <- exp(-a*deltat)-1

#Probabilité de transition pour un branchage 1


pu1 <- function(j)  { 1/6 + ((j*M)^2 + j*M)/2 }
pm1<-  function(j) { 2/3 - (j*M)^2 }
pd1 <- function(j)  {1/6 + ((j*M)^2 - j*M)/2 }


#proba de transition pour le Branchage 2 

puu2 <- function(j) { 1/6 + ((j*M)^2 - j*M)/2}
pu2 <- function(j) { -1/3 - (j*M)^2 + 2*j*M }
pm2 <- function(j) {7/6 + ((j*M)^2 - 3*j*M)/2 }

#proba de transition pour le Branchage 3

pm3<- function(j){7/6 + ((j*M)^2 + 3*j*M)/2}
pd3<- function(j){ -1/3 - (j*M)^2 - 2*j*M }
pdd3 <- function(j) { 1/6 + ((j*M)^2 + j*M)/2}



vardeltaR1 <- volatilité2 * (1  - exp(-2*a*deltat))/2*a
#deltaR1bis <- (3*(vardeltaR1))^1/2 
deltaR <- volatilité2 * (3*deltat)^0.5
  
R1 <- function(j){j*deltaR} 
t <- i*deltat 
#eltaR1 <- (-a*R1)*deltat + ecarttype*deltaz

jmax= floor(0.184/a*deltat) + 1.0
IndiceJ=c()
IndiceI=c()
IndiceI<- 1:T
IndiceJ <- jmax:-jmax 

###étape 1 de la méthode de HULL AND WHITE 

arbreEtape1 <- function(){
  Tab=matrix(nrow=2*jmax + 1,ncol= T,dimnames=list(IndiceJ,IndiceI))
  
for(i in 1:T){Tab[jmax + 1,i] = R1(0)}
  for (i in 1:T)
    {
    for(j in 1:jmax)
        {
      
      Tab[j,i] = R1((jmax + 1) - j)
      if(jmax+1-i>=1)
      {
        k <- 1:(jmax+1-i)
        Tab[k,i] = 0
      }
      
      
        }
    }
T2 <- (jmax*2)+1
T1 <- jmax + 2
  
     for(i in 1 :T)
     { 
       for(j in T1:T2)
        {
     Tab[ j ,i] = R1(-(j-(jmax+1)))
      if(jmax+1+i<=T2)
     {
       k <- (jmax+1+i):T2
        Tab[k,i] = 0
     }
       }
   }
    
		return(Tab)
}
arbreEtape1()

###Proba de passage entre les noeuds de deux périodes qui se suivent 

ProbaTransition <- function(l,j){ 
  
  if (-jmax < l && l<jmax)
  {
    if (l -1 == j)
    {
      proba = pd1(l)
      return(proba)
    }
    if (l == j) 
    {
      proba = pm1(l)
      return(proba)
    }
    if (l + 1 == j) 
    {
      proba =  pu1 (l)
      return(proba)
    }
    return(0)
  }
  if (l == -jmax)
  {
    if (l == j )
    {
      proba = pm2(l)
      return(proba)
    }
    if (l+1 == j) 
    {
      proba = pu2(l)
      return(proba)
    }
    if (l + 2 == j) 
    { 
      proba =  puu2 (l)
      return(proba)
    }
    return(0)
  }
  if (l == jmax)
  {
    if (l == j )
    {
      proba = pm3(l)
      return(proba)
    }
    if (l-1 == j) 
    {
      proba = pd3(l)
      return(proba)
    }
    if (l - 2 == j) 
    { 
      proba =  pdd3 (l)
      return(proba)
    }
    return(0)
  }
}


tot <- 2*jmax + 1 
probatran=matrix(nrow = tot,ncol = tot,dimnames = list(IndiceJ, IndiceJ))
for(l in 1:tot){ for (j in 1:tot ){probatran[l,j] <- ProbaTransition(jmax+ 1 -l ,jmax + 1 - j)}}
probatran

 ### Etape 2 de la méthdode de HULL AND WHITE 

 #alpha= matrix(0,nrow = T,ncol=1)
 alpha=c(0.03824,0,0)
 vecteurnoeud=-jmax:jmax
 #CalculerLesAlpha <-function()
 #{
 tailleq=T
 Q=matrix(0,nrow = tailleq,ncol=1)
 for(i in 2:T)
 {
   if (i==2)
   {
     #calcul de alpha 2, Q, puis P, puis Alfpa
     Q[1]=probatran[3,2]*exp(-(alpha[1]*deltat)) #proba up a partir du oeud B
     Q[2]=probatran[3,3]*exp(-(alpha[1]*deltat))
     Q[3]=probatran[3,4]*exp(-(alpha[1]*deltat))
     Qformule=matrix(0,nrow = tailleq,ncol=1)
     Qformule=c(Q[1]*exp(-deltaR*deltat),Q[2],Q[3]*exp(deltaR*deltat))
     alpha[2]=log(sum(Qformule))-log(exp(-P[i]*i))
   }
   #On determine la taille de Q A chaque etape de l'iteration

   if(i>=T)
   {
     tailleq=jmax*2+1
     tailleN=tailleq
   }
   if (tailleq==jmax*2+1)
   {
     #On initialise Q et on calcul de vecteur q suivant qm+1,j
     Qformule2=matrix(0,nrow = tailleq,ncol=1)
     QQ=matrix(0,nrow = tailleq,ncol=1)
     Qformule2=matrix(0,nrow = tailleq,ncol=1)
     for (qi in -jmax:jmax) 
     {
       for (qii in 1:length(Q)) #taille de q precedente
       {
         QQ[qi+jmax+1]=QQ[qi+jmax+1]+Q[qii]*probatran[qii,qi+jmax+1]*exp(-(alpha[i-1]+(qii-length(Q)+1)*deltaR)*deltat) #calcul de Qm+1,j
       }
     }
     alpha[i]=log(sum(QQ*(exp(-vecteurnoeud*deltaR*deltat))))-log(exp(-P[i]*i))
     Q=QQ
   }#fin if2
   
 }
 #  return(alpha)
 #}
 #CalculerLesAlpha()
 alpha
 
###Arbre final 

 ###On opère une translation des noeuds obtenus à l'étape 1 en ajoutant les alphas 
 Arbrefinal <- function(){
   R=matrix(nrow = 2*jmax + 1, ncol = T)
   for (i  in 1:T) {
     for(j in 1:tot){
       R[j,i]<- arbreEtape1()[j,i] + alpha[i]
       if(jmax+1-i>=1)
       {
         k <- 1:(jmax+1-i)
         R[k,i] = 0
       }
       if(jmax+1+i<=2*jmax+1)
       {
         h <- (jmax+1+i):(2*jmax+1)
         R[h,i] = 0
         
       }
     }
     
   }
   
   return(R)
 }
 Arbrefinal()

 