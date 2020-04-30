library(cluster)
install.packages("EnvStats")
library(EnvStats)
library(car)
#!change working directory!
#Data Inladen


wijn <- read.csv("C:/Users/Gebruiker/Desktop/wijn.csv", header=FALSE, sep = ",")
##schaling
for (i in 2:13){
  boxplot(wijn[,i]~wijn$V1)
}

###clustering ongeschaald
chemwijn=wijn[,2:14]
pclst = pam(chemwijn,2)
print(pclst$silinfo)
table(wijn[,1],pclst$clustering)

pclst = pam(chemwijn,4)
table(wijn[,1],pclst$clustering)
##clustering geschaald
chemwijnsc=scale(chemwijn)
for (i in 1:7){
  print(pam(chemwijnsc,i)$silinfo$avg.width)
}
pclstsc=pam(chemwijnsc,3)
table(wijn$V1,pclstsc$clustering)

##fuzzy clustering
fclstsc=fanny(chemwijnsc,3,memb.exp=1.5)
fclstsc$silinfo$avg.width
table(fclstsc$clustering,wijn$V1)

##hyrarchical clustering
aclstsc=agnes(chemwijnsc)
aclstsc
pltree(aclstsc)
dclstsc=diana(chemwijnsc)
pltree(dclstsc,main="dendrogram van hyrarchishe clustering",xlab="wijnen",ylab="hoogte",sub="")
table(wijn$V1,cutree(dclstsc,3))

###noraliteit
for (i in 2:14){
  hist(wijn[,i])
}
for (i in 2:14){
  qqnorm(wijn[,i])
}
isnorm = rep(FALSE,13)
for (i in 2:14){
  isnorm[i-1]=shapiro.test(wijn[,i])$p.value>0.01
}
islognorm = rep(FALSE,13)
for (i in 2:14){
  islognorm[i-1]=shapiro.test(log(wijn[,i]))$p.value>0.01
}
lambda=rep(0,13)
isbcnorm = rep(FALSE,13)
for (i in 2:14){
  lambda[i-1]=boxcox(wijn[,i],optimize=TRUE)$lambda
  lmda=lambda[i-1]
  isbcnorm[i-1]=0.01<shapiro.test((wijn[,i]^lmda-1)/lmda)$p.value
}
print(isnorm)
print(islognorm)
print(isbcnorm)
print(lambda)
#we transformeren de data afhankelijk van voorgaande testresultaten
wijntrans=data.frame(wijn$V1)
for (i in 2:14){
  if (isnorm[i-1]) {
    wijntrans[,i]=wijn[,i]
  }
  else if(islognorm[i-1]){
    wijntrans[,i]=log(wijn[,i])
  }
  else{
    wijntrans[,i]=((wijn[,i])^lambda[i-1]-1)/lambda[i-1]
  }
}
implement=c(T)
implement=append(implement,isbcnorm)
scatterplotMatrix(wijntrans[,implement],diagonal="boxplot",smooth=F,regLine=F)
proteïne = wijn[,13]
densityPlot(proteïnes[soort],main="dichtheid van proteinegehalte",xlab="hoeveelheid proteïne",ylab="benaderde dichtheid")

##normaliteit van verdeelde data
soort=array(c(F),dim=c(178,3))
for (j in 1:3){
  soort[,j]=wijn$V1==j
}
wijn1=wijn[soort[,1],]
wijn2=wijn[soort[,2],]
wijn3=wijn[soort[,3],]
isnorm=array(c(FALSE),dim = c(3,13))
for (j in 1:3){
  for (i in 2:14){
    isnorm[j,i-1]=0.01<shapiro.test(wijn[soort[,j],i])$p.value
  }
}
print(isnorm)  
scatterplotMatrix(wijn3,diagonal="boxplot",smooth=F,regLine=F)
#verdeelde druiven op getransformeerde data
for (j in 1:3){
  for (i in 2:14){
    isnorm[j,i-1]=0.01<shapiro.test(wijntrans[soort[,j],i])$p.value
  }
}
print(isnorm)
print(isbcnorm)

#genereren van voorbeeld plot
plot(density(wijntrans[,7]),main="Dichtheden van fenolen",xlab="fenol gehalte",ylab="dichtheid")
temp=density(wijntrans[soort[,1],7])
lines(temp$x,temp$y*length(wijn[soort[,1],1])/length(wijn[,1]),col='red')
for (j in 2:14){
for (i in 1:3){
plot(density(wijntrans[soort[,i],j]))
}}
#in een tweede poging nemen we een gewogen gemiddelde van de lambda's means = array(c(0),dim=c(3,13))
#merk op dat de validiteit van deze methode onzeker is
means = array(c(0),dim=c(3,13))
vars = array(c(0),dim=c(3,13))
for (j in 1:3){
  for (i in 2:14){
    means[j,i-1]=mean(wijn[soort[,j],i])
    vars[j,i-1]=var(wijn[soort[,j],i])
  }
}
lambda = array(c(0),dim=c(3,13))
for(j in 1:3){
  for (i in 2:14){
    lambda[j,i-1]=boxcox(wijn[soort[,j],i],lambda=c(-3,3),optimize=T)$lambda
  }
}
lambda
weights=array(c(0),dim=c(3,13))
#als gewichten word er rekening gehouden met de verandering van de kromming voor de gegeven data, hoe meer de kromming veranderd,
#en hoe groter de sprijding van de variabele daar, hoe resistenter
weights=abs(means^(lambda-2)*(1+(lambda-1)*log(means)))*vars
weights
for (i in 1:13){
  weights[,i]=weights[,i]/sum(weights[,i])
}
weights
newlambda=array(c(0),dim=c(1,13))
for (i in 1:13){
  newlambda[,i]=sum(weights[,i]*lambda[,i])
}
newlambda  
isbcnorm=array(c(F),dim=c(3,13))
for (i in 2:14){
  for (j in 1:3){
    isbcnorm[j,i-1]=0.05<shapiro.test(((wijn[soort[,j],i])^newlambda[1,i-1]-1)/newlambda[1,i-1])$p.value
  }
}
isbcnorm
transdata=(wijn[,2:14]^newlambda-1)/newlambda
scatterplotMatrix(transdata,diagonal="boxplot",smooth=F,regLine=F)
for (j in 1:3){
  scatterplotMatrix(transdata[soort[,j],],smooth=F,regLine=F)
}
for (j in 1:3){
  scatterplotMatrix(transdata[soort[,j],useful],smooth=F,regLine=F)
}


  