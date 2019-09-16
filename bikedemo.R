set.seed(1)
library(MASS)
library(mgcv)
library(refund) 
library(fda)
#data<-read.csv("hour.csv")
#data<-data[data$weekday==6,]
#id<-factor(data$dteday)
#levels(id)<-c(1:105)
#bikedatanp<-data.frame(id=id,hour=(data$hr/23),temp=data$temp,atemp=data$atemp,hum=data$hum,windspeed=data$windspeed,casual=data$casual)
#write.csv(bikedatanp,file="bikenp.csv",row.names = FALSE)

data<-read.csv("bikenp.csv")
attach(data)
n=length(unique(id))
m=length(unique(hour))
T<-list()
for(i in 1:n)
{ind<-which(id==i)
T[[i]]<-hour[ind]
}
UT<-c()
for( i in 1:n)
{UT<-union(UT,T[[i]])}
UT<-sort(UT)
Y<-list()
for(i in 1:n)
{ind<-which(id==i)
Y[[i]]<-casual[ind]
}
GrandY<-NULL
for(i in 1:n)
{GrandY<-c(GrandY,Y[[i]])
}
####log transform#########
GrandY<-log(1+GrandY)
y<-GrandY
#######################
mydata<-data[,c(1,2,3,4,5,6)]
names(mydata)[2]<-c("time")
#######preprocess############
source('NPFCMselection.R')
mydata<- pprocess(mydata)
################add random covariates#############
p=16
A<-matrix(0,n,p)
B<-matrix(0,n,p)
for(i in 1:n)
{for(j in 1:p)
{
  A[i,j]<-rnorm(1,0,2) #mean 50 in parametric
}
}
for(i in 1:n)
{for(j in 1:p)
{
  B[i,j]<-rnorm(1,0,2) 
}
}
X<-function(i,j,t){A[i,j]*sqrt(2)*sin(2*pi*j*t)+B[i,j]*sqrt(2)*cos(2*pi*j*t)}

for(i in 1:2512)
  for(j in 7:22)
  {
    {mydata[i,j]<- X(mydata$id[i],(j-6),mydata$time[i])}}

names(mydata)[7:22]<-paste("pseudo",1:16)
names(mydata)[-c(1:2)]
#####################################################################

### Perform variable selection in NPFCM#################

NPFCM.select(y,mydata,basist=c(4:6),basisx=c(4:6))