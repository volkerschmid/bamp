# bamp.R v1.3.0.1 (08 aug 2006)

# edit first age of first age group
startingage <- 45
# edit first period of first period group
startingperiod <- 1945
# edit years per period group
yearsperperiod <- 1
# edit name of your .ini file
inifile <- "bamp.ini"
# edit root folder
bamp.dir <- "c://bamp//"
# set psfile to 0 if you don't want a postscript file
psfile <- 0

# O.K. 
# Unix users: save and start with "Splus < bamp.ssc"
# Windows users: select "select all" in "Edit" menu and let "Run" 

inifile <- paste(bamp.dir, inifile, sep="")
first.ini <- read.table(file=inifile, sep=":")
ini <- as.character((first.ini[[2]]))
name <- (first.ini[[1]])
name0 <- rep("",length=length(name))
for (i in 1:length(name)){
name1<-unlist(strsplit(as.character(name[i])," "))
name1<-unlist(strsplit(as.character(name1),"\t"))
name1<-unlist(strsplit(as.character(name1),"\n"))
name2<-""
for (j in 1:length(name1))
name2<-paste(name2,name1[j],sep="")
name0[i]<-name2
}
names(ini)<-name0
rm(first.ini)

cancerdata <- paste(bamp.dir, ini["cancerdata"], sep="")
popdata <- paste(bamp.dir, ini["populationdata"], sep="")
outputfolder <- paste(bamp.dir, ini["outputfolder"], sep="")

dataorder <- as.numeric(ini["dataorder"])
numberofagegroups <- as.numeric(ini["numberofagegroups"])
numberofperiods <- as.numeric(ini["numberofperiods"])
cols <- numberofperiods
ppa <- as.numeric(ini["periodsperagegroup"])
numberofcohorts <- (ppa*(numberofagegroups - 1))+numberofperiods;
numberofpredictions <- as.numeric(ini["numberofpredictions"])
prognosis <- as.numeric(ini["prognosis"])
if (prognosis==0) numberofpredictions <- 0
quantile <- vector(length=5)
quantile[1] <- as.numeric(ini["quantile1"])
quantile[2] <- as.numeric(ini["quantile2"])
quantile[3] <- as.numeric(ini["quantile3"])
quantile[4] <- as.numeric(ini["quantile4"])
quantile[5] <- as.numeric(ini["quantile5"])
numberofquantiles <- 5
qumed <- 0
for (i in 5:1)
{
	if (quantile[i]==-1){numberofquantiles <- numberofquantiles-1}	
	if (quantile[i]==0.5){qumed <- i}
}
perplus <- 0
cohplus <- 0

ageblock <- as.numeric(ini["ageblock"])
periodblock <- as.numeric(ini["periodblock"])
if (periodblock >7)
{
	periodblock <- periodblock - 7
	perplus <- 1
}
cohortblock <- as.numeric(ini["cohortblock"])
if (cohortblock >7)
{
	cohortblock <- cohortblock - 7
	cohplus <- 1
}
rw <- max(ageblock,max(periodblock,cohortblock))

if (is.na(numberofpredictions)){numberofpredictions<-0}

# labels for axis
periodaxis<-vector(length=as.integer((numberofperiods+numberofpredictions)/3))
temp<-0
if (yearsperperiod==1)
{
for (i in seq(1,(numberofperiods+numberofpredictions),by=3))
{
temp<-temp+1
periodaxis[temp]<-as.character(startingperiod+i-1)
}
}
if (yearsperperiod>1)
{
for (i in seq(1,(numberofperiods+numberofpredictions),by=3))
{
temp<-temp+1
periodaxis[temp]<-paste((startingperiod+((i-1)*yearsperperiod)),"-",(startingperiod+(i*yearsperperiod)-1))
}
}

ageaxis<-vector(length=as.integer(numberofagegroups/3))
temp<-0
if ((ppa*yearsperperiod)==1)
{
for (i in seq(1,(numberofagegroups),by=3))
{
temp<-temp+1
ageaxis[temp]<-as.character(startingage+i-1)
}
}
if ((ppa*yearsperperiod)>1)
{
for (i in seq(1,numberofagegroups,by=3))
{
temp<-temp+1
ageaxis[temp]<-paste((startingage+((i-1)*ppa*yearsperperiod)),"-",(startingage+(i*ppa*yearsperperiod)-1))
}
}
yearsperc <- (ppa+1)*yearsperperiod
startingcohort<-startingperiod-(startingage+(ppa*yearsperperiod*numberofagegroups))
cohortaxis<-vector(length=as.integer(numberofcohorts/10))
temp<-0
for (i in seq(1,numberofcohorts,by=10))
{
temp<-temp+1
cohortaxis[temp]<-paste(startingcohort+((i-1)*yearsperperiod),"-",startingcohort+((i-1)*yearsperperiod)+yearsperc-1)
}


if (psfile==1) postscript(file=paste(outputfolder,"bamp.ps",sep=""), horizontal=F)

linien <- c(1,1,1,1)
if (numberofquantiles==3) linien <- c(3,1,3)
if (numberofquantiles==5) linien <- c(3,5,1,5,3)

if (rw==1){
# plots of age period and cohort effects, only when rw=1
theta <- matrix(scan(paste(outputfolder,"theta.txt",sep="")),ncol=numberofquantiles)
phi <- matrix(scan(paste(outputfolder,"phi.txt",sep="")),ncol=numberofquantiles)
psi <- matrix(scan(paste(outputfolder,"psi.txt",sep="")),ncol=numberofquantiles)
if (perplus==1) phibeta <- matrix(scan(paste(outputfolder,"beta-per.txt",sep="")),ncol=numberofquantiles)
if (cohplus==1) psibeta <- matrix(scan(paste(outputfolder,"beta-coh.txt",sep="")),ncol=numberofquantiles)
par(mfrow=c((3+perplus+cohplus),1))

plot.ts(theta,lty=linien,col=c(1),axes=FALSE,plot.type="single",ann=FALSE)
title(main="Age Period and Cohort effects")
axis(2)
title(sub="Age")
axis(1, at=seq(1,numberofagegroups,by=3),labels=ageaxis)
box() 
 

plot.ts(phi,lty=linien,col=c(1),plot.type="single",ann=FALSE,axes=FALSE)
if (perplus==1)
{
	title(sub="Period x Covariate")
}
else
{
	title(sub="Period")
}
axis(2)
axis(1, at=seq(1,numberofperiods+numberofpredictions,by=3),labels=periodaxis)
if (numberofpredictions>0)abline(v=numberofperiods)
box() 

if (perplus==1)
{
	plot.ts(phibeta,lty=linien,col=c(1),plot.type="single",ann=FALSE,axes=FALSE)
	title(sub="Period")
	axis(2)
	axis(1, at=seq(1,numberofperiods+numberofpredictions,by=3),labels=periodaxis)
	if (numberofpredictions>0)abline(v=numberofperiods)
	box()
} 

plot.ts(psi,lty=linien,,col=c(1),plot.type="single",ann=FALSE,axes=FALSE)
if (numberofpredictions>0)abline(v=dim(psi)[1]-numberofpredictions)
axis(2)
if (cohplus==1)
{
	title(sub="Cohort x Covariate")
}
else
{
	title(sub="Cohort")
}
axis(1, at=seq(1,numberofcohorts,by=10),labels=cohortaxis)
box() 

if (cohplus==1)
{
	plot.ts(psibeta,lty=linien,,col=c(1),plot.type="single",ann=FALSE,axes=FALSE)
	title(sub="Cohort")
	if (numberofpredictions>0)abline(v=dim(psi)[1]-numberofpredictions)
	axis(2)
	axis(1, at=seq(1,numberofcohorts,by=10),labels=cohortaxis)
	box() 
}

}



if (dataorder==0)
{
y <- t(matrix(scan(cancerdata),ncol=numberofagegroups,nrow=cols))
n <- t(matrix(scan(popdata),ncol=numberofagegroups,nrow=cols))
if (prognosis==3) 
	{
	nn <- t(matrix(scan(popdata),ncol=numberofagegroups,nrow=cols+numberofpredictions))
	}
}
if (dataorder==1)
{
y <- matrix(scan(cancerdata),nrow=numberofagegroups,ncol=cols)
n <- matrix(scan(popdata),nrow=numberofagegroups,ncol=cols)
if (prognosis==3) 
	{
	nn <- matrix(scan(popdata),nrow=numberofagegroups,ncol=cols+numberofpredictions)
	}
}

# plot predicted mortality propability for each agegroup
pr <- array(scan(paste(outputfolder,"pr.txt",sep="")),c(numberofperiods+numberofpredictions,numberofagegroups,numberofquantiles))
par(mfrow=c(4,2))
if (prognosis==3)
{
ypred<-array(scan(paste(outputfolder,"y<-pred.txt",sep="")),c(numberofpredictions,numberofagegroups,numberofquantiles))
}
for (i in 1:numberofagegroups)
{
plot.ts(pr[,i,],lty=linien,ylim=c(0,max(pr[,i,],(y[i,]/n[i,]))),,col=c(1),ann=FALSE,plot.type="single",axes=FALSE)
title(sub=paste("agegroup",(startingage+((i-1)*ppa)),"-",(startingage+(i*ppa)-1)))
for (j in 1:cols)
points(j,(y[i,j]/n[i,j]),pch="X")
if (prognosis==3)
{
	for (j in 1:numberofpredictions)
		points(j+numberofperiods,ypred[j,i,qumed]/nn[i,numberofperiods+j],pch="o")
}
title(xlab="period",ylab="pr")
if (numberofpredictions>0) abline(v=numberofperiods)
axis(2)
axis(1, at=seq(1,numberofperiods+numberofpredictions,by=3),labels=periodaxis)
box() 

}

par(mfrow=c(2,1))

# plot ysumpred
if (prognosis==3)
{
 ysum0 <- apply(y,2,sum)
 ysum <- matrix(ncol=numberofperiods+numberofpredictions,nrow=numberofquantiles)
 for (i in 1:numberofquantiles)
	{
		ysum[i,1:numberofperiods] <- ysum0
	}
 ysum[,(numberofperiods+1):(numberofperiods+numberofpredictions)] <- as.matrix(read.table(paste(outputfolder,"ysum<-pred.txt",sep="")))

plot.ts(t(ysum),lty=linien,sub="predicted cases",col=c(1),ann=FALSE,plot.type="single",axes=F)
abline(v=numberofperiods)
axis(2)
axis(1, at=seq(1,numberofperiods+numberofpredictions,by=3),labels=periodaxis)
box() 
}


# plots of z 
z<-array(scan(paste(outputfolder,"z.txt",sep="")),c(numberofperiods,numberofagegroups,numberofquantiles))

schalter<-0
for (i in 1:numberofperiods)
for (j in 1:numberofagegroups)
{
if (z[i,j,1]>0){
 if (schalter==1){points(j,i,pch="-")}
 if (schalter==0){plot(j,i,pch="-",ylim=c(1,numberofperiods),xlim=c(1,numberofagegroups),ylab="period",xlab="agegroup",axes=FALSE);schalter<-1}
}
}
for (i in 1:numberofperiods)
for (j in 1:numberofagegroups)
{
if (z[i,j,numberofquantiles]<0){
 if (schalter==1){points(j,i,pch="+")}
 if (schalter==0){plot(j,i,pch="+",ylim=c(1,numberofperiods),xlim=c(1,numberofagegroups),ylab="period",xlab="agegroup",axes=FALSE);schalter<-1}
}
}
if (schalter==1)
{
axis(1, at=seq(1,numberofagegroups,by=3),labels=ageaxis)
axis(2, at=seq(1,numberofperiods+numberofpredictions,by=3),labels=periodaxis)
box() 
title(sub="significant overdispersion")
}

if (psfile) dev.off()
#!lp bamp.ps
