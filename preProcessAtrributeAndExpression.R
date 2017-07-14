#This program is Free for Non-Commercial Use.
#Author: Herty Liany. email: lianyh@gis.a-star.edu.sg, e0146315@u.nus.edu
#Copyright Year(2017)
#############################################################################

args<-commandArgs(TRUE)
data_attr<-args[1]  #your gene expression data's sample ID attributes, refer to the example folder for the format
dataPath<-args[2] #your gene expression data, first column is geneName/Symbol
out<-args[3]   #this is your output folder, do end with slash /

attr=read.table(data_attr,check.names=FALSE)
mdat<-as.matrix(read.table(dataPath,header=T,check.names=FALSE))
mdat <-na.omit(mdat)

sampleName=colnames(mdat)
geneName=mdat[,1]
nrowData=nrow(mdat)
nAttr=ncol(attr)
colAttributeheader=paste(matrix(colnames(attr)),collapse=",")

combs=2^nAttr
un=unique(attr)

groups_idx<-matrix(rep(c(NA)),nrow(un),2)
groups<-vector("list",nrow(un))

a<-array(as.matrix(attr),dim=c(nrow(attr),ncol(attr)))
count=1

indexes<-NULL
newSampleAttr=NULL
newExpDataMatrix=NULL
newExpDataMatrix=cbind(newExpDataMatrix,geneName)
newGroupInfo=matrix(rep(c(NA)),nrow(un),2)
for (i in 1:(nrow(un)))
{
	cond=un[i,]
	x<-which(apply(a,1,function(x)all(x==cond)))
	groups[[i]]=rownames(attr[x,])
	
	newGroupInfo[i,1]=paste("Group",i,":",collapse=" ")
	tempVal=paste(matrix(cond),collapse=",")
	newGroupInfo[i,2]=paste(colAttributeheader,"[",tempVal,"]",sep="")
	
	groups_idx[i,1]=count
	groups_idx[i,2]=(count -1) + length(x)
	count = count + length(x)
	
	newSampleAttr=rbind(newSampleAttr,attr[x,])
	y=match(groups[[i]],sampleName)
	newExpDataMatrix=cbind(newExpDataMatrix,mdat[,y])
}

f=paste(out, "processedSampleAttributeMDCX.txt", sep = "")
write.table(t(newSampleAttr),f,sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
f=paste(out, "processedExpressionMDCX.txt", sep = "")
write.table(newExpDataMatrix,f,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
f=paste(out, "processedGroupsMDCX.txt", sep = "")
write.table(groups_idx,f,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
f=paste(out, "GroupsInfoMDCX.txt", sep = "")
write.table(newGroupInfo,f,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

print("DONE")

