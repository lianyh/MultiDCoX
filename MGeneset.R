#This program is Free for Non-Commercial Use.
#Author: Herty Liany. email: lianyh@gis.a-star.edu.sg, e0146315@u.nus.edu
#Copyright Year(2017)
#############################################################################

args<-commandArgs(TRUE)


if(length(args) < 5) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --arg1 = row number idx or genelist   - a file
      --arg2 = dataset attribute   - a file
      --arg3 = gene expression dataset   - a file
      --arg4 = dataset groups(columns) index    - a file
      --arg5 = outputFolder/    - folder path ended with /

      --help             
 
      Example:
      ./MGeneset.R --arg1=genelist --arg2=processedSampleAttributeMDCX.txt --arg3=processedExpressionMDCX.txt --arg4=processedGroupsMDCX.txt --arg5=outputFolder/ \n\n")
 
  q(save="no")
}

genelistFile=args[1]
data_attr=args[2]
dataPath=args[3]
groups=args[4]
out=args[5]

################################################## FUNCTION #######################################################################


checkIfNumeric <- function (x) {
  out <- tryCatch(scan(x,what=numeric()), error = function(e) NULL)
  return(out)
}
################################################## END FUNCTION #######################################################################


################################################## MAIN ###########################################################################

	attr=read.table(data_attr,header=T,check.names=FALSE)

	g_attr<-as.matrix(attr)
	
	mydat<-read.table(dataPath,header=T,check.names=FALSE)
	mydat <-na.omit(mydat)
	gnlist<-mydat[,1]
	mydat<-mydat[,-1]
	myGeneList<-NULL
	gidx<-NULL

	status=checkIfNumeric(genelistFile)
	if (is.null(status)) {
		myGeneList<-read.table(genelistFile,header=F,check.names=FALSE)
		gidx=match(as.matrix(myGeneList),as.matrix(gnlist))
	}else {
	gidx=p=as.vector(scan(genelistFile,what=numeric()))
	myGeneList=gnlist[gidx]}

	dat<-as.matrix(mydat)
	dat<-dat[gidx,]
	num_row = nrow(dat)
	exprDataCol=colnames(dat)


	#nrAttr refers to number of covariates you have generated in newSampleAttributeMDCX.txt (first column)
	nrAttr=nrow(g_attr)
	rAttrName=row.names(g_attr)

	group=read.table(groups,header=F,row.names=NULL,sep="\t",quote="")
	groups_idx <- matrix(as.matrix(group),ncol=ncol(group))
	groups_comb<-matrix(rep(c(NA)),nrow(groups_idx),ncol(groups_idx))

	myTimestamp=format(Sys.time(), "%H_%M_%S_%d_%m_%Y")
	log_summary=paste(out,"GenesetResults_",myTimestamp,".txt", sep = "")
	str=paste("Method","Run_Number","GIDX","Genes","Total_Genes",sep="\t")
	

	All_coeff <- vector()
	fm<-matrix()
	lm_sampling=NULL
	str_lm2=NULL

	print("START")
	colCombinations=mapply(x=groups_idx[,1],y=groups_idx[,2],function(x,y){ ((y-x)+1) * (((y-x)+1)-1) /2 })
	length_lm=sum(colCombinations)
	z<-matrix(rep(c(0)),num_row,length_lm)
	
	for( h in 1:nrAttr)
	{
		# initialiaze variables for all possible covariates, i.e.: ER, P53,etc
		assign(paste(rAttrName[h],"_lm",sep=""),NULL)
		assign(paste(rAttrName[h],"_max_seeds",sep=""),matrix())
		assign(paste(rAttrName[h],"_min_seeds",sep=""),matrix())
		assign(paste(rAttrName[h],"_max_coeff",sep=""),NULL)
		assign(paste(rAttrName[h],"_min_coeff",sep=""),NULL)
		assign(paste(rAttrName[h],"_max_idx",sep=""),NULL)
		assign(paste(rAttrName[h],"_min_idx",sep=""),NULL)
		
		#initialize header results for all covariates (coeff and pvalues columns)
		str=paste(str,paste(rAttrName[h],"_coefficient",sep=""),paste(rAttrName[h],"_pvalue",sep=""),sep="\t")
		
	}
	
	#append permutation columns header after all coeff and pvalue columns
	for( h in 1:nrAttr)
	{
		str=paste(str,paste("Count_",rAttrName[h],sep=""),paste("pvalue_",rAttrName[h],sep=""),sep="\t")
		str_lm2=paste(str_lm2,paste(rAttrName[h],"_lm",sep=""),sep="+")
	}
	
	#write out the header results to be output
	write(str,file=log_summary)

	idx=0
	#loop through all possible combinations of all covariates (in all individuals) generated in newGroupsMDCX.txt
	for(k in 1:nrow(groups_idx))
	{
		start=groups_idx[k,1]
		end=groups_idx[k,2]
		groups_comb[k,1]=idx+1
		
		#loop through each sample in the group (newGroupsMDCX.txt) and do pair-wise calculation of eilm and stored in z matrix
		for(i in start:(end-1))
		{
			for(j in (i+1):end) 
			{
				idx=idx+1	
				#nrAttr refers to the number of covariates you have generated in newSampleAttributeMDCX.txt (first column)
				#g_attr[h,i] refers to attribute, i.e. -1,0 or 1.
				for( h in 1:nrAttr)   
				{
					f=paste(rAttrName[h],"_lm",sep="")
					if (g_attr[h,i] == g_attr[h,j])
					{
						eq<-paste(f,"[",idx,"] <- g_attr[h,i]",sep="")
						eval(parse(text=eq))
					}
				}
				#here is for the z matrix(which is the main eilm matrix)
				x =dat[,i]
				y =dat[,j]
				z[,idx]=(x-y)
			}
		}

		groups_comb[k,2]=idx
		#append all the individuals in the group (newGroupsMDCX.txt)
		lm_sampling<-append(lm_sampling,groups_comb[k,1]:groups_comb[k,2])
	}
	
	#assign all lm_sampling to the attributes, i.e: ER_lm=ER_lm[lm_sampling], etc
	for( h in 1:nrAttr)
	{
		f=paste(rAttrName[h],"_lm",sep="")
		eq<-paste(f,"<-",f,"[",'lm_sampling',"]",sep="")
		eval(parse(text=eq))
	}
	#to z matrix as well
	z=z[,lm_sampling]

        alm = (apply(z,2,mean))^2
	f=paste("o=lm(alm~",str_lm2,")",sep="")
	eval(parse(text=f))
	b=summary(o)

	#permutation lm i.e permute_ER_lm + permute_P53_lm + permute_Grade_lm
	str_plm=""
	for (h in 1:nrAttr)
	{
		#initialize count for all covariates =0
		f=paste("count_",rAttrName[h],"<- 0",sep="")
		eval(parse(text=f))
	}

	#do permutation 100 times to get the permutation results
	for(m in 1:100)
	{
		for( h in 1:nrAttr)
		{
			f=paste("permute_",rAttrName[h],"_lm <- sample(",rAttrName[h],"_lm)",sep="")
			eval(parse(text=f))
			str_plm=paste(str_plm,paste("permute_",rAttrName[h],"_lm",sep=""),sep="+")
		}
		
		f=paste("new_o=lm(alm~",str_plm,")",sep="")
		eval(parse(text=f))
		new_b=summary(new_o)
		
		for( h in 1:nrAttr)
		{
			if(o$coeff[h+1] > 0 && new_o$coeff[h+1] > o$coeff[h+1] && b$coefficients[h+1,4] < new_b$coefficients[h+1,4])
			{
				f=paste("count_",rAttrName[h],"<- count_",rAttrName[h]," + 1",sep="")
				print("f pos")
				print(f)
				eval(parse(text=f))
			}
			if(o$coeff[h+1] < 0 && new_o$coeff[h+1] < o$coeff[h+1] && b$coefficients[h+1,4] < new_b$coefficients[h+1,4])
			{
				f=paste("count_",rAttrName[h],"<- count_",rAttrName[h]," + 1",sep="")
				
				eval(parse(text=f))
			}
		}
	}

	str=paste("GenesetStats","-",paste(gidx,collapse=","),paste(as.matrix(myGeneList),collapse=","),length(gidx),sep="\t")
	for( h in 1:nrAttr) { str=paste(str,eval(parse(text=paste("o$coeff[h+1]",sep=""))),eval(parse(text=paste("b$coefficients[h+1,4]",sep=""))),sep="\t") }
	
	#permutation count and pvalue for each covariate
	for( h in 1:nrAttr) {str=paste(str,eval(parse(text=paste("count_",rAttrName[h],sep=""))),eval(parse(text=paste("count_",rAttrName[h],"/100",sep=""))),sep="\t")}
	
	#############################################
	#WRITING THE RESULTS OUT TO THE OUTPUT FILE
	#############################################
	write(str,file=log_summary,append=TRUE)
	
	print("PRINTING COVARIATE/FACTORWISE CO-EXPRESSION PLOT...")
	#############################################
	#PRINT COVARIATE/FACTORWISE CO-EXPRESSION PLOT
	#############################################
	tdp=paste(out, "COVARIATES_Plot_",myTimestamp,".png", sep = "")
	png(tdp,width=1200, height=700)
	par( mfrow = c( nrAttr, 3 ) )
	for( h in 1:nrAttr)
	{
		covariate=g_attr[h,]
		c=unique(covariate)
		covar_name=row.names(g_attr)[[h]]
	
		for( k in 1:length(c))
		{

			colnames=names(covariate[covariate==c[k]])
			w=match(as.matrix(colnames),as.matrix(exprDataCol))

			y=dat[,w]
			plot(y[1,],type="l",col="black",ylim=c(-12,12),main=paste(covar_name," [",c[[k]],"]",sep=""),xlab="Samples",ylab="gene Expression")
			for(i in 2:length(gidx)){lines(y[i,],col="black")}
			p=as.matrix(apply(y,2,mean))
			lines(p,col="gray",lwd=1)
		}
	}
	dev.off()

	print("PRINTING STRATUMWISE CO-EXPRESSION PLOT...")
	#############################################
	#PRINT STRATUMWISE CO-EXPRESSION PLOT
	#############################################
	tdp=paste(out, "STRATUM_Plot_",myTimestamp,".png", sep = "")
	png(tdp,width=1200, height=700)
	d=nrow(group)
	d=ifelse(d%%3==0,(d/3),as.integer(d/3)+1)

	par( mfrow = c(d, 3 ) )
	for( h in 1:nrow(group))
	{
		startg=groups_idx[h,1]
		endg=groups_idx[h,2]
		y=dat[,startg:endg]
		plot(y[1,],type="l",col="black",ylim=c(-12,12),main=paste("Group ",h," [col ",startg,"-",endg,"]",sep=""),xlab="Samples",ylab="gene Expression")
		for(i in 2:length(gidx)){lines(y[i,],col="black")}
		p=as.matrix(apply(y,2,mean))
		lines(p,col="gray",lwd=1)
	}
	dev.off()

print("COMPLETED")

################################################## END MAIN ###########################################################################
