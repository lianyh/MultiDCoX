#This program is Free for Non-Commercial Use.
#Author: Herty Liany. email: lianyh@gis.a-star.edu.sg, e0146315@u.nus.edu
#Copyright Year(2017)
#############################################################################

rAttrName <<- NULL
gnlist  <<- NULL
gidxLength  <<- 5  #default (the bigger the number of genes to be in a geneset, it will render to stricter/fewer output or no results. Higher number will render to higher computational time)
numGenesets <<- 6  #default (this is the maximum number of iterations of top seed genes if the previous seed gene cannot find the co-expressed geneset. Higher number will render to higher computational time)
AllCovariatesCoeff <<- vector()
log_summary <<- NULL

args<-commandArgs(TRUE)
if(length(args) < 4) {
  args <- c("--help")
}


## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --arg1 = dataset attribute   - a file
      --arg2 = gene expression dataset   - a file
      --arg3 = dataset groups(columns) index    - a file
      --arg4 = outputFolder/    - folder path ended with /

      --arg5 = minimum number of co-expressed genes(OPTIONAL)  - numeric (5 default)
      --arg6 = maximum number of iterations of top seed genes(OPTIONAL)   - numeric (6 default)

      --help             
 
      Example:
      ./MDCoX.R --arg1=processedSampleAttributeMDCX.txt --arg2=processedExpressionMDCX.txt --arg3=processedGroupsMDCX.txt --arg4=outputFolder/ \n\n")
 
  q(save="no")
}


data_attr=args[1]
dataPath=args[2]
groups=args[3]
out=args[4]

if(length(args) ==5) {
	gidxLength=args[5]
}
if(length(args) ==6 ){
	gidxLength=args[5]
	numGenesets=args[6]
}




####################################################### FUNCTION #################################################################################
fxGetCovariateGeneSet<-function(prefix,lm_seeds,initGene,num_row,gnlist,z,rAttrName,threshold_all,fm)
{
	nt=nrow(lm_seeds)
	nrAttr=length(rAttrName)
	str_lm2=""
	
	for( h in 1:nrAttr)
	{
		#retrieve all lm for each covariate, i.e. ER_lm, P53_lm, Grade_lm
		f=paste(rAttrName[h],"_lm <- fm[,",h+1,"]",sep="")
		eval(parse(text=f))

		#init string i.e ER_lm,P53_lm,Grade_lm
		str_lm2=paste(str_lm2,paste(rAttrName[h],"_lm",sep=""),sep="+")
	}

	t <- 0
	count <- 0
	while((!is.null(nt)) && t < nt && count < numGenesets) { 
		t=t+1
		g_idx=c(lm_seeds[t,1],lm_seeds[t,2])
		previous_coeff=0
		for(c in 1:num_row)
		{
			m=initGene[c]
			d=gnlist[g_idx]
			currGene=gnlist[m]
			if ( (!(any(g_idx==m)))  )
			{
			   if ((!(any(d==currGene))) && (!is.na(currGene)) && (currGene!="NA") )
			   {
				temp2_g_idx<-append(g_idx,m)
				subMat=z[temp2_g_idx,]
				alm = (apply(subMat,2,mean))^2
				f=paste("o=lm(alm~",str_lm2,")",sep="")
				eval(parse(text=f))
				b=summary(o)

				#check the coefficient for that particular covariate max or min (need to loop through all covariates, 
				#to make sure the profile (i.e. Min for MinGrade, MinP53, and Max for MaxGrade, MaxP53) is matching, pvalue is significant and  
				#the current coeff is greater than previous one.
				#for MinGrade, you can assume is grade1, while MaxGrade is grade 2 or 3.
				for (h in 1:nrAttr)
				{
					minMeth <- paste("Min",rAttrName[h],sep="")
					maxMeth <- paste("Max",rAttrName[h],sep="")
					if (prefix == minMeth || prefix == maxMeth)
					{
						f=paste("var_coeff=(o$coeff[h+1])",sep="")
						eval(parse(text=f))
						f=paste("pval=(b$coefficients[h+1,4])",sep="")
						eval(parse(text=f))
						if ((var_coeff / threshold_all) > 1)
						{
							if(pval < 0.05 && abs(var_coeff) > abs(previous_coeff))
							{
								g_idx<-append(g_idx,m)
								previous_coeff=var_coeff
								break
							}	
						}
					}
				}
			   }
			}
		}
		
		#reshuffle the indices for every next new gene pair seed
		initGene<-sample(num_row,num_row)
		subMat=z[g_idx,]
		alm = (apply(subMat,2,mean))^2
		f=paste("o=lm(alm~",str_lm2,")",sep="")
		eval(parse(text=f))
		b=summary(o)

		############################################################################
		#do Augmentation,Filtering for each covariate, lastly print out the results#
		############################################################################
		g_idx<-fxMainAugmentationAndFiltering(z,gnlist,prefix,t,rAttrName,g_idx,num_row,o,b,threshold_all,str_lm2,fm)
			

		if(length(g_idx) >= gidxLength) count = count +1

		lm_seeds=matrix(lm_seeds[!lm_seeds[,1] %in% g_idx],ncol=3)
		lm_seeds=matrix(lm_seeds[!lm_seeds[,2] %in% g_idx],ncol=3)
	
		nt=nrow(lm_seeds)
	}
}

fxMainAugmentationAndFiltering<-function(z,gnlist,prefix,t,rAttrName,gidx,num_row,o,b,threshold_all,str_lm2,fm)
{
	nrAttr=length(rAttrName)
	profileMatrix=matrix(rep(c(NA)),nrAttr,3)

	#permutation lm i.e permute_ER_lm + permute_P53_lm + permute_Grade_lm
	str_plm=""
	
	for (h in 1:nrAttr)
	{
		profileMatrix[h,1]=0 #initialize, 0 means the whole profile is insignificant, 1 means significant ( greater than coeffcient threshold and pvalue <0.01)
		profileMatrix[h,2]=as.numeric(o$coeff[h+1]) #coeff value
		profileMatrix[h,3]=as.numeric(b$coefficients[h+1,4])  #pvalue
		#status =indicate it is max or min
		status <- ifelse(profileMatrix[h,2] > 0, 1, 0)
		
		f=paste(rAttrName[h],"_maximization <- status",sep="")
		eval(parse(text=f))
		
		if((profileMatrix[h,2] / threshold_all > 1) && profileMatrix[h,3] <= 0.01) 
		{
			if(status==1){ profileMatrix[h,2] = max(mean(profileMatrix[h,2] + threshold_all), threshold_all) }
			else {profileMatrix[h,2] = min( mean(profileMatrix[h,2] + threshold_all), threshold_all)}
			profileMatrix[h,1]=1
		}

	}
	eblm=(apply(z[gidx,],2,mean))
	#agumentation (adding of all possible genes to the geneset)
	gidx <- augmentCoefficient(z,num_row,gidx,gnlist,eblm,rAttrName,profileMatrix,prefix,fm,str_lm2)

	subMat=z[gidx,]
	alm = (apply(subMat,2,mean))^2
	f=paste("o=lm(alm~",str_lm2,")",sep="")
	eval(parse(text=f))
	b=summary(o)
	signFlag=1

	#check and update profile after augmentation
	for (h in 1:nrAttr)
	{
		#To double confirm the final gidx set is coherent with the run max/min method
		minMeth <- paste("Min",rAttrName[h],sep="")
		maxMeth <- paste("Max",rAttrName[h],sep="")
		if (prefix == minMeth && as.numeric(o$coeff[h+1]) > 0 ){signFlag=0}
		if (prefix == maxMeth && as.numeric(o$coeff[h+1]) < 0 ){signFlag=0}
		
		#update profile after augmentation
		profileMatrix[h,2]=as.numeric(o$coeff[h+1])
		profileMatrix[h,3]=as.numeric(b$coefficients[h+1,4])
	}
	#discard all geneset with length <  gidxLength (default =5)
	if(length(gidx) >= gidxLength && signFlag==1) #prior filtering, #genes should greater than equal to default settings
	{
		#do filtering after augmentation
		gidx<- filterCoefficient(z,length(gidx),gidx,rAttrName,profileMatrix,prefix,fm,str_lm2)
		if(length(gidx) >= gidxLength) #make sure the final geneset is greater than the default #genes
		{
			subMat=z[gidx,]
			alm = (apply(subMat,2,mean))^2
			f=paste("o=lm(alm~",str_lm2,")",sep="")
			eval(parse(text=f))
			b=summary(o)
		
			for (h in 1:nrAttr){
				#initialize count for all covariates =0
				f=paste("count_",rAttrName[h],"<- 0",sep="")
				eval(parse(text=f))

				#store all coefficients for density plot
				AllCovariatesCoeff <<- append(AllCovariatesCoeff,b$coefficients[h+1,1])
			}
			#do permutation 100 times to get the permutation results
			for(m in 1:100)
			{
				for( h in 1:nrAttr)
				{
					f=paste("permute_",rAttrName[h],"_lm <- sample(",rAttrName[h],"_lm)",sep="")
					eval(parse(text=f))
					#permutation lm i.e permute_ER_lm + permute_P53_lm + permute_Grade_lm
					#assign all covariates's lm model from fm
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
						eval(parse(text=f))
					}
					if(o$coeff[h+1] < 0 && new_o$coeff[h+1] < o$coeff[h+1] && b$coefficients[h+1,4] < new_b$coefficients[h+1,4])
					{
						f=paste("count_",rAttrName[h],"<- count_",rAttrName[h]," + 1",sep="")
						eval(parse(text=f))
					}
				}
			}

			str=paste(prefix,t,paste(gidx,collapse=","),paste(gnlist[gidx],collapse=","),length(gidx),sep="\t")
			for( h in 1:nrAttr) { str=paste(str,eval(parse(text=paste("o$coeff[h+1]",sep=""))),eval(parse(text=paste("b$coefficients[h+1,4]",sep=""))),sep="\t") }
			#permutation count and pvalue for each covariate
			for( h in 1:nrAttr) {str=paste(str,eval(parse(text=paste("count_",rAttrName[h],sep=""))),eval(parse(text=paste("count_",rAttrName[h],"/100",sep=""))),sep="\t")}
			
			#############################################
			#WRITING THE RESULTS OUT TO THE OUTPUT FILE
			#############################################
			write(str,file=log_summary,append=TRUE)

		}
	}
	return(gidx)
}


#the idea of augmentcoeff is gene to be added (augmented) to the geneset must satisfy both critera:
#1) (The average of original geneset's diff expression from Samplei to Samplej + a gene k's diff expression)/2 should be greater than the previous genetset profile prior adding of a gene k
#2) the mean of (the original geneset's diff expression from Samplei to Samplej  + a gene k's diff expression) should be greater than the previous genetset profile prior adding of a gene k
augmentCoefficient<- function(z,num_row,gidx,gnlist,eblm,rAttrName,profileMatrix,prefix,fm,str_lm2)
{
	nrAttr=length(rAttrName)
	#assign all covariates's lm model from fm
	for( h in 1:nrAttr)
	{
		f=paste(rAttrName[h],"_lm <- fm[,",h+1,"]",sep="")
		eval(parse(text=f))
	}
	#to get all the genelist index other than the gidx
	k=1:num_row
	k=k[!(k %in% gidx)]
	g=1:length(k)  #the idea to add this is for t(z[k,]), alm[final_all,]  is correnpond to the number of genes minus gidx
	
	#criteria #1
	#eblm is the average of the Geneset's expression value 
	#subMat below is the average of  gidx alm (eblm)+ a particular gene k
	subMat=t((eblm + t(z[k,]))/2)
	alm = (subMat)^2
	qr_coeff=qr.coef(qr(fm),y=t(alm))
	
	#criteria #2
	#eblm_gidx below apply "sum", and average the gene length + 1, this is the default way to find the alm
	eblm_gidx=(apply(z[gidx,],2,sum))
	#the original method to find alm of the geneset's eilm + a gene k's eilm
	subMat=t((eblm_gidx + t(z[k,]))/(length(gidx)+1)) #change from eblm_gidx to gidx
	alm_all=(subMat)^2
		
	qr_coeff_main_all=qr.coef(qr(fm),y=t(alm_all))

	final_all<-vector()

	#loop through to find the subset within the subset of all covariates profile with status=1, ex:er_status=1,p53_status=1, find subset results of er and then p53 that matched the initial profile
	for (h in 1:nrAttr)
	{
		status=profileMatrix[h,1]  #if profile is significant for that particular covariate then status 1 otherwise is 0
		coeff=profileMatrix[h,2]   #coeff value
	
		##################
		# the profile that is significant only will be chosen for augmentation (status =1)
		if(status ==1)
		{
			xx_coeff = qr_coeff[h+1,]
			#(this is to satisfy criteria #1)make sure the new coefficient is greater than the previous profile #criteria 1
			final_idx1 = g[(xx_coeff / coeff) > 1] 

			#(this is to satisfy criteria #2) make sure the new coefficient is greater than the previous profile #criteria 2
			xx_coeff_all=qr_coeff_main_all[h+1,]
			final_idx_all = g[(xx_coeff_all / coeff) > 1] 

			#make sure the added gene should satisfy both criteria 1 and 2
			final_idx= intersect(final_idx1,final_idx_all)

			final_all<-append(final_all,final_idx)
		}
		
	}
	final_all=unique(final_all)
	
	#after getting all the genes with coefficients value higher than previous one, 
	#we check for PVALUE and profile significance after adding that particular gene into the geneset
	if(length(final_all) > 0)
	{
		b <- NULL
		#only 1 gene found to be added 
		if(length(final_all) == 1) 
		{
			f=paste("o=summary(lm(alm[final_all,] ~ ",str_lm2,"))",sep="") 
			eval(parse(text=f))
			b=o
		}else{ #more than 1 gene 
			f=paste("o=summary(lm(t(alm[final_all,]) ~ ",str_lm2,"))",sep="") 
			eval(parse(text=f))
		}
		
		#check the geneset profile again and make sure the new pvalue (after the added genes) is more significant than the previous profile's pvalue
		for(a in 1: length(final_all))
		{
			
			if(length(final_all) > 1) {
				b=o[[a]]
			}
			
			currGene=gnlist[final_all[a]]
			d=gnlist[gidx]
			count=0
			if((!(currGene %in% d)) && (!is.na(currGene)) && (currGene != "NA") ) # we exclude probes which come from the same gene (we chose only one of it)
			{
				for (h in 1:nrAttr)
				{
					status=profileMatrix[h,1]  #status profile is significance =1, non significant=0
					if(status ==1)
					{
						pval_xx=(b$coefficients[h+1,4])
						coeff_xx=(b$coefficients[h+1,1])
						if(!is.null(coeff_xx) && (coeff_xx / profileMatrix[h,2] > 1) && pval_xx < profileMatrix[h,3])
						{
							count=count+1
						}
					}
				}
				#to make sure the adding of a particular gene still maintain the profile significance
				if(count == sum(profileMatrix[1:nrAttr,1]) )
				{
					 #assign back the seq g to its former genelist index k
					gidx<-append(gidx,k[final_all[a]])

					
				}
			}
		}
	}
	return(gidx)
}


#to filter away some genes in the geneset if possible, the filtering will be applicable only if the geneset's overall coefficient and pvalue are getting improved.
filterCoefficient<- function(z,num_row,gidx,rAttrName,profileMatrix,prefix,fm,str_lm2)
{
	previous_coeff=0
	var_coeff_er=0
	pval_er=0
	var_coeff_p53=0
	pval_p53=0
	var_coeff_grade=0
	pval_grade=0

	nrAttr=length(rAttrName)
	#assign all covariates's lm model from fm
	for( h in 1:nrAttr)
	{
		f=paste(rAttrName[h],"_lm <- fm[,",h+1,"]",sep="")
		eval(parse(text=f))
	}
	
	#filter or remove all genes in the geneset that when those genes are removed, the overall geneset's coefficient increase and pvalue remain significant
	for(a in 1:length(gidx))
	{
		temp_gidx<-gidx[-a]
		eblm=(apply(z[temp_gidx,],2,mean))
		alm = (eblm)^2
		
		f=paste("o=lm(alm~",str_lm2,")",sep="")
		eval(parse(text=f))
		b=summary(o)

		count=0
		for (h in 1:nrAttr)
		{
			status=profileMatrix[h,1]  #status significant=1, insignificant =0
			if(status ==1)
			{
				var_coeff=(o$coeff[h+1])
				pval=(b$coefficients[h+1,4])
			
				minMeth <- paste("Min",rAttrName[h],sep="")
				maxMeth <- paste("Max",rAttrName[h],sep="")

				#prefix refers to currently optimized covariate, some covariate might not in the current optimization 
				#but still the whole profile (all covariates regardless if is in current optimization or not) shall be improved
				if (prefix == minMeth || prefix == maxMeth)
				{
					if((var_coeff / profileMatrix[h,2]) > 1.2 && pval < 0.01 && pval < profileMatrix[h,3])
					{
						count=count+1
					}
				} #else if other covariates (less stringent) which is not under current covariate optimization run
				else 
				{
					if((var_coeff / profileMatrix[h,2]) > 1 && pval < 0.01 && pval < profileMatrix[h,3])
					{
						count=count+1
					}
				}
				
			}
		}
		#to make sure the removal of a particular gene still maintain the profile significance
		if(count == sum(profileMatrix[1:nrAttr,1]))
		{
			gidx<-gidx[-a]		
		}
	}
	return(gidx)
}

################################################ END FUNCTION #######################################################################







################################################## MAIN ###########################################################################
	attr=read.table(data_attr)
	g_attr<-as.matrix(attr)

	mydat<-read.table(dataPath,header=T,check.names=FALSE)

	mydat <-na.omit(mydat)
	gnlist<-mydat[,1]
	mydat<-mydat[,-1]

	dat<-as.matrix(mydat)
	num_row = nrow(dat)

	#nrAttr refers to number of covariates you have generated in newSampleAttributeMDCX.txt (first column)
	nrAttr=nrow(g_attr)
	rAttrName=row.names(g_attr)

	group=read.table(groups,header=F,row.names=NULL,sep="\t",quote="")
	groups_idx <- matrix(as.matrix(group),ncol=ncol(group))
	groups_comb<-matrix(rep(c(NA)),nrow(groups_idx),ncol(groups_idx))

	myTimestamp=format(Sys.time(), "%H_%M_%S_%d_%m_%Y")
	log_summary=paste(out,"MultiDCoxGenesetResults_",myTimestamp,".txt", sep = "")
	str=paste("Method","Run_Number","GIDX","Genes","Total_Genes",sep="\t")
	if (!(file.exists(out)))
	{
	     dir.create(out)
	} 

	All_coeff <- vector()
	fm<-matrix()
	lm_sampling=NULL
	
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
	
	#i.e str_lm = ER_lm,P53_lm,Grade_lm
	str_lm=""
	for (h in 1:nrAttr){str_lm=paste(str_lm,paste(rAttrName[h],"_lm",sep=""),sep=",")}
	#fm is all covariates' lm in vertical order (column bind), usable for all subsequent qr.coef function
	fm = paste("fm <- cbind(rep(1,length_lm)",str_lm,")",sep="")
	eval(parse(text=fm))

	####################################################################################
	#TO FIND A THRESHOLD of COEFFICIENT cut-off for gene seed selection in later phase.
	#The threshold cut-off is based on 2% sampling of total probes
	####################################################################################
	sampling_size=as.integer(0.02 * num_row) #0.02 is 2% out of total probes
	#shuffle genes indices
	initGene=sample(num_row,num_row) #randomly shuffle the index genes
	for(k in 1:sampling_size){
		g_idx=c(initGene[k])

		#ignore gene with NA name, move to next iteration
		if (length(g_idx) == 0 || is.na(gnlist[g_idx]) || gnlist[g_idx]=="NA" ) { next }
	
		#to exclude the g_idx which paired with its own partner, i.e. gene index 24 and gene index 24 pairing should be excluded
		c_num_row=1:num_row
		c_num_row=c_num_row[c_num_row!=g_idx]
		
		#randomly select 2% of all genes
		s_size=as.integer(0.02 * length(c_num_row))
		sampling<-sample(c_num_row,s_size)
		alm=(t(z[sampling,] + z[g_idx,])/2)^2

		qr_coeff=qr.coef(qr(fm),y=alm)
		genes=gnlist[sampling]
		#append all coefficient value into All_coeff variable, which is used to derive the threshold coefficient cut-off.
		for (h in 1:nrAttr){
			eq<-paste("All_coeff <- append(All_coeff,qr_coeff[",h+1,",])",sep="")
			eval(parse(text=eq))
		}
		print("#################")
	}
	
	coeff_max=sort(All_coeff,TRUE)[10] # to choose the top 10th of the max coefficient
	coeff_min=sort(All_coeff,FALSE)[10] # to choose the top 10th of the min coefficient
	threshold_all=abs(max(coeff_max,abs(coeff_min)))/2 #ignore sign (+ or -), both ends are included for selection

	# for (-)  covariates, assign negative sign (threshold_all_min) to the threshold_all above
	threshold_all_min=threshold_all * -1
	print("threshold_all")
	print(threshold_all)
	
	########################################################################################################################################
	#To get paired gene seeds  whereby the covariates coeffient values is greater than the threshold cut-off (threshold_all)
	########################################################################################################################################
	sampling_size=as.integer(0.2 * num_row) #0.2 is 20% out of total probes
	#again, shuffle the genes indices
	initGene=sample(num_row,num_row)
	
	#####START FOR LOOP
	for(k in 1:sampling_size){
		g_idx=c(initGene[k])
		
		#ignore gene with NA name
		if (is.na(gnlist[g_idx]) || gnlist[g_idx]=="NA") { next }

		#to exclude the g_idx that pair with its own partner, i.e. 24,24 (to be excluded)
		c_num_row=1:num_row
		c_num_row=c_num_row[c_num_row!=g_idx]

		#randomly select 2% of all genes
		s_size=as.integer(0.2 * length(c_num_row))
		sampling<-sample(c_num_row,s_size)
		alm=(t(z[sampling,] + z[g_idx,])/2)^2

		qr_coeff=qr.coef(qr(fm),y=alm)
		
		#FOR EACH COVARIATE, find gene with its coefficient value (from qr.coef function) that greater than the threshold_all cut-off 
		for (h in 1:nrAttr)
		{
			eq = paste("xx_coeff = qr_coeff[",h+1,",]",sep="")
			eval(parse(text=eq))

			#################For Max (positive covariate)
			f=paste(rAttrName[h],"_max_coeff = xx_coeff[xx_coeff > threshold_all]",sep="") 
			eval(parse(text=f))
			
			f=paste(rAttrName[h],"_max_idx = sampling[xx_coeff > threshold_all]",sep="") 
			eval(parse(text=f))
			
			f=paste("genes=gnlist[",rAttrName[h],"_max_idx]",sep="")
			eval(parse(text=f))
			
			f=paste(rAttrName[h],"_max_idx=",rAttrName[h],"_max_idx[genes != 'NA']",sep="")
			eval(parse(text=f))
			
			f=paste(rAttrName[h],"_max_coeff=",rAttrName[h],"_max_coeff[genes != 'NA']",sep="")
			eval(parse(text=f))

			f=paste(rAttrName[h],"_max_seed=matrix(cbind(rep(g_idx),",rAttrName[h],"_max_idx,",rAttrName[h],"_max_coeff),length(",rAttrName[h],"_max_idx),3)",sep="")
			eval(parse(text=f))

			
			#################For Min (negative covariate) 
			f=paste(rAttrName[h],"_min_coeff = xx_coeff[xx_coeff < threshold_all_min]",sep="") 
			eval(parse(text=f))

			f=paste(rAttrName[h],"_min_idx = sampling[xx_coeff < threshold_all_min]",sep="") 
			eval(parse(text=f))

			f=paste("genes=gnlist[",rAttrName[h],"_min_idx]",sep="")
			eval(parse(text=f))

			f=paste(rAttrName[h],"_min_idx=",rAttrName[h],"_min_idx[genes != 'NA']",sep="")
			eval(parse(text=f))
			
			f=paste(rAttrName[h],"_min_coeff=",rAttrName[h],"_min_coeff[genes != 'NA']",sep="")
			eval(parse(text=f))
			f=paste(rAttrName[h],"_min_seed=matrix(cbind(rep(g_idx),",rAttrName[h],"_min_idx,",rAttrName[h],"_min_coeff),length(",rAttrName[h],"_min_idx),3)",sep="")
			eval(parse(text=f))

			################# Final first seeds for each covariate (both positive and negative covariates)
			f=paste("if (length(",rAttrName[h],"_max_idx) > 0 ) { if (!is.na(",rAttrName[h],"_max_seeds[1])) {",rAttrName[h],"_max_seeds=rbind(",rAttrName[h],"_max_seeds,",rAttrName[h],"_max_seed)} else {",rAttrName[h],"_max_seeds=",rAttrName[h],"_max_seed} }",sep="")
			eval(parse(text=f))

			f=paste("if (length(",rAttrName[h],"_min_idx) > 0 ) { if (!is.na(",rAttrName[h],"_min_seeds[1])) {",rAttrName[h],"_min_seeds=rbind(",rAttrName[h],"_min_seeds,",rAttrName[h],"_min_seed)} else {",rAttrName[h],"_min_seeds=",rAttrName[h],"_min_seed} }",sep="")
			eval(parse(text=f))
		}
		#END FOR LOOP FOR EACH COVARIATE
		print("#################")
	}
	#####END FOR LOOP
	
	##############################################################################################################################################################
	#GET THE COVARIATE "GeneSet" FROM A LIST OF SEED GENES (PAIRED GENE), THIS INCLUDDE: AUGMENTATION, FILTERING AND LASTLY PRINTING OUT THE FINAL GENESET RESULTS
	##############################################################################################################################################################
	for (h in 1:nrAttr)
	{
		if((!is.na(eval(parse(text=paste(rAttrName[h],"_max_seeds[1]",sep=""))))) &&
(!is.na(eval(parse(text=paste(rAttrName[h],"_max_seeds[1,2]",sep=""))))))
{



			f=paste(rAttrName[h],"_max_seeds=",rAttrName[h],"_max_seeds[order(",rAttrName[h],"_max_seeds[,3],decreasing=T),] ",sep="")
			eval(parse(text=f))
			
			f=paste("fxGetCovariateGeneSet('Max",rAttrName[h],"',",rAttrName[h],"_max_seeds,initGene,num_row,gnlist,z,rAttrName,threshold_all,fm) ",sep="")
			eval(parse(text=f))
		}
		if((!is.na(eval(parse(text=paste(rAttrName[h],"_min_seeds[1]",sep=""))))) &&
(!is.na(eval(parse(text=paste(rAttrName[h],"_min_seeds[1,2]",sep=""))))))
{
			f=paste(rAttrName[h],"_min_seeds=",rAttrName[h],"_min_seeds[order(",rAttrName[h],"_min_seeds[,3],decreasing=F),] ",sep="")
			eval(parse(text=f))
			
			f=paste("fxGetCovariateGeneSet('Min",rAttrName[h],"',",rAttrName[h],"_min_seeds,initGene,num_row,gnlist,z,rAttrName,threshold_all_min,fm)",sep="")
			eval(parse(text=f))
		}
	}
	
	#########################################
	#PRINT OUT A PLOT OF THRESHOLD DENSITY
	#########################################
	if(length(AllCovariatesCoeff) > 0)
	{
		tdp=paste(out, "Threshold_Density_Plot.png", sep = "")
		png(tdp,width=1200, height=700)
		plot(density(AllCovariatesCoeff),main="Density Coeff Plot",xlab="coefficient")
		dev.off()
	}
	print("COMPLETED")

####################################################################### END MAIN #######################################################################################


