################################################
##
##NAME: ASE.analysis_auto_V1.R
##AUTHOR: James Koltes
##DATE: March 4, 2015
##EXPECTED DATA FORMAT:
##column1: chr_position
##columns2+3 = reference_allele_count & alternate_allele_count
##after the first column, it is expected that each animal has a pair of alleles for the proceeding 2 columns
##OUTPUT FILE FORMAT:
## 6 columns, including: SNP.ID, alt.allele.prob, CI.95.upper, CI.95.lower, pvalue, qvalue)
##TODO: Need to add annotation data to the results and format data as GFF3 or VCF to create data tracks for viz
##POTENTIAL PROBLEMS: If qvalue library is not installed, this code will crash prior to writing the results to file!
##
##COMMAND LINE USAGE: R CMD BATCH ASE.analysis_auto_V1.R &
##
################################################
options(warn=-1) #turn warning messages off
commandLine = "TRUE"

##Instantiate variables
args <- commandArgs(TRUE)
outFileName = args[4]
dataDirectory = args[1]
outputFileDirectory = args[3]
infile = args[2]
minorAllele.obj = 0
minorAllele.ind	= 0 #This is the index of SNPs that have too low freq, need to be removed.
errors <- 0
outfile = data.frame(SNP.ID=NA, alt.allele.prob=NA, CI.95.upper=NA, CI.95.lower=NA , pvalue=NA)
snpName = NULL
ds.counter = NULL
#Initiate variables- rep for total number SNPs later
props <- NA
## Initialize pvals.props
pvals.props <- NA
ci.ul <- NA
ci.ll <- NA

#Set SNP filtering paramters
percentZerosTollerated = 0.2 #i.e. allow up to 20% of the individuals to have 0 counts for the alternate and reference alleles.
#require at least 50X coverage on average across SNPs
depthFilter = 50
#Require minor allele to be at least 0.20
filter.minorAlleleFreq = 0.2

##Read data and parse important info
setwd(dataDirectory)
dat = read.table(file=infile,header=TRUE,sep="",na.strings = "NA")
row.names(dat) = dat[,1] #this assumes the first column is the SNP ID
bcknames = dat[,1]
dat = dat[,-1] #dat ONLY contains allele counts from here on out!
#This arg defines the first column with allele count data
dat.col1 = 1
#Pick out the allele counts
ref.alleles.all  = dat[,(rep(seq(dat.col1,ncol(dat),by=2),times=1))]
alt.alleles.all = dat[,(rep(seq(dat.col1+1,ncol(dat),by=2),times=1))]
aid.loc = as.factor(seq(1:(ncol(dat)/2)))

##ERROR HANDLING FUNCTION
tryCatch.W.E <- function(expr)
 {
     W <- NULL
     # warning handler
     w.handler <- function(w){ 
 	W <<- w
 	invokeRestart("muffleWarning")
     }
     #normally next line warning = warning = w.handler
     #list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
 	 list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
 	 #print(paste("Error at position",i,sep=" "))
 	 #print()
  }
  
#create is.empty() function
  is.empty <- function(x, mode=NULL){
     if (is.null(mode)) mode <- class(x)
     identical(vector(mode,1),c(x,vector(class(x),1)))
 }


########### Start the data filtering process ########### 
##First- filter on # zeros (zero Filter)
zeroCutoff = floor(percentZerosTollerated*(ncol(dat)/2)) #specifies the max number of zeros allowed across all samples (ncol/2 = # samples since two alleles per sample), floor rounds down to be more strict on # zeros allowed

##Zero filter: limit the # of zero counts
#First, filter the reference allele counts
numZeros <- apply(ref.alleles.all == 0,1,sum) #get the total number of zeros in a row here
nzIndex <- which(numZeros>zeroCutoff) #get index of rows with more zeros than the filter here
filter1.ref <- ref.alleles.all[-nzIndex,]
filter1.alt <- alt.alleles.all[-nzIndex,]
##Now filter the alternate allele counts
numZeros <- apply(filter1.alt ==0,1,sum)
nzIndex.alt <- which(numZeros>zeroCutoff)
##Next two lines contain the zero filtered allele counts for clinical animals
ref.alleles.all  <- filter1.ref[-nzIndex.alt,]
alt.alleles.all <- filter1.alt[-nzIndex.alt,]
sel.snp.clin = row.names(ref.alleles.all) #these are just row names- not sure if still needed
ref.tot = as.numeric(rowSums(ref.alleles.all))

##SNP Filtering Strategy one: Number of counts allowed to have zeros
ind3 = which((as.numeric(ref.tot)>=depthFilter)=="TRUE")
ref.allele1.1 = ref.alleles.all[ind3, ]
alt.allele1.1 = alt.alleles.all[ind3, ]
alt.tot = as.numeric(rowSums(alt.allele1.1))
ind4 = which((as.numeric(alt.tot)>=depthFilter)=="TRUE")
ref.alleles.all = ref.allele1.1[ind4, ]
alt.alleles.all = alt.allele1.1[ind4, ]

##SNP Filtering Strategy two: Read Depth filtering
##Sum reference and alt alleles to get total count per sample (i.e. read depth at a SNP)
SNPdepth = ref.alleles.all + alt.alleles.all
SNPdepthAVE = as.numeric( rowSums(SNPdepth)/ncol(SNPdepth) )
ind3 = which((as.numeric(SNPdepthAVE)>=depthFilter)=="TRUE")
ref.allele1.1 = ref.alleles.all[ind3, ]
alt.allele1.1 = alt.alleles.all[ind3, ]
alt.tot = as.numeric(rowSums(alt.allele1.1))
ind4 = which((as.numeric(alt.tot)>=depthFilter)=="TRUE")
ref.alleles.cl = ref.allele1.1[ind4, ]
alt.alleles.cl = alt.allele1.1[ind4, ]


##SNP Filtering Strategy three: minor allele frequency filter
counter3 = 1
for (i in 1:dim(alt.alleles.cl)[1])
{
	rowCts.alt = alt.alleles.cl[i, ]
	rowCts.ref = ref.alleles.cl[i, ]
	for(j in 1:dim(rowCts.alt)[2])
	{
		#sum each allele pair
		total = sum( (rowCts.alt[ ,j]), (rowCts.ref[ ,j]) )
		#get the minor allele freq
		if(total != 0)
		{
			alt.freq = (rowCts.alt[ ,j])/total
			ref.freq = (rowCts.ref[ ,j])/total
		} else {alt.freq=0; ref.freq=0}
		#save the minor allele freqs across row
		if( (alt.freq < filter.minorAlleleFreq ) || (alt.freq > (1-filter.minorAlleleFreq) ) )
		{
			minorAllele.obj[j] =  0
		} else{minorAllele.obj[j] =  1}
		if( j==(dim(rowCts.alt)[2]) )
		{
			#Count # alleles with freq less than the filter
			cl.minor = as.matrix(minorAllele.obj[1:dim(rowCts.alt)[2]]) #I think dim(rowCts.alt)[2] is ncol?
			numZeros1 <- apply(cl.minor == 0,2,sum) #2 indicates by column
			#require 4 of 5 SNPs to meet the filter, thus each has 1
			if( numZeros1 > 1 )
				{
				minorAllele.ind[counter3] = i
				counter3 = counter3 + 1
				}
		} #END: if(j==(dim(rowCts.alt)[2]))
	} #END: for(j in 1:dim(rowCts.alt)[2])	
} #END: for (i in 1:dim(alt.alleles.all)[1])
##remove the rows that do not meet the allele frequency filter
ref.alleles.cl = ref.alleles.cl[-minorAllele.ind, ]
alt.alleles.cl = alt.alleles.cl[-minorAllele.ind, ]


########### Run the ASE analysis ###########
snp.ID = row.names(alt.alleles.cl)
require(lme4)
for(m in 1:length(snp.ID))
{
	##grab counts from full data
	cnt.ref = as.integer(ref.alleles.cl[m,])
	cnt.alt = as.integer(alt.alleles.cl[m,])
	##Test for errors to catch them to keep the loop running to analyze all SNPs
	error.test <- tryCatch.W.E(glmer(cbind(cnt.alt, cnt.ref) ~ (1 | aid.loc), family = binomial, nAGQ=3))
	if( (isS4(error.test$value)) || (is.null(grepRaw("Downdated",error.test$value))) && (is.null(error.test$warning)) && (is.null(grepRaw("maxstephalfit",error.test$value)))  )
		{ 
			glm.ase <- glmer(cbind(cnt.alt, cnt.ref) ~  (1 | aid.loc), family = binomial, nAGQ=3)
			fe <- fixef(glm.ase)
			num.fe <- length(fe)
			error.test2 <- tryCatch.W.E(vcov(glm.ase))
			if((is.empty(grepRaw("Lapackk",error.test2$value))))
			{
				vc <- vcov(glm.ase)
				D <- diag(num.fe)
				D[1, ] <- 1
				lmean <- fe[c("(Intercept)")] #with fe- need the names of the effects for sub, non, clinc 							from the glm above 
				props[m] <- plogis(lmean)
				## Use the column of D for the delta method
  				d <- t(D)
  				se <- diag(sqrt(d %*% vc %*% t(d)))
  				## Get the p-value by comparing the lmean to a normal with our
 				## standard error. It's a two-tailed test so we multiply times two.
  				pvals.props[m] <- 2 * pnorm(abs(lmean), mean=0, sd=se, lower.tail=FALSE)  
  				ci.ul[m] <- plogis(lmean + 1.96 * se)
  				ci.ll[m] <- plogis(lmean - 1.96 * se)
				snpName[m] = snp.ID[m]
			} #END: if((is.null(grepRaw("Lapackk",error.test2$value))))
		} #END: if( (isS4(error.test$value)) || (is.null(grepRaw("Downdated",error.test$value))) )
		if(commandLine == "FALSE")
		{
			if(m==1){print(paste("Total # Tests to run: ",length(snp.ID)))}
			if(m%%50 == 0){print(paste("Currently on SNP test #: ", m))}
		}
} #END: for(m in 1:length(snp.ID))

##Create outfile data object
outfile = data.frame(SNP.ID=snpName, alt.allele.prob=props, CI.95.upper=ci.ul, CI.95.lower=ci.ll , pvalue=pvals.props)

##Determine qvalues
library(qvalue)
p = outfile$pvalue
qobj = qvalue(p)
outfile$qvalue = qobj$qvalues

#Write results to file
setwd(outputFileDirectory)
write.table(outfile, file=outFileName, row.names=FALSE, quote=FALSE,sep="\t")
