### Founder and parents 

#1a get genetic map 
Gdat_c=read.csv("Genetic map_for_Simulation_editedto0_sept17.csv",header = T)

Gdat_c$MarkerID=as.character(Gdat_c$MarkerID) ### Make sure that MarkerID column is character

#2a get genetic data from vcf file using Gformat package 
library(Gformat)
Gdat0=read.vcf("Genotype_Simu.vcf") ##### read vcf file using readvcf function from Gformat (kaler's package)
colnames(Gdat0)=Gdat0[1,] ### set up column names 
Gdat1=Gdat0[-1,] ### removed 1st row as we dont need that 
Gdat1$chr=as.numeric(Gdat1$`#CHROM`) ### make sure chr are numeric 
Gdat2=Gdat1[,c(2100,10:2099)]  ### subset the columns that we need. Remeber this is may vary depending upon  the number of genotype 
rownames(Gdat2)=Gdat1$ID ### set rownames to marker ID


Gdat3=Gdat2[Gdat_c$MarkerID,] ## Sort the Gdat2 (Genotype data) according to marker order of Gdat_c (genetic map) 
### lets make sure marker order matches in two file 
sum(!(rownames(Gdat3)==as.character(Gdat_c$MarkerID))) ### this should give 0 



###### Now lets use both files to make genemap and haplotype file 
### Here we are using reduced datsets to reduce the memory usase. 


#1c for making genemap file, you need chromosome and position id 
### genMapc0=Gdat_c[c(6,3)]
genMapc0=as.data.frame(Gdat_c[c("chrom","Genetic.Position")]) ### Lets subset the file to 2 columns that we need 
genMapc1=split(genMapc0[,"Genetic.Position"],genMapc0$chrom) ### Lets split the position by chromosome the file to 2 columns that we 

#1d for making haplotype file, you need chromosome, SNP info and genotypes. We will run on subset of genotypes 
## for this we will use Gdat2 file from above 

haplotypes=as.data.frame(Gdat3)
haplotypes0=haplotypes[,-1] ## romove first row 

### create function `foo`` to pick first (before /) or second number (after /) from each vcf column. 
# here , x is genetic data file and g can be 1 (pick first value ie before /) or 2 (pick second value ie after /)
foo=function(x,g){
  sp=strsplit(x,"/")
  as.numeric(unlist(lapply(sp, `[[`, g)))
  
}

haplotypes0.1=apply(haplotypes0, 2,function(x)foo(x,1)) ### haplotype 1
haplotypes0.2=apply(haplotypes0, 2,function(x)foo(x,2)) ### haplotype 2
haplotypes.a=as.data.frame(cbind(haplotypes0.1,haplotypes0.2))
haplotypes.a1=haplotypes.a[,order(colnames(haplotypes.a))] ### This will sort and put bith haps next to each other 
rownames(haplotypes.a1)=rownames(haplotypes0) ### setting row names as marker ID

#### next we split the haplotype by chromosome and transponse, this is required format for this package to make founderpop
haplotypes1=split(haplotypes.a1,haplotypes$chr)
haplotypes2=lapply(haplotypes1, function(x){
  y=t(x)
  y
})


##1 e maker founder pop
founderPop1 = newMapPop(genMap=genMapc1,haplotypes=haplotypes2,ploidy = 2,inbred = F)

