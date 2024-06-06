# Set working directory to where both the current file and 
# 'geneDrive.power.source.R' are located
setwd("")

source('TARE_source.R')

# Set any set of parameters, each column represents a set
#   x - rate of wt/D -> D/D in germline
#   y - rate of wt/D -> r/D in germline
#   z - rate of wt -> r in any embryo where mom has at least one D allele
x <- c(0.78,0.8)
y <- c(0.1,0.1)
z <- c(0.31,0.31)

# Set the name of the phenotype file you wish to read from
# See "PhenoEx.txt" for accepted format
X.haplo <- read.delim2("PhenoExample.txt",sep="")

# model
# mode==7: m only
# mode==8: f only
# mode==9: v only
# mode==10: m=f only
# mode==11: m=f=v
modes <- c(10)

# effective population size. 
# avg=average of current and previous generation
# abs=an absolute value
ne.mode <- 'avg'

# prewritten parameter modes
# 'new' if customized logic is implemented
para.modes <- c('codom.haplo')

# lower bound of parameter 1,2,3...
para.lower <- c(0.01,0.001)

# upper bound of parameter 1,2,3...
para.upper <- c(10,1.0)    ################### What is "10" for?

# start point of parameter 1,2,3...
para.start <- c(0.5,0.05)

# file name of the output result file
filename <- 'test.csv'

run(modes,para.modes,ne.mode,para.lower,para.upper,para.start,filename)
