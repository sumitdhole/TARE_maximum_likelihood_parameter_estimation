# Each generation comprises of 12 genotypes: 
#   Wt/Wt, Wt/r, Wt/D, D/D, D/r, r/r x 6 male and 6 female. 
# Three drive parameters are: 
#   x - rate of wt/D -> D/D in germline
#   y - rate of wt/D -> r/D in germline
#   z - rate of wt -> r in any embryo where mom has at least one D allele

init.matrix <- function() {
  trans.matrix <- read.table("trans.matrix.txt", header=T)
  trans.matrix[,c(3:8)] <- trans.matrix[,c(3:8)]/4.0
  return(trans.matrix)
}
  
# germline conversion parent (drive conversion and resistance formation)
# geno.freq is an array of 6 genotype frequencies that add up to 1
germ <- function(geno.freq, x, y, z) {
  new4 = geno.freq[4] + geno.freq[3]*x 
  new5 = geno.freq[5] + geno.freq[3]*y
  new3 = geno.freq[3]*(1-x-y)
  return (c(geno.freq[1],geno.freq[2],new3,new4,new5,geno.freq[6]))
}

# germline conversion (drive conversion and resistance formation)
# geno.freq is an array of 6 genotype frequencies that add up to 1
germ <- function(geno.freq, x, y, z) {
  new4 = geno.freq[4] + geno.freq[3]*x 
  new5 = geno.freq[5] + geno.freq[3]*y
  new3 = geno.freq[3]*(1-x-y)
  return (c(geno.freq[1],geno.freq[2],new3,new4,new5,geno.freq[6]))
}

# embryo conversion (only resistance formation)
# geno.freq is an array of 6 genotype frequencies that add up to 1
embryo <- function(geno.freq, x, y, z) {
  new5 = geno.freq[5] + geno.freq[3]*z
  new6 = geno.freq[6] + geno.freq[1]*z^2 + geno.freq[2]*z ########nothing here that identifies only wt.wt or wt.r with D mother
  new2 = geno.freq[2]*(1-z) + 2*geno.freq[1]*z*(1-z) ###BUGGGG - the last term should be 2*genofreq[1]*z(1-z). I have fixed this in this file
  new1 = geno.freq[1]*(1-z)^2
  new3 = geno.freq[3]*(1-z)
  return (c(new1,new2,new3,geno.freq[4],new5,new6))
}

# generate transition matrix given parameter x, y, z
gen.matrix <- function(og.matrix, x, y, z) {
  # wt/D convert wt->D or wt->r in germline
  new.matrix <- og.matrix
  germ.drive <- germ(c(0,0,1,0,0,0),x,y,z)
  for(i in which((new.matrix$male-3)*(new.matrix$female-3)==0) ) {
      new.male <- rep(0,6)
      new.female <- rep(0,6)
      if(new.matrix[i,]$male==3&&new.matrix[i,]$female==3) {
        new.male <- germ.drive
        new.female <- germ.drive
      } else if(new.matrix[i,]$male==3) {
        new.male <- germ.drive
        new.female[new.matrix[i,]$female]<-1
      } else if(new.matrix[i,]$female==3) {  #########this else if is unnecessary. could just be else
        new.female <- germ.drive
        new.male[new.matrix[i,]$male]<-1
      }
      new.pair <- as.vector(outer(new.female,new.male))
      new.progeny <- new.pair*og.matrix[,3:8]
      new.matrix[i,3:8] <- colSums(new.progeny)/sum(new.progeny)
  }

  # */D female convert wt->r in embryo
  for (k in which(new.matrix$female%in%c(3,4,5))) {
    new.matrix[k,3:8] <- embryo(new.matrix[k,c(3:8)],x,y,z)
  }
  
  return(new.matrix)
}

# og.matrix <- init.matrix()
# new.matrix <- gen.matrix(og.matrix, 0.76, 0.06, 0.32)
 for(i in c(1:36)) {
   if(new.matrix[i,]$male==1){new.matrix[i,]$male<-'wt.wt'}
   if(new.matrix[i,]$male==2){new.matrix[i,]$male<-'wt.r'}
   if(new.matrix[i,]$male==3){new.matrix[i,]$male<-'wt.D'}
   if(new.matrix[i,]$male==4){new.matrix[i,]$male<-'D.D'}
   if(new.matrix[i,]$male==5){new.matrix[i,]$male<-'D.r'}
   if(new.matrix[i,]$male==6){new.matrix[i,]$male<-'r.r'}
   if(new.matrix[i,]$female==1){new.matrix[i,]$female<-'wt.wt'}
   if(new.matrix[i,]$female==2){new.matrix[i,]$female<-'wt.r'}
   if(new.matrix[i,]$female==3){new.matrix[i,]$female<-'wt.D'}
   if(new.matrix[i,]$female==4){new.matrix[i,]$female<-'D.D'}
   if(new.matrix[i,]$female==5){new.matrix[i,]$female<-'D.r'}
   if(new.matrix[i,]$female==6){new.matrix[i,]$female<-'r.r'}
 }
 View(new.matrix)
# write.table(new.matrix,'TRANSITION.txt',row.names = F)


######## prewritten fitness parameter logic #######

# Generate 6 mate choice fitness parameters for 6 male phenotypes
gen.m.codom.para <- function(m) {
  # assume female phenotypes don't matter in mate choice
  # assume m.Dr = m.Dwt = sqrt(m.DD)
  m.para <- rep(1,6)  
  m.para[4] = m
  # D codominant
  m.para[3] = sqrt(m)
  m.para[5] = sqrt(m)
  return(m.para)
}

# Generate 6 mate choice fitness parameters under DOMINANT model
gen.m.dom.para <- function(m) {
  # assume female phenotypes don't matter in mate choice
  # assume m.Dr = m.Dwt = m.DD
  m.para <- rep(1,6)  
  m.para[4] = m
  # D dominant
  m.para[3] = m
  m.para[5] = m
  return(m.para)
}

# Generate 36 mate choice fitness parameters for 36 mating pairs
# 36 parameters are 6 replicates of a set of 6 parameters since 
# male phenotypes don't matter in a mating pair
gen.f.codom.para <- function(f) {
  # assume male phenotypes don't matter in fecundity
  # assume f.Dr = f.Dwt = sqrt(f.DD)
  f.para <- rep(1,6)  
  f.para[4] = f
  # D codominant                      ############# Multiplicative fecundity cost of D to females
  f.para[3] = sqrt(f)
  f.para[5] = sqrt(f)
  return(rep(f.para,6))
}

# Generate 36 mate choice fitness parameters for 36 mating pairs
# under a DOMINANT model
gen.f.dom.para <- function(f) {
  # assume male phenotypes don't matter in fecundity
  # assume f.Dr = f.Dwt = f.DD
  f.para <- rep(1,6)  
  f.para[4] = f
  # D dominant
  f.para[3] = f
  f.para[5] = f
  return(rep(f.para,6))
}

# Generate 12 viablity parameters for 12 phenotypes (sex dependent)  ################ Does this mean 12 GENOtypes?
# under HAPLOLETHAL model
gen.v.haplo.para <- function(v) {
  # assume v.Dwt = sqrt(v.DD)
  # assume male and female have same viability
  # assume v.rr = v.wtr = v.Dr = 0
  v.para <- rep(1,6)
  v.para[6] = 0
  v.para[4] = v
  v.para[3] = sqrt(v)
  v.para[2] = 0
  v.para[5] = 0
  return(rep(v.para,2))
}

# Generate 12 viablity parameters for 12 phenotypes (sex dependent)
# under HAPLOLETHAL DOMINANT model
gen.v.haplo.dom.para <- function(v) {
  # assume v.Dwt = v.DD
  # assume male and female have same viability
  # assume v.rr = v.wtr = v.Dr = 0
  v.para <- rep(1,6)
  v.para[6] = 0
  v.para[4] = v
  v.para[3] = v
  v.para[2] = 0
  v.para[5] = 0
  return(rep(v.para,2))
}

# Generate 12 viablity parameters for 12 phenotypes (sex dependent)
# under TARE model
gen.v.tare.para <- function(v) {
  # assume v.Dr = v.Dwt = sqrt(v.DD)
  # assume male and female have same viability
  # assume v.rr = 0
  v.para <- rep(1,6)
  v.para[6] = 0
  v.para[4] = v
  # TARE codominant
  v.para[3] = sqrt(v)
  #v.para[5] = sqrt(v)
  # TARE dominant
  #v.para[3] = v
  #v.para[5] = v
  return(rep(v.para,2))
}

######## CUSTOMIZE fitness parameter logic here #######
# mate choice parameters
gen.m.para <- function(m){
  # return an array of 6 parameters
  return(rep(1,6))
}

# fecundity parameters
gen.f.para <- function(f){
  # return an array of 36 parameters
  return(rep(1,36))
}

# viability parameters
gen.v.para <- function(v){
  # return an array of 6 parameters
  return(rep(1,12))
}

# gn0, gn1 are two arrays of 12 genotype frequency
# trans.matrix is 36 x 6 (assuming autosomal locus)
# m.para is an array of length 6 **assume female phenotypes don't matter in mate choice
# f.para is an array of length 36
# v.para is an array of length 12
step.help <- function(gn0, trans.matrix, m.para, f.para, v.para) {
  male <- gn0[1:6]/sum(gn0[1:6])
  female <- gn0[7:12]/sum(gn0[7:12])
  
  # MATE CHOICE TAKES PLACE HERE
  male <- male*m.para/sum(male*m.para)          ########## Male freqs are weighted by mating siccess
  pairs <- as.vector(outer(female, male))
  matrix0 <- trans.matrix[,3:8]
  matrix1 <- pairs*matrix0
  
  # FECUNDITY SELECTION TAKES PLACE HERE
  matrix1 <- f.para*matrix1
  
  # VIABILITY SELECTION TAKES PLACE HERE
  matrix1 <- cbind(matrix1, matrix1) #assuming autosomal locus
  matrix1 <- t(v.para*t(matrix1))
  gn1 <- colSums(matrix1)/sum(matrix1)
  return(gn1)
}

# given the frequency for parent generation, predict progeny frequency
step <- function(gn0, x,y,z, m, f, v, para.mode){
  og.matrix <- init.matrix()
  trans.matrix <- gen.matrix(og.matrix, x, y, z)
  if(para.mode=='codom.haplo'){
    m.para <- gen.m.codom.para(m)
    f.para <- gen.f.codom.para(f)
    v.para <- gen.v.haplo.para(v)
  } else if(para.mode=='dom.haplo'){
    m.para <- gen.m.dom.para(m)
    f.para <- gen.f.dom.para(f)
    v.para <- gen.v.haplo.dom.para(v)
  } else if(para.mode=='dom.tare'){
    m.para <- gen.m.dom.para(m)
    f.para <- gen.f.dom.para(f)
    v.para <- gen.v.tare.dom.para(v)
  } else if(para.mode=='new'){
    m.para <- gen.m.para(m)
    f.para <- gen.f.para(f)
    v.para <- gen.v.para(v)
  }
  gn1 <- step.help(gn0, trans.matrix, m.para, f.para, v.para)
  return(gn1)
}

# plot frequency trajectory 
plot.freq <- function(freq) {
  generation <- nrow(freq)-1
  col.1 <- 'black'
  col.2 <- 'dark red'
  col.3 <- 'red'
  
  x.range <- c(0:generation)
  plot(x.range,freq[,1], ylim=c(0,1), xlim=c(0,generation+1), type="o", col=col.1, pch=17, axes=F)
  points(x.range,freq[,2], type="o",col=col.1, pch=19, lty=3)
  points(x.range,freq[,3], type="o",col=col.2, pch=17, lty=3)
  points(x.range,freq[,4], type="o",col=col.3, pch=17)
  points(x.range,freq[,5], type="o",col=col.2, pch=19, lty=3)
  points(x.range,freq[,6], type="o",col=col.1, pch=19)
  
  axis(1,padj=-0.5)
  axis(2,padj=0.5)
  
  legend(0.58+0.33,y=0.99, c('wt/wt','wt/r','r/r'), pch=c(17,19,19), lty=c(1,3,1), col=c(col.1,col.1,col.1), bty="n", y.intersp=1.04)
  legend(generation/2+0.08-0.33,y=0.99, c('D/wt','D/r','D/D'), pch=c(17,19,17), lty=c(3,3,1), col=c(col.2,col.2,col.3), bty="n", y.intersp=1.04)
}

# plot phenotype (red/wt) frequency trajectory    ############## Expected phenotype freq.
plot.pheno.freq <- function(freq,real='') {
  generation <- nrow(freq)-1
  col.1 <- 'black'
  col.2 <- 'red'
  col.3 <- 'green'
  
  x.range <- c(0:generation)
  if(length(freq[1,])==12) {
    plot(x.range,freq[,3]+freq[,4]+freq[,5]+freq[,9]+freq[,10]+freq[,11], ylim=c(0,1), xlim=c(0,generation+1), type="o", col=col.2, pch=17, axes=F)
    points(x.range,freq[,1]+freq[,2]+freq[,6]+freq[,7]+freq[,8]+freq[,12], type="o",col=col.1, pch=17, lty=3)
  } else if(length(freq[1,])==2) {
    plot(x.range,freq[,1], ylim=c(0,1), xlim=c(0,generation+1), type="o", col=col.2, pch=17, axes=F)
    points(x.range,freq[,2], type="o",col=col.1, pch=17, lty=3)
  } else{
    print("WRONG FREQUENCY LIST")
  }
  if(real=='haplo') {
    X.haplo = read.delim2("haplo.txt",sep="")
    points(c(1:nrow(X.haplo))-1,X.haplo[,1]/(X.haplo[,1]+X.haplo[,2]), type="o",col=col.3, pch=17, lty=3)
  } else if(real=='tare1') {
    X = read.delim2("TARE.txt",sep="")
    points(c(1:nrow(X))-1,X[,1]/(X[,1]+X[,2]), type="o",col=col.3, pch=17, lty=3)
  } else if(real=='tare2') {
    X = read.delim2("TARE2.txt",sep="")
    points(c(1:nrow(X))-1,X[,1]/(X[,1]+X[,2]), type="o",col=col.3, pch=17, lty=3)
  }
  
  axis(1,padj=-0.5)
  axis(2,padj=0.5)
  
  #legend(0.58+0.33,y=0.99, c('wt/wt','wt/r','r/r'), pch=c(17,19,19), lty=c(1,3,1), col=c(col.1,col.1,col.1), bty="n", y.intersp=1.04)
  #legend(generation/2+0.08-0.33,y=0.99, c('D/wt','D/r','D/D'), pch=c(17,19,17), lty=c(3,3,1), col=c(col.2,col.2,col.3), bty="n", y.intersp=1.04)
}

## Simulate genotype trajectories for given parameters and vector x of starting 
#  frequencies over g generations. If n.e=0, deterministic trajectories are simulated. 
evolve <- function(start,n.e,m,f,v,x,y,z,g,para.mode) {
  freq <- rbind(start/sum(start),NULL)
  for(i in (1:g)) {
    new.gen <- step(freq[i,], x,y,z, m, f, v,para,mode)
    if(n.e > 0) {
      X.m <- t(rmultinom(1,n.e/2,new.gen[1:6]))
      X.m <- X.m/sum(X.m)
      X.f <- t(rmultinom(1,n.e/2,new.gen[7:12]))
      X.f <- X.f/sum(X.f)
    } else {
      X.m <- new.gen[1:6]
      X.f <- new.gen[7:12]
    }
    X <- c(X.m,X.f)
    freq <- rbind(freq,X,deparse.level=0)
  }
  return(freq)
}

#simulation test
# start <- rep(c(0.695, 0, 0, 0.305, 0, 0),2)
# 
# m <- 2
# f <- 1
# v <- 1
# x <- 0
# y <- 1
# z <- 1
# 
# g <- 6
# n.e <- 0
# 
# sim.freq <- evolve(start,n.e,m,f,v,x,y,z,g,'codom.haplo')
# sim.pheno.freq <- data.frame()
# for(i in c(1:nrow(sim.freq))) {
#   pheno <- twelve.to.two(sim.freq[i,])
#   print(pheno)
#   sim.pheno.freq <- rbind.data.frame(sim.pheno.freq,pheno)
# }
# plot.pheno.freq(sim.pheno.freq,'tare1')

# Collapse 12 genotype frequency to 2 phenotype frequency
twelve.to.two <- function(g) {
  if (length(g)!=12) {stop("wrong freq length")}
  red <- g[3]+g[4]+g[5]+g[9]+g[10]+g[11]
  wt  <- g[1]+g[2]+g[6]+g[7]+g[8]+g[12]
  return(c(red,wt)/(red+wt))
}

# Expend 2 phenotype frequency to 12 genotype frequency based on predicted 
# genotype frequency
two.to.twelve <- function(g12,p2) {
  p2 <- p2/sum(p2)
  g2 <- twelve.to.two(g12)
  for(gi in c(3,4,5,9,10,11)) {g12[gi]=g12[gi]/g2[1]*p2[,1]} #######If g2[1] is ever 0, this would give an infinite value. Should include a condition for checking it's positive. If it is 0, then g12[gi] will be 0 to begin with and need not be changed.
  for(gi in c(1,2,6,7,8,12)) {g12[gi]=g12[gi]/g2[2]*p2[,2]}
  return(g12)
}

## Functions for calculating logL of parameter set and one generation transition.
# both pred.freq and real.freq have two numbers: red & wt
rho <- function(pred.freq,real.freq,n) {
  if (length(pred.freq)!=2||length(real.freq)!=2) {stop("wrong freq length")}
  pred.freq <- pred.freq/sum(pred.freq)
  real.freq <- real.freq/sum(real.freq)
  real.num <- real.freq*n
  
  log.a <- log(n)-log(1-(pred.freq[1]^n+pred.freq[2]^n)/2)  ########### This is log(eq.7) in Liu, Champer et al 2019
  r <- log.a+lgamma(n+1)-lgamma(real.num[1]+1)-lgamma(real.num[2]+1) ########### This (+few lines below) is log(eq.5) in Liu, Champer et al 2019
  
  if(pred.freq[1]>0){
    r <- r+real.num[1]*log(pred.freq[1])
  }
  if(pred.freq[2]>0){
    r <- r+real.num[2]*log(pred.freq[2])
  }
  return(r)
}

## logL generator
logL.gen <- function(x,y,z) {
  ## Functions for calculating logL of parameter set and phenotype trajectory.
  logL <- function(para,X.list,mode,para.mode) {
    l <- 0
    # para: optimization parameters
    # X: table with observed phenotype numbers from single cage
    m <- 1
    f <- 1
    v <- 1
    ne <- 1
    
    # full model
    if(mode==1) {
      m <-  para[1]
      f <-  para[2]
      v <-  para[3]
      ne <- para[4]
    }
    
    # test all parameter model
    if(mode==2) {
      m <-  para[1]
      f <-  para[2]
      v <-  para[3]
      ne <- para[4]
      x <- para[5]
      y <- para[6]
      z <- para[7]
    }
    
    # mate choice model w/o ne
    if(mode==3) {
      m <-  para[1]
    }
    
    # fecundity model w/o ne
    if(mode==4) {
      f <-  para[1]
    }
    
    # viability model w/o ne
    if(mode==5) {
      v <-  para[1]
    }
    
    # full model w/o ne
    if(mode==6) {
      m <-  para[1]
      f <-  para[2]
      v <-  para[3]
    }
    
    # mate choice model
    if(mode==7) {
      m <-  para[1]
      ne <- para[2]
    }
    
    # fecundity model
    if(mode==8) {
      f <-  para[1]
      ne <- para[2]
    }
    
    # viability model
    if(mode==9) {
      v <-  para[1]
      ne <- para[2]
    }
    
    # mate choice=fecundity model
    if(mode==10) {
      m <-  para[1]
      f <-  para[1]
      ne <- para[2]
    }
    
    # mate choice=fecundity=viability model
    if(mode==11) {
      m <-  para[1]
      f <-  para[1]
      v <-  para[1]
      ne <- para[2]
    }
    
    for(cage in c(1:length(X.list))) {  ########## X.list creates a list of datasets 'X.haplo' for different cages
      
      X <- as.data.frame(X.list[cage])  ######## X converts an element of X.list (data of a single cage) into useable data
      
      parent.g <- rep(0,12)   ###492-497 starting generation of genotypes in a cage
      init <- X[1,]/sum(X[1,])
      parent.g[4] <- init[,1]/2
      parent.g[10] <- init[,1]/2
      parent.g[1] <- init[,2]/2
      parent.g[7] <- init[,2]/2
      
      for(t in (2:nrow(X))) {
        # parent: observed phenotype frequencies in generation t-1
        parent.p  <- X[t-1,]/sum(X[t-1,])
        
        # progeny: observed phenotype frequencies in generation t
        progeny.p <- X[t,]/sum(X[t,])
        
        if (ne.mode=='avg') {     ######## ne.mode='avg' uses the mean of 2 subsequent gens as the ne
          Ne <- ne*(sum(X[t-1,])+sum(X[t,]))/2
        } else if (ne.mode=='abs') {
          Ne <- ne
        } else{
          print('INVALID NE MODE')
        }
        if(mode%in%c(6,5,4,3)) {Ne <- 0.05*sum(X[t-1,])}  ######## not sure how 0.05 was picked
        
        # expected genotype frequencies in generation t
        progeny.g.pred <- step(parent.g, x,y,z,m,f,v,para.mode)
        progeny.p.pred <- twelve.to.two(progeny.g.pred)
        parent.g <- two.to.twelve(progeny.g.pred, progeny.p) ###predicted parental genotypes in generations after the first
        
        l <- l + rho(progeny.p.pred, progeny.p, Ne)
      } 
    }
    return(l)  
  }
  
  return(logL)
}

# CI generator
ci <- function(logL,X.list,df,para,mode,para.mode){
  #dL <- c(3.841,5.991,7.815,9.488,11.07,12.59,14.07)
  dL <- rep(3.841,7)      ####### chi sq value that needs to be exceeded for sign diff for dof=2 & p<0.05
  Delta <- 0.001
  
  para.left <- para
  para.right <- para
  high.logL <- logL(para,X.list,mode,para.mode)
  for(i in c(1:df)){
    new.para <- para
    # start delta
    delta <- para[i]/2
    while(T){
      para.left[i] <- max(0.001,para.left[i]-delta) # avoid 0
      new.para[i] <- para.left[i]
      new.logL <- logL(new.para,X.list,mode,para.mode)
      if(2*(high.logL-new.logL)>=dL[df]){
        if(2*delta<Delta) {break}
        else {
          para.left[i] <- para.left[i]+delta
          delta <- delta/2
        }
      }
    }
    new.para <- para
    # start delta
    delta <- para[i]/2
    while(T){
      para.right[i] <- para.right[i]+delta
      new.para[i] <- para.right[i]
      new.logL <- logL(new.para,X.list,mode,para.mode)
      if(2*(high.logL-new.logL)>=dL[df]){
        if(delta<Delta) {break}
        else {
          para.right[i] <- para.right[i]-delta
          delta <- delta/2
        }
      }
    }
  }
  true.interval <- data.frame(para)
  true.interval <- cbind.data.frame(true.interval,para.left)
  true.interval <- cbind.data.frame(true.interval,para.right)
  return(true.interval)
}

# parameter model mode number to text name
mode.to.name <- function(mode) {
  if(mode==7){return('m')}
  if(mode==8){return('f')}
  if(mode==9){return('v')}
  if(mode==10){return('m=f')}
  if(mode==11){return('m=f=v')}
}

run <- function(modes,para.modes,ne.mode,para.lower,para.upper,para.start,filename) {
  all.results <- data.frame()
  for(i in c(1:length(x))) {
    set.string <- paste0('x:',x[i],' y:',y[i],' z:',z[i])
    print(set.string)
    for(test.para.mode in para.modes) {
      for(test.mode in modes) {
        print(mode.to.name(test.mode))
        logL <- logL.gen(x[i],y[i],z[i])
        MLE.all3  <- optim(para.start,logL, method="L-BFGS-B",control=list(fnscale=-1),
                           lower=para.lower, upper=para.upper, mode=test.mode, X.list=list(X.haplo),
                           para.mode=test.para.mode, hessian=T)
        print(MLE.all3$par)
        true.interval <- ci(logL,list(X.haplo),length(para.start),MLE.all3$par,test.mode,test.para.mode)
        l <- logL(MLE.all3$par,list(X.haplo),mode=test.mode,test.para.mode)
        result <- data.frame(c(test.para.mode,l,mode.to.name(test.mode),true.interval[1,],true.interval[2,],x[i],y[i],z[i]))
        colnames(result) <- c('mode','logL','para','para.value','p.l','p.r','ne','ne.l','ne.r','x','y','z')
        all.results <- rbind.data.frame(all.results,result)
      }
    }
  }
  write.csv(all.results,filename,row.names = F)
}