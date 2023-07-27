toep.to.circ <- function(toep,n1,n2){
    g <- c(t( (row(toep) == 1 ) * row(toep) )) * rep(rep(1:n2,each=n1),n1*n2)
    g <- factor(g)
    row.split <- split(toep,g)[-1]
    tmp.fnc <- function(x) {return(c(x,0,x[length(x):2]))}
    circ.first.row <- lapply(row.split,tmp.fnc)
    circ.first.row <-  c(circ.first.row, list(rep(0,2*n1)),
            circ.first.row[length(circ.first.row):2])
    return(circ.first.row)
}

get.circ.mat <- function(circ.row){
    len <- length(circ.row)
    n = length(circ.row[[1]])
    make.mat <- function(mat.row){
        m <- length(mat.row)
        x <- matrix(NA,nrow=m,ncol=m)
        x[1,] <- mat.row
        for(i in 2:m){
            x[i,] <- mat.row[c((m-i+2):m,1:(m-i+1))]
        }
        return(x)
    }
    circ.mat.list <- lapply(circ.row,make.mat)
    circ.mat <- do.call(cbind,circ.mat.list)
    for(i in 2:len){
        circ.mat <- rbind(circ.mat,do.call(cbind,circ.mat.list[c((len-i+2):len,1:(len-i+1))]))
    }
   return(circ.mat)
}

library(mvtnorm)
D = 5
loc <- as.matrix(expand.grid((1:D)/(D+1),(1:D)/(D+1) ))
colnames(loc) <- c("x", "y")
loc_dist <- as.matrix(dist(loc,method="euclidean",upper=F,diag=F))
cov_mat <-  exp(- (loc_dist/0.5)^2) + diag(1e-6,nrow=D^2,ncol=D^2)
circ_first_row = toep.to.circ(cov_mat,D,D)
circ_cov_mat = get.circ.mat(circ_first_row)
isSymmetric(circ_cov_mat)
eignvalues <- eigen(circ_cov_mat)
eignvalues$values

set.seed(10324)
x = rnorm(D^2,mean=0,sd=1)
log_d = dmvnorm(x=x,mean=rep(0,D^2),sigma=cov_mat,log=TRUE)
x_circ = rnorm(D^2*4,0,1)
log_d = dmvnorm(x_circ,mean=rep(0,D^2*4),sigma=circ_cov_mat,log=TRUE)
