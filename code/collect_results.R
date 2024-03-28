load("data/precip.RData")
marginal.fit.files <- list.files("data/",pattern = "marginal_fit_.*.RData",full.names = TRUE)

for(i in 1:length(marginal.fit.files)){
    load(marginal.fit.files[i],e<-new.env())
    print(is.null(e$data.df.gpd))
}

### collect results for bootstrap ### 
library(stringr)
boot.files <- list.files(path="data/bootstrap/",pattern="fit_bootstrap_",full.names = TRUE)
region.idx <- as.numeric(str_extract(str_extract(boot.files,"_[1-8]\\."),"[1-8]"))
boot.idx <- as.numeric(str_extract(str_extract(boot.files,"_\\d+_"),"\\d+"))

results.boot.list <- list()
for(i in 1:length(boot.files)){
    load(boot.files[i],e<-new.env())
    results.boot.list[[i]] <- e$result.list
}

idx.valid = !unlist(lapply(results.boot.list,function(x){any(!unlist(lapply(x,function(y){length(y)==4})))}))

system(paste("rm",boot.files[!idx.valid]))
results.boot.list <- results.boot.list[idx.valid]
boot.files <- boot.files[idx.valid]
region.idx <- region.idx[idx.valid]
boot.idx <- boot.idx[idx.valid]
boot.collect <- function(data,idx,norm.idx,season.idx,region.idx,boot.idx){
    estimates.boot = matrix(unlist(lapply(data[region.idx==idx & boot.idx != 301],function(x){y=x[[norm.idx]][[season.idx]]$par;y[1] <- 2/(1+exp(-y[1]));y})),ncol=3,byrow=TRUE)
    est.true = results.boot.list[region.idx == idx  & boot.idx == 301][[1]][[norm.idx]][[season.idx]]$par; est.true[1] <- 2/(1+exp(-est.true[1]))
    n = nrow(estimates.boot)
    sd.jackknife = apply(estimates.boot,2,function(x){return(sd(x))})
    
    return(list(sd=sd.jackknife,true=est.true,jack=estimates.boot))
}

boot.result.list <- list(list(),list())

for(norm.idx in 1:2){
    for(season.idx in 1:4){
        result.list <- list()
        for(idx in 1:8){
            result.list[[idx]] <- boot.collect(data=results.boot.list,idx=idx,norm.idx=norm.idx,season.idx=season.idx,region.idx=region.idx,boot.idx=boot.idx)
        }
        boot.result.list[[norm.idx]][[season.idx]] <- result.list
    }
}

save(boot.result.list, results.boot.list, region.idx, boot.files, boot.idx,file="data/dep.fit.boot.results.RData")



