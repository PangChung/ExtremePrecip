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

save(results.boot.list,region.idx, boot.files, boot.idx,file="data/dep.fit.boot.results.RData")


load("data/dep.fit.boot.results.RData")

jackknife <- function(data,idx,norm.idx,season.idx,region.idx,boot.idx){
    estimates.boot = matrix(unlist(lapply(data[region.idx==idx & boot.idx != 301],function(x){x[[norm.idx]][[season.idx]]$par})),ncol=3,byrow=TRUE)
    est.true = results.boot.list[region.idx == idx  & boot.idx == 301][[1]][[norm.idx]][[season.idx]]$par
    n = nrow(estimates.boot)
    sd.jackknife = apply(estimates.boot,2,function(x){return(sd(x)*(n-1)/sqrt(n))})
    return(list(sd=sd.jackknife,true=est.true,jack=estimates.boot))
}

jack.result.list <- list(list(),list())
for(norm.idx in 1:2){
    for(season.idx in 1:4){
        result.list <- list()
        for(idx in 1:8){
            result.list[[idx]] <- jackknife(data=results.boot.list,idx=idx,norm.idx=norm.idx,season.idx=season.idx,region.idx=region.idx,boot.idx=boot.idx)
        }
        jack.result.list[[norm.idx]][[season.idx]] <- result.list
    }
}
str(jack.result.list)
