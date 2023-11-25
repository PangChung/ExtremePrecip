load("data/precip.RData")
marginal.fit.files <- list.files("data/",pattern = "marginal_fit_.*.RData",full.names = TRUE)
dep.fit.files <- list.files("data/",pattern="fit_pot_ST_*",full.names = TRUE)

for(i in 1:length(marginal.fit.files)){
    load(marginal.fit.files[i],e<-new.env())
    print(is.null(e$data.df.gpd))
}


results.list <- list()
for(i in 1:length(dep.fit.files)){
    load(dep.fit.files[i],e1<-new.env())
    e1$file <- dep.fit.files[i]
    results.list[[i]] <- e1
}

save(results.list,file="data/dep.fit.results.RData")
