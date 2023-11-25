marginal.fit.files <- list.files("data/",pattern = "marginal_fit_.*.RData",full.names = TRUE)
dep.fit.files <- list.files("data/",pattern="fit_pot_ST_*",full.names = TRUE)

for(i in 1:length(marginal.fit.files)){
    load(marginal.fit.files[i],e<-new.env())
    print(is.null(e$data.df.gpd))
}

time.vec <- c()
time.elapse <- c()
estimates <- c()

for(i in 1:length(dep.fit.files)){
    load(dep.fit.files[i],e1<-new.env())
    time.vec <- c(time.vec,e1$result$time[5])
    time.elapse <- c(time.elapse,e1$result$time[3])
    estimates <- rbind(estimates,e1$result$par)
}
