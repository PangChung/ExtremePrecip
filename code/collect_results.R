### collect results for bootstrap ### 
rm(list=ls())
library(stringr)
library(tidyr)
library(knitr)
library(kableExtra)
boot.files <- list.files(path="data/boot4/",pattern="fit_bootstrap_0_\\d+_\\d+.RData",full.names = TRUE)
numbers.idx = str_extract_all(boot.files,"\\d+")
region.idx <- sapply(numbers.idx,function(x){as.numeric(x[4])})
boot.idx <- sapply(numbers.idx,function(x){as.numeric(x[3])})

# for(i in 0:300){
#     for(j in 1:8){
#         file0 = paste0("data/boot4/fit_bootstrap_",i,"_",j,".RData")
#         file1 = paste0("data/boot4/fit_bootstrap_1_",i,"_",j,".RData")
#         if(file.exists(file0) & file.exists(file1)){
#             load(file0,e0<-new.env());load(file1,e1<-new.env())
#             result.list <- e0$result.list
#             result.list[[1]] <- e1$result.list[[1]]
#             print("done")
#             save(result.list,file=paste0("data/boot4/fit_bootstrap_0_",i,"_",j,".RData"))
#         }
#     }
#     print(i)
# }

results.boot.list <- list()
for(i in 1:length(boot.files)){
    load(boot.files[i],e<-new.env())
    results.boot.list[[i]] <- e$result.list
}

idx.valid = !unlist(lapply(results.boot.list,function(x){any(!unlist(lapply(x,function(y){length(y)==4})))}))

system(paste("rm",boot.files[!idx.valid]))
ist <- results.boot.list[idx.valid]
boot.files <- boot.files[idx.valid]
region.idx <- region.idx[idx.valid]
boot.idx <- boot.idx[idx.valid]
boot.collect <- function(data,idx,norm.idx,season.idx,region.idx,boot.idx){
    estimates.boot = matrix(unlist(lapply(data[region.idx==idx & boot.idx != 0],function(x){y=x[[norm.idx]][[season.idx]]$par;y[1] <- 2/(1+exp(-y[1]));y})),ncol=3,byrow=TRUE)
    est.true = results.boot.list[region.idx == idx  & boot.idx == 0][[1]][[norm.idx]][[season.idx]]$par; est.true[1] <- 2/(1+exp(-est.true[1]))
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

idx.grid = expand.grid(region=1:8,season=1:4,risk=1:2)
boot.result.df <- data.frame(
    region = idx.grid[,1],
    season = idx.grid[,2],
    risk = idx.grid[,3],
    shape = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[1]}),
    lambda0 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[2]}),
    lambda1 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[3]}),
    sd.shape = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$sd[1]}),
    sd.lambda0 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$sd[2]}),
    sd.lambda1 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$sd[3]}),
    low.shape = apply(idx.grid,1,function(x){quantile(boot.result.list[[x[3]]][[x[2]]][[x[1]]]$jack[,1],0.025)}),
    high.shape = apply(idx.grid,1,function(x){quantile(boot.result.list[[x[3]]][[x[2]]][[x[1]]]$jack[,1],0.975)}),
    low.lambda0 = apply(idx.grid,1,function(x){quantile(boot.result.list[[x[3]]][[x[2]]][[x[1]]]$jack[,2],0.025)}),
    high.lambda0 = apply(idx.grid,1,function(x){quantile(boot.result.list[[x[3]]][[x[2]]][[x[1]]]$jack[,2],0.975)}),
    low.lambda1 = apply(idx.grid,1,function(x){quantile(boot.result.list[[x[3]]][[x[2]]][[x[1]]]$jack[,3],0.025)}),
    high.lambda1 = apply(idx.grid,1,function(x){quantile(boot.result.list[[x[3]]][[x[2]]][[x[1]]]$jack[,3],0.975)})
)

boot.result.df[which(boot.result.df$low.lambda1 * boot.result.df$high.lambda1 > 0),1:3]
boot.result.df[which(abs(boot.result.df$lambda1) > boot.result.df$sd.lambda1*1.96),1:3]

save(boot.result.df,boot.result.list, results.boot.list, region.idx, boot.files, boot.idx,file="data/dep.fit.boot.results3.RData")

boot.result.df$shape.latex <- paste0("(",round(boot.result.df$low.shape,3),", ",round(boot.result.df$high.shape,3),")")

boot.result.df$lambda0.latex <- paste0("(",round(boot.result.df$low.lambda0,3),", ",round(boot.result.df$high.lambda0,3),")") 

boot.result.df$lambda1.latex <- paste0("(",round(boot.result.df$low.lambda1,3),", ",round(boot.result.df$high.lambda1,3),")") 

df.latex <- boot.result.df[,c(1:6,16,17,18)]
df.latex[,4:6] <- round(df.latex[,4:6],3)
df.latex$season = season[df.latex$season]
df.latex$region = factor(region.name[df.latex$region],levels=region.name)
df.latex <- rbind(df.latex,df.latex)
df.latex$type = rep(1:2,each=nrow(df.latex)/2)
df.latex$shape[df.latex$type==2] = df.latex$shape.latex[df.latex$type==2]
df.latex$lambda0[df.latex$type==2] = df.latex$lambda0.latex[df.latex$type==2]
df.latex$lambda1[df.latex$type==2] = df.latex$lambda1.latex[df.latex$type==2]
df.latex <- df.latex[,c(1:6,10)]

df.latex <- rbind(pivot_wider(df.latex[df.latex$type==1,],names_from=season,values_from=c(shape,lambda0,lambda1)),pivot_wider(df.latex[df.latex$type==2,],names_from=season,values_from=c(shape,lambda0,lambda1)))

df.latex = df.latex[order(df.latex$region,df.latex$risk),]
df.latex = df.latex[,c(1:3,order(rep(1:4,3)) + 3)]
# Use kable to create a basic table
table <- df.latex[,-c(2,3)] %>%
  kable(format = "latex", booktabs = TRUE)

# Use kableExtra to rotate the table
table <- table %>%
  kable_styling(latex_options = c("hold_position","scale_down")) %>%  landscape()
# Print the table
print(table)
