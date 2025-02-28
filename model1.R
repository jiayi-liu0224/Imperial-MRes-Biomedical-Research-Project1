### Stimulate the data first
library(rstan)
library(phytools)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

##read the tree and traits when sig2=0.001
N_c<-100
N_t<-100
z0<-0
tree<-readRDS("/rds/general/user/jl1924/home/tree.rds")
y<-readRDS("/rds/general/user/jl1924/home/trait1.rds")
C<-vcv.phylo(tree)

###run the stan model
data=list(N_t=N_t,N_c=N_c,z0=z0,C=C,Y=y)
stanmodel<-stan_model(file='/rds/general/user/jl1924/home/non_hier_rscript_fit/non_hier1.stan',model_name = 'non_hier_tree1')
stanmodel.fit<-sampling(stanmodel,data=data, chains = 4, cores = 8,
                         iter = 4000, warmup = 500)

saveRDS(y,file='/rds/general/user/jl1924/home/non_hier_rscript_fit/non_hier1_trait.rds')
saveRDS(stanmodel.fit,file = '/rds/general/user/jl1924/home/non_hier_rscript_fit/non_hier1_fit.rds')