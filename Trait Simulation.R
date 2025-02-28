library(phytools)
library(ggplot2)
library(tidyr)
library(reshape2)

##Number of Cells
N_c<-1000
##Number of traits
N_t<-100

##Simulate a pure-birth phenotype tree
tree<-pbtree(n=1000,b=1,d=0)


trait_stimulation<-function(tree,N_c,N_t,sig2,fig,title){
  y<- matrix(0, nrow = N_c, ncol = N_t)
  ##Simulate the traits under BM
  for(i in 1:N_t){
    y[,i]<-fastBM(tree,mu=0, sig2=sig2[i],bounds=c(-1,1),nsim=1)
  }
  # order the tips
  tip_labels <- as.numeric(gsub("t", "", tree$tip.label))
  y_new<-y[order(tip_labels),]
  gene_expression_matrix<-t(y_new)
  rownames(gene_expression_matrix) <- paste("trait", 1:N_t, sep="_") 
  colnames(gene_expression_matrix) <- paste("Cell", 1:N_c, sep="_")  
  gene_expression_long <- melt(gene_expression_matrix)
  colnames(gene_expression_long) <- c("Trait", "Cell", "Expression") 
  
  # Plot the heatmap
  p<-ggplot(gene_expression_long, aes(x = Cell, y = Trait, fill = Expression)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +  
    theme_minimal() + 
    scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-1, 1),  
    oob = scales::squish 
  ) +
    labs(title = title, x = "Cells", y = "Traits") +theme(
      title =  element_text(family = "Calibri", size = 18,face = 'bold'),
    axis.text.x = element_blank(),
    axis.title.x = element_text(family = "Calibri", size = 18,face = 'bold'), 
    axis.title.y = element_text(family = "Calibri", size = 18,face = 'bold'),
    legend.title   = element_text(family = "Calibri", size = 18,face = 'bold'),
    legend.text =  element_text(family = "Calibri", size = 16),
    plot.title = element_text(hjust = 0.5)) 
  ggsave(fig,plot=p,width = 10,height=8)
  return(y_new)
}

##Simulate the trait evolution when the mutation rate is fixed across all the traits
set.seed(123)
y1<-trait_stimulation(tree,N_c,N_t=100,rep(0.001,100),'trait_value_heatmap_001.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " =0.001"))))
y2<-trait_stimulation(tree,N_c,N_t=100,rep(0.01,100),'trait_value_heatmap_01.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " =0.01"))))
y3<-trait_stimulation(tree,N_c,N_t=100,rep(0.05,100),'trait_value_heatmap_05.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " =0.05"))))
y4<-trait_stimulation(tree,N_c,N_t=100,rep(0.1,100),'trait_value_heatmap_10.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " =0.1"))))
####Simulate the trait evolution when the mutation rate is varying and generated from a given normal distribution
sig1<-abs(rnorm(N_t,mean=0.001,sd=0.1))
sig2<-abs(rnorm(N_t,mean=0.01,sd=0.1))
sig3<-abs(rnorm(N_t,mean=0.05,sd=0.1))
sig4<-abs(rnorm(N_t,mean=0.1,sd=0.1))
y5<-trait_stimulation(tree,N_c,N_t,sig1,'trait_value_vary_heatmap_001.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " ~ ", N(0.001, 0.01)))))
y6<-trait_stimulation(tree,N_c,N_t,sig2,'trait_value_vary_heatmap_01.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " ~ ", N(0.01, 0.01)))))
y7<-trait_stimulation(tree,N_c,N_t,sig3,'trait_value_vary_heatmap_05.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " ~ ", N(0.05, 0.01)))))
y8<-trait_stimulation(tree,N_c,N_t,sig4,'trait_value_vary_heatmap_10.jpg',expression(bold(paste("Trait Value Heatmap when ",sigma^2, " ~ ", N(0.1, 0.01)))))
