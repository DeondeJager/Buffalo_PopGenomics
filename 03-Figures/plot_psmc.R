# Function to plot PSMC using .psmc output
# Code from https://figshare.com/articles/Plot_PSMC_results/3996156/1
# Also available at: https://datadryad.org/stash/dataset/doi:10.5061/dryad.0618v
# which is from this paper: https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12606

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  # extract data for each iteration (30)
  
  START<-grep("^RD",X) # line numbers
  END<-grep("^//",X) # line numbers
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE) # \theta_0 \rho_0
  RS<-grep("^RS",X,value=TRUE) # k t_k \lambda_k \pi_k \sum_{l\not=k}A_{kl} A_{kk}
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s # scale 

  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  a$Generation<-as.numeric(2*N0*a[,3])
  a$Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne,mu,g)
  #plot(Ne~YearsAgo)
}
