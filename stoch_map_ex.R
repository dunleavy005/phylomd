rm(list = ls(all = TRUE))
setwd("C:/Users/Amrit/Documents/UW Research/Stochastic Mapping")
library(ape)
library(manipulate)
library(Matrix)
library(Rcpp)
library(rphast)
source("code/stoch_mapping_functions.R")
sourceCpp('C:/Users/Amrit/Downloads/expmat.cpp')
set.seed(0)


moments.ctmcjumps.ex = function(t, rate.mat, label.mat) {
  rate.mat.L=rate.mat*label.mat
  outp_list = vector("list", length=length(t))
  for (i in 1:length(t)) {
    tmp = matrix(0, nrow=12, ncol=12)
    tmp[1:3,1:3] = rate.mat
    tmp[4:6,4:6] = rate.mat
    tmp[7:9,7:9] = rate.mat
    tmp[10:12,10:12] = rate.mat
    tmp[1:3,4:6] = rate.mat.L
    tmp[4:6,7:9] = rate.mat.L
    tmp[7:9,10:12] = rate.mat.L
    tmp_array = array(expm(tmp*t[i])[1:nrow(rate.mat),], dim=c(nrow(rate.mat), ncol(rate.mat), 4))
    outp_array = array(0, dim=c(nrow(rate.mat), ncol(rate.mat), 4))
    
    for (j in 1:dim(tmp_array)[3]) {
      tmp_array[,,j] = tmp_array[,,j] * factorial(j-1)
    }
    
    stirling_num = stirling_num_table(3);
    
    for (j in 1:dim(tmp_array)[3]) {
      for (k in 1:j) {
        outp_array[,,j] = outp_array[,,j] + stirling_num[j,k] * tmp_array[,,k]
      }
    }
    
    outp_list[[i]] = outp_array
  }
  return(outp_list)
}


###### initialize tree objects

N = 3
tree = rtree(N, rooted=T)
tree = reorder(tree, order = "pruningwise")
rate.mat = matrix(c(-0.8,0.3,0.5,0.2,-0.7,0.5,0.2,0.3,-0.5), nrow=3, ncol=3, byrow=T)
label.mat = matrix(1,nrow=3,ncol=3) - diag(1,3)
rate.mat.L = rate.mat * label.mat
root.dist = c(0.2,0.3,0.5)
edge.set = 1:nrow(tree$edge)
##edge.moments = ctmc_moments_wrapper(tree$edge.length, rate.mat, rate.mat.L, 3)
edge.moments = moments.ctmcjumps.ex(tree$edge.length, rate.mat, label.mat)
num.edges = nrow(tree$edge)
num.states = nrow(rate.mat)
num.term.nodes = length(tree$tip.label)
root.node = length(tree$tip.label) + 1
root_branches = which(tree$edge[,1] == root.node)
root1 = root_branches[1]
root2 = root_branches[2]

I = function(cond) if(cond) 1 else 0

###### E(\omega_1^2 \omega_2)

edge.set1 = edge.set
edge.set2 = edge.set

lnames = c("{}", "{1}", "{2}", "{11}", "{12}", "{1,1}", "{1,2}", "{1,1,2}", "{{11},{2}}", "{{12},{1}}", "{112}")
V = replicate(length(lnames), matrix(0, nrow=num.states, ncol=num.edges), simplify=F)
names(V) = lnames

tip.01 = cbind(diag(1,3), diag(1,3), diag(1,3), diag(1,3), diag(1,3), diag(1,3))[,1:N]

for (i in 1:nrow(tree$edge)) {
  
  child_node = tree$edge[i,2]
  
  if (!any(tree$edge[,1] == child_node)) {
    ## terminal branch - initialize
    V[["{}"]][, i] = edge.moments[[i]][,,1] %*% tip.01[,child_node]
    V[["{1}"]][,i] = edge.moments[[i]][,,2] %*% tip.01[,child_node] * I(i %in% edge.set1)
    V[["{2}"]][,i] = edge.moments[[i]][,,2] %*% tip.01[,child_node] * I(i %in% edge.set2)
    V[["{11}"]][,i] = edge.moments[[i]][,,3] %*% tip.01[,child_node] * I(i %in% edge.set1)
    V[["{12}"]][,i] = edge.moments[[i]][,,3] %*% tip.01[,child_node] * I((i %in% edge.set1) & (i %in% edge.set2))
    V[["{112}"]][,i] = edge.moments[[i]][,,4] %*% tip.01[,child_node] * I((i %in% edge.set1) & (i %in% edge.set2))
  } else {
    ## internal branch - recurse
    child_branches = which(tree$edge[,1] == child_node)
    ch1_edge = child_branches[1]
    ch2_edge = child_branches[2]
    
    #### "{}", "{1}", "{2}", "{11}", "{12}", "{1,1}", "{1,2}", "{1,1,2}", "{{11},{2}}", "{{12},{1}}", "{112}"
    
    V[["{}"]][,i] = edge.moments[[i]][,,1] %*% (V[["{}"]][,ch1_edge] * V[["{}"]][,ch2_edge])
    
    V[["{1}"]][,i] = edge.moments[[i]][,,2] %*% (V[["{}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) * I(i %in% edge.set1) +
      edge.moments[[i]][,,1] %*% ((V[["{1}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]))
    
    V[["{2}"]][,i] = edge.moments[[i]][,,2] %*% (V[["{}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) * I(i %in% edge.set2) +
      edge.moments[[i]][,,1] %*% ((V[["{2}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]))
    
    V[["{11}"]][,i] = edge.moments[[i]][,,3] %*% (V[["{}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) * I(i %in% edge.set1) +
      edge.moments[[i]][,,1] %*% ((V[["{11}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{11}"]][,ch2_edge]))
    
    V[["{12}"]][,i] = edge.moments[[i]][,,3] %*% (V[["{}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) * I((i %in% edge.set1) & (i %in% edge.set2)) +
      edge.moments[[i]][,,1] %*% ((V[["{12}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{12}"]][,ch2_edge]))
    
    #bad
    V[["{112}"]][,i] = edge.moments[[i]][,,4] %*% (V[["{}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) * I((i %in% edge.set1) & (i %in% edge.set2)) +
      edge.moments[[i]][,,1] %*% ((V[["{112}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{112}"]][,ch2_edge]))
    
    V[["{1,1}"]][,i] = 2 * edge.moments[[i]][,,2] %*% ( (V[["{1}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) ) * I(i %in% edge.set1) +
      edge.moments[[i]][,,1] %*% ( (V[["{1,1}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1,1}"]][,ch2_edge]) +
                                     (2 * V[["{1}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) )
    
    V[["{1,2}"]][,i] = edge.moments[[i]][,,2] %*% ( (V[["{2}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]) ) * I(i %in% edge.set1) +
      edge.moments[[i]][,,2] %*% ( (V[["{1}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) ) * I(i %in% edge.set2) +
      edge.moments[[i]][,,1] %*% ( (V[["{1,2}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1,2}"]][,ch2_edge]) +
                                     (V[["{1}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]) + (V[["{2}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) )
    
    #bad
    V[["{{11},{2}}"]][,i] = edge.moments[[i]][,,3] %*% ( (V[["{2}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]) ) * I(i %in% edge.set1) +
      edge.moments[[i]][,,2] %*% ( (V[["{11}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{11}"]][,ch2_edge]) ) * I(i %in% edge.set2) +
      edge.moments[[i]][,,1] %*% ( (V[["{{11},{2}}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{{11},{2}}"]][,ch2_edge]) +
                                     (V[["{11}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]) + (V[["{2}"]][,ch1_edge] * V[["{11}"]][,ch2_edge]) )
    
    #bad
    V[["{{12},{1}}"]][,i] = edge.moments[[i]][,,3] %*% ( (V[["{1}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) ) * I((i %in% edge.set1) & (i %in% edge.set2)) +
      edge.moments[[i]][,,2] %*% ( (V[["{12}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{12}"]][,ch2_edge]) ) * I(i %in% edge.set1) +
      edge.moments[[i]][,,1] %*% ( (V[["{{12},{1}}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{{12},{1}}"]][,ch2_edge]) +
                                     2 * ((V[["{12}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) + (V[["{1}"]][,ch1_edge] * V[["{12}"]][,ch2_edge])) )
    
    #bad
    V[["{1,1,2}"]][,i] = 2 * edge.moments[[i]][,,2] %*% ( (V[["{1,2}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1,2}"]][,ch2_edge]) +
                                                            (V[["{1}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]) + (V[["{2}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) ) * I(i %in% edge.set1) +
      edge.moments[[i]][,,2] %*% ( (V[["{1,1}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1,1}"]][,ch2_edge]) +
                                     2 * V[["{1}"]][,ch1_edge] * V[["{1}"]][,ch2_edge] ) * I(i %in% edge.set2) +
      edge.moments[[i]][,,1] %*% ( (V[["{1,1,2}"]][,ch1_edge] * V[["{}"]][,ch2_edge]) + (V[["{}"]][,ch1_edge] * V[["{1,1,2}"]][,ch2_edge]) +
                                     (V[["{1,1}"]][,ch1_edge] * V[["{2}"]][,ch2_edge]) + (V[["{2}"]][,ch1_edge] * V[["{1,1}"]][,ch2_edge]) +
                                     2 * ((V[["{1,2}"]][,ch1_edge] * V[["{1}"]][,ch2_edge]) + (V[["{1}"]][,ch1_edge] * V[["{1,2}"]][,ch2_edge])) )
    
    
  }
  
  
  
}


exp112 = sum(root.dist %*% ( (V[["{112}"]][,root1] * V[["{}"]][,root2]) + (V[["{}"]][,root1] * V[["{112}"]][,root2]) +
                           (V[["{11}"]][,root1] * V[["{2}"]][,root2]) + (V[["{2}"]][,root1] * V[["{11}"]][,root2]) +    
                           (V[["{{11},{2}}"]][,root1] * V[["{}"]][,root2]) + (V[["{}"]][,root1] * V[["{{11},{2}}"]][,root2]) +  #### im here
                           (V[["{{12},{1}}"]][,root1] * V[["{}"]][,root2]) + (V[["{}"]][,root1] * V[["{{12},{1}}"]][,root2]) +
                           2 * ((V[["{12}"]][,root1] * V[["{1}"]][,root2]) + (V[["{1}"]][,root1] * V[["{12}"]][,root2])) +
                           (V[["{1,1}"]][,root1] * V[["{2}"]][,root2]) + (V[["{2}"]][,root1] * V[["{1,1}"]][,root2]) +
                           (V[["{1,1,2}"]][,root1] * V[["{}"]][,root2]) + (V[["{}"]][,root1] * V[["{1,1,2}"]][,root2]) +
                           2 * ((V[["{1,2}"]][,root1] * V[["{1}"]][,root2]) + (V[["{1}"]][,root1] * V[["{1,2}"]][,root2])) )) /
  sum(root.dist %*% (V[["{}"]][,root1] * V[["{}"]][,root2]))

exp112

x=phylojumps.sim(tree, rate.mat, root.dist, scale=F, states=0:2, seq.data=matrix(c(0,1,2)), N=100000)
mean(x^3)
