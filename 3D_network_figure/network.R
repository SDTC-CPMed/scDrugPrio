library(shiny)
library(shinyjs)
library(igraph)
library(plotly)




#Read the data
network = read.table('input_for_network_updated.txt', sep = "\t", header = TRUE, stringsAsFactors = FALSE)
network = network[order(network$Molecule_A),]
rownames(network) = 1:length(rownames(network))

#Delete those connections where Molecule B is not available in "molecule_A" columns --> no information about what
#kind of protein it is
network = network[-703063,]
network = network[-703095,]
network = network[-703127,]
network = network[-703303,]
network = network[-703341,]
network = network[-703370,]
network = network[-703376,]
network = network[-703530,]
network = network[-703616,]
network = network[-703807,]
network = network[-703814,]
network = network[-704002,]
rownames(network) = 1:length(rownames(network))

#Scale the size
network$Drug_ranking = network$Drug_ranking*1.7

#names - all nodes
names = unique(network$Molecule_A)

#set the colors
col = rep(0, length(names))
for (i in 1:length(names)){
  idx = which(network$Molecule_A == names[i])[1]
  col[i] = network$Color_indication[idx]
}

#Include only interaction where one of the nodes is DEG
DEGs = which(col > 0)
ind_DEGs = rep(0, length(network$Molecule_A))
j = 1
for(i in 1:length(network$Molecule_A)){
  if(network$Molecule_A[i] %in% names[DEGs] && network$Molecule_B[i] %in% names[DEGs] && network$Molecule_A[i] != network$Molecule_B[i]){
    ind_DEGs[j] = i
    j = j+1
  }
}
ind_DEGs = ind_DEGs[1:(j-1)]
network_DEGs = network[ind_DEGs,]
rownames(network_DEGs) = 1:length(rownames(network_DEGs))

#names_DEGs - nodes in network of only DEGs
names_DEGs = unique(network_DEGs$Molecule_A)


#Get all the variables that are needed 
rank = network_DEGs$Drug_ranking
rank[is.na(rank)] = 0
Score = rep(0, length(names_DEGs))
col = rep(0, length(names_DEGs))
node_size = rep(0, length(names_DEGs))
hover_name = rep(0, length(names_DEGs))
Mtype = rep(0, length(names_DEGs))
for (i in 1:length(names_DEGs)){
  idx = which(network_DEGs$Molecule_A == names_DEGs[i])[1]

  Score[i] = rank[idx]
  col[i] = network_DEGs$Color_indication[idx]
  if (network_DEGs$Molecule_type[idx] == 'drug'){
    node_size[i] = 3.0
    hover_name[i] = network_DEGs$Molecule_A[idx]
  }
  if (network_DEGs$Molecule_type[idx] == 'protein'){
    node_size[i] = network_DEGs$Average_FC[idx]
    hover_name[i] = network_DEGs$Molecule_description[idx]
  }
  Mtype[i] = network_DEGs$Molecule_type[idx]
}





#Only keep one connection if there are two in both directions
ind_DEGs = rep(0, length(network_DEGs$Molecule_A))
j = 1
for(i in 1:length(network_DEGs$Molecule_A)){
  if(network_DEGs$Molecule_A[i] > network_DEGs$Molecule_B[i]){
    ind_DEGs[j] = i
    j = j+1
  }
}
ind_DEGs = ind_DEGs[1:j]
network_DEGs = network_DEGs[ind_DEGs,]
rownames(network_DEGs) = 1:length(rownames(network_DEGs))


# Set how the things should look like
col_DEGs = factor(col, labels = c("blue", "red", "yellow"))
score_DEGs  = Score
node_size_DEGs = (0.15+node_size)*2
hover_name_DEGs = hover_name


# Set sources and targets as needed for the network
Source = rep(0, length(network_DEGs$Molecule_A))
Target = rep(0, length(network_DEGs$Molecule_A))
for (i in 1:length(network_DEGs$Molecule_A)){
  Source[i] = which(network_DEGs$Molecule_A[i] == names_DEGs)
  Target[i] = which(network_DEGs$Molecule_B[i] == names_DEGs)
}


# Get the Targets - This will be needed for lighting up the neighboring nodes after the click
Targets = list()
for (i in 1:length(names_DEGs)){
  idx = which(network_DEGs$Molecule_A == names_DEGs[i])
  t = list()
  k = 1
  
  
  for(j in idx){
    t[k] = network_DEGs$Molecule_B[j]
    k = k+1
  }
  idx = which(network_DEGs$Molecule_B == names_DEGs[i])
  for (j in idx){
    t[k] = network_DEGs$Molecule_A[j]
    k = k+1
  }
  Targets <- append(Targets, list(t))
}


#Get the neighboors. This will be needed for lighting up the neighboring nodes after the click
Neighboors = list()
for(i in 1:length(names_DEGs)){
  N = list()
  for(j in 1:length(Targets[[i]])){
    N[j] = which(Targets[[i]][[j]] == names_DEGs)
  }
  N[j+1] = i
  Neighboors <- append(Neighboors, list(N))
}

#Get the connections. This will be needed for lighting up the neighboring nodes after the click
Connections = list()
for(i in 1:length(names_DEGs)){
  C1 = which(names_DEGs[i] == network_DEGs$Molecule_A)
  C2 = which(names_DEGs[i] == network_DEGs$Molecule_B)
  C = c(C1, C2)
  Connections = append(Connections, list(C))
}




# Get the edges as needed for the network
Edges = rep(0, 2*length(Source))
for (i in 1:length(Source)){
  Edges[2*i-1] = Source[i]
  Edges[2*i] = Target[i]
}
length(unique(Edges))

#Create the graph
G = graph(Edges, directed = FALSE)
set.seed(42)
layt=layout.drl(G, options = list(cooldown.attraction = 0.4, cooldown.temperature=1000, cooldown.damping.mult = 1, crunch.temperature = 15000))


#Get the coordinates
Xn = rep(0, length(col_DEGs))
Yn = rep(0, length(col_DEGs))
N = length(col_DEGs)

for(k in 1:N){
  Xn[k]=round(layt[k,1], 3)# x-coordinates of nodes
  Yn[k]=round(layt[k,2], 3)# y-coordinates
}

Zn = score_DEGs
Xe= rep(0, 3*length(Source))
Ye= rep(0, 3*length(Source))
Ze= rep(0, 3*length(Source))

library(rt3)

for (i in 1:length(Source)){
  Xe[3*i-2] = round(layt[Edges[2*i-1],1], 3)# x-coordinates of edge ends
  Xe[3*i-1] = round(layt[Edges[2*i],1], 3)# x-coordinates of edge ends
  Xe[3*i] = NONE# x-coordinates of edge ends
  Ye[3*i-2] = round(layt[Edges[2*i-1],2], 3)# x-coordinates of edge ends
  Ye[3*i-1] = round(layt[Edges[2*i],2], 3)# x-coordinates of edge ends
  Ye[3*i] = NONE# x-coordinates of edge ends
  
  Ze[3*i-2] = round(Zn[which(names_DEGs == network_DEGs$Molecule_A[i])[1]], 3)
  Ze[3*i-1] = round(Zn[which(names_DEGs == network_DEGs$Molecule_B[i])[1]], 3)
  Ze[3*i] = NONE
}

line_data = data.frame(cbind(Xe,Ye,Ze))

point_data = data.frame(cbind(Xn,Yn,Zn,col_DEGs))




#save everything that will be needed 
save(fig,Neighboors,point_data,col_DEGs,hover_name_DEGs,names_DEGs,node_size_DEGs,score_DEGs, Connections, line_data, cam,  file="shiny.RData")

