getwd()
# Save the data file to a location on your hard drive and specify the path here (Windows systems use forward slashes)
dir_path <-"C:/Users/ABHI/Desktop/sem1/social_media/adv"
#getwd()
setwd(dir_path)
# clear everything out of memory
rm(list=ls())  

# Load primary school data, contact data
infile_edges<-"edges.csv"
infile_nodes<-"nodes.csv"

## Load package
library(igraph)
edge_frame=read.csv(infile_edges, header = TRUE, sep = ",")
node_frame=read.csv(infile_nodes, header = TRUE, sep = ",")


## Create bins required for Part 2
# points
node_frame$ln_points_bins <- cut(node_frame$ln_points,5)
# fdi
node_frame$ln_fdi_bins <- cut(node_frame$ln_fdi,5)
# fdi_pergdp
node_frame$fdi_pergdp_bins <- cut(node_frame$fdi_pergdp,5)
# ICT_goods
node_frame$ICT_goods_bins <- cut(node_frame$ICT_goods_import_percent,5)
# internet
node_frame$internet_bins <- cut(node_frame$internet_users_percent,5)
# immigration
node_frame$immigration_bins <- cut(node_frame$immigration_pct,5)




# Create the graphs
SAP_graph_directed <- graph.data.frame(edge_frame, directed = TRUE, vertices= node_frame)
E(SAP_graph_directed)$weight <-1
SAP_graph<-simplify(SAP_graph_directed, edge.attr.comb="sum")



# Decomposing the graph to get main components - undirected
set.seed(1000)
temp <- graph.data.frame(edge_frame,  directed = FALSE, vertices= node_frame)
decompose_SAP <- decompose(temp, mode = c("weak", "strong"), max.comps = NA, min.vertices = 4)


# Recreate graph by combining the above main components
combined <- decompose_SAP[[1]]
SAP_main <- for(i in 2:length(decompose_SAP))
{
  temp <- decompose_SAP[[i]]
  combined <- union(combined,temp)
  print(i)
}

# Louvain algorithm
name <- "Louvain algorithm"
set.seed(1000)
SAP_comm_louvain <- cluster_louvain(combined, weights = E(combined)$weight)
c.m <- membership(SAP_comm_louvain)
nclust <- max(c.m)
# zzz <- table(c.m, emp.country, useNA = c("no"))

# Fast Greedy algorithm
name <- "Fast Greedy algorithm"
set.seed(1000)
SAP_comm_fast <- fastgreedy.community(combined, weights=E(combined)$weight)
c.m <- membership(SAP_comm_fast)
nclust <- max(c.m)
# zzz <- table(c.m, emp.country, useNA = c("no"))


# Walktrap algorithm
name <- "Walktrap algorithm"
set.seed(1000)
SAP_comm_fast <- walktrap.community(combined, weights=E(combined)$weight)
c.m <- membership(SAP_comm_fast)
nclust <- max(c.m)
# zzz <- table(c.m, emp.country, useNA = c("no"))


# Spinglass algorithm
name <- "Spinglass algorithm"
set.seed(1000)
SAP_comm_fast <- spinglass.community(combined, weights=E(combined)$weight)
c.m <- membership(SAP_comm_fast)
nclust <- max(c.m)
# zzz <- table(c.m, emp.country, useNA = c("no"))

# Label Propogation algorithm
name <- "Label Propogation algorithm"
set.seed(1000)
SAP_comm_fast <- label.propagation.community(combined, weights=E(combined)$weight)
c.m <- membership(SAP_comm_fast)
nclust <- max(c.m)


## Create a function for extracting the giant component
giant.component <- function(x) 
{
  cl <- clusters(x)
  # induced.subgraph(x, which(cl$membership == which.max(cl$csize)-1)-1)
  induced.subgraph(x, which(cl$membership == which.max(cl$csize)))
}

# Getting the giant component
SAP_giant_dir <- giant.component(SAP_graph)
is.simple(SAP_giant_dir)
is.directed(SAP_giant_dir)
is.weighted(SAP_giant_dir)


## After running multiple algorithms, we find that choosing the Giant Component that was created by 
## using an undirected, non-simplified input graph with edge weights = 1 
## and the Louvian algorithm for community detection
## gives the lowest number of communities (91) --> If you use a simplified graph, you get 91 communities
## The same has been created below

# Louvain algorithm
name <- "Louvain algorithm"
set.seed(1000)
# E(SAP_giant_dir)$weight <- 1
louvain_giant <- cluster_louvain(as.undirected(SAP_giant_dir), weights = E(as.undirected(SAP_giant_dir))$weight)
# louvain_giant <- cluster_louvain(SAP_giant_dir, weights = E(SAP_giant_dir)$weight)
c.m <- membership(louvain_giant)
nclust <- max(c.m)      # This gives no. of communities(nclust) = 91


## Plotting the giant component and the communities
plot(louvain_giant, SAP_giant_dir, 
     vertex.label= NA, 
     vertex.size=2, edge.arrow.size=.2,
     main="Louvain algorithm",
     layout= layout.fruchterman.reingold, 
     colors=c.m)


################################################################################################################
### Find the network metrics to decribe the network structure


## ASSORTATIVITY - DIRECTED GRAPH -------> Part1 (For giant component)

# Degree assortativity
is.directed(SAP_giant_dir)
giant_degr<-degree(SAP_giant_dir)
assortativity(SAP_giant_dir, types1=giant_degr, directed=TRUE)
# Conclusion
# -0.01697888 - negative assortativity shows presense of a hub and spoke structure

# Assortativity for ln_points
emp.ln_points <- get.vertex.attribute(SAP_giant_dir, "ln_points")
emp.ln_points[is.na(emp.ln_points)] <- 0
assortativity(SAP_giant_dir, types1=emp.ln_points, directed=TRUE)
# Conclusion
# [1] -0.03240836 - shows the 2 immediate points which are connected have
# opposite level of intelligence coefficient- negative homophily
modularity(louvain_giant)
# 0.7945445 // This number is changed to 0.7858535
class(louvain_giant)
vcount(SAP_giant_dir)
comm <- communities(louvain_giant)


# Extarcing sub graph based on communities from giant graph - individual assortativity
# Degree assortativity of individual communities
individual.assortativity<-NULL
for (i in seq(1,91,1)) 
{
  gr <- induced.subgraph(SAP_giant_dir, louvain_giant$membership == i)
  dg <- degree(gr)
  community <- i
  as <- assortativity(gr,types1 = dg)
  a <- data.frame(community, as)
  individual.assortativity <- rbind(individual.assortativity,a)
}

# ln_points assortativity of individual communities
individual.assortativity_ln<-NULL
for (i in seq(1,91,1)) {
  gr <- induced.subgraph(SAP_giant_dir, louvain_giant$membership == i)
  emp.ln_points <- get.vertex.attribute(gr, "ln_points")
  emp.ln_points[is.na(emp.ln_points)] <- 0
  bb <- assortativity(gr, types1 = emp.ln_points, directed=TRUE)
  community <- i
  ln <- data.frame(community,bb )
  individual.assortativity_ln<-rbind(individual.assortativity_ln,ln)
}


## Other features 

# In the giant network
# OVerall (global) clustering = 0.005204738
overall_clustering <- transitivity(SAP_giant_dir)
# Individual clustering
clustering_i <- transitivity(SAP_giant_dir, "local")
# Average clustering = 0.0352601
avg_clustering <-mean(clustering_i, na.rm = TRUE)


# Clustering in sub communities

sub_comm<-NULL
for (i in seq(1,91,1))
{
  
  gr<-induced.subgraph(SAP_giant_dir,louvain_giant$membership == i)
  
  #individual vs overall clustering in subgraphs
  overall_clustering <- transitivity(gr)
  
  # Individual clustering
  clustering_i <- transitivity(gr, "local")
  # Average clustering
  avg_clustering <-mean(clustering_i, na.rm = TRUE)
  cc<-data.frame(i, overall_clustering, avg_clustering)
  sub_comm<-rbind(sub_comm,cc)
}


## ASSORTATIVITY - DIRECTED GRAPH -------> Part2 (For the entire graph)

# Degree assortativity
is.directed(SAP_graph)
giant_degr <- degree(SAP_graph)
assortativity(SAP_graph, types1=giant_degr, directed=TRUE)
# Conclusion
#-0.01096796 - positive assortativity shows presense of a hub and spoke structure

# Assortativity for ln_points
emp.ln_points <- get.vertex.attribute(SAP_graph, "ln_points")
emp.ln_points[is.na(emp.ln_points)] <- 0
assortativity(SAP_graph, types1=emp.ln_points, directed=TRUE)
# Conclusion
# [1] 0.02234436 - shows the 2 immediate points which are connected have
# Similar level of intelligence coefficient- positive homophily


# Individual vs overall clustering in the entire graph
# Overall clustering = 0.005203974
overall_clustering_g <- transitivity(SAP_graph)         
# Individual clustering : 
clustering_i_g <- transitivity(SAP_graph, "local")
# Average clustering = 0.02770246
avg_clustering_g <- mean(clustering_i_g, na.rm = TRUE)


# Avg. path length and diameter
average.path.length(SAP_giant_dir, directed=TRUE)
#7.855268
diameter(SAP_giant_dir)
#29 ---> This is coming out to be 33
# Summarize the graph structure
summary(SAP_giant_dir)


## Use the inverse of log weight for some of the network measure calculations
inv_weight <- 1/log(E(SAP_giant_dir)$weight  + 1)
num_weight <- E(SAP_giant_dir)$weight 

# Embeddedness/ inverse of structural hole access (see Burt 2004)
constraints_SAP <- round(constraint(SAP_giant_dir, nodes=V(SAP_giant_dir)), digits=4)
# Degree centrality
degree_sap <- degree(SAP_giant_dir)
# Node betweenness
betweens_SAP <- round(betweenness(SAP_giant_dir, v=V(SAP_giant_dir), directed = TRUE, nobigint =TRUE, normalized = FALSE))
# Edge betwenness
edgebetweens_SAP<-edge.betweenness(SAP_giant_dir, e=E(SAP_giant_dir), directed = TRUE)
# Local clustering coefficients
clustering_SAP <- transitivity(SAP_giant_dir, type="local", vids=V(SAP_giant_dir)) 


### Plots for inference

# PLOTTING CHARTS FOR Q1
plot1<-hist(giant_degr,col="blue",
            xlab="Degree", ylab="Frequency",
            main="Degree distribution: SAP Network- Giant Component")

# Plots 1 and 2: Can run them together
par(mfrow=c(1, 2))
edge_frame<-data.frame(edgebetweens_SAP, num_weight, inv_weight)
a_edge<-aggregate(edgebetweens_SAP ~ inv_weight, data=edge_frame, mean)
plot(a_edge, col="blue", log="xy", xlab="Weight of edge", ylab="Average Betweenness of edges")
node_frame<-data.frame(betweens_SAP, constraints_SAP, clustering_SAP, degree_sap)
a_node<-aggregate(betweens_SAP ~ clustering_SAP, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Clustering", ylab="Average Betweenness of nodes")


# Plot set 2: Four plots 
par(mfrow=c(2, 2))
a_node<-aggregate(betweens_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Betweenness")
a_edge<-aggregate(edgebetweens_SAP ~ num_weight, data=edge_frame, mean)
plot(a_edge, col="blue", log="xy", xlab="Weight of edge", ylab="Average Betweenness of edges")
a_node<-aggregate(clustering_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Clustering")
a_node<-aggregate(constraints_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Constraint (Embeddedness)")

# Log-log degree distributino
par(mfrow=c(1, 2))
d.net <-degree(SAP_giant_dir)
dd.net <- degree.distribution(SAP_giant_dir)
d <- 1:max(d.net)-1
ind <- (dd.net != 0)
plot(d[ind], dd.net[ind], log="xy", col="blue",
     xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
     main="Log-Log Degree Distribution")

# CHUNK 8# Average neighbor degree versus vertex degree
a.nn.deg <- graph.knn(SAP_giant_dir,V(SAP_giant_dir))$knn
plot(d.net, a.nn.deg, log="xy", 
     col="goldenrod", xlab=c("Log Vertex Degree"),
     ylab=c("Log Average Neighbor Degree"))



#---------------------------------------------------------------------------------------------------------------
# PART B
#---------------------------------------------------------------------------------------------------------------

# Plotting the distribution of nodes across all the countries
c.m <- membership(louvain_giant)
nclust <- max(c.m)      # This gives no. of communities(nclust) = 91
emp.country <- get.vertex.attribute(SAP_giant_dir, "country")


# Checking for country 
t <- table(c.m, emp.country, useNA = c("no"))   # sum(!is.na(t))
data_country <- as.data.frame(t)    # country <- as.data.frame.matrix(t)
colnames(data_country)[3] <- "Num_nodes"
data_country$emp.country <- as.character(data_country$emp.country)
data_country$emp.country[data_country$emp.country==''] <- "Unknown"
temp <- aggregate(Num_nodes ~ emp.country, data_country, sum)
colnames(temp)[2] <- "Num_nodes_country"
data_country <- merge(x = data_country, y = temp, by.x = "emp.country", by.y = "emp.country", all.x = TRUE)
data_country <- data_country[order(-data_country$Num_nodes_country, -data_country$Num_nodes),]

par(mfrow=c(1, 1))
barplot(height=data_country$Num_nodes,data_country$Num_nodes_country,names.arg=data_country$emp.country)


# Set n/w attributes for running ANOVA tests
sapnew.country <- get.vertex.attribute(SAP_giant_dir, "country")
sapnew.points_bins <- get.vertex.attribute(SAP_giant_dir, "ln_points_bins")
sapnew.fdi_bins <- get.vertex.attribute(SAP_giant_dir, "ln_fdi_bins")
sapnew.pergdp_bins <- get.vertex.attribute(SAP_giant_dir, "fdi_pergdp_bins")
sapnew.goods_bins <- get.vertex.attribute(SAP_giant_dir, "ICT_goods_bins")
sapnew.internet_bins <- get.vertex.attribute(SAP_giant_dir, "internet_bins")
sapnew.immigration_bins <- get.vertex.attribute(SAP_giant_dir, "immigration_bins")
sapnew.lat_bins <- get.vertex.attribute(SAP_giant_dir, "lat_bins")
sapnew.lng_bins <- get.vertex.attribute(SAP_giant_dir, "lng_bins")


#checking for country 
a1<-table(c.m, sapnew.country, useNA = c("no"))
c1<-as.data.frame(a1)
anova1<-aov(Freq~sapnew.country, data=c1)
summary(anova1)

#anova for points
a2<-table(c.m, sapnew.points_bins, useNA = c("no"))
c2<-as.data.frame(a2)
anova2<-aov(Freq~sapnew.points_bins, data=c2)
summary(anova2)

#anova for fdi
a3<-table(c.m, sapnew.fdi_bins, useNA = c("no"))
c3<-as.data.frame(a3)
anova3<-aov(Freq~sapnew.fdi_bins, data=c3)
summary(anova3)

#anova for fdi_pergdp
a4<-table(c.m, sapnew.pergdp_bins, useNA = c("no"))
c4<-as.data.frame(a4)
anova4<-aov(Freq~sapnew.pergdp_bins, data=c4)
summary(anova4)

#anova for ICT_goods_bins
a5<-table(c.m, sapnew.goods_bins, useNA = c("no"))
c5<-as.data.frame(a5)
anova5<-aov(Freq~sapnew.goods_bins, data=c5)
summary(anova5)

#anova for internet_bins
a6<-table(c.m, sapnew.internet_bins, useNA = c("no"))
c6<-as.data.frame(a6)
anova6<-aov(Freq~sapnew.internet_bins, data=c6)
summary(anova6)

#anova for immigration_bins
a7<-table(c.m, sapnew.immigration_bins, useNA = c("no"))
c7<-as.data.frame(a7)
anova7<-aov(Freq~sapnew.immigration_bins, data=c7)
summary(anova7)

# Binning for latitudes and longitues does not make sense
# #anova for lat_bins
# a8<-table(c.m, sapnew.lat_bins, useNA = c("no"))
# c8<-as.data.frame(a8)
# anova8<-aov(Freq~sapnew.lat_bins, data=c8)
# summary(anova8)
# 
# #anova for lng_bins
# a9<-table(c.m, sapnew.lng_bins, useNA = c("no"))
# c9<-as.data.frame(a9)
# anova9<-aov(Freq~sapnew.lng_bins, data=c9)
# summary(anova9)


## Conclusion
# All the above node attributes have a SIGNIFICANT(based on p-value) effect 
# on the distribution  of nodes among the communities


#---------------------------------------------------------------------------------------------------------------
# PART C
#---------------------------------------------------------------------------------------------------------------

# Generate attributes required for finding influencers on ln_points
c.m <- membership(louvain_giant)
nclust <- max(c.m)      # This gives no. of communities(nclust) = 78
emp.country <- get.vertex.attribute(SAP_giant_dir, "country")
emp.ln_points <- get.vertex.attribute(SAP_giant_dir, "ln_points")
emp.ln_fdi <- get.vertex.attribute(SAP_giant_dir, "ln_fdi")
emp.fdi_pergdp <- get.vertex.attribute(SAP_giant_dir, "fdi_pergdp")
emp.ICT_goods_import_percent <- get.vertex.attribute(SAP_giant_dir, "ICT_goods_import_percent")
emp.internet_users_percent <- get.vertex.attribute(SAP_giant_dir, "internet_users_percent")
emp.immigration_pct <- get.vertex.attribute(SAP_giant_dir, "immigration_pct")
emp.lat <- get.vertex.attribute(SAP_giant_dir, "lat")
emp.lng <- get.vertex.attribute(SAP_giant_dir, "lng")


degree_sap_dir <- degree(SAP_giant_dir)   # Degree centrality
betweens_SAP_dir <- betweenness(SAP_giant_dir, v=V(SAP_giant_dir), directed = TRUE, 
                                nobigint =TRUE, normalized = FALSE)       # Node betweenness
edgebetweens_SAP_dir<-edge.betweenness(SAP_giant_dir, e=E(SAP_giant_dir), directed = TRUE)    # Edge betwenness
individual_clustering_SAP_dir <- transitivity(SAP_giant_dir, type="local", vids=V(SAP_giant_dir))   # Local clustering coefficients
SAP_giant_dir$hub_scr<-hub.score(SAP_giant_dir)$vector      # Hub score
SAP_giant_dir$auth_scr<-authority.score(SAP_giant_dir)$vector     # Authority score


# Setting the node attributes
SAP_giant_dir <- set.vertex.attribute(SAP_giant_dir, "community", value=c.m)
SAP_giant_dir <- set.vertex.attribute(SAP_giant_dir, "degree_sap_dir", value=degree_sap_dir)
SAP_giant_dir <- set.vertex.attribute(SAP_giant_dir, "betweens_SAP_dir", value=betweens_SAP_dir)
SAP_giant_dir <- set.vertex.attribute(SAP_giant_dir, "individual_clustering_SAP_dir", value=individual_clustering_SAP_dir)
SAP_giant_dir <- set.vertex.attribute(SAP_giant_dir, "hub_scr", value=hub.score(SAP_giant_dir)$vector)
SAP_giant_dir <- set.vertex.attribute(SAP_giant_dir, "auth_scr", value=authority.score(SAP_giant_dir)$vector)

# Setting the edge attributes
SAP_giant_dir <- set.edge.attribute(SAP_giant_dir, "edgebetweens_SAP_dir", value=edgebetweens_SAP_dir)


## In order to do the regression, we need to extraxt the data from the graph and convert it into a data frame
SAP_giant_dir.edgeframe <- NULL
SAP_giant_dir.edgeframe <- as.data.frame(get.edgelist(SAP_giant_dir))
SAP_giant_dir.nodeframe <- NULL
SAP_giant_dir.nodeframe <- as.data.frame(lapply(list.vertex.attributes(SAP_giant_dir),function(x) get.vertex.attribute(SAP_giant_dir,x)))
# Rename columns to appropriate mnemonic names
colnames(SAP_giant_dir.nodeframe) <- list.vertex.attributes(SAP_giant_dir)



## Aggregate attributes to a community level to create the dataset for regression
SAP1 <- aggregate(cbind(ln_points, ln_fdi, fdi_pergdp, ICT_goods_import_percent, 
                        internet_users_percent, immigration_pct) 
                  ~ community, data = SAP_giant_dir.nodeframe, sum)


# Merge SAP1 with the community attributes

comm_attr <- NULL
for (i in seq(1,91,1))
{
  
  gr<-induced.subgraph(SAP_giant_dir,louvain_giant$membership == i)
  
  # # Overall clustering
  # overall_clustering_comm <- transitivity(gr)
  # 
  # # Average individual clustering
  # clustering_i_comm <- transitivity(gr, "local")
  # avg_clustering_comm <- mean(clustering_i, na.rm = TRUE)
  
  # Degree centrality
  degree_sap_dir_comm <- degree(gr)
  avg_degree <- mean(degree_sap_dir_comm)
  
  # Node betweenness
  betweens_SAP_dir_comm <- betweenness(gr, v=V(gr), directed = TRUE, 
                                       nobigint =TRUE, normalized = FALSE)
  avg_betweens <- mean(betweens_SAP_dir_comm)
  
  # Hub score
  gr$hub_scr <- hub.score(gr)$vector      
  avg_hub_scr <- mean(gr$hub_scr)
  
  # Authority score
  gr$auth_scr <- authority.score(gr)$vector
  avg_auth_scr <- mean(gr$auth_scr)
  
  cc <- data.frame(i, 
                   # overall_clustering, avg_clustering, 
                   avg_degree, avg_betweens, avg_hub_scr, avg_auth_scr)
  
  comm_attr<-rbind(comm_attr,cc)
  
}
colnames(comm_attr)[1] <- "community"

# Creating the analytical dataset for regression
ads <- cbind(SAP1, sub_comm[-1], comm_attr[-1])


## Perform the regression

# Check for linearity in parameters
pairs(ads[-1])

fit <- lm(ln_points ~ .-community, data=ads)
summary(fit)
fit
## Conclusion
# The results of the regression show that there is a negative relation betwwen the degree and betweenness
# centralities and the knowledge contained within the the network. This does not make sense. 
# Such a result might be because the effects of two or more of the independent variables are 
# confounded with each other. This can be checked by running the regression using just the network parameters
# This has been done below

fit <- lm(ln_points ~ overall_clustering + avg_clustering + avg_degree + 
            avg_betweens + avg_hub_scr + avg_auth_scr, data=ads)
summary(fit)
fit
# As hypothesized we now see a positive relation between the knowledge contained within a community and the  
# degreee and betweenness centralities within the network

