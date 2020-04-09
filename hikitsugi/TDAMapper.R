generate.mapper <- function(X,slice,X.dist=dist(X),num_interval=10,num_bins=10,dont.plot=F){
  require(TDAmapper)
  require(igraph)
  require(TDA)
  if (missing(slice)) {
    m <- apply(X, 2, mean)
    slice <- apply(X,1,function(z)sum(mapply(function(x,y)abs(x^2-y^2), m,z))^(1/2))
    print("sliced by distance from centroid of X")
  }
  mpr <<- mapper1D(X.dist,filter_values = slice,num_intervals = num_interval,num_bins_when_clustering = num_bins)
  grh <<- graph.adjacency(mpr$adjacency,"undirected")
  if(!dont.plot)plot(grh)
  
  adj.dist <- mpr$adjacency
  adj.dist[adj.dist==0] <- 1000
  diagram <<- ripsDiag(adj.dist,2,100,"arbitrary")
  return(c(ints = num_interval,
           bins = num_bins,
           vers = mpr$num_vertices,
           levs = max(mpr$level_of_vertex),
           nholes=length(diagram[[1]][diagram[[1]][,1]==1,])/3))
}

mapper.holes.estimation <- function(X,slice,range_interval=c(5,20),range_bins=c(5,20)) {
  if (missing(slice)) {
    m <- apply(X, 2, mean)
    slice <- apply(X,1,function(z)sum(mapply(function(x,y)abs(x^2-y^2), m,z))^(1/2))
    print("sliced by distance from centroid of X")
  }
  X.dist <- dist(X)
  
  mapper.holes.stats <<- array(0,dim = c(range_interval[2]-range_interval[1]+1,range_bins[2]-range_bins[1]+1,5))
  for (ints in range_interval[1]:range_interval[2]) {
    for (bins in range_bins[1]:range_bins[2]) {
      mapper.holes.stats[ints-range_interval[1]+1,bins-range_bins[1]+1,] <<- generate.mapper(X,slice,X.dist,num_interval = ints,num_bins = bins,dont.plot = T)
      print(paste("finished",ints,"int,",bins,"bins"))
    }
  }
  dimnames(mapper.holes.stats) <<- list(paste0("int",range_interval[1]:range_interval[2]),paste0("bin",range_bins[1]:range_bins[2]),c("ints","bins","vers","levs","holes"))
}

show.overlap <- function(X,from.level=F){
  require(Rtsne)
  plot(1,type = "n",xlim = c(1,11),ylim = c(1,9),ylab = "y",xlab = "x")
  node <- if(from.level){mpr$level_of_vertex}else{1:mpr$num_vertices}
  belong <- if(from.level){mpr$points_in_level}else{mpr$points_in_vertex}
  X.tsne <- Rtsne(X)
  plot(X.tsne,type = "n")
  for (n in node) {
    text(X.tsne[belong[[n]],1],X.tsne[belong[[n]],2],labels = n,xlim = c(1,11),ylim = c(1,9),col=rainbow(length(node))[n])
  }
}