library(mcgraph)

for (n in c(6,12,18)) {
    ang = mcg.angie(nodes=n,edges=n*1.5)

    if (nrow(ang) != n) {
        stop(sprintf("mcg.angie(nodes=%i,edges=%i) does not produce %i nodes",
                     n,n*1.5,n))
    }

    if (sum(ang)/2 != n*1.5) {
        stop(sprintf("mcg.angie(nodes=%i,edges=%i) does not produce %i edges",
                     n,n*1.5,n))
    }
}

# check mcg.shortest.path
ang=mcg.angie(nodes=3,edges=3)
if (length(which(mcg.shortest.paths(ang)==1)) != 6) {
    stop("Error: mcg.shortest.paths for 3 clique not correctly calculated!")
}

if(!max(mcg.shortest.paths(mcg.cluster(nodes=18,edges=12,cluster=2)))==Inf) {
   stop("Error: mcg.shortest.paths for clusters not correctly Inf!")
}

# check mcg.components
if (length(table(mcg.components(mcg.cluster(nodes=10,edges=15,cluster=2)))) != 2) {
    stop("Error: mcg.components not giving two components for two clusters!")
}

if (length(table(mcg.components(mcg.cluster(nodes=12,edges=15,cluster=3)))) != 3) {
    stop("Error: mcg.components not giving two components for three clusters!")
}
