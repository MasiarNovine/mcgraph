#!/usr/bin/env Rscript

usage <- function () {
    cat("graph2data version 0.5\n")
    cat("Usage: graph2data.R [-i 10] [-n 100] adj.tab [out.tab]\n")
    cat("    -h, --help: display this help page\n")
    cat("    -i val (optional): number of iterations (default: 15)\n")
    cat("    -d (optional): debug mode, displays final correlations (default: off)\n")    
    cat("    -n val (optional): number of values (default: 100)\n")
    cat("    -s val (optional): random seed, default 0, which means no seed\n")
    cat("    adj.tab an adjacency matrix of a directed graph, example see below.\n")
    cat("    out.tab (optional) output filename (default: stdout)\n")
    cat("\nExample matrix:\n\n")
    cat("	A	B	C	D	E	F	G
A	0	0	1	0	0	0	0
B	0	0	1	0	0	0	0
C	0	0	0	1	0	0	0
D	0	0	0	0	1	1	0
E	0	0	0	0	0	1	0
F	0	0	0	0	0	0	0
G	-1	0	0	0	0	0	0
\n")
    cat("\nAuthor:  Detlef Groth, University of Potsdam, 2020\n")
    cat("License: MIT\n")
    
    q()
}
# arguments:
# A: adjacency matrix for directed graph, possibly weighted
# n: number of samples for each variable/node
# iter: number of iterations, default: 15 which produces on overage |r|~=0.5, 
#       the higher iter the higher the strength of the associations
# sd: standard deviation for the rnorm function, default: 2
# val: normal distribution mean for the start value, default: 100
# prop: proportion of the target node value take from the source node
# noise: sd for the noise value added after each iteration using rnorm 
#        function with mean 0
# Hint: play mainly with the iter value to increase and decrease the 
#       strength of the associations
#
mcg.graph2data <- function (A,n=100,iter=15,val=100,sd=2,prop=0.1,noise=1) {
    res=matrix(0,ncol=n,nrow=nrow(A))
    rownames(res)=rownames(A)
    for (i in 1:n) {
        units=rnorm(nrow(A),mean=val,sd=sd)
        names(units)=rownames(A)
        for (j in 1:iter) {
            for (node in sample(rownames(A))) {
                targets=colnames(A)[which(A[node,]!=0)]
                for (target in sample(targets)) {
                    if (TRUE) {
                        P=abs(A[node,target])
                        nval=units[[node]]*(prop*P)
                        nval=nval+units[[target]]*(1-(prop*P))
                        if (A[node,target]<0) {
                            diff=nval-units[[target]]
                            nval=units[[target]]-diff
                        }
                        units[[target]]=nval;
                    } else {
                        units[[target]]=units[[target]]+(units[node]*prop*A[node,target])
                    }
                                        
                }
            }
            units=units+rnorm(length(units),sd=noise)
        }
        
        res[,i]=units
    }
    return(res)
}
# wrapp function for the library method above
graph2data <-  function (infile,iter=15,n=100,noise=1) {
    M=read.table(infile,sep="\t",header=TRUE,row.names=1)
    res=mcg.graph2data(M,n=n,iter=iter,noise=noise)
    return(t(res))
}

main <- function (argv) {
    iter=10
    infile=""
    outfile=""
    n=100
    iter=15
    noise=1
    debug=FALSE
    seed=0
    if (length(argv) == 0) {
        usage()
    }
    while (length(argv) > 0) {
        if (argv[1] %in% c("-h","--help")) {
            usage()
        }
        else if (argv[1] %in% c("-d","--debug")) {
            debug=TRUE
            argv=argv[2:length(argv)]
        }
        else if (argv[1] == "-s") {
            if (length(argv)> 1 && grepl("^[0-9]+$",argv[2])) {
                seed=as.numeric(argv[2])
                argv=argv[3:length(argv)]
            } else {
                cat("Error: missing numerical operand for argument -s\n")
                usage()
            }
        }
        else if (argv[1] == "-i") {
            if (length(argv)> 1 && grepl("^[0-9]+$",argv[2])) {
                iter=as.numeric(argv[2])
                argv=argv[3:length(argv)]
            } else {
                cat("Error: missing numerical operand for argument -i\n")
                usage()
            }
        }
        else if (argv[1] == "-n") {
            if (length(argv)> 1 && grepl("^[0-9]+$",argv[2])) {
                n=as.numeric(argv[2])
                argv=argv[3:length(argv)]
            } else {
                cat("Error: missing numerical operand for argument -n\n")
                usage()
            }
        }
        else if (argv[1] == "-r") {
            if (length(argv)> 1 && grepl("^[\\.0-9]+$",argv[2])) {
                noise=as.numeric(argv[2])
                argv=argv[3:length(argv)]
            } else {
                cat("Error: missing numerical operand for argument -r\n")
                usage()
            }
        }

        else if (infile == "" & !file.exists(argv[1])) {
            cat(paste("Error: not existing file",argv[1],"!\n",sep=""))
            usage()
        } 
        else if (infile == "") {
            infile = argv[1]
            if (length(argv) ==2) {
                outfile=argv[2]
            } else {
                outfile = ""
            }
            break
        }
    }
    if (seed>0) {
        set.seed(seed)
    }
    data=graph2data(infile,iter=iter,n=n,noise=noise)
    data=round(data,4)
    if (outfile == "") {
        print(data)
    } else {
        write.table(data,file=outfile,sep="\t",quote=FALSE)
        cat(paste("Outfile ",outfile,"written!\n",sep=" "))
        if (debug) {
            print(summary(data))
            print(round(cor(data),1))
        }
            
    }
}
if (sys.nframe() == 0L && !interactive()) {
    main(commandArgs(trailingOnly=TRUE))
}
