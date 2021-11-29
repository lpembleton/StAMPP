#' Genetic Distance Calculation
#'
#' This function calculates Nei's genetic distance (Nei 1972) or Nei's Da distance 1983 between populations or individuals
#'
#' Updated syntax and transition of StAMPP to tidyverse
#'
#' @param geno a data frame containing allele frequency data generated from stamppConvert, or a genlight object containing genotype data, individual IDs, population IDs and ploidy levels
#' @param pop logical. True if genetic distance should be calculated between populations, false if it should be calculated between individual
#' @param measure a character string defining the distance measure to use: "standard" for the Neis standard genetic distance 1972 or "DA" for Neis DA distance 1983.
#' @return A object of class matrix which contains the genetic distance between each population or individual
#' @examples
#' # import genotype data and convert to allele frequecies
#' data(potato.mini, package="StAMPP")
#' potato.freq <- stamppConvert(potato.mini, "r")
#' # Calculate genetic distance between individuals
#' potato.D.ind <- stamppNeisD(potato.freq, FALSE, "standard")
#' # Calculate genetic distance between populations
#' potato.D.pop <- stamppNeisD(potato.freq, TRUE, "standard")
#' @author Luke Pembleton <lpembleton at barenbrug.com>
#' @references Nei M (1972) Genetic Distance between Populations. The American Naturalist 106, 283-292.
#' @import adegenet
#' @importFrom utils combn
#' @export
stamppNeisD <- function(geno, pop=TRUE, measure="standard"){

  if(class(geno)=="genlight"){  #if input file is a genlight object convert to a data.frame

    geno2 <- geno

    geno <- as.matrix(geno2) #extract genotype data from genlight object
    sample <- row.names(geno) #individual names
    pop.names <- pop(geno2) #population names
    ploidy <- ploidy(geno2) #ploidy level
    geno=geno*(1/ploidy) #convert genotype data (number of allele 2) to precentage allele frequency
    geno[is.na(geno)]=NaN
    format <- vector(length=length(geno[,1]))
    format[1:length(geno[,1])]="genlight"


    pops <- unique(pop.names) #population names

    pop.num <- vector(length=length(geno[,1])) #create vector of population ID numbers

    for (i in 1:length(geno[,1])){
      pop.num[i]=which(pop.names[i]==pops) #assign population ID numbers to individuals
    }

    genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, ploidy, format), stringsAsFactors=FALSE)

    geno <- cbind(genoLHS, geno) #combine genotype data with labels to form stampp geno file

    geno[,2]=as.character(pop.names)

    geno[,4]=as.numeric(as.character(geno[,4]))

    row.names(geno)=NULL

  }

  if(pop==TRUE){ #if calculating genetic distance between populations

    p <- split(geno[,-c(1:5)], f=geno[,3])
    p <- matrix(unlist(lapply(p, colMeans, na.rm=T)), ncol=(ncol(geno)-5), byrow=T)
    row.names(p) <- unique(geno[,2])

  }else{ #individual pairwise distances

    p <- as.matrix(geno[,-c(1:5)])
    row.names(p) <- geno[,1]

  }


  if(measure=="standard"){

    idx <- unlist(combn(seq_len(nrow(p)), 2, simplify=F))
    grp_id <- rep(1:(length(idx)/2), each=2)

    p0 <- p
    p0[is.na(p0)] <- 0 #p matrix with NAs converted to zeros as subsequent matrix multiplication will negate NA comparisons within pairs.

    jxy <- tcrossprod(p0, p0) + tcrossprod(1-p0, 1-p0)

    j_u1 <- p^2
    j_u2 <- (1-p)^2 #alt allele
    j_u1 <- j_u1[idx,] #expand to pair combinations
    j_u2 <- j_u2[idx,]

    idx <- order(rep(1:(nrow(j_u1)/2), each=2), rep(2:1, nrow(j_u1)/2))
    j_u1 <- j_u1*!is.na(j_u1[idx,]) #where there is a NA at a locus in a pairwise comparison replace with 0
    j_u2 <- j_u2*!is.na(j_u2[idx,])

    jxjy <- rowSums(j_u1, na.rm=T) + rowSums(j_u2, na.rm=T)
    jxjy <- jxjy[seq(1, length(jxjy), 2)]*jxjy[seq(2, length(jxjy), 2)]


    I <- jxy[lower.tri(jxy)]/sqrt(jxjy)
    D <- -log(I)

    neis_d <- matrix(0, nrow=nrow(p), ncol=nrow(p), dimnames=list(row.names(p), row.names(p)))
    neis_d[lower.tri(neis_d)] <- as.numeric(sprintf("%.6f", D))
    neis_d[upper.tri(neis_d)] <- t(neis_d)[upper.tri(neis_d)]

    return(neis_d)

  }

  if(measure=="DA"){

    idx <- unlist(combn(seq_len(nrow(p)), 2, simplify=F))
    grp_id <- rep(1:(length(idx)/2), each=2)

    p0 <- p
    p0[is.na(p0)] <- 0 #p matrix with NAs converted to zeros as subsequent matrix multiplication will negate NA comparisons within pairs.

    p_u1 <- p
    p_u2 <- (1-p) #freq of the alternative allele
    p_u1 <- p_u1[idx,] #expand to pair combinations
    p_u2 <- p_u2[idx,]

    sqrt_xy_u1 <- rowSums(apply(p_u1, 2, function(x) { sqrt(x[seq(1,length(x), by=2)]*x[seq(2,length(x), by=2)]) }), na.rm=T) #sqrt of Xu*Yu
    sqrt_xy_u2 <- rowSums(apply(p_u2, 2, function(x) { sqrt(x[seq(1,length(x), by=2)]*x[seq(2,length(x), by=2)]) }), na.rm=T)

    p1 <- p
    p1[!is.na(p1)] <- 1 #where not a NA use 1's to sum up number of loci in pairwise comparison
    p1[is.na(p1)] <- 0

    L <- tcrossprod(p1, p1)*2 #number of loci in the pairwise comparison *2 alleles

    DA <- 1-((sqrt_xy_u1+sqrt_xy_u2)/L[lower.tri(L)])

    neis_DA <- matrix(0, nrow=nrow(p), ncol=nrow(p), dimnames=list(row.names(p), row.names(p)))
    neis_DA[lower.tri(neis_DA)] <- as.numeric(sprintf("%.6f", DA))
    neis_DA[upper.tri(neis_DA)] <- t(neis_DA)[upper.tri(neis_DA)]

    return(neis_DA)


  }

}
