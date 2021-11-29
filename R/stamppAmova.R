#' Analysis of Molecular Variance
#'
#' Calculates an AMOVA based on the genetic distance matrix from stamppNeisD() using the amova() function from the package PEGAS for exploring within and between population variation
#'
#' @param dist.mat the matrix of genetic distances between individuals generated from stamppNeisD()
#' @param geno a data frame containing allele frequency data generated from stamppConvert, or a genlight object containing genotype data, individual IDs, population IDs and ploidy levels
#' @param perm the number of permutations for the tests of hypotheses
#' @details Uses the formula distance ~ populations, to calculate an AMOVA for population differentiation and within & between population variation.
#' This function uses the amova function from the PEGAS package.
#' @return An object of class "amova" which is a list containing a table of sum of square deviations (SSD), mean square deviations (MSD) and the number of degrees of freedom as well as the variance components
#' @examples
#' # import genotype data and convert to allele frequecies
#' data(potato.mini, package="StAMPP")
#' potato.freq <- stamppConvert(potato.mini, "r")
#' # Calculate genetic distance between individuals
#' potato.D.ind <- stamppNeisD(potato.freq, FALSE, "standard")
#' # Calculate AMOVA
#' stamppAmova(potato.D.ind, potato.freq, 100)
#' @author Luke Pembleton <lpembleton at barenbrug.com>
#' @references Paradis E (2010) pegas: an R package for population genetics with an integrated-modular approach. Bioinformatics 26, 419-420. <doi:10.1093/bioinformatics/btp696>
#' @import adegenet pegas
#' @importFrom stats as.dist
#' @export
stamppAmova <- function(dist.mat, geno, perm=100){

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

  pop.names <- geno[,2]

  pop.names <- factor(pop.names) #updated line for compatibility with pegas 0.6

  dist.mat <- as.dist(dist.mat)

  res <- amova(dist.mat ~ pop.names, nperm=perm)

  return(res)

}
