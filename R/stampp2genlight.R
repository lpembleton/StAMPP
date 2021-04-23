#####################################################
#
#' Convert StAMPP genotype data to genlight object
#'
#' Converts a StAMPP formated allele frequency data frame generated from the stamppConvert function to a genlight object for use in other packages
#'
#' @param geno a data frame containing allele frequency data generated from stamppConvert
#' @param pop logical. True if population IDs are present in the StAMPP genotype data, False if population IDs are absent.
#' @details StAMPP only exports to genlight objects as they are able to handle mixed ploidy datasets unlike genpop and genloci objects.
#' The genlight object allows the intergration between StAMPP and other common R packages such as ADEGENET
#' @return A object of class genlight which contains genotype data, individual IDs, population IDs (if present) and ploidy levels
#' @examples
#' # import genotype data and convert to allele frequecies
#' data(potato.mini, package="StAMPP")
#' potato.freq <- stamppConvert(potato.mini, "r")
#' # Convert the StAMPP formatted allele frequency data frame to a genlight object
#' potato.genlight <- stampp2genlight(potato.freq, TRUE)
#' @author Luke Pembleton <luke.pembleton at agriculture.vic.gov.au>
#' @import adegenet
#' @importFrom methods new
#' @export
stampp2genlight <- function(geno, pop=TRUE){


      data <- geno[,-(1:5)]   #matrix of allele frequencies
      ind <- geno[,1] #individual ids
      ploidy.levels <- geno[,4] #ploidy
      pop.names <- geno[,2] #population ids

      data=data*ploidy.levels #convert percentage of allele A to number of allele A based on ploidy level

      data <- new("genlight", data, parallel=FALSE) #convert genotype data to genlight object
      indNames(data)=ind #add individual ids to genlight object
      ploidy(data)=ploidy.levels #add ploidy levels to genlight object

      if(pop==TRUE){

        pop(data)=pop.names #if population ids are present, add to genlight object

      }


    return(data)

  }
