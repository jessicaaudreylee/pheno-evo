### these functions process the lists of values
### spit out by NetLogo at each timepoint

#' Histogram from the numeric values from a single cell of a NetLogo output
#' table
#'
#' In Pheno-Evo experiments, it is common to store all the values from all
#' agents for a particular variable at each given timepoint. In the resulting
#' Netlogo data table, those values all go into a single entry of the table at
#' each timepoint. [parsenumeric()] can be used to split up those
#' values and export them into a new dataframe - specifically, it returns a
#' histogram of the values.
#'
#' In general, the normal Pheno-Evo user should not have to use this function.
#' It is called by several other functions used for analyzing Pheno-Evo output
#' data, such as [calc.fun.div()] and
#' [phenotype.histogram()]. For a function that returns the vector of
#' numbers and does not calculate density, see [extractnumeric()]. For
#' a function that returns the values as factors rather than numeric, see
#' [parsefactors.norm()].
#'
#' @param charstring The full character string from one entry of a NetLogo data
#'   table, including brackets and quotation marks.
#'
#' @return A dataframe with a single column named s4, giving population density
#'   estimates across the range of the trait value from 0 and 1.00 in bins of
#'   width 0.01.
#'
#' @export
#'
#' @examples
#' data(PE.ends)
#' PE.ends$degrade.rate[4]
#' parsenumeric(PE.ends$degrade.rate[4])
parsenumeric<-function(charstring){
  s1<-stringr::str_replace(as.character(charstring),'\\[','')
  s2<-stringr::str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]])
  s3[which(s3>1)]<-1
  s3[which(s3<0)]<-0
  s4<-hist(s3, breaks=seq(0, 1, 0.01), plot=F)$density
  #s4<-hist(s3, breaks=seq(0, 1, 0.01), plot=F)$density
  return(data.frame(s4))
}




#' Factors from the character values from a single cell of a NetLogo output
#' table
#'
#' In Pheno-Evo experiments, it is common to store all the values from all
#' agents for a particular variable at each given timepoint. In the resulting
#' Netlogo data table, those values all go into a single entry of the table at
#' each timepoint. [parsenfactors.norm()] can be used to split up
#' those values and export them into a new dataframe as factors, along with
#' information about their frequency in the population. As an example: the list
#' of individual barcodes carried by cells in a Pheno-Evo population.
#'
#' In general, the normal Pheno-Evo user should not have to use this function.
#' It is called by several other functions used for analyzing Pheno-Evo output
#' data, such as [barcode.timecourse()]. For a function that
#' calculates density and returns a dataframe, see [parsenumeric()].
#' For a function that returns a vector of numeric values, see
#' [extractnumeric()].
#'
#' @param charstring The full character string from one entry of a NetLogo data
#'   table, including brackets and quotation marks.
#'
#' @return A dataframe with three columns: s3 = factor levels; Freq = number of
#'   individuals at each factor level; Prop = proportion of the total community.
#'
#' @export
#'
#' @examples
#' data(PE.ends)
#' PE.ends$barcode[4]
#' parsefactors.norm(PE.ends$barcode[4])
parsefactors.norm<-function(charstring){
  s1<-stringr::str_replace(as.character(charstring),'\\[','')
  s2<-stringr::str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]]) # a list of all the barcodes
  s4<-data.frame(table(s3)) # a table of frequencies of all the barcodes (s3 = barcode name)
  s4$s3<-as.character(s4$s3) # make the barcode names into character strings
  s4$prop<-s4$Freq/sum(s4$Freq)
  return(data.frame(s4))
}




#' Extract a set of numeric values from a single cell of a NetLogo output table
#'
#' In Pheno-Evo experiments, it is common to store all the values from all
#' agents for a particular variable at each given timepoint. In the resulting
#' Netlogo data table, those values all go into a single entry of the table at
#' each timepoint. [extractnumeric()] can be used to split up those
#' values and export them as a vector.
#'
#' In general, the normal Pheno-Evo user should not have to use this function.
#' It is called by several other functions used for analyzing Pheno-Evo output
#' data, such as [summarize.endpoint()]. For a function that
#' calculates density and returns a dataframe, see [parsenumeric()].
#' For a function that returns the values as factors rather than numeric, see
#' [parsefactors.norm()].
#'
#' @param charstring The full character string from one entry of a NetLogo data
#'   table, including brackets and quotation marks.
#'
#' @return A vector containing all the numbers that were in the cell of the
#'   table.
#'
#' @export
#'
#' @examples
#' data(PE.ends)
#' PE.ends$degrade.rate[4]
#' extractnumeric(PE.ends$degrade.rate[4])
extractnumeric<-function(charstring){
  s1<-stringr::str_replace(as.character(charstring),'\\[','')
  s2<-stringr::str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]])
  return(s3)
}




#' Interpret the "xydr" variable from Pheno-Evo model runs
#'
#' In Pheno-Evo, xydr is a property of a cell that is a list of the
#' x-coordinate, the y-coordinate, and the degrade rate. These are reported
#' together in order to enable spatial analysis of degrade rate values, which
#' would not otherwise be possible because coordinate data would be separated from degrade rate data.
#' [extractXYdr()] takes all the xydr values reported in one entry
#' (one timepoint of one experiment) of a Pheno-Evo data table and separates
#' them out into a new dataframe with appropriately labeled columns.
#'
#' In general, the normal Pheno-Evo user should not have to use this function.
#' It is called by several other functions used for analyzing Pheno-Evo
#' output data, such as [degraderate.variogram()].
#'
#' @param charstring The full character string from one entry in the xydir column
#'   of a NetLogo data table, including brackets and quotation marks.
#'
#' @return A dataframe with three columns: xcor, ycor, and degrade.rate.
#' Each row represents the data from one cell.
#' @export
#'
#' @examples
#' data(PE.ends)
#' PE.ends$xydr[4]
#' extractXYdr(PE.ends$xydr[4])
extractXYdr<-function(charstring){
  s1<-stringr::str_replace(as.character(charstring),'\\[\\[','')
  s2<-stringr::str_replace(as.character(s1),'\\]\\]','')
  s3<-strsplit(s2, split='\\] \\[')[[1]]
  s4<-data.frame(s3)
  s5<-tidyr::separate(s4, 1, into=c('xcor','ycor','degrade.rate'),
                      sep=' ', remove=T, convert=T)
  return(s5)
}


