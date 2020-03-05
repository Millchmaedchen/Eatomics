# Experimental setup functions to remove non-overlapping subjects in clinical and in protein abundance
matchedSamples<- function(expDesign, proteinAbundance){
  #expDesign = expDesign %>% rownames_to_column("PatientID")
  l = proteinAbundance %>% colnames()
  r = as.character(expDesign$PatientID)
  Reduce(intersect, list(l, r))
}

matchedExpDesign <- function(expDesign, proteinAbundance){
  expDesign = expDesign %>% rownames_to_column("PatientID")
  proteinAbundance = proteinAbundance %>% column_to_rownames("Gene names") %>% as.data.frame()
  matchedExpDesign <- expDesign %>% filter(PatientID %in% matchedSamples(expDesign, proteinAbundance))
  matchedExpDesign = matchedExpDesign %>% column_to_rownames("PatientID")
  return(matchedExpDesign)
}
