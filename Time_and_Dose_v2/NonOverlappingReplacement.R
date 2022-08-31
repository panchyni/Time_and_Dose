NonOverlappingReplacement <- function(check){
  # Go through suggested symbols and replaced the suggested
  # symbol with the exist symbol if the suggested symbol already exists
  # because duplicate name are not allowed and visiaul inspection suggests
  # at least some of these same symbol genes have different expression patterns
  
  # Make a table of all Suggested Symbols looking for dupliates
  check_suggested_table <- table(check$Suggested.Symbol)
  
  for (i in 1:dim(check)[1]){
    #print(i)
    
    # If this is a replacement row
    if(check[i,]$Approved == 'FALSE'){
      
      # If its actually an new symbol
      if (check[i,]$Suggested.Symbol != check[i,]$x){
        # If the suggested symbols is duplicated
        if (check_suggested_table[check[i,]$Suggested.Symbol] > 1){
          # Replace the suggested symbols with the existing symbol
          check[i,]$Suggested.Symbol <- check[i,]$x
          #print(check[i,]$x)
        }
      }
    }
  }
  return(check)
}

### Check code ###

check_NonOverlappingReplacement <- function(){
  
  # Test 1
  dummy_check <- rbind( c("A","TRUE","A"),    
                        c("B","TRUE","B"),
                        c("C","FALSE","D"),
                        c("E","FALSE","F"),
                        c("G","FALSE","H"),  # This shoule be replace as a dup
                        c("G2","FALSE","H"), # This shoule be replace as a dup
                        c("I","TRUE","I"),
                        c("J","FALSE","I"), # This shoule be replace as a dup
                        c("K","TRUE","K"),
                        c("L","FALSE","K"), # This shoule be replace as a dup
                        c("M","FALSE","K"),
                        c("N","FALSE","A"),
                        c("O","FALSE","O"), # This should be kept, even if false, because not replacement
                        c("P","FALSE","P"), # This should be kept, even if false, because not replacement
                        c("Q","FALSE","P")) # This shoule be replace as a dup
  dummy_check <- as.data.frame(dummy_check)
  colnames(dummy_check) <- c("x","Approved","Suggested.Symbol")
  dummy_check_NoOverlap <- NonOverlappingReplacement(dummy_check)
  #dummy_check_NoOverlap$x
  print(dummy_check_NoOverlap$Suggested.Symbol)
  # Should read A, B, D, F, G, G2, I, J, K, L, M, N, O, P, Q
  
  # Check 2 --- Real Gene Names + Dummies
  # Update gene map
  library(HGNChelper)
  currentmamp <- getCurrentHumanMap()
  
  MarkerGenes <- readRDS("MarkerGenes.rds")
  check_markers <- checkGeneSymbols(MarkerGenes,unmapped.as.na=FALSE,map=currentmamp)
  table(check_markers$Approved)
  
  print(check_markers[check_markers$Approved=="FALSE",])
  
  check_markers_NoOverlap <- NonOverlappingReplacement(check_markers)
  print(check_markers_NoOverlap[check_markers_NoOverlap$Approved=="FALSE",]) # SHould be same as above
  
  # Add in dummy genes named
  # ATP5MK "TRUE" ATP5MK
  # ATP5MPL2 "FALSE" ATP5MJ
  
  check_markers_dummy <- rbind(check_markers,c("ATP5MK",TRUE,"ATP5MK"),c("ATP5MPL2",FALSE,"ATP5MJ"))
  colnames(check_markers_dummy) <- c("x","Approved","Suggested.Symbol")
  print(check_markers_dummy[check_markers_dummy$Approved=="FALSE",])
  
  check_markers_dummy_NoOverlap <-NonOverlappingReplacement(check_markers_dummy)
  print(check_markers_dummy_NoOverlap[check_markers_dummy_NoOverlap$Approved=="FALSE",]) # Should revert the bottom 3 genes
}
