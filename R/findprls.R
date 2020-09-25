#' @title MS Lesion Identification and PRL Classification with Pre-trained Model
#' @description This function takes in a lesion probability map, a lesion segmentation mask, and a T2*-phase image for a single subject and identifies + classifies lesions as PRL or not.
#' @param probmap Lesion probability map. We recommend using lesion segmentation algorithm MIMoSA.
#' @param lesmask Lesion segmentation mask. Given a probability threshold, automatically binarizes lesion probability map into a segmentation mask.
#' @param phasefile Location of the T2*-phase image
#' @param predmodel Model used to generate classifications. This automatically uses a pre-trained model ported with the package.
#' @param disc Calculate discretized versions of first order radiomic features? 
#' @export
#' @import lesiontools
#' @import RIA
#' @import Rfast
#' @return A list with 3 objects: a NIfTI file of a lesion label map, a dataframe containing all radiomic features, and a vector containing lesion-wise probability of being a PRL

findprls = function(probmap, lesmask = probmap > thresh, phasefile, 
                    predmodel = pretrainedmodel, disc = T){
  # run lesion identification code
  lesident <- lesion_identification(probmap = probmap,
                                    lesmask = lesmask)
  
  leslabels.out = lesident$leslabels
  print("lesion identification done!")
  
  # run ria feature extraction
  ria.obj = extract_ria(phasefile = phasefile, leslabels = leslabels.out, disc = disc)
  ria.df = as.data.frame(ria.obj)
  print("radiomic feature extraction done!")
  
  # rename variable names to match the ones saved in the model 
  names.temp = as.character(names(ria.df))
  #names.temp = sapply(X = names.temp, function(X){gsub(pattern = "orig.", replacement = "", x = X)})
  names.temp = sapply(X = names.temp, function(X){gsub(pattern = "%", replacement = ".", x = X)})
  names(ria.df) = names.temp
  
  return(list(leslabels = leslabels.out, ria.df = ria.df, preds = predict(pretrainedmodel, newdata = ria.df, type = "prob")))
  
}
