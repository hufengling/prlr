#' @title MS Lesion Identification and PRL Classification with Pre-trained Model
#' @description This function takes in a lesion probability map, a lesion segmentation mask, and a T2*-phase image for a single subject and identifies + classifies lesions as PRL or not.
#' @param phase_path Filepath of T2*-phase image with appropriate preprocessing already done. Cannot be antsImage or nifti object.
#' @param prob_map Lesion probability map. We recommend using lesion segmentation algorithm MIMoSA. Can be path to .nii.gz, antsImage or nifti object
#' @param bin_map Lesion segmentation mask obtained by thresholding prob_map. Can be path to .nii.gz, antsImage or nifti object
#' @param register_maps Logical indicating whether prob_map and bin_map need to be registered to the phase image. (Default = F)
#' @param prob_map_image If register_maps == T, an image in the prob_map/bin_map space must be provided
#' @param epi_image If register_maps == T, an EPI image in the same space as the phase image must be provided. Registration to phase does not work well due to noise.
#' @param verbose Logical indicating whether print statements should be made (Default = F)
#'
#' @export
#'
#' @import RIA
#' @import Rfast
#' @import ANTsRCore
#' @import extrantsr
#' 
#' @return A list with 2 objects: a NIfTI file of a lesions where each lesion is labeled with its probability of being a PRL and a dataframe containing all radiomic features.

findprls <- function(phase_path, prob_map, bin_map,
                     register_maps = F,
                     prob_map_image = NULL, epi_image = NULL,
                     verbose = F) {
  prob_map <- check_ants(prob_map)
  bin_map <- check_ants(bin_map)
  
  if (sum(bin_map) == 0) {
    warning("No lesions detected")
    return(NULL)
  }
  if (class(phase_path) != "character") {
    stop("phase_path must be a string path to the phase image.")
  }
  
  if (register_maps) {
    if (is.null(prob_map_image) | is.null(epi_image)) {
      stop("If you need APRL to register the prob_map/bin_map, both raw_flair and raw_epi need to be provided.")
    }
  }
  
  lesion_mask <- label_lesion(prob_map, bin_map, mincluster = 30)
  
  if (register_maps) {
    if (verbose) {
      print("Registering prob_mask to epi space")
    }
    prob_map_image <- check_ants(prob_map_image)
    epi_image <- check_ants(epi_image)
    
    transform <- antsRegistration(epi_image, prob_map_image,
                                  typeofTransform = "Rigid"
    )
    lesion_mask <- antsApplyTransforms(
      fixed = epi_image, moving = lesion_mask,
      transformlist = c(transform$fwdtransforms),
      interpolator = "nearestNeighbor"
    )
    lesion_mask <- antsCopyImageInfo(check_ants(phase_path), lesion_mask)
    lesion_mask <- ants2oro(lesion_mask)
  }
  
  # run ria feature extraction
  ria.obj <- extract_ria(
    phase_path = phase_path,
    lesion_labels = lesion_mask,
    disc = F
  )
  ria.df <- as.data.frame(ria.obj)
  
  # rename variable names to match the ones saved in the model
  names.temp <- as.character(names(ria.df))
  names.temp <- sapply(X = names.temp, function(X) {
    gsub(pattern = "%", replacement = ".", x = X)
  })
  names(ria.df) <- names.temp
  
  pretrainedmodel <- prlr::prlmodel_orig
  
  predictions <- stats::predict(pretrainedmodel,
                                newdata = ria.df,
                                type = "prob")
  for (j in 1:nrow(predictions)) {
    lesion_mask[lesion_mask == j] <- predictions["rimpos"][j]
  }
  
  return(list(
    leslabels = lesion_mask,
    ria.df = ria.df
  ))
}
