# prlr
This R package implements the method described in Lou (2020): https://www.biorxiv.org/content/10.1101/2020.08.31.276238v1. This package takes in a lesion probability map, a lesion segmentation mask, and the filepath to a phase image. It then identifies distinct lesions and provides PRL classifications for each identified lesion based on a pretrained model.

The primary function in this package is the function `findprls`. A sample script is provided below.

```
library(fslr)

probmap = readnii('/probmap.nii.gz')
lesmask = readnii('/lesmask.nii.gz')
phasefile = '/phasefile'
pretrainedmodel = prlr::prlmodel_orig

findprls_out = findprls(probmap = probmap, lesmask = lesmask, phasefile = phasefile, disc = T)

leslabels_img = findprls_out$leslabels
ria.df = findprls_out$ria.df
preds = findprls_out$preds

saveRDS(findprls_out, "/out/findprls_out_disc.RDS")
```
