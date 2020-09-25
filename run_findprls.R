
library(fslr)

source("/lesion_identification.R")
source("/extract_ria.R")
source("/findprls.R")
source("/getleslabels.R")

probmap = readnii('/probmap.nii.gz')
lesmask = readnii('/lesmask.nii.gz')
phasefile = '/phasefile'
pretrainedmodel = readRDS('/pretrainedmodel.RDS')

findprls_out = findprls(probmap = probmap, lesmask = lesmask, phasefile = phasefile, predmodel = pretrainedmodel)
print(unique(findprls_out$preds))

saveRDS(findprls_out, "/out/findprls_out_disc.RDS")