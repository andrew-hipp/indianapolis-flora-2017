## do it!
## ah 2016-05-06
## for Indianapolis prairie

require(ape)
require(phylobase)

prairieStep <- 3

if(prairieStep == 1) {

# tree.zanne <- read.tree('../DATA/phylo.zanne.tre')
spp.mat <- read.delim('../DATA/spp.2016-05-06.tsv', as.is = TRUE)
which.use <- which(!as.character(spp.mat$exclude) %in% '1')

source('../SCRIPTS/zzz.makeMat.R')

indi.tree.v1 <- make.matAndTree(tree.zanne, mat = spp.mat[which.use, ], name.column = 'species',matName = 'indianapolis.matrix.2016-05-06b.csv')

stop('Done with step 1; edit tree matrix that was output, line 5 of 00.treeMaking.R, and go on to step 2!')

}

if(prairieStep == 2) {
  source('../SCRIPTS/zzz.weldTaxa.R')
  indi.mat.edited <- read.csv('indianapolis.matrix.2016-05-06b.csv', as.is = T, row.names = 1)
  indi.tree.edited <-weldTaxa(indi.tree.v1, indi.mat.edited)
  indi.tree.v2 <- drop.tip(indi.tree.edited, which(!indi.tree.edited$tip.label %in% row.names(indi.mat.edited)))
  write.tree(indi.tree.v2, 'indianapolis.tree.pruned.2016-05-06.tre')
}

if(prairieStep == 3) {
  indi.tree.v3 <- drop.tip(indi.tree.v2, c('Quercus_muehlenbergii', 'Quercus_imbricaria'))
  indi.mat.fixers <- read.csv('indianapolis.matrix.2016-05-06c-fixes.csv', as.is = T, row.names = 1)
  indi.tree.v3 <- weldTaxa(indi.tree.v3, indi.mat.fixers, tr = indi.tree.v3)
  indi.tree.v3 <- drop.tip(indi.tree.v3, which(!indi.tree.v3$tip.label %in% row.names(indi.mat.fixers)))
  write.tree(indi.tree.v3, '../DATA/indianapolis.tree.pruned.v3.2016-05-06.tre')
}
