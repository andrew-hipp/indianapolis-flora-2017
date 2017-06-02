## visualizing and stats

library(phytools)
library(caper)
library(ggtree)

read.dat <- TRUE
plot.supp <- TRUE
d.tests <- TRUE
plot.figs <- TRUE

if(read.dat) {
spp <- read.delim('../DATA/spp.2016-09-09.tsv', as.is = TRUE)
which.use <- which(!as.character(spp$exclude) %in% '1')

traits <- read.delim('../DATA/traits.2016-09-11.tsv', as.is = T, row.names = 4)
row.names(traits) <- gsub(" ", "_", row.names(traits))
traits$spp <- row.names(traits)
ind.tree <- read.tree('../RESULTS/R.OUT/indianapolis.tree.pruned.v3.2016-05-06.tre')

ind.tree <- drop.tip(ind.tree, which(!ind.tree$tip.label %in% row.names(traits)))
ind.tree$node.label[ind.tree$node.label %in% c('NA', '')] <- paste('node', seq(sum(ind.tree$node.label %in% c('NA', ''))))
orders <- c('Poales', 'Asparagales', 'Asterales', 'Dipsacales', 'Gentianales',
            'Ericales', 'Caryophyllales', 'Rosales', 'Fabales', 'Brassicales', 'Malvales', 'Magnoliales')
ind.tree$node.label.orders <- ind.tree$node.label
ind.tree$node.label.orders[which(!ind.tree$node.label.orders %in% orders)] <- NA

traits <- traits[which(row.names(traits) %in% ind.tree$tip.label), ]
traits$Status <- ifelse(traits$Status == "0", 0, 1)

traits$losers <- as.integer(traits$MC.Hist) == 1 & as.integer(traits$MC.Now) == 0
traits$winners <- as.integer(traits$MC.Hist) == 1 & as.integer(traits$MC.Now) == 1
traits$incomers <- as.integer(traits$MC.Hist) == 0 & as.integer(traits$MC.Now) == 1

growth.form <- c('forb', 'fern', 'graminoid', 'shrub', 'tree',
                 'succulent', 'epiphyte', 'vine', 'parasite',
				 'aquatic', 'other', 'missing')
growth.lumped <- c('Other', 'Other', 'Graminoid', 'Woody', 'Woody',
                   'Other', 'Other', 'Other', 'Other',
				   'Other', 'Other', 'Other')
growth.col <- structure(c('yellow', 'brown', 'gray'), names = c('Graminoid', 'Woody', 'Other'))

traits$growthLumped <- traits$Growth.form
traits$growthLumped[grep('4', traits$Growth.form)] <- '4'
traits$growthLumped[grep(',', traits$Growth.form, fixed = TRUE)] <- '1'
traits$growthLumped <- growth.lumped[as.integer(traits$growthLumped)]
traits$growthLumped[traits$growthLumped == 'Other'] <- NA

traits$mc.text <- ifelse(traits$MC.Now == 1, "Currently in Indianapolis", NA)

traits$pollinationText <- ''
traits$pollinationText[traits$Pollination.System == 1] <- NA
traits$pollinationText[!is.na(traits$pollinationText)] <- 'Abiotic pollination'

traits$graminoid <- ifelse(traits[, 'growthLumped'] == "Graminoid", 1, 0)
traits$graminoid[is.na(traits$graminoid)] <- 0

traits$woody <- ifelse(traits[, 'growthLumped'] == "Woody", 1, 0)
traits$woody[is.na(traits$woody)] <- 0

traits$abioticPollination <- ifelse(is.na(traits[, 'pollinationText']), 0, 1)



ind.native <- drop.tip(ind.tree, row.names(traits)[which(traits$Status == 1)])
traits.native <- traits[ind.native$tip.label, ]

}

if(plot.supp) {

pdf('SUPPL.1.indianapolis.tree.pdf', 30, 30)
plot(ind.tree, 'fan', cex = 0.4)
nodelabels(grep('ales', ind.tree$node.label, value = T), grep('ales', ind.tree$node.label, value = F) + length(ind.tree$tip.label), cex = 1, bg = 'white')
dev.off()

pdf('SUPPL.2.trial.plots.pdf', 11, 8.5)
layout(matrix(1:4, 2,2))
plotTree.wBars(ind.tree, structure(traits$MC.Hist, names = row.names(traits)), lwd = 0.5, type = 'fan', scale = 10,  mar = c(1,1,2,1))
nodelabels(grep('ales', ind.tree$node.label, value = T), grep('ales', ind.tree$node.label, value = F) + length(ind.tree$tip.label), cex = 0.5, bg = 'white')
title(main = "MC historical")
plotTree.wBars(ind.tree, structure(traits$MC.Now, names = row.names(traits)), lwd = 0.5, type = 'fan', scale = 10,  mar = c(1,1,2,1))
nodelabels(grep('ales', ind.tree$node.label, value = T), grep('ales', ind.tree$node.label, value = F) + length(ind.tree$tip.label), cex = 0.5, bg = 'white')
title(main = "MC now")
plotTree.wBars(ind.tree, structure(traits$region, names = row.names(traits)), lwd = 0.5, type = 'fan', scale = 10, mar = c(1,1,2,1))
nodelabels(grep('ales', ind.tree$node.label, value = T), grep('ales', ind.tree$node.label, value = F) + length(ind.tree$tip.label), cex = 0.5, bg = 'white')
title(main = "Region")
plotTree.wBars(ind.tree, structure(traits$Status, names = row.names(traits)), lwd = 0.5, type = 'fan', scale = 10, mar = c(1,1,2,1))
nodelabels(grep('ales', ind.tree$node.label, value = T), grep('ales', ind.tree$node.label, value = F) + length(ind.tree$tip.label), cex = 0.5, bg = 'white')
title(main = "Native status")
dev.off()

}

if(d.tests) {

ind.native.compv1 <- comparative.data(ind.native, traits.native[c('region', 'MC.Hist', 'MC.Now', 'spp', 'losers', 'winners', 'incomers', 'graminoid', 'woody', 'abioticPollination')], 'spp')
logfile <- file(format(Sys.time(), 'd.tests.%Y-%m-%d.txt'), 'wa')
writeLines('Historical flora, natives', con = logfile)
out.native.hist.d <- phylo.d(ind.native.compv1, binvar = MC.Hist)

writeLines('\n\nCurrent flora, natives', con = logfile)
out.native.now.d <- phylo.d(ind.native.compv1, binvar = MC.Now)
close(logfile)

ind.all.compv1 <- comparative.data(ind.tree, traits[c('region', 'MC.Hist', 'MC.Now', 'winners', 'losers', 'incomers', 'spp', 'graminoid', 'woody', 'abioticPollination')], 'spp')
ind.all.compv2 <- comparative.data(ind.tree, traits[c('region', 'MC.Hist', 'MC.Now', 'winners', 'losers', 'incomers', 'Status', 'spp')], 'spp')
out.all.losers <- phylo.d(ind.all.compv1, binvar = losers)
out.all.winners <- phylo.d(ind.all.compv1, binvar = winners)
out.all.incomers <- phylo.d(ind.all.compv1, binvar = incomers)
out.all.native <- phylo.d(ind.all.compv2, binvar = Status)

out.native.losers <- phylo.d(ind.native.compv1, binvar = losers)
out.native.winners <- phylo.d(ind.native.compv1, binvar = winners)
out.native.incomers <- phylo.d(ind.native.compv1, binvar = incomers)

out.all.hist.d <- phylo.d(ind.all.compv1, binvar = MC.Hist)
out.all.now.d <- phylo.d(ind.all.compv1, binvar = MC.Now)

out.native.woody.d <- phylo.d(ind.native.compv1, binvar = woody)
out.native.graminoid.d <- phylo.d(ind.native.compv1, binvar = graminoid)
out.native.abioticPollination.d <- phylo.d(ind.native.compv1, binvar = abioticPollination)

out.all.woody.d <- phylo.d(ind.all.compv1, binvar = woody)
out.all.graminoid.d <- phylo.d(ind.all.compv1, binvar = graminoid)
out.all.abioticPollination.d <- phylo.d(ind.all.compv1, binvar = abioticPollination)

write.csv(traits[traits$incomers == 1 & traits$Status == 0, ], 'newcomers.2016.05.06.csv')
write.csv(traits[traits$incomers == 1 & traits$Status == 0, ], 'newcomers.native.2016.05.06.csv')
write.csv(traits[traits$incomers == 1, ], 'newcomers.all.2016.05.06.csv')

}

ind.tree.dat <- match.phylo.data(ind.tree, traits)
ind.tree.mpd <- ses.mpd(t(ind.tree.dat$data[, c('MC.Hist', 'MC.Now')]), cophenetic(ind.tree.dat$phy))
ind.tree.mntd <- ses.mntd(t(ind.tree.dat$data[, c('MC.Hist', 'MC.Now')]), cophenetic(ind.tree.dat$phy))

ind.tree.native.dat <- match.phylo.data(ind.native, traits)
ind.tree.native.mpd <- ses.mpd(t(ind.tree.native.dat$data[, c('MC.Hist', 'MC.Now')]), cophenetic(ind.tree.native.dat$phy))
ind.tree.native.mntd <- ses.mntd(t(ind.tree.native.dat$data[, c('MC.Hist', 'MC.Now')]), cophenetic(ind.tree.native.dat$phy))



if(plot.figs) {
  pdf(format(Sys.time(), 'fig.x.treeWithTraits.%Y-%m-%d.pdf'), width = 7.5, height = 7.5)
  plotColors <- c('Abiotic pollination' = 'burlywood3', 'Currently in Indianapolis' = 'blue', Woody = 'black', Graminoid = 'orange', 'x - All other character states' = 'grey90')
  trt2 <- traits[ind.tree$tip.label, c('growthLumped', 'mc.text', 'pollinationText')]
  for(i in names(trt2)) trt2[[i]][is.na(trt2[[i]])] <- 'x - All other character states'
  p <- ggtree(ind.tree, layout='circular', colour = 'gray')
  p <- gheatmap(p, trt2, width = 0.1, font.size = 0)
  p <- p + geom_label(aes(x=branch, label = c(rep(NA, length(ind.tree$tip.label)), ind.tree$node.label.orders)), fill = 'white', colour = 'black', size = 2)
  p <- p + scale_fill_manual('Traits', values = plotColors)
  p <- p + theme(legend.position = c(0.5, 0.4),
               legend.key.size = unit(0.3, 'cm'))
  print(p)
  dev.off()
}

## rarefy by genus
ind.tree.genera <- sapply(ind.tree$tip.label, function(x) strsplit(x, "_", fixed = T)[[1]][1])
genus.sets <- lapply(1:100, function(w) sapply(unique(ind.tree.genera), function(x) {
  if(sum(ind.tree.genera == x) == 1) out <- which(ind.tree.genera == x)
  else out <- sample(which(ind.tree.genera == x), 1)
  out
  })
)
ind.tree.rarefied <- lapply(genus.sets, function(x) drop.tip(ind.tree, setdiff(seq(length(ind.tree$tip.label)), x)))
ind.tree.rarefied.comp <- lapply(ind.tree.rarefied, function(x) comparative.data(x, traits[c('region', 'MC.Hist', 'MC.Now', 'winners', 'losers', 'incomers', 'spp', 'graminoid', 'woody', 'abioticPollination')], 'spp'))
ind.tree.rarefied.nowD <- lapply(ind.tree.rarefied.comp, phylo.d, binvar = MC.Now)
