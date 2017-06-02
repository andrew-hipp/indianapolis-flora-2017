## Make matrix for binding taxa to a known tree
## ahipp@mortonarb.org, 2016-05-06

require(ape)

make.matAndTree <- function(tr, mat, rows.include = 'all', columns.include = 'all', name.column = 'species', trim2sp = T, matName = NA) {
  ## format the matrix
  if(identical(rows.include, 'all')) rows.include <- seq(dim(mat)[1])
  if(identical(columns.include, 'all')) columns.include <- seq(dim(mat)[2])
  mat.temp <- mat[rows.include, columns.include]

  ## do the species tree
  tr$tip.label <- gsub(" ", "_", tr$tip.label)
  taxa <- gsub(" ", "_", mat.temp[[name.column]])
  if(trim2sp) taxa <- label.elements(taxa, "_", 1:2, "_", fixed = T)
  taxa <- unique(taxa)
  tr.out.spp <- drop.tip(tr, tr$tip.label[!tr$tip.label %in% taxa])
  taxa.matched <- taxa[taxa %in% tr.out.spp$tip.label]
  taxa.not.matched <- taxa[!taxa %in% tr.out.spp$tip.label]

  ## do the genus tree
  genera.missing <- unique(sapply(taxa.not.matched, function(x) strsplit(x, "_")[[1]][1]))
  genera.tr.tips <- sapply(tr$tip.label, function(x) strsplit(x, "_")[[1]][1])
  include.genus.tree <- unique(c(taxa.matched, tr$tip.label[genera.tr.tips %in% genera.missing]))
  tr.out.genera <- drop.tip(tr, tr$tip.label[!tr$tip.label %in% include.genus.tree])

  ## make the renaming matrix
  mat.rename <- matrix('', nrow = length(taxa), ncol = 3, dimnames = list(taxa, c('seqOrSplice', 'spliceRule', 'spliceTaxa')))
  mat.rename[taxa.matched, 'seqOrSplice'] <- 'seq'
  mat.rename[taxa.not.matched, 'seqOrSplice'] <- '** SPLICE **'
  mat.rename[taxa.not.matched, 'spliceRule'][sapply(taxa.not.matched, function(x) strsplit(x, "_")[[1]][1]) %in% genera.tr.tips] <- 'genus'
  mat.rename[taxa.not.matched, 'spliceRule'][mat.rename[taxa.not.matched, 'spliceRule'] != 'genus'] <- '???'

  ## then spit everything back
  if(!is.na(matName)) write.csv(mat.rename, matName)
  out <- list(tr.spp = tr.out.spp, tr.genera = tr.out.genera, matrix.full = mat, matrix.renaming = as.data.frame(mat.rename), missing.spp = taxa.not.matched, matrix.rows.used = rows.include, name.used = name.column, timestamp = date())
  class(out) <- 'prairie.match.step1'
  message("You'll probably want to save the 'mat.rename' object\nof the list you've just been handed so you can edit\nin your favorite spreadsheet software")
  return(out)
  }

label.elements <- function(x, delim = '|', returnNum = 1, returnDelim = ' ', ...) {
  ## finds any label at the tips; default assumes pipe delimitation, and the element of interest is b/f the first pipe
    if('phylo' %in% class(x)) labelVector <- x$tip.label
    else labelVector <- x
    out <- sapply(labelVector, function(x) paste(strsplit(x, delim, ...)[[1]][returnNum], collapse = returnDelim))
    out
    }
