## code for binding missing taxa
## A Hipp, for Larkin et al. paper
## updated 2015-06-24 for prairie experiment, AH

require(phytools)
require(ape)
require(phylobase)

weldTaxa <- function(dat = NA, taxa = NA, tr = NA, branchScalar = 0.50, logFileBase = 'prairie.matchv2.weldingLog') {
## tr is a phylo-style tree
## taxa is either a vector of names, with genus separated from sp. by "_", or a list of length three, in which:
##     $vect is a vector of taxa to weld in, and
##     $spliceRule is a vector of binding methods, as described under $cladeMembers below
##     $spliceTaxa is a vector of:
##        "mrca" pipe-delimited taxa whose MRCA is the node to which the corresponding taxon in $vect should be welded
##        "sister" or "substitute" a single taxon to which the taxon will be welded at branchScalar proportional distance from subtending node
##        "genus" nothing; the genus will be used to identify the attachment node
##		  "drop" means to drop the taxon

  if(!class(dat) %in% 'prairie.match.step1') warning('This function expects output from do.matAndTree.step1')
  if(is.na(taxa[1])) {
    warning('no taxa handed in; using taxon matrix from dat')
	taxa <- dat$matrix.renaming
	}
  nTaxa <- 0
  if(class(taxa) %in% c("matrix", "data.frame", "list")) {
    whichTaxa <- taxa$seqOrSplice == "** SPLICE **"
    taxa <- list(vect = as.character(row.names(taxa)[whichTaxa]),
                 bindMethod = as.character(taxa$spliceRule[whichTaxa]),
                 cladeMembers = as.character(taxa$spliceTaxa[whichTaxa])
    )
  }
  if(class(taxa) == "list") nTaxa <- length(taxa[[1]])
  if(class(taxa) == "character") nTaxa <- length(taxa)
  if(nTaxa == 0) stop("Something is awry with your taxon list")

  if(is.na(tr[[1]])) tr <- dat$tr.genera
  trList <- vector('list', nTaxa)
  weldNodeVector <- numeric(nTaxa)
  ## make a tree for each taxon
  lastCount <- 0
  out <- file(paste(logFileBase,'_', format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), '.txt', sep = ""), open = 'a')

  for(i in 1:nTaxa) {
    #trDists <- dist.nodes(tr)
    newCount <- length(tr$tip.label)
    if(newCount == lastCount) {
      writeOut <- paste("**", taxa$vect[i-1], "failed with method", taxa$bindMethod[i-1])
      message(writeOut)
      writeLines(writeOut, con = out)
    }
    message(paste("Doing", taxa$vect[i], "using spliceRule =", taxa$bindMethod[i], "-- number of taxa in current tree =", newCount))
    lastCount <- newCount
    weldDist <- 0
    trDepth <- -1

	if(taxa$bindMethod[i] == 'genus') {
      tempTaxa <- grep(strsplit(taxa$vect[i],"_")[[1]][1], tr$tip.label, fixed = TRUE, value = TRUE)
      if(length(tempTaxa) == 0) message(paste(taxa$vect[i], "-- returned no matches and was skipped"))
      else {
        taxa$cladeMembers[i] <- paste(tempTaxa, collapse = "|")
        taxa$bindMethod[i] <- "mrca"
      }
    }

    if(taxa$bindMethod[i] == 'mrca') {
      taxaForMRCA <- unlist(strsplit(taxa$cladeMembers[i], "|", fixed = TRUE))
      if(length(taxaForMRCA) == 1) taxa$bindMethod[i] <- "sister"
      else {
        nHeights <- nh(tr)
        weldNodeVector[i] <- findMRCA(tr, taxaForMRCA)
        tr4 <- as(tr, 'phylo4')
        allCladeTaxa <- descendants(tr4, weldNodeVector[i])
        trDepth <- mean(nHeights[as.character(allCladeTaxa), 'height'] - nHeights[as.character(weldNodeVector[i]), 'height'])
      }
    }

    if(taxa$bindMethod[i] %in% c('substitute', 'sister')) { ## currently treats these the same
      if(!taxa$cladeMembers[i] %in% tr$tip.label) message(paste(taxa$vect[i], "skipped -- taxon", taxa$cladeMembers[i], "is missing from your tree."))
      else {
        trEdge <- try(which(tr$edge[, 2] == which(tr$tip.label == taxa$cladeMembers[i])))
        if(class(trEdge) == "try-error") message("***", taxa$vect[i], "failed")
        weldNodeVector[i] <- tr$edge[trEdge, 2]
        weldDist <- trDepth <- tr$edge.length[trEdge] * branchScalar
        }
    }

	if(taxa$bindMethod[i] == 'drop') tr <- drop.tip(tr, taxa$vect[i])

    if(trDepth != -1) {
      trList[[i]] <- list(edge = matrix(c(2L, 1L), 1, 2),
                          tip.label = taxa$vect[i],
                          edge.length = trDepth,
                          Nnode = 1L)
      class(trList[[i]]) <- "phylo"
      ## if(taxa$bindMethod[i] == "sister") save(trList, tr, weldNodeVector, weldDist, file = 'debug.Rdata')
      bindResults = try(bind.tree(tr, trList[[i]], weldNodeVector[i], weldDist))
      if(class(bindResults) == "try-error") message(paste("No matches to", weldNodeVector[i]))
      else {
        message("------------------------------------")
        message(paste("I just made a tree with", length(bindResults$tip.label), "taxa"))
        message(paste("New tip:", taxa$vect[i]))
        message(paste("Node:", weldNodeVector[i]))
        message(paste("Taxon corresponding to that node:",tr$tip.label[weldNodeVector[i]]))
        message(paste("New tree depth and position:", trDepth, weldDist))
        message("")
        tr <- bindResults
      }
    }
  }
  close(out)
  return(tr)
  }

nh <- function(tr) {
    nh <- cbind(as.numeric(tr$edge), as.numeric(nodeHeights(tr)))
    nh <- nh[!duplicated(nh[, 1]), ]
    nh <- nh[order(nh[, 1]), ]
    dimnames(nh) <- list(as.character(nh[, 1]), c("node", "height")) # this matrix gives the distance of each node to the root
    nh
  }

tidyName <- function(x) tolower(gsub(" ", "", gsub(".", "", x, fixed = TRUE), fixed = TRUE))
