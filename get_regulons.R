getRegulon <- function(reg, tf, d) {
  #d: "neg"/"pos"
  digitpattern <- "[[:digit:]]+\\.[[:digit:]e-]+"
  alltargets <- lapply(1:reg[TF == tf & dir == d, .N], function(i) {
    targets <- strsplit(
      gsub(digitpattern, "",
           gsub(" ", "",
                gsub("'", "",
                     gsub("]", "",
                          gsub("[", "",
                               gsub("[()]", "", reg[TF == tf & dir == d, TargetGenes][i]),
                               fixed = TRUE),
                          fixed = TRUE),
                )
           )
      )
      , ",")[[1]]
    targets[sapply(targets, function(p) {p != ""})]
  })
  Reduce(union, alltargets)
}


library(data.table)
tab <- fread("./SCENIC/test.motifs.csv", sep = ",")


tab[, dir := ifelse(grepl("activating", Context), "pos", "neg")]
regulons <- list()
d <- "neg"
regulons[[d]] <- sapply(unique(tab[dir == d, TF]), function(tf) {  
  getRegulon(tab, tf, d)})  
d <- "pos"
regulons[[d]] <- sapply(unique(tab[dir == d, TF]), function(tf) {  
  getRegulon(tab, tf, d)})

print("final regulons:")
print(sapply(regulons, length))
saveRDS(regulons,"./SCENIC/regulons.rds")