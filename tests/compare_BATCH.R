library(jsonlite)

setwd('~/git/sierra-local')
directory <- './hivdb/hivdb-data/'
gene <- 'PR'

sierra_file_list <- list.files(path=directory, pattern=paste0(gene, '-[0-9]+-sierra'))
#local_file_list <- list.files(path=directory, pattern=paste0(gene, '-[0-9]+-local'))
local_file_list <- gsub("sierra", "refactor-v8.5", sierra_file_list)


compare <- function(sierra, local) {
  complete.matches <- 0
  unmatching.headers <- c()
  unmatched.drugs <- 0
  
  #figure out which ones match
  for (sierra_index in 1:nrow(sierra)) {
    # extracts a data frame summarizing HIVdb score by drug
    s_scores <- sierra$drugResistance[[sierra_index]]$drugScores[[1]]
    header <- sierra$inputSequence[1]$header[sierra_index]
    
    # find corresponding sequence label in <local> JSON
    local_index <- match(header, local$inputSequence[1]$header)
    if (is.na(local_index)) { # might be due to header unique count difference
      local_index <- match(paste0(header, '.', sierra_index-1), local$inputSequence[1]$header)
    }
    if (is.na(local_index)) { # last resort
      warning("Failed to match sequence label ", header, ", assuming linear order")
      local_index <- sierra_index
    }
    l_scores <- local$drugResistance[[local_index]]$drugScores[[1]]
    
    if(identical(s_scores$score, l_scores$score)) {
      complete.matches <- complete.matches + 1
    } else {
      unmatching.headers <- c(unmatching.headers, header)
      unmatched.drugs <- unmatched.drugs + sum(s_scores$score %in% l_scores$score)
    }
  }
  
  return (list(
    complete.matches=complete.matches,
    unmatching.headers=unmatching.headers,
    unmatched.drugs=unmatched.drugs
  ))
}



require(parallel)
res <- mclapply(1:length(sierra_file_list), 
  function(i) {
    sfile <- sierra_file_list[i]
    lfile <- local_file_list[i]
    sierra <- fromJSON(paste0(directory, sfile))
    local <- fromJSON(paste0(directory, lfile))
    compare(sierra, local)    
  }, mc.cores=10
)

# post processing
results_unmatched <- sapply(res, function(x) length(x$unmatching.headers))
results_p <- sapply(res, function(x) x$complete.matches / 
                      (x$complete.matches+length(x$unmatching.headers)))
results_unmatched_headers <- sapply(res, function(x) paste(x$unmatching.headers, collapse=' '))

## previous serial code below

#results <- c()
#results_p <- c()
#results_unmatched <- c()
#results_unmatched_headers <- c()
#for (i in 1:length(sierra_file_list)) {
#  sfile <- sierra_file_list[i]
#  lfile <- local_file_list[i]
#  print(sfile)
#  
#  sierra <- fromJSON(paste0(directory, sfile))
#  local <- fromJSON(paste0(directory, lfile))
#  
#  # assertthat::are_equal(nrow(sierra), nrow(local))
#  res.list <- compare(sierra, local)
#  # unpack results
#  complete.matches <- res.list$complete.matches
#  unmatching.headers <- res.list$unmatching.headers
#  unmatched.drugs <- res.list$unmatched.drugs
#  
#  results <- c(results, length(unmatching.headers))
#  results_p <- c(results_p, complete.matches / (complete.matches + length(unmatching.headers)))
#  results_unmatched <- c(results_unmatched, length(unmatching.headers))
#  results_unmatched_headers <- c(results_unmatched_headers, paste(unmatching.headers, collapse=' '))
#}

out <- data.frame(f=sierra_file_list, p=results_p, unmatched=results_unmatched, 
                  headers=results_unmatched_headers)
#View(out)
write.table(out, file=paste(gene, '-analysis-v5.tsv', sep=''), 
            sep='\t', col.names = NA, quote=FALSE)


# RT-13-sierra.json is missing!  I had to type with v8.6.1 using Sierra web service
rt13 <- read.table('hivdb/hivdb-data/RT-13-ResistanceSummary.tsv', sep='\t', header=T)

local <- fromJSON('hivdb/hivdb-data/RT-13-refactor-v8.5.json')
nrow(rt13) == nrow(local)
header <- local$inputSequence[[1]]

# ordering was preserved
all(match(header, rt13$Sequence.Name) == 1:1000)

# extract drug scores
drug.scores <- rt13[, grepl('\\.Score', names(rt13))]
row.names(drug.scores) <- rt13$Sequence.Name
names(drug.scores) <- sapply(names(drug.scores), function(s) strsplit(s, '\\.')[[1]][1])

local.scores <- sapply(local$drugResistance, function(x) {
  x$drugScores[[1]]$score
})
local.scores <- as.data.frame(t(local.scores))
row.names(local.scores) <- header
local.drug.names <- local$drugResistance[[1]]$drugScores[[1]]$drug

# these are in the same order
match(names(drug.scores), local.drug.names$name)
# 3TC is also LMV

res <- sapply(1:nrow(drug.scores), function(i) {
  all(drug.scores[i, 9:19] == local.scores[i,])
})

table(res)  # they all match

