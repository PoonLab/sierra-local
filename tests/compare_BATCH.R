library(jsonlite)

directory <- './hivdb/hivdb-data/'
gene <- 'RT'

sierra_file_list <- list.files(path=directory, pattern=paste0(gene, '-[0-9]+-sierra'))
local_file_list <- list.files(path=directory, pattern=paste0(gene, '-[0-9]+-local'))

results <- c()
results_p <- c()
results_unmatched <- c()
results_unmatched_headers <- c()

for (file in sierra_file_list) {
  sierra <- fromJSON(paste0(directory,file))
  local <- fromJSON(paste0(directory,gsub('sierra', 'local', file)))
  
  # assertthat::are_equal(nrow(sierra), nrow(local))
  
  complete.matches <- 0
  unmatching.headers <- c()
  unmatched.drugs <- 0
  
  #figure out which ones match
  for (sierra_index in 1:nrow(sierra)) {
    s_scores <- sierra$drugResistance[[sierra_index]]$drugScores[[1]]
    header <- sierra$inputSequence[1]$header[sierra_index]
    
    local_index <- match(header, local$inputSequence[1]$header)
    if (is.na(local_index)) { # might be due to header unique count difference
      local_index <- match(paste0(header, '.', sierra_index-1), local$inputSequence[1]$header)
    }
    if (is.na(local_index)) { # last resort
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
  
  complete.matches
  unmatching.headers
  unmatched.drugs
  
  results <- c(results, length(unmatching.headers))
  results_p <- c(results_p, complete.matches / (complete.matches + length(unmatching.headers)))
  results_unmatched <- c(results_unmatched, length(unmatching.headers))
  results_unmatched_headers <- c(results_unmatched_headers, paste(unmatching.headers, collapse=' '))
}

out <- data.frame(f=sierra_file_list, p=results_p, unmatched=results_unmatched, headers=results_unmatched_headers)
View(out)
write.table(out, file=paste(gene, 'analysisv4.tsv'), sep='\t', col.names = NA, quote=FALSE)
