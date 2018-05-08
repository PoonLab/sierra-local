library(jsonlite)

sierra <- fromJSON('./hivdb/sierrapy.json')
local <- fromJSON('./hivdb/local.json')

assertthat::are_equal(nrow(sierra), nrow(local))

complete.matches <- 0
unmatching.headers <- c()
unmatched.drugs <- 0
has.mixtures <- c()

#figure out which ones match
for (sequence in 1:nrow(sierra)) {
  s_scores <- sierra$drugResistance[[sequence]]$drugScores[[1]]
  l_scores <- local$drugResistance[[sequence]]$drugScores[[1]]
  if(identical(s_scores$score, l_scores$score)) {
    complete.matches <- complete.matches + 1
  } else {
    unmatching.headers <- c(unmatching.headers, sierra$inputSequence[1]$header[sequence])
    unmatched.drugs <- unmatched.drugs + sum(s_scores$score != l_scores$score)
  }
  
  if(any(nchar(sierra$alignedGeneSequences[[sequence]]$mutations[[1]]$AAs) > 1)) {
    has.mixtures <- c(has.mixtures, sequence)
  }
}
complete.matches
unmatching.headers
unmatched.drugs

#examine the unmatched
headers <- c()
s_score <- c()
l_score <- c()
drug <- c()
s_partialscore <- c()
l_partialscore <- c()

for (header in unmatching.headers) {
  index <- match(header, sierra$inputSequence[1]$header)
  s_mutations <- sierra$alignedGeneSequences[[index]]$mutations[[1]]
  l_mutations <- local$alignedGeneSequences[[index]]$mutations[[1]]
  s_scores <- sierra$drugResistance[[index]]$drugScores[[1]]
  l_scores <- local$drugResistance[[index]]$drugScores[[1]]
  
  for(i in 1:nrow(s_scores)) {
    headers <- c(headers, header)
    s_score <- c(s_score, s_scores$score[i])
    l_score <- c(l_score, l_scores$score[i])
    drug <- c(drug, s_scores$drug$name[i])
    
    s_p <- ""
    s_l <- ""
    
    for(j in length(s_scores$partialScores[[i]]$mutations)) {
      s_p <- paste(s_p, paste0(s_scores$partialScores[[i]]$mutations[[j]]$text, collapse = "."))
    }
    
    for(j in length(l_scores$partialScores[[i]]$mutations)) {
      s_l <- paste(s_l, paste0(l_scores$partialScores[[i]]$mutations[[j]]$text, collapse = "."))
    }

    s_partialscore <- c(s_partialscore, s_p)
    l_partialscore <- c(l_partialscore, s_l)
  }
}

unmatched <- as.data.frame(cbind(headers, s_score, l_score, drug, s_partialscore, l_partialscore), stringsAsFactors = F)
View(unmatched)

sum(as.numeric(unmatched$l_score) - as.numeric(unmatched$s_score))