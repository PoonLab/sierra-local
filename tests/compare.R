library(jsonlite)

setwd('~/git/sierra-local')

g <- 'RT'
num <- 98  # which subset to analyze
sierra <- fromJSON(paste0('./hivdb/hivdb-data/', g, '-', num, '-sierra.json'))
local <- fromJSON(paste0('./hivdb/hivdb-data/', g, '-', num, '-local.json'))

# sierra <- fromJSON(paste0('./hivdb/hivdb-data/DQ-local.json'))
# local <- fromJSON(paste0('./hivdb/hivdb-data/DQ-sierra.json'))


assertthat::are_equal(nrow(sierra), nrow(local))

complete.matches <- 0
# unmatching.headers <- c("ABC", "AZT", "FTC", "3TC", "TDF", "D4T", "DDI",
#                         "BIC", "DTG", "EVG", "RAL",
#                         "EFV", "ETR", "NVP", "RPV",
#                         "ATV/r", "DRV/r", "LPV/r", "FPV/r",
#                         "IDV/r", "NFV", "SQV/r", "TPV/r")
unmatching.headers <- c()
unmatched.drugs <- 0

#figure out which ones match
for (sierra_index in 1:nrow(sierra)) {
  s_scores <- sierra$drugResistance[[sierra_index]]$drugScores[[1]]
  header <- sierra$inputSequence[1]$header[sierra_index]
  local_index <- sierra_index
  l_scores <- local$drugResistance[[local_index]]$drugScores[[1]]
  
  if(identical(s_scores$score, l_scores$score)) {
    complete.matches <- complete.matches + 1
  } else {
    unmatching.headers <- c(unmatching.headers, header)
    unmatched.drugs <- unmatched.drugs + sum(s_scores$score %in% l_scores$score)
  }
}

# sierra.table <- data.frame(matrix(nrow=nrow(sierra), ncol=length(unmatching.headers)))
# colnames(sierra.table) <- unmatching.headers
# rownames(sierra.table) <- sierra$inputSequence[1]$header

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
  sierra_index <- match(header, sierra$inputSequence[1]$header)
  # local_index <- match(header, local$inputSequence[1]$header)
  # if (is.na(local_index)) { # might be due to header unique count difference
  #   local_index <- match(paste0(header, '.', sierra_index), local$inputSequence[1]$header)
  # }
  local_index <- sierra_index
  s_mutations <- sierra$alignedGeneSequences[[sierra_index]]$mutations[[1]]
  l_mutations <- local$alignedGeneSequences[[local_index]]$mutations[[1]]
  s_scores <- sierra$drugResistance[[sierra_index]]$drugScores[[1]]
  l_scores <- local$drugResistance[[local_index]]$drugScores[[1]]

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
unmatched$s_level <- cut(as.numeric(unmatched$s_score), breaks=c(-Inf, 10, 15, 30, 60, Inf), labels=c(1, 2, 3, 4, 5))
unmatched$l_level <- cut(as.numeric(unmatched$l_score), breaks=c(-Inf, 10, 15, 30, 60, Inf), labels=c(1, 2, 3, 4, 5))

View(unmatched)

25
50
AZT
K219R
M41L.L210W
3
4
3
KC221011.35963.B.939
35
60