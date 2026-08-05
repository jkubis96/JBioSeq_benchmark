library(iCodon)
library(readxl)

table <- read_xlsx('random_CDS_jbst_extended.xlsx')

score_base <- iCodon::predict_stability('human')

score_base(sequence)

table$native_iCodon <- NA
table$jbst_iCodon <- NA
table$vectorbuilder_iCodon <- NA
table$codontransformer_iCodon <- NA



for (s in 1:nrow(table)) {
  table$native_iCodon[s] <- score_base(table$seq[s])
  table$jbst_iCodon[s] <- score_base(table$jbst_sequence[s])
  table$vectorbuilder_iCodon[s] <- score_base(table$vectorbuilder_sequence[s])
  table$codontransformer_iCodon[s] <- score_base(table$codontransformer_sequence[s])


}


check <- read_xlsx('sources/check_list_mrna.xlsx')

count_matches <- function(pattern, sequences) {
  sum(sapply(sequences, function(x) {
    m <- gregexpr(pattern, x, fixed = TRUE)[[1]]
    if (m[1] == -1) 0 else length(m)
  }))
}

check_df <- data.frame(
  motif = check$Sequence,
  seq = sapply(check$Sequence, count_matches, sequences = table$seq),
  jbst_sequence = sapply(check$Sequence, count_matches, sequences = table$jbst_sequence),
  vectorbuilder_sequence = sapply(check$Sequence, count_matches, sequences = table$vectorbuilder_sequence),
  codontransformer_sequence = sapply(check$Sequence, count_matches, sequences = table$codontransformer_sequence)
)

check_df
