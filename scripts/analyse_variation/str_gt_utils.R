# library(tidyverse)

load_str_panel <- function(filepath) {
  df_panel <- read.delim(filepath, header=FALSE, na.strings=".") %>% 
    cbind(., str_split_fixed(.$V7, "[:]", n=Inf)) %>% 
    dplyr::select(-c(V5, V6, V7, `3`)) %>% 
    mutate(
      `1` = ifelse(`1` == ".", NA, `1`),
      `4` = as.integer(`4`)
    )
  
  colnames(df_panel) <- c("chr", "start", "end", "period", "repeat_id", "msa", "max_p_stretch")
  
  df_panel <- df_panel %>% 
    mutate(
      tmp_id = paste0(chr, "_", start),
      ref = (end - start + 1) / period
    )
  
  return(df_panel)
}

add_panel_info <- function(df, panel, new_colname) {
  df <- df %>% 
    left_join(
      panel %>% select(tmp_id) %>% distinct() %>% mutate(tmp = TRUE),
      by="tmp_id"
    ) %>%
    mutate(tmp = if_else(is.na(tmp), FALSE, TRUE))
  
  colnames(df)[ncol(df)] <- new_colname
  
  return(df)
}

calc_patient_len_dif <- function(a0, b0, a1, b1) {
  option_a = abs(a0 - a1) + abs(b0 - b1)
  option_b = abs(a0 - b1) + abs(b0 - a1)
  return(min(option_a, option_b))
}
