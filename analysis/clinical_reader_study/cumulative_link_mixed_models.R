library(ordinal)
library(emmeans)

infile  <- "~/Downloads/collated_reader_study_results.csv"
out_dir <- "~/Downloads"
exclude_ids <- c(29, 73, 105, 122, 141, 143, 152, 166, 184)

df <- read.csv(infile, stringsAsFactors = FALSE)
df <- df[!df$global_set_id %in% exclude_ids, ]

df$global_set_id <- factor(df$global_set_id)
df$user_id       <- factor(df$user_id)
df$source        <- factor(df$source, levels = c(
  "original_history",
  "additional_history",
  "llm_indication_claude",
  "llm_indication_qwen"
))
df$source <- relevel(df$source, ref = "original_history")

source_lookup <- c(
  "original_history"      = "Referring Clinician",
  "additional_history"    = "Radiologist",
  "llm_indication_claude" = "Claude 3.5 Sonnet",
  "llm_indication_qwen"   = "Qwen 2.5-7B Instruct"
)

outcomes <- c("comprehensiveness", "factuality", "conciseness")

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("< .001")
  formatC(p, digits = 3, format = "f")
}

for (out in outcomes) {
  df[[paste0(out, "_ord")]] <- ordered(df[[out]])
}

models <- lapply(setNames(outcomes, outcomes), function(out) {
  form <- as.formula(paste0(out, "_ord ~ source + (1|global_set_id) + (1|user_id)"))
  clmm(form, data = df, link = "logit", Hess = TRUE)
})

table4_list <- lapply(names(models), function(name) {
  emm <- emmeans(models[[name]], ~ source)
  res <- as.data.frame(contrast(emm, method = "trt.vs.ctrl", ref = 1, infer = c(TRUE, TRUE)))
  
  nonref_raw <- gsub(" - .*", "", res$contrast)
  
  data.frame(
    Criterion = tools::toTitleCase(name),
    Source = source_lookup[nonref_raw],
    `Odds Ratio` = round(exp(res$estimate), 3),
    `Confidence Interval` = paste0("(", round(exp(res$asymp.LCL), 3), ", ", 
                                   round(exp(res$asymp.UCL), 3), ")"),
    `P value` = sapply(res$p.value, fmt_p),
    check.names = FALSE
  )
})

table4 <- do.call(rbind, table4_list)
write.csv(table4, file.path(out_dir, "table4_combined_by_criterion.csv"), row.names = FALSE)

table_s6_list <- lapply(names(models), function(name) {
  emm <- emmeans(models[[name]], ~ source)
  res <- as.data.frame(pairs(emm, adjust = "tukey", infer = c(TRUE, TRUE)))
  
  source_a_raw <- gsub(" - .*", "", res$contrast)
  source_b_raw <- gsub(".* - ", "", res$contrast)
  
  data.frame(
    Criterion = tools::toTitleCase(name),
    `Source A VS Source B` = paste0(source_lookup[source_a_raw], " vs ", source_lookup[source_b_raw]),
    `Logit Estimate` = round(res$estimate, 4),
    `Standard Error` = round(res$SE, 4),
    `P Value` = sapply(res$p.value, fmt_p),
    check.names = FALSE
  )
})

table_s6 <- do.call(rbind, table_s6_list)

src_names <- names(source_lookup)
order_indices <- which(lower.tri(matrix(0, length(src_names), length(src_names))), arr.ind = TRUE)
desired_levels <- paste0(source_lookup[src_names[order_indices[,2]]], " vs ", source_lookup[src_names[order_indices[,1]]])

table_s6$`Source A VS Source B` <- factor(table_s6$`Source A VS Source B`, levels = unique(desired_levels))
table_s6 <- table_s6[order(factor(table_s6$Criterion, levels = tools::toTitleCase(outcomes)), 
                           table_s6$`Source A VS Source B`), ]

write.csv(table_s6, file.path(out_dir, "s6_combined_table_S6_format.csv"), row.names = FALSE)

