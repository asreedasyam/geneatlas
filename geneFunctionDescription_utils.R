library(data.table)
library(dplyr)
library(magrittr)
library(sentimentr)
library(lexicon)
library(tidytext)
library(tibble)
library(qdapRegex)


#' Title generate lexicon for gene function descriptions
#'
#' @param annotation.files 
#' @param description.column 
#' @param file.suffix 
#' @param out.dir 
#'
#' @return
#' @export
#'
#' @examples
generate.lexicon_for_geneDescriptions <- function(annotation.files, 
                                                  description.column = "Description",
                                                  file.suffix, 
                                                  out.dir = "."){
  
  ## get function descriptions or deflines 
  deflines <- unique(unlist(lapply(annotation.files, function(i) {
    s <- data.table::fread(i)[[description.column]]
    return(s)
    
  })))
  
  ## -- Remove strings bounded between a left and right marker.
  stripBetween <- function(x, start = "<", end = ">", fixed = TRUE){
    qdapRegex::rm_between(x, start, end, fixed = fixed)
  }
  
  
  deflines <- stripBetween(deflines , 
                           start = c("\\[ co-ortho", "\\[co-ortho", "\\[ ortho", "\\[ortho", "\\[ORG"),
                           end = c("$|\\]", "$|\\]", "$|\\]", "$|\\]", "$|\\]"),
                           fixed = F)
  
  
  text_df <- tibble(line = c(1:length(deflines)), text = deflines)
  
  defl_words <- text_df %>% unnest_tokens(word, text)
  
  ## remove stop words
  data(stop_words)
  defl_words %<>%
    anti_join(stop_words)
  defl_words <- unique(defl_words$word)
  
  ## remove punctuations
  defl_words_filtered <- sapply(defl_words, function(x) get_uniqueWords(x), USE.NAMES = F)
  defl_words_filtered <- unlist(strsplit(defl_words_filtered, " "), use.names = F)
  ## remove elements which are just numbers after removing punctuations
  defl_words_filtered <- sapply(defl_words_filtered, function(x) {if(grepl("[a-zA-Z]+",x)) x else gsub("[0-9]","",x)}, USE.NAMES = FALSE)
  defl_words_filtered <- defl_words_filtered[defl_words_filtered!=""]
  ## remove single alphabet
  defl_words_filtered <- defl_words_filtered[!grepl("\\b[aA-zZ]\\b",defl_words_filtered)]
  ## remove words like DUF630, UCP030210, PD694200, UPF0497
  
  rm_unfunctions <- c("^(pd)([0-9]+)", "^(sf)([0-9]+$)", "^(map)([0-9]+$)", "^(duf)([0-9]+)", "(upf)([0-9]+)", "(ucp)([0-9]+)",
                      "^(pf)([0-9]+)", "pfam", "^(osfbx)([0-9]+)", "^(caaa)([0-9]+)", "^[0-9](caaa)([0-9]+)", "^[0-9]([a-zA-Z])([0-9]{6,7})",
                      "^[0-9]([a-zA-Z]{1,5})([0-9]{4,8})",
                      "^[0-9]([a-zA-Z]+)([0-9]+)([a-zA-Z])([0-9]+)", "^([0-9]{6,7})([a-zA-Z]$)",
                      "^(f2g)([0-9]+)",
                      "^([0-9]+)(of)([0-9]+)", "^([a-zA-Z]+)([0-9]{4,10}$)", "^([a-zA-Z]+)([0-9]{1,3})([a-zA-Z])([0-9]{3,8}$)")
  
  defl_words_filtered <- defl_words_filtered[!grepl(paste(rm_unfunctions, collapse = "|"),defl_words_filtered)]
  
  
  
  ## -- Go over scores of a few organisms to add/delete keywords/sentences in add_defline_priority function
  # Example:
  # similar to expressed protein in Arabidopsis thaliana; [ co-ortholog (2of2) of At5g56520, ]
  # Protein of Unknown Function (DUF239); similar to expressed protein in Arabidopsis thaliana; [ co-ortholog (2of3) of At3g13510, At1g55360, At5g56530, ]
  # similar to expressed protein in Arabidopsis thaliana; [ ortholog of At5g12900,]
  # Uncharacterised protein family UPF0090; similar to expressed protein in Arabidopsis thaliana; [ ortholog of At1g77122,]
  ## add priority info for each word/sentence
  ## NOTE: we can manually include words/sentences and their priorities here
  add_defline_priority <- function(s, 
                                   priority_low= NULL, priority_vlow = NULL,
                                   priority_vvlow = NULL, priority_lowest = NULL){
    
    if(is.null(priority_low))
      priority_low <- c("carrier","conserved","peptide")
    #priority_low <- c("binding","carrier","conserved","peptide")
    
    if(is.null(priority_vlow))
      priority_vlow <- c("domain","binding","function","proteins","protein","expressed","family","similar",
                         "plant","containing","transcription","superfamily","hypothetical","catalytics")
    
    if(is.null(priority_vvlow))
      priority_vvlow <- c("arabidopsis","thaliana","populus","eukaryotic","duplicated", "putativ", "putative","putatively", "ortholog", "co-ortholog", "factor")
    
    if(is.null(priority_lowest))
      priority_lowest <- c("uncharacterized","uncharacterised","unknown","decoy", "duf", "ucp", "upf")
    
    df <- data.table::data.table(words = s)
    df <- df %>% mutate(priority = dplyr::case_when(words %in% priority_low ~ 0.2,
                                                    words %in% priority_vlow ~ 0.05,
                                                    words %in% priority_vvlow ~ 0.01,
                                                    words %in% priority_lowest ~ -2,
                                                    !words %in% priority_low | !words %in% priority_vlow | !words %in% priority_vvlow | !words %in% priority_lowest ~ 1))
    
    ### NOTE: Add additional words or sentences and their priorities here ###
    additional_words <- c("open reading frame","expressed protein", "always early", "blast2go", "org")
    additional_priority <- c(0, 0, 1, 0, 0)
    df <- dplyr::bind_rows(df, data.table(words = additional_words, priority = additional_priority))
    
    df <- unique(df)
    return(df)
  }
  defl_words_filtered_df <- add_defline_priority(defl_words_filtered)
  
  ## modify valence shifters
  remove_valence_shifters <- c("no", "not", "high", "highly", "little", "heavy", "extra", "never")
  
  ## Check if any valence shifters are present in the built-in list
  check_valenceShifters_in_deflines = FALSE
  if(check_valenceShifters_in_deflines){
    def_wrds <- sapply(deflines, function(x) strsplit(x, "\\s+"), USE.NAMES = F)
    for(i in check_valence_shifters){
      message(rep("* = ",20),"\n\t", i, "\n", rep("* = ",20))
      chks <- grep(paste0("\\b",i,"\\b"), def_wrds)
      message("Total found : ",length(chks), " / ", length(def_wrds))
      if(length(chks) > 3)
        chks <- sample(chks, size = 3, replace = F)
      print(deflines[chks])
    }
  }
  hash_valence_shifters_defl <- lexicon::hash_valence_shifters
  hash_valence_shifters_defl <- hash_valence_shifters_defl  %>% filter(!x %in% remove_valence_shifters) %>% data.table
  write.csv(hash_valence_shifters_defl, file = file.path(out.dir, 
                                                         file = sprintf("hash_valence_shifters_defl_PolarityDT_%s.csv", file.suffix)))
  
  write.csv(defl_words_filtered_df, file = file.path(out.dir,
                                                     file = sprintf("hash_sentiment_functionDescription_words_PolarityDT_%s.csv",file.suffix)))
}




#' Title Categorize gene function descriptions based on sentiment score
#'
#' @param annotation.file 
#' @param polarity.file 
#' @param file.suffix 
#' @param description.column 
#' @param sentiment_threshold 
#' @param out.dir 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
categorize.functionDescriptions <- function(annotation.file, 
                                            polarity.file,
                                            file.suffix,
                                            description.column = "Description",
                                            sentiment_threshold = 0.1999,
                                            out.dir = ".",
                                            verbose = TRUE){
  
  
  fd <- fread(annotation.file)
  if(any(!c("GeneID", description.column) %in% colnames(fd)))
    stop(sprintf("annotation.file must contain GeneID and %s columns", description.column))
  if(!file.exists(polarity.file))
    stop("** provide polarity file **")
  polarity_data <- data.table::fread(polarity.file, drop = 1)
  
  remove_valence_shifters <- c("no", "not", "high", "highly", "little", "heavy", "extra", "never")
  valence_shifters_data <- lexicon::hash_valence_shifters
  valence_shifters_data <- valence_shifters_data  %>% filter(!x %in% remove_valence_shifters) %>% data.table
  
  polarity_data <- polarity_data %>% data.frame
  if(!is_key(polarity_data))
    polarity_data <- as_key(polarity_data, comparison = valence_shifters_data)
  
  ## REMOVE DUPLICATE WORDS IN FUNCTION DESCRIPTION
  ss <- lapply(fd[[description.column]], function(x) get_uniqueWords(x))
  
  ss <- sentimentr:::make_class(ss, "get_sentences", "get_sentences_character")
  
  make_sentence_df2 <- function(sents, retention_regex = NULL){
    indx <- wc <- NULL
    ids <- sentimentr:::add_row_id(sents)
    
    ## ORIGINAL    
    ##    text.var <- gsub("[^a-z',;: ]|\\\\d:\\\\d|\\\\d ", " ", stringi::stri_trans_tolower(gsub("(\\\\s*)([;:,]+)",
    ##        " \\\\2", unlist(sents))))
    text.var <- gsub("[^a-z',;: 0-9]", " ", stringi::stri_trans_tolower(gsub("(\\\\s*)([;:,]+)"," \\\\2", unlist(sents))))
    dat <- data.frame(id = ids, sentences = text.var, wc = sentimentr:::count_words(text.var),
                      stringsAsFactors = FALSE)
    data.table::setDT(dat)
    dat[, `:=`(indx, wc < 1), by = c("id", "sentences", "wc")][(indx),
                                                               `:=`(c("sentences", "wc"), NA)][, `:=`(indx, NULL)][]
  }
  ## Change the original function in sentimentr to retain say ycf9 as ycf9 and b12d as b12d instead of (ycf and b d)
  assignInNamespace("make_sentence_df2", make_sentence_df2, ns = "sentimentr")
  
  senti_res <- sentiment_by(ss, polarity_dt = polarity_data, valence_shifters_dt = valence_shifters_data)
  
  senti_res$GeneID <- fd$GeneID
  
  senti_res[[description.column]] <- fd[[description.column]]
  
  if(verbose)
    message(sprintf("Number of Deflines with average sentiment below %.3f = %d / %d", 
                    sentiment_threshold, 
                    (senti_res %>% filter(ave_sentiment < sentiment_threshold) %>% nrow),
                    nrow(fd)))
  senti_res %<>% mutate(Category = case_when(ave_sentiment > sentiment_threshold ~ "Good",
                                             ave_sentiment == sentiment_threshold ~ "Good",
                                             ave_sentiment < sentiment_threshold ~ "Poor") )
  
  senti_res %<>% dplyr::select(c('GeneID', all_of(description.column), 'Category', 'ave_sentiment', 'word_count'))
  
  data.table::setnames(senti_res, old = c("ave_sentiment", "word_count"), new = c("Sentiment_score", "Word_count"))
  if(verbose)
    message(sprintf("Median sentiment score = %.2f", median(senti_res$Sentiment_score)))
  
  write.csv(senti_res, file = file.path(out.dir,sprintf("FunctionDescription_categories_%s.csv",file.suffix)))
  return(senti_res)
}
