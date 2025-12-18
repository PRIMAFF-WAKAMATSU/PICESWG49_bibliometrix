# **************************************************************
# Filename:              Bibliometrics.R
# Created:               Sys.Date()
# Last modified:         Sys.Date() by Hiroki Wakamatsu
# Program Description:   PICES WG49 bibliometric analysis
# Authors:               hwakamatsu@affrc.go.jp
# **************************************************************

# --------------------------------------------------
# Setup
# --------------------------------------------------
#setwd("C:/Users/Wakamatsu/Dropbox/Biblio")

# --------------------------------------------------
# Load required packages
# --------------------------------------------------
required_packages <- c(
  "dplyr","stringr","quanteda","tidyr","igraph"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages, dependencies = TRUE)
invisible(lapply(required_packages, library, character.only = TRUE))

# --------------------------------------------------
# Load Data
# --------------------------------------------------
review <- read.csv("included.csv", stringsAsFactors = FALSE) %>% as.data.frame()

# --------------------------------------------------
# Text Cleaning
# --------------------------------------------------
clean_text <- function(x) {
  x %>%
    str_remove_all("<[^>]+>") %>%
    str_remove_all("[\\(\\)\"':/\\?]") %>%
    str_replace_all("\\s{2,}", " ") %>%
    str_trim()
}

review <- review %>%
  mutate(
    Title    = clean_text(Title),
    Abstract = clean_text(Abstract)
  )

texts <- paste(review$Title, review$Abstract, sep = ". ")

corp <- corpus(texts, docnames = paste0("doc_", seq_along(texts)))

# --------------------------------------------------
# Tokenization and vocabulary reduction
# --------------------------------------------------
custom_stop <- c(
  "coast","relat","area","differ","indic","signific","condit","observ",
  "studi","result","data","using","used","use","based","analysis","larg","pattern",
  "show","shows","showed","can","may","also","within","among","caus","trend",
  "observed","however","here","paper", "event","high", "ecolog","negat",
  "year","years", "due", "associated", "condition", "conditions","two","time","low",
  "northeast", "period","northern","southern","degre","across","increase","increased","sea",
  "similar","temperature", "reduc","temperatures","change","changes","2014-2016","occur",
  "decad","import","declin","shift","environment","suggest","associ","frequenc","mean","howev","recent",
  "persist","strong", "forag","mix","posit","potenti","analysi","decreas","contribut", "understand",
  "rate","includ","examin","spatial","heat","climate","water","product","chang","level","find",
  "compar","provid","affect","temperatur","record","higher","intens","tropic","influenc","averag","assess",
  "season","along","state","biolog","variabl","region","rang","variat","climat","follow","like","current",
  "found","structur","investig","site","system","earli","near","reveal","day","china","forc",
  "valu","three","identifi","base","scale","develop","role","limit","process","first",
  "western","combin","central","relationship","collect","continu","experienc","histor","interannu",
  "total","unpreced","remain","addit","organ","lead","factor","domin","distribut","major","large-scal",
  "summer","winter","spring","north","south","west","east", #checked if season or location is associated with network graph
  "long-term","subtrop","lower","model"
)


# Abbreviation and variant normalization dictionary
abbr_dict <- list(
  ocean_acidification = c("oa", "o.a."),
  temperature = c("sst"),
  mhw = c("mhw","mhws"),
  el_nino = c("enso", "el nino"),
  anomali = c("anomalies","anomaly")
)

toks <- tokens(
  corp,
  remove_punct   = TRUE,
  remove_numbers = TRUE
) %>%
  tokens_tolower() %>%
  tokens_remove(stopwords("en")) %>%
  # Normalize abbreviations and spelling variants
  tokens_replace(
    pattern = unlist(abbr_dict),
    replacement = rep(names(abbr_dict), lengths(abbr_dict)),
    valuetype = "fixed"
  ) %>%
  # Ensure mhw and mhws are unified
  tokens_replace(
    pattern = c("mhws"),
    replacement = c("mhw"),
    valuetype = "fixed"
  ) %>%
  # Apply English stemming
  tokens_wordstem(language = "english") %>%
  # Force unification of anomal / anomali after stemming
  tokens_replace(
    pattern = c("anomal", "anomali"),
    replacement = c("anomali", "anomali"),
    valuetype = "fixed"
  ) %>%
  # Remove custom stopwords
  tokens_remove(custom_stop) %>%
  # Remove tokens shorter than three characters
  tokens_keep(min_nchar = 3)

# --------------------------------------------------
# Construct and trim document–feature matrix
# --------------------------------------------------
dfm0 <- dfm(toks) %>%
  dfm_trim(min_termfreq = 3) %>%                 # Remove low-frequency terms
  dfm_trim(max_docfreq  = 0.60, docfreq_type="prop")  # Remove overly common terms

# Convert DFM to long format (doc_id, token)
token_ta <- quanteda::convert(dfm0, to="data.frame") %>%
  pivot_longer(cols = -doc_id, names_to = "token", values_to = "freq") %>%
  filter(freq > 0) %>%
  select(doc_id, token)

# --------------------------------------------------
# Dictionary-based semantic standardization
# --------------------------------------------------
token_dict <- list(
  mhw = c("marine heatwave","marine heatwaves","heatwave","heat wave","heat waves",
          "north pacific marine heatwave","northeast pacific marine heatwave",
          "pacific marine heatwave","extreme heatwave","heatwaves"),
  co2 = c("carbon dioxide","anthropogenic co2","co2 emissions"),
  sst = c("sea surface temperature","sea surface temperatures","surface temperature",
          "water temperature","temperature extremes","marine extreme temperature","sst"),
  enso = c("el nino","enso","nino"),
  climate_change = c("climate change","climate","climate event","climate events",
                     "global warming","warm"),
  acidification = c("ocean acidification")
)

standardize_token <- function(x, dict) {
  out <- x
  for (key in names(dict)) out[x %in% dict[[key]]] <- key
  out
}

token_ta <- token_ta %>%
  mutate(token = standardize_token(token, token_dict))

# --------------------------------------------------
# Limit terms used for visualization (Top N by frequency)
# --------------------------------------------------
top_n_words <- 70

term_freq_all <- colSums(dfm0)
top_terms <- names(sort(term_freq_all, decreasing = TRUE))[1:min(top_n_words, length(term_freq_all))]

token_ta_top <- token_ta %>%
  mutate(token = as.character(token)) %>%
  filter(token %in% top_terms)

# --------------------------------------------------
# Construct co-occurrence network (Top N terms)
# --------------------------------------------------
token_pairs <- token_ta_top %>%
  distinct(doc_id, token) %>%
  group_by(doc_id) %>%
  summarise(
    pairs = list({
      toks2 <- sort(unique(token))
      if (length(toks2) >= 2) {
        as.data.frame(t(combn(toks2, 2)), stringsAsFactors = FALSE) %>%
          setNames(c("from","to"))
      } else NULL
    }),
    .groups = "drop"
  ) %>%
  filter(!sapply(pairs, is.null)) %>%
  tidyr::unnest(pairs)

cooccur_df <- token_pairs %>%
  count(from, to, sort = TRUE)

min_edge <- 2
filtered_pairs <- cooccur_df %>%
  filter(n >= min_edge)

g <- graph_from_data_frame(filtered_pairs, directed = FALSE)
E(g)$weight <- filtered_pairs$n

# Remove isolated nodes
g <- induced_subgraph(g, vids = V(g)[degree(g) > 0])

# --------------------------------------------------
# Node size based on term frequency
# --------------------------------------------------
V(g)$tf <- as.numeric(term_freq_all[as.character(V(g)$name)])
V(g)$tf[is.na(V(g)$tf)] <- 0

tf_min <- min(V(g)$tf)
tf_max <- max(V(g)$tf)
vsize <- 6 + 34 * (V(g)$tf - tf_min) / (tf_max - tf_min + 1e-9)

# Edge width based on co-occurrence strength
ew <- E(g)$weight
ew <- 1 + 4 * (ew - min(ew)) / (max(ew) - min(ew) + 1e-9)

# --------------------------------------------------
# Community detection and color assignment
# --------------------------------------------------
clusters <- cluster_louvain(g, weights = E(g)$weight)
V(g)$cluster <- as.character(membership(clusters))

cluster_ids <- sort(unique(V(g)$cluster))
cluster_colors_map <- setNames(rainbow(length(cluster_ids)), cluster_ids)
V(g)$color <- adjustcolor(cluster_colors_map[V(g)$cluster], alpha.f = 0.55)

# --------------------------------------------------
# Plot the co-occurrence network
# --------------------------------------------------
set.seed(123)
layout_fr <- layout_with_fr(g, weights = E(g)$weight, niter = 2000)

par(mar = c(0, 0, 3, 0))
plot(
  g,
  layout = layout_fr,
  vertex.label = V(g)$name,
  vertex.label.cex = 0.8,
  vertex.label.color = "black",
  vertex.color = V(g)$color,
  vertex.frame.color = "black",
  vertex.size = vsize,
  edge.width = ew,
  edge.color = "gray70",
  main = paste0("Token Co-occurrence Network (Top ", top_n_words, " terms)")
)

legend(
  "bottomleft",
  legend = paste0("Cluster ", cluster_ids),
  col = cluster_colors_map[cluster_ids],
  pch = 19,
  pt.cex = 1.2,
  bty = "n"
)



