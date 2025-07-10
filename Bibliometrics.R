# **************************************************************
# Filename:              Bibliometrics.R
# Created:               Sys.Date()
# Last modified:         Sys.Date() by Hiroki Wakamatsu
# Program Description:   PICES WG49 bibliometric analysis
# Authors:               hwakamatsu@affrc.go.jp
# Data In:               path: "C:\Users\Wakamatsu\Documents\Research\PICES\WG49\Biblio\Bibliometrics.R"
# Data Out:              path: "C:\Users\Wakamatsu\Documents\Research\PICES\WG49\Biblio\Bibliometrics.R"
# **************************************************************

# --------------------------------------------------
# Program Setup
# --------------------------------------------------

# Set working directory (update this path as needed)
setwd("~/Research/PICES/WG49/Biblio")


# --------------------------------------------------
# Load Required Packages
# --------------------------------------------------

# Install and load the 'haven' package to read Stata .dta files
#install.packages("haven")
library(haven)
#install.packages("bibliometrix")
library(bibliometrix)
#install.packages("quanteda")
library(quanteda)
library(ggplot2)
1library(dplyr)
library(stringr)
library(ggplot2)
library(tidygraph)
#install.packages("ggraph")
library(ggraph)

# --------------------------------------------------
# Load Data
# --------------------------------------------------

# Read csv files of raw data (pacific) and screened data (review_data)
pacific <- read.csv("Pacific.csv", stringsAsFactors = FALSE)
review <- read.csv("review_data.csv", stringsAsFactors = FALSE)

# decapitalize, remove spaces and convert to common format
pacific <- pacific %>%
  mutate(
    Title_lower = tolower(str_trim(Title))
  )

review <- review %>%
  mutate(
    Title_lower = tolower(str_trim(Title))
  )

# Matching with titles
# Make titles decap
matched_all <- inner_join(pacific, review, by = "Title_lower")
# Rename var name of keywords from Automatic.Tags.x to Manual.Tags.x because papers contain keywords in either variable
matched_all <- matched_all %>%
  mutate(
    Manual.Tags.x = if_else(
      is.na(Manual.Tags.x) | str_trim(Manual.Tags.x) == "",
      Automatic.Tags.x,
      Manual.Tags.x
    )
  )

# Remove duplicates and remove one with less information
matched_unique <- matched_all %>%
  group_by(Title_lower) %>%
  slice_max(order_by = rowSums(!is.na(across(everything()))), n = 1, with_ties = FALSE) %>%
  ungroup()

# See the result of matching
cat("Matched ", nrow(matched_unique), " / ", nrow(review), "Total\n")

# 3 mismatched papers --> matched data
keys_to_add <- c("FUP3RWJK", "GNLUIS6T", "JGEQLCZ4")

# Detect 3 papers in Pacific.csv
rows_to_add <- pacific %>% filter(Key %in% keys_to_add)

# Add these to matched_unique
matched_unique <- bind_rows(matched_unique, rows_to_add)

# See the result of matching
cat("Matched：", nrow(matched_unique), " / ", nrow(review), "Total\n")

# save to csv
write.csv(matched_unique, "matched_pacific_output.csv", row.names = FALSE)

# Analyze using Bibliometrix
# Convert variable names for bibliometrix
M <- matched_unique %>%
  rename(
    AU = Author.x,            # Authors
    TI = Title.x,             # Title
    SO = `Publication.Title.x`, # Source (Journal or Book Title)
    DT = Item.Type.x,         # Document Type
    PY = Publication.Year.x,  # Publication Year
    DE = Manual.Tags.x,       # Author Keywords (manual tags here)
    AB = Abstract.Note.x,     # Abstract
    DI = DOI.x,               # DOI
    SN = ISSN.x,              # ISSN
    IS = Issue.x,             # Issue
    VL = Volume.x,            # Volume
    BP = Pages.x,             # Beginning page (if page numbers)
    UT = `Archive.Location.x`,# Accession Number / WOS ID
    DB = Library.Catalog.x,    # Database (if applicable)
    ID = Manual.Tags.x
  )

M$DB <- "wos"

write.csv(M, "biblio_ready.csv", row.names = FALSE)
M_csv <- read.csv("biblio_ready.csv", stringsAsFactors = FALSE)
data <- convert2df(as.data.frame(M_csv_raw), dbsource = "isi", format = "csv")



# Display summary of the dataset
results <- biblioAnalysis(data, sep = ";")
summary <- summary(object = results, k = 10)

# Authors' Production over Time
authors <- authorProdOverTime(data)



# Plot collaboration network (authors)
netmatrix <- biblioNetwork(data, analysis = "collaboration", network = "authors", sep = ";")
networkPlot(netmatrix, n = 30, Title = "Collaboration Network", type = "fruchterman", size = T)

keyword_pairs <- keyword_pairs %>%
  mutate(
    from = tolower(str_trim(from)),
    to   = tolower(str_trim(to)),
    
    # california current → Cal. current
    from = ifelse(from == "california current", "Cal. current", from),
    to   = ifelse(to   == "california current", "Cal. current", to),
    
    # ocean acidification → acidification
    from = ifelse(from == "ocean acidification", "acidification", from),
    to   = ifelse(to   == "ocean acidification", "acidification", to),
    
    # marine heatwave variants → MHW
    from = ifelse(from %in% c("marine heatwave", "marine heatwaves", "marine heat wave", "marine heat waves"),
                  "MHW", from),
    to   = ifelse(to %in% c("marine heatwave", "marine heatwaves", "marine heat wave", "marine heat waves"),
                  "MHW", to),
    
    # anthropogenic co2 → co2
    from = ifelse(from == "anthropogenic co2", "co2", from),
    to   = ifelse(to   == "anthropogenic co2", "co2", to),
    
    # climate-change → climate change
    from = ifelse(from == "climate-change", "climate change", from),
    to   = ifelse(to   == "climate-change", "climate change", to),
    
    # climate → climate change
    from = ifelse(from == "climate", "climate change", from),
    to   = ifelse(to   == "climate", "climate change", to),
    
    # sea surface temperature variants → sst
    from = ifelse(from %in% c("sea surface temperature", "sea-surface temperature", "temperature"),
                  "SST", from),
    to   = ifelse(to %in% c("sea surface temperature", "sea-surface temperature", "temperature"),
                  "SST", to),
    
    # el nino variants → enso
    from = ifelse(from %in% c("el nino", "el-nino", "el.nino"), "ENSO", from),
    to   = ifelse(to   %in% c("el nino", "el-nino", "el.nino"), "ENSO", to)
  ) %>%
  filter(
    !from %in% c("pacific", "north pacific", "water", "ocean"),
    !to   %in% c("pacific", "north pacific", "water", "ocean")
  )


# Optional: keyword co-occurrence network
# 1. from側のキーワード出現頻度で上位30を抽出
top_keywords <- keyword_pairs %>%
  count(from, sort = TRUE) %>%
  top_n(40, n) %>%
  pull(from)

# 2. 上位キーワード同士に絞って共起ペア抽出
filtered_pairs <- keyword_pairs %>%
  filter(from %in% top_keywords & to %in% top_keywords)


# グラフ再作成
g_filtered <- graph_from_data_frame(filtered_pairs, directed = FALSE)
layout_matrix <- layout_with_fr(g_filtered)

# ノードを半透明の青系に（alpha = 0.4 で透け感を出す）
node_color <- rgb(0.2, 0.4, 0.8, alpha = 0.4)  # R,G,B = 青系, alpha = 透明度

# エッジも少し透かす（任意）
edge_color <- rgb(0.3, 0.3, 0.3, alpha = 0.5)  # 灰色半透明

# グラフ描画
plot(
  g_filtered,
  layout = layout_matrix,
  vertex.color = node_color,        # ★ 透明ノード
  vertex.label.cex = 1.1,
  vertex.size = degree(g_filtered)^0.5 + 5,
  vertex.label.color = "black",
  edge.color = edge_color,          # 半透明線（任意）
  edge.width = E(g_filtered)$weight,
  main = "Keyword Co-occurrence"
)



#Extract Title + Abstract for quanteda
texts <- paste(data$TI, data$AB, sep = ". ")
docnames <- make.unique(data$TI)

# Create corpus
corp <- corpus(texts, docnames = docnames)

# Tokenize
tokens_all <- tokens(corp, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE)

# Define a topic-specific dictionary
topic_dict <- dictionary(list(
  MHW = c("marine heatwave", "heatwaves", "MHW", "marine heat wave", "marine heat waves", "ENSO", "el-nino", "el nino"),
  Hypoxia = c("hypoxia", "low oxygen", "deoxygenation"),
  OA = c("acidification", "ocean acidification", "ph", "calcification", "california current", "co2"),
  Climate = c("climate change", "warming", "global warming"),
  Social = c("fishery", "fisheries", "aquaculture", "community", "communities", "resilience"),
  Ecosystem = c("ecosystem","zooplankton","fish", "mortality")
))

# Build document-feature matrix and lookup dictionary
dfm_topic <- dfm(tokens_all)
dfm_topics <- dfm_lookup(dfm_topic, dictionary = topic_dict)
topic_counts <- convert(dfm_topics, to = "data.frame")


# Classify documents by dominant topic
topic_cols <- c("MHW", "Hypoxia", "OA", "Climate", "Social", "Ecosystem")

topic_class <- topic_counts %>%
  mutate(doc = doc_id) %>%
  rowwise() %>%
  mutate(
    Top_Topic = topic_cols[which.max(c_across(all_of(topic_cols)))]
  )


# View classification result
head(topic_class[, c("doc", "Top_Topic")])

topic_class %>%
  count(Top_Topic) %>%
  ggplot(aes(x = reorder(Top_Topic, -n), y = n, fill = Top_Topic)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Documents by topic", x = "Topics", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none")

# 1. Join topic classification with publication year
topic_with_year <- topic_class %>%
  left_join(data %>% select(TI, PY), by = c("doc" = "TI"))  # assuming doc names match titles

# 2. Count documents per topic per year, excluding year 2025
topic_trend <- topic_with_year %>%
  filter(PY != 2025) %>%         # Exclude 2025
  count(PY, Top_Topic)

# 3. Plot
ggplot(topic_trend, aes(x = PY, y = n, color = Top_Topic)) +
  geom_line(size = 1) +
  labs(
    title = "Topic Trends by Year",
    x = "Year",
    y = "Number of Documents",
    color = "Topic"
  ) +
  scale_x_continuous(
    breaks = seq(min(topic_trend$PY), max(topic_trend$PY), by = 2) 
  ) +
  theme_minimal()

# Other Plot of Keyword Co-occurrence Network
ggraph(g_tbl, layout = "fr") +  # Fruchterman-Reingold layout
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "gray40") +
  geom_node_point(aes(size = degree(g)), color = "steelblue", alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  theme_void() +
  labs(title = "Keyword Co-occurrence Network") +
  theme(plot.title = element_text(hjust = 0.5))

# Show only labels at higher degree
g_tbl <- g_tbl %>%
  mutate(label = ifelse(degree(g) >= 5, name, ""))

ggraph(g_tbl, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "gray40") +
  geom_node_point(aes(size = degree(g)), color = "steelblue", alpha = 0.7) +
  geom_node_text(aes(label = label), repel = TRUE, size = 4, max.overlaps = Inf) +
  theme_void() +
  labs(title = "Keyword Co-occurrence Network") +
  theme(plot.title = element_text(hjust = 0.5))


