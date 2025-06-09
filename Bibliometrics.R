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

# --------------------------------------------------
# Load Data
# --------------------------------------------------

# Load Stata data file (update the file path as needed) 
file <- "savedrecs.txt"
data <- convert2df(file, dbsource = "wos", format = "plaintext")
# Display summary of the dataset
summary(data)

# Use the 'AB' (abstract) column
texts <- data$AB  
names(texts) <- data$TI  # Optionally name the texts by their titles

# Create a dictionary for three climate-related topics
topic_dict <- dictionary(list(
  MHW = c("marine heatwave", "heatwaves", "MHW"),
  Hypoxia = c("hypoxia", "low oxygen", "deoxygenation"),
  OA = c("acidification", "ocean acidification", "pH")
))
# Tokenize the texts (remove punctuation and numbers)
docnames <- make.unique(as.character(data$TI))
corp <- corpus(texts, docnames = docnames)
# Create a document-feature matrix (DFM)
dfm_topic <- dfm(tokens_all)

# Apply the dictionary to count topic-specific terms
dfm_lookup <- dfm_lookup(dfm_topic, dictionary = topic_dict)
dfm_df <- convert(dfm_lookup, to = "data.frame")

# View all 50 rows
print(dfm_df, row.names = FALSE) 

# Convert to data frame for plotting or aggregation
topic_counts <- convert(dfm_lookup, to = "data.frame")

# Total frequency by topic
colSums(dfm_lookup)

# Plot the topic frequency
topic_sums <- data.frame(
  topic = names(colSums(dfm_lookup)),
  count = as.numeric(colSums(dfm_lookup))
)

ggplot(topic_sums, aes(x = topic, y = count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Topic Frequency in Abstracts", x = "Topic", y = "Keyword Matches")




