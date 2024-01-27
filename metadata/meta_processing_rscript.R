#import library
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#import metadata extracted from json
geo <- read_csv("geo.csv", col_names = FALSE)
year <- read_csv("year.csv", col_names = FALSE)
hosts <- read_csv("hosts.csv", col_names = FALSE)
completeness <- read.csv("completeness", col_names = FALSE)
spc <- read_csv("spc.csv", col_names = FALSE)

#import checkm data
combined <- read.csv("E:/zyf/ucl/mb/ncbicli/checkm_report.csv", header=FALSE)


#processing to combine columns and remove uninformative value
hosts <- hosts %>%
  unite("combined", c("X2", "X3", "X4", "X5", "X6",), sep = ",", na.rm = TRUE)
hosts$host <- hosts$combined %>%
  str_replace_all("N/A|missing|not applicable|Unknown|,", "") %>%
  str_replace_all("^,|,$|,,+", "") %>%
  str_trim()
hosts <- subset(hosts, select = -combined)

year <- year %>%
  unite("combined", c("X2", "X3", "X4", "X5", "X6",), sep = ",", na.rm = TRUE)
year$date <- year$combined %>%
  str_replace_all("N/A|not collected|,|missing|not applicable", "") %>%
  str_replace_all("^,|,$|,,+", "") %>%
  str_trim()
year <- subset(year, select = -combined)

geo <- geo %>%
  unite("combined", c("X2", "X3", "X4", "X5", "X6",), sep = ",", na.rm = TRUE)
geo$geo <- geo$combined %>%
  str_replace_all("N/A|,|missing|not applicable|not collected", "") %>%
  str_replace_all("^,|,$|,,+", "") %>%
  str_trim()
geo <- subset(geo, select = -combined)


#rename column and parse metadata using genome accession
colnames(hosts) <- c("genome", "host")
colnames(year) <- c("genome", "date")
colnames(completeness) <- c("genome", "completeness")
colnames(geo) <- c("genome", "geo")
colnames(spc) <- c("genome", "spc")
df <- merge(hosts, year, by = "genome", all.x = TRUE)
df <- merge(df, completeness, by = "genome", all.x = TRUE)
df <- merge(df, geo, by = "genome", all.x = TRUE)
df <- merge(df, spc, by = "genome", all.x = TRUE)
df <- merge(df, contamination, by = "genome", all.x = TRUE)

#remove genus name mycobacterium to process in species level
df$spc <- df$spc %>%
  str_replace_all("Mycobacterium", "") %>%
  str_replace_all("^,|,$|,,+", "") %>%
  str_trim()

#processing date and geological data to remove missing data and retain only year and country value
df$date <- df$date %>%
  str_replace_all("Unknown", "") %>%
  str_replace_all("^,|,$|,,+", "") %>%
  str_trim()
df$geo <- gsub("Missing|not determined|^$", "", filtered$geo, ignore.case = TRUE)
df$geo <- gsub(":.*$", "", filtered$geo)


#parsing own checkm report because some genome from metadata miss contamination/completeness data
combined <- combined[, -c(2,3,4,5,6,7,8,9,10,11,14)]
colnames(combined) <- c("genome", "comp", "cont")
combined$genome <- sub("([^_]*_[^_]*)_.*", "\\1", combined$genome)
combined$genome <- sub("_.*","", combined$genome)
df <- merge(df, combined, by = "genome", all.x = TRUE)
df$contamination[is.na(df$contamination) & !is.na(df$cont)] <- df$cont[is.na(df$contamination) & !is.na(df$cont)]
df$completeness[is.na(df$completeness) & !is.na(df$comp)] <- df$comp[is.na(df$completeness) & !is.na(df$comp)]

#extract data with >95% completeness and <5% contamination
filtered <- subset(df, completeness > 95)
filtered <- subset(filtered, contamination < 5)



#human MTB subsampling
#avoid overrepresentation of MTB species with human host by subsampling
df2 <- filtered %>% filter(spc == "tuberculosis")
df2$date <- gsub(".*(\\d{4}).*", "\\1", df2$date)
df2$completeness <- as.numeric(df2$completeness)
df2$contamination <- as.numeric(df2$contamination)
df3 <- subset(df2, is.na(host) | host != "Homo sapiens")

non_mtb <- filtered %>%
  filter(spc != 'tuberculosis')
non_hmtb <- rbind(df3, non_mtb) #Mycobacterium species which are not MTB with human hosts
df2 <- subset(df2, !is.na(host) & host == "Homo sapiens")

subhmtb <- df2 %>%  #MTB with human host is subsampled by collection date and country
  distinct(geo, date, .keep_all = TRUE)
sub <- rbind(non_hmtb, subhmtb)

#write the list of genome accession for genome annotation (prokka)
write.table(sub$genome, "acc_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)





















