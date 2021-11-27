library(rvest)
library(dplyr)
library(purrr)
library(Matrix)
bands_source <- read_html("http://www.redhotjazz.com/band_body.html")
bands <- bands_source %>%
  html_node(xpath = "/html/body/center/table") %>%
  html_table(header = TRUE)
bands_url <- bands_source %>%
  html_node(xpath = "/html/body/center/table") %>%
  html_nodes("a") %>%
  html_attr("href")
bands_name <- bands_source %>%
  html_node(xpath = "/html/body/center/table") %>%
  html_nodes("a") %>%
  html_text()
bands_full <- tibble::tibble(name = bands_name, url = bands_url) %>% left_join(bands, by = c("name" = "Name of the Band"))
names(bands_full) <- c("Band", "url", "L", "Years", "Recorded", "Live")
bands_full <- bands_full %>%
  mutate(Band = gsub("[\n,\r]", "", Band)) %>%
  mutate(Band = stringi::stri_unescape_unicode(Band)) %>%
  group_by(Band) %>%
  summarise(url = paste(url, collapse = ","), Years = paste(Years, collapse = ","), Recorded = paste(Recorded, collapse = ","), Live = paste(Live, collapse = ","))
base_url <- "http://www.redhotjazz.com/"
bands_details <- list()
for (u in bands_url) {
  print(u)
  try({
    source <- read_html(paste(base_url, u, sep = "/"))
    musicians <- source %>%
      html_node(xpath = "/html/body/center[6]/table") %>%
      html_table(header = TRUE, fill = TRUE)
    musicians$band <- u
    bands_details[[u]] <- musicians
  })
}

bands_details <- bands_details %>% keep(function(b) {
  all(colnames(b) == c("Artist", "Instrument", "band"))
})
musicians_bands <- do.call(rbind, bands_details) %>%
  left_join(bands_full, by = c("band" = "url")) %>%
  rename(url = band)
G <- musicians_bands %>%
  select(Artist, Band) %>%
  group_by(Artist, Band) %>%
  summarise(n = 1) %>%
  ungroup() %>%
  filter(!is.na(Artist), !grepl("\\?", Artist)) %>%
  mutate(Band = gsub("[\n,\r]", "", Band)) %>%
  mutate(Band = stringi::stri_unescape_unicode(Band)) %>%
  mutate(Artist = gsub('[\n,\r,\\",\\,\\.,_,-]', "", Artist)) %>%
  mutate(Artist = stringi::stri_unescape_unicode(Artist)) %>%
  tidyr::spread(Band, n, fill = 0)

Xns <- as.matrix(G[, -1])
rownames(Xns) <- G$Artist
ij <- which(X > 0, arr.ind = TRUE)
X <- sparseMatrix(ij[, 1], ij[, 2], x = Xns[ij])
colnames(X) <- colnames(Xns)
rownames(X) <- rownames(Xns)

col_meta <- tibble(Band = colnames(X)) %>% left_join(bands_full)

Jazz_full <- list(X = X, col_meta = col_meta)
devtools::use_data(Jazz_full)

Xb <- t(X) %*% X
diag(Xb) <- 0
Xb[Xb > 1] <- 1

sol <- greed(Xb)


Xm <- X %*% t(X)
diag(Xm) <- 0
Xm[Xm > 1] <- 1

sol <- greed(Xm)
