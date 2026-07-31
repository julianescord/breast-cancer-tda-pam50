library(base64enc)
h <- readLines("docs/report.html", warn=FALSE)
h <- paste(h, collapse="\n")
imgs <- c(
  "__IMG_MAPPER__"  = "results/figures/mapper_graph.png",
  "__IMG_KM__"      = "results/figures/km_within_dominant.png",
  "__IMG_VOLCANO__" = "results/figures/c2_volcano.png",
  "__IMG_METABRIC__"= "results/figures/metabric_km_c2.png"
)
for (ph in names(imgs)) {
  uri <- dataURI(file = imgs[[ph]], mime = "image/png")
  h <- gsub(ph, uri, h, fixed = TRUE)
}
writeLines(h, "docs/report_built.html")
cat("OK - tamaño final:", round(file.size("docs/report_built.html")/1024), "KB\n")
