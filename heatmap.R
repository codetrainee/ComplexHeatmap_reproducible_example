## Prepare R environment
rvst <- R.version.string
file_prefix <- strsplit(rvst, split = " ")[[1]][[3]]
path <- sprintf("./%s", file_prefix)
if (!dir.exists(path)) dir.create(path)
repos = "https://cran.rstudio.com"
.libPaths(path)
install.packages("BiocManager", repos = repos)
library("BiocManager")
BiocManager::version()
BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
RNGversion("3.5") # R3.6 set.seed: major change

## First heatmap result ====================================
all <- readRDS("all.rds")

## Same with heatmap for sig proteins
all_ha_df <- data.frame(Comparison = colnames(all))
all_ha_df$Comparison <- gsub("(.+)[0-9]{1}$", "\\1", all_ha_df$Comparison)

## construct the customised annotation
all_ha <-
  columnAnnotation(
    df = all_ha_df,
    col = list(Comparison = setNames(
      rainbow(5),
      sort(unique(all_ha_df$Comparison))
    )),
    annotation_legend_param = list(
      Comparison = list(
        nrow = 3,
        title = "Comparison",
        title_gp = gpar(fontsize = 8),
        labels_gp = gpar(fontsize = 6)
      )
    )
  )

## set the seed for reproducibility
set.seed(123)

all_heatmap <-
  Heatmap(
    all,
    top_annotation = all_ha,
    heatmap_legend_param = list(
      title = "Log2 Fold Change",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 6),
      legend_direction = "horizontal",
      legend_width = unit(1.5, "in")
    ),
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE
  )

png(file = sprintf("%s_all.png", file_prefix), res = 300, width = 6, height = 6, units = "in")

draw(all_heatmap,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")

for (an in colnames(all_ha_df)) {
    decorate_annotation(an, {
      ## annotation names on the left
      grid.text(
        an,
        unit(1, "npc") + unit(2, "mm"),
        0.6,
        default.units = "npc",
        just = "left",
        gp = gpar(fontsize = 8)
      )
    })
}

dev.off()

## Second heatmap result ===============================================
sig <- readRDS("sig.rds")

set.seed(123)

sig_matrix <- as.matrix(sig)
rownames(sig_matrix) <- toupper(sig_matrix[, "Gene"])
sig_matrix <- sig_matrix[, !colnames(sig_matrix) %in% c("Gene")]
mode(sig_matrix) <- "double"

sig_ha_df <- as.data.frame(sig$cluster)
colnames(sig_ha_df) <- "cluster"

sig_ha_row <-
  rowAnnotation(
    df = sig_ha_df,
    col = list(cluster = setNames(
      topo.colors(9),
      sort(unique(sig$cluster))
    )),
    width = unit(0.2, "in"),
    annotation_legend_param = list(
      nrow = 1,
      title = "Cluster",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      legend_width = unit(1, "in")
    )
  )

sig_heatmap <-
  Heatmap(
    sig_matrix[, !colnames(sig_matrix) %in% c("cluster")],
    top_annotation = all_ha,
    heatmap_legend_param = list(
      title = "Log2 Fold Change",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 6),
      legend_direction = "horizontal",
      legend_width = unit(1.5, "inch")
    ),
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = FALSE,
    split = sig_matrix[, "cluster"],
    gap = unit(1.5, "mm")
  )

png(file = sprintf("%s_sig.png", file_prefix), res = 300, width = 6, height = 6, units = "in")

draw(
  sig_ha_row + sig_heatmap,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)

for (an in colnames(all_ha_df)) {
    decorate_annotation(an, {
      ## annotation names on the left
      grid.text(
        an,
        unit(1, "npc") + unit(2, "mm"),
        0.6,
        default.units = "npc",
        just = "left",
        gp = gpar(fontsize = 8)
      )
    })
  }

dev.off()