#' QC_plot Function performing QC by PCA analysis.
#'
#' @param gnmsn_obj path to the matrix (n * m having features i.e. genes, proteins, chipseq peaks in the raws and individual samples in the columns)
#' @param design_path design matrix describing the samples. design matrix must be composed by:  Sample_ID, Sample_Condition and Sample_Replicate columns (layers in Sample_conditions will be defined by "_" separator)
#' @param custom_table vector of Sample_ID names to be used as initial reference for normalisation. If random specified it select randomly 30% of the samples (To be improoved considering global coverage and intensity in order to avoid failed samples) (default NULL)
#' @param saving_path full path where results are saved . default to the working directory
#' @param norm_methods normalisation methods to compare the GNMSN : can be one of "quantile_normalisation", "median_of_ratios_normalization", "uqua", "tmm"
#'
#' @returns Save plots and table to the "saving_pth"/Norm_quality/
#' @export
#'
#' @import grDevices
#' @import graphics
#' @import FactoMineR
#' @import circlize
#' @import ComplexHeatmap
#' @import grid
#' @import methods
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Typical usage
#' param <- BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)
#'
#' Result <- RunNorm(mat,
#'   deg_design,
#'   fix_reference = "random",
#'   row_name_index = 1,
#'   saving_path = out_dir,
#'   n_pop = 1, BiocParam = param
#' )
#'
#' QC_plot(Result, deg_design, saving_path = out_dir)
#' }
#'
QC_plot <- function(gnmsn_obj,
                    design_path,
                    custom_table = NULL,
                    saving_path = NULL,
                    norm_methods = "all") {
  if (!is.null(saving_path)) {
    save <- paste0(saving_path, "/Norm_quality/")
    dir.create(save)
  } else {
    stop("ERROR: to plot the Normalisation pair plots, a path to a general output folder is required")
  }

  norm_mat <- gnmsn_obj[["norm_mat"]]
  ori_mat <- gnmsn_obj[["ori_mat"]]
  scaling_factors <- gnmsn_obj[["scaling_factors"]]

  norm_mat <- norm_mat[rowMeans(norm_mat) > 0, ]
  if (max(norm_mat) > .Machine$integer.max) {
    norm_mat <- norm_mat * ((.Machine$integer.max - 1) / max(norm_mat))
  }

  ori_mat <- ori_mat[rowMeans(ori_mat) > 0, ]
  if (max(ori_mat) > .Machine$integer.max) {
    ori_mat <- ori_mat * ((.Machine$integer.max - 1) / max(ori_mat))
  }

  norm_all <- c("quantile_normalisation", "median_of_ratios_normalization", "uqua", "tmm")
  if (is.null(norm_methods)) {
    cat("|| No additional normalisation to compare\n")
    norm_methods <- NULL
  } else {
    if (norm_methods == "all") {
      norm_methods <- norm_all
    } else {
      if (all(norm_methods %in% norm_all)) {
        norm_methods <- norm_methods
      }
    }
  }

  list_tables <- list()
  list_tables[["GNMSN"]] <- norm_mat

  for (ii in norm_methods) {
    if (ii == "quantile_normalisation") {
      cat("|| Compute quantile_normalisation\n")
      list_tables[[ii]] <- quantile_normalisation(ori_mat)
    }
    if (ii == "median_of_ratios_normalization") {
      cat("|| Compute median_of_ratios_normalization\n")
      list_tables[[ii]] <- median_of_ratios_normalization(ori_mat)
    }
    if (ii == "uqua") {
      cat("|| Compute uqua\n")
      list_tables[[ii]] <- uqua(ori_mat)
    }
    if (ii == "tmm") {
      cat("|| Compute tmm\n")
      list_tables[[ii]] <- tmm(ori_mat)
    }
  }

  if (!is.null(custom_table)) {
    list_tables[["custom_table"]] <- custom_table
  }


  lapply(1:length(list_tables), function(x) {
    nm <- paste0(names(list_tables)[x], ".txt")
    write.table(list_tables[[x]], file = paste0(save, nm), sep = "\t", col.names = NA)
  })


  # To assess quality we compute the PCA() of:
  # rows-centered effects
  # column-centeres effects
  # compute a vst with design ~1 exploiting DESeq2 w/o replicate
  cat("|| VST\n")
  mat_vsts <- lapply(1:length(list_tables), function(y) {
    x <- list_tables[[y]]
    x <- x[rowMeans(x) > 0, ]
    if (max(x) > .Machine$integer.max) {
      x <- x * ((.Machine$integer.max - 1) / max(x))
    }
    vs_mat <- VarianceStabilized(x)
    return(vs_mat)
  })
  names(mat_vsts) <- names(list_tables)
  pca_vsts <- lapply(mat_vsts, function(x) {
    rowwise <- t(scale(t(x), scale = FALSE, center = TRUE))
    # Add the selection of top varying genes:
    set.seed(123456)
    rpca <- FactoMineR::PCA(rowwise, graph = FALSE, scale.unit = FALSE, ncp = ncol(rowwise))
    return(rpca)
  })
  names(pca_vsts) <- names(mat_vsts)
  # Assemble a grip plot for ggplot
  if (is.character(design_path)) {
    design <- read.delim(design_path, sep = "\t", stringsAsFactors = FALSE)
  } else {
    design <- as.data.frame(as.matrix(design_path), stringsAsFactors = FALSE)
  }
  design <- within(design, Condition <- data.frame(do.call("rbind", strsplit(as.character(design$Sample_Condition), "_", fixed = TRUE))))
  design <- design[order(design[, "Condition"][, tail(colnames(design[, "Condition"]), n = 1)]), ]
  design <- design[order(design[, "Condition"][, head(colnames(design[, "Condition"]), n = 1)]), ]
  w <- which(design$Sample_ID %in% colnames(ori_mat))
  design <- design[w, ]
  colnames(design) <- gsub("\\.", "", colnames(design))
  design <- as.data.frame(as.matrix(design), stringsAsFactors = FALSE)
  cn <- c(colnames(design)[grep("^Condition", colnames(design))], colnames(design)[grep("Replicate", colnames(design))])
  for (col in rev(cn)) {
    design <- design[order(design[, col]), ]
  }
  rownames(design) <- paste0(design[rownames(design), ]$Sample_Condition, "_", design[rownames(design), ]$Sample_Replicate)
  # Annotation
  d <- design[, c(grep("^Condition", colnames(design)), c(grep("Sample_Replicate", colnames(design))))]
  d <- as.matrix(d)
  colnames(d) <- gsub("\\.", "", colnames(d))
  col_list <- apply(d, 2, function(x) {
    cc <- circlize::rand_color(length(unique(x)), hue = NULL, luminosity = "bright")
    names(cc) <- unique(x)
    return(cc)
  })

  dd <- 100 / nrow(d)
  row_ha <- ComplexHeatmap::HeatmapAnnotation(df = d, col = col_list, annotation_name_gp = grid::gpar(fontsize = 2), border = TRUE, simple_anno_size = grid::unit(dd, "cm"), show_legend = FALSE)

  # Plots:
  pp_comp <- lapply(1:length(pca_vsts), function(ii) {
    x <- pca_vsts[[ii]]
    eig <- x$eig[1:10, ]
    column_ha <- ComplexHeatmap::rowAnnotation(PCA = ComplexHeatmap::anno_barplot(eig[, 2]))
    ppc <- x$var$coord[(design$Sample_ID), 1:10]
    # add dimension of dots - redundant:
    ppc <- ppc[design$Sample_ID, ]
    rownames(ppc) <- rownames(design)
    ppc_size <- abs(ppc)
    ppc_size <- abs(ppc)
    ppc_size <- sqrt(ppc_size / pi)
    ppc_size <- ppc_size / max(ppc_size)
    write.table(ppc, paste0(save, names(pca_vsts)[ii], "_VST_all_PCA.txt"), sep = "\t", col.names = NA)
    ##
    ppc <- t(ppc)
    ppc_size <- t(ppc_size)
    ##
    col_pal <- c("#1127CC", "white", "#FFB800")
    col_fun <- circlize::colorRamp2(c(min(ppc), 0, max(ppc)), col_pal)
    PCA <- ComplexHeatmap::Heatmap(ppc,
      width = ncol(ppc_size) * grid::unit(dd, "cm"),
      height = nrow(ppc_size) * grid::unit(dd, "cm"),
      col = col_fun, rect_gp = grid::gpar(type = "none"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = NA, fill = NA))
        grid::grid.circle(
          x = x, y = y, r = abs(ppc_size[i, j]) / 2 * grid::unit(dd, "cm"),
          gp = grid::gpar(fill = col_fun(ppc[i, j]), col = NA)
        )
      },
      cluster_rows = FALSE, cluster_columns = FALSE, row_dend_reorder = FALSE, column_dend_reorder = FALSE,
      row_names_gp = grid::gpar(fontsize = 3),
      column_names_gp = grid::gpar(fontsize = 3),
      row_names_side = "left", border = TRUE, top_annotation = row_ha, right_anno = column_ha,
      heatmap_legend_param = list(title = "PC score", just = c("right", "top")), show_heatmap_legend = FALSE,
      show_row_names = TRUE, show_column_names = TRUE
    )

    pdf(paste0(save, names(pca_vsts)[ii], "_VST_all_PCA.pdf"), width = 50, height = 25)
    ComplexHeatmap::draw(PCA, heatmap_legend_side = "bottom")
    dev.off()
  })
}
