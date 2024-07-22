#' @import grid
#' @import circlize
#' @import ComplexHeatmap
#' @import grDevices
#' 
NULL

#' Plot expression of a receptor's ligands by other cell types as a chord plot
#'
#' Creates a chord plot of expression of ligands that can activate a specified
#' receptor where chord widths correspond to mean ligand expression by the cluster.
#'
#' @param dom Domino object that has undergone network building with [build_domino()]
#' @param receptor Name of a receptor active in at least one cell type in the domino object
#' @param ligand_expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to be rendered between the cell type and the receptor
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @param cell_colors Named vector of color names or hex codes where names correspond to the plotted cell types and the color values
#' @return Renders a circos plot to the active graphics device
#' @export new_circos_ligand_receptor
#' 
new_circos_ligand_receptor <- function(
    dom, receptor, ligand_expression_threshold = 0.01, cell_idents = NULL,
    cell_colors = NULL) {
  arc_df <- circos_lr_shape_data(
    dom = dom, ligand_expression_threshold = ligand_expression_threshold, 
    cell_idents = cell_idents
  )
  circos_lr_plot(
    arc_df = arc_df, ligand_expression_threshold = ligand_expression_threshold, 
    cell_idents = cell_idents, cell_colors = cell_colors
  )
}

#' Create plotting data.frame for a circos plot from a domino object
#'
#' lorem ipsum
#'
#' @param dom Domino object that has undergone network building with [build_domino()]
#' @param receptor Name of a receptor active in at least one cell type in the domino object
#' @param ligand_expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to be rendered between the cell type and the receptor
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @return Data fram describing arcs and chords rendered in the receptor's circos plot
#' @export circos_lr_shape_data
circos_lr_shape_data <- function(
    dom, receptor, ligand_expression_threshold = 0.01, cell_idents = NULL){
  ligands <- dom@linkages$rec_lig[[receptor]]
  signaling_df <- NULL
  if (is.null(cell_idents)) {
    # default to all cluster labels in domino object in alphabetical order
    cell_idents <- sort(unique(dom@clusters))
  }
  # obtain expression values from cl_signaling matrices
  active_chk <- vapply(dom@linkages$clust_rec, FUN.VALUE = logical(1), FUN = function(x) {
    receptor %in% x
  })
  if (sum(active_chk)) {
    # obtain a signaling matrix where receptor is active
    active_cell <- names(active_chk[active_chk == TRUE])
    sig <- dom@cl_signaling_matrices[active_cell][[1]]
    cell_names <- gsub("^L_", "", colnames(sig))
    for (l in ligands) {
      df <- data.frame(origin = paste0(cell_names, "-", l), destination = receptor, mean.expression = unname(sig[rownames(sig) ==
                                                                                                                   l, ]))
      signaling_df <- rbind(signaling_df, df)
    }
  } else {
    stop("No clusters have active ", receptor, " signaling")
  }
  signaling_df$mean.expression[is.na(signaling_df$mean.expression)] <- 0
  # create a scaled mean expression plot for coord widths greater than 1 by dividing by the max
  # expression [range (0-1)] scaled.mean will only be used when the max expression is > 1
  signaling_df$scaled.mean.expression <- signaling_df$mean.expression / max(signaling_df$mean.expression)
  # exit function if no ligands are expressed above ligand expression threshold
  if (sum(signaling_df[["mean.expression"]] > ligand_expression_threshold) == 0) {
    stop("No ligands of ", receptor, " exceed ligand expression threshold.")
  }
  # initialize chord diagram with even ligand arcs
  arc_df <- signaling_df[, c("origin", "destination")]
  arc_df["ligand.arc"] <- 1
  # receptor arc will always sum to 4 no matter how many ligands and cell idents are plotted
  arc_df["receptor.arc"] <- 4 / (nrow(signaling_df))
  return(arc_df)
}



#' Plot expression of a receptor's ligands by other cell types as a chord plot
#'
#' lorem ipsum
#'
#' @param arc_df Data frame descibing arcs and chords
#' @param ligand_expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to be rendered between the cell type and the receptor
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @param cell_colors Named vector of color names or hex codes where names correspond to the plotted cell types and the color values
#' @return Renders a circos plot to the active graphics device
#' @export circos_lr_plot
circos_lr_plot <- function(arc_df, ligand_expression_threshold = 0.01, cell_idents = NULL,
                           cell_colors = NULL){
  # name grouping based on [cell_ident]
  nm <- c(receptor, arc_df$origin)
  group <- structure(c(nm[1], gsub("-.*", "", nm[-1])), names = nm)
  # order group as a factor with the receptor coming first
  group <- factor(group, levels = c(
    receptor, sort(unique(gsub("-.*", "", nm))[-1]) # alphabetical order of the other cell idents
  ))
  # colors for ligand chords
  lig_colors <- ggplot_col_gen(length(ligands))
  names(lig_colors) <- ligands
  # colors for [cell_ident] arcs
  if (is.null(cell_colors)) {
    cell_colors <- ggplot_col_gen(length(cell_idents))
    names(cell_colors) <- cell_idents
  }
  grid_col <- c("#FFFFFF") # hide the arc corresponding to the receptor by coloring white
  for (i in seq_along(ligands)) {
    grid_col <- c(grid_col, rep(lig_colors[i], length(cell_idents)))
  }
  names(grid_col) <- c(receptor, signaling_df$origin)
  circlize::circos.clear()
  circlize::circos.par(start.degree = 0)
  circlize::chordDiagram(arc_df,
                         group = group, grid.col = grid_col, link.visible = FALSE, annotationTrack = c("grid"),
                         preAllocateTracks = list(track.height = circlize::mm_h(4), track.margin = c(circlize::mm_h(2), 0)), big.gap = 2
  )
  for (send in signaling_df$origin) {
    if (signaling_df[signaling_df$origin == send, ][["mean.expression"]] > ligand_expression_threshold) {
      if (max(signaling_df[["mean.expression"]]) > 1) {
        expr <- signaling_df[signaling_df$origin == send, ][["scaled.mean.expression"]]
        max_width <- signif(max(signaling_df[["mean.expression"]]), 2)
      } else {
        expr <- signaling_df[signaling_df$origin == send, ][["mean.expression"]]
        max_width <- 1
      }
      circlize::circos.link(send, c(0.5 - (expr / 2), 0.5 + (expr / 2)), receptor, 2, col = paste0(
        grid_col[[send]],
        "88"
      ))
    }
  }
  sector_names <- circlize::get.all.sector.index()
  cell_sectors <- cell_idents[cell_idents %in% gsub("-.*", "", sector_names)]
  for (cell in cell_sectors) {
    row_pick <- sector_names[grepl(paste0("^", cell), sector_names)]
    if (length(row_pick)) {
      circlize::highlight.sector(sector_names[grepl(paste0("^", cell, "-"), sector_names)],
                                 track.index = 1,
                                 col = cell_colors[cell], text = cell, cex = 1, facing = "inside", text.col = "black",
                                 niceFacing = FALSE, text.vjust = -1.5
      )
    }
  }
  # highlight receptor sector
  circlize::highlight.sector(sector_names[grepl(paste0("^", receptor, "$"), sector_names)],
                             track.index = 1,
                             col = "#FFFFFF", text = receptor, cex = 1.5, facing = "clockwise", text.col = "black", niceFacing = TRUE,
                             pos = 4
  )
  # create legends
  lgd_cells <- ComplexHeatmap::Legend(
    at = as.character(cell_idents), type = "grid", legend_gp = grid::gpar(fill = cell_colors),
    title_position = "topleft", title = "cell identity"
  )
  lgd_ligands <- ComplexHeatmap::Legend(
    at = ligands, type = "grid", legend_gp = grid::gpar(fill = lig_colors), title_position = "topleft",
    title = "ligand"
  )
  chord_width <- 10 / (4 + length(cell_idents) * length(ligands))
  lgd_chord <- ComplexHeatmap::Legend(
    at = c(ligand_expression_threshold, max_width), col_fun = circlize::colorRamp2(c(
      ligand_expression_threshold,
      max_width
    ), c("#DDDDDD", "#DDDDDD")), legend_height = grid::unit(chord_width, "in"), title_position = "topleft",
    title = "ligand expression"
  )
  lgd_list_vertical <- ComplexHeatmap::packLegend(lgd_cells, lgd_ligands, lgd_chord)
  ComplexHeatmap::draw(lgd_list_vertical, x = grid::unit(0.02, "npc"), y = grid::unit(0.98, "npc"), just = c("left", "top"))
}

