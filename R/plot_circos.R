#' Plot expression of a receptor's ligands by other cell types as a chord plot
#'
#' Creates a chord plot of expression of ligands that can activate a specified
#' receptor where chord widths correspond to mean ligand expression by the cluster.
#'
#' @param dom Domino object that has undergone network building with [build_domino()]
#' @param receptor Name of a receptor active in at least one cell type in the domino object
#' @param ligand_expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to
#'   be rendered between the cell type and the receptor
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @param cell_colors Named vector of color names or hex codes where names correspond to the plotted cell
#'   types and the color values
#' @return Renders a circos plot to the active graphics device
#' @export circos_ligand_receptor
#' @examples
#' data(DominoObjects)
#' dom <- DominoObjects$built_dom_tiny
#' # basic usage
#' circos_ligand_receptor(dom, receptor = "CXCR3")
#' # specify colors
#' cols <- c("red", "orange", "green")
#' names(cols) <- dom_clusters(dom)
#' circos_ligand_receptor(dom, receptor = "CXCR3", cell_colors = cols)
#'
circos_ligand_receptor <- function(
    dom, receptor, ligand_expression_threshold = 0.01, cell_idents = NULL,
    cell_colors = NULL
) {
    check_arg(dom, allow_class = "domino", allow_len = 1)
    check_arg(receptor, allow_class = "character", allow_len = 1)
    check_arg(ligand_expression_threshold, allow_class = "numeric", allow_len = 1, allow_range = c(0, Inf))
    check_arg(cell_idents, allow_class = c("character", "NULL"))
    if (!is.null(cell_colors)) {
        check_arg(cell_colors, allow_class = "character", need_names = TRUE)
    }

    # pull signaling information from the domino result
    ligands <- dom@linkages$rec_lig[[receptor]]

    if (is.null(cell_idents)) {
        # default to all cluster labels in domino object in alphabetical order
        cell_idents <- sort(unique(dom@clusters))
    }

    signaling_df <- obtain_circos_expression(
        dom = dom, receptor = receptor, ligands = ligands,
        ligand_expression_threshold = ligand_expression_threshold,
        cell_idents = cell_idents
    )
    # render circos plot
    render_circos_ligand_receptor(
        signaling_df = signaling_df, receptor = receptor,
        cell_colors = cell_colors,
        ligand_expression_threshold = ligand_expression_threshold
    )
}


#' Obtain Circos Expression
#'
#' Pull expression data from a domino object and format for plotting as a receptor-oriented circos plot.
#'
#' @param dom Domino object that has undergone network building with build_domino()
#' @param receptor Name of a receptor active in at least one cell type in the domino object
#' @param ligands Character vector of ligands capable of interaction with the receptor
#' @param ligand_expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to
#'   be rendered between the cell type and the receptor
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @return a data frame where each row describes plotting parameters of ligand-receptor interactions to pass
#'   to render_circos_ligand_receptor()
#' @keywords internal

obtain_circos_expression <- function(dom, receptor, ligands, ligand_expression_threshold = 0.01, cell_idents = NULL) {
    signaling_df <- NULL
    # obtain expression values from cl_signaling matrices
    active_chk <- vapply(
        dom@linkages$clust_rec,
        FUN.VALUE = logical(1), FUN = function(x) {
            receptor %in% x
        }
    )
    if (!any(active_chk)) {
        stop("No clusters have active ", receptor, " signaling")
    }
    # obtain a signaling matrix where receptor is active
    active_cell <- names(active_chk[active_chk])
    sig <- dom@cl_signaling_matrices[active_cell][[1]]
    cell_names <- gsub("^L_", "", colnames(sig))

    # ensure only ligands present in the signaling matrix are considered
    lig_check <- ligands %in% rownames(sig)
    if (!all(lig_check)) {
        message(
            "Ligands: ", toString(ligands[!lig_check]),
            " of receptor ", receptor, " are listed in the rl_map, but not present in the signaling matrix."
        )
        if (!any(lig_check)) {
            stop("No ligands of receptor ", receptor, " present in signaling matrix.")
        }
        message("Only ligands: ", paste(ligands[lig_check], collapse = ","), " will be considered.")
    }
    ligands <- ligands[lig_check]

    lig_signal_ls <- lapply(
        setNames(ligands, nm = ligands),
        function(l) {
            return(data.frame(
                origin = paste0(cell_names, "-", l),
                destination = receptor,
                mean.expression = unname(sig[rownames(sig) == l, ]),
                sender = cell_names,
                ligand = l,
                receptor = receptor
            ))
        }
    )
    signaling_df <- purrr::list_rbind(lig_signal_ls)


    if (!is.null(cell_idents)) {
        signaling_df <- signaling_df[signaling_df$sender %in% cell_idents, ]
    }

    signaling_df$mean.expression[is.na(signaling_df$mean.expression)] <- 0
    # create a scaled mean expression plot for coord widths greater than 1 by dividing by the max
    # expression [range (0-1)] scaled.mean will only be used when the max expression is > 1
    signaling_df$scaled.mean.expression <- signaling_df$mean.expression / max(signaling_df$mean.expression)
    # exit function if no ligands are expressed above ligand expression threshold
    if (!any(signaling_df[["mean.expression"]] > ligand_expression_threshold)) {
        stop("No ligands of ", receptor, " exceed ligand expression threshold.")
    }
    signaling_df["ligand.arc"] <- 1
    # receptor arc will always sum to 4 no matter how many ligands and cell idents are plotted
    signaling_df["receptor.arc"] <- 4 / (nrow(signaling_df))

    return(signaling_df)
}

#' Render Circos Ligand Receptor Plot
#'
#' Renders a circos plot from the output of [obtain_circos_expression()] to the active graphics device
#'
#' @param signaling_df Data frame output from [obtain_circos_expression()]
#' @param receptor Name of a receptor active in at least one cell type in the domino object
#' @param ligand_expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to be
#'   rendered between the cell type and the receptor
#' @param cell_colors Named vector of color names or hex codes where names correspond to the plotted cell types and
#'   the color values
#' @return a circlize plot is rendered to the active graphics device
#' @keywords internal

render_circos_ligand_receptor <- function(
    signaling_df, receptor, cell_colors = NULL, ligand_expression_threshold = 0.01
) {
    ligands <- sort(unique(signaling_df$ligand))

    # colors for [cell_ident] arcs
    cell_idents <- sort(unique(signaling_df$sender))
    if (is.null(cell_colors)) {
        cell_colors <- ggplot_col_gen(length(cell_idents))
        names(cell_colors) <- cell_idents
    }
    # ensure the vector cell_ident colors is in alphabetical order so that the legend matches the plot
    cell_colors <- cell_colors[sort(names(cell_colors))]

    # chords colored by ligand type
    lig_colors <- ggplot_col_gen(length(ligands))
    names(lig_colors) <- ligands
    origin_cols <- vapply(
        signaling_df$ligand,
        FUN.VALUE = character(1), FUN = function(l) {
            return(lig_colors[l])
        }
    )

    # first index of color vector set to white to hid receptor arc
    grid_col <- c("#FFFFFF", origin_cols)
    names(grid_col) <- c(receptor, signaling_df$origin)

    # name grouping based on [cell_ident]
    l_name_mask <- paste0(paste(paste0("-", ligands), collapse = "|"), "$")
    arc_name <- c(receptor, gsub(l_name_mask, "", signaling_df$origin))
    group <- setNames(arc_name, c(receptor, signaling_df$origin))

    circlize::circos.clear()
    circlize::circos.par(start.degree = 0)
    circlize::chordDiagram(
        signaling_df[, c("origin", "destination", "ligand.arc", "receptor.arc")],
        group = group,
        grid.col = grid_col, link.visible = FALSE,
        annotationTrack = "grid",
        preAllocateTracks = list(
            track.height = circlize::mm_h(4),
            track.margin = c(circlize::mm_h(2), 0)
        ),
        big.gap = 2
    )
    # pre-compute global maxima to avoid repeated calls inside loop
    max_expr <- max(signaling_df[["mean.expression"]])
    for (send in signaling_df$origin) {
        me <- signaling_df[signaling_df$origin == send, ][["mean.expression"]]
        # skip links below threshold
        if (me <= ligand_expression_threshold) next
        if (max_expr > 1) {
            expr <- signaling_df[signaling_df$origin == send, ][["scaled.mean.expression"]]
            max_width <- signif(max_expr, 2)
        } else {
            expr <- me
            max_width <- 1
        }
        circlize::circos.link(send, c(0.5 - (expr / 2), 0.5 + (expr / 2)), receptor, 2, col = paste0(
            grid_col[[send]],
            "88"
        ))
    }
    
    sector_names <- circlize::get.all.sector.index()
    cell_sectors <- cell_idents[cell_idents %in% signaling_df$sender]

    # pick cell sectors based on the start of the sector name being the cell type
    for (cell in cell_sectors) {
        row_pick <- sector_names[startsWith(sector_names, cell)]
        if (length(row_pick)) {
            circlize::highlight.sector(
                sector_names[startsWith(sector_names, cell)],
                track.index = 1, col = cell_colors[cell],
                text = cell, cex = 1, facing = "inside", text.col = "black",
                niceFacing = FALSE, text.vjust = -1.5
            )
        }
    }

    # highlight receptor sector
    circlize::highlight.sector(
        sector_names[startsWith(sector_names, receptor)],
        track.index = 1, col = "#FFFFFF",
        text = receptor, cex = 1.5, facing = "clockwise", text.col = "black",
        niceFacing = TRUE, pos = 4
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
    ComplexHeatmap::draw(lgd_list_vertical, x = grid::unit(0.02, "npc"), y = grid::unit(0.98, "npc"),
        just = c("left", "top"))
}
