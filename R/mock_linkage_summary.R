#' Create a mock linkage summary object
#'
#' @return obj a linkage summary object
#' @export
mock_linkage_summary <- function() {
    linkage_sum_tiny <- new("linkage_summary",
        subject_meta = data.frame(
            "subject_names" = paste0("P", 1:6),
            "group" = c(rep("G1", 3), rep("G2", 3))
        ),
        subject_names = factor(
            paste0("P", 1:6),
            levels = paste0("P", 1:6)
        ),
        subject_linkages = list(
            "P1" = list(
                "C1" = list(
                    "tfs" = c("TF1", "TF2", "TF3", "TF4"),
                    "rec" = c("R1", "R2", "R3", "R4"),
                    "incoming_lig" = c("L1", "L2", "L3", "L4"),
                    "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                    "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
                ),
                "C2" = list(
                    "tfs" = c("TF2", "TF3", "TF4"),
                    "rec" = c("R2", "R3", "R4"),
                    "incoming_lig" = c("L2", "L3", "L4"),
                    "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                    "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
                )
            ),
            "P2" = list(
                "C1" = list(
                    "tfs" = c("TF1", "TF2"),
                    "rec" = c("R1", "R2"),
                    "incoming_lig" = c("L1", "L2"),
                    "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
                    "rec_lig" = c("R1 <- L1", "R2 <- L2")
                ),
                "C2" = list(
                    "tfs" = c("TF3", "TF4"),
                    "rec" = c("R3", "R4"),
                    "incoming_lig" = c("L3", "L4"),
                    "tfs_rec" = c("TF3 <- R3", "TF4 <- R4"),
                    "rec_lig" = c("R3 <- L3", "R4 <- L4")
                )
            ),
            "P3" = list(
                "C1" = list(
                    "tfs" = c("TF1", "TF2"),
                    "rec" = c("R1", "R2"),
                    "incoming_lig" = c("L1", "L2"),
                    "tfs_rec" = c("TF1 <- R1", "TF2 <- R2"),
                    "rec_lig" = c("R1 <- L1", "R2 <- L2")
                ),
                "C2" = list(
                    "tfs" = c("TF3"),
                    "rec" = c("R3"),
                    "incoming_lig" = c("L3"),
                    "tfs_rec" = c("TF3 <- R3"),
                    "rec_lig" = c("R3 <- L3")
                )
            ),
            "P4" = list(
                "C1" = list(
                    "tfs" = c("TF2", "TF3", "TF4"),
                    "rec" = c("R2", "R3", "R4"),
                    "incoming_lig" = c("L2", "L3", "L4"),
                    "tfs_rec" = c("TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                    "rec_lig" = c("R2 <- L2", "R3 <- L3", "R4 <- L4")
                ),
                "C2" = list(
                    "tfs" = c("TF1", "TF2", "TF3", "TF4"),
                    "rec" = c("R1", "R2", "R3", "R4"),
                    "incoming_lig" = c("L1", "L2", "L3", "L4"),
                    "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                    "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
                )
            ),
            "P5" = list(
                "C1" = list(
                    "tfs" = c("TF3"),
                    "rec" = c("R3"),
                    "incoming_lig" = c("L3"),
                    "tfs_rec" = c("TF3 <- R3"),
                    "rec_lig" = c("R3 <- L3")
                ),
                "C2" = list(
                    "tfs" = c("TF1", "TF2", "TF3", "TF4"),
                    "rec" = c("R1", "R2", "R3", "R4"),
                    "incoming_lig" = c("L1", "L2", "L3", "L4"),
                    "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3", "TF4 <- R4"),
                    "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3", "R4 <- L4")
                )
            ),
            "P6" = list(
                "C1" = list(
                    "tfs" = c(),
                    "rec" = c(),
                    "incoming_lig" = c(),
                    "tfs_rec" = c(),
                    "rec_lig" = c()
                ),
                "C2" = list(
                    "tfs" = c("TF1", "TF2", "TF3"),
                    "rec" = c("R1", "R2", "R3"),
                    "incoming_lig" = c("L1", "L2", "L3"),
                    "tfs_rec" = c("TF1 <- R1", "TF2 <- R2", "TF3 <- R3"),
                    "rec_lig" = c("R1 <- L1", "R2 <- L2", "R3 <- L3")
                )
            )
        )
    )
    return(linkage_sum_tiny)
}
