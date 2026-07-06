quantify_cells <- function(mask_file, marker_files) {
  mask <- tiff::readTIFF(mask_file, as.is = TRUE)
  lab  <- as.vector(mask); keep <- lab != 0
  g    <- lab[keep]

  area <- rowsum(rep(1L, length(g)), g)                        # 面积
  x    <- rowsum(as.vector(col(mask))[keep], g) / area         # 质心 x
  y    <- rowsum(as.vector(row(mask))[keep], g) / area         # 质心 y

  out <- data.frame(cell_id = as.integer(rownames(area)),
                    x = x[,1], y = y[,1], area = area[,1],
                    row.names = NULL)

  for (m in names(marker_files)) {
    img <- as.matrix(RBioFormats::read.image(marker_files[[m]], normalize = FALSE))
    stopifnot(identical(dim(img), dim(mask)))
    out[[paste0(m, "_mean")]] <- (rowsum(as.vector(img)[keep], g) / area)[, 1]
  }
  out
}

## 用法
res <- quantify_cells(
  cellpose_mask_tif_file,
  list(INS = INS_ome_file, GLU = GLU_ome_file,Delta = SST_ome_file,PP = PP_ome_file, epsilon =Ghrline_ome_file)
)