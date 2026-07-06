
library(tiff)
library(EBImage)

fill_hold <- function(x, radius = 5) {
  bw <- x > 0.1                                        # 阈值二值化
  k  <- makeBrush(2 * radius + 1, shape = "disc")      # r=5 -> size 11
  bw <- closing(bw, k)                                 # 二值闭运算(先膨胀后腐蚀)
  bw <- fillHull(bw)                                   # 填补孔洞
  (bw > 0) * 255L                                      # 返回 0/255 整数矩阵
}
islet_marker_addghrelin <- function(sample_id,roi_id,
                                    channel_dir= "Age_channel_img_260104",
                                    out_dir    = "Islet_mask") {


  ## 匹配文件名:sample_id 加 "-" 前缀,避免 "P1" 命中 "P10"


  INS <- readTIFF(file.path(channel_dir, sample_id, roi_id,
                                "141Pr_INS_Mask.tif"), as.is = TRUE)
    GLU <- readTIFF(file.path(channel_dir, sample_id, roi_id,
    "151Eu_GLU_Mask.tif"), as.is = TRUE)
    PP <- readTIFF(file.path(channel_dir, sample_id, roi_id,
    "153Eu_PP_Mask.tif"), as.is = TRUE)
    SST <- readTIFF(file.path(channel_dir, sample_id, roi_id,
    "159Tb_SST_Mask.tif"), as.is = TRUE)

  ## 用 fill_hold 把 ghrelin 通道扩成实心 mask,再叠加到 islet mask
  merged <-fill_hold( INS+GLU+PP+SST)
  merged[merged > 255] <- 255                          # clip 到 uint8 范围

  ## uint8 写出:先归一化到 [0,1] 再指定 8-bit
fname  <- paste(sample_id,roi_id,'islet_mask.tiff',collapse ='-')
  writeImage(Image(merged / 255), file.path(out_dir, fname),
             type = "tiff", bits.per.sample = 8L)
  invisible(fname)
}