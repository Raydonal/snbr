## Arctic lake sediment data (Aitchison, 1986): sand, silt and clay
## proportions of 39 sediment samples at different water depths.
## Source used here: Klaus/artigodaskewnormal/geolog.mat (the file used by
## the original Ox programs), checked on 05/10/2026 against the data set
## ArcticLake of the CRAN package compositions (2.0-9): identical values.
G <- as.matrix(read.table("../../artigodaskewnormal/geolog.mat", skip = 1))
arcticlake <- data.frame(depth = G[, 1] * 100, sand = G[, 2], silt = G[, 3], clay = G[, 4])
arcticlake$logit_clay <- log(arcticlake$clay / (1 - arcticlake$clay))
save(arcticlake, file = "data/arcticlake.rda", compress = "xz")
