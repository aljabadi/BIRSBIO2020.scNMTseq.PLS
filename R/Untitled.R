## workflow commands in R

remotes::install_deps(dependencies = TRUE, repos = BiocManager::repositories())

devtools::load_all()
mae1 <- readRDS('savedata/scnmtseq_gastrulation_mae-sce.rds')
mae2 <- scNMT(dataType = 'mouse_gastrulation', modes = '*', version = '2.0.0', dry.run = FALSE, verbose = FALSE)
