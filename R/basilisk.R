##### Using mofapy2 in R via basilisk #####

env_mofapy2 <- basilisk::BasiliskEnvironment(
  envname = "env_mofapy2",
  pkgname = "BIRSBIO2020.scNMTseq.PLS",
  packages = c("mofapy2==0.5.6")
)
