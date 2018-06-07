##summarize NF1 results

library(fendR)
require(synapser)
synLogin()

files=synapser::synTableQuery("SELECT distinct `PCSF Result` FROM syn12334021 WHERE ( ( \"NF1 WT\" = '+/+' ) AND ( \"NF1 KO\" = '-/-' ) ) AND mu > 0.004")$asDataFrame()

sapply(unlist(files),function(x){fendR::doAllNetworkAssess(synId=x,tableId="syn12334021")})


files=synapser::synTableQuery("SELECT distinct `PCSF Result` FROM syn12334021 WHERE ( ( \"NF1 WT\" = '+/+' ) AND ( \"NF1 KO\" = '+/-,-/-' ) ) AND mu > 0.004")$asDataFrame()

sapply(unlist(files),function(x){fendR::doAllNetworkAssess(synId=x,tableId="syn12334021")})
