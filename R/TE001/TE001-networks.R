##
# TE001 Networks Generation
# -------------------------

# x <- readRDS("~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/TE001_unPruned.rds")
# length(x)
x <- readRDS("~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/TE001-cpm.rds")
.isc_network.tfs <- aracne2regulon( afile = "~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/List-derived/finalNetwork_4col_TFs.tsv" , eset = x , format = "3col" , verbose = TRUE )
length(.isc_network.tfs)
.isc_network.cotfs <- aracne2regulon( afile = "~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/List-derived/finalNetwork_4col_CoTF.tsv" , eset = x , format = "3col" , verbose = TRUE )
length(.isc_network.cotfs)
.isc_network.nucleus <- aracne2regulon( afile = "~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/List-derived/finalNetwork_4col_Nucleus.tsv" , eset = x , format = "3col" , verbose = TRUE )
length(.isc_network.nucleus)
.isc_network.cyto <- aracne2regulon( afile = "~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/List-derived/finalNetwork_4col_Cyto.tsv" , eset = x , format = "3col" , verbose = TRUE )
length(.isc_network.cyto)
.isc_network.mem <- aracne2regulon( afile = "~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/List-derived/finalNetwork_4col_Membrane.tsv" , eset = x , format = "3col" , verbose = TRUE )
length(.isc_network.mem)
.isc_network.ligands <- aracne2regulon( afile = "~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/List-derived/finalNetwork_4col_Ligands.tsv" , eset = x , format = "3col" , verbose = TRUE )
length(.isc_network.ligands)

.isc_network <- c(.isc_network.tfs,.isc_network.cotfs,.isc_network.nucleus,.isc_network.cyto,.isc_network.mem,.isc_network.ligands)
length(.isc_network)

# dim(my_counts)
# TE001.viper <- viper( eset = my_counts , regulon = pruneRegulon(.isc_network,50) , method = "mad" , minsize = 50 , verbose = TRUE )
# saveRDS(TE001.viper,"~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE001_networks/TE001-viper.rds")

# my_counts <- readRDS( file.path("~/Clouds/Dropbox/Data/isc/TE002-cpm.rds") )
# dim(my_counts)
# TE002.viper <- viper( eset = my_counts , regulon = pruneRegulon(.isc_network,50) , method = "mad" , minsize = 50 , verbose = TRUE )
# saveRDS(TE002.viper,"~/Clouds/Dropbox/Califano Wang Collaboration/Data/te-samples/TE002-analysis/TE002-viper.rds")

Olfm4
Cenpf
Krt19