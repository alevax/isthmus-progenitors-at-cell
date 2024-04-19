## 
# Ayazz Data Analysis
# -------------------

ayyaz_cpm <- readRDS("/Volumes/ac_lab_scratch/lv2395/wang_intestine_sc/C05/c05_cpm.rds")
dim(ayyaz_cpm)
ayyaz_cpm <- convertFeaturesOnMatrix(ayyaz_cpm,"SYMBOL","ENSEMBL","mouse")
rownames(ayyaz_cpm) <- getHumanOrthologousFromMouseEnsemblIds(rownames(ayyaz_cpm))
ayyaz_cpm <- convertFeaturesOnMatrix(ayyaz_cpm,"ENSEMBL","SYMBOL")
ayyaz_s_i <- getStemnessIndex(ayyaz_cpm)
summary(ayyaz_s_i)

kelley_cpm <- readRDS("/Volumes/ac_lab_scratch/lv2395/wang_intestine_sc/KR001_redo/kr001-redo_cpm.rds")
dim(kelley_cpm)
kelley_cpm <- convertFeaturesOnMatrix(kelley_cpm,"SYMBOL","ENSEMBL","mouse")
rownames(kelley_cpm) <- getHumanOrthologousFromMouseEnsemblIds(rownames(kelley_cpm))
kelley_cpm <- convertFeaturesOnMatrix(kelley_cpm,"ENSEMBL","SYMBOL")
kelley_s_i <- getStemnessIndex(kelley_cpm)
kelley_s_i
summary(kelley_s_i)

TE001_cpm <- readRDS("~/Clouds/Dropbox/Data/isc/TE001-cpm.rds")
dim(TE001_cpm)
TE001_cpm <- convertFeaturesOnMatrix(TE001_cpm,"SYMBOL","ENSEMBL","mouse",verbose = T)
rownames(TE001_cpm) <- getHumanOrthologousFromMouseEnsemblIds(rownames(TE001_cpm))
TE001_cpm <- convertFeaturesOnMatrix(TE001_cpm,"ENSEMBL","SYMBOL")
TE001_s_i <- getStemnessIndex(TE001_cpm)
summary(TE001_s_i)

TE002_cpm <- readRDS("~/Clouds/Dropbox/Data/isc/TE002-cpm.rds")
dim(TE002_cpm)
TE002_cpm <- convertFeaturesOnMatrix(TE002_cpm,"SYMBOL","ENSEMBL","mouse",verbose = T)
rownames(TE002_cpm) <- getHumanOrthologousFromMouseEnsemblIds(rownames(TE002_cpm))
TE002_cpm <- convertFeaturesOnMatrix(TE002_cpm,"ENSEMBL","SYMBOL")
TE002_s_i <- getStemnessIndex(TE002_cpm)
summary(TE002_s_i)

N <- 1000
df <- data.frame( s_i = c( sample(ayyaz_s_i,N,replace = T) , 
													 # sample(kelley_s_i,N,replace = T) , 
													 sample(TE001_s_i,N,replace = T) 
													 # sample(TE002_s_i,N,replace = T) 
													 ) ,
									dataset = c( rep("ayyaz",N) , 
															 # rep("kelley",N) , 
															 rep("TE001",N) 
															 # rep("TE002",N) 
															 ) 
									)

p <- ggplot( df , aes(x=s_i, fill=dataset) ) + 
	geom_density(alpha=.3) +
	scale_fill_brewer(palette="Accent",direction = -1) +
	ylab( paste0( "Stemness Index Density") ) +
	# ylab( paste0( "SU2C Patients") ) +
	xlab( paste0( "Stemness Index Score") ) +
	# coord_equal() +
	theme_minimal() +
	theme( legend.position = "right",
		text = element_text(face="italic", colour="black" , size=10 ) ,
		axis.title.x = element_text(face="bold", colour="black" , size=10 ) ,
		axis.title.y = element_text(face="bold", colour="black" , size=10 ) ,
		axis.text.x = element_text(face="italic", colour="black" , angle = 0 ,size = 10 ) ,
		axis.text.y = element_text(face="italic", colour="black" , angle = 0 , size = 10 ) )	

pdf(file.path(reports.dir,"TE001-vs-Ayyazz-stemness-index-on-gene-expression.pdf"))
	print(p)
dev.off()
	
	
	




