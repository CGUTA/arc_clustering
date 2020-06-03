

Channel.fromPath("${params.matrix}")
.into{ input_mtx ; ordering_mtx}



process distance_matrix {
	
	input:
	file(matrix) from input_mtx

	output:
	file("distance.tsv") into distance

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	essentialome <- fread("$matrix")

	if("$params.cols" != "null"){
		essentialome = essentialome %>% melt(id.vars="V1") %>% dcast(variable ~ V1)
	}

	essentialome[,-1] %>% dist %>% as.matrix %>% as.data.table(keep.rownames=TRUE) %>% fwrite("distance.tsv", sep="\t")

	"""
}

process clustering {
	
	input:
	file(distance_matrix) from distance

	output:
	file("clustering.rds") into clustering

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	clusters <- fread("$distance_matrix") %>% 
	.[, -1] %>%
	as.matrix %>%
	as.dist %>%
	hclust(method = "ward.D2")


	saveRDS(clusters, "clustering.rds")
	"""
}


process ordering{
	publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${params.id}_$filename" }
	
	input:
	file(cluster_data) from clustering
	file matrix from ordering_mtx

	output:
	file("*_ordered_matrix.tsv")

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	clusters <- readRDS("$cluster_data")

	if("$params.cols" != "null"){
		fread("$matrix") %>% .[,c(1, clusters[["order"]] + 1)] %>%
  	fwrite("cols_ordered_matrix.tsv", sep="\t")
	} else{
		fread("$matrix") %>% .[clusters[["order"]],] %>%
  	fwrite("rows_ordered_matrix.tsv", sep="\t")
	}

	

	"""
}
