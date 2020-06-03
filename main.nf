

Channel.fromPath("${params.matrix}")
.into{ input_mtx ; ordering_mtx}

process transpose {
	
	input:
	file(matrix) from input_mtx

	output:
	file "oriented_matrix.rds" into oriented_mtx

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	essentialome <- fread("$matrix")

	if("$params.cols" != "null"){
		essentialome = essentialome %>% melt(id.vars="V1") %>% dcast(variable ~ V1)
	}

	saveRDS(essentialome, "oriented_matrix.rds")	

	"""
}

oriented_mtx.into{ hclust_mtx; umap_mtx}


process distance_matrix {
	
	input:
	file(matrix) from hclust_mtx

	output:
	file("distance.tsv") into distance

	when:
	params.hclust

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	data <- readRDS("$matrix")

	data[,-1] %>% dist %>% as.matrix %>% as.data.table(keep.rownames=TRUE) %>% fwrite("distance.tsv", sep="\t")

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
	file cluster_data

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


///UMAP
process umap {
	publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${params.id}_$filename" }
	
	input:
	file(matrix) from umap_mtx

	output:
	file "umap.tsv"

	when:
	params.umap

	"""
	#!/usr/bin/env Rscript

	library(uwot)
	library(magrittr)
	library(data.table)

	data <- readRDS("$matrix")

	umap_result <- uwot::umap(data[,-1] %>% as.matrix, $params.uwot_args)
	
	rownames(umap_result) <- data[[1]]

	fwrite(umap_result %>% as.data.table(keep.rownames=TRUE), "umap.tsv")

	"""
}