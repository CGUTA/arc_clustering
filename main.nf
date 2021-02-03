

Channel.fromPath("${params.matrix}")
.into{ input_mtx ; ordering_mtx}

process transpose {
	
	input:
	file(matrix) from input_mtx

	output:
	file "oriented_matrix.txt" into oriented_mtx

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	essentialome <- fread("$matrix")

	if("$params.orientation" == "column"){
		setnames(essentialome, colnames(essentialome), c("V1", colnames(essentialome)[-1]))
		essentialome = essentialome %>% melt(id.vars="V1") %>% dcast(variable ~ V1)
	}

	fwrite(essentialome, "oriented_matrix.txt", sep = "\t")

	"""
}

oriented_mtx.into{ dist_mtx; umap_mtx}


process distance_matrix {
	publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${params.id}_${params.distance}_$filename" }
	cpus params.cpus
	
	input:
	file(matrix) from dist_mtx

	output:
	file("distance.tsv") into distance

	script:
	if( ["euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"].contains(params.distance) )
		"""
		#!/usr/bin/env Rscript

		library(data.table)
		library(magrittr)

		data <- fread("$matrix", header = TRUE)

		data[,-1] %>% as.matrix %>% dist(method = "${params.distance}") %>% as.matrix %>% as.data.table(keep.rownames=TRUE) %>% fwrite("distance.tsv", sep="\t")

		"""

    else if( ["oneMinusdCor", "oneMinusdCov", "dCor", "dCov"].contains(params.distance) )
        """
        DistanceMatrix -s "$matrix" -o row -d dCor --od . --of distance.tsv --rf false --lf distanceMatrixLog.txt --tc $params.cpus
        """
    else if( ["pearson", "spearman", "kendall"].contains(params.distance) )
        """
        #!/usr/bin/env Rscript

		library(data.table)
		library(magrittr)

		data <- fread("$matrix", header = TRUE)

		cor_matrix = data[,-1] %>% as.matrix %>% cor(method = "${params.distance}")
		colnames(cor_matrix) <- colnames(data[, -1])
		rownames(cor_matrix) <- colnames(data[, -1])
		as.data.table(cor_matrix, keep.rownames=TRUE) %>% fwrite("distance.tsv", sep="\t")
        """

    else
        error "Invalid distance parameter: ${params.distance}"
}

distance.into{ distance_cluster ; distance_print}

process print_distance {
	publishDir "$params.outdir/", mode: 'copy', saveAs: { filename -> "${params.id}_${params.distance}_$filename" }

	input:
	file(distance_matrix) from distance_print

	output:
	file distance_matrix

	when:
	params.print_distance

	"""
	echo printing

	"""
}




process clustering {
	
	input:
	file(distance_matrix) from distance_cluster

	output:
	file("clustering.rds") into clustering

	when:
	params.hclust

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	invert_if_cor <- function(x) {
		if("$params.distance" %in% c("spearman", "pearson", "kendall")){
			1 - x
		} else{
			x
		}	
	}

	clusters <-  fread("$distance_matrix") %>% 
	.[, -1] %>%
	as.matrix %>%
	invert_if_cor %>%
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

	if("$params.orientation" == "column"){
		fread("$matrix") %>% .[,c(1, clusters[["order"]] + 1), with=FALSE] %>%
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

	data <- fread("$matrix")

	umap_result <- uwot::umap(data[,-1] %>% as.matrix, $params.uwot_args)
	
	rownames(umap_result) <- data[[1]]

	fwrite(umap_result %>% as.data.table(keep.rownames=TRUE), "umap.tsv")

	"""
}
