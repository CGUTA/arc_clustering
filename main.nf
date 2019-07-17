to_transpose = Channel.create()
ordering_rows_mtx = Channel.create()
ordering_columns_mtx = Channel.create()

input_mtx = Channel.fromPath("${params.matrix}")
				.tap(to_transpose)
				.tap(ordering_rows_mtx)
				.tap(ordering_columns_mtx)
				.map { x -> tuple('cols', x) }

process transpose {
  input:
  file matrix from to_transpose

  output:
  set val("rows"), file('transposed_matrix.tsv') into transposed_mtx

  when:
  params.rows

  script:
  """
  #!/usr/bin/env Rscript

  library(data.table)
  library(magrittr)

  fread("$matrix") %>% melt(id.vars="V1") %>% dcast(variable ~ V1) %>%
  fwrite("transposed_matrix.tsv", sep="\t")

  """
}


process distance_matrix {
	
	input:
	set mode, file(matrix) from input_mtx.mix(transposed_mtx)

	output:
	set mode, file("distance.tsv") into distance

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	essentialome <- fread("$matrix")
	essentialome[,-1] %>% dist %>% as.matrix %>% as.data.table %>% fwrite("distance.tsv", sep="\t")

	"""
}

process clustering {
	
	input:
	set mode, file(distance_matrix) from distance

	output:
	set mode, file("clustering.rds") into clustering

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	clusters <- fread("$distance_matrix") %>% 
	as.matrix %>%
	as.dist %>%
	hclust(method = "ward.D2")


	saveRDS(clusters, "clustering.rds")
	"""
}

ordering_rows = Channel.create()
ordering_columns = Channel.create()

clustering.choice(ordering_rows, ordering_columns) { it[0] === "rows" ? 0 : 1 }

process ordering_rows {
	
	input:
	set mode, file(cluster_data) from ordering_rows
	file matrix from ordering_rows_mtx

	output:
	set mode, file(ordered_matrix)

	when:
	params.rows | params.cols

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	clusters <- readRDS($cluster_data)

	fread("$matrix") %>% .[clusters["order"],] %>%
  	fwrite("ordered_matrix.tsv", sep="\t")

	"""
}

process ordering_columns {
	
	input:
	set mode, file(cluster_data) from ordering_columns
	file matrix from ordering_columns_mtx

	output:
	file ordered_matrix

	when:
	params.rows | params.cols

	"""
	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)

	clusters <- readRDS($cluster_data)

	fread("$matrix") %>% .[,c(1, clusters["order"] + 1)] %>%
  	fwrite("ordered_matrix.tsv", sep="\t")

	"""
}