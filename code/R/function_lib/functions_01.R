
get_btx_palette <- function() {
	BTX_COL <-
		list(
			blue_d_in = c(189, 197, 255),
			blue_d_out = c(0, 62, 193),
			blue_l_in = c(179, 223, 255),
			blue_l_out = c(0, 151, 229),
			green_in = c(195, 249, 195),
			green_out = c(0, 169, 0),
			pink_in = c(254, 195, 252),
			pink_out = c(184, 56, 180),
			purple_in = c(212, 178, 255),
			purple_out = c(151, 73, 254),
			gray_in = c(214, 214, 214),
			gray_out = c(0, 0, 0)	
		)
	tmp <- rep(NA, length(BTX_COL))	
	for (i in 1:length(BTX_COL)) {
	  tmp[i] <- 
		rgb(
			red = (BTX_COL[[i]][1]), 
			green = (BTX_COL[[i]][2]), 
			blue = (BTX_COL[[i]][3]),
			maxColorValue = 255
		)
	}
	names(tmp) <- names(BTX_COL)
	BTX_COL <- as.list(tmp)	
	return(BTX_COL)

}

get_palette <- function(annot, max_val = NULL) {
	#paletteer::palettes_d_names |> dplyr::filter(length == length(phase_cols) & type == "sequential") |> View()

	## Target Type --------------------------------------------------------
	BTX_COL <- get_btx_palette()
	color_type <-
		c(
			"amplicon" = BTX_COL$green_out,             
			"portfolio" = BTX_COL$blue_d_in,             
			"positive_control_adc" = BTX_COL$blue_d_out, 
			"positive_control_radio" = BTX_COL$blue_l_out,
			"additional_targets" = BTX_COL$pink_out 
		)
	## Target Phase -------------------------------------------------------
	color_phase <- levels(annot$`Highest Phase of Development`)
	color_phase <- paletteer::paletteer_d("MoMAColors::Ernst", n = length(color_phase))
	names(color_phase) <- levels(annot$`Highest Phase of Development`)
	## Target Level -------------------------------------------------------
	if (!is.null(max_val)) {
		hokusai_heatmap <- as.character(paletteer::paletteer_d("MetBrewer::Hokusai3"))		
		breaks <- seq(0, round(max_val), length.out = length(hokusai_heatmap))
		hokusai_heatmap <- circlize::colorRamp2(breaks = breaks, colors = hokusai_heatmap)		
		
		tealrose_heatmap <- as.character(paletteer::paletteer_c(palette = "grDevices::TealRose", direction = -1, n = 300))
		breaks_1 <- seq(0, round(max_val), length.out = length(tealrose_heatmap))
		tealrose_heatmap <- circlize::colorRamp2(breaks = breaks_1, colors = tealrose_heatmap)

		# grDevices::Blues 3
		# grDevices::Purples 2
		# grDevices::Purple-Blue
		purple_heatmap <- as.character(paletteer::paletteer_c(palette = "grDevices::Purples 3", direction = -1, n = 300))		
		breaks <- seq(0, round(max_val), length.out = length(purple_heatmap))
		purple_heatmap <- circlize::colorRamp2(breaks = breaks, colors = purple_heatmap)				
	}else{
		hokusai_heatmap = NULL
		tealrose_heatmap = NULL
		purple_heatmap = NULL
	}
	## Amplicon incidence
	amp_cols <-	circlize::colorRamp2(breaks = c(0, 100), colors = c(BTX_COL$gray_in, BTX_COL$gray_out))

	return(list(color_type = color_type, color_phase = color_phase, hokusai_heatmap = hokusai_heatmap, tealrose_heatmap = tealrose_heatmap, purple_heatmap = purple_heatmap, amp_cols = amp_cols))	
}
	
write_heatmap <- function(x, file_out) {
	svglite::svglite(file_out, width = 32.47, height = 13.1)
	ComplexHeatmap::draw(
		x, 
		heatmap_legend_side = "top", 
		annotation_legend_side = "top"
	)
	invisible(dev.off())
	return(file_out)
}

write_svg <- function(x, file_out, w = 32.47, h = 13.1) {
	svglite::svglite(file_out, width = w, height = h)
	plot(x)
	invisible(dev.off())
	return(file_out)
}

get_target_list <- function(filepath) {
	out <- 
		readxl::read_excel(filepath) |>
		dplyr::filter(Status != "3-Reserved (AMP)") |>
		tidyr::separate_rows(ENSGID, sep = ";")
	return(out)

}


prepare_ensembl_ids <- function(x) {

	# Convert from gene.symbol to ensembl.gene
	library(EnsDb.Hsapiens.v86)
	output <- 
		ensembldb::select(
			EnsDb.Hsapiens.v86, 
			keys = x, 
			keytype = "SYMBOL", 
			columns = columns(EnsDb.Hsapiens.v86)
		) |>
		dplyr::filter(grepl("ENSG", GENEID)) |>
		dplyr::filter(TXBIOTYPE == "protein_coding") |>
		dplyr::select(GENEID, SYMBOL, UNIPROTID) |>
		dplyr::distinct()

	# Test outupt
	if(all(x %in% output$SYMBOL)) {
		return(output)
	}else{
		logger::log_warn("Missing ensembl ID.")
		return(output)
	}
}	

merge_targets <- function(x, y) {

	# tar_load(current_targets); x = current_targets; tar_load(references_ensgid); y = references_ensgid
	out <- 
		dplyr::full_join(x, y, by = c("SYMBOL", "ENSGID" = "GENEID")) |>
		dplyr::select(SYMBOL, ENSGID) |>
		dplyr::distinct(.keep_all = TRUE)
	return(out)

}

get_hpa_filter <- function(x, y, p) {

	# tar_load(normal_ihc_data); x = normal_ihc_data; tar_load(param_ls); p = param_ls; tar_load(current_targets); y = current_targets
	out <- 
		x |>
		dplyr::filter(grepl(paste0(p$core_tissue, collapse = "|"), Tissue)) |>
		dplyr::filter(Gene %in% y$ENSGID)
	return(out)
}

hpa_download <- function(url_str){
	#' Download and load Human Protein Atlas data.
	#'
	#' This function downloads a data file from the Human Protein Atlas (HPA)
	#' and loads it into an R data frame. For details, see:
	#' https://www.proteinatlas.org/about/download
	#'
	#' @param url_str (string) The name of the file to download from the HPA.
	#'   Defaults to "rna_cancer_sample.tsv.gz".
	#' @param target_ids (data.frame) Output of get_ensembl_id()
	#'
	#' @return A data frame containing the loaded data, or NULL if the download
	#'   or loading fails.
	#'
	#' @examples
	#' \dontrun{
	#'   hpa_data <- hpa_download("rna_cancer_sample.tsv.gz")
	#'   if (!is.null(hpa_data)) {
	#'     print(head(hpa_data))
	#'   }
	#' }
	#'
	#' @export

	# Check input arguments
	checkmate::assertCharacter(
		url_str, 
		len = 1, 
		any.missing = FALSE, 
		null.ok = FALSE
	)	
		
	# Define paths
	url_id <- 
		stringr::str_glue("https://www.proteinatlas.org/download/tsv/{url_str}")
	dest_id <- 
		stringr::str_glue("../../input/hpa/{url_str}")
	
	# Download file if it doesn't exist
	if(file.exists(dest_id)) {
		logger::log_info("Download exists:", dest_id)
	}else{
		tryCatch({
		  download.file(url_id, destfile = dest_id)
		  logger::log_info("Downloaded data to:", dest_id)
		}, error = function(e) {
		  logger::log_error("Failed to download data:", e$message)
		  return(NULL) # Return NULL if download fails
		})		
	}
	
	# Load files
	# TODO filter to genes of interest on the fly to save memory 
	# TODO Save to rds and remove zip file, and check for rds
	tryCatch({
		input_hpa <- readr::read_tsv(dest_id)
	}, error = function(e) {
		logger::log_error("Failed to read data:", e$message)
		return(NULL) # Return NULL if reading fails
	})
	
	# Test output
	if (!exists("input_hpa") || is.null(input_hpa)) {
		logger::log_warn("No data to return.")
		return(NULL)
	}else{
		return(input_hpa)
	}


}

get_hpa_overview <- function() {
	library(HPAanalyze) # devtools version
	# Retrieve all genes/data available in HPA
	histology_tb <- 
		HPAanalyze::hpaDownload(downloadList = 'histology')
	return(histology_tb)
}

get_ensgids <- function(x, y) {
	library(HPAanalyze) # devtools version
	get_me <- 
		y$normal_tissue |>
		dplyr::filter(tissue == "Kidney")

	# Filter to those with membrane expression (for comp. tractability)
	subcellular <- 
		y$subcellular_location |>
		dplyr::filter(
			grepl("plasma membrane|junction|Focal adhesion sites", main_location, ignore.case = TRUE) | 
			grepl("plasma membrane|junction|Focal adhesion sites", additional_location, ignore.case = TRUE)
		)
	get_me <- 
		get_me |>
		dplyr::filter(ensembl %in% subcellular$ensembl)
	# Filter for pilot study
	if (FALSE) {
		get_me <-
			get_me |>
			dplyr::filter(gene %in% c("DDR1", "ST14", "ITGA3", "NECTIN4", "ERBB2", "FOLH1", "SSTR2"))	
	}

	ensgids <- unique(c(get_me$ensembl, x$ENSGID)) # Add in any from internal portfolio
	ensgids <- as.list(ensgids) # NOTE: Subset for debugging [1:10]
	return(ensgids)
}

get_image_path <- function(ensgid) {
	library(HPAanalyze)
	library(dplyr)
	logger::log_info(ensgid)
	xml_tmp <- 
		tryCatch({
		  HPAanalyze::hpaXmlGet(ensgid) # XML nodes
		}, error = function(e) {
		  logger::log_error("Failed to download data:", e$message)
		  return(NULL) # Return NULL if download fails
		})
	if(!is.null(xml_tmp)) {
		ab_ls <- HPAanalyze::hpaXmlAntibody(xml_tmp) # Antibody info
		expr_ls <- HPAanalyze::hpaXmlTissueExpr(xml_tmp) # expression details; URLs; both healthy and cancer 
		names(expr_ls) <- ab_ls$id
		if(sum(unlist(lapply(X = expr_ls, FUN = nrow))) > 0) {
			ensgid_tb <-
				expr_ls |>
				dplyr::bind_rows(.id = "ab_id") |>
				dplyr::mutate(ENSGID = ensgid) |>
				dplyr::filter(snomedCode1 == "M-00100") |> 
				dplyr::filter(tissueDescription2 %in% "Kidney") # filter to all 'Normal tissue' 				
			return(ensgid_tb)		
		}else{
			return(NULL)
		}
	}else{
		return(NULL)
	}
}

get_img_dir <- function(x) {

	logger::log_info("Image files:", length(x))
	return("../../input/hpa/ihc/")

}

get_url <- function(x) {
	# tar_load(ensgids_paths_tb); x = ensgids_paths_tb$imageUrl[1]
	return(as.list(ensgids_paths_tb$imageUrl))
}

get_ihc <- function(x) {
	dest_id <- file.path("../../input/hpa/ihc/", basename(unlist(x)))
	# Download file if it doesn't exist
	if(file.exists(dest_id)) {
		logger::log_info("Download exists:", dest_id)
	}else{
		tryCatch({
		  download.file(x, destfile = dest_id, mode = "wb")
		  logger::log_info("Downloaded data to:", dest_id)
		}, error = function(e) {
		  logger::log_error("Failed to download data:", e$message)
		  return(NULL) # Return NULL if download fails
		})		
	}
	# Process the file
	#img_hist_ls[[i]] <- process_ihc_img(path = dest_id)
	return(dest_id)
}

get_binned_hist <- function(path) {
	library(colordistance)
	# Process image
	# Aim: Quantize -> kmeans -> count
	# Load an image
	# img <- 
		# magick::image_read(dest_id) |> 
		# image_quantize(max = 30)
	# First, off-the-shelf
	# clusters_km <- 
		# colordistance::getKMeanColors(
			# dest_id, 
			# n = 6, 
			# color.space = "hsv",
			# plotting = TRUE
		# ) # Not super automatable downstream
	#cluster_tb <- colordistance::extractClusters(clusters_km)
	#colordistance::plotClusters(cluster_tb)
	# This should work irrespective of expected number of clusters
	# Computes a histogram in RGB/HSV colorspace by sorting pixels into a specified number of bins.
	binned_hist <- 
		colordistance::getImageHist(
			path, 
			bins = 3, 
			hsv = TRUE, 
			plotting = FALSE
		) |>
	dplyr::mutate(path = basename(path))
	return(binned_hist)
}

clean_img_transparent <- function(img_path, dest_dir = "../../input/hpa/ihc/"){
  	# Aim: given an image, select colour in top right-hand corner and remove
	library(magick) 

	# Ensure the data directory exists
	if(!dir.exists(dest_dir)) {
		dir.create(dest_dir, showWarnings = FALSE)
	}
	dest_id <- file.path("../../input/hpa/ihc/", basename(unlist(img_path)))
	# Process file if it doesn't exist
	if(file.exists(dest_id)) {
		logger::log_info("Download exists:", dest_id)
		return(dest_id)
	}else{
		tryCatch({
			#download.file(x, destfile = dest_id, mode = "wb")
			img <-  
				magick::image_read(img_path) |>
				magick::image_contrast(sharpen = 1) |>
				magick::image_fill(
					color = "transparent", 
					refcolor = "white", 
					fuzz = 4,
					point = "+1+1" # start at top left 1 pixel 
				) 
			# A future iteration will sample color from first pixel; hard to translate
			# top_left_pixel <- image_crop(img, "1x1+0+0")
			# pixel_data <- image_data(top_left_pixel, channels = "rgb")[ , , 1]
			# test2 = image_transparent(tmp, color = pixel_data) 
			magick::image_write(img, path = dest_id, format = "jpg") # TODO: JPEG does not handle transparency! Use PNG or set black background later (opted for this)
			image_destroy(img)
			return(dest_id)		  
		  logger::log_info("Downloaded data to:", dest_id)
		}, error = function(e) {
		  logger::log_error("Failed to download data:", e$message)
		  return(NULL) # Return NULL if download fails
		})		
	}
}

get_hist_dist <- function(x, y) {

	# tar_load(image_hist_ls); x = image_hist_ls; tar_load(ensgids_paths); y = ensgids_paths
	
	# Aim: Average histograms by antibody, and then calculate earth mover distance between antibodies
	
	mappings <- 
		y |>
		dplyr::select(imageUrl, ENSGID, ab_id) |>
		dplyr::mutate(file_n = basename(imageUrl))
	### --- combineClusters at the antibody level ---		
	ab_ids <- unique(mappings$ab_id)
	out_ls <- vector("list", length(ab_ids))
	names(out_ls) <- ab_ids
	for (id in ab_ids) {
		find_me <-
			mappings |>
			dplyr::filter(ab_id == id) |>
			dplyr::pull(file_n)
		test <- all(unlist(lapply(X = x[find_me], FUN = is.null)))
		if (!test) {
			out_ls[[id]] <- 
				colordistance::combineList(x[find_me], method = "mean")
		}
		
	}
	
	cdm <- 
		colordistance::getColorDistanceMatrix(
			Filter(Negate(is.null), out_ls), 
			method = "emd", 
			plotting = FALSE
		)
	return(cdm)
}


get_heatmap <- function(x, y, p, g, a) {
	
	# tar_load(img_dist); x = img_dist; tar_load(ensgids_paths); y = ensgids_paths; tar_load(param_ls); p = param_ls; tar_load(hpa_overview_ls); g = hpa_overview_ls; tar_load(all_targets); a = all_targets
	
	library(ggplot2)
	library(ComplexHeatmap)
	library(dplyr)

	# Note cell_type annotation is at gene level, not antibody level
	# This means the heatmap may flag anitbodies within a protein target that do not align with the single label (in heatmap). On website, annotations are better, but not programatically accesible.
	normal_annotations <-
		g$normal_tissue |> 
		dplyr::filter(tissue == "Kidney") |>
		dplyr::mutate(level = na_if(level, "N/A")) |>
		dplyr::mutate(level = factor(level, levels = c("Not representative", "Not detected", "Descending", "Ascending", "Low", "Medium", "High"))) |>
		dplyr::group_by(ensembl) |>
		tidyr::pivot_wider(names_from = cell_type, values_from = level) |>
		dplyr::select(!gene:reliability)

	mappings <- 
		y |>
		dplyr::select(ENSGID, ab_id) |>
		tidyr::unnest_longer(ENSGID) |>
		dplyr::distinct() |>
		dplyr::left_join(g$subcellular_location, by = c("ENSGID" = "ensembl")) |>
		dplyr::mutate(highlight = gene %in% p$focus_targets) |>
		dplyr::left_join(normal_annotations, by = c("ENSGID" = "ensembl"))

	## Extract matrix ---------------------------------------------------------
	m <- x
	
	## --- Colours ------------------------------------------------------------
	BTX_COL <- get_btx_palette()

	## --- Row annotations ---
	# Kidney architecture staining
	cell_type <- unique(g$normal_tissue$cell_type)
	cell_type <-
		mappings |>
		dplyr::select(ab_id, contains(cell_type))
	cell_type <-
		cell_type[match(colnames(m), cell_type$ab_id), ] |>
		dplyr::select(!ab_id) |>
		as.matrix()
	tealrose_heatmap <- as.character(paletteer::paletteer_c(palette = "grDevices::TealRose", direction = 1, n = 4))
	names(tealrose_heatmap) <- c("Not detected", "Low", "Medium", "High")
	right_ha <- 
		ComplexHeatmap::rowAnnotation(
			cell_type = cell_type,
			col = list(cell_type = tealrose_heatmap)
		)
			
	## --- Top Column annotations ---
	## Target charateristics
	top_ha <- 
		ComplexHeatmap::HeatmapAnnotation(
			portfolio = as.character(mappings$highlight[match(colnames(m), mappings$ab_id)]),
			annotation_name_side = "left",
			col = list(portfolio = c("TRUE" = BTX_COL$blue_d_out, "FALSE" = BTX_COL$gray_in))
		)
		
	## --- Bottom Column annotations ---
	## Antibody charateristics	
	reliability_col <- as.character(paletteer::paletteer_d("lisa::GeorgiaOKeeffe", direction = 1, n = 4))
	names(reliability_col) <- c("Enhanced", "Supported", "Approved", "Uncertain")
	bottom_ha <- 
		ComplexHeatmap::HeatmapAnnotation(
			reliability = mappings$reliability[match(colnames(m), mappings$ab_id)], 
			single_cell_var_intensity = mappings$single_cell_var_intensity[match(colnames(m), mappings$ab_id)], 
			single_cell_var_spatial = mappings$single_cell_var_spatial[match(colnames(m), mappings$ab_id)],
			annotation_name_side = "left",
			col = list(reliability = reliability_col)			
		)

	## --- Update names ---
	# These are not unique, so need to keep antibody identifiers till end
	colnames(m) <- mappings$gene[match(colnames(m), mappings$ab_id)]
	rownames(m) <- mappings$gene[match(rownames(m), mappings$ab_id)]

	fh = function(x) fastcluster::hclust(dist(x))
	a_plt <-
		ComplexHeatmap::Heatmap(
			m, 
			col = colorRampPalette(c("royalblue4", "ghostwhite", "violetred2"))(299),
			# col = paletteer::paletteer_c("grDevices::Purples 3", 299),
			cluster_columns = fh, 
			cluster_rows = fh, 
			row_dend_reorder = TRUE,
			column_dend_reorder = TRUE,
			# name = "Confidence",
			# heatmap_legend_param = list(
				# title = "IHC level", at = c(0:4), 
				# labels = c("NA", "ND", "L", "M", "H"), 
				# direction = "horizontal"
			# ),
			na_col = "white",
			top_annotation = top_ha, 
			bottom_annotation = bottom_ha, 
			right_annotation = right_ha,
			cluster_row_slices = TRUE,
			row_gap = unit(5, "mm"),
			row_title = NULL,
			column_title = NULL,
			row_names_gp = gpar(fontsize = 18),
			column_names_rot = 45,
			# row_names_rot = 45,
			column_names_gp = gpar(fontsize = 18, hjust = 1),
			rect_gp = gpar(col = "white", lwd = 1),			
			show_row_dend = TRUE,
			show_column_dend = TRUE,
			show_row_names = TRUE,
			show_column_names = TRUE
		)

	return(a_plt)
}


get_stacked_barplot <- function(x, y) {

	# tar_load(image_hist_ls); x = image_hist_ls; tar_load(ensgids_paths); y = ensgids_paths
	
	# Aim: Average histograms by antibody, and then calculate earth mover distance between antibodies
	
	mappings <- 
		y |>
		dplyr::select(imageUrl, ENSGID, ab_id) |>
		dplyr::mutate(file_n = basename(imageUrl))
	### --- combineClusters at the antibody level ---		
	ab_ids <- unique(mappings$ab_id)
	out_ls <- vector("list", length(ab_ids))
	names(out_ls) <- ab_ids
	for (id in ab_ids) {
		find_me <-
			mappings |>
			dplyr::filter(ab_id == id) |>
			dplyr::pull(file_n)
		test <- all(unlist(lapply(X = x[find_me], FUN = is.null)))
		if (!test) {
			out_ls[[id]] <- 
				colordistance::combineList(x[find_me], method = "mean")
		}
		
	}
	
    # Generate hex colors for each pixel
    lambda <- function(img) {
		grDevices::rgb(
			suppressMessages(
				colordistance::convertColorSpace(
					from = "Lab",
					to = "sRGB", 
					color.coordinate.matrix = img, 
					sample.size = "all", 
					from.ref.white = "D65"
				)
			)
		)
		out <- 
			img |>
			dplyr::mutate(rgb_cols = rgb_cols) 
		return(out)		
	}
	to_plot <- lapply(X = out_ls, FUN = lambda)
	names(to_plot) <- names(out_ls)
	to_plot <-
		to_plot |>
		dplyr::bind_rows(.id = "antibody")
	to_plot <-
		to_plot |>
		dplyr::mutate(fill_key = interaction(antibody, rgb_cols))
	
	# Rank order by non-staining colour
	od <- 
		to_plot |>
		dplyr::filter(rgb_cols == "#D8CCCD") |>
		dplyr::arrange(Pct) |>
		dplyr::pull(antibody)
	to_plot <-
		to_plot |>
		dplyr::mutate(antibody = factor(antibody, levels = od))
	
	# IHC colors
	color_palette <- setNames(to_plot$rgb_cols, to_plot$fill_key)
	
	# TODO: Flip coordinates (Label: Target | Ab)
	# TODO: only label targets of interest
	a_plot <- 
		to_plot |>
			ggplot(aes(fill = fill_key, y = Pct, x = antibody)) + 
			geom_bar(position="stack", stat="identity")	+
			scale_fill_manual(values = color_palette) +	
			theme(legend.position = "none")	  
	return(a_plot)
}




###############################################################################
# Experiments

# Define the image transformations
preprocess <- function(img) {
  img |>
    # Convert the image to a tensor
    transform_to_tensor() |>
    # Resize to the expected input size of the model (e.g., 224x224 for ResNet)
    transform_resize(size = c(224, 224)) |>
    # Normalize with the mean and standard deviation of the ImageNet dataset
    transform_normalize(mean = c(0.485, 0.456, 0.406), std = c(0.229, 0.224, 0.225)) |>
    # Add a batch dimension
    `$`unsqueeze(1)
}


get_torch <- function() {

	library(torch)
	library(torchvision)

	# Set the device to CPU
	device <- torch_device("cpu")
	# Load a pre-trained ResNet-18 model
	model <- models_resnet18(pretrained = TRUE)
	# Remove the final classification layer (fc) to get embeddings
	model$fc <- nn_identity()
	# Move the model to the CPU and set it to evaluation mode
	model$to(device = device)
	model$eval()
	# Example list of image file paths
	image_files <- list.files("./tmp/", full.names = TRUE, pattern = "\\.jpg$|\\.png$")
	# A list to store the embeddings
	embeddings_list <- list()
	# Loop through each image file
	for (file_path in image_files) {
	  # Load the image
	  img <- magick::image_read(file_path)
	  # Preprocess the image
	  img_tensor <- preprocess(img)
	  # Generate the embedding (no gradient calculation needed)
	  with_no_grad({
		embedding <- model(img_tensor$to(device = device))
	  })
	  # Store the embedding (move it back to the CPU if it were on a GPU)
	  embeddings_list[[file_path]] <- as.array(embedding$cpu())
	}
	# Combine all embeddings into a single matrix
	embeddings_matrix <- do.call(rbind, embeddings_list)



}


get_hpa_xml <- function(x) {
	library(HPAanalyze) # devtools version
	library(progress) # devtools version

	# tar_load(all_targets); x = all_targets

	# Retrieve all genes/data available in HPA
	# histology_tb <- 
		# HPAanalyze::hpaDownload(downloadList = 'histology')
	
	## Loop through each gene (must parallelise to save ~2d)
	# parallelisation at function level not working in targets
	# out_ls <- vector("list", length(ensgids)) 
	# names(out_ls) <- get_me$ensgids
	logger::log_info("Processing ENSGIDS:", length(ensgids))
		  
	lambda <- function(ensgid) {
		library(HPAanalyze)
		library(dplyr)
		logger::log_info(ensgid)
		xml_tmp <- 
			tryCatch({
			  HPAanalyze::hpaXmlGet(ensgid) # XML nodes
			}, error = function(e) {
			  logger::log_error("Failed to download data:", e$message)
			  return(NULL) # Return NULL if download fails
			})
		if(!is.null(xml_tmp)) {
			ab_ls <- HPAanalyze::hpaXmlAntibody(xml_tmp) # Antibody info
			expr_ls <- HPAanalyze::hpaXmlTissueExpr(xml_tmp) # expression details; URLs; both healthy and cancer (verified with NECTIN4) 
			names(expr_ls) <- ab_ls$id
			ensgid_tb <-
				expr_ls |>
				dplyr::bind_rows(.id = "ab_id") |>
				dplyr::mutate(ENSGID = ensgid) |>
				dplyr::filter(snomedCode1 == "M-00100") |> 
				dplyr::filter(tissueDescription2 %in% "Kidney") # filter to all 'Normal tissue' 				
			# expr_sum <- HPAanalyze::hpaXmlTissueExprSum(xml_tmp)	
		}
		return(ensgid_tb)
	}
	
	library(parallel)
	numCores <- parallel::detectCores()

	if(FALSE) {
		library(parabar)
		backend <- 
			parabar::start_backend(
				cores = numCores/2, 
				cluster_type = "psock", 
				backend_type = "async"
			)
		out_ls <- parabar::par_lapply(backend, x = ensgids, fun = lambda)
		parabar::stop_backend(backend)
		# parabar::clear(backend)	
	}else{
		cl <- parallel::makeCluster(numCores/2, type = "PSOCK")
		# 2. Load necessary packages on each worker
		parallel::clusterEvalQ(cl, {
			library(targets)
			library(HPAanalyze)
			library(dplyr)
			library(logger)
		})
		# 3. Export the specific function you want to run to each worker
		parallel::clusterExport(cl, varlist = "lambda", envir = environment())
		out_ls <- parallel::parLapply(cl, X = ensgids, fun = lambda)
		parallel::stopCluster(cl)	
	}

	out_tb <-
		out_ls |>
		dplyr::bind_rows() 

  return(out_tb)
}


get_hpa_healthy_img <- function(x) {

	# tar_load(hpa_ihc_url_tb); x = hpa_ihc_url_tb
	library(HPAanalyze) # devtools version

	pb <- 
		progress::progress_bar$new(
			format = "  downloading [:bar] :percent eta: :eta",
			total = nrow(img_tb), 
			clear = FALSE, 
			width= 60
		  )
	img_hist_ls <- list("vector", nrow(img_tb))
	for (i in 1:nrow(img_tb)) {
		pb$tick()
		url_id <- img_tb$imageUrl[i]
		dest_id <- file.path("../../input/hpa/ihc/", basename(url_id))
		# Download file if it doesn't exist
		if(file.exists(dest_id)) {
			logger::log_info("Download exists:", dest_id)
		}else{
			tryCatch({
			  download.file(url_id, destfile = dest_id, mode = "wb")
			  logger::log_info("Downloaded data to:", dest_id)
			}, error = function(e) {
			  logger::log_error("Failed to download data:", e$message)
			  return(NULL) # Return NULL if download fails
			})		
		}

		img_hist_ls[[i]] <- process_ihc_img(path = dest_id)		
	}

  return(img_hist_ls)
}


sort_hpa_healthy_img <- function(x, y) {

	# TODO: PLot stacked barplots for each IHC image
	# Return Gini 

	# tar_load(hpa_ihc_hist_ls); x = hpa_ihc_hist_ls; tar_load(hpa_ihc_url_tb); y = hpa_ihc_url_tb
	for(i in 1:length(x)) {
	
		to_plot <-
			x[[i]] |> 
			dplyr::rowwise() |>
			dplyr::mutate(hex = hsv(h = h, s = s, v = v)) |>
			dplyr::ungroup() |>
			dplyr::arrange(h, s, v) |>
			dplyr::mutate(bin = as.character(row_number()))
	
		colour_mapping <- setNames(unique(to_plot$hex), unique(to_plot$bin))
		
		# TODO: add in all associated metadata for original IHC image
		
		a_plot <-
			to_plot |>
			ggplot2::ggplot(aes(fill = bin, y = Pct, x = '')) + 
			geom_bar(position = "stack", stat = "identity") +
			scale_fill_manual(values = colour_mapping) +
			guides(fill = "none") +
			theme_bw() 
		gini <- DescTools::Gini(to_plot$Pct, na.rm = TRUE)
	
	}
	

}


plot_ihc_montage <- function(y, p, metadata, target) {

	## This function will create a montage of all HPA IHC images for a specific target, and a pre-specified list of healthy tissues, and write this to a pdf

	# target = "CEACAM5"; metadata = all_targets; tar_load(hpa_ihc_url_tb); y = hpa_ihc_url_tb	
	# Load IHC images for a specific Antibody
	# Annotate images: for tumour, use a coloured border or a symbol (L-H|W-S|%)
	# Tile into a montage (sort by pct brown)

	ensgid <- 
		metadata |>
		dplyr::filter(SYMBOL == target) |>
		dplyr::pull(ENSGID)

	tissue__ids <- c(p$core_tissue, "Bone marrow", "Pancreas", "Salivary gland")

	filepath_ls <- 
		y |>
		dplyr::filter(ENSGID %in% ensgid) |>
		dplyr::filter(snomedCode1 == "M-00100") |> # filter to all 'Normal tissue'
		dplyr::filter(tissueDescription2 %in% tissue__ids) |>
		dplyr::arrange(snomedCode1, snomedCode2, snomedCode3) |>
		dplyr::mutate(label = tissueDescription2)
	
	# --- Read the images directly from URLs ---
	all_img <- magick::image_read(path = filepath_ls$imageUrl) 
	#  --- Annotate each image with its label ---
	annotated_images <- 
		magick::image_annotate(
			all_img, 
			text = filepath_ls$label, 
            size = 300, 
			color = "#D6D6D6", 
			boxcolor = "#003EC1",
            gravity = "southwest"
		) |>
		magick::image_annotate(
			text = filepath_ls$ab_id, 
            size = 100, 
			color = "#D6D6D6", 
			boxcolor = "#B838B4",
            gravity = "northwest"
		) 		

	# --- Create a montage of the annotated images ---
	montage <- 
		magick::image_montage(
			annotated_images, 
			geometry = 'x200+10+10', 
			tile = '3x3', 
			bg = 'white', 
			shadow = FALSE
		)
	
	# Display the montage
	# print(montage)
	magick::image_write(montage, path = glue::glue("IMG/tma_healthy_{target}.pdf"), format = "pdf")

}

get_ihc_fet <- function(y, metadata, target) {

	## This function will apply Fishers Exact Test to contingency table of HPA IHC image annotations for a target 

	# target = "ST14"; metadata = all_targets; tar_load(hpa_ihc_url_tb); y = hpa_ihc_url_tb	
	# Load IHC images for a specific Antibody
	# Annotate images: for tumour, use a coloured border or a symbol (L-H|W-S|%)
	# Tile into a montage (sort by pct brown)

	ensgid <- 
		metadata |>
		dplyr::filter(SYMBOL == target) |>
		dplyr::pull(ENSGID)

	filepath_ls <- 
		y |>
		dplyr::filter(ENSGID %in% ensgid) |>
		dplyr::filter(snomedCode1 != "M-00100") |> # filter out all 'Normal tissue'
		dplyr::arrange(snomedCode1, snomedCode2, snomedCode3) |>
		dplyr::mutate(label = tissueDescription2) |>
		dplyr::mutate(staining = factor(staining, levels = c("Not detected", "Low", "Medium", "High"))) |>
		dplyr::mutate(intensity = factor(intensity, levels = c("Negative", "Weak", "Moderate", "Strong")))
	
	# --- Summarise annotation data ---
	# For tumour data, create contingency tables of staining x intensity per indication -> mosaic plot, suitable for Fishers Exact Test
	# TODO: compare global to local (indication-specific) table
	global_tb <- 
		filepath_ls |>
		# dplyr::group_by(tissueDescription1, tissueDescription2) |>
		dplyr::count(staining, intensity) |>
		tidyr::pivot_wider(names_from = intensity, values_from = n)


}


