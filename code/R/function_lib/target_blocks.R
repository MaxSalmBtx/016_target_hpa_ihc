
core_data_ls <- 
	list(
		tar_target(param_file, "target_params.json", format = "file"),
		tar_target(param_ls, jsonlite::read_json(param_file, simplifyVector = TRUE)),
		tar_target(references_ensgid, prepare_ensembl_ids(x = c(param_ls$references))),
		tar_target(current_targets_f, "../../input/Bicycle_Merged_Target_List_Reserved_targets_06_02_26_FINAL.xlsx", format = "file"),
		tar_target(current_targets, get_target_list(current_targets_f)),
		tar_target(all_targets, merge_targets(x = current_targets, y = references_ensgid))
	)

hpa_ihc_image_module <- 
	list(
		tar_target(hpa_overview_ls, get_hpa_overview()),
		tar_target(ensgids, get_ensgids(x = all_targets, y = hpa_overview_ls)),
		tar_target(
			name = ensgids_paths,
			command = get_image_path(ensgids),
			pattern = map(ensgids)
		),
		tar_target(
			name = ensgids_paths_tb,
			command = dplyr::bind_rows(ensgids_paths) 
		),
		tar_target(urls, as.list(ensgids_paths_tb$imageUrl)),
		tar_target(
			name = remote_data_files,
			command = as.character(urls),
			pattern = map(urls),
			format = "url" # watch for changes at the source (slow).
		),
		# tar_target(
			# name = image_paths,
			# command = get_ihc(urls),
			# pattern = map(urls),
			# cue = tar_cue(mode = "always") # for some reason, this is not always run first
		# ),
		tar_target(
			name = image_clean,
			command = clean_img_transparent(remote_data_files), 
			pattern = map(remote_data_files)
		),
		tar_target(ihc_f, get_img_dir(x = image_clean), format = "file"),
		tar_target(
			name = image_hist_ls, 
			command = 	
				colordistance::getLabHistList(
					ihc_f, 
					bins = 3, 
					color.space = "hsv",
					lower = c(0, 0, 0),
					upper = c(0.1, 0.1, 0.1), # Filter out black background
					ref.white = "D65",
					plotting = FALSE
				)
			),
		tar_target(img_dist, get_hist_dist(x = image_hist_ls, y = ensgids_paths)),
		tar_target(ihc_hm, get_heatmap(x = img_dist, y = ensgids_paths, p = param_ls, g = hpa_overview_ls, a = all_targets)), 
		tar_target(ihc_bp, get_stacked_barplot(x = image_hist_ls, y = ensgids_paths)),
		# Specific report for this output
		tar_target(report_ihc_healthy_t, class(ihc_hm)[1] == "tibble"),		
	    tarchetypes::tar_quarto(
			report_ihc_healthy,
			path = "report_ihc_healthy.qmd",
			quiet = FALSE,
			execute_params = list(dummy_param = report_ihc_healthy_t)
		)		
	)

