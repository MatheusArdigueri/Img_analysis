# Shared functions. Run either numbered script in RStudio, after editing settings.
# Install ONCE, then restart R:
# install.packages('remotes')
# remotes::install_local('C:/your/path/ppe_0.2.3 (1).tar.gz',
#                        dependencies = NA, upgrade = 'never')
# install.packages('sf') # optional, for plot polygon summaries
library(terra)
library(ppe)
stopifnot(as.character(packageVersion('ppe')) == '0.2.3')
ppe_options(sample_n = 30000, threads = 2, seed = 42)

save_tif <- function(x, path) {
  writeRaster(x, path, overwrite = TRUE, datatype = 'FLT4S',
              gdal = c('COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=IF_SAFER'))
  rast(path)
}

# Band stacking: explicit order, single-band inputs, identical grids.
# Alignment is optional because resampling cannot correct poor coregistration.
build_stack <- function(paths, out, divisor = 1, align = FALSE) {
  stopifnot(length(paths) > 1, all(file.exists(paths)),
            is.finite(divisor), divisor > 0)
  layers <- lapply(paths, rast)
  if (any(vapply(layers, terra::nlyr, numeric(1)) != 1))    
    stop('Each input must contain exactly one reflectance band.')
  template <- layers[[1]]
  if (!nzchar(crs(template))) stop('Input has no CRS. Fix georeferencing first.')
  for (i in seq_along(layers)) {
    if (!compareGeom(template, layers[[i]], stopOnError = FALSE)) {
      if (!align) stop('Grid mismatch in ', paths[i],
                      '. Inspect alignment, or explicitly set align = TRUE.')
      layers[[i]] <- if (!same.crs(layers[[i]], template)) {
        project(layers[[i]], template, method = 'bilinear')
      } else resample(layers[[i]], template, method = 'bilinear')
    }
  }
  # Combine bands into one SpatRaster
  x <- layers[[1]]
  
  for (i in seq.int(2, length(layers))) {
    x <- c(x, layers[[i]])
  }
  
  # Convert stored values to reflectance
  x <- x / divisor
  names(x) <- names(paths)
  # Input nodata must already be correctly encoded. Never blanket-mask zeros:
  # a zero in one band can be a valid reflectance observation.
  save_tif(x, out)
}

# Same package ensemble for both sensor types. No XGBoost installation needed.
# All available indices are computed, but only selected layers segment the image.
run_ppe <- function(imp, full_stack, output_dir, plant_threshold = 0.4,
                    plots_file = NULL, plot_id = 'plot') {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  catalog <- ppe_list_vis(imp)
  write.csv(catalog, file.path(output_dir, 'index_formulas_and_availability.csv'),
            row.names = FALSE)
  available <- catalog$vi[catalog$available]
  vi <- ppe_compute_vis(imp, available,
                        outfile = file.path(output_dir, 'all_indices.tif'))
  spread <- ppe_spread(imp, include_vis = TRUE)
  write.csv(data.frame(layer = names(spread$spreads),
                       spread = unname(spread$spreads),
                       spread_class = unname(spread$spread_class)),
            file.path(output_dir, 'sample_spread.csv'), row.names = FALSE)
  write.csv(spread$correlation, file.path(output_dir, 'sample_correlation.csv'))
  selection <- ppe_select(imp, spread = spread, vi = TRUE, bands = TRUE,
                          top_n = 6, max_correlation = 0.7)
  # Supply only selected indices: avoids duplicated/unneeded segmentation layers.
  selected_vis <- intersect(selection$chosen, available)
  selected_stack <- if (length(selected_vis)) {
    ppe_compute_vis(imp, selected_vis,
                   outfile = file.path(output_dir, 'segmentation_indices.tif'))
  } else NULL
  seg <- ppe_segment(imp, selection, vi_stack = selected_stack,
                     strict_bimodality = FALSE)
  votes <- ppe_vote_fractions(seg)
  ppe_vote_threshold_guide(votes) # inspect before choosing a threshold
  mask <- ppe_create_mask(votes, plant_thresh = plant_threshold,
                          import = imp, ambiguous = 'na')
  # Conservative handling: uncertain pixels remain NA, rather than forced to soil.
  # Package mask: 1=plant, 0=soil, NA=uncertain/nodata.
  save_tif(mask$mask, file.path(output_dir, 'plant_soil_labels.tif'))
  save_tif(votes$fraction, file.path(output_dir, 'plant_vote_fraction.tif'))
  # Retain ORIGINAL reflectance layers, including EVERY HSI band.
  # The condition is applied to the whole stack; no soil is filled with zero.
  plant_path <- file.path(output_dir, 'plant_only_reflectance.tif')
  soil_path <- file.path(output_dir, 'soil_only_reflectance.tif')
  plant <- save_tif(ifel(mask$mask == 1, full_stack, NA), plant_path)
  save_tif(ifel(mask$mask == 0, full_stack, NA), soil_path)
  plant_vis <- ppe_apply_mask(imp, mask, vi_names = available,
                 outfile = file.path(output_dir, 'plant_only_all_indices.tif'))
  # Image-wide descriptive statistics (not classification accuracy).
  tab <- global(plant_vis$raster, c('mean', 'sd', 'min', 'max'), na.rm = TRUE)
  tab$index <- rownames(tab)
  write.csv(tab, file.path(output_dir, 'plant_index_summary.csv'), row.names = FALSE)
  write.csv(as.data.frame(freq(mask$mask)),
            file.path(output_dir, 'classified_pixel_counts.csv'), row.names = FALSE)
  qc <- ppe_qc(imp, votes = votes)
  write.csv(qc$summary, file.path(output_dir, 'brightness_vote_QC.csv'), row.names = FALSE)
  # Optional: all SIX package plot statistics + canopy cover.
  if (!is.null(plots_file)) {
    if (!file.exists(plots_file)) stop('Plot polygon file does not exist.')
    plot_summary <- ppe_summarize_plots(plots_file, mask = mask, vi_stack = vi,
      id_field = plot_id, vi_scope = 'plant',
      stats = c('mean', 'median', 'max', 'sd', 'cv', 'min'),
      outfile = file.path(output_dir, 'plot_metrics.csv'))
  }
  # Preserve all package diagnostics, including segmentation and selection.
  saveRDS(list(spread = spread, selection = selection, segmentation = seg,
               votes = votes, mask = mask, qc = qc),
          file.path(output_dir, 'package_diagnostics.rds'))
  capture.output(sessionInfo(), file = file.path(output_dir, 'sessionInfo.txt'))
  ppe_plot_mask(imp, mask)
  message('Plant-only GeoTIFF saved to: ', normalizePath(plant_path))
  invisible(list(import = imp, mask = mask, plant = plant, indices = plant_vis,
                 votes = votes, diagnostics = qc, plant_tiff = plant_path))
}
