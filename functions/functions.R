# Helper Functions #

# * add_taxa_to_ranges
# * appearance_plot_save
# * arrange_pbdb_diversity_data
# * build_bins_df
# * count_bin_flapp
# * count_pbdb_taxa
# * count_pbdb_taxon_diversity
# * count_taxa_in_slice
# * diversity_plot_axes
# * diversity_plot_eot_label
# * diversity_plot_facet_taxon
# * diversity_plot_join
# * diversity_plot_save
# * diversity_plot_theme
# * filter_slice_taxa
# * pretty_num
# * process_pyrate_results
# * pyrate_plot_axes
# * pyrate_plot_curves
# * pyrate_plot_eotlabel
# * pyrate_plot_save
# * pyrate_plot_shifts
# * pyrate_plot_theme
# * read_pyrate_bayesfactors
# * read_pyrate_results
# * tidy_palaeogeog_polygons
# * write_for_sqs

add_taxa_to_ranges <-
  function(ranges_table = pbdb_ranges, class_list = classification) {
    # Match a vector of taxon names to taxa in the ranges table, create a vector
    # and join to the table so that taxa in the ranges have a higher taxon
    # assigned to them.
    #
    # Args:
    #   ranges_table: a table of taxon ranges with an "accepted_name" column
    #   class_list: a list of accepted names to match in the ranges table.
    #
    # Returns:
    #   A tibble with "taxon" column assigning the taxonomic group for each
    #   accepted name.
    taxon_col <- vector("character", length = nrow(ranges_table))
    for (i in seq_along(class_list)) {
      rows <-
        ranges_table$accepted_name %in% class_list[[i]]
      taxon_col[rows] <- names(class_list)[i]
    }
    ranges_table %>%
      add_column(taxon = taxon_col)
  }

appearance_add_silhouettes <- function() {
  ggdraw() +
  draw_plot(appearance_plot) +
  pmap(
    appearance_silhouettes,
    ~ draw_image(
      image_data(..1, size = 256)[[1]],
      x = ..2, y = ..3,
      height = 0.08, width = 0.08
    )
  )
}

appearance_plot_save <- function(plot_obj = last_plot()) {
  ggsave(
    paste0(dirs$figs, "fig_shelf_appearances.pdf"),
    plot = plot_obj,
    device = cairo_pdf,
    width = 80,
    height = 160,
    units = "mm"
  )
}

appearance_plot_theme <- function() {
  # A custom theme that compresses the legend.
  list(
    diversity_plot_theme(),
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "mm"),
      legend.spacing = unit(1, "mm"),
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 8)
    )
  )
}

arrange_app_data <- function(appearance_data) {
  # Arrange counts of first and last appearances into long form ready for
  # plotting. Convert the appearance curves to a facet with labels to plot
  # nicely.
  #
  # Args:
  #   appearance_data: (data.frame) table with bin ages and appearance data in
  #     'first_app' and 'last_app' columns.
  #
  # Returns:
  #   A long-form tibble that can be passed to e.g. gpplot. NB appearance counts
  #   are included in a column called 'Diversity' as this is used in other
  #   functions called for plotting.
  appearance_data %>%
    pivot_longer(
      c(first_app, last_app),
      names_to = "appearance_curve",
      # I use this to more easily reuse code from the diversity plot above –
      # these are counts of first and last appearances within that bin.
      values_to = "Diversity"
    ) %>%
    mutate(
      appearance_curve =
        factor(appearance_curve, labels = c("First", "Last"))
    )
}

arrange_map_occurrences <- function(occurrences) {
  # arrange PBDB occurrence data and assign taxon and latitude band for map
  # plotting
  occurrences %>%
    bind_rows() %>%
    filter(!is.na(paleolat)) %>%
    mutate(
      lat_band =
        case_when(
          paleolat > 23.5 | paleolat < -23.5 ~ "High latitude",
          paleolat <= 23.5 & paleolat >= -23.5 ~ "Tropics"
        ) %>%
        as_factor(),
      taxon =
        case_when(
          class  == "Anthozoa" ~ "Anthozoa",
          class  == "Bivalvia" ~ "Bivalvia",
          phylum == "Bryozoa" ~ "Bryozoa",
          phylum == "Echinodermata" ~ "Echinodermata",
          class  == "Gastropoda" ~ "Gastropoda"
        ),
      age_name =
        fct_relevel(age_name, pg_ages$unit_name)
    )
}

arrange_pbdb_diversity_data <-
  function(wide_diversity_data) {
    wide_diversity_data %>%
      pivot_longer(
        cols = -any_of(c("midpoint", "collections")),
        names_to = "Taxon",
        values_to = "Diversity",
        values_drop_na = TRUE
      ) %>%
      mutate(
        Group = factor(Taxon, levels = pbdb_diversity_labels)
      )
  }

build_bins_df <- function (start, end, interval = NULL, n = NULL) {
  # Construct a data.frame of bin start and end dates from a time range and
  # interval length or number required.
  #
  # Args:
  #   start, end: the start and end points/time/depth of the series
  #   interval:   length of each time interval/bin
  #   n:          number of bins to split the series into
  #
  #   NB a warning will be given if both the interval length and number of bins
  #   are specified
  #
  # Returns:
  #   A data.frame with bins in rows and start and end dates and mid points in
  #   each column.
  require(tibble)
  # Warn if both interval length and number of intervals are specified together.
  if (!is.null(interval) & !is.null(n)) {
    warning("Interval length and number of bins should not be used together.")
  }

  # convert number of bins to bin length
  if (is.null(interval) & !is.null(n)) {
    interval <- abs((start - end) / n)
  }

  # correct the interval sign if given in 'before present' style
  if (start > end) {
    interval <- -abs(interval)
  }

  # get a vector of bin boundary times
  bin_bound_seq <- seq(start, end, interval)

  # construct data.frame
  # offset start and end times by 1 to produce bins of required length
  bins <- tibble(
    start = bin_bound_seq[-length(bin_bound_seq)],
    end   = bin_bound_seq[-1]
  )

  # add bin midpoints
  bins$midpoint <- rowMeans(bins)
  bins
}

count_bin_flapp <- function(ranges_table, bin_data, name = "genus") {
  # Counts the number of taxa that have their first and last occurrences within
  # a set of bins based on a set of taxon first and last appearance ages.
  #
  # Args:
  #   ranges_table: (data.frame) first and last appearance dates for a set of
  #     taxa.
  #   bin_data: (data.frame) bin boundaries with columns start_Ma and end_Ma.
  #   name: (character) the name of the column that includes taxon names.
  #
  # Returns:
  #   A tibble of per-bin data with counts of first and last appearances.
  # must order bins with oldest first to get counting direction correct
  bin_data <-
    bin_data %>% 
    arrange(desc(start_Ma))
  bin_taxa <-
    bin_data %>%
    select(start_Ma, end_Ma) %>%
    pmap(
      function(start_Ma, end_Ma) {
        ranges_table %>%
          filter(max_ma >= end_Ma & min_ma < start_Ma) %>%
          pull(name)
      }
    )
  appearances <- vector("list", nrow(bin_data) - 2)
  for(i in seq(2, nrow(bin_data) - 1)) {
    curr_bin <- bin_data[i, ]
    taxa_curr <- bin_taxa[[i]]
    taxa_prev <- bin_taxa[[i - 1]]
    taxa_next <- bin_taxa[[i + 1]]
    appearances[[i - 1]] <-
      curr_bin %>%
      bind_cols(
        tibble(
          first_app = length(which(!taxa_curr %in% taxa_prev)),
          last_app = length(which(!taxa_curr %in% taxa_next))
        )
      )
  }
  bind_rows(appearances) %>%
    mutate(
      midpoint = (start_Ma + end_Ma) / 2
    )
}

count_lbf_in_slice <- function(bin_margin, occurrence_table) {
  # Counts the number of taxa between defined margins at the specified rank
  # from PBDB occurrence data.
  #
  # Args: occurrence_table: a table of PBDB data with accepted_rank,
  # accepted_name, min_ma, and max_ma columns. rank: taxonomic level at which
  # to count unique names.
  #
  # Returns: The number of unique taxonomic names.
  occurrence_table[
    occurrence_table$min_ma <= bin_margin[["midpoint"]] &
      occurrence_table$max_ma >= bin_margin[["midpoint"]],
    ] %>%
      nrow()
}

count_pbdb_taxa <- function(occurrence_table, rank = "species") {
  # Counts the number of taxa between defined margins at the specified rank
  # from PBDB occurrence data.
  #
  # Args: occurrence_table: a table of PBDB data with accepted_rank,
  # accepted_name, min_ma, and max_ma columns. rank: taxonomic level at which
  # to count unique names.
  #
  # Returns: The number of unique taxonomic names.
  length(
    unique(
      occurrence_table[
        occurrence_table$accepted_rank == rank,
        "accepted_name"
      ]
    )
  )
}

count_pbdb_taxon_diversity <- function() {
  # Counts the slice diversity at species rank (by default) for taxa indentified
  # in taxon_matching.
  pblapply(
    taxon_matching,
    function(tax) {
      apply(
        bins, 1, count_taxa_in_slice,
        pbdb_ranges %>%
          filter(taxon == tax$taxon_name)
      )
    }
  )
}

count_taxa_in_slice <-
  function(bin_margin, ranges_table = pbdb_ranges) {
  # Count the number of taxa within a bin from PBDB occurrence data.
  #
  # Args:
  #   bin_margin: a row/vector of slice timepoint indicated by "midpoint".
  #   ranges_table: a data.frame or tibble of PBDB data with taxon ranges in
  #     min_ma, and max_ma columns.
  #
  # Returns:
  #   The number of taxa in the specified bin.
  ranges_table %>%
    filter(
      min_ma <= bin_margin["midpoint"] &
        max_ma >= bin_margin["midpoint"]
    ) %>%
    nrow()
}

create_taxon_vector <-
  function (taxon) {
    # Create a vector of row numbers from the ranges table with assigned taxa
    # that can be joined later.
    #
    # Args:
    #   taxon: the taxon from taxon_matching
    #
    # Returns:
    #   A vector of taxon names matching the filtered occurrences.
    pbdb_occurrences %>%
      filter(taxon$sp_occurrences) %>%
      pull(accepted_name) %>%
      unique()
  }

diversity_plot_axes <- function() {
  # Add reversed x-axis (decreasing towards the present, a custom colour scale,
  # geological time scale at the bottom, and x-axis label.
  list(
    scale_x_reverse(limits = c(66, 23)),
    expand_limits(y = 0),
    xlab("Age (Ma)"),
    scale_discrete_manual(
      values =
        c("#2F9599", "#EC2049", "#4378CB", "black", "#F26B38", "#A7226E"),
      aesthetics = c("colour")
    ),
    coord_geo(
      expand = TRUE,
      pos = as.list(rep("bottom", 2)),
      dat = list("stages", "epochs"),
      height = list(unit(1.5, "mm"), unit(4, "mm")),
      size = list(3, 4),
      skip = c("Late Cretaceous", "Miocene", stages$name),
      center_end_labels = TRUE
    )
  )
}

diversity_plot_label <- function(events = event_ages, facet_name = "Diversity") {
  # Add rectangle and text label indication the temporal range of the
  # Eocene–Oligocene transition (34–35 Ma).
  #
  # Args:
  #   events: (data.frame) a table of event start and end ages with columns:
  #     * event: name of the event (label string to plot)
  #     * xmax, xmin: start and end ages (in that order) of the event.
  #   facet_name: (character) the name of the facet on which to show the label.
  list(
    geom_rect(
      data = events,
      mapping =
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        alpha = 1.0, fill = "#F7DB4F", inherit.aes = FALSE
    ),
    geom_text(
      data =
        . %>%
          filter(Group == facet_name) %>%
          summarise(yval = rep(max(Diversity), nrow(events))) %>%
          add_column(
            event_label = events %>% pull(event),
            Group = factor(facet_name),
            midpoint = events %>%
              transmute((xmin + xmax) / 2) %>%
              pull(1),
            hjust = events %>% pull(hjust),
            vjust = events %>% pull(vjust)
          ),
      mapping =
        aes(x = midpoint, y = yval, label = event_label, hjust =
          hjust, vjust = vjust),
      size = 2, inherit.aes = FALSE
    )
  )
}

diversity_plot_facet_taxon <- function() {
  # Organise facets for curves splitting by taxon.
  list(
    facet_wrap(
      vars(Group), ncol = 1, scales = "free_y",
      strip.position = "right"
    )
  )
}

diversity_plot_join <- function() {
  p1 <-
    plot_grid(
      plot_raw_diversity, plot_sqs_diversity,
      labels = "AUTO"
    )
  ggdraw() +
    draw_plot(p1) +
    apply(diversity_annotations, 1, function(taxon) {
      draw_image(
        rphylopic::image_data(taxon[["phylopic"]], size = 256)[[1]],
        x = as.numeric(taxon[["x"]]), y = as.numeric(taxon[["y"]]),
        height = 0.05, width = 0.05
      )
    })
}

diversity_plot_save <- function(plot_obj = last_plot()) {
  ggsave(
    paste0(dirs$figs, "fig_pbdb_diversity.pdf"),
    plot = plot_obj,
    device = cairo_pdf,
    width = 166,
    height = 150,
    units = "mm"
  )
}

diversity_plot_theme <- function() {
  # A custom theme that removes the legend
  list(
    theme_light(),
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8)
    )
  )
}

filter_slice_taxa <- function(occurrence_table, margin) {
  # Filters unique taxon names at the specified margin.
  #
  # Args:
  #   occurrence_table: a table of PBDB data with accepted_rank,
  #                     accepted_name, min_ma, and max_ma columns.
  #   margins: a two element vector
  #            of maximum and minimum ages.
  #
  # Returns: Row-occurrences at the defined slice margin.
  occurrence_table[
    occurrence_table$min_ma <= margin[["midpoint"]] &
      occurrence_table$max_ma >= margin[["midpoint"]],
  ]
}

lbfdiversity_add_silhouette <- function() {
  ggdraw() +
    draw_plot(lbf_diversity_plot) +
    draw_image(
      image_data(lbf_phylopic$phylopic_uuid, size = 256)[[1]],
      x = 0.15, y = 0.80,
      height = 0.15, width = 0.15
    )
}

lbfdiversity_plot_save <- function(plot_obj = last_plot()) {
  ggsave(
    paste0(dirs$figs, "fig_lbf_diversity.pdf"),
    plot = plot_obj, 
    device = cairo_pdf,
    width = 80,
    height = 60,
    units = "mm"
  )
}

map_plot_palaeogeography <- function() {
  list(
    geom_polygon(
      data = palaeogeog_data,
      aes(
        x     = long,
        y     = lat,
        fill  = layer,
        group = polygon_id
      ),
      colour = NA
    )
  )
}

map_plot_occurrences <- function() {
  list(
    geom_point(
      data = pg_occs,
      aes(x = paleolng, y = paleolat, colour = lat_band),
      size = 0.25, na.rm = TRUE
    ),
    scale_discrete_manual(
      values =
        c(
          "High latitude" = "#003DFF",
          "Tropics"       = "#FF004D"
        ),
        aesthetics = c("colour"),
        name = "Occurrence latitude"
    ),
    guides(colour =
      guide_legend(
        nrow = 2, byrow = TRUE,
        override.aes = list(size = 3)
      )
    )
  )
}

map_plot_prettify <- function() {
  # Adds ggplot layers to make palaeogeographical maps prettier:
  #
  # * Mollweide projection
  # * Equator, tropics lines, outline border
  # * Use ggthemes::theme_map to tidy the plot area
  # * legend at the bottom
  # * use palaeogeographical colour scheme
  # * two row legend
  list(
    coord_map("mollweide"),
    geom_segment(
      data = tibble(
        x = rep(-180, 3),
        xend = rep(180, 3),
        y = c(-23.5, 0, 23.5)
      ),
      aes(x = x, xend = xend, y = y, yend = y),
      colour = "grey70", size = 0.3
    ),
    geom_rect(
        data = data.frame(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        color = 1, fill = NA, size = 0.3
    ),
    theme_map(),
    theme(
      legend.position  = c(1.00, 0.05),
      legend.justification = c("right", "bottom"),
      legend.box = "horizontal",
      panel.spacing    = unit(1, "mm"),
      strip.background = element_blank()
    ),
    scale_discrete_manual(
      values =
        c(
          "Shallow marine" = "#45D8FF",
          "Land"           = "#FFD23A",
          "Mountain"       = "#FF8D51",
          "Ice cap"        = "#DAD3FF"
        ),
      aesthetics = c("fill"),
      name = "Palaeogeography"
    ),
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  )
}

map_save_base <- function(plot_obj = last_plot()) {
  ggsave(
    file = paste0(dirs$figs, "fig_base_maps.pdf"),
    plot = plot_obj,
    device = cairo_pdf,
    width = 166,
    height = 220,
    units = "mm"
  )
}

map_save_occurrences <- function(plot_obj = last_plot()) {
  ggsave(
    file = paste0(dirs$figs, "fig_occurrence_maps.pdf"),
    plot = plot_obj,
    device = cairo_pdf,
    width = 166,
    height = 220,
    units = "mm"
  )
}

pretty_num <- function(number) {
  # Prints a large number with thousands separated by spaces. If the number is
  # greater than 10^6 then it will be printed in scientific notation.
  #
  # Args:
  #   number: number to print
  #
  # Returns:
  #   A formatted pretty number with thousands separator.
  if (number < 1e6) {
    as_scientific <- FALSE
  } else {
    as_scientific <- TRUE
  }
  prettyNum(number, scientific = as_scientific, big.mark = " ")
}

process_pyrate_results <- function() {
  # Reads the results of the PyRate analyses for each taxon and arranges the
  # data into a long-form table ready for plotting with ggplot2. Results are
  # output to CSV files.
  data_rtt <-
    read_lines(paste0(dirs$data, "t260-200RTT_plots.r"))
  # read in relevant lines for each subplot
  source(textConnection(data_rtt[7:10]))
  sp_rates <- tibble(
    time = time,
    rate = rate,
    minHPD = minHPD,
    maxHPD = maxHPD,
    plt_curve = rep("sp_rates", length(time))
  )
  source(textConnection(data_rtt[14:15]))
  sp_shifts <- tibble(
    time = mids,
    rate = counts,
    plt_curve = rep("sp_shifts", length(mids))
  )
  source(textConnection(data_rtt[22:25]))
  ex_rates <- tibble(
    time = time,
    rate = rate,
    minHPD = minHPD,
    maxHPD = maxHPD,
    plt_curve = rep("ex_rates", length(time))
  )
  source(textConnection(data_rtt[29:30]))
  ex_shifts <- tibble(
    time = mids,
    rate = counts,
    plt_curve = rep("ex_shifts", length(mids))
  )
  source(textConnection(data_rtt[35:38]))
  net_rates <- tibble(
    time = time,
    rate = net_rate,
    minHPD = net_minHPD,
    maxHPD = net_maxHPD,
    plt_curve = rep("net_rates", length(time))
  )
  source(textConnection(data_rtt[45:48]))
  longevity <- tibble(
    time = time,
    rate = rate,
    minHPD = minHPD,
    maxHPD = maxHPD,
    plt_curve = rep("longevity", length(time))
  )
  source(textConnection(data_rtt[17:18]))
  bf <- tibble(
    low = bf2,
    high = bf6,
  )
  # join the results into a single long tibble
  data_file <-
    sp_rates %>%
      bind_rows(sp_shifts) %>%
      bind_rows(ex_rates) %>%
      bind_rows(ex_shifts) %>%
      bind_rows(net_rates) %>%
      bind_rows(longevity)
  # write files
  write_csv(
    data_file,
    paste0(dirs$pyrate, "ichthyosaur-pyrate-results.csv")
  )
  write_csv(
    bf,
    paste0(dirs$pyrate, "ichthyosaur-pyrate-bf.csv")
  )
}

pyrate_facet_taxon_curve <- function() {
  # Organise grid of facets by taxon (columns) and curves plotted (rows)
  # To be added to a ggplot object.
  facet_grid(
    rows = vars(plt_curve),# cols = vars(taxon),
    scales = "free_y", shrink = FALSE, switch = "y",
    labeller = label_wrap_gen(width = 20)
  )
}

pyrate_plot_curves <- function() {
  # Add line and uncertainty ribbons to a ggplot, for diversification rates and
  # longevity
  list(
    geom_ribbon(
      data = . %>% filter(!str_detect(plt_curve, "shifts")),
      alpha = 0.2, colour = NA
    ),
    geom_line(
      data = . %>% filter(!str_detect(plt_curve, "shifts"))
    )
  )
}

pyrate_plot_axes <- function() {
  # Add reversed x-axis (decreasing towards the present), a custom colour scale,
  # geological time scale at the bottom, and x-axis labels.
  list(
    scale_x_reverse(limits = c(250, 214)),
    scale_discrete_manual(
      values = plot_colours,
      aesthetics = c("colour", "fill")
    ),
    coord_geo(
      expand = TRUE,
      pos = as.list(rep("bottom", 2)),
      dat = list("stages", "epochs"),
      height = list(unit(1.5, "mm"), unit(4.0, "mm")),
      size = list(3, 3),#unit(1, "mm"), unit(2.5, "mm")),
      skip = c("Late Cretaceous", "Miocene", stages$name),
      center_end_labels = TRUE
    ),
    labs(x = "Age (Ma)", y = NULL)
  )
}

pyrate_plot_eotlabel <- function() {
  # Add rectangle and text label indicating the temporal range of the
  # Eocene–Oligocene transition (34–33.5 Ma).
  # To be added to a ggplot object.
  list(
    geom_rect(

      mapping =
        aes(xmin = 33.5, xmax = 34, ymin = -Inf, ymax = Inf),
      alpha = 1.0, fill = "#F7DB4F", inherit.aes = FALSE,
    ),
    geom_text(
      data = . %>%
        filter(plt_curve == "Speciation rate") %>%
        summarise(yval = max(maxHPD), plt_curve = factor("Speciation rate")),
      mapping = aes(x = 33.75, y = yval, label = "EOT"),
      size = 2, inherit.aes = FALSE
    )
  )
}

#diversity_plot_label <- function(events = event_ages, facet_name = "Diversity") {
#  # Add rectangle and text label indication the temporal range of the
#  # Eocene–Oligocene transition (34–35 Ma).
#  #
#  # Args:
#  #   events: (data.frame) a table of event start and end ages with columns:
#  #     * event: name of the event (label string to plot)
#  #     * xmax, xmin: start and end ages (in that order) of the event.
#  #   facet_name: (character) the name of the facet on which to show the label.
#  list(
#    geom_rect(
#      data = events,
#      mapping =
#        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
#        alpha = 1.0, fill = "#F7DB4F", inherit.aes = FALSE
#    ),
#    geom_text(
#      data =
#        . %>%
#          filter(Group == facet_name) %>%
#          summarise(yval = rep(max(Diversity), nrow(events))) %>%
#          add_column(
#            event_label = events %>% pull(event),
#            Group = factor(facet_name),
#            midpoint = events %>%
#              transmute((xmin + xmax) / 2) %>%
#              pull(1),
#            hjust = events %>% pull(hjust)
#          ),
#      mapping =
#        aes(x = midpoint, y = yval, label = event_label, hjust = hjust),
#      size = 2, inherit.aes = FALSE
#    )
#  )
#}

pyrate_plot_save <- function(plot_obj = last_plot()) {
  # Save the plot as a PDF file
  ggsave(
    paste0(dirs$figs, "fig_pyrate_results.pdf"),
    plot = plot_obj,
    device = cairo_pdf,
    width = 226,
    height = 166,
    units = "mm"
  )
}

pyrate_plot_shifts <- function() {
  # Add frequency shifts as segments with horizontal lines to indicate Bayes
  # Factor levels
  list(
    geom_segment(
      data = . %>%
        filter(str_detect(plt_curve, "shifts")),
      size = 0.5
    ),
    geom_hline(
      mapping = aes(yintercept = value),
      colour = "grey70", inherit.aes = FALSE
    )#,
    # geom_text(
    #   # data = . %>% filter(taxon == "Anthozoa"),
    #   mapping = aes(x = 66, y = value, label = line_label),
    #   hjust = 0, vjust = -0.2, size = 2,
    #   inherit.aes = FALSE
    # )
  )
}

pyrate_plot_theme <- function() {
  # A custom theme that removes minor horizontal gird lines and all vertical
  # grid lines. Also hides the legend.
  list(
    theme_light(),
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text.y = element_text(size = 8)
    )
  )
}

read_palaeogeog_data <-
  function(
    data_dirs,
    subset_names = names(palaeogeog_colours)
  ) {
    # Read palaeogeographical reconstructions in OGR-RMT file format from a
    # directory. The data are processed to assign subsets and join multiple files.
    data_dirs %>% 
      map2(
        subset_names,
        ~ tidy_palaeogeog_polygons(polygon_dir = .x, polygon_desc = .y)
      ) %>%
      bind_rows() %>%
      mutate(
        polygon_id = str_c(id, group, age_name, layer) %>% as_factor(),
        age_name = factor(age_name, levels = pg_ages$unit_name)
      ) %>%
      sort_polugonid_order(layer_names = subset_names)
}

read_pyrate_bayesfactors <- function() {
  # Read and combine bayes factor CSV files output from Pyrate.
  read_csv(paste0(dirs$pyrate, "ichthyosaur-pyrate-bf.csv")) %>%
  # Pivot to a long form table ready for plotting
  # pivot_longer(cols = -taxon) %>%
  # Repeat each BF value ready to assign curve names
  slice(rep(1:n(), each = 2)) %>%
  # Add curve names to plot in sp_shifts and ex_shifts facet rows
  add_column(
    plt_curve =
      rep(
        c("sp_shifts", "ex_shifts"),
        nrow(.) / 2
      ),
    line_label =
      rep(
        c("BF = 2", "BF = 6"),
        nrow(.) / 2
      )
  )
}

read_pyrate_results <- function() {
  # Read in main results of the Pyrate analyses from CSV files, combine with
  # pyrate_bayesfactors, and format the curve labels for plotting in ggplot2.
  read_csv(paste0(dirs$pyrate, "ichthyosaur-pyrate-results.csv")) %>%
  bind_rows(pyrate_bayesfactors) %>%
  mutate(
    # Convert the plot names columns to factors with labels for easy facet
    # plotting
    plt_curve =
      factor(
        plt_curve,
        levels = names(plot_labels),
        labels = plot_labels
      ),
    # Convert time to positive format
    time = abs(time)
  )
}

sort_polugonid_order <-
  function(
    polygon_data = palaeogeog_data,
    layer_names = palaeogeog_colours
    )
  {
    # Order the levels of the polygon_id column by the layers represented.
    #
    # Args:
    #   polygon_data: (data.frame) square data describing the polygons with a
    #     polygon_id _factor_ column built from layer names.
    #   layer_names: (character) names used to build the polygon_id in the order
    #     of plotting (first to last = bottom to top).
    #
    # Returns:
    #   A tibble with the levels of polygon_id column ordered according to
    #   layer_names.
    id_level_order <-
      layer_names %>%
        map(
          function(layr) {
            polygon_data$polygon_id %>%
              levels() %>%
              str_which(layr)
          }
        ) %>%
        unlist()
    polygon_data %>%
      mutate(
        polygon_id = fct_relevel(polygon_id, levels(polygon_id)[id_level_order])
      )
  }

tidy_palaeogeog_polygons <-
  function(
    polygon_dir,
    polygon_desc,
    strat_data = pg_ages
  ) {
    # Read and tidy sequential OGR-GMT palaeogeographic map polygon outlines.
    #
    # 'Tidy' data in this sense means a tibble (data.frame) with vertices for
    # the polygons contained in rows. age_Ma, age_name, and layer name columns
    # are added for plotting in facets and subsetting.
    #
    # The steps are:
    #
    # 1. list files
    # 2. extract reconstruction age from filename
    # 3. match to and assign chronostratigraphic names to times
    # 4. read and tidy data
    #
    # Args:
    #   polygon_dir: (character) the directory in which to list and read OGR
    #     data
    #   polygon_desc: (character) description of the polygons that will be used
    #     for the 'layer' column
    #   strat_data: (data.frame) data of charonostratigraphic bin ranges with
    #     columns named 'midpoint' and 'unit_name'
    #
    # Returns:
    #   A tibble of multiple polygons that can be plotted with, e.g., geom_map.
    poly_files <-
      list.files(polygon_dir, full.names = TRUE)
    poly_times <-
      poly_files %>%
      str_extract("[:digit:]+\\.[:digit:]+") %>%
      as.double()
    names(poly_times) <-
      pg_ages %>%
        slice(which(floor(pg_ages$midpoint) %in% floor(poly_times))) %>%
        pluck("unit_name")
    poly_data <-
      poly_files %>%
        map(readOGR) %>%
        map(tidy) %>%
        map2(poly_times, ~ add_column(.x, age_Ma = .y)) %>%
        map2(names(poly_times), ~ add_column(.x, age_name = .y)) %>%
        bind_rows() %>%
        add_column(layer = polygon_desc)
  }

write_for_sqs <- function(data_for_sqs, filepath = "output/data_for_sqs.csv") {
  # Write data formatted for SQS perl script to a CSV file. This file should
  # have no row names, NAs indicated with blank cells ("") and no quotes.
  #
  # Args:
  #   data_for_sqs: (data.frame) formatted for SQS with following columns:
  #                 - collection_no (leftmost column)
  #                 - occurrence.genus_name/species_name depending on taxonomic
  #                   level
  #                 - max_ma, min_ma or interval
  #                 - collection.reference_no
  #
  # Returns:
  #   Writes a CSV file with the relevant columns.
  write.csv(
    data_for_sqs,
    file = filepath,
    na = "",
    quote = FALSE,
    row.names = FALSE
  )
}
