#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, dependencies=TRUE, type="source", repos="http://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}

check_package("reshape2")
check_package("ggplot2")

plot_polydiv <- function(polyfile, divfile, scaffold_list, prefix, ref, plotfile, window_size, ymax) {
   poly <- read.table(polyfile,
      colClasses=c("character", "integer", "numeric"),
      col.names=c("Scaffold", "Position", "Polymorphism"));
   div <- read.table(divfile,
      colClasses=c("character", "integer", "numeric"),
      col.names=c("Scaffold", "Position", "Divergence"));
   polydiv <- poly;
   polydiv$Divergence <- div$Divergence;
   polydiv_melted <- melt(polydiv,
      id.vars=c("Scaffold", "Position"),
      variable.name="Statistic",
      value.name="Value");
   sampleid <- strsplit(prefix, "_", fixed=TRUE)[[1]];
   ggplot(subset(polydiv_melted, Scaffold %in% scaffold_list),
      aes(x=Position, y=Value, group=Statistic, colour=Statistic)) +
      geom_line() +
      facet_wrap(~ Scaffold, ncol=2) +
      theme_bw() +
      ylim(0, ymax) +
      ggtitle(paste0("Polymorphism and Divergence to ", ref, " in ", round(window_size/1000), "kb windows for ", sampleid));
   ggsave(plotfile, width=10.5, height=8.0, units="in");
}
drosophila_plot_polydiv <- function(dir, prefix, ref, polystem, divstem, window_size, ymax) {
   suffix <- paste0("_w", round(window_size/1000), "kb");
   tsvsuffix <- paste0(suffix, ".tsv");
   pdfsuffix <- paste0(suffix, ".pdf");
   polyfile <- paste0(dir, prefix, "_", polystem, tsvsuffix);
   divfile <- paste0(dir, prefix, "_", divstem, tsvsuffix);
   scaffold_list <- c("X", "2L", "2R", "3L", "3R", "4");
   plotfile <- paste0(dir, prefix, "_DxyPi", pdfsuffix);
   plot_polydiv(polyfile, divfile, scaffold_list, prefix, ref, plotfile, window_size, ymax)
}

dir <- options[1]
prefix <- options[2]
ref <- options[3]
polystem <- options[4]
divstem <- options[5]
window_size <- as.integer(options[6])
ymax <- as.numeric(options[7])

drosophila_plot_polydiv(dir, prefix, ref, polystem, divstem, window_size, ymax)
