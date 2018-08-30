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

plot_polydiv_faceted <- function(polyfile, divfile, scaffold_list, prefix, ref, plotfile, window_size, ymax) {
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
   ggplot(subset(polydiv_melted, Scaffold %in% scaffold_list),
      aes(x=Position, y=Value, group=Statistic, colour=Statistic)) +
      geom_line() +
      facet_wrap(~ Scaffold, ncol=2) +
      theme_bw() +
      ylim(0, ymax) +
      ggtitle(paste0("Polymorphism and Divergence in ", round(window_size/1000), "kb windows for ", prefix, " against ", ref));
   ggsave(plotfile, width=10.5, height=8.0, units="in");
}

plot_polydiv_separate <- function(polyfile, divfile, scaffold_list, prefix, ref, plotfile, window_size, ymax) {
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
   pdf(plotfile, width=10.5, height=8.0)
   for (scaffold in scaffold_list) {
      print(ggplot(subset(polydiv_melted, Scaffold == scaffold),
         aes(x=Position, y=Value, group=Statistic, colour=Statistic)) +
         geom_line() +
         theme_bw() +
         ylim(0, ymax) +
         ggtitle(paste0("Polymorphism and Divergence in ", round(window_size/1000), "kb windows on scaffold ", scaffold, " for ", prefix, " against ", ref)))
   }
   dev.off()
}

individual_plot_polydiv <- function(dir, prefix, ref, window_size, ymax, scaffold_list) {
   suffix <- paste0("_w", round(window_size/1000), "kb");
   tsvsuffix <- paste0(suffix, ".tsv");
   pdfsuffix <- paste0(suffix, ".pdf");
   polyfile <- paste0(dir, prefix, "_poly", tsvsuffix);
   divfile <- paste0(dir, prefix, "_div", tsvsuffix);
   plotfile <- paste0(dir, prefix, "_polydiv", pdfsuffix);
   if (length(scaffold_list) > 6) {
      plot_polydiv_separate(polyfile, divfile, scaffold_list, prefix, ref, plotfile, window_size, ymax)
   } else {
      plot_polydiv_faceted(polyfile, divfile, scaffold_list, prefix, ref, plotfile, window_size, ymax)
   }
}

dir <- options[1]
prefix <- options[2]
ref <- options[3]
window_size <- as.integer(options[4])
ymax <- as.numeric(options[5])
scaflist <- options[6]

if (scaflist == "allArms") {
   scaffold_list <- c("X", "2L", "2R", "3L", "3R", "4");
} else if (scaflist == "arms") {
   scaffold_list <- c("X", "2L", "2R", "3L", "3R");
} else if (scaflist == "chroms") {
   scaffold_list <- c("X", "2", "3");
} else if (scaflist == "allChroms") {
   scaffold_list <- c("X", "2", "3", "4");
} else {
   scaffold_list <- unlist(strsplit(scaflist, ","));
}

individual_plot_polydiv(dir, prefix, ref, window_size, ymax, scaffold_list)
