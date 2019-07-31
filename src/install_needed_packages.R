#!/usr/bin/env Rscript

#####################################
## name: install_needed_packages.R
## date: 2019-06-21
## version 0.1.0
## author: Dylan Morris <dhmorris@princeton.edu>
##
## installs needed packages for
## reproducing Australia influenza
## study
##
####################################

install_if_absent <- function(package_name){
    if (!suppressPackageStartupMessages(
             require(package_name, character.only=TRUE))){
      install.packages(pkgs=package_name,
                       repos="http://cloud.r-project.org")
  }
  else
      cat(sprintf("Package %s already installed\n", package_name))
}

needed_packages <- c(
    "rstan",
    "dplyr",
    "readr",
    "ggplot2",
    "tidybayes",
    "cowplot")

for (package in needed_packages)
    install_if_absent(package)
