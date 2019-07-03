#####################################
# name: Makefile
# date: 2019-06-21
# version: 0.1.0
# author: Dylan Morris <dhmorris@princeton.edu>
#
# Makefile to generate analyses
# for Lam et al 2019 study of
# Australian influenza epidemiology
#####################################


#####################################
# Expected bash settings
#
# Check these vs your local
# machine setup if you are having
# difficulty reproducing the
# analysis
#####################################

CXX := g++-8
CXXFLAGS := -O3 -std=c++11
MKDIR := @mkdir -p
RM := rm -rf

R_OPTIONS = --vanilla
R_COMMAND := Rscript $(R_OPTIONS)

#################################
# directory structure
##################################
SRC := src
DAT := dat
OUT := out

MCMC_CHAINS := $(OUT)/mcmc_chains
RAW_DATA_DIR := $(DAT)/raw
CLEANED_DATA_DIR := $(DAT)/cleaned
PLOT_PATH = $(OUT)/figures

FIGEXT := png
SRCEXT := py
CHAINS_SUFFIX = _chains.Rds

#####################################
# Installation / dependencies
#
# Rules for prepping analysis
#####################################

.PHONY: depend

depend:
	$(R_COMMAND) $(SRC)/install_needed_packages.R


#####################################
# ANALYSIS LOCATIONS
#
# where to find scripts and data for
# analyses
#####################################

#######
#data
#######

RAW_CLIMATE_DATA = $(RAW_DATA_DIR)/mean_fortnightly_climate_30years.csv

CLEAN_EPI_DATA = $(CLEANED_DATA_DIR)/epi_table.csv
CLEAN_STAN_DATA = $(CLEANED_DATA_DIR)/clean_stan_data.csv

################################
# for Bayesian analyses in Stan
################################
MODEL_FITTING_SCRIPT = $(SRC)/stan_fit.R
INCIDENCE_MODEL_NAME = multilevel_incidence_model
MODEL_DATA = $(CLEAN_STAN_DATA)
MODELS = $(INCIDENCE_MODEL_NAME)
CHAINS = $(addsuffix $(CHAINS_SUFFIX), \
   $(addprefix $(MCMC_CHAINS)/, $(MODELS))) 


#################################
# RULES
#################################

##########################
# rules for data cleaning
##########################
DATA_PROCESSING_SCRIPT = $(SRC)/make_epi_table.R
STAN_DATA_CLEANING_SCRIPT = $(SRC)/clean_data_for_stan.R

$(CLEAN_STAN_DATA): $(STAN_DATA_CLEANING_SCRIPT) $(CLEAN_EPI_DATA) $(RAW_CLIMATE_DATA)
	$(MKDIR) $(CLEANED_DATA_DIR)
	$(R_COMMAND) $^ $@


##########################
# rules for Stan models
##########################

# generic recipe for model chains output
$(MCMC_CHAINS)/%$(CHAINS_SUFFIX): $(SRC)/%.stan $(MODEL_FITTING_SCRIPT) $(STAN_DATA)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $(MODEL_FITTING_SCRIPT) $< $(MODEL_DATA) $@


#################################
# figures
##################################

FIGURES := figure_final_size_difference.$(FIGEXT) figure_susceptibility_distribution.$(FIGEXT)

FIG_PATHS := $(addprefix $(OUT)/, $(FIGURES))

FIG_DEPS = flu_final_size_model.py model_parameters.py
FIG_DEP_PATHS = $(addprefix $(SRC)/, $(FIG_DEPS)) 

$(OUT)/figure%.$(FIGEXT): $(SRC)/figure%.$(SRCEXT) $(FIG_DEP_PATHS)
	$(MKDIR) $(dir $@)
	./$< $@

.PHONY: figs
figs: $(FIG_PATHS)
	@echo $(FIG_PATHS)

#################################
# cleanup
##################################

GARBAGE := $(SRC)/*~ $(SRC)/__pycache__/ 

.PHONY: clean
clean:
	$(RM) -r $(OUT)
	$(RM) -r $(GARBAGE)
