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
SRCEXT := R
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
RAW_EPI_DATA = $(RAW_DATA_DIR)/epi_table.csv

CLEAN_EPI_DATA = $(CLEANED_DATA_DIR)/epi_table_with_climate.csv
CLEAN_STAN_DATA = $(CLEANED_DATA_DIR)/clean_stan_data.csv

################################
# for Bayesian analyses in Stan
################################
MODEL_FITTING_SCRIPT = $(SRC)/fit_stan_model.R
NORMED_INCIDENCE_MODEL_NAME = normed_multilevel_incidence_regression
MODEL_DATA = $(CLEAN_STAN_DATA)
MODELS = $(NORMED_INCIDENCE_MODEL_NAME)
CHAINS = $(addsuffix $(CHAINS_SUFFIX), \
   $(addprefix $(MCMC_CHAINS)/, $(MODELS)))

.PHONY: chains

chains: $(CHAINS)


#################################
# RULES
#################################

##########################
# rules for data cleaning
##########################
DATA_PROCESSING_SCRIPT = $(SRC)/make_epi_table.R
EPI_CLEANING_SCRIPT = $(SRC)/add_climate_to_epi_table.R
STAN_DATA_CLEANING_SCRIPT = $(SRC)/clean_data_for_stan.R

$(CLEAN_STAN_DATA): $(STAN_DATA_CLEANING_SCRIPT) $(CLEAN_EPI_DATA)
	$(MKDIR) $(CLEANED_DATA_DIR)
	$(R_COMMAND) $^ $@

$(CLEAN_EPI_DATA): $(EPI_CLEANING_SCRIPT) $(RAW_EPI_DATA) $(RAW_CLIMATE_DATA)
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

PLOTTING_STYLE = $(SRC)/plotting_style.R

FIGURES := figure_final_size_difference.$(FIGEXT) figure_susceptibility_distribution.$(FIGEXT) figure_posterior_effects.$(FIGEXT) figure_posterior_by_subtype.$(FIGEXT)

FIG_PATHS := $(addprefix $(OUT)/, $(FIGURES))

FIG_DEPS = flu_final_size_model.py model_parameters.py
FIG_DEP_PATHS = $(addprefix $(SRC)/, $(FIG_DEPS)) 

$(OUT)/figure_posterior_effects.$(FIGEXT): $(SRC)/figure_posterior_effects.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_INCIDENCE_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@

$(OUT)/figure_posterior_by_subtype.$(FIGEXT): $(SRC)/figure_posterior_by_subtype.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_INCIDENCE_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@

.PHONY: figs
figs: $(FIG_PATHS)
	@echo $(FIG_PATHS)

##################################
# reported quantities and diagnostics
##################################

$(OUT)/chain_diagnostics.csv: $(SRC)/chain_diagnostics.R $(CHAINS) $(CLEAN_STAN_DATA)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@

$(OUT)/reported_quantities.csv: $(SRC)/extract_reported_quantities.R $(CHAINS) $(CLEAN_STAN_DATA)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@

#################################
# cleanup
##################################

GARBAGE := $(SRC)/*~ $(SRC)/__pycache__/ 

.PHONY: clean
clean:
	$(RM) -r $(OUT)
	$(RM) -r $(GARBAGE)

