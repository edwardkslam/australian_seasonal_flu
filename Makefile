#####################################
# name: Makefile
# author: Dylan Morris
# <dhmorris@princeton.edu>
#
# Makefile to generate analyses
# for Lam et al 2020 study of
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

#####################################
# directory structure
#####################################
SRC := src
DAT := dat
OUT := out

MCMC_CHAINS := $(OUT)/mcmc_chains
RAW_DATA_DIR := $(DAT)/raw
CLEANED_DATA_DIR := $(DAT)/cleaned
PLOT_PATH = $(OUT)/figures

FIGEXT := .png .eps
SRCEXT := R
CHAINS_SUFFIX = _chains.Rds

#####################################
## all targets
#####################################

all: depend chains summary_stats figs


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
NORMED_WITH_TEMP_MODEL_NAME = normed_multilevel_incidence_regression_temp
EPIDEMIC_PROBABILITY_MODEL_NAME = multilevel_epi_probability_model

MODEL_DATA = $(CLEAN_STAN_DATA)

MODELS = $(NORMED_INCIDENCE_MODEL_NAME) $(NORMED_WITH_TEMP_MODEL_NAME)

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

## generic recipe for model chains output
$(MCMC_CHAINS)/%$(CHAINS_SUFFIX): $(SRC)/%.stan $(MODEL_FITTING_SCRIPT) $(MODEL_DATA)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $(MODEL_FITTING_SCRIPT) $< $(MODEL_DATA) $@


#################################
# figures
##################################

PLOTTING_STYLE = $(SRC)/plotting_style.R

FIG_CLEANUP = @$(RM) Rplots.pdf

FIGURES := figure_posterior_effects figure_posterior_by_subtype figure_posterior_sds figure_posterior_effects_temp figure_posterior_by_subtype_temp figure_posterior_sds_temp

FIG_PATHS := $(foreach EXT, $(FIGEXT), $(addsuffix $(EXT), $(addprefix $(OUT)/, $(FIGURES))))

.PHONY: figs
figs: $(FIG_PATHS)
	@echo $(FIG_PATHS)


$(OUT)/figure_posterior_effects.%: $(SRC)/figure_posterior_effects.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_INCIDENCE_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(OUT)/figure_posterior_effects_temp.%: $(SRC)/figure_posterior_effects.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_WITH_TEMP_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(OUT)/figure_posterior_by_subtype.%: $(SRC)/figure_posterior_by_subtype.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_INCIDENCE_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(OUT)/figure_posterior_by_subtype_temp.%: $(SRC)/figure_posterior_by_subtype.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_WITH_TEMP_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


$(OUT)/figure_posterior_sds.%: $(SRC)/figure_posterior_sds.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_INCIDENCE_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(OUT)/figure_posterior_sds_temp.%: $(SRC)/figure_posterior_sds.$(SRCEXT) $(MCMC_CHAINS)/$(NORMED_WITH_TEMP_MODEL_NAME)$(CHAINS_SUFFIX) $(CLEAN_STAN_DATA) $(PLOTTING_STYLE)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)



##################################
# reported quantities and diagnostics
##################################

$(OUT)/chain_diagnostics.csv: $(SRC)/chain_diagnostics.R $(CHAINS) $(CLEAN_STAN_DATA)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@

$(OUT)/reported_quantities.csv: $(SRC)/extract_reported_quantities.R $(CHAINS) $(CLEAN_STAN_DATA)
	$(MKDIR) $(dir $@)
	$(R_COMMAND) $^ $@

.PHONY: summary_stats

summary_stats: $(OUT)/chain_diagnostics.csv $(OUT)/reported_quantities.csv

#################################
# cleanup
##################################

GARBAGE := $(SRC)/*~ $(SRC)/__pycache__/ 

.PHONY: clean
clean:
	$(RM) -r $(OUT)
	$(RM) -r $(GARBAGE)

