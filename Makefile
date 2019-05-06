#################################
# commands
##################################
RM = rm -f
CXXFLAGS = -O3 -std=c++11
MKDIR = @mkdir -p

#################################
# directory structure
##################################
SRCDIR := src
DATA := dat
OUT := out

FIGEXT := png
SRCEXT := py

#################################
# figures
##################################

FIGURES := figure_final_size_difference.$(FIGEXT) figure_susceptibility_distribution.$(FIGEXT)

FIG_PATHS := $(addprefix $(OUT)/, $(FIGURES))

FIG_DEPS = flu_final_size_model.py model_parameters.py
FIG_DEP_PATHS = $(addprefix $(SRCDIR)/, $(FIG_DEPS)) 

$(OUT)/figure%.$(FIGEXT): $(SRCDIR)/figure%.$(SRCEXT) $(FIG_DEP_PATHS)
	$(MKDIR) $(dir $@)
	./$< $@

.PHONY: figs
figs: $(FIG_PATHS)
	@echo $(FIG_PATHS)

#################################
# cleanup
##################################

GARBAGE := $(SRCDIR)/*~ $(SRCDIR)/__pycache__/ 

.PHONY: clean
clean:
	$(RM) -r $(OUT)
	$(RM) -r $(GARBAGE)

