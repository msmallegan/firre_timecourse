######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)

# Treatment vs. Control
rmarkdown::render(file.path("01_firre_induction_vs_control/",
                            "firre_responder_differential_expression.Rmd"), 
                  md_document(variant = "markdown_github"))
