# load libs
require(knitr)
require(rmarkdown)
require(ape)
require(phangorn)
require(spider)

# first need to run 'sudo apt-get install pandoc-citeproc'

# cd
setwd("/home/rupert/Projects/pseudolithoxus-novsp/manuscript")

# run the R code in knitr and render the document in HTML format
render("manuscript-master.Rmd", "html_document") #?render
