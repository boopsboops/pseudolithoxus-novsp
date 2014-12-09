# load libs
require(knitr)
require(markdown)
require(rmarkdown)
# first need to run 'sudo apt-get install pandoc-citeproc'

#cd
setwd("/home/rupert/Projects/pseudolithoxus-novsp/manuscript")

# run the R code in knitr and render the document in HTML format
knit("manuscript-master.Rmd") # ?knit
render("manuscript-master.Rmd", "html_document") #?render
