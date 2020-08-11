all: 
	
clean:
	rm -r output
	$(info you have a clean slate)

housekeeping:
	mkdir -p output/html
	mkdir -p output/figs

### Process raw MiSeq data ###
MiSeq: output/html/MiSeqProcessing.html
output/html/MiSeqProcessing.html output/rds/phy.rds output/rds/phy.erm.rds: code/MiSeqProcessing.R | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'

### Analyses ###
analysis: output/rds/phy.rds output/rds/phy.erm.rds \
	output/figs/Sankey.pdf \
	output/html/ChronosequenceComposition.html output/figs/ordinations.pdf \
	output/html/ChronosequenceDiversity.html output/figs/divPlot.pdf \
	output/html/FertilizationPlots.html output/figs/fertDivPlot.pdf \

output/figs/Sankey.pdf: code/SankeyPlot.R output/rds/phy.rds | housekeeping
	Rscript $<
	
output/html/ChronosequenceComposition.html output/figs/ordinations.pdf: code/ChronosequenceComposition.R output/rds/phy.rds output/rds/phy.erm.rds | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'

output/html/ChronosequenceDiversity.html output/figs/divPlot.pdf: code/ChronosequenceDiversity.R output/rds/phy.rds output/rds/phy.erm.rds | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'

output/html/FertilizationPlots.html output/figs/fertDivPlot.pdf: code/FertilizationPlots.R output/rds/phy.rds output/rds/phy.erm.rds | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'
	
	