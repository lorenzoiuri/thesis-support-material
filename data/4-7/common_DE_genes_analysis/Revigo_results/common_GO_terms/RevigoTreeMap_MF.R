# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");
setwd("/mnt/hdd0/univ/lab2/04-26/PANTHER/common/full_noFDR/revigo")

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003674","molecular_function",100.000,7.480,1.000,0.000,"molecular_function"),
c("GO:0003777","microtubule motor activity",0.370,3.240,0.759,0.000,"microtubule motor activity"),
c("GO:0050512","lactosylceramide 4-alpha-galactosyltransferase activity",0.006,1.340,0.957,0.106,"microtubule motor activity"),
c("GO:0005488","binding",89.584,7.804,1.000,0.000,"binding"),
c("GO:0015643","toxic substance binding",0.055,1.455,0.997,0.000,"toxic substance binding"),
c("GO:0016887","ATPase",2.615,1.358,1.000,0.000,"ATPase"),
c("GO:0045237","CXCR1 chemokine receptor binding",0.006,1.340,0.950,0.001,"CXCR1 chemokine receptor binding"),
c("GO:0005515","protein binding",76.443,10.059,0.990,0.003,"protein binding"),
c("GO:0044877","protein-containing complex binding",6.861,1.991,0.994,0.008,"protein-containing complex binding"),
c("GO:0035639","purine ribonucleoside triphosphate binding",10.045,2.418,0.554,0.009,"purine ribonucleoside triphosphate binding"),
c("GO:0043168","anion binding",13.219,1.983,0.942,0.339,"purine ribonucleoside triphosphate binding"),
c("GO:1901265","nucleoside phosphate binding",11.864,2.445,0.792,0.341,"purine ribonucleoside triphosphate binding"),
c("GO:0097367","carbohydrate derivative binding",12.400,1.467,0.994,0.010,"carbohydrate derivative binding"),
c("GO:0036094","small molecule binding",13.755,1.613,0.994,0.010,"small molecule binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_MF.pdf", width=10.5, height=10.5 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Molecular function GO terms enriched in TP53 silencing and in LINC01605 knockout experiments",  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

