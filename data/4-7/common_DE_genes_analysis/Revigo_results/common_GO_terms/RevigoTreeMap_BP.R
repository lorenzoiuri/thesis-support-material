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
revigo.data <- rbind(c("GO:0008150","biological_process",100.000,2.745,1.000,0.000,"biological_process"),
c("GO:0009987","cellular process",86.803,2.939,1.000,0.000,"cellular process"),
c("GO:0051179","localization",32.339,3.275,1.000,0.000,"localization"),
c("GO:0051606","detection of stimulus",3.999,3.462,0.979,0.000,"detection of stimulus"),
c("GO:0006979","response to oxidative stress",2.131,2.629,0.954,0.111,"detection of stimulus"),
c("GO:0010243","response to organonitrogen compound",5.597,1.374,0.950,0.129,"detection of stimulus"),
c("GO:0042542","response to hydrogen peroxide",0.679,2.545,0.905,0.291,"detection of stimulus"),
c("GO:0062197","cellular response to chemical stress",1.526,1.752,0.934,0.320,"detection of stimulus"),
c("GO:0071549","cellular response to dexamethasone stimulus",0.163,1.365,0.941,0.388,"detection of stimulus"),
c("GO:0090316","positive regulation of intracellular protein transport",1.010,2.378,0.693,0.000,"positive regulation of intracellular protein transport"),
c("GO:0001954","positive regulation of cell-matrix adhesion",0.314,1.876,0.869,0.152,"positive regulation of intracellular protein transport"),
c("GO:2000379","positive regulation of reactive oxygen species metabolic process",0.387,1.423,0.885,0.156,"positive regulation of intracellular protein transport"),
c("GO:0031323","regulation of cellular metabolic process",34.139,1.857,0.899,0.179,"positive regulation of intracellular protein transport"),
c("GO:0043065","positive regulation of apoptotic process",3.017,1.654,0.849,0.195,"positive regulation of intracellular protein transport"),
c("GO:0048522","positive regulation of cellular process",31.548,2.255,0.816,0.295,"positive regulation of intracellular protein transport"),
c("GO:1903954","positive regulation of voltage-gated potassium channel activity involved in atrial cardiac muscle cell action potential repolarization",0.006,1.340,0.791,0.390,"positive regulation of intracellular protein transport"),
c("GO:0030335","positive regulation of cell migration",2.838,1.401,0.759,0.395,"positive regulation of intracellular protein transport"),
c("GO:0051041","positive regulation of calcium-independent cell-cell adhesion",0.006,1.340,0.908,0.469,"positive regulation of intracellular protein transport"),
c("GO:0031325","positive regulation of cellular metabolic process",18.132,1.559,0.805,0.480,"positive regulation of intracellular protein transport"),
c("GO:0060341","regulation of cellular localization",4.565,1.318,0.804,0.484,"positive regulation of intracellular protein transport"),
c("GO:0140014","mitotic nuclear division",1.139,2.496,0.899,0.000,"mitotic nuclear division"),
c("GO:0030199","collagen fibril organization",0.566,1.389,0.968,0.160,"mitotic nuclear division"),
c("GO:0016043","cellular component organization",32.036,3.539,0.950,0.262,"mitotic nuclear division"),
c("GO:0048285","organelle fission",2.075,1.442,0.947,0.312,"mitotic nuclear division"),
c("GO:0007010","cytoskeleton organization",6.803,2.073,0.939,0.397,"mitotic nuclear division"),
c("GO:0008608","attachment of spindle microtubules to kinetochore",0.129,1.668,0.952,0.497,"mitotic nuclear division"),
c("GO:0016322","neuron remodeling",0.062,2.344,0.978,0.003,"neuron remodeling"),
c("GO:0034501","protein localization to kinetochore",0.084,1.397,0.970,0.003,"protein localization to kinetochore"),
c("GO:0051234","establishment of localization",25.782,2.197,0.954,0.181,"protein localization to kinetochore"),
c("GO:0008104","protein localization",11.374,1.327,0.955,0.433,"protein localization to kinetochore"),
c("GO:0033036","macromolecule localization",13.735,2.052,0.959,0.457,"protein localization to kinetochore"),
c("GO:0050673","epithelial cell proliferation",0.494,1.485,0.998,0.004,"epithelial cell proliferation"),
c("GO:0051301","cell division",2.737,2.186,0.998,0.005,"cell division"),
c("GO:0071840","cellular component organization or biogenesis",33.281,3.636,0.997,0.008,"cellular component organization or biogenesis"),
c("GO:0006928","movement of cell or subcellular component",8.828,1.886,0.997,0.011,"movement of cell or subcellular component"),
c("GO:0060668","regulation of branching involved in salivary gland morphogenesis by extracellular matrix-epithelial cell signaling",0.006,1.340,0.976,0.019,"regulation of branching involved in salivary gland morphogenesis by extracellular matrix-epithelial cell signaling"),
c("GO:1900450","negative regulation of glutamate receptor signaling pathway",0.006,1.340,0.974,0.023,"negative regulation of glutamate receptor signaling pathway"),
c("GO:0043086","negative regulation of catalytic activity",4.324,1.467,0.916,0.027,"negative regulation of catalytic activity"),
c("GO:0048872","homeostasis of number of cells",1.251,1.391,0.963,0.028,"homeostasis of number of cells"),
c("GO:0007088","regulation of mitotic nuclear division",0.595,1.721,0.892,0.034,"regulation of mitotic nuclear division"),
c("GO:0033043","regulation of organelle organization",6.741,2.103,0.878,0.449,"regulation of mitotic nuclear division"),
c("GO:0050678","regulation of epithelial cell proliferation",1.896,1.314,0.955,0.039,"regulation of epithelial cell proliferation"),
c("GO:0065009","regulation of molecular function",17.257,2.194,0.955,0.043,"regulation of molecular function"),
c("GO:0051128","regulation of cellular component organization",13.214,1.672,0.939,0.056,"regulation of cellular component organization"),
c("GO:0051726","regulation of cell cycle",6.405,1.309,0.946,0.070,"regulation of cell cycle"),
c("GO:0007605","sensory perception of sound",0.847,1.485,0.995,0.096,"sensory perception of sound"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_BP.pdf", width=10.5, height=15 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Biological process GO terms enriched in TP53 silencing and in LINC01605 knockout experiments",
  inflate.labels = F,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

