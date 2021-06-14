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
revigo.data <- rbind(c("GO:0005575","cellular_component",100.000,2.199,1.000,0.000,"cellular_component"),
c("GO:0005871","kinesin complex",0.271,3.197,0.776,0.000,"kinesin complex"),
c("GO:0044599","AP-5 adaptor complex",0.005,1.340,0.946,0.113,"kinesin complex"),
c("GO:0043235","receptor complex",2.602,1.682,0.973,0.187,"kinesin complex"),
c("GO:0043232","intracellular non-membrane-bounded organelle",26.812,2.327,0.778,0.266,"kinesin complex"),
c("GO:0072686","mitotic spindle",0.830,1.848,0.752,0.465,"kinesin complex"),
c("GO:0015630","microtubule cytoskeleton",6.901,2.597,0.746,0.491,"kinesin complex"),
c("GO:0110165","cellular anatomical entity",99.079,2.260,1.000,0.000,"cellular anatomical entity"),
c("GO:0062023","collagen-containing extracellular matrix",2.261,1.631,1.000,0.000,"collagen-containing extracellular matrix"),
c("GO:0005622","intracellular anatomical structure",78.328,6.363,1.000,0.000,"intracellular anatomical structure"),
c("GO:0012505","endomembrane system",24.572,1.670,1.000,0.000,"endomembrane system"),
c("GO:0031974","membrane-enclosed lumen",29.281,1.590,1.000,0.000,"membrane-enclosed lumen"),
c("GO:0043226","organelle",73.523,5.717,1.000,0.000,"organelle"),
c("GO:0005737","cytoplasm",63.105,8.035,0.954,0.017,"cytoplasm"),
c("GO:0005741","mitochondrial outer membrane",1.064,2.386,0.788,0.022,"mitochondrial outer membrane"),
c("GO:0030659","cytoplasmic vesicle membrane",4.203,1.775,0.792,0.345,"mitochondrial outer membrane"),
c("GO:0005765","lysosomal membrane",2.096,1.391,0.778,0.457,"mitochondrial outer membrane"),
c("GO:0098588","bounding membrane of organelle",11.424,2.049,0.795,0.499,"mitochondrial outer membrane"),
c("GO:0043227","membrane-bounded organelle",69.128,6.339,0.914,0.033,"membrane-bounded organelle"),
c("GO:0070013","intracellular organelle lumen",29.281,1.747,0.887,0.102,"membrane-bounded organelle"),
c("GO:0043231","intracellular membrane-bounded organelle",63.542,3.544,0.886,0.198,"membrane-bounded organelle"),
c("GO:0043229","intracellular organelle",69.272,4.604,0.902,0.221,"membrane-bounded organelle"),
c("GO:0048471","perinuclear region of cytoplasm",3.858,1.724,0.962,0.044,"perinuclear region of cytoplasm"),
c("GO:0019867","outer membrane",1.203,1.860,0.959,0.072,"outer membrane"),
c("GO:0031090","organelle membrane",19.230,2.411,0.898,0.107,"outer membrane"),
c("GO:0005829","cytosol",28.349,2.030,0.944,0.075,"cytosol"),
c("GO:0043228","non-membrane-bounded organelle",26.854,2.323,0.940,0.097,"non-membrane-bounded organelle"),
c("GO:0005783","endoplasmic reticulum",10.562,1.783,0.891,0.097,"endoplasmic reticulum"),
c("GO:0005794","Golgi apparatus",8.434,1.472,0.895,0.449,"endoplasmic reticulum"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap_CC.pdf", width=10.5, height=15 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Cellular component GO terms enriched in TP53 silencing and in LINC01605 knockout experiments",  
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

