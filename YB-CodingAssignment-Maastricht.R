################################################################################

# Yasmijn Balder
# 25-07-2021
# Coding assignment for PhD position AT2021.240 in Maastricht

################################################################################

# Import packages
library(biomaRt)
library(RCy3)
library(stringr)
library(rWikiPathways)

# Import data
variants <- read.csv("./variant_list.txt", col.names = 'rsID', header = FALSE)

################################################################################

## Step 1: Retrieve the associated gene for each variant using BioMart

# Create empty column to store the associated genes
variants$gene <- rep('', nrow(variants))

# Create bioMart database
# listEnsembl()
# variation = useEnsembl(biomart = "snps")
# listDatasets(variation)
variation = useEnsembl(biomart = "snps", dataset = "hsapiens_snp", 
                       mirror = 'uswest')

# Create an empty data frame to store hits of rs with multiple genes associated
multiple.hits <- data.frame('','')
# Give this dataframe the same names, in order to combine then later
names(multiple.hits) <- c('rsID', 'gene')

# Obtain the associated gene for each variant and store in the data.frame
for (rsnum in 1:nrow(variants)) {
  # If there is 1 hit
  if (dim(getBM(attributes = 'associated_gene', 
                       filters = 'snp_filter', mart = variation, 
                       values = variants$rsID[rsnum]))[1] == 1){
  # Add this gene to the dataframe
  variants$gene[rsnum] <- getBM(attributes = 'associated_gene', 
                                 filters = 'snp_filter', mart = variation, 
                                 values = variants$rsID[rsnum])[,1]
  # If there are multiple hits
  } else if (dim(getBM(attributes = 'associated_gene', 
                       filters = 'snp_filter', mart = variation, 
                       values = variants$rsID[rsnum]))[1] > 1){
    # Find all unique hits
    hits <-  getBM(attributes = 'associated_gene',
           filters = 'snp_filter', mart = variation, 
           values = variants$rsID[rsnum])[,1]
    separated <- str_split(hits, ',')
    combined <- unique(unlist(separated))
    # Add the first hit to the data frame
    variants$gene[rsnum] <- combined[1]
    # Save the other hits to add after this loop
    for (hit in 2:length(combined)) {
      additional <- c(variants$rsID[rsnum], combined[hit])
      names(additional) <- c('rsID', 'gene')
      multiple.hits <- rbind(multiple.hits, additional)
    }
  }
}

# Combine the first hits and additional hits to obtain them all
edge.list <- rbind(variants, multiple.hits)
# Turn NA gene entry into empty entry
edge.list[which(is.na(edge.list$gene)),] <- ""
# Remove empty row
edge.list <- edge.list[-which(edge.list$rsID == ""),]

################################################################################

## Step 2: Create a network in Cytoscape linking the SCPs to the genes

# Connect to cytoscape (should be opened!)
cytoscapePing()

# Create a dataframe with all rsIDs and genes
nodes <- data.frame(id = c(unique(edge.list$rsID), 
                           # Don't create empty nodes for genes
                           unique(edge.list$gene)[-which(unique(edge.list$gene) == "")]),
                    group = c(rep('rs', length(unique(edge.list$rsID))),
                              # Don't create empty nodes for genes
                              rep('gene', length(unique(edge.list$gene)[- which(unique(edge.list$gene) == "")]))),
                    stringsAsFactors = FALSE)
# Create a dataframe to explain the edges
edge.list <- edge.list[-which(edge.list$gene == ""),]
edges <- data.frame(source = edge.list$rsID, 
                    target = edge.list$gene,
                    interaction = rep("interacts", nrow(edge.list)),
                    stringsAsFactors = FALSE)

# Open the network in cytoscape
createNetworkFromDataFrames(nodes, edges, title = "var and genes", 
                            collection = "YB_AT2021.240")

# Create style
style.name = "myStyle"
defaults <- list(NODE_SHAPE = "circle",
                 NODE_SIZE = 70,
                 EDGE_TRANSPARENCY = 120,
                 NODE_LABEL_COLOR="#FFFFFF",
                 NODE_BORDER_PAINT="#FFFFFF")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','group','d',c("rs","gene"), 
                               c("#56B4E9","#E69F00"))
# Apply style
createVisualStyle(style.name, defaults, list(nodeLabels, nodeFills))
setVisualStyle(style.name)

################################################################################

## Step 3: Add known pathways

# Make sure the cytoscape app is installed
installApp('WikiPathways')

# Extend the network in cytoscape with the pathways (this step doesn't work)
extend.cmd = paste('cytargetlinker extend idAttribute="shared name" linkSetDirectory="', 
                   paste(getwd(), 'LinkSets', sep = '/'),
                   '" network=current direction=TARGETS', sep="")
commandsRun(extend.cmd)
layoutNetwork()

# Filter out genes without any pathway information
filter1.cmd = "network select edgeList=all"
filter2.cmd = "network select extendEdges=true"
filter3.cmd = "network create nodeList=selected edgeList=selected networkName=selection source=current"
commandsRun(filter1.cmd)
commandsRun(filter2.cmd)
commandsRun(filter3.cmd)

################################################################################

## Step 4: Save session and image

saveSession('YB_AT2021240_session')
png.path = paste(getwd(),'YB_AT2021240_image', sep = '/')
exportImage(png.path, 'PNG')
