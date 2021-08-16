Contamination_genes <- c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1")
Prolif_genes <- c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4")
RBC_genes <- c("HBB", "HBA1", "HBA2", "HBD")

Endo_genes_SN_SC <- c("CTGF", "FCGR2B", "S100A13", "FCN2", "FCN3", "LYVE1", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
				"F8", "CALCRL", "SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF",
"ENG", "PECAM1", "RAMP3", "INMT", "DNASE1L3", "LIFR", "TIMP3", "C7"
				)
Endo_genes <- c("EPCAM", "PECAM1", "VWF", "CLEC4G", "CLEC4M", 
			"CD34", "CD14", "LYVE1", "RSPO3", "WNT2", "COL1A2", "TFF3",
			"ENG", "STAB2", "CLDN5", "SPARCL1", "RBP7")
Endo_cvLSEC_vs_ppLSEC <- c("CTSD", "CTSL", "CLEC1B", "MS4A6A", 
				"STAB1", "CLEC4G", "CRHBP", "DNASE1L3", "FCN2", "FCN3")
Endo_ppLSEC_vs_cvLSEC_portEndo <- c("MGP", "VIM", "ADIRF", "SPARCL1", "CLU", "S100A6", 
							"CD9", "CLEC14A", "AQP1", "TM4DF1")

Endo_genes_dot <- unique(c("FCN2", "FCN3", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
			"SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF", 
			"ENG", "PECAM1", "RAMP3", "DNASE1L3", "LIFR", "C7", "TIMP3",
			"ACKR1", "WNT2", "RSPO3", "INMT", "PLAC8",
			"PODXL", "RBP7", "PLVAP", "GSN", "CD34",
			"CCL21", "S100A6", "CST3", "ALB", "APOA1"))

Endo_genes_SN_SC, "ACKR1", "WNT2", "RSPO3", 
		Endo_cvLSEC_vs_ppLSEC, Endo_ppLSEC_vs_cvLSEC_portEndo, 
		"CCL21", "MYL9", "FABP4", "IGFBP5", "COL1A2", "LGALS1", "CRYAB", 
		"PODXL", "JAG2", "RBP7", "EBF`", "PLVAP", "GSN", "CD34", "ADAM15", "TIMP3"))

Macrophage_genes <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", 
				"CD163", "HLA-DRA", "HLA-DPA1", "CD74", "FABP5",
				"PLAUR", "LYZ", "S100A4", "S100A8", "VCAN", "FCN1")

Macrophage_genes_dot <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", "FCGR3A", "FCER1G",
				"CD163", "HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "VSIG4", "CD74", 
				"GPNMB", "ACTP5", "FABP5", "SPP1", "FABP4", "TREM2", "LGALS3", "CTSB", "PSAP", "APOE", "APOC1", "NPC2",
				 "LYZ",  "S100A4", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1", "FTH1", "CD68", 
				"PLAUR", "SRGN", "AREG", "THBS1", "CXCL3", "IL1B", "CCL3", "LILRA5", "ILRB2", "PLAC8", "CD52", 
				 "LST1", Contamination_genes)

Macrophage_genes_dot <- c(
		"MARCO", "CD5L", "LYVE1", "SLC40A1", "FTL", "CD163", "SEPP1", "C1QC", "C1QB", "C1QA", "CTSB", "HMOX1", "VCAM1", 
		"HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "CD74", "VSIG4", 
		"LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1",
		"FABP5", "ACP5", "PLD3", "FTH1", "CD68", "APOE", "PSAP", "CSTB", "LGMN",
		"RBP7", "FOLR2", "FCER1G", "MS4A7", "TIMP1",
		"JUND", "FOS", "NFKBIA", "ACTG1", "CD14", "CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B",
		"PLAUR", "SRGN", "AREG", "THBS1", "CXCL3", "IL1B", "CCL3",
		"PLAC8", "CD54", "LST1", "IFITM3", "AIF1", "COTL1",
		"DNASE1L3", "FCN2", "CCL14", "FCN3", "SPARC", "CLEC1B", "ENG",
		"ALB", "SERPINA1", "APOA1", "HP", "FGA")

NKT_genes <- c("CD8A", "CD3D", "TRAC", "TRBC2", "TRDC", "GNLY", "GZMB", "GZMA", "CCL5", 
			"NKG7", "FCGR3A", "FGFBP2", "CD8B", "IL7R", "CD74", "HLA-DRB1")
NKT_exhaustion = c("PDCD1") 
NKT_Treg = c("ITGAE", "FOXP3")
NKT_genes_dot_old <- c("CD52","NCR1", "CXCR6", "CXCR3", "CD69", "PTPRC", "CD8A", "CD3D", "CD3E", "CD8B", "CD4", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", 
				"IL32", "IL7R", "LTB", "IL17A", "IL18R1", "CD44", "KLRB1", "KLRC1", "KLRF1", "KLRK1", "CCL4", "CCL5", "NKG7",
				"FCGR3A", "FGFBP2", "IL2RB", "GZMA", "GZMB", "CSF2", "GNLY", "KLRD1",
				"CD74", "CD79A", "CD79B", "HLA-DRB1", "HLA-DRA", "AIF1", "PDCD1",
				Prolif_genes, Contamination_genes, RBC_genes)

NKT_genes_dot <- c("CD3D", "CD3E", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", 
				"CD8A", "CD8B", "CCL3", "CCL4", "CCL5", "IL7R", "LTB", "KLRB1", "TPT1",
				"KLRC1", "KLRF1", "GZMK", "CMC1", "XCL1", "XCL2", "GZMB", "FCGR3A", "GNLY", "CXCR6", "CD69", "EOMES", "TBX21",
				"CD79A", "CD79B", "CD74", "HLA-DRB1", "HLA-DRA", "HLA-DPA1", "HLA-DQB1",
				"AIF1", "VSIG4", "LYZ", "VCAN", "CD163", "C1QC", "ITGAM", "ITGAE", "CST3",
				"MKI67", "BRIC5", "TOP2A", "CDK1", "HMGB2", "HBB", "HBA1", "HBA2", "HBD",
				"ALB", "SERPINA1", "APOA1", "FGA" 
				)

require(Seurat)
require(ggplot2)

obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Macrophage_harmony_Subcluster.rds")
cluster_col = "Coarse_clusters"
manual_cluster_anno <- c(
		"InfSynap", "NonInf", "NonInf", "Inflam", "RetNonInf", "InflamActiv", 
		"Inflam", "MHCII", "Activated", "MHCII", "PhagoNonInf", 
		"Debris", "Inflam", "Doublet", "NonInf")
# Changed 25 May 2021 so that Cluster 2 (third in list) is "NonInf" rather than "Debris", 
# and Cluster 8 is "Activated" instead of "Repair"
obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]


FeaturePlot(obj, Macrophage_genes)

DotPlot(obj, features=unique(Macrophage_genes_dot), group.by="Coarse_clusters", split.by="assay_type", cols=c("red", "blue"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/Colour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/SubColour_Scheme.R")

flow_cyto_genes <- c("CD68", "MRC1", "CD14", "CD274", "PTPRC", "CD86", "CD163", "IL2RA", "FOXP3",
				"NCAM1", "NCR1", "CD3E", "PTPRC", "CD8A", "IL7R", "ICOS", "KLRD1", "IL2RA",
				"CD4", "PDCD1", "SELL", "HAVCR2", "CCR3", "CD3E", "IL7R", "CD8A", "PTPRC", 
				"LAG3", "CTLA4", "CD27", "CD4")

#### Integrate Blair's Macropahge Data ####
set.seed(101)
blair <- Read10X("C:/Users/tandrews/Documents/UHNSonya/RawData/BlairMacs/Blair_Macs_cellranger_out/filtered_feature_bc_matrix")
blair <- blair[,Matrix::colSums(blair) > 875 & Matrix::colSums(blair > 0) > 500]
blair <- CreateSeuratObject(blair, project="Blair")

#Match genes
blair <- blair[rownames(blair) %in% rownames(obj),]
obj_merge <- obj[rownames(obj) %in% rownames(blair),]
blair <- blair_merge[match( rownames(obj_merge), rownames(blair)),]

#Normalize
orig_norm_factor = 1500
blair_norm_factor = orig_norm_factor/median(Matrix::colSums(obj_merge@assays$RNA@counts)) * median(Matrix::colSums(blair@assays$RNA@counts))

blair <- NormalizeData(blair, normalizeation.method="LogNormalize", scale.factor=blair_norm_factor)

dim(blair)

#Scale -> UMAP
set.seed(189)
blair <- ScaleData(blair, features=rownames(blair))
blair <- RunPCA(blair, features=VariableFeatures(obj), dims=1:20)
blair <- RunUMAP(blair, dims=1:20)
DimPlot(blair)

## Cluster then match clusters ##
blair <- Seurat::FindNeighbors(blair, dims=1:20)
blair <- Seurat::FindClusters(blair)
DimPlot(blair)

saveRDS("Clustered_Blair.rds")

# Plan 1 -> merge together = Fail! too much background RNA

#merged_obj <- merge(obj_merge[ , obj_merge@meta.data$assay_type == "3pr"], blair_merge, add.cell.ids=c("Map2", "Blair"))
##Add on
#merged_obj@meta.data$Coarse_clusters[is.na(merged_obj@meta.data$Coarse_clusters)] <- "Blair"
#DotPlot(merged_obj, features=unique(Macrophage_genes_dot), group.by="Coarse_clusters")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#merged_obj <- ScaleData(merged_obj)
#merged_obj <- RunPCA(merged_obj, dims=20, features=VariableFeatures(obj))
#merged_obj <- RunUMAP(merged_obj, dims=1:20)
#DimPlot(merged_obj, group.by="sample") ## Mainly by Sample! - Plan Fail! - Scaling is necessary...

# Plan 2 -> map cells to nearest neighbours using spearman correlations.
require(Hmisc)
# Use Scaled data for reference.
# Subset to key markers from dot plot
match_features = Macrophage_genes_dot # VariableFeatures(obj_merge) # Macrophage_genes_dot
mat <- obj_merge@assays$RNA@scale.data[rownames(obj_merge) %in% match_features,]
mat_blair <- blair@assays$RNA@data[rownames(blair) %in% match_features,]
require(qlcMatrix)

cor_mat <- corSparse(mat, mat_blair)
match <- apply(cor_mat, 2, function(x){  obj_merge@meta.data$Subcluster_Manual[ which(x == max(x)) ]  })

blair@meta.data$matched_annotation <- match

table(match)

blair@meta.data$annotation_HVG <- blair@meta.data$matched_annotation
blair@meta.data$annotation_Key <- blair@meta.data$matched_annotation

table(blair@meta.data$annotation_HVG[blair@meta.data$annotation_HVG == blair@meta.data$annotation_Key])

DimPlot(blair, group.by="annotation_Key")

saveRDS("Clustered_Anno_Blair.rds")

# Plan 3 - Markers of clusters & match across.

#Renorm: blair@assays$RNA@data <- Matrix::t(Matrix::t(blair@assays$RNA@data)/colSums(blair@assays$RNA@data)*2500)


for (cluster in 0:13) {
	blair_markers <- FindMarkers(blair, group.by="seurat_clusters", ident.1=cluster, logfc.threshold=-Inf)
	head(blair_markers, 20)
	write.table(blair_markers, paste("DE_Cluster_Markers_c", cluster, ".txt", sep=""), quote=TRUE, row.names=TRUE, col.names=TRUE, sep=",")
}

DotPlot(blair, features=unique(Macrophage_genes_dot), group.by="seurat_clusters")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

X11()
DimPlot(blair, group.by="seurat_clusters")
FeaturePlot(blair, features="nFeature_RNA")

# Manual Annotation: 
# 0 = Activated
# 1
# 2
# 3
# 4
# 5 = Activated
## side clusters
# 6 = low expression / high MT
# 7 = markers = c("CADM1", "C6", "DNAJC15", "FOLR3", "OXCT1", "CXCL14")
# 8 = LSEC, markers = c("DNASE1L3", "FCN3", "CLEC1B")
# 9 = heading towards Inflammatory, markers = c("DNASE1L3", "FCN3", "CLEC1B")
# 10 = RBC, markers = c("ALAS2", "HBA2", "HBA1", "HBD")
# 11
# 12
# 13 = has hepatocyte genes


pseudobulks_20Livers <- get_pseudobulk_means(obj_merge@assays$RNA@data, obj@meta.data$Subcluster_Manual, obj@meta.data$sample)
pseudobulks_blair <- group_rowmeans(blair@assays$RNA@data, blair@meta.data$seurat_clusters)
profile_20Livers <- group_rowmeans(obj_merge@assays$RNA@scale.data, obj@meta.data$Subcluster_Manual)

pseudo_sim <- rcorr(pseudobulks_blair, pseudobulks_20Livers)
profile_sim <- rcorr(pseudobulks_blair[rownames(pseudobulks_blair) %in% Macrophage_genes_dot,], profile_20Livers[rownames(pseudobulks_blair) %in% Macrophage_genes_dot,], type="spearman")
require(gplots)
heatmap.2(profile_sim$r[1:14,15:24], trace="none")

### DE Blair vs LiverMap2.0


## Blair vs Devo_map

# source: https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/

# First lets get the gene expression profiles of all macrophages from the Devo map.
dir = "C:/Users/tandrews/Documents/UHNSonya/ExternalData/DevelopmentalCellAtlas"

#Which cells are macrophages?
devo_meta <- readRDS(paste(dir, "/df_cell.RDS", sep=""))

table(devo_meta[grepl("Myeloid", devo_meta[,15]),15])

objs <- Sys.glob(paste(dir, "/*count.RDS", sep=""))
genes_meta <- readRDS(paste(dir, "df_gene.RDS", sep="/"))


all_profiles <- c();
organ <- c();
best_genes <- c();
# Gather Myeloid data
for (this_i in 1:length(objs)) {

	this_obj <- readRDS(objs[this_i])
	this_counts <- this_obj[,colnames(this_obj) %in% devo_meta[grepl("Myeloid", devo_meta[,15]), 1] ]
	this_name <- unique(devo_meta[ devo_meta[,1] %in% colnames(this_obj) , "Organ"])

	# select genes >> myeloid than any other cell type
	# mean expression across all cell-types
	this_types <- devo_meta[match( colnames(this_obj), devo_meta[,1] ), 15] 
	type_profiles <- group_rowmeans(this_obj, this_types)
	
	if (ncol(this_counts) > 0) {
		profile <- Matrix::rowMeans(this_counts)
		all_profiles <- cbind(all_profiles, profile)

		# Store myeloid-specific genes
		macro <- grep("Myeloid", colnames(type_profiles))
		macro_genes <- apply(type_profiles, 1, function(x) {x[macro]/max(x[-macro]) > 10 & x[macro] == max(x) & max(x) > 0.05})
		best_genes <- c(as.character(best_genes), as.character(genes_meta[macro_genes,3]))
	
		organ <- c(organ, this_name)
	} else {
		next;
	}
}

rownames(all_profiles) <- genes_meta[,3]
colnames(all_profiles) <- organ

# unique genes
macrophage_genes <- unique(best_genes)

key_genes <- all_profiles[rownames(all_profiles) %in% macrophage_genes,]

#fix weird duplicates
key_genes <- key_genes[rowSums(key_genes) > 0,]

require(gplots)
heatmap.2(key_genes, dendrogram="none", trace="none", scale="row")

# match blair data
key_genes_blair <- pseudobulks_blair[match(rownames(key_genes), rownames(pseudobulks_blair)),]
key_genes_blair <- key_genes_blair[!is.na(key_genes_blair[,1]),]
key_genes <- key_genes[match(rownames(key_genes_blair), rownames(key_genes)),]
identical(rownames(key_genes_blair), rownames(key_genes))
# combine
full_key_genes <- cbind(key_genes, key_genes_blair)
rownames(full_key_genes) <- rownames(key_genes)
full_key_genes <- apply(full_key_genes, 2, scale) # This is much better than straight normalization.

# heatmap
heatmap.2(full_key_genes, dendrogram="none", trace="none", scale="none")

# correlations
require(Seurat)
heatmap.2( cor(full_key_genes, method="spearman")[ 1:ncol(key_genes), ( ncol(key_genes)+1 ):ncol( full_key_genes )] , 
		trace="none", mar=c(4,6), col=PurpleAndYellow(25), key.title="", key.xlab="Spearman correlation")

### Liver-specific macrophage-markers

# make table of which genes are markers in which tissues
macrophage_gene_table <- matrix(0, nrow=length(macrophage_genes), ncol=length(organ))
colnames(macrophage_gene_table) <- organ
rownames(macrophage_gene_table) <- macrophage_genes

for (this_i in 1:length(objs)) {

	this_obj <- readRDS(objs[this_i])
	this_name <- unique(devo_meta[ devo_meta[,1] %in% colnames(this_obj) , "Organ"])

	# select genes >> myeloid than any other cell type
	# mean expression across all cell-types
	this_types <- devo_meta[match( colnames(this_obj), devo_meta[,1] ), 15] 
	type_profiles <- group_rowmeans(this_obj, this_types)
	
	if (length(grep("Myeloid", colnames(type_profiles))) > 0 ) {
		# Store myeloid-specific genes
		macro <- grep("Myeloid", colnames(type_profiles))
		macro_genes <- apply(type_profiles, 1, function(x) {x[macro]/max(x[-macro]) > 10 & x[macro] == max(x) & max(x) > 0.05})
		macro_genes <- as.character(genes_meta[macro_genes,3])
		macrophage_gene_table[,this_name] <- as.numeric(rownames(macrophage_gene_table) %in% macro_genes);
	} else {
		next;
	}
}

# Add expression data
macrophage_gene_expr = rowMeans(pseudobulks_blair)
macrophage_gene_expr <- macrophage_gene_expr[match(rownames(macrophage_gene_table), names(macrophage_gene_expr))]
liver_mac_gene_expr <- all_profiles[,"Liver"]
liver_mac_gene_expr <- liver_mac_gene_expr[match(rownames(macrophage_gene_table), names(liver_mac_gene_expr))]
Adrenal_mac_gene_expr <- all_profiles[,"Adrenal"]
Adrenal_mac_gene_expr <- Adrenal_mac_gene_expr[match(rownames(macrophage_gene_table), names(Adrenal_mac_gene_expr))]


data <- data.frame(macrophage_gene_table, expression_BLAIR=macrophage_gene_expr, expression_Liver=liver_mac_gene_expr, expression_Adrenal=Adrenal_mac_gene_expr)
write.table(data, "Reference_fetal_tissue-specific.csv", row.names=T, col.names=T, sep=",", quote=TRUE)

# Plots with UpSetR
require(UpSetR)

# Only upto 5 sets it looks like? - subset to most similar based on correlations.
upset(data)

data <- data[, c("Liver", "Spleen", "Adrenal", "Placenta", "Pancreas", "expression_BLAIR", "expression_Liver", "expression_Adrenal")]
upset( data, 
	attribute.plots=list(gridrows=100, ncols=1,
		plots=list(list(plot=scatter_plot, x="expression_Liver", y="expression_BLAIR", queries=T),
			     list(plot=scatter_plot, x="expression_Adrenal", y="expression_BLAIR", queries=T))
		), 
	order.by="freq",
	queries=list(list(query = intersects, params=list("Liver"), active=T), list(query = intersects, params=list("Adrenal"))) )


tmp <- data[data[,"Liver"] == 1 & rowSums(data[,1:5]) == 1,]
tmp[order(tmp[,"expression_Liver"], decreasing=TRUE),]