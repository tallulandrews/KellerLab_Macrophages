Contamination_genes <- c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1")
Prolif_genes <- c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4")
RBC_genes <- c("HBB", "HBA1", "HBA2", "HBD")

Endo_genes_SN_SC <- c("CTGF", "FCGR2B", "S100A13", "FCN2", "FCN3", "LYVE1", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
				"F8", "CALCRL", "SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF",
"ENG", "PECAM1", "RAMP3", "INMT", "DNASE1L3", "LIFR", "TIMP3", "C7"
				)
Endo_genes <- c("EPCAM", "PContamination_genes <- c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1")
Prolif_genes <- c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4")

Hepato_Portal <- c("AGT", "ALB", "ALDOB", "APOA1", "CYP1A2", "CYP2A6", "CYP2A7", "CYP2D6", "FGA", "FGB", "GLS2", "HAL", "HAMP", "SDS")
Hepato_Central <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "CES1", "CYP2E1", "CYP3A4", "CYP3A5", "DCXR", "GSTA1", "OAT")

RBC_genes <- c("HBB", "HBA1", "HBA2", "HBD")

Endo_genes_SN_SC <- c("CTGF", "FCGR2B", "S100A13", "FCN2", "FCN3", "LYVE1", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
				"F8", "CALCRL", "SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF",
"ENG", "PECAM1", "RAMP3", "INMT", "DNASE1L3", "LIFR", "TIMP3", "C7"
				)
Endo_genes <- c("PECAM1","CLEC4M", "CD34", "COL1A2", "ENG", "STAB2", "CLDN5")
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

label_genes <- function(genes, label) {names(genes) <- rep(label, length(genes)); return(genes)}
Final_Macrophage_genes <- c(
	label_genes(c("AREG", "CCL3", "CD83", "CXCL2", "CXCL3", "IL1B", "NAMPT", "PLAUR", "SRGN", "THBS1"), "Activated"),
	label_genes(c("C1QA", "C1QB", "C1QC", "CD163", "CD5L", "CTSB", "FTL", "HMOX1", "MARCO", "SEPP1", "SLC40A1", "VCAM1"), "Kupffer"),
	label_genes(c("ACP5", "CSTB", "FABP5", "LGMN", "PLD3", "PSAP"), "LAM-like"),
	label_genes(c("CD74", "HLA-DPA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1"), "MHCII"),
	label_genes(c("FCN1", "LYZ", "MNDA", "S100A12", "S100A4", "S100A6", "S100A8", "S100A9", "VCAN"), "Monocyte"),
	label_genes(c("FCER1G", "FOLR2", "MS4A7", "TIMD4", "TIMP1"), "Resident"))

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
		"Synapse", "Kupffer", "Kupffer", "Monocyte", "Resident", "Mono-Act", 
		"Monocyte", "MHCII", "Activated", "MHCII", "LAM-like", 
		"Debris", "Monocyte", "LSEC-Doublet", "Kupffer")

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

outname="Redo_Blair"
set.seed(101)
#blair <- Read10X("C:/Users/tandrews/Documents/UHNSonya/RawData/BlairMacs/Blair_Macs_cellranger_out/filtered_feature_bc_matrix")
blair <- Read10X("C:/Users/tandrews/Documents/UHNSonya/RawData/BlairMacs/STARSolo_mapping/BlairMacs_Human_STARsolo_filtered_gzs")
blair <- CreateSeuratObject(counts = blair, project = "HumanGenome", min.cells = 3, min.features = 200)
#Normalize
orig_norm_factor = 1500
blair <- NormalizeData(blair, normalizeation.method="LogNormalize", scale.factor=orig_norm_factor)

humanmap <- blair
humanmap[["percent.mt"]] <- PercentageFeatureSet(humanmap, pattern = "^MT-")
VlnPlot(humanmap, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(humanmap, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(humanmap, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2plot1 + plot2
humanmap <- subset(humanmap, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)
humanmap <- FindVariableFeatures(humanmap, selection.method = "vst", nfeatures = 2000)
grep("^MT-", VariableFeatures(humanmap))
humanmap <- ScaleData(humanmap, vars.to.regress = "percent.mt")

blair <- humanmap
#Scale -> UMAP
set.seed(189)
blair <- ScaleData(blair, features=rownames(blair))
blair <- RunPCA(blair, features=VariableFeatures(obj), dims=1:20)
blair <- RunUMAP(blair, dims=1:20)

## Cluster then match clusters ##
blair <- Seurat::FindNeighbors(blair, dims=1:10)
blair <- Seurat::FindClusters(blair, res=0.5)

blair <- readRDS(paste(outname, "Clustered_redo.rds", sep="_"))
DimPlot(blair)

##Marker Gene Table##
blair_markers <- FindAllMarkers(blair, logfc.threshold=-Inf)
blair_markers <- blair_markers[,blair_markers$p_val_adj < 0.001]
write.table(blair_markers, "TableS2_our_human_mac_cluster_markers.csv", sep=",")


DotPlot(blair, feature=c(Final_Macrophage_genes, label_genes(Endo_genes, "Endo"), label_genes(Hepato_Portal, "P-Hepato"), label_genes(Hepato_Central, "C-Hepato"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png("Figure_3A_UMAP.png", width=6*2.1/2, height=2*(2+2.4)/2, units="in", res=300)
DimPlot(blair)
dev.off()

png("Figure_3B_Dotplot.png", width=6*2.1, height=2*2, units="in", res=300)
DotPlot(blair, feature=c(Final_Macrophage_genes, label_genes(Endo_genes, "Endo"), "APOE", "IPMK", "COX5A", "H3F3B"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png("Figure_S3B_FeaturePlots.png", width=6*2.1, height=2*2.4, units="in", res=300)
FeaturePlot(blair, c("PECAM1", "DNASE1L3", "FCGR2B", "STAB2", "OIT3", "GPR182", "RAMP3", "F8"), ncol=4)
dev.off()

png("Figure_S3A_FeaturePlots.png", width=6*2.1, height=2*2.4, units="in", res=300)
FeaturePlot(blair, c("C1QC", "CD163", "FOLR2", "LYVE1", "CD38", "FCGR3A", "IL10", "VCAN"), ncol=4)
dev.off()

png("Figure_S3C_FeaturePlots.png", width=6*2.1, height=2*2.4/2, units="in", res=300)
FeaturePlot(blair, c("SERPINA1", "ALB", "CYP2D6", "ASS1"), ncol=4)
dev.off()

png("Figure_3C_FeaturePlots.png", width=6*2.1, height=2*2.4, units="in", res=300)
FeaturePlot(blair, c("PTPRC", "MARCO", "CD5L", "TIMD4", "CD68", "S100A6", "FCN1", "MNDA"), ncol=4)
dev.off()

# Plan 2 -> map cells to nearest neighbours using spearman correlations.
require(Hmisc)
# Use Scaled data for reference.
# Subset to key markers from dot plot
match_features = Final_Macrophage_genes[Final_Macrophage_genes %in% intersect(rownames(obj), rownames(blair))]
mat <- obj@assays$RNA@scale.data[rownames(obj) %in% match_features,]
mat_blair <- blair@assays$RNA@data[rownames(blair) %in% match_features,]
require(qlcMatrix)

cor_mat <- corSparse(mat, mat_blair)
match <- apply(cor_mat, 2, function(x){  obj@meta.data$Subcluster_Manual[ which(x == max(x)) ]  })
blair@meta.data$matched_annotation <- match

DimPlot(blair, group.by="matched_annotation")
table(match)

# Plan 3 - Markers of clusters & match across.

DE <- FindAllMarkers(blair)
write.table(DE, "Supplementary_Table_HumanDE.csv", sep=",")

## ---- Pseudobulk correlations --- ##

# --- Blair vs 20 Liver Map --- #

# Match Genes
obj <- obj[rownames(obj) %in% rownames(blair),]
blair <- blair[rownames(blair) %in% rownames(obj),]

endo_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Endo_harmony_Subcluster.rds")
cluster_col = "Coarse_clusters"
manual_cluster_anno <- c("cvLSEC", "cvLSEC", "cvLSEC", "ppLSEC", "cvLSEC", "cvEndo", 
	"ppLSEC", "interLSEC", "HepContam", "Arterial")
endo_obj@meta.data$Subcluster_Manual <- manual_cluster_anno[endo_obj@meta.data[,cluster_col]]
endo_obj <- endo_obj[rownames(endo_obj) %in% rownames(blair)]

#All Macrophage Markers
Idents(obj) <- (obj@meta.data$Subcluster_Manual)
mac_marks <- FindAllMarkers(obj, logfc.threshold=0.25, pos=TRUE)
mac_marks <- mac_marks[mac_marks[,2] > 0.3 & mac_marks$cluster %in% names(Final_Macrophage_genes) & !grepl("^MT-", mac_marks$gene) & !grepl("^RP[SL]", mac_marks$gene),]
write.table(mac_marks, "prototype_TableS3.csv", sep=",")

#All Endo Markers
Idents(endo_obj) <- (endo_obj@meta.data$Subcluster_Manual)
endo_marks <- FindAllMarkers(endo_obj, logfc.threshold=0.25, pos=TRUE)
endo_marks <- endo_marks[endo_marks[,2] > 0.3 & endo_marks$cluster %in% c("cvLSEC", "ppLSEC") & !grepl("^MT-", endo_marks$gene) & !grepl("^RP[SL]", endo_marks$gene),]
write.table(endo_marks, "prototype_TableS3_endo.csv", sep=",")



require(Hmisc)

profile_all_blair <- group_rowmeans(blair@assays$RNA@data, blair@meta.data$seurat_clusters)
profile_all_20Livers <- group_rowmeans(obj@assays$RNA@data, obj@meta.data$Subcluster_Manual)
profile_all_endo <- group_rowmeans(endo_obj@assays$RNA@data, endo_obj@meta.data$Subcluster_Manual)

common_genes <- intersect(rownames(profile_all_blair), intersect(rownames(profile_all_endo), rownames(profile_all_20Livers)))
key_genes <- unique(c(Final_Macrophage_genes, Endo_genes, mac_marks$gene, endo_marks$gene))
key_genes <- key_genes[key_genes %in% common_genes]

profile_20Livers <- profile_all_20Livers[match(key_genes, rownames(profile_all_20Livers)),]
profile_endo <- profile_all_endo[match(key_genes, rownames(profile_all_endo)),]
profile_blair <- profile_all_blair[match(key_genes, rownames(profile_all_blair)),]

profile_sim <- rcorr(profile_blair, cbind(profile_endo, profile_20Livers), type="spearman")
heatdata <- t(profile_sim$r[c("Kupffer", "LAM-like", "MHCII",  "Activated", "Monocyte", "cvLSEC", "ppLSEC"),c("0","1","2","3","4","5")])
heatdata[heatdata<0] <-0
require(gplots)
png(paste("Figure3D_20LiverMap_corr_heatmap_fancy.png", sep="_"), width=6, height=8, units="in", res=300)
heatmap.2(round(heatdata, digits=2), 
		trace="none", symbreaks=FALSE, col=colorRampPalette(c("black","firebrick", "darkorange","yellow"))(30), mar=c(8,4),
		key.title="", key.xlab="Spearman Cor", Colv=FALSE, Rowv=FALSE)
dev.off()



## ----- Blair vs Devo_map ----- ##

# source: https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")
require(Seurat)
# First lets get the gene expression profiles of all macrophages from the Devo map.
dir = "C:/Users/tandrews/Documents/UHNSonya/ExternalData/DevelopmentalCellAtlas"

#Which cells are macrophages?
devo_meta <- readRDS(paste(dir, "/df_cell.RDS", sep=""))

table(devo_meta[grepl("Myeloid", devo_meta[,15]),15])

objs <- Sys.glob(paste(dir, "/*count.RDS", sep=""))
genes_meta <- readRDS(paste(dir, "df_gene.RDS", sep="/"))

#objs <- objs[c(1,4,5,6,7,8,10,11,12,13)]

all_profiles <- c();
organ <- c();
best_genes <- c();
gene_2_organ_table <- c();
# Gather Myeloid data
for (this_i in 1:length(objs)) {

	#if (this_i == 3) {next;}

	#this_counts <- this_obj[,colnames(this_obj) %in% devo_meta[grepl("Myeloid", devo_meta[,15]), 1] ]
	#if (ncol(this_counts) == 0) {this_counts <- this_obj[,colnames(this_obj) %in% devo_meta[grepl("Microglia", devo_meta[,15]), 1] ]}

	if (this_i == 3) { 

		type_profiles <- readRDS("Cerebrum_type_profiles_profile.rds")
		this_name <- "Cerebrum"
	} else {
		# select genes >> myeloid than any other cell type
		# mean expression across all cell-types
		this_obj <- readRDS(objs[this_i])
		this_name <- unique(devo_meta[ devo_meta[,1] %in% colnames(this_obj) , "Organ"])		
		this_types <- devo_meta[match( colnames(this_obj), devo_meta[,1] ), 15] 
		type_profiles <- group_rowmeans(this_obj, this_types)
	}	

	print(colnames(type_profiles))
	macro <- which(grepl("Microglia", colnames(type_profiles)) | grepl("Myeloid", colnames(type_profiles)))
	myeloid_profile <- type_profiles[,macro]


	#print(objs[this_i])
	#print(table(this_types))
	#print(ncol(this_counts))
	
	if (length(macro) > 0) {
		#profile <- Matrix::rowMeans(this_counts)
		all_profiles <- cbind(all_profiles, myeloid_profile)

		# Store myeloid-specific genes
		#macro <- grep("Myeloid", colnames(type_profiles))
		macro_genes <- apply(type_profiles, 1, function(x) {x[macro]/max(x[-macro]) > 10 & x[macro] == max(x) & max(x) > 0.05})
		macro_genes_symbol <- as.character(genes_meta[macro_genes,3])
		best_genes <- c(as.character(best_genes), macro_genes_symbol)
		gene_2_organ_table <- rbind(gene_2_organ_table, cbind(macro_genes_symbol, as.character(rep(this_name, length(macro_genes_symbol))) ))
	
		organ <- c(organ, rep(this_name, length(macro)))
	} else {
		print(this_name)
		print(colnames(type_profiles))
		next;
	}
}

write.table(gene_2_organ_table, "prototype_SupplementaryTable4.csv", sep=",")

rownames(all_profiles) <- genes_meta[,3]
colnames(all_profiles) <- organ

saveRDS(list(profiles = all_profiles,best_genes=best_genes, macro_genes=macro_genes), "full_refAtlas_Myeloid_profiles.rds")

Devo_Data <- readRDS("full_refAtlas_Myeloid_profiles.rds")
all_profiles <- Devo_Data$profiles
macro_genes <- Devo_Data$macro_genes
best_genes <- Devo_Data$best_genes

# unique genes
macrophage_genes <- unique(best_genes)
write.table(macrophage_genes, file="Supplementary_Table_devo_macrophage_genes.csv", sep=",")

key_profiles <- all_profiles[rownames(all_profiles) %in% macrophage_genes,]

#fix weird duplicates
key_profiles <- key_profiles[rowSums(key_profiles) > 0,]

require(gplots)
png(paste(outname, "reference_atlas_key_genes_heatmap.png", sep="_"), width=6, height=8, units="in", res=150)
heatmap.2(key_profiles, dendrogram="none", trace="none", scale="row", mar=c(6,4))
dev.off()

# match blair data
common_genes <- intersect(rownames(key_profiles), rownames(profile_all_blair))

key_genes_blair <- profile_all_blair[match(common_genes, rownames(profile_all_blair)),]
key_profiles <- key_profiles[match(common_genes, rownames(key_profiles)),]

identical(rownames(key_genes_blair), rownames(key_profiles))

# combine
full_key_genes <- cbind(key_profiles, key_genes_blair)
full_key_genes <- apply(full_key_genes, 2, scale) # This is much better than straight normalization.
rownames(full_key_genes) <- rownames(key_profiles)

# make fancy table
# gene x organ 1/0s + mean expression in blair & reference
fancy_table <- matrix(0, ncol=ncol(key_profiles), nrow=length(common_genes))
rownames(fancy_table) <- common_genes
colnames(fancy_table) <- sort(colnames(key_profiles))
for ( i in 1:nrow(fancy_table) ) {
	g <- rownames(fancy_table)[i];
	organs <- gene_2_organ_table[gene_2_organ_table[,1] == g,2]
	fancy_table[g,] <- colnames(fancy_table) %in% organs
}
full_fancy_table <- cbind(fancy_table, key_profiles, key_genes_blair)
write.table(full_fancy_table, file="Supplementary_table_devo_mac_genes.csv", sep=",")

# heatmap
heatmap.2(full_key_genes, dendrogram="none", trace="none", scale="none")

png(paste(outname, "vs_Fetal_profiles_full_heatmap.png", sep="_"), width=8, height=10, units="in", res=150)
heatmap.2(full_key_genes, dendrogram="none", trace="none", scale="row")
dev.off()

# correlations
require(Seurat)
png("Figure3E_Fetal_profile_full_correlations_fancy.png", width=10, height=8, units="in", res=300)
heatmap.2( t( cor(full_key_genes, method="spearman")[ 1:ncol(key_profiles), ( ncol(key_profiles)+1 ):ncol( full_key_genes )]) , cexCol=0.2+1/log10(7),
		trace="none", mar=c(8,4), symbreaks=FALSE, col=colorRampPalette(c("black","firebrick", "darkorange","yellow"))(30), key.title="", key.xlab="Spearman correlation")
dev.off()
#col=colorRampPalette(c("white", "black"))(20)

### Extra Stuff Blair vs Devo Map ###

### Liver-specific macrophage-markers

# make table of which genes are markers in which tissues
macrophage_gene_table <- matrix(0, nrow=length(macrophage_genes), ncol=length(organ))
colnames(macrophage_gene_table) <- organ
rownames(macrophage_gene_table) <- macrophage_genes

for (this_i in 1:length(objs)) {

	if (this_i == 3) {
		type_profiles <- readRDS("Cerebrum_type_profiles_profile.rds")
		this_name <- "Cerebrum"

	} else {

		this_obj <- readRDS(objs[this_i])
		this_name <- unique(devo_meta[ devo_meta[,1] %in% colnames(this_obj) , "Organ"])

		# select genes >> myeloid than any other cell type
		# mean expression across all cell-types
		this_types <- devo_meta[match( colnames(this_obj), devo_meta[,1] ), 15] 
		type_profiles <- group_rowmeans(this_obj, this_types)
	}

	macro <- which(grepl("Microglia", colnames(type_profiles)) | grepl("Myeloid", colnames(type_profiles)))
	
	if (length(macro) > 0 ) {
		# Store myeloid-specific genes
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
Cerebellum_mac_gene_expr <- all_profiles[,"Cerebellum"]
Cerebellum_mac_gene_expr <- Cerebellum_mac_gene_expr[match(rownames(macrophage_gene_table), names(Cerebellum_mac_gene_expr))]



data <- data.frame(macrophage_gene_table, expression_BLAIR=macrophage_gene_expr, 
				expression_Liver=liver_mac_gene_expr, 
				expression_Cerebellum=Cerebellum_mac_gene_expr,
				expression_Adrenal=Adrenal_mac_gene_expr)
write.table(data, paste(outname, "Reference_fetal_tissue-specific_full.csv", sep="_"), row.names=T, col.names=T, sep=",", quote=TRUE)

# Plots with UpSetR
require(UpSetR)

# Only upto 5 sets it looks like? - subset to most similar based on correlations.

data <- data[, c("Liver", "Spleen", "Adrenal", "Placenta", "Pancreas", "expression_BLAIR", "expression_Liver", "expression_Adrenal", "expression_Cerebellum")]
png(paste(outname, "Tissue_specific_full_dataset.png", sep="_"), width=10, height=8, units="in", res=300)
upset( data, 
	attribute.plots=list(gridrows=100, ncols=1,
		plots=list(list(plot=scatter_plot, x="expression_Liver", y="expression_BLAIR", queries=T),
			     list(plot=scatter_plot, x="expression_Adrenal", y="expression_BLAIR", queries=T))
		), 
	order.by="freq",
	queries=list(list(query = intersects, params=list("Liver"), active=T), 
			 list(query = intersects, params=list("Adrenal"))) )

dev.off()


tmp <- data[data[,"Liver"] == 1 & rowSums(data[,1:5]) == 1,]
tmp[order(tmp[,"expression_Liver"], decreasing=TRUE),]

### Blair Mouse ###
set.seed(101)
mousemap <-Read10X(data.dir = "C:/Users/tandrews/Documents/UHNSonya/RawData/BlairMacs/STARSolo_mapping/BlairMacs_Mouse_STARsolo_filtered")
mousemap <- CreateSeuratObject(counts = mousemap, project = "MouseGenome", min.cells = 3, min.features = 200)
mousemap[["percent.mt"]] <- PercentageFeatureSet(mousemap, pattern = "mt-")
mousemap <- subset(mousemap, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)
mousemap <- NormalizeData(mousemap, normalization.method = "LogNormalize", scale.factor = 10000)
mousemap <- FindVariableFeatures(mousemap, selection.method = "vst", nfeatures = 2000)
mousemap <- ScaleData(mousemap, vars.to.regress = "percent.mt", features=rownames(mousemap))
#mousemap <- ScaleData(mousemap)
mousemap <- RunPCA(mousemap, features = VariableFeatures(object = mousemap))
mousemap <- FindNeighbors(mousemap, dims = 1:10)
mousemap <- FindClusters(mousemap, resolution = 0.8)
mousemap <- RunUMAP(mousemap, dims = 1:15)
DimPlot(mousemap, reduction = "umap")
saveRDS(mousemap, "new_mouse_map_obj.rds")
#mousemap <- readRDS("new_mouse_map_obj.rds")

##Marker Gene Table##
Idents(mousemap) <- mousemap@meta.data$seurat_clusters
mouse_markers <- FindAllMarkers(mousemap, logfc.threshold=-Inf)
mouse_markers <- mouse_markers[mouse_markers$p_val_adj < 0.001,]
write.table(mouse_markers, "TableS2_our_mouse_eaten_cluster_markers.csv", sep=",")

#Mouse Clusters
#COLOUR
colour <-c("Degraded Cells"="#ef1a1c","Neutrophils"="#377eb8","Kupffer Cells"="#4daf4a","Erythrocytes"="#984ea3","Endothelial Cells"="#ff7f00","LSECs"="#ffff33","Hepatocytes"="#f781bf","NA"="#969696")

new.cluster.ids <- c("1","1","1","1","2","3","4","5","6","7") 
names(new.cluster.ids) <- levels(mousemap)
mousemap <- RenameIdents(mousemap, new.cluster.ids)
new.cluster.ids <- c("Degraded Cells","Neutrophils","Kupffer Cells","Erythrocytes","Endothelial Cells","LSECs","Hepatocytes") 
names(new.cluster.ids) <- levels(mousemap)
mousemap <- RenameIdents(mousemap, new.cluster.ids)
#png("Figure7A_MouseUMAP.png", width=6, height=6, units="in", res=300)
png("Figure_7A_MouseUMAP_dim.png", width=6*2.1/2, height=2*(2+2.4)/2, units="in", res=300)
#DimPlot(mousemap, reduction = "umap")+NoLegend()+annotate("text", x = 0, y = 1.5, label = "Degraded Cells", fontface=2)+annotate("text", x = -1, y = 5, label = "Endothelial Cells", fontface=2)+annotate("text", x = 18.5, y = 2.5, label = "Hepatocytes", fontface=2)+annotate("text", x = 7.5, y = -9.5, label = "LSECs", fontface=2)+annotate("text", x = -1, y = -12.5, label = "Erythrocytes", fontface=2)+annotate("text", x = 0, y = -8, label = "Kupffer Cells", fontface=2)+annotate("text", x = -2.5, y = -7, label = "Neutrophils", fontface=2)+scale_color_manual(values=colour)
DimPlot(mousemap, reduction = "umap")+NoLegend() +annotate("text", x = -3, y = -0.5, label = "Degraded Cells", fontface=2)+annotate("text", x = -2, y = 4, label = "Endothelial Cells", fontface=2)+annotate("text", x = -10.5, y = -3.5, label = "Hepatocytes", fontface=2)+annotate("text", x = -10, y = 10.5, label = "LSECs", fontface=2)+annotate("text", x = 2, y = 16, label = "Erythrocytes", fontface=2)+annotate("text", x = -1, y = 11, label = "Kupffer Cells", fontface=2)+annotate("text", x = 0, y = 6.5, label = "Neutrophils", fontface=2)+scale_color_manual(values=colour)
dev.off()

#Determining which macrophages have eaten which cell type.
blair@meta.data$"hPSC-Macrophage phagocytosed Last Meal Analysis"<- mousemap@active.ident[match(rownames(blair@meta.data),rownames(mousemap@meta.data))] 
#png("Figure7C_MouseOntoHuman.png", width=6, height=6, units="in", res=300)
png("Figure_7A_MouseOntoHuman_dim.png", width=6*2.1/2, height=2*(2+2.4)/2, units="in", res=300)
DimPlot(blair, group.by = "hPSC-Macrophage phagocytosed Last Meal Analysis")+scale_color_manual(values=colour)+ggtitle("Mouse Cell Phagocytosed")
dev.off()

my_barplot <- function(tab) {
	tab <- tab/sum(tab)*100;
	coords <- barplot(tab, las=2, col=unlist(colour), ylim=c(0,100), ylab="Cells (%)")
	text(coords, tab, paste(round(tab, digits=1), "%"), pos=3)
}

par(mar=c(8,4,1,1))
my_barplot( table(blair@meta.data$"hPSC-Macrophage phagocytosed Last Meal Analysis") )
png("SupplFigure_TotalEaten.png", width=6, height=6, units="in", res=300)
par(mar=c(8,4,1,1))
my_barplot( table(blair@meta.data$"hPSC-Macrophage phagocytosed Last Meal Analysis"[blair@meta.data$seurat_clusters %in% 0:3]) )
dev.off()

MouseDE <- FindAllMarkers(mousemap)
write.table(MouseDE, "Supplementary_MouseDE.csv", sep=",")


label_genes <- function(genes, label) {names(genes) <- rep(label, length(genes)); return(genes)}
neutro <- label_genes(c("Cd52", "S100a4", "S100a11", "S100a8", "S100a9"), "Neutrophil") # , "Chil1"
rbc <- label_genes(c("Alas2", "Hbb-bs", "Hba-a1", "Hbb-bt"), "Erythro") # "Hbq1b", "Car2", "Hemgn", 
kupffer <- label_genes(c("C1qa", "C1qc", "Cd68", "Lilra5", "Acp5"), "Kupffer")
degraded <- label_genes(c("Gm5644", "Gm11517", "Rpl39l", "Rpl19-ps12"), "Unkn")
endo <- label_genes(c("Col1a2", "Myh10",  "Pcdh17", "Hmgb1", "Plk2" ), "Endo") # "Dpp10",
lsec <- label_genes(c("Stab1", "Stab2", "Clec4g", "Pecam1", "Sparc"), "LSEC") #Igfbp7
hepato <- label_genes(c("Alb", "Aldob", "Aldh8a1", "Sult1d1", "Adh4"), "Hepato")

tmp <- c(degraded, neutro, kupffer, rbc, endo, lsec, hepato); names(tmp) <- NULL
png("Figure7B_MouseMarkers.png", width=6, height=6.5, units="in", res=300)
DotPlot(mousemap, features=tmp)+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## Shiny Objects ##
blair_shiny <- blair
blair_shiny@meta.data <- blair_shiny@meta.data[,-6]
mouse_shiny <- mousemap
mouse_shiny@meta.data$annotation <- Idents(mouse_shiny)
mouse_shiny@meta.data <- mouse_shiny@meta.data[,-5]

mouse_shiny@reductions <- mouse_shiny@reductions[2]
blair_shiny@reductions <- blair_shiny@reductions[2]

write.table(blair_shiny@meta.data, "ForUpload_Human_metadata.csv", sep=",")
write.table(mouse_shiny@meta.data, "ForUpload_Mouse_metadata.csv", sep=",")
saveRDS(blair_shiny, "ForShiny_Human_obj.rds")
saveRDS(mouse_shiny, "ForShiny_Mouse_obj.rds")
