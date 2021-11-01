#### Packages

library(shiny)

library(Seurat)
blair_obj <- readRDS("data/Clustered_Blair.rds")
all_genes <- sort(unique(rownames(blair_obj)))
all_features <- colnames(blair_obj@meta.data)[2:ncol(blair_obj@meta.data)]

not_gene_token <- "(None)"

#### UI
ui <- fluidPage(
	# App title ----
	titlePanel("Keller Lab iPSC Derived Macrophages"),

	# Sidebar Layout with input and output definitions ----
	sidebarLayout (

		# Sidebar panel for inputs ----
		sidebarPanel(
		

			# Input: Select Gene to display
			selectInput(inputId = "gene",
	                  label = "Choose a gene:",
                  	choices = c(not_gene_token, all_genes), 
				selected = not_gene_token),

			# Input: Select Gene to display
			selectInput(inputId = "feature",
	                  label = paste('Choose a feature (Set gene to:"', not_gene_token, '")', sep =""),
                  	choices = all_features, 
				selected = "zonation_score"),

			# Input: Set resolution of image to download
			sliderInput("res",
				label = "Download Image Resolution (dpi)",
                       	min = 50, 
				max = 500, 
				step = 50,
				value = 300),

			# Download Button
			downloadButton("downloadPlot", "Download Plot")
			
		),
	
		# Main panel for displaying outputs ----
		mainPanel(
			# Formatted Text for explanation & Citation information
			p(strong("Creator:"), tags$a(href="https://uk.linkedin.com/in/tallulah-andrews-60169122", "Tallulah Andrews")),
			p(strong("Data by:"), 
				tags$a(href="https://ca.linkedin.com/in/blair-gage-6053b744", "Blair Gage"), "&",
				tags$a(href="https://www.uhn.ca/Research/Research_Institutes/McEwen_Stem_Cell_Institute/Research/Pages/Gordon_Keller.aspx#lab", "Keller Lab")
			),
			p(strong("Code:"), tags$a(href="https://github.com/tallulandrews/KellerLab_Macrophages", "https://github.com/tallulandrews/KellerLab_Macrophages")),


			# Output: SpatialFeaturePlot
			plotOutput(outputId = "umapFeaturePlot")

		)
	)
)

#### Server
server <- function(input, output) {

	# Seurat Spatial Feature Plot of the specified gene and slice.
	plotImage <- function() {
		toplot <- input$gene
		if (input$gene == not_gene_token) {
			toplot <- input$feature
		}
		print(suppressWarnings(Seurat::FeaturePlot( blair_obj, features=toplot)))
	}


	# renderPlot = react to input & output is a plot
	output$umapFeaturePlot <- renderPlot({
		plotImage()
	})

	output$downloadPlot <- downloadHandler(
		filename = paste("UMAP", input$gene, input$slice, "plot.png", sep="_"),
		content = function(file) {
		png(file, width=6, height=6, units="in", res=input$res)
		plotImage()
		dev.off()
    })  

}

#### shinyApp Call
shinyApp(ui = ui, server = server)