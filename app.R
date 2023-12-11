#### Packages
# Works but uses too much memory #

library(shiny)
library(scales)
library(Seurat)
library(ggplot2)
human_obj <- readRDS("data/ForShiny_Human_obj_minimal.rds")
mouse_obj <- readRDS("data/ForShiny_Mouse_obj_minimal.rds")
all_hgenes <- sort(unique(rownames(human_obj$expr)))
all_mgenes <- sort(unique(rownames(mouse_obj$expr)))

not_gene_token <- "(None)"

#tmp <- colnames(human_obj$metadata); tmp[5] <- "clusters"; colnames(human_obj$metadata)<-tmp
#tmp <- colnames(mouse_obj$metadata); tmp[5] <- "clusters"; colnames(mouse_obj$metadata)<-tmp

all_mfeatures <- colnames(mouse_obj$metadata)[2:ncol(mouse_obj$metadata)]
all_hfeatures <- colnames(human_obj$metadata)[2:ncol(human_obj$metadata)]
all_features <- unique(c(all_hfeatures, all_mfeatures))
all_features <- all_features[all_features != "orig.ident"]

get_colours <- function(identities) {
        require(scales)
        if (class(identities) != "factor") {
                identities <- factor(identities)
        }
        identities <- levels(identities)
        my_color_palette <- hue_pal()(length(identities))
	  names(my_color_palette) <- identities
        return(my_color_palette)
}

#ReDef Seurat Functions

UMAP_plot <- function(obj, values, value_name, is.human=TRUE, is.discrete=FALSE) {
	coords <- obj$metadata[,c("UMAP1", "UMAP2")]
	#coords <- data.frame(UMAP1=as.numeric(coords[,1]), UMAP2=as.numeric(coords[,2]))
	coords[,value_name] <- values
	colnames(coords) <- c("UMAP1", "UMAP2", value_name)

	if (is.human) {
		mytitle <- paste(value_name, "(human)")
	} else {
		mytitle <- paste(value_name, "(mouse)")
	}

	if ( !is.discrete ){
		# Continuous
		thisplot <- ggplot(coords, aes(UMAP1, UMAP2))+geom_point(aes(color=.data[[value_name]]), size=0.5)+
				ggtitle(mytitle)+theme_classic()+scale_colour_gradient2(low="cadetblue1", mid="deepskyblue", high="blue3", 
				midpoint=max(0.01,mean(values,na.rm=TRUE)))

	} else {
		# Discrete
		this_palette <- get_colours(unique(values))
		
		coords[,3] <- factor(coords[,3])
		thisplot <- ggplot(coords, aes(UMAP1, UMAP2))+ geom_point(aes(color=.data[[value_name]]), size=0.5) + 
				theme_classic() +ggtitle(mytitle)+ scale_color_manual(values=this_palette)

	}
	print(thisplot+theme(plot.title = element_text(face=2, size=15, hjust=0.5)))
}




#### UI
ui <- fluidPage(
	# App title ----
	titlePanel("Keller Lab iPSC Derived Macrophages"),

	# Sidebar Layout with input and output definitions ----
	sidebarLayout (

		# Sidebar panel for inputs ----
		sidebarPanel(
		

			# Input: Select Human Gene to display
			selectInput(inputId = "hgene",
	                  label = "Choose a Human gene:",
                  	choices = c(not_gene_token, all_hgenes), 
				selected = not_gene_token),

			# Input: Select Mouse Gene to display
			selectInput(inputId = "mgene",
	                  label = "Choose a Mouse gene:",
                  	choices = c(not_gene_token, all_mgenes), 
				selected = not_gene_token),

			# Plot M gene on H data
			checkboxInput(inputId="m_on_h", label="Plot Mouse expression on Human cells", value = FALSE, width = NULL),

			# Input: Select Gene to display
			selectInput(inputId = "feature",
	                  label = paste('Choose a feature (Set gene to:"', not_gene_token, '")', sep =""),
                  	choices = all_features, 
				selected = "clusters"),

			# Input: Set resolution of image to download
			sliderInput("res",
				label = "Download Image Resolution (dpi)",
                       	min = 50, 
				max = 500, 
				step = 50,
				value = 300),

			# Download Button
			downloadButton("downloadMPlot", "Download Mouse Plot"),
			# Download Button
			downloadButton("downloadHPlot", "Download Human Plot")

			
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
			plotOutput(outputId = "umapFeaturePlot1"),
			plotOutput(outputId = "umapFeaturePlot2")

		)
	)
)





#### Server
server <- function(input, output, session) {

	# Seurat Spatial Feature Plot of the specified gene and slice.
	plotHImage <- function() {
		is.discrete=FALSE
		if(input$m_on_h) {
			toplot <- input$mgene
			values <- mouse_obj$expr[toplot,]
			values <- values[match(rownames(human_obj$metadata), names(values))]
			value_name <- input$mgene
			print(suppressWarnings(UMAP_plot( human_obj, values=values, value_name=value_name, is.human=TRUE)))
		} else {
			if (input$hgene == not_gene_token) {
				toplot <- input$feature
				if (toplot %in% c("clusters", "annotation", "Last Meal")) {
					is.discrete = TRUE
				}
				values <- human_obj$metadata[,toplot]
				value_name <- input$feature
			} else {
				toplot <- input$hgene
				values <- human_obj$expr[toplot,]
				value_name <- input$hgene
			}
		}
		print(suppressWarnings(UMAP_plot( human_obj, values=values, value_name=value_name, is.human=TRUE, is.discrete=is.discrete)))
	}
	plotMImage <- function() {
		is.discrete=FALSE
		if (input$mgene == not_gene_token) {
			toplot <- input$feature
			if (toplot %in% c("clusters", "annotation", "Last Meal")) {
				is.discrete = TRUE
			}
			values <- mouse_obj$metadata[,toplot]
			value_name <- input$feature
		} else {
			toplot <- input$mgene
			values <- mouse_obj$expr[toplot,]
			value_name <- input$mgene
		}
		print(suppressWarnings(UMAP_plot( mouse_obj, values=values, value_name=value_name, is.human=FALSE, is.discrete=is.discrete)))
	}



	# renderPlot = react to input & output is a plot
	output$umapFeaturePlot1 <- renderPlot({
		plotHImage()
	})
	output$umapFeaturePlot2 <- renderPlot({
		plotMImage()
	})


	output$downloadHPlot <- downloadHandler(
		filename = paste("UMAP", input$hgene, "human_plot.png", sep="_"),
		content = function(file) {
		png(file, width=6, height=6, units="in", res=input$res)
		plotHImage()
		dev.off()
    })  
	output$downloadMPlot <- downloadHandler(
		filename = paste("UMAP", input$mgene, "mouse_plot.png", sep="_"),
		content = function(file) {
		png(file, width=6, height=6, units="in", res=input$res)
		plotMImage()
		dev.off()
    })  

	# Stop when the browser tab is closed.
	session$onSessionEnded(stopApp)

}

#### shinyApp Call
shinyApp(ui = ui, server = server)
 

}

#### shinyApp Call
shinyApp(ui = ui, server = server)
