library(shiny)
library(rgdal)
library(rgeos)
library(adehabitatHR)
library(ggmap)
library(ggplot2)

# Thresholds of IUCN
EOO <- c(100, 2000, 20000)
AOO <- c(10, 500, 2000)


shinyUI <- fluidPage(

  titlePanel("A Theshold Calculator for IUCN's Threat Levels using eBird Maps"),
  
  sidebarPanel(
    
    helpText(h4("Download a file for a single species from eBird. Upload your data in the same format (zip) that you got it from eBird.")),
    fileInput('zipfile', h3('Upload',accept = c('.zip','.csv'))),
    
    helpText(h4("To use recent data and ignore historical data.")),
    sliderInput('startYear', 'Ignore data till Year', min=2000, max=2016,
                value=2000, step=1, round=0),
    
    helpText(h4("Shorter the distance traveled, more the precision.")),
    sliderInput('maxDistance', 'Ignore lists with distance (in km) travelled is greater than...', min=0, max=100,
                value=50, step=1, round=0),
    
    helpText(h4("Short duration lists will inherently be confined to small area.")),
    sliderInput('maxDuration', 'Ignore lists with duration of list (in hours) is greater than...', min=0, max=24,
                value=8, step=0.5, round=0),
    
    helpText(h4("Used only for AOO. IUCN suggests one way to calculate AOO is by summing constant sized-squares where the taxon is present.One Degree is roughly 111km.")),
    selectInput('squareSize', 'Size of Square', c("One Degree"=1, "Half Degree"=0.5, "Quarter Degree"=0.25, "8th Degree"=0.125, "16th Degree"=0.0625), selected=0.125)
    
  ),
  
  
  mainPanel(
    headerPanel(textOutput ("species")),
    tabsetPanel(
      tabPanel("MCP", 
               fluidRow(
                 column (3, 
                         textOutput ("mcparea")        
                 ),
                 br(),
                 column (1, 
                         plotOutput('MCP')))),
      tabPanel("Area of Occupancy", 
               fluidRow(
                 column (3, 
                         textOutput ("aooarea")        
                 ),
                 br(),
                 column (1, 
                         plotOutput('AOO')))),
      tabPanel("About", 
               br(), h1("Minimum Convex Polygon Calculator"), 
               br(), p("eBird has data on the range of species."), 
               br(), p("IUCN uses MCP to calculate Extend of Occurrence (EOO) of each species."), 
               br(), p("This app helps calculate MCP from eBird data output."),
               br(), p("IUCN provides the option to count constant sized squares of arbitrary size to calculate Are of Occupancy (AOO) of each species."), 
               br(), p("This app helps calculate AOO from eBird data output by defining the square size to be used."),
               br(), a("Understanding eBird Data Quality", href = "http://help.ebird.org/customer/portal/articles/1055676-understanding-the-ebird-review-and-data-quality-process"),
               br(), br(), a("Understanding Extent of Occurance and Area of Occupancy", href = "http://www.iucnredlist.org/technical-documents/categories-and-criteria/2001-categories-criteria"),
               br(), br(), p("Fore any comments/queries. Contact Praveen J: paintedstork@gmail.com")
               
      )
    )
  )
)


options(shiny.maxRequestSize=30*1024^2) 


shinyServer <- function(input, output) {
  
  output$species <- renderText ( {
    if(is.null(input$zipfile)) return(NULL)
    
    ebd <- read.delim(unz(input$zipfile$datapath,gsub('zip','txt',input$zipfile$name)), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    ebd$COMMON.NAME[1]
  })
  
  output$mcparea <- renderText( { 
    if(is.null(input$zipfile)) return(NULL)
    
    ebd_raw <- read.delim(unz(input$zipfile$datapath,gsub('zip','txt',input$zipfile$name)), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    
    ebd <- processEbdFiles (ebd_raw, 
                            input$startYear,
                            input$maxDistance,
                            input$maxDuration)
    
    
    ebd   <- subset(ebd, select = c("SCIENTIFIC.NAME", "LATITUDE", "LONGITUDE"))
    sp::coordinates(ebd) <- ~LONGITUDE+LATITUDE
    
    sp::proj4string(ebd) <- CRS("+init=epsg:4326")
    
    mcp_species_area <- mcp(spTransform(ebd, CRS("+init=epsg:3857")), percent=100, unout="km2")
    
    iucn_msg <- ifelse (mcp_species_area$area < EOO[1], "Potentially Critically Endangered", 
                  ifelse (mcp_species_area$area < EOO[2], "Potentially Endangered",
                    ifelse (mcp_species_area$area < EOO[3], "Potentially Vulnerable", "")))
    paste ("MCP= ",prettyNum(as.integer(mcp_species_area$area), big.mark=",",scientific=FALSE), "sq.km.", iucn_msg)
  })
  
  
  output$MCP <-   renderPlot ( {
    
    if(is.null(input$zipfile)) return(NULL)
    
    ebd_raw <- read.delim(unz(input$zipfile$datapath,gsub('zip','txt',input$zipfile$name)), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    ebd_file_name <- 'ebd_wispet_relFeb-2017'
    
    unzip(paste('..\\data\\',ebd_file_name,'.zip',sep=''))
#    ebd_raw <- read.delim(paste(ebd_file_name,'.txt',sep=''), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")


    ebd <- processEbdFiles (ebd_raw, 
                            input$startYear,
                            input$maxDistance,
                            input$maxDuration)
                            
    
    #Subset to useful data and convert to sp
    ebd   <- subset(ebd, select = c("SCIENTIFIC.NAME", "LATITUDE", "LONGITUDE"))
    sp::coordinates(ebd) <- ~LONGITUDE+LATITUDE
    
    #Calculate MCP
    mcp_species <- mcp(ebd, percent=100, unout="km2")
    

    #Pick bounding square and convert to form usable in google
    bounds <- bbox(mcp_species)
    sbbox <- make_bbox(lon = c(bounds[1], bounds[3]), lat = c(bounds[2], bounds[4]), f = 2)
    
#    sbbox <- unlist(geocode('Bengaluru, Karnataka')) + c(0,.02)
#    sbbox <- unlist(geocode('San Francisco, California')) + c(0,.02)
    
    gmap <- get_map (location=sbbox, maptype = "terrain", source = "google", col="color")

    ggmap(gmap) +
      geom_polygon(data=fortify(mcp_species), aes(x=long, y=lat, group=group), color="red", fill = "red", alpha=0.5) +
      coord_map()  +
      theme(legend.position = "none")
  }, height = 700, width = 700)

  output$aooarea <- renderText( { 
    if(is.null(input$zipfile)) return(NULL)
    
    ebd_raw <- read.delim(unz(input$zipfile$datapath,gsub('zip','txt',input$zipfile$name)), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    ebd_file_name <- 'ebd_barfly1_relFeb-2017'
    
    unzip(paste('..\\data\\',ebd_file_name,'.zip',sep=''))
#    ebd_raw <- read.delim(paste(ebd_file_name,'.txt',sep=''), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    
    ebd <- processEbdFiles (ebd_raw, 
                            input$startYear,
                            input$maxDistance,
                            input$maxDuration)
    
    
    squareSize <- as.numeric (input$squareSize)
    
    #Making coordinates coarse
    ebd <- within (ebd, 
                   LONGITUDE <-  squareSize * as.integer(LONGITUDE/squareSize + 0.5))
    ebd <- within (ebd, 
                   LATITUDE  <-  squareSize * as.integer(LATITUDE/squareSize  + 0.5))
    
    ebd   <- subset(ebd, select = c("LATITUDE", "LONGITUDE"))
    
    #Make it one entry per square
    ebd <- ebd[!duplicated(ebd[c("LONGITUDE","LATITUDE")]),]

    aoo_species_area <- 111 * 111 * as.numeric (input$squareSize) * as.numeric (input$squareSize) * nrow(ebd)
    
    iucn_msg <- ifelse (aoo_species_area < AOO[1], "Potentially Critically Endangered", 
                        ifelse (aoo_species_area < AOO[2], "Potentially Endangered",
                                ifelse (aoo_species_area < AOO[3], "Potentially Vulnerable", "")))
    
    paste ("AOO=",prettyNum(aoo_species_area, big.mark=",",scientific=FALSE), "sq.km.",iucn_msg)
  })
  
  
  output$AOO <-   renderPlot ( {
    
    if(is.null(input$zipfile)) return(NULL)
    
    ebd_raw <- read.delim(unz(input$zipfile$datapath,gsub('zip','txt',input$zipfile$name)), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    ebd_file_name <- 'ebd_barfly1_relFeb-2017'
    
    unzip(paste('..\\data\\',ebd_file_name,'.zip',sep=''))
#    ebd_raw <- read.delim(paste(ebd_file_name,'.txt',sep=''), na.strings = c("NA", "", "null"), as.is=TRUE, quote="")
    
    ebd <- processEbdFiles (ebd_raw, 
                            input$startYear,
                            input$maxDistance,
                            input$maxDuration)
            
    # Moving into sp domain for calculating the bounding polygon by doing mcp first
    sp_ebd <- ebd
    sp::coordinates(sp_ebd) <- ~LONGITUDE+LATITUDE

    #Calculate MCP
    mcp_species <- mcp(sp_ebd, percent=100, unout="km2")
    
    #Pick bounding square and convert to form usable in google
    bounds <- bbox(mcp_species)
    sbbox <- make_bbox(lon = c(bounds[1], bounds[3]), lat = c(bounds[2], bounds[4]), f = 2)
    
        
    squareSize <- as.numeric (input$squareSize)
    
    #Making coordinates coarse
    ebd <- within (ebd, 
                   LONGITUDE <-  squareSize * as.integer(LONGITUDE/squareSize + 0.5))
    ebd <- within (ebd, 
                   LATITUDE  <-  squareSize * as.integer(LATITUDE/squareSize  + 0.5))
    
    ebd   <- subset(ebd, select = c("LATITUDE", "LONGITUDE"))
    
    #Make it one entry per square
    ebd <- ebd[!duplicated(ebd[c("LONGITUDE","LATITUDE")]),]

#    sbbox <- unlist(geocode('Bengaluru, Karnataka')) + c(0,.02)
#    sbbox <- unlist(geocode('San Francisco, California')) + c(0,.02)
    
    gmap <- get_map (location=sbbox, maptype = "terrain", source = "google", col="color")
    ebd$color <- 1 
    ggmap(gmap) + geom_tile(data=ebd, aes(x=LONGITUDE,y=LATITUDE, color = "black", fill = "blue", alpha=0.1)) +
      coord_map()  +
      theme(legend.position = "none")
  }, height = 700, width = 700)
  
} 

processEbdFiles <- function (ebd, startYear, maxDistance, maxDuration) {
  
      #Add unique list identifier for removing duplicates
      ebd <- within (ebd, UNIQUE_SAMPLING_ID <-  ifelse(is.na(GROUP.IDENTIFIER),SAMPLING.EVENT.IDENTIFIER,GROUP.IDENTIFIER))
      
      #If subspecies, copy subspecies common name
      ebd <- within (ebd, COMMON.NAME <-  ifelse(CATEGORY=='issf',SUBSPECIES.COMMON.NAME,COMMON.NAME))

      #Remove entries from shared lists
      ebd <- ebd[!duplicated(ebd[c("UNIQUE_SAMPLING_ID","COMMON.NAME")]),]

      #Add Year
      ebd <- within (ebd, 
                     YEAR <-  as.numeric(format(as.Date(OBSERVATION.DATE),"%Y")))

      #Filter out historical lists 
      if(!is.null(startYear))
      {
        ebd <- ebd [which(ebd$YEAR >= startYear), ]
      }

      #Filter out long traveling lists
      if(!is.null(maxDistance))
      {
        ebd <- ebd [which(ebd$EFFORT.DISTANCE.KM <= maxDistance), ]
      }

      #Filter out long duration lists 
      if(!is.null(maxDuration))
      {
        ebd <- ebd [which(ebd$DURATION.MINUTES <= 60 * maxDuration), ]
      }

      return (ebd)
}

shinyApp(ui = shinyUI, server = shinyServer)

