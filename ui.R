#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(data.table)
library(DT)
library(plotly)
library(EnsDb.Hsapiens.v86)
ensdb <- EnsDb.Hsapiens.v86
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(shinycssloaders)
#library(biomaRt)


#dataTrack_nano1 <- DataTrack(range = "/home/smaegol/analyses/ONT/RNA/wgs2/NA12878-DirectRNA.pass.dedup.NoU.fastq.hg38.minimap2.sorted.bam.bw")

#polyA_deduplicated_summary_hg38 <- readRDS("/home/smaegol/analyses/ONT/RNA/wgs2/R/polyA_deduplicated_summary_hg38.rds")
#polyA_deduplicated_all_hg38_dt <- readRDS("/home/smaegol/analyses/ONT/RNA/wgs2/R/polyA_deduplicated_all_hg38_dt.rds")
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  
  # Application title
  titlePanel("polyA_tails"),
  
  
  #   sidebarLayout(
  #     
  #   sidebarPanel(
  #  
  #      DT::dataTableOutput('x1',width=800),
  #     plotlyOutput('x2', height = 500, width=800),
  #     checkboxInput("show_alignment","Show reads alignment",FALSE),
  #     checkboxInput("show_only_selected_transcripts","Show Only transcripts with alignments",FALSE),width = 800),
  #  
  # mainPanel(
  #  tabsetPanel(type="tabs",position = "right",
  #    tabPanel("browser",plotOutput('x3', height = 500, width=700)),
  #    tabPanel("browser2",plotOutput('x4', height = 500, width=700)))
  # )
  # )
  # )
  
  # # sidebarPanel(
  fluidRow(
    column(5,
           DT::dataTableOutput('table')),
    column(5,plotlyOutput('boxplot', height = 500, width=1000) %>% withSpinner(type = 4))
  ),
  fluidRow(
    column(8,
           tabsetPanel(type="tabs",position = "left",
                      # tabPanel("ENCODE",plotOutput('encode', height = 500, width=1000)),
                       tabPanel("distribution",plotlyOutput("distribution")),
                       tabPanel("alignment",plotOutput('alignment', height = 500, width=1000)),
                       tabPanel("detailed_summary",verbatimTextOutput("detailed_summary")))),
    column(3,checkboxInput("show_alignment","Show reads alignment",FALSE),
           checkboxInput("show_only_selected_transcripts","Show Only transcripts with alignments",FALSE),
           checkboxInput("histogram_distribution","Show distribution as histogram",FALSE),
           actionButton("reset_click","Reset to all isoforms per gene")))
  
  # ),
  # mainPanel(
  #  tabsetPanel(type="tabs",position = "right",
  #    tabPanel("isoforms",plotlyOutput('x2', height = 500, width=1000)),
  #    tabPanel("browser",plotOutput('x3', height = 500, width=1000),
  # plotOutput('x4', height = 500, width=1000)))
  # #)
  
  
)