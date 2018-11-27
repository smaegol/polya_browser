#######################################################################################
###                                                                                 ###
###     Copyright (C) 2018  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################

library(shiny)
library(dplyr)
library(ggplot2)
library(data.table)
library(DT)
library(plotly)
library(EnsDb.Mmusculus.v79)
ensdb <- EnsDb.Mmusculus.v79
library(Gviz)
library(rtracklayer)
library(GenomicRanges)

library(BSgenome.Mmusculus.UCSC.mm10)
library(shinycssloaders)

#Set GVIZ and ensembldb compatibility options
seqlevelsStyle(ensdb) <- "UCSC"
options(ucscChromosomeNames=FALSE)

#paths to files to read (mapping and coverage)
bigwig_file_wt<-BigWigFile("wt_0h_LPS_all.fastq.mm18.mapped.sorted.bam.bw")
bigwig_file_mut<-BigWigFile("AC_0h_LPS_all.fastq.mm18.mapped.sorted.bam.bw")
bam_file_wt <- "wt_0h_LPS_all.fastq.mm18.mapped.sorted.bam"
bam_file_mut <- "AC_0h_LPS_all.fastq.mm18.mapped.sorted.bam"

axis_track <- GenomeAxisTrack() #initialize genome axis track for Gviz

selected_transcript = "Gapdh" #default transcript to show

get_location <- function(location2,show_only_selected_transcripts = FALSE) {
  locationz<-GenomicRanges::reduce(location2)
  end(locationz)[1]<-end(locationz)[length(locationz)]
  locationz<-locationz[1]
  location_from=start(locationz)
  location_to = end(locationz)
  
  length_location = location_to - location_from
  location_from = location_from - 0.1*length_location
  location_to = location_to + 0.1*length_location
  end(locationz)[1] <- location_to
  start(locationz)[1] <- location_from 
  return(locationz)
}

#read alignment track for given location
get_alignment_track <- function(locationz) {
  location_from=start(locationz)
  location_to = end(locationz)
  align_Track <- AlignmentsTrack(bam_file_wt,isPaired = FALSE,allow.nonnarrowing=TRUE,start=location_from,end=location_to)
  return(align_Track)
}

#get gene track for given location
get_gene_track <- function(location,show_only_selected_transcripts = FALSE,selected_isoform=NA,found_isoforms=c()) {
  location_from <- start(location)
  location_to <- end(location)
  #print(location)
  if (show_only_selected_transcripts == TRUE) {
    gene_track <- GeneRegionTrack(location,col.line = NULL, col = NULL)
  }
  else {
    gene_track <- GeneRegionTrack(getGeneRegionTrackForGviz(ensdb,chromosome = as.character(seqnames(location))[1],start=location_from,end=location_to),thinBoxFeature=c("utr", "ncRNA", "utr3", "utr5", "miRNA", "lincRNA"))
  }
  # print(gene_track)
  # print("gen")
  feature(gene_track)[transcript(gene_track) %in% found_isoforms]<-"found"
  if (selected_isoform!='NA') {
    feature(gene_track)[transcript(gene_track)==selected_isoform]<-"good"
    interestcolor <- list("found"="green", "good"="red")
  }
  else {
    interestcolor <- list("found"="green")
  }
  # print(feature(gene_track))
  #print(interestcolor)
  Gviz::displayPars(gene_track) <- interestcolor
  return(gene_track)
}

plot_genome <- function(locationz,gene_track,align_track=NA,show_only_selected_transcripts = FALSE) {
  
  #print(locationz)
  # Load data, at a reasonable level of detail for width(location)
  n <- min(width(locationz), 1000)
  location_from <- start(locationz)
  location_to <- end(locationz)
  coverage_wt <- rtracklayer::summary(
    bigwig_file_wt, locationz, n, "max")[[1]]
  
  coverage_mut <- rtracklayer::summary(
    bigwig_file_mut, locationz, n, "max")[[1]]
  
  data_track_wt <- DataTrack(
    coverage_wt, data=coverage_wt$score,  
    name="wt_0h_LPS coverage", type="l", col="#000000", legend=FALSE)
  
  data_track_mut<- DataTrack(
    coverage_mut, data=coverage_mut$score,  
    name="AC_0h_LPS coverage", type="l", col="#00AA00", legend=FALSE)
  
  
  ideoTrack <- IdeogramTrack(genome = "mm10", chromosome = as.character(seqnames(locationz))[1])
  
  sTrack <- SequenceTrack(Mmusculus)
  
  tracks_list <- list(ideoTrack,sTrack,gene_track,data_track_wt,data_track_mut)
  
  if (!is.na(align_track)) {
    tracks_list <- c(tracks_list,align_track)
  }
  
  plotTracks(
    tracks_list,
    chromosome=as.character(seqnames(locationz)), 
    from=location_from, to=location_to, transcriptAnnotation = "transcript")
}

plot_genome2 <- function(location2,show_only_selected_transcripts = FALSE) {
  locationz<-GenomicRanges::reduce(location2)
  end(locationz)[1]<-end(locationz)[length(locationz)]
  locationz<-locationz[1]
  location_from=start(locationz)
  location_to = end(locationz)
  
  length_location = location_to - location_from
  #location_from = location_from - 0.1*length_location
  #location_to = location_to + 0.1*length_location
  end(locationz)[1] <- location_to
  start(locationz)[1] <- location_from
  #print(locationz)
  # Load data, at a reasonable level of detail for width(location)
  n <- min(width(locationz), 1000)
  
  if (show_only_selected_transcripts == TRUE) {
    gene_track <- GeneRegionTrack(location2,col.line = NULL, col = NULL)
  }
  else {
    gene_track <- GeneRegionTrack(getGeneRegionTrackForGviz(ensdb,chromosome = as.character(seqnames(locationz))[1],start=location_from,end=location_to))
  }
  align_Track <- AlignmentsTrack(bam_file,isPaired = FALSE,allow.nonnarrowing=TRUE)
  
  ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = as.character(seqnames(locationz))[1])
  
  #biomTrack <- BiomartGeneRegionTrack(genome = "hg38",chromosome = as.character(seqnames(locationz))[1], start = start(locationz), end = end(locationz), name = "ENSEMBL",biomart = mart)
  sTrack <- SequenceTrack(Hsapiens)
  plotTracks(
    list(sTrack,gene_track,align_Track),
    chromosome=as.character(seqnames(locationz)), 
    from=location_from, to=location_to, transcriptAnnotation = "transcript")
}


plot_genome3 <- function(location2,show_only_selected_transcripts = FALSE,selected_isoform = '',found_isoforms = c()) {
  locationz<-GenomicRanges::reduce(location2)
  end(locationz)[1]<-end(locationz)[length(locationz)]
  locationz<-locationz[1]
  location_from=start(locationz)
  location_to = end(locationz)
  
  length_location = location_to - location_from
  location_from = location_from - 0.1*length_location
  location_to = location_to + 0.1*length_location
  end(locationz)[1] <- location_to
  start(locationz)[1] <- location_from
  #print(locationz)
  # Load data, at a reasonable level of detail for width(location)
  n <- min(width(locationz), 1000)
  
  
  ck <- AlignmentsTrack(bam_file,isPaired = FALSE)
  
  
  if (selected_isoform=='NA') {

    feature(gene_track)[transcript(gene_track)==selected_isoform]<-"good"
    interestcolor <- list("found"="green", "good"="red")

  }
  else {
    interestcolor <- list("found"="green")
  }
 
  displayPars(gene_track) <- interestcolor
  
  
  ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = as.character(seqnames(locationz))[1])
  
  sTrack <- SequenceTrack(Hsapiens)
  plotTracks(
    list(ideoTrack,sTrack,gene_track),
    chromosome=as.character(seqnames(locationz)), 
    from=location_from, to=location_to, transcriptAnnotation = "transcript")
}


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  selected_isoform <- reactiveValues(isoform_name='NA')
  
  output$table = DT::renderDataTable(macrophages_merged_summary, server = TRUE, selection='single')
  output$detailed_summary = renderPrint({
    selected_row <- input$table_rows_selected
    selected_transcript = macrophages_merged_summary[selected_row,]$transcript
    macrophages_merged_stat_wilcoxon_summary %>% filter(transcript==selected_transcript) %>% glimpse()
    })
  
  #a scatterplot with certain points highlighted
  output$boxplot = renderPlotly({
    selected_row <- input$table_rows_selected
    
    data_transcript = data_transcript()
    isoforms_plot <- ggplot(data_transcript,aes(x=ensembl_transcript_id,y=polya_length)) + geom_boxplot()  + ggtitle(selected_transcript)
    if (nrow(data_transcript)<10000) {
      isoforms_plot <- isoforms_plot + geom_jitter(aes(color=sample_name))
    }
    ggplotly(isoforms_plot)
    
  })
  
  
  output$distribution = renderPlotly({
    selected_row <- input$table_rows_selected
    selected_transcript = macrophages_merged_summary[selected_row,]$transcript
    data_transcript = data_transcript()
    
    if (selected_isoform$isoform_name!='NA') {
      if (input$histogram_distribution==FALSE) {
        distribution_plot <- ggplot(data_transcript %>% filter(ensembl_transcript_id==selected_isoform$isoform_name),aes(x=polya_length,colour=group)) + geom_density()  + ggtitle(selected_transcript)
      }
      else {
        distribution_plot <- ggplot(data_transcript %>% filter(ensembl_transcript_id==selected_isoform$isoform_name),aes(x=polya_length,colour=group)) + geom_histogram() + facet_wrap(. ~ group)  + ggtitle(selected_transcript)
      }
          }
    else {
      if (input$histogram_distribution==FALSE) {
        distribution_plot <- ggplot(data_transcript,aes(x=polya_length,colour=group)) + geom_density()  + ggtitle(selected_transcript)
      }
      else {
        distribution_plot <- ggplot(data_transcript,aes(x=polya_length,colour=group)) + geom_histogram() + facet_wrap(. ~ group)  + ggtitle(selected_transcript)
      }
    }
    ggplotly(distribution_plot)
    
  })
  
  
  
  
  data_transcript <- reactive({
    selected_row <- input$table_rows_selected 
    selected_transcript = macrophages_merged_summary[selected_row,]$transcript
    data_transcript = macrophages_merged_dt[transcript == selected_transcript]
    data_transcript
  })
  
  
  
  # get genomic ranges of transcripts of interest
  gene_model <- reactive({
    getGeneRegionTrackForGviz(ensdb, filter = AnnotationFilter(~ tx_id %in% transcript_isoforms() ))
  },label="gene_model")
  
  
  observeEvent(input$reset_click, {
    selected_isoform$isoform_name <- "NA"
  })  
  
  
  clicked_isoform <- reactive({
    event_data("plotly_click")
  })
  
  selected_isoform_event <- observeEvent(clicked_isoform(), {
    s = clicked_isoform()
    x_plot = as.numeric(s$x)
    if (length(x_plot)!=0) {
      transcript_number = round(x_plot[1])
      data_transcript = data_transcript()
      transcript_isoforms <- levels(as.factor(data_transcript$ensembl_transcript_id))
      selected_isoform$isoform_name = transcript_isoforms[transcript_number]
    }
    else {
      selected_isoform$isoform_name = "NA"
    }
  },label="selected_iso")
  
  locationz <- reactive({
    get_location(gene_model())
  },label="locationz")
  
  gene_track <- reactive({
    locationz = locationz()
    selected_transcript_name = selected_isoform$isoform_name
    #print(selected_transcript_name)
    get_gene_track(locationz,found_isoforms = transcript_isoforms(),selected_isoform = selected_transcript_name)
    
  },label="gene_track")
  
  align_track <- reactive({
    locationz = locationz()
    get_alignment_track(locationz)
    
  },label="align_track")
  
  # get vector or names of analyzed transcript isoforms
  transcript_isoforms <- reactive({
    data_transcript = data_transcript()
    unique(data_transcript$ensembl_transcript_id)
  },label="transcript_isoforms")
  
  output$encode = renderPlot({
    
    locationz = locationz()
    gene_track = gene_track()
    
        plot_genome(locationz,gene_track)
    
    })
  
  output$alignment = renderPlot({
    
    if (input$show_alignment==TRUE) {
      locationz = get_location(gene_model())
      gene_track = gene_track()
      align_track = align_track()
      print(gene_track)
      plot_genome(locationz,gene_track,align_track = align_track(),input$show_only_selected_transcripts)
      
    }
    else {
      print("Alignment plot disabled")
    }
  })
  
  output$clicks = renderPrint({
    s = event_data("plotly_click")
    x_plot = as.numeric(s$x)
    transcript_number = round(x_plot[1])
    data_transcript = data_transcript()
    transcript_isoforms <- levels(as.factor(data_transcript$ensembl_transcript_id))
    transcript_name = transcript_isoforms[transcript_number]
    if (length(s) == 0) {
      "Click on plot"
    }
    else {
      cat(paste0("Selected:",transcript_name," \n\n"))
      
      as.list(s)
      
    }
    
  })
  
}
