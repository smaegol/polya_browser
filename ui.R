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


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  
  # Application title
  titlePanel("polyA_tails macrophages"),
  
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
  
  
  
)