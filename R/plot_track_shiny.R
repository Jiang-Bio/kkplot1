

#' Title
#'
#' @param bwfile
#'
#' @return
#' @export
#'
#' @examples
plot_track_shiny <- function(bwfile){

  library(shiny)
  library(JBrowseR)
  ui <- fluidPage(
    titlePanel("JBrowseR"),
    # this adds to the browser to the UI, and specifies the output ID in the server
    JBrowseROutput("browserOutput")
  )

  server <- function(input, output, session) {
    # create the necessary JB2 assembly configuration
    assembly <- assembly(
      "/wrk/data/genome/yz_genome_data/aragenome/ara.fa",
      bgzip = F
    )

    # create configuration for a JB2 GFF FeatureTrack
    annotations_track <- track_feature(
      "/wrk/data/genome/yz_genome_data/aragenome/athaliana.gff3.gz",
      assembly
    )
    wiggle_track <- track_wiggle(bwfile,assembly,bgzip=F)
    # create the tracks array to pass to browser
    tracks <- tracks(
      annotations_track,
      wiggle_track
    )

    # set up the default session for the browser
    default_session <- default_session(
      assembly,
      c(annotations_track)
    )

    theme <- theme("#5da8a3", "#333")

    # link the UI with the browser widget
    output$browserOutput <- renderJBrowseR(
      JBrowseR(
        "View",
        assembly = assembly,
        tracks = tracks,
        location = "Chr3:1..50000",
        defaultSession = default_session,
        theme = theme
      )
    )
  }

  shinyApp(ui, server)


}
