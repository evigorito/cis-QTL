library(rmarkdown)

#' Render report
#'
#' @param script full path to R script to render
#' @return saves a pdf
#' @export
#' main()

main <- function(script){
    render(script)
}

main(snakemake@input[['script']])
