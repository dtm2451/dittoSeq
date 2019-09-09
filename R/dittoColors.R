#' Extracts the DittoSeq default colors
#' @author Daniel Bunis
#' @param get.names Logical, whether only the names of the default DittoSeq color panel should be returned instead
#' @description Creates a colors vector of 24 colors in hexadecimal form,
#' a modification of the protanope and deuteranope friendly colors from Wong, B. Nature Methods, 2011.
#'
#' The colors are:
#'
#' 1-7 = Suggested color panel from Wong, B. Nature Methods, 2011, minus black
#'
#' 8 = gray40
#'
#' 9-16 = 25\% darker versions of colors 1-8
#'
#' 17-24 = 25\% lighter versions of colors 1-8
#'
#' @export
#' @examples
#' dittoColors()
#'
#' #To retrieve names:
#' dittoColors(get.names = TRUE)
dittoColors <- function(get.names = FALSE) {
    if (get.names) {
        return(c(
            "orange", "skyBlue", "bluishGreen", "yellow", "blue",
            "vermillion", "reddishPurple", "gray40", "orange-25%", "skyBlue-25%",
            "bluishGreen-25%", "yellow-25%", "blue-25%", "vermillion-25%", "reddishPurple-25%",
            "gray40-25%", "orange+25%", "skyBlue+25%", "bluishGreen+25%", "yellow+25%",
            "blue+25%", "vermillion+25%", "reddishPurple+25%", "gray40+25%"))
    } else {
        return(c(  # DittoSeq-v0.2.10 Colors
            "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
            "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
            "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
            "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"))
    }
}
