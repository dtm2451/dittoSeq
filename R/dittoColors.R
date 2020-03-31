#' Extracts the dittoSeq default colors
#' @author Daniel Bunis
#' @param reps Integer which sets how many times the original set of colors should be repeated
#' @param get.names Logical, whether only the names of the default dittoSeq color panel should be returned instead
#' @return A string vector with length = 24.
#' @description Creates a string vector of 40 unique colors, in hexadecimal form, repeated 100 times.
#' Or, if \code{get.names} is set to \code{TRUE}, outputs the names of the colors which can be helpful as reference when adjusting how colors get used.
#'
#' These colors are a modification of the protanope and deuteranope friendly colors from Wong, B. Nature Methods, 2011.
#'
#' Truly, only the first 1-7 are maximally (red-green) color-blindness friendly, but the lightened and darkened versions (plus grey) in slots 8-40 still work releatively well at extending their utility further.
#' Note that past 40, the colors simply repeat in order to most easily allow dittoSeq visualizations to handle situations requiring even more colors.
#'
#' The colors are:
#'
#' 1-7 = Suggested color panel from Wong, B. Nature Methods, 2011, minus black
#' \itemize{
#' \item 1- orange = "#E69F00"
#' \item 2- skyBlue = "#56B4E9"
#' \item 3- bluishGreen = "#009E73"
#' \item 4- yellow = "#F0E442"
#' \item 5- blue = "#0072B2"
#' \item 6- vermillion = "#D55E00"
#' \item 7- reddishPurple = "#CC79A7"
#' }
#'
#' 8 = gray40
#'
#' 9-16 = 25\% darker versions of colors 1-8
#'
#' 17-24 = 25\% lighter versions of colors 1-8
#'
#' 25-32 = 40\% lighter versions of colors 1-8
#'
#' 33-40 = 40\% darker versions of colors 1-8
#'
#' @export
#' @examples
#' dittoColors()
#'
#' #To retrieve names:
#' dittoColors(get.names = TRUE)
dittoColors <- function(reps = 100, get.names = FALSE) {
    if (get.names) {
        cols <- c(
            "orange", "skyBlue", "bluishGreen", "yellow",
            "blue", "vermillion", "reddishPurple", "gray40",
            "orange-25%", "skyBlue-25%", "bluishGreen-25%", "yellow-25%",
            "blue-25%", "vermillion-25%", "reddishPurple-25%", "gray40-25%",
            "orange+25%", "skyBlue+25%", "bluishGreen+25%", "yellow+25%",
            "blue+25%", "vermillion+25%", "reddishPurple+25%", "gray40+25%",
            "orange+40%", "skyBlue+40%", "bluishGreen+40%", "yellow+40%",
            "blue+40%", "vermillion+40%", "reddishPurple+40%", "gray40+40%",
            "orange-40%", "skyBlue-40%", "bluishGreen-40%", "yellow-40%",
            "blue-40%", "vermillion-40%", "reddishPurple-40%", "gray40-40%")
    } else {
        cols <- c(  # Colors as of dittoSeq-v0.2.10
            "#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#666666",
            "#AD7700", "#1C91D4", "#007756", "#D5C711",
            "#005685", "#A04700", "#B14380", "#4D4D4D",
            "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
            "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
            "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
            "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3",
            "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
            "#00446B", "#803800", "#8D3666", "#3D3D3D")
    }
    return(rep(cols, reps))
}
