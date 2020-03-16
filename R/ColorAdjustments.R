#### Darken: For darkening colors ####
#' Darkens input colors by a set amount
#'
#' @description A wrapper for the darken function of the colorspace package.
#' @param colors the color(s) input. Can be a list of colors, for example, /code{dittoColors()}.
#' @param percent.change # between 0 and 1. the percentage to darken by. Defaults to 0.25 if not given.
#' @param relative TRUE/FALSE. Whether the percentage should be a relative change versus an absolute one. Default = TRUE.
#' @return Return a darkened version of the color in hexadecimal color form (="#RRGGBB" in base 16)
#' @examples
#' Darken("blue") #"blue" = "#0000FF"
#' #Output: "#0000BF"
#' Darken(dittoColors()[1:8]) #Works for multiple color inputs as well.
#'
#' @author Daniel Bunis
#' @export
Darken <- function(colors, percent.change = 0.25, relative = TRUE) {
    colorspace::darken(
        colors, amount = percent.change, space = "HLS", fixup = TRUE,
        method = ifelse(relative, "relative", "absolute"))
}

#### Lighten: For lightening colors ####
#' Lightens input colors by a set amount
#'
#' @description A wrapper for the lighten function of the colorspace package.
#' @param colors the color(s) input. Can be a list of colors, for example, /code{dittoColors()}.
#' @param percent.change # between 0 and 1. the percentage to darken by. Defaults to 0.25 if not given.
#' @param relative TRUE/FALSE. Whether the percentage should be a relative change versus an absolute one. Default = TRUE.
#' @return Return a lighter version of the color in hexadecimal color form (="#RRGGBB" in base 16)
#' @examples
#' Lighten("blue") #"blue" = "#0000FF"
#' #Output: "#4040FF"
#' Lighten(dittoColors()[1:8]) #Works for multiple color inputs as well.
#'
#' @author Daniel Bunis
#' @export
Lighten <- function(colors, percent.change = 0.25, relative = TRUE) {
    colorspace::lighten(
        colors, amount = percent.change, space = "HLS", fixup = TRUE,
        method = ifelse(relative, "relative", "absolute"))
}

#### Simulate: For simulating what a plot would look like as seen by a colorblind person ####
#' Simulates what a colorblind person would see for any dittoSeq plot!
#'
#' @description Essentially a wrapper function for colorspace's deutan(), protan(), and tritan() functions. This function will output any dittoSeq plot as it might look to an individual with one of the common forms of colorblindness: deutanopia/deutanomaly, the most common, is when the cones mainly responsible for red vision are defective. Protanopia/protanomaly is when the cones mainly responsible for green vision are defective. In tritanopia/tritanomaly, the defective cones are responsible for blue vision. Note: there are more severe color deficiencies that are even more rare. Unfortunately, for these types of color vision deficiency, only non-color methods, like lettering or shapes, will do much to help.
#' @param type The type of colorblindness that you want to simulate for. Options: "deutan", "protan", "tritan". Anything else, and you will get an error.
#' @param plot.function The plotting function that you want to use/simulate. not quoted. and make sure to remove the () that R will try to add.
#' @param ... other paramters that can be given to dittoSeq plotting functions, including color.panel, used in exactly the same way they are used for those functions. (contrary to the look of this documentation, color.panel will still default to dittoColors() when not provided.)
#' @param color.panel,min.color,max.color The set of colors to be used.
#' @return Outputs a dittoSeq plot with the color.panel / min.color & max.color updated as it might look to a colorblind individual.
#'
#' Note: Does not currently adjust dittoHeatmap.
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' Simulate("deutan", dittoDimPlot, object=myRNA, var="clustering", size = 2)
#' Simulate("protan", dittoDimPlot, myRNA, "clustering", size = 2)
#' Simulate("tritan", dittoDimPlot, myRNA, "clustering", size = 2)
#'
#' @author Daniel Bunis
#' @importFrom colorspace deutan protan tritan
#' @export
Simulate <- function(
    type = c("deutan","protan","tritan"),
    plot.function, ..., color.panel = dittoColors(),
    min.color = "#F0E442", max.color = "#0072B2") {

    type <- match.arg(type)

    color.adj <- eval(expr = parse(text = paste0(
        "colorspace::",type,"(color.panel)")))
    min.adj <- eval(expr = parse(text = paste0(
        "colorspace::",type,"(min.color)")))
    max.adj <- eval(expr = parse(text = paste0(
        "colorspace::",type,"(max.color)")))

    #Make the plot!
    plot.function(color.panel = color.adj, min.color = min.adj,
                  max.color = max.adj, ... )
}
