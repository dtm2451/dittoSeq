#### Darken: For darkening colors ####
#' Darkens input colors by a set amount
#'
#' @description A wrapper for the darken function of the colorspace package.
#' @param colors the color(s) input. Can be a list of colors, for example, MYcolors.
#' @param percent.change # between 0 and 1. the percentage to darken by. Defaults to 0.25 if not given.
#' @param relative TRUE/FALSE. Whether the percentage should be a relative change versus an absolute one. Default = TRUE.
#' @return Return a darkened version of the color in hexadecimal color form (="#RRGGBB" in base 16)
#' @examples
#' Darken("blue") #"blue" = "#0000FF"
#' #Output: "#0000BF"
#' Darken(MYcolors[1:8]) #Works for multiple color inputs as well.

Darken <- function(colors, percent.change = 0.25, relative = TRUE){

  colorspace::darken(colors, amount = percent.change, space = "HLS", fixup = TRUE, method = ifelse(relative,"relative","absolute"))
}

#### Lighten: For lightening colors ####
#' Lightens input colors by a set amount
#'
#' @description A wrapper for the lighten function of the colorspace package.
#' @param colors the color(s) input. Can be a list of colors, for example, MYcolors.
#' @param percent.change # between 0 and 1. the percentage to darken by. Defaults to 0.25 if not given.
#' @param relative TRUE/FALSE. Whether the percentage should be a relative change versus an absolute one. Default = TRUE.
#' @return Return a lighter version of the color in hexadecimal color form (="#RRGGBB" in base 16)
#' @examples
#' Lighten("blue") #"blue" = "#0000FF"
#' #Output: "#4040FF"
#' Lighten(MYcolors[1:8]) #Works for multiple color inputs as well.

Lighten <- function(colors, percent.change = 0.25, relative = TRUE){

  colorspace::lighten(colors, amount = percent.change, space = "HLS", fixup = TRUE, method = ifelse(relative,"relative","absolute"))
}

#### Simulate: For simulating what a plot would look like as seen by a colorblind person ####
#' Simulates what a colorblind person would see for any DittoSeq plot!
#'
#' @description Essentially a wrapper function for colorspace's deutan(), protan(), and tritan() functions. This function will output any DittoSeq plot as it might look to an individual with one of the common forms of colorblindness: deutanopia/deutanomaly, the most common, is when the cones mainly responsible for red vision are defective. Protanopia/protanomaly is when the cones mainly responsible for green vision are defective. In tritanopia/tritanomaly, the defective cones are responsible for blue vision. Note: there are more severe color deficiencies that are even more rare. Unfortunately, for these types of color vision deficiency, only non-color methods, like lettering or shapes, will do much to help.
#' @param type The type of colorblindness that you want to simulate for. Options: "deutan", "protan", "tritan". Anything else, and you will get an error.
#' @param plot.function The plotting function that you want to use/simulate. not quoted. and make sure to remove the () that R will try to add.
#' @param ... other paramters that can be given to DittoSeq plotting functions, including color.panel, used in exactly the same way they are used for those functions. (contrary to the look of this documentation, color.panel will still default to MYcolors when not provided.)
#' @param color.panel The set of colors to be used.  Not required to be given, as contrary to the look of this documentation, it will still default to MYcolors when not provided.
#' @return Outputs a DittoSeq plot with the color.panel updated as it might look to a colorblind individual. Note: Does not currently work for DBHeatmap or for continuous variable plotting in DBDimPlot.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' Simulate("deutan", DBDimPlot, var = "RNA_snn_res.1", object = "pbmc", size = 2)
#' Simulate("protan", DBDimPlot, "RNA_snn_res.1", "pbmc", size = 2)
#' Simulate("tritan", DBDimPlot, "RNA_snn_res.1", "pbmc", size = 2)

Simulate <- function(type = "deutan", plot.function, ..., color.panel = MYcolors){

  #Check that type was given properly
  if(!(type=="deutan"|type=="protan"|type=="tritan")){
    return("Error: type must be 'deutan', 'protan', or 'tritan'")
  }

  #Simulate the color panel for the given color blindness type.
  color.p <- eval(expr = parse(text = paste0("colorspace::",type,"(color.panel)")))

  #Make the plot!
  plot.function(color.panel = color.p, ... )
}
