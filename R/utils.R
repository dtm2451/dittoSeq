.remove_legend <- function(ggplot) {
    ggplot + theme(legend.position = "none")
}

.grab_legend <- function(ggplot) {
    cowplot::ggdraw(cowplot::get_legend(ggplot))
}

.all_cells <- function(object = DEFAULT) {
    target <- data.frame(
        use = c("@samples", "@cell.names"),
        row.names = c("RNAseq", "Seurat.v2"))
    if (.class_of(object) %in% c("Seurat.v3", "SingleCellExperiment")) {
        return(colnames(object))
    } else {
        eval(expr = parse(text = paste0("object", target[.class_of(object),])))
    }
}

.which_cells <- function(cells.use, object = DEFAULT) {
    all.cells <- .all_cells(object)
    if (is.null(cells.use)) {
        return(all.cells)
    }
    if (is.logical(cells.use)) {
        return(all.cells[cells.use])
    }
    cells.use
}

.class_of <- function (object = DEFAULT) {
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    class <- class(object)

    #If a Seurat, need to add what version
    if (grepl("Seurat|seurat",class)) {
        if(object@version >= '3.0.0') {
            #Then needs to be
            class <- "Seurat.v3"
        } else {
            class <- "Seurat.v2"
        }
    }

    class
}

.var_OR_get_meta_or_gene <- function(var, object = DEFAULT, data.type) {
    OUT <- var
    if (length(var)==1 && typeof(var)=="character") {
        if (is.meta(var, object)) {
            OUT <- meta(var, object)
        }
        if (is.gene(var, object)) {
            OUT <- gene(var, object, data.type)
        }
    }
    names(OUT) <- .all_cells(object)
    OUT
}

.which_data <- function(data.type, object=DEFAULT) {
    #Set up data frame for establishing how to deal with different input object types
    target <- data.frame(
        RNAseq = c(
            "@data",
            "@counts",
            "error_Do_not_use_scaled_for_RNAseq_objects"),
        Seurat.v2 = c("@data", "@raw.data", "@scale.data"),
        Seurat.v3 = c("nope", "counts", "scale.data"),
        SingleCellExperiment = c(
            "SingleCellExperiment::logcounts(",
            "SingleCellExperiment::counts(",
            "error_do_not_use_scaled_for_SCE_objects("),
        stringsAsFactors = FALSE,
        row.names = c("normalized","raw","scaled"))
    if (.class_of(object) == "SingleCellExperiment") {
        OUT <- as.matrix(
            eval(expr = parse(text = paste0(
                target[data.type,.class_of(object)],
                "object)"))))
    } else {
        if (.class_of(object)=="Seurat.v3") {
            if(data.type == "normalized"){
                OUT <- Seurat::GetAssayData(object)
            } else {
                OUT <- Seurat::GetAssayData(
                    object,
                    slot = target[data.type,.class_of(object)])
            }
        } else {
            OUT <- eval(expr = parse(text = paste0(
                "object",
                target[data.type,.class_of(object)]
            )))
        }
    }
    OUT
}

.extract_Reduced_Dim <- function(reduction.use, dim=1, object=DEFAULT) {
    # If object is a Seurat object
    if (.class_of(object)=="Seurat.v2") {
        OUT <- list(eval(expr = parse(text = paste0(
            "object@dr$",reduction.use,"@cell.embeddings[,",dim,"]"))))
        OUT[2] <- paste0(eval(expr = parse(text = paste0(
            "object@dr$",reduction.use,"@key"))),dim)
    }
    if (.class_of(object)=="Seurat.v3") {
        OUT <- list(eval(expr = parse(text = paste0(
            "object@reductions$",reduction.use,"@cell.embeddings[,",dim,"]"))))
        OUT[2] <- paste0(eval(expr = parse(text = paste0(
            "object@reductions$",reduction.use,"@key"))),dim)
    }
    if (.class_of(object)=="RNAseq"){
        OUT <- list(eval(expr = parse(text = paste0(
            "object@reductions$",reduction.use,"$embeddings[,",dim,"]"))))
        OUT[2] <- paste0(eval(expr = parse(text = paste0(
            "object@reductions$",reduction.use,"$key"))),dim)
    }
    if (.class_of(object)=="SingleCellExperiment"){
        OUT <- list(eval(expr = parse(text = paste0(
            "SingleCellExperiment::reducedDim(",
            "object, type = '",reduction.use,"')[,",dim,"]"))))
        OUT[2] <- paste0(.gen_key(reduction.use),dim)
    }

    names(OUT) <- c("embeddings","name")
    OUT
}

.gen_key <- function (reduction.use){
    key <- reduction.use
    if (grepl("pca|PCA", reduction.use)) {
        key <- "PC"
    }
    if (grepl("cca|CCA", reduction.use)) {
        key <- "CC"
    }
    if (grepl("cca.aligned", reduction.use)) {
        key <- "aligned.CC"
    }
    if (grepl("ica|ICA", reduction.use)) {
        key <- "IC"
    }
    if (grepl("tsne|tSNE|TSNE", reduction.use)) {
        key <- "tSNE_"
    }
    key
}

.make_hover_strings_from_vars <- function(data.hover, object, data.type) {
    # Overall: if do.hover=TRUE and data.hover has a list of genes / metas called
    #   c(var1, var2, var3, ...), then for all cells, make a string:
    #   "var1: var1-value\nvar2: var2-value\nvar3: var3-value\n..."
    #   vars that are not genes of metadata are ignored.
    fillable <- vapply(
        seq_along(data.hover),
        function(i)
            (is.meta(data.hover[i],object) |
                is.gene(data.hover[i],object) |
                (data.hover[i]=="ident")),
        logical(1))
    data.hover <- data.hover[fillable]
    if (is.null(data.hover)) {
        stop("No genes or metadata names added to `hover.data`")
    }

    # Create dataframe to contain the hover.info
    features.info <- data.frame(row.names = .all_cells(object))
    features.info <- vapply(
        data.hover,
        function(this.data)
            as.character(.var_OR_get_meta_or_gene(this.data,object, data.type)),
        character(nrow(features.info)))
    names(features.info) <- data.hover[fillable]

    # Convert each row of dataframe to 'colname1: data1\ncolname2: data2\n...'
    hover.strings <- .make_hover_strings_from_df(features.info)
}

.make_hover_strings_from_df <- function(df){
    vapply(
        seq_len(nrow(df)),
        function(row){
            paste(as.character(vapply(
                seq_len(ncol(df)),
                function(col){
                    paste0(names(df)[col],": ",df[row,col])
                }, FUN.VALUE = character(1))
                ),collapse = "\n")
        }, FUN.VALUE = character(1))
}

.rename_and_or_reorder <- function(orig.data, reorder = NULL, relabels = NULL) {
    if (is.numeric(orig.data)) {
        return(orig.data)
    }
    rename.args <- list(x = orig.data)
    if (!(is.null(reorder))) {
        rename.args$levels <- levels(factor(rename.args$x))[reorder]
    }
    if (!(is.null(relabels))) {
        rename.args$labels <- relabels
    }
    do.call(factor, args = rename.args)
}

.leave_default_or_null <- function(
    target, default, null.if = FALSE, default.when = "make") {
    if (!is.null(target)) {
        if (target==default.when) {
            if (null.if) {
                NULL
            } else {
                target <- default
            }
        }
    }
    target
}
