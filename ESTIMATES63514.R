library(utils)
rforge <- "http://r-forge.r-project.org"
library(estimate)
library(tidyverse)
help(package="estimate")

#This code were used#
CervicalCancerExpr <- system.file("extdata", "63514_raw.txt", package="estimate")
filterCommonGenes(input.f= CervicalCancerExpr, output.f="gene_expression63514.gct", id ="GeneSymbol")
gene_expression63514 <- system.file("extdata", "gene_expression63514.gct", package="estimate")
estimateScore(gene_expression63514,"EstimateScore63514.gct", platform = c("affymetrix", "agilent", "illumina"))
Estimatescore63514 <- system.file("extdata", "EstimateScore63514.gct", package="estimate")
plotPurity(Estimatescore63514, samples="all_samples", platform = c("affymetrix", "agilent", "illumina"),
output.dir="estimated_purity_plots63514")

#######################################################################################################
Full code version
#Part1
# Filter Common genes
CervicalCancerExpr <- system.file("extdata", "63514_raw.txt", package="estimate")
filterCommonGenes(input.f= CervicalCancerExpr, output.f="gene_expression63514.gct", id ="GeneSymbol")


## Check arguments
stopifnot((is.character(CervicalCancerExpr) && length(CervicalCancerExpr) == 1 && nzchar(CervicalCancerExpr)) ||
              (inherits(CervicalCancerExpr, "connection") && isOpen(CervicalCancerExpr, "r")))
stopifnot((is.character("gene_expression63514.gct") && length("gene_expression63514.gct") == 1 && nzchar("gene_expression63514.gct")))
possible_values <- c("GeneSymbol", "EntrezID")
id <- "GeneSymbol"
id <- match.arg(id, choices = possible_values)

## Load estimate.RData
data(common_genes)

## Read input data
input.data <- read.table(CervicalCancerExpr, header=TRUE, row.names=1, sep="\t",
                           quote="", stringsAsFactors=FALSE)

if(id == "GeneSymbol")
  input.f1 <- merge (common_genes,input.data, by.x="GeneSymbol", by.y="row.names")
  rownames(input.f1) <- input.f1$GeneSymbol
  input.f1 <- input.f1[,-1:-ncol(common_genes)]
  x.mismatched.genes <- nrow(common_genes) - nrow(input.f1)
  print(paste("This dataset includes ", nrow(input.f1), "genes.", sep=""))
  print(paste(x.mismatched.genes, "genes were mismatched.", sep=""))
  outputGCT(input.f1,output.f  = "gene_expression63514.gct")

#Part2 ให้ผลเหมือนกับการทำfunction outputGCT ในPart1
    if(is.data.frame(input.f1)==TRUE)
      exp.data <- input.f1 else
      exp.data <- read.table(input.f1, header=TRUE, row.names=1, sep="\t", quote="" )


    exp.data1 <- data.frame(NAME=rownames(exp.data),Description=rownames(exp.data),exp.data)
    column1 <- colnames(exp.data1)
    column1[1] <- "NAME"
    column1[2] <- "Description"
    exp.data1$NAME <- factor(exp.data1$NAME)
    exp.data1$Description <- factor(exp.data1$Description)
    levels(exp.data1[,1]) <- c(levels(exp.data1[,1]),"NAME")
    levels(exp.data1[,2]) <- c(levels(exp.data1[,2]),"Description")
    exp.data2 <- rbind(column1,exp.data1)

    row1 <- rep("",length(1:ncol(exp.data)))
    row1_2 <- data.frame(row1,row1)
    row1_2 <- t(row1_2)
    No_gene <- nrow(exp.data1)
    No_sample <- (ncol(exp.data1)-2)
    GCT <- matrix(c("#1.2",No_gene,"",No_sample),nrow=2,ncol=2)
    gct <- cbind(GCT,row1_2)
    colnames(gct) <- colnames(exp.data2)
    tmp <- rbind(gct,exp.data2)
    write.table(tmp,"Data.gene.gct", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    invisible(NULL)


#Part3
    gene_expression63514 <- system.file("extdata", "gene_expression63514.gct", package="estimate")
    estimateScore(gene_expression63514,"EstimateScore63514.gct", platform = c("affymetrix", "agilent", "illumina"))
    ?estimateScore
      ## Check arguments
      stopifnot(is.character("gene_expression63514.gct") && length("gene_expression63514.gct") == 1 && nzchar("gene_expression63514.gct"))
      stopifnot(is.character("EstimateScore63514") && length("EstimateScore63514") == 1 && nzchar("EstimateScore63514"))
      possible_values1 <- c("affymetrix", "agilent", "illumina")
      platform <- "affymetrix"
      platform <- match.arg(platform, choices = possible_values1)


      ## Load estimate.RData
      data(SI_geneset)

      ## Read input dataset(GCT format)
      ds <- read.delim(gene_expression63514, header=TRUE, sep="\t", skip=2, row.names=1, blank.lines.skip=TRUE, as.is=TRUE, na.strings = "")
      descs <- ds[,1]
      ds <- ds[-1]
      row.names <- row.names(ds)
      names <- names(ds)
      dataset <- list(ds = ds, row.names = row.names, descs = descs, names = names)

      m <- data.matrix(dataset$ds)
      gene.names <- dataset$row.names
      sample.names <- dataset$names
      Ns <- length(m[1, ]) # Number of genes
      Ng <- length(m[, 1]) # Number of samples
      temp <- strsplit(gene_expression63514, split="/")
      s <- length(temp[[1]])
      input.file.name <- temp[[1]][s]
      temp <- strsplit(input.file.name, split=".gct")
      input.file.prefix <-  temp[[1]][1]

      ## Sample rank normalization
      for (j in 1:Ns) {
        m[, j] <- rank(m[, j], ties.method="average")
      }
      m <- 10000*m/Ng

      ## SI_geneset
      gs <- as.matrix(SI_geneset[,-1],dimnames=NULL)
      N.gs <- 2
      gs.names <- row.names(SI_geneset)

      ## Loop over gene sets
      score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
      for (gs.i in 1:N.gs) {
        gene.set <- gs[gs.i,]
        gene.overlap <- intersect(gene.set, gene.names)
        print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
        if (length(gene.overlap) == 0) {
          score.matrix[gs.i, ] <- rep(NA, Ns)
          next
        } else {
          ES.vector <- vector(length=Ns)

          ## Enrichment score
          for (S.index in 1:Ns) {
            gene.list <- order(m[, S.index], decreasing=TRUE)
            gene.set2 <- match(gene.overlap, gene.names)
            correl.vector <- m[gene.list, S.index]

            TAG <- sign(match(gene.list, gene.set2, nomatch=0))    # 1 (TAG) & 0 (no.TAG)
            no.TAG <- 1 - TAG
            N <- length(gene.list)
            Nh <- length(gene.set2)
            Nm <-  N - Nh
            correl.vector <- abs(correl.vector)^0.25
            sum.correl  <- sum(correl.vector[TAG == 1])
            P0 <- no.TAG / Nm
            F0 <- cumsum(P0)
            Pn <- TAG * correl.vector / sum.correl
            Fn <- cumsum(Pn)
            RES <- Fn - F0
            max.ES <- max(RES)
            min.ES <- min(RES)
            if (max.ES > - min.ES) {
              arg.ES <- which.max(RES)
            } else {
              arg.ES <- which.min(RES)
            }
            ES <- sum(RES)
            EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = TAG)
            ES.vector[S.index] <- EnrichmentScore$ES
          }
          score.matrix[gs.i, ] <- ES.vector
        }
      }

      score.data <- data.frame(score.matrix)
      names(score.data) <- sample.names
      row.names(score.data) <- gs.names
      estimate.score <- apply(score.data, 2, sum)

      if (platform != "affymetrix"){
        score.data <- rbind(score.data, estimate.score)
        rownames(score.data) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
      } else {
        ## Calculate ESTIMATE-based tumor purity (This function is limited in the Affymetrix platform)
        convert_row_estimate_score_to_tumor_purity <- function(x) {
          stopifnot(is.numeric(x))
          cos(0.6049872018 + 0.0001467884*x)
        }

        est.new <- NULL
        for (i in 1:length(estimate.score)) {
          est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
          est.new <- rbind(est.new, est_i)
          if(est_i >= 0){
            next
          } else {
            message(paste(sample.names[i],": out of bounds", sep=""))
          }
        }
        colnames(est.new) <- c("TumorPurity")
        estimate.t1 <- cbind(estimate.score, est.new)
        x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 0
        estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
        score.data <- rbind(score.data, t(estimate.t1))
        rownames(score.data) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
      }
      outputGCT(score.data, "EstimateScore63514(2).gct")

      Score63514 <- system.file("extdata", "EstimateScore63514(1).gct", package="estimate")
      scores <- read.table(Score63514 , skip = 2, header = TRUE, sep = "\t", row.names = 1)
      head(scores)

#Part4
      plotPurity <- function(scores,
                             samples="all_samples",
                             platform = c("affymetrix", "agilent", "illumina"),
                             output.dir="estimated_purity_plots") {

        ## Check arguments
        stopifnot((is.character(scores) && length(scores) == 1 && nzchar(scores)) ||
                    (inherits(scores, "connection") && isOpen(scores, "r")))
        stopifnot(is.character(output.dir) && length(output.dir) == 1 && nzchar(output.dir))
        platform <- match.arg(platform)

        if (platform != "affymetrix"){
          stop("not implemented")
        }
        ## Begin processing

        ##---------------------------------------------------------------------------
        get_estimates_df <- function(scores) {
          #estimate <- read.table(scores, skip=2, header=TRUE, row.names=1, sep="\t")
          estimate <- read.delim(scores, skip=2, row.names=1)
          as.data.frame(t(estimate[, -1]))
        }

        ##---------------------------------------------------------------------------
        convert_row_estimate_score_to_tumor_purity <- function(x) {
          stopifnot(is.numeric(x))
          cos(0.6049872018 + 0.0001467884*x)
        }

        ## Load estimate.RData
        data(PurityDataAffy)

        ## Read ESTIMATE data file
        estimate.df <- get_estimates_df(scores)
        samplenames <- rownames(estimate.df)
        Affy.model <- PurityDataAffy
        pred.p <- Affy.model[, 5:7]
        est <- estimate.df[, 3]
        est.new <- estimate.df[, 4]

        ## Create output directory
        dir.create(output.dir)

        ## ESTIMATE based tumor purity in scatterplot with prediction interval
        message("Plotting tumor purity based on ESTIMATE score")

        max.af <- max(Affy.model$ESTIMATEScore)
        min.af <- min(Affy.model$ESTIMATEScore)

        if (samples[1] == "all_samples"){
          Num.S <- nrow(estimate.df)
        } else {
          Num.S <- as.numeric(length(samples))
        }

        for (i in 1:Num.S) {
          if(samples[1] =="all_samples"){
            samplename <- samplenames[i]
          } else {
            samplename <- samples[i]
          }

          png.filename <- file.path(output.dir, sprintf("%s.png", samplename))
          png(filename=png.filename, width=480, height=480)

          geMin <- est[i] >= min.af
          leMax <- est[i] <= max.af
          withinMinMax <- geMin && leMax

          xlim <- if (!withinMinMax) {
            ## Expands plot boundary
            adjustment <- 500    # Arbitrary
            if (geMin) {
              from <- min.af
              to   <- est[i] + adjustment
            } else {
              from <- est[i] - adjustment
              to   <- max.af
            }
            c(from, to)
          } else {
            NULL
          }

          plot(Affy.model$tumor.purity~Affy.model$ESTIMATEScore, Affy.model,
               main=samplename,
               type="n",
               xlab="ESTIMATE score",
               xlim=xlim,
               ylab="Tumor purity",
               ylim=c(0, 1))
          par(new=TRUE)
          points(Affy.model$ESTIMATEScore, Affy.model$tumor.purity, cex=0.75, col="lightgrey")
          if (withinMinMax) {
            ## Prediction interval
            matlines(Affy.model$ESTIMATEScore, pred.p, lty=c(1, 2, 2), col="darkgrey")
          } else {
            matlines(Affy.model$ESTIMATEScore, pred.p, lty=c(1, 2, 2), col="darkgrey")
            par(new=TRUE)
            curve(convert_row_estimate_score_to_tumor_purity,
                  from, to, n=10000, col="grey", ylim=c(0, 1), xlab="", ylab="")
          }
          points(est[i], est.new[i], pch=19, cex=1.25)
          abline(h=est.new[i], col="black", lty=2)
          abline(v=est[i], col="black", lty=2)

          dev.off()
        }
      }
