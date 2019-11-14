# original library
library(rusda)

# test parameters
x = "Protomyces inouyei"
database = "FH"
spec_type <- "fungus"

# dependencies
library(foreach)
library(XML)
library(httr)
library(RCurl)
library(testthat)

# other tweaked function
source("./code/getHF_function.R")

# tweaked function
associations2 <- function (x, database = c("FH", "SP", "both"), spec_type = c("plant", 
                                                                             "fungus"), clean = TRUE, syn_include = TRUE, process = TRUE) 
{
  if (!url.exists("r-project.org") == TRUE) 
    stop("Not connected to the internet. Please create a stable connection and try again.")
  if (!is.character(getURL("https://nt.ars-grin.gov/fungaldatabases/index.cfm"))) 
    stop(" Database is not available : http://nt.ars-grin.gov/fungaldatabases/index.cfm")
  expect_match(spec_type, ("fungus|plant"))
  expect_match(database, ("FH|SP|both"))
  if (length(grep("_", x)) > 0) {
    x <- gsub("_", " ", x)
  }
  words <- vapply(strsplit(x, "\\W+"), length, integer(1))
  if (any(words == 1) & any(words == 2)) {
    stop(paste(" check if you specified ONLY genus names or ONLY species names \n", 
               "AFAICS you provided:  \n", sum(words == 1), "  genus name(s)  \n", 
               sum(words == 2), "  species name(s) ", sep = ""))
  }
  if (all(words == 1)) {
    x <- lapply(x, ncbiSpecies, clean = TRUE, sub = FALSE)
    x <- unlist(x)
  }
  if (length(grep("\\sx\\s", x)) > 0) 
    stop(" no hybrids allowed as input ")
  tax <- strsplit(x, " ")
  if (process == TRUE) {
    message("... retrieving data ... for:")
  }
  p <- foreach(i = seq_along(tax)) %do% getHF(tax[[i]], process, 
                                              spec_type = spec_type)
  taxa <- lapply(tax, function(x) {
    paste(as.character(x[1]), as.character(x[2]))
  })
  co <- lapply(p, getCOND)
  if (process == TRUE) {
    message("... extracting Synonyms ...")
  }
  syns <- lapply(p, getSYNS, process = process, taxa = taxa)
  names(syns) <- taxa
  if (process == TRUE & database == "FH" | database == "both") {
    message("... extracting Fungus-Hosts DB ...")
  }
  i <- NULL
  hosts_hf <- foreach(i = seq_along(taxa)) %do% {
    if (length(co[[i]]$hfu) == 0 | length(co[[i]]$hf.st) == 
        0) {
      hf <- "nodata"
    }
    if (length(co[[i]]$hf.st) > 0) {
      hf.c <- grep("The Literature database has", p[[i]])
      ifelse(length(hf.c) > 0, hf.sp <- hf.c, hf.sp <- (grep("No records were found in the Literature database", 
                                                             p[[i]])))
      if (length(hf.sp) == 0) {
        hf.sp <- (grep(paste("There are no records for ", 
                             taxa[[i]], " in the Literature database", sep = ""), 
                       p[[i]]))
      }
      p[[i]][(co[[i]]$hf.st + 1):(hf.sp - 1)]
    }
  }
  names(hosts_hf) <- unlist(taxa)
  if (process == TRUE & database == "SP" | database == "both") {
    message("... extracting Specimens DB ...")
  }
  i <- NULL
  hosts_sp <- foreach(i = seq_along(taxa)) %do% {
    if (length(co[[i]]$sp) == 0 | length(co[[i]]$spe.st) == 
        0) {
      specim <- "nodata"
    }
    if (length(co[[i]]$spe.st) > 0) {
      spe.sp <- grep("Systematic Mycology and Microbiology Laboratory[.]", 
                     p[[i]])
      spe.sp <- spe.sp[length(spe.sp)]
      specim <- p[[i]][(co[[i]]$spe.st + 1):(spe.sp - 1)]
    }
  }
  names(hosts_sp) <- unlist(taxa)
  if (syn_include == FALSE) {
    if (process == TRUE) {
      message("... excluding synonyms ...")
    }
    no_syns <- function(x) {
      st <- foreach(i = seq_along(taxa)) %do% grep(taxa[[i]], 
                                                   x[[i]])
      sp <- foreach(i = seq_along(taxa)) %do% {
        sy <- paste(syns[[i]][!syns[[i]] == taxa[[i]]], 
                    collapse = "|")
        grep(sy, x[[i]], value = FALSE)
      }
      spp <- list()
      for (i in seq_along(taxa)) {
        if (is.integer(sp[[i]]) && length(sp[[i]]) == 
            0L) {
          spp[[i]] <- integer(0)
        }
        if (length(st[[i]]) > 0 & length(sp[[i]]) > 0) 
          spp[[i]] <- sp[[i]][(sp[[i]] > st[[i]][1]) & 
                                sp[[i]] == ((min(st[[i]][1] - sp[[i]]) * 
                                               -1) + st[[i]][1])]
      }
      res <- list()
      for (i in seq_along(taxa)) {
        if (length(st[[i]]) == 0 & is.integer(spp[[i]]) && 
            length(spp[[i]]) == 0L) {
          res[[i]] <- x
        }
        if (length(st[[i]]) > 0 & length(spp[[i]]) == 
            0L) 
          res[[i]] <- x[[i]][st[[i]][1]:length(x[[i]])]
        if (length(st[[i]]) > 0 & length(spp[[i]]) > 
            0) 
          if (length(st[[i]]) > 0) 
            res[[i]] <- x[[i]][st[[i]][1]:(spp[[i]] - 
                                             1)]
          else res[[i]] <- x[[i]][st[[i]]:(spp[[i]] - 
                                             1)]
      }
      return(res)
    }
    res <- lapply(list(hosts_hf, hosts_sp), no_syns)
    hosts_hf <- res[[1]]
    hosts_sp <- res[[2]]
  }
  if (database == "FH") {
    res <- hosts_hf
  }
  if (database == "SP") {
    res <- hosts_sp
  }
  if (database == "both") {
    res <- foreach(i = seq_along(hosts_hf)) %do% c(hosts_hf[[i]], 
                                                   hosts_sp[[i]])
    names(res) <- names(hosts_hf)
    res <- lapply(res, function(x) {
      if (length(grep("nodata", x)) == 2) {
        x <- "nodata"
      }
      if (!length(grep("nodata", x)) == 2) {
        x
      }
    })
  }
  if (clean == TRUE) {
    if (process == TRUE) {
      message("... cleaning step ...")
    }
    res <- lapply(res, clean_step, taxa = taxa, syns = syns, 
                  spec_type = spec_type, synonyms_incl = TRUE)
  }
  res <- lapply(res, unique)
  names(res) <- taxa
  return(list(synonyms = syns, associations = res))
}
