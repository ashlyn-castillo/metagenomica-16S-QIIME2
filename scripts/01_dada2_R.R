# Adaptación de la práctica Flyteek de Alex Andrés Jarrín Jurado.
# Versión reorganizada para el portafolio; no ejecutada de extremo a extremo.
# Declaracion de la libreria
library(dada2)
packageVersion("dada2")
# Ruta definida
path <- Sys.getenv("RAW_DIR")
if (!dir.exists(path)) stop("Configure RAW_DIR con la carpeta FASTQ original")
list.files(path)
# Procesamiento de nombres
fnFs <- sort(list.files(path,pattern="_R1_001\\.fastq$", full.names = TRUE))
fnRs <- sort(list.files(path,pattern="_R2_001\\.fastq$", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs),"_"), `[`,1)
# Selecciona solo FASTQ sin comprimir para evitar duplicados .fastq.gz
stopifnot(length(fnFs) == 5L, length(fnRs) == 5L)
stopifnot(identical(sub("_R1_001.fastq$", "", basename(fnFs)),
                    sub("_R2_001.fastq$", "", basename(fnRs))))
# Graficos de calidad, inspeccion
plotQualityProfile(fnFs[1:1])
plotQualityProfile(fnRs[1:1])
# Corte y filtrado
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)


# MODELO Y TAZAS DE ERROR

errF <- learnErrors(filtFs, multithread=  TRUE) 
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# INFERENCIA DE LA MUESTRA - APLICACION DEL  ALGORITMO DADA2

dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)
dadaFs[[1]]

# FUSIONAR LAS LECTURAS EMPAREJADAS

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) # Mantener los valores predeterminados de mergePairs
head(mergers[[1]])

# CONSTRUIR TABLAS ASVs

seqtab <- makeSequenceTable(mergers)
dim(seqtab)



table(nchar(getSequences(seqtab)))

# ELIMINAR QUIMNERAS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# SEGUIMIENTO DE LAS LECTURAS
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("entrada", "filtrados", "SinRuidoF", "SinRuidoR", "Fusionaddas", "NoQuimericas")
rownames(track) <- sample.names
head(track)

# ASIGNAR TAXONOMIA

taxa <- assignTaxonomy(
  seqtab.nochim,
  Sys.getenv("SILVA_TRAIN"),
  multithread = TRUE
)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

taxa <- addSpecies(
  taxa,
  Sys.getenv("SILVA_SPECIES")
)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

table(!is.na(taxa[, "Family"]))
table(!is.na(taxa[, "Species"]))
sort(table(taxa[, "Genus"]), decreasing = TRUE)[1:20]
