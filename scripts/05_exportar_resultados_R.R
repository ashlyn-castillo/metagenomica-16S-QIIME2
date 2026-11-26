# Exportar resultados existentes SIN repetir DADA2 ni la optimización del árbol.
# Ejecutar en la sesión de RStudio que contiene los resultados del taller.
# Si se reabrió RStudio, cargar antes Analisis_16S.RData y las bibliotecas.
library(phyloseq)
library(ggplot2)
destino <- choose.dir(caption = "Seleccione la carpeta portafolio_16S extraída")
if (is.na(destino)) stop("Exportación cancelada")
salida <- file.path(destino, "resultados", "R")
dir.create(salida, recursive = TRUE, showWarnings = FALSE)
requeridos <- c("track", "alfa", "seqtab.nochim", "taxa")
faltantes <- requeridos[!vapply(requeridos, exists, logical(1), inherits = TRUE)]
if (length(faltantes)) stop("Faltan objetos: ", paste(faltantes, collapse = ", "))
write.csv(track, file.path(salida, "seguimiento_lecturas.csv"))
write.csv(alfa, file.path(salida, "diversidad_alfa.csv"))
write.csv(taxa, file.path(salida, "taxonomia_SILVA138_2.csv"))
write.csv(data.frame(Observed = rowSums(seqtab.nochim > 0)),
          file.path(salida, "riqueza_observada.csv"))
for (nombre in c("jaccard", "bray", "d_unifrac", "d_wunifrac")) {
 if (exists(nombre, inherits = TRUE)) {
  write.csv(as.matrix(get(nombre)), file.path(salida, paste0(nombre, ".csv")))
 }
}
for (nombre in c("g_filos", "g_generos", "g_top10", "g_heatmap")) {
 if (exists(nombre, inherits = TRUE)) {
  ggsave(file.path(salida, paste0(nombre, ".png")), plot = get(nombre),
         width = 11, height = 7, dpi = 300)
 }
}
if (exists("coord") && all(c("Eje1", "Eje2", "Muestra") %in% names(coord))) {
 g <- ggplot(coord, aes(Eje1, Eje2)) + geom_point(size = 4) +
  geom_text(aes(label = Muestra), vjust = -0.8) + theme_minimal() +
  labs(title = "PCoA de Bray-Curtis", x = "Eje 1", y = "Eje 2")
 ggsave(file.path(salida, "PCoA_Bray_Curtis.png"), plot = g,
        width = 9, height = 7, dpi = 300)
}
if (exists("arbol_vista")) {
 pdf(file.path(salida, "arbol_20_ASVs.pdf"), width = 12, height = 8)
 tryCatch({plot(arbol_vista, show.node.label = FALSE, cex = 0.8,
                main = "Árbol de las 20 ASVs con más lecturas")
           ape::add.scale.bar()}, finally = dev.off())
}
writeLines(capture.output(sessionInfo()), file.path(salida, "versiones_R.txt"))
cat("Exportación terminada en:", salida, "\n")
