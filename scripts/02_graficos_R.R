# Adaptación de la práctica Flyteek de Alex Andrés Jarrín Jurado.
# Versión reorganizada para el portafolio; no ejecutada de extremo a extremo.
# Ejecutar después de 01_dada2_R.R en la misma sesión.

#  PREPARACION Y OBJETO PHYLOSEQ

library(phyloseq)                  # Integra abundancia, taxonomia y metadatos
library(ggplot2)                   # Crea graficos

# Los objetos seqtab.nochim, taxa y track vienen del pipeline DADA2.

# Una fila de metadata por muestra. 
samdf <- data.frame(
  SampleID = rownames(seqtab.nochim),
  row.names = rownames(seqtab.nochim)
)

# Comprobar que todas las ASVs tengan una fila en la tabla taxonomica
stopifnot(all(colnames(seqtab.nochim) %in% rownames(taxa)))

# Ordenar la taxonomia para hacerla coincidir con la tabla de ASVs
taxa2 <- taxa[colnames(seqtab.nochim), , drop = FALSE]

# Crear objeto unificado: muestras x ASVs, taxonomia y metadata
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa2),
  sample_data(samdf)
)

# Verificacion del objeto
ps
ntaxa(ps)                         # Cantidad total de ASVs
nsamples(ps)                      # Cantidad de muestras
sample_sums(ps)                   # Lecturas finales de cada muestra

# Convertir la tabla a formato largo para usar ggplot2
df <- psmelt(ps)

# Marcar la falta de clasificacion sin borrar las ASVs
df$Phylum <- as.character(df$Phylum)
df$Genus <- as.character(df$Genus)
df$Phylum[is.na(df$Phylum) | df$Phylum == ""] <- "Sin asignar"
df$Genus[is.na(df$Genus) | df$Genus == ""] <- "Sin asignar"


# RETENCION DE LECTURAS
# Lecturas finales / lecturas iniciales, para cada muestra
retencion <- 100 * track[, "NoQuimericas"] / track[, "entrada"]
round(retencion, 2)

barplot(
  retencion,
  col = "#20A39E",
  ylim = c(0, 100),
  main = "Retención de lecturas por muestra",
  ylab = "Retención (%)",
  xlab = "Muestras"
)



# COMPOSICION RELATIVA A NIVEL FILO
# Agrupar (sumar) todas las lecturas de las ASVs de cada filo
filos <- aggregate(Abundance ~ Sample + Phylum, data = df, sum)

# En cada muestra, calcular el porcentaje que representa cada filo
filos$Porcentaje <- 100 * filos$Abundance /
  ave(filos$Abundance, filos$Sample, FUN = sum)

g_filos <- ggplot(
  filos,
  aes(x = Sample, y = Porcentaje, fill = Phylum)
) +
  geom_col() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Composición bacteriana por filo",
    x = "Muestras", y = "Abundancia relativa (%)", fill = "Filo"
  )

g_filos



# COMPOSICION RELATIVA POR GENERO
# Ranking global: solo los generos con una asignacion disponible
ranking <- aggregate(
  Abundance ~ Genus,
  data = df[df$Genus != "Sin asignar", ],
  sum
)
ranking <- ranking[order(ranking$Abundance, decreasing = TRUE), ]
top10 <- head(ranking$Genus, 10)

# Mantener los 10 generos principales; los demas quedan como "Otros"
# Las ASVs sin genero identificado permanecen como "Sin asignar"
df$GrupoGenero <- ifelse(
  df$Genus == "Sin asignar", "Sin asignar",
  ifelse(df$Genus %in% top10, df$Genus, "Otros")
)

# Sumar lecturas por muestra y categoria
generos <- aggregate(
  Abundance ~ Sample + GrupoGenero, data = df, sum
)

# Convertir a porcentaje dentro de cada muestra
generos$Porcentaje <- 100 * generos$Abundance /
  ave(generos$Abundance, generos$Sample, FUN = sum)

g_generos <- ggplot(
  generos,
  aes(x = Sample, y = Porcentaje, fill = GrupoGenero)
) +
  geom_col() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Top 10 géneros bacterianos por muestra",
    x = "Muestras", y = "Abundancia relativa (%)", fill = "Género"
  )

g_generos


# RANKING GLOBAL DE GENEROS
# Mostrar el ranking exacto (suma de lecturas de las cinco muestras)
head(ranking, 10)

g_top10 <- ggplot(
  head(ranking, 10),
  aes(x = reorder(Genus, Abundance), y = Abundance)
) +
  geom_col(fill = "#168AAD") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Los 10 géneros más abundantes",
    x = "Género bacteriano", y = "Número total de lecturas"
  )

g_top10


# HEATMAP DE GENEROS PRINCIPALES
# Usar los porcentajes calculados arriba para los 10 generos
datos_heatmap <- generos[generos$GrupoGenero %in% top10, ]

g_heatmap <- ggplot(
  datos_heatmap,
  aes(x = Sample, y = GrupoGenero, fill = Porcentaje)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#006D77") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Heatmap de géneros bacterianos",
    x = "Muestras", y = "Género", fill = "Abundancia (%)"
  )

g_heatmap


# RIQUEZA OBSERVADA DE ASVS
# Contar ASVs con al menos una lectura en cada muestra
riqueza <- rowSums(seqtab.nochim > 0)
riqueza

barplot(
  riqueza,
  col = "#52B69A",
  main = "Riqueza observada por muestra",
  xlab = "Muestras", ylab = "Número de ASVs"
)

# ASVS COMPARTIDAS ENTRE MUESTRAS
# Para cada ASV: contar en cuantas de las cinco muestras aparece
compartidas <- colSums(seqtab.nochim > 0)

# Numero de ASVs que aparecen en 1, 2, 3, 4 o 5 muestras
resumen_compartidas <- table(
  factor(compartidas, levels = 1:nrow(seqtab.nochim))
)
resumen_compartidas

barplot(
  resumen_compartidas,
  col = "#4361EE",
  main = "ASVs compartidas entre muestras",
  xlab = "Número de muestras donde aparece la ASV",
  ylab = "Número de ASVs"
)


# RESOLUCION DE LA ASIGNACION TAXONOMICA
# Cuantas ASVs tienen informacion en cada rango taxonomico
clasificadas <- c(
  Familia = sum(!is.na(taxa2[, "Family"])),
  Genero  = sum(!is.na(taxa2[, "Genus"])),
  Especie = sum(!is.na(taxa2[, "Species"]))
)
clasificadas

barplot(
  clasificadas,
  col = c("#48CAE4", "#0096C7", "#023E8A"),
  main = "Resolución de la asignación taxonómica",
  ylab = "Número de ASVs clasificadas"
)

# RESUMEN 

cat("ASVs:", ntaxa(ps), "| Muestras:", nsamples(ps), "\n")
print(round(retencion, 2))
print(head(ranking, 10))
print(riqueza)
print(resumen_compartidas)
print(clasificadas)


