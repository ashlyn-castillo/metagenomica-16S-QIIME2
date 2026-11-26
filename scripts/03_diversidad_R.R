# Adaptación de la práctica Flyteek de Alex Andrés Jarrín Jurado.
# Versión reorganizada para el portafolio; no ejecutada de extremo a extremo.
library(phyloseq)
library(ggplot2)
alfa <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon", "Simpson")
)

alfa$ExpShannon <- exp(alfa$Shannon)
alfa$InvSimpson <- 1 / (1 - alfa$Simpson)

round(alfa, 2)
library(vegan)

mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) mat <- t(mat)

rarecurve(
  mat, step = 100, label = TRUE,
  col = 1:nrow(mat), lwd = 2,
  xlab = "Lecturas acumuladas",
  ylab = "ASVs esperadas"
)
prof_comun <- min(rowSums(mat))
prof_comun

riqueza_igual <- rarefy(mat, sample = prof_comun)

cbind(
  Observada = rowSums(mat > 0),
  Rarefaccion = round(riqueza_igual, 1)
)
# Jaccard: presencia y ausencia
ps_pa <- transform_sample_counts(
  ps, function(x) {
    x[] <- as.numeric(x > 0)
    x
  }
)

jaccard <- phyloseq::distance(
  ps_pa, method = "jaccard"
)

round(as.matrix(jaccard), 2)

# Bray-Curtis: abundancias relativas
ps_rel <- transform_sample_counts(
  ps, function(x) x / sum(x)
)
bray <- phyloseq::distance(ps_rel, method = "bray")
round(as.matrix(bray), 2)
ord <- ordinate(
  ps_rel,
  method = "PCoA",
  distance = bray
)

coord <- as.data.frame(ord$vectors[, 1:2])
colnames(coord) <- c("Eje1", "Eje2")
coord$Muestra <- rownames(coord)

ggplot(coord, aes(x = Eje1, y = Eje2)) +
  geom_point(size = 4) +
  geom_text(aes(label = Muestra), vjust = -0.8) +
  theme_minimal() +
  labs(
    title = "PCoA — Bray-Curtis",
    x = "Eje 1",
    y = "Eje 2"
  )