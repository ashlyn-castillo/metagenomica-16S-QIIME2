# Adaptación de la práctica Flyteek de Alex Andrés Jarrín Jurado.
# Versión reorganizada para el portafolio; no ejecutada de extremo a extremo.
# Ejecutar después de 02_graficos_R.R. Conserva ps y sus identificadores.
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(phyloseq)
seqs <- colnames(seqtab.nochim)
names(seqs) <- seqs
stopifnot(setequal(names(seqs), taxa_names(ps)))
alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor = NA)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
# La ejecución documentada utilizó k=1: optGamma fue ignorado.
# Se explicita esa configuración para no afirmar que se ajustó gamma.
fit <- pml(treeNJ, data = phang.align, k = 1)
set.seed(123)
fitGTR <- optim.pml(fit, model = "GTR", optInv = TRUE,
                   optGamma = FALSE, rearrangement = "stochastic",
                   control = pml.control(trace = 1))
tree <- fitGTR$tree
stopifnot(setequal(tree$tip.label, taxa_names(ps)))
tree_rooted <- phangorn::midpoint(tree)
ps_tree <- merge_phyloseq(ps, tree_rooted)
stopifnot(ape::is.rooted(phy_tree(ps_tree)), ntaxa(ps_tree) == ntaxa(ps))
d_unifrac <- UniFrac(ps_tree, weighted = FALSE)
d_wunifrac <- UniFrac(ps_tree, weighted = TRUE, normalized = TRUE)
round(as.matrix(d_unifrac), 3)
round(as.matrix(d_wunifrac), 3)
top20 <- head(names(sort(taxa_sums(ps_tree), decreasing = TRUE)), 20)
ps_top20 <- prune_taxa(top20, ps_tree)
arbol_vista <- phy_tree(ps_top20)
arbol_vista$tip.label <- paste0("ASV", match(arbol_vista$tip.label, taxa_names(ps_tree)))
plot(arbol_vista, show.node.label = FALSE, cex = 0.7,
     main = "Árbol de las 20 ASVs con más lecturas")
ape::add.scale.bar()
