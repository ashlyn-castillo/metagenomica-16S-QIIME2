#!/usr/bin/env bash
# Adaptación de la práctica de Flyteek. Ejecutar en un directorio nuevo.
# Entorno: conda activate rachis-qiime2-2026.10
# Uso: RAW_DIR=/ruta/fastq CLASSIFIER=/ruta/clasificador.qza bash scripts/qiime2_v4.sh
set -euo pipefail
: "${RAW_DIR:?Indique el directorio con los diez FASTQ comprimidos originales}"
: "${CLASSIFIER:?Indique la ruta absoluta del clasificador SILVA 144 V4}"
export RAW_DIR
python - <<'MANIFEST'
import os,csv
from pathlib import Path
p=Path(os.environ['RAW_DIR']).resolve()
with open('manifest.tsv','w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 w.writerow(['sample-id','forward-absolute-filepath','reverse-absolute-filepath'])
 for i in range(1,6):
  a=p/f'M{i}_L001_R1_001.fastq.gz'; b=p/f'M{i}_L001_R2_001.fastq.gz'
  assert a.is_file() and b.is_file(), f'Falta el par M{i}'
  w.writerow([f'M{i}',str(a),str(b)])
MANIFEST
printf 'sample-id\nM1\nM2\nM3\nM4\nM5\n' > sample-metadata.tsv
qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path manifest.tsv   --output-path demux-paired-end.qza   --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize   --i-data demux-paired-end.qza   --o-visualization demux-summary.qzv

qiime dada2 denoise-paired   --i-demultiplexed-seqs demux-paired-end.qza   --p-trunc-len-f 240   --p-trunc-len-r 160   --p-trim-left-f 0   --p-trim-left-r 0   --p-max-ee-f 2   --p-max-ee-r 2   --p-n-threads 4   --o-table table-dada2.qza   --o-representative-sequences rep-seqs-dada2.qza   --o-denoising-stats denoising-stats.qza   --o-base-transition-stats base-transition-stats.qza   --verbose

qiime metadata tabulate   --m-input-file denoising-stats.qza   --o-visualization denoising-stats.qzv

qiime feature-table summarize   --i-table table-dada2.qza   --o-feature-frequencies feature-frequencies.qza   --o-sample-frequencies sample-frequencies.qza   --o-summary table-dada2.qzv

qiime feature-table tabulate-seqs   --i-data rep-seqs-dada2.qza   --o-visualization rep-seqs-dada2.qzv

qiime tools validate   "${CLASSIFIER}"

qiime feature-classifier classify-sklearn   --i-classifier "${CLASSIFIER}"   --i-reads rep-seqs-dada2.qza   --p-n-jobs 1   --o-classification taxonomy.qza   --verbose

qiime metadata tabulate   --m-input-file taxonomy.qza   --o-visualization taxonomy.qzv

qiime taxa barplot   --i-table table-dada2.qza   --i-taxonomy taxonomy.qza   --o-visualization taxa-barplot.qzv

qiime phylogeny align-to-tree-mafft-fasttree   --i-sequences rep-seqs-dada2.qza   --o-alignment aligned-rep-seqs.qza   --o-masked-alignment masked-aligned-rep-seqs.qza   --o-tree unrooted-tree.qza   --o-rooted-tree rooted-tree.qza   --verbose

qiime diversity alpha-rarefaction   --i-table table-dada2.qza   --i-phylogeny rooted-tree.qza   --p-max-depth 3700   --p-min-depth 100   --p-steps 20   --p-iterations 10   --p-metrics observed_features shannon faith_pd   --o-visualization alpha-rarefaction.qzv

qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table-dada2.qza   --p-sampling-depth 3000   --m-metadata-file sample-metadata.tsv   --output-dir core-metrics-results

qiime metadata tabulate   --m-input-file core-metrics-results/shannon_vector.qza   --o-visualization core-metrics-results/shannon-table.qzv

qiime metadata tabulate   --m-input-file core-metrics-results/faith_pd_vector.qza   --o-visualization core-metrics-results/faith-pd-table.qzv

qiime metadata tabulate   --m-input-file core-metrics-results/evenness_vector.qza   --o-visualization core-metrics-results/pielou-table.qzv

qiime metadata tabulate   --m-input-file core-metrics-results/observed_features_vector.qza   --o-visualization core-metrics-results/observed-features-table.qzv

qiime tools export   --input-path core-metrics-results/bray_curtis_distance_matrix.qza   --output-path core-metrics-results/bray-curtis-export
