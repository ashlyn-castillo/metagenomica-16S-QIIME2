# Análisis de amplicones bacterianos 16S V4 con R y QIIME 2

**Diana Ashlyn Castillo Abad — proyecto de aprendizaje, noviembre de 2026.**

La práctica se realizó durante el curso: Análisis bioinformático de metagenómica bacteriana 16S, impartido por Alex Andrés Jarrín Jurado en Flyteek. En este trabajo se describe la ejecución y adaptación de los procedimientos del curso en R y QIIME 2, este último en Ubuntu mediante WSL en Windows. Los datos corresponden a una actividad del curso y no a un estudio experimental propio.

## Objetivo y alcance

El objetivo fue procesar cinco pares de archivos FASTQ de las muestras didácticas M1–M5 mediante control de calidad, inferencia de ASVs, asignación taxonómica, reconstrucción filogenética y análisis de diversidad. La procedencia biológica de las muestras y las condiciones de redistribución del conjunto no están documentadas en los archivos disponibles. El repositorio no incluye los FASTQ ni los materiales completos del curso.

## Flujo ejecutado

FASTQ → manifiesto → importación → calidad → DADA2 → tabla y secuencias de ASVs → clasificación SILVA → barras taxonómicas → MAFFT/FastTree → rarefacción → diversidad alfa y beta → PCoA.

Parámetros DADA2: truncamiento forward/reverse 240/160; trim-left 0/0; maxEE 2/2; cuatro hilos. Clasificador recibido: `SILVA_144_SSURef_NR99_uniform_classifier_V4-515f-806r.qza`. Rarefacción exploratoria de 100 a 3700 lecturas, 20 profundidades y 10 iteraciones. Core metrics con 3000 lecturas por muestra.

## Resultados documentados

- Cinco muestras, **158 ASVs y 37 713 lecturas finales**.
- Lecturas finales M1–M5: 5015, 16822, 5482, 3715 y 6679.
- Retención global: aproximadamente 84,40 %.
- Riqueza en la tabla completa: 101, 135, 69, 87 y 91 ASVs.
- Los conteos finales coinciden con la ejecución previa de DADA2 en R. La coincidencia de los conteos no establece la equivalencia de todos los objetos o métodos.
- Bray–Curtis tras rarefacción: menor distancia M3–M5 (~0,277), mayor M1–M3 (~0,483).

`resultados/seguimiento_lecturas.csv` transcribe los conteos verificados en las capturas de la ejecución; no es una exportación automática del artefacto. Las exportaciones realizadas desde R se encuentran en `resultados/R/`: seguimiento de lecturas, diversidad alfa, riqueza, taxonomía, cuatro matrices de distancia, cinco figuras PNG, un árbol PDF y las versiones de R. Los artefactos QIIME 2 no están incluidos.

## Reproducción

1. Se crea o activa el entorno QIIME 2 2026.10. El archivo YAML registra el entorno Linux utilizado; la instalación puede variar entre sistemas.
2. Se obtienen, por una vía autorizada, los pares originales comprimidos, con nombres `M1_L001_R1_001.fastq.gz` hasta `M5_L001_R2_001.fastq.gz` y el clasificador indicado.
3. Se crea un directorio de trabajo vacío y, desde él, se ejecuta el script con las variables `RAW_DIR` y `CLASSIFIER` definidas mediante rutas absolutas. Los archivos de entrada se conservan.

```bash
conda activate rachis-qiime2-2026.10
RAW_DIR=/ruta/a/fastq CLASSIFIER=/ruta/al/clasificador.qza bash /ruta/al/repositorio/scripts/qiime2_v4.sh
```

El script organiza comandos observados en el historial. Se comprobó su sintaxis con `bash -n`; **esta versión reorganizada no se ejecutó de extremo a extremo**. La rarefacción aleatoria puede producir cifras diferentes al repetirla. Las versiones se registran en `reproducibilidad/`.

## Interpretación y límites

Una ASV no equivale necesariamente a una especie. Los conteos de lecturas no representan cantidades de células ni permiten determinar viabilidad, actividad o función. No se dispone de grupos experimentales documentados para evaluar efectos de tratamientos. La taxonomía de esta ejecución usa SILVA 144 y no es directamente intercambiable con SILVA 138.2 utilizada en R. El árbol de QIIME 2 se construyó con MAFFT/FastTree y difiere del flujo filogenético de R.

## Estado y próximos pasos

La práctica V4 se ejecutó y las exportaciones de R se incorporaron al repositorio. Queda pendiente comprobar los scripts reorganizados mediante una ejecución completa desde una sesión limpia. Las condiciones de reutilización del material docente no están documentadas. El ejercicio independiente V3–V4 no forma parte de estos resultados. Los scripts de R se organizaron en cinco etapas (01–05). Los bloques de árbol y UniFrac se reconstruyeron a partir de las instrucciones y salidas de la sesión. El modelo de R utilizó una categoría de tasas (sin gamma); no se calculó bootstrap.

## Reconocimientos y referencias

- Curso Flyteek y guía docente `QIIME2_16S_V4_pipeline.pdf`.
- QIIME 2: https://qiime2.org/
- DADA2: https://benjjneb.github.io/dada2/
- SILVA: https://www.arb-silva.de/

No se asigna una licencia de redistribución a material de terceros cuyo permiso no esté documentado.

## Ruta en R

Se definen `RAW_DIR` (FASTQ sin comprimir), `SILVA_TRAIN` y `SILVA_SPECIES` mediante `Sys.setenv()`. Los scripts 01, 02 y 03 se ejecutan en ese orden dentro de una misma sesión. El bloque 04 reconstruye el árbol con el modelo GTR y una categoría de tasas, y calcula las distancias UniFrac. El bloque 05 exporta objetos ya calculados sin repetir ese procesamiento. Las etiquetas ASV abreviadas del árbol son identificadores de presentación, no nombres de especies.

Se ajustaron los patrones FASTQ para excluir las copias `.gz`, se retiró el bloque que usaba `ps` antes de construirlo y se separó la filogenia del pipeline inicial. No se incorporan categorías Control/Tratado, dado que no se dispone de un diseño experimental que las respalde. Estas adaptaciones requieren validación en una ejecución limpia; la evidencia disponible corresponde a los resultados observados durante la práctica guiada.

## Figuras exportadas en R

![Composición por filos](resultados/R/g_filos.png)

![Ordenación PCoA con Bray–Curtis](resultados/R/PCoA_Bray_Curtis.png)

Otras figuras: [géneros](resultados/R/g_generos.png), [taxones principales](resultados/R/g_top10.png), [mapa de calor](resultados/R/g_heatmap.png) y [árbol de 20 ASVs](resultados/R/arbol_20_ASVs.pdf).

La diversidad alfa exportada desde R se calculó sobre la tabla completa; el índice de Shannon utiliza el logaritmo natural. QIIME 2 utiliza aquí una tabla rarefactada a 3000 lecturas y Shannon en base 2; los valores no deben compararse directamente. El árbol en R utiliza GTR, una categoría de tasas, sin gamma ni bootstrap, y enraizamiento por punto medio.
