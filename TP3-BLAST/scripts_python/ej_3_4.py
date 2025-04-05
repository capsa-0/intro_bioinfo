# El objetivo es analizar el resultado de un blastp del proteoma de Staphylococcus aureus N315 contra una base de
# enterotoxinas curada. Se quiere obtener una tabla con información sobre los mejores hits, eligiéndolos en función 
# de la cobertura de la secuencia query y la identidad. Además de algunos parámetros que devuelve blastp, se quieren 
# encontrar variantes/polimorfismos usando como referencia la secuencia de los hits

from Bio import SearchIO
from auxiliary_functions import summarize_blastp
import sys

blast_file = sys.argv[1]
blast_results = list(SearchIO.parse(blast_file, "blast-xml"))

# La base de datos de enterotoxinas curadas contiene secuencias distintas con el mismo ID e.g.

#>Gr11/seh-2p - truncated NC_017351.1
#NIILFFTLLFVLTSYAKAEDLHDKSELSDLALANAYGQYNHPFIKENVISSEIIGEKDLIFRNQGDNGNDLR
#VKFASAGLAQNFKNKNVDIYGASFYYKCEKVSENISECLYGGTTLNNEKLEQERVIGANVWVNGIQKET
#ELIRTNKKNVTLQELDIKMRKNYLINIEFIIKTVK
#>Gr11/seh-2p - fully repaired NC_017351.1
#MRKKIRILFXFFTLLFVLTSYAKAEDLHDKSELSDLALANAYGQYNHPFIKENVISSEIIGEKDLIFRNQGD
#NGNDLRVKFASAGLAQNFKNKNVDIYGASFYYKCEKVSENISECLYGGTTLNNEKLEQERVIGANVWV
#NGIQKETELIRTNKKNVTLQELDIKMRKXLSNKYRIYYKDSEIRKGLIEFDMKTPRDYSFDIYDLKGENDY
#EIDKIYEDNKTLKSEDISHIDIYLYTK

# En este caso se trata de la misma proteína pero una es la versión truncada. Como "truncated" está separada de seh-2p,
# al parsear, se toma como ID "seh-2p", lo que eleva un warning; Igualmente SearchIO simplemente asigna una nueva ID. 
# Se podría usar un sufijo de la forma "_truncated" en la base de datos (seh-2p_truncated) para solucionar ésto.

# Uso include_variants=1, para que el dataframe incluya las variantes encontradas en cada alineamiento
# La función para encontrar variables está definida en auxiliary_functions.py
summary = summarize_blastp(blast_results,include_variants=1)

# Parámetros de cobertura e identidad mínimos establecidos desde la terminal
min_cov_query=float(sys.argv[3])
min_ident=float(sys.argv[4])

# Filtro el dataframe, quedándome solamente con aquellos hits que cumplan los criterios establecidos
summary_filtered = summary.loc[
    (summary["cobertura_query"] > min_cov_query) &
    (summary["porcentaje_identidad"] > min_ident)
]
print("Proteome_N315 vs curated_enterotoxins")
print(summary_filtered)

# Exporto el dataframe como .tsv
summary_filtered.to_csv(sys.argv[2] + "_resumen_filtrado.tsv", sep="\t", index=False)