# El objetivo es analizar el resultado de un blastp de VFDB () contra el proteoma de Staphylococcus aureus N315. Se
# quiere obtener una tabla que muestre los resultados para cada secuencia query, incluyendo datos importantes del hit:
# id query, id hit, cobertura de la query, cobertura del hit, identidad


from Bio import SearchIO
# summarize_blastp es una función definida en el script auxiliary_functions.py,
# devuelve un dataframe de pandas con los parámetros de interés de cada hit
from auxiliary_functions import summarize_blastp
import sys

# El formato de parseo "blast-text" no funciona bien con nuevas versiones de blast (se sobre el formato por defecto)
# (>2.2.26+, https://biopython.org/wiki/SearchIO)
#blast_results = SearchIO.parse(blast_file, "blast-text")

blast_file = sys.argv[1]
blast_results = list(SearchIO.parse(blast_file, "blast-xml"))
summary=summarize_blastp(blast_results)

print("VFDB_setB_pro vs Proteome_N315")
print(summary)

# Exporto el dataframe como .tsv
summary.to_csv(sys.argv[2] + "_resumen.tsv", sep="\t", index=False)
