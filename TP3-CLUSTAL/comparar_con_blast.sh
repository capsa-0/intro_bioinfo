#!/bin/bash

# Crear archivos de salida con encabezado (celda vacía al principio)
echo -ne "\t" > bit_score_matrix.tsv
echo -ne "\t" > positives_matrix.tsv

# Obtener todos los archivos .fasta en fastas/
files=(fastas/*.fasta)
names=()

# Extraer nombres sin ruta ni extensión
for f in "${files[@]}"; do
    name=$(basename "$f" .fasta)
    names+=("$name")
done

# Escribir nombres de las columnas
for name in "${names[@]}"; do
    echo -ne "${name}\t" >> bit_score_matrix.tsv
    echo -ne "${name}\t" >> positives_matrix.tsv
done

echo >> bit_score_matrix.tsv
echo >> positives_matrix.tsv

# Comparaciones por pares
for i in "${!files[@]}"; do
    query="${files[i]}"
    row_name="${names[i]}"
    
    echo -ne "${row_name}\t" >> bit_score_matrix.tsv
    echo -ne "${row_name}\t" >> positives_matrix.tsv

    for j in "${!files[@]}"; do
        subject="${files[j]}"

        if [ "$i" -eq "$j" ]; then
            echo -ne "-\t" >> bit_score_matrix.tsv
            echo -ne "-\t" >> positives_matrix.tsv
        else
            result=$(blastp -query "$query" -subject "$subject" -matrix BLOSUM45 \
                     -outfmt "6 bitscore positive length" 2>/dev/null)

            # Extraer valores (primera línea)
            bit=$(echo "$result" | awk '{print $1}' | head -n1)
            positives=$(echo "$result" | awk '{print $2}' | head -n1)
            length=$(echo "$result" | awk '{print $3}' | head -n1)

            # Calcular porcentaje de positivos
            if [[ -n "$positives" && -n "$length" && "$length" -ne 0 ]]; then
                percent=$(awk -v p="$positives" -v l="$length" 'BEGIN { printf "%.1f", (p / l) * 100 }')
            else
                percent=0
            fi

            echo -ne "${bit:-0}\t" >> bit_score_matrix.tsv
            echo -ne "${percent:-0}\t" >> positives_matrix.tsv
        fi
    done

    echo >> bit_score_matrix.tsv
    echo >> positives_matrix.tsv
done
