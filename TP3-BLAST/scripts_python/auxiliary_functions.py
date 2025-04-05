import pandas as pd

def summarize_blastp(blast_results,include_variants=0):
    
    # Inicializo listas que contendrán los parámetros importantes del blast y serán las columnas del dataframe. 
    id_querys=[]
    id_hits=[]
    coverage_querys=[]
    coverage_hits=[]
    p_ident=[]
    variants=[]

    # Itero sobre cada secuencia query de VFDB
    for query in blast_results:
        query_len=query.seq_len
        # Compruebo que la query tenga hits
        if query.hits:
            # Itero sobre cada hit, extrayendo la ID de la secuencia subject y el nombre de la proteína.
            for hit in query.hits:
                hit_len=hit.seq_len
                hit_id=hit.id.split('|')[1] if '|' in hit.id else hit.id
                # Itero sobre cada hsp (high scoring-sequence pair) y añado a las listas los datos del
                # hsp que quiero en las columnas del dataframe. (Lo que llamo coverage_hits es en realidad
                # la covertura de cada hsp)
                for hsp in hit.hsps:
                    id_querys.append(query.id.split('|')[1] if '|' in query.id else query.id)
                    id_hits.append(hit_id)
                    coverage_querys.append((hsp.aln_span/query_len)*100)
                    coverage_hits.append((hsp.aln_span/hit_len)*100)
                    p_ident.append((hsp.ident_num / hsp.aln_span) * 100)

                    if include_variants:
                        variants.append(find_variants(hsp.query.seq,hsp.hit.seq))

        # Si no se encuentran hits para la query añado su ID y agrego 0 al resto de listas.         
        else:  
            id_querys.append(query.id)
            id_hits.append(0) 
            coverage_querys.append(0)
            coverage_hits.append(0)
            p_ident.append(0)
            if include_variants:
                variants.append(0)

    # Genero el dataframe 
    summary_blastp = pd.DataFrame({
        "query_id": id_querys,
        "hit_id": id_hits,
        "cobertura_query": coverage_querys,
        "cobertura_hit": coverage_hits,
        "porcentaje_identidad": p_ident
    })

    if include_variants:
        summary_blastp['variantes'] = variants

    return summary_blastp

aa_code = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }

def find_variants(q,s):
    global aa_code
    variants=''
    last_aa,pos='',0
    for i,aa in enumerate(s):
        if aa in ['X','*'] or q[i] in ['X','*']:
            continue
        if aa!='-':
            last_aa,pos=aa,i+1
        if aa != q[i]:
            if aa == '-':
                prev_aa=s[:i].replace('-','')[-1]
                variants+=f'| {aa_code[last_aa]}{pos}_{i+1}ins{aa_code[q[i]]} '
            elif q[i] == '-':
                variants+=f'| {aa_code[aa]}{i+1}del '
            else:
                variants+=f'| {aa_code[aa]}{i+1}{aa_code[q[i]]} '
    return variants[2:-1]
