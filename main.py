#*****************************************************************
#------- Nos sirve para leer el archivo del genoma completo ------
#*****************************************************************
from typing import List, Dict, Tuple
#Funcion que nos ayuda a normalizar nuestro genomas y genes en una sola linea sin espacios
def leer_archivo(ruta):
    seq = []
    with open(ruta, "r", encoding="utf-8") as f: #Se abre el archivo fasta y se cierra con with
        for line in f: #Para cada linea del archivo 
            if line.startswith(">"): # Si empieza con > 
                continue  # Se ignora la cabecera
            seq.append(line.strip())  # Se quitan saltos de línea y espacios
    return "".join(seq) # Se unen todos los arreglos en uno solo sin espacios


#*****************************************************************
#------- Nos sirve para leer el archivo de la proteinas completo 
#*****************************************************************

#Funcion que nos ayuda a formatear de manear correcta el nombre de la proteinas y su secuencia en aminoacidos
def leer_archivo_proteina(ruta):
    lista_proteinas = []
    secuencia_actual = ""
    with open(ruta, "r", encoding="utf-8") as f: #Se abre el archivo fasta y se cierra con with
        for line in f: #Para cada linea del archivo 
            if line.startswith(">"): # Si empieza con > 
                nombre_actual = line[1:].strip()  # Se toma el nombre
            else:
                secuencia_actual = line.strip()# Se quitan saltos de línea y espacios
                diccionario_proteina = {
                    "nombre": nombre_actual,
                    "secuencia": secuencia_actual
                }
                lista_proteinas.append(diccionario_proteina)
    return lista_proteinas # Se unen todos los arreglos en uno solo sin espacios

#*****************************************************************
#------- Nos sirve para crear el arreglo LPS --------------------- 
#*****************************************************************

def computeLPSArray(pat, M, lps):
    #Se inicializo la variable que guarda la longitud del prefijo mas largo
    len = 0
    i = 1
    lps[0] = 0 #El primero elemento siempre empieza en 0
    #Se recorre el patron para calcular el arreglo lps
    while i < M:
        if pat[i] == pat[len]:
            #Si hay coincidencias entonces aumento la longitud y guardo el valor en el lps
            lps[i] = len + 1
            len += 1
            i += 1
        else:
            if len != 0:
                #Si no hay coincidencia y len no es cero, se actualiza len usando el valor anteirior del lps
                len = lps[len - 1]
            else:
                #Si len es cero entonces pongo cero en el lps y sigo avanzando
                lps[i] = 0
                i += 1 

#********************************************************************
#--- Nos sirve para utilizar el algoritmo de string matching KMP --- 
#********************************************************************
def KMPSearch(pat, txt, nombre):
    #Se guarda la longitud del texto en N
    N = len(txt) 
    #Se guarda la longitud del patron
    M = len(pat)
    #Se crea un lps del tamaño del patron
    lps = [0]*M
    #Se consigue el lps
    computeLPSArray(pat, M, lps)
    #Indice que nos sirve para el texto
    i = 0
    #Indice que nos sirve para el patron
    j = 0
    #Se recorre el texto buscando el patron
    while i < N:
        if txt[i] == pat[j]:
            #Si hay coincidencia se avanza al siguiente caracter en ambos indices
            i += 1
            j += 1
        if j == M:
            #Si se encuentra la misma longitud de j que la del patron entonces se retorna el indice inicio, final, secuencia y nombre
            return i-j, i-1, txt[i-j:i-j + 12], nombre
        elif i < N and txt[i] != pat[j]:
            if j != 0:
                #Si hay un error y j no es cero, se actualiza j usando el lps
                j = lps[j-1]
            else:
                #Si j es cero, solo se avanza en el texto
                i += 1
    #Si no se encontro el patron, se regresa -1 por seguridad
    return -1, -1, -1, nombre

#**************************************************************************
#--- Nos sirve para estructurar el genoma a partir de diferentes frames --- 
#**************************************************************************   

def frames(genoma):
    #Arreglo que nos sirve para guardar las 3 secuencias del genoma traducido en los diferentes tipos de marcos
    proteinas_genoma = []
    #Se almacena la secuencia de cada frame aqui
    secuencia=""
    #Los 3 diferentes de frames se representan con un for
    for frame in range(3):
        #Se va desde el indice 0 al fin del genoma en pasos de 3
        for i in range(frame, len(genoma), 3):
            #El codon actual seria conformado por 3 nucleotidos
            codon_actual = genoma[i:i+3]
            #Si el codon es menor a 3 entonces se rompe el ciclo ahi, esto se hace por irregularidades que los marcos provoquen
            if len(codon_actual) < 3:
                break
            #Si no entonces se siguen añadiendo codones a la secuencia 
            else:
                secuencia = secuencia + codones[codon_actual]
        #Una vez se termina un frame se agrega al arreglo y se continua con el siguiente frame
        proteinas_genoma.append(secuencia)
        #Se limpia la secuencia
        secuencia = ""
    #Al terminar todos los frames se retorna el arreglo
    return proteinas_genoma

#**************************************************************************
#------------------- Funciones para el algoritmo de Manacher --------------
#**************************************************************************   

def format_manchester(s: str) -> str:
        """
        Function: Format list
        Pourpose: Give format to the list in order to be processed by
        Manchester alglrithm in algorihms.py file.
        Complexity: O(n)
        """
        simbol_1 = '@'
        simbol_2 = '#'
        simbol_3 = '$'
        char = 0

        result = [simbol_1]
        new_size = 2 * len(s) + 3

        for i in range(1, new_size - 1):
            if i % 2 != 0:
                result.append(simbol_3)
            else:
                if char < len(s):
                    result.append(s[char])
                    char += 1

        result.append(simbol_2)
        return result

def mancher_algorithm(s: str) -> str:
        """
        Function: Manacher's Algorithm
        Purpose: Find the largest palindrome substring in a given string.
        Complexity: O(n), where n is the length of the input string.
        """
        simbol_1 = '@'
        simbol_2 = '#'
        simbol_3 = '$'

        new_list = format_manchester(s)  # formateador
        n = len(new_list)
        if n == 0:
            return ""

        P = [0] * n
        center = 0
        limit = 0

        for i in range(1, n - 1):
            if i < limit:
                simetric = 2 * center - i
                P[i] = min(limit - i, P[simetric])

            gap = P[i] + 1
            # expandir cuidando los límites
            while i - gap >= 0 and i + gap < n and new_list[i - gap] == new_list[i + gap]:
                P[i] += 1
                gap += 1

            if i + P[i] > limit:
                limit = i + P[i]
                center = i

        # obtener el mejor centro
        best_center = max(range(n), key=lambda k: P[k])
        best_radius = P[best_center]
        L = best_center - best_radius
        R = best_center + best_radius

        # reconstruir
        out = []
        for k in range(L, R + 1):
            ch = new_list[k]
            if ch not in (simbol_1, simbol_2, simbol_3):
                out.append(ch)

        return "".join(out)


def run_palindrome_for_file(label: str, file_path: str) -> None:
    try:
        seq = leer_archivo(file_path)
        longest = mancher_algorithm(seq)
        print(f"\n[{label}] {file_path}")
        print(f"- Longitud de la secuencia: {len(seq)}")
        print(f"- Longitud del palindromo mas largo: {len(longest)}")
        print(f"- Palindromo mas largo: {longest}")
    except Exception as e:
        print(f"\n[{label}] {file_path} -> ERROR: {e}")

#******************************************************************************************
#--- Nos sirve para inicializar y mandar a llamar a todas nuestras funciones y archivos --- 
#******************************************************************************************
        
#********************************************************************
#---------------------- LLamadas a problema 1 ----------------------
#********************************************************************

print("********************************************************************")
print("---------------------- Comienzo a problema 1 ----------------------")
print("********************************************************************")
# -------- Archivos de prueba -----------
# #Este seria el texto 
genoma = leer_archivo("archivos/SARS-COV-2-MN908947.3.txt") # Archivo del genoma completo

# #Estos serian los patrones
genes = [
    (leer_archivo("archivos/gen-M.txt"), "genM"),
    (leer_archivo("archivos/gen-S.txt"), "genS"),
    (leer_archivo("archivos/gen-ORF1AB.txt"), "genORF")
]

for gen, nombre in genes:
    
    inicio, final, secuencia, nombre = KMPSearch(gen,genoma,nombre)

    print("El nombre del gen es: ", nombre)
    print("El inicio es:", inicio)
    print("El final es:", final)
    print("El inicio es:", secuencia)


#********************************************************************
#---------------------- LLamadas a problema 2 ----------------------
#********************************************************************

print("********************************************************************")
print("---------------------- Comienzo a problema 2 ----------------------")
print("********************************************************************")
         
# Archivos para analizar palindromos
palindrome_files = {
    "Gen M": "archivos/gen-M.txt",
    "Gen ORF1ab": "archivos/gen-ORF1AB.txt",
    "Gen S": "archivos/gen-S.txt",
}

for label, file_path in palindrome_files.items():
    run_palindrome_for_file(label, file_path)

#********************************************************************
#---------------------- LLamadas a problema 3 ----------------------
#********************************************************************

print("********************************************************************")
print("---------------------- Comienzo problema 3 ----------------------")
print("********************************************************************")
#Tabla de aminoacidos, aqui vemos las equivalencias de tripletes de nucleotidos a aminoacidos
codones = {
    "ATG": "M", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "AAT": "N", "AAC": "N", "GAT": "D", "GAC": "D",
    "TGT": "C", "TGC": "C", "CAA": "Q", "CAG": "Q",
    "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H", "ATT": "I", "ATC": "I", "ATA": "I",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L",
    "AAA": "K", "AAG": "K", "TTT": "F", "TTC": "F",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGG": "W", "TAT": "Y", "TAC": "Y",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TAA": "*", "TGA": "*", "TAG": "*"
}

resultado = frames(genoma)
proteinas = leer_archivo_proteina("archivos/seq-proteins.txt")
for frame, prot_res in enumerate(resultado):
    for prot_ref in proteinas:
        inicio, final, secuencia, nombre= KMPSearch(prot_ref["secuencia"],prot_res,prot_ref["nombre"])
        if inicio != -1:
            print ("El nombre de la protiena es:", nombre)
            print("Frame:", frame)
            print("El inicio del aminoacido en el genoma es en el indice: ", inicio * 3 + frame)
            print("El final del aminoacido en el genoma es en el indice: ", final * 3 + frame)
            primeros_4_codones = genoma[inicio*3: inicio*3 + 12]
            print("La secuencia de los 4 codones de la proteina es: ", primeros_4_codones)
            primeros_4_aminoacidos = prot_ref["secuencia"][:4]
            print("La secuencia de los 4 aminoacidos de la proteina es: ", primeros_4_aminoacidos)
            print("**************************")
        else:
            pass
            # print("No se encontro la proteina con el nombre: ", nombre)

#Llamada slippery sequence
# slippery_seq = "TTTAAAC"
# shift = +1
# proteinas_afectada
# print()

# ======== PROBLEMA 4: Wuhan (2019) vs Texas (2020) ==========================
# ============================================================================

#********************************************************************
#---------------------- Comienzo a problema 4 ----------------------
#********************************************************************
# Wuhan (2019) vs Texas (2020): iguales, diferencias y efecto en AA
from typing import List, Dict, Tuple

# ---------- Utilidades específicas del P4 (reutiliza 'codones') ----------
def aa_from_codon_p4(codon: str) -> str:
    """Traduce 1 codón usando tu tabla global 'codones'."""
    return codones.get(codon, "X") if len(codon) == 3 else "X"

def comparar_nt_p4(a: str, b: str) -> Tuple[List[Dict], List[Tuple[int, int]]]:
    """
    Devuelve:
      - difs: [{'pos': i, 'a': 'A', 'b': 'G'}, ...]
      - iguales_runs: tramos contiguos iguales [(ini, fin_inclusivo), ...]
    """
    n = min(len(a), len(b))
    difs: List[Dict] = []
    iguales: List[Tuple[int, int]] = []
    i = 0
    while i < n:
        if a[i] == b[i]:
            j = i
            while j < n and a[j] == b[j]:
                j += 1
            iguales.append((i, j - 1))
            i = j
        else:
            difs.append({"pos": i, "a": a[i], "b": b[i]})
            i += 1
    return difs, iguales

def cambios_codones_p4(difs: List[Dict], ref: str, alt: str, frame: int = 0) -> List[Dict]:
    """
    Anota cada diferencia con el codón y AA, en el marco indicado (0/1/2).
    - Si la diferencia cae en un codón incompleto (al final), se omite.
    """
    out = []
    for d in difs:
        pos = d["pos"]
        # Ajuste de codón según frame (marco de lectura)
        c0 = ((pos - frame) // 3) * 3 + frame
        if c0 < frame or c0 + 2 >= len(ref) or c0 + 2 >= len(alt):
            continue
        codon_ref = ref[c0:c0 + 3]
        codon_alt = alt[c0:c0 + 3]
        aa_ref = aa_from_codon_p4(codon_ref)
        aa_alt = aa_from_codon_p4(codon_alt)
        out.append({
            "pos_nt": pos,
            "nt_ref": d["a"],
            "nt_alt": d["b"],
            "codon_start": c0,
            "frame": frame,
            "codon_ref": codon_ref,
            "codon_alt": codon_alt,
            "aa_ref": aa_ref,
            "aa_alt": aa_alt,
            "tipo": "no-sinonimo" if aa_ref != aa_alt else "sinonimo"
        })
    return out

# ---------- Impresiones bonitas (muestras para no inundar consola) ----------
def print_tramos_iguales_p4(tramos: List[Tuple[int, int]], max_show: int = 10):
    print("\nTramos IGUALES (muestras):")
    if not tramos:
        print("  (no se encontraron tramos iguales)")
        return
    for k, (s, e) in enumerate(tramos[:max_show], 1):
        print(f"  {k:02d}. nt[{s}:{e}] (len={e - s + 1})")
    if len(tramos) > max_show:
        print(f"  ... {len(tramos) - max_show} tramos adicionales no mostrados")

def print_difs_p4(difs: List[Dict], max_show: int = 20):
    print("\nDIFERENCIAS nucleotídicas (muestras):")
    for k, d in enumerate(difs[:max_show], 1):
        print(f"  {k:02d}. nt {d['pos']}: {d['a']}→{d['b']}")
    if len(difs) > max_show:
        print(f"  ... {len(difs) - max_show} diferencias adicionales no mostradas")

def print_cambios_codones_p4(rows: List[Dict], max_show: int = 50):
    print("\nCambios de CODÓN / AMINOÁCIDO (muestras):")
    for k, r in enumerate(rows[:max_show], 1):
        print(
            f"  {k:02d}. [frame {r['frame']}] nt {r['pos_nt']}: {r['nt_ref']}→{r['nt_alt']} | "
            f"{r['codon_ref']}({r['aa_ref']}) → {r['codon_alt']}({r['aa_alt']}) "
            f"| codón@{r['codon_start']}-{r['codon_start']+2} | {r['tipo']}"
        )
    if len(rows) > max_show:
        print(f"  ... {len(rows) - max_show} cambios adicionales no mostrados")

# ---------- Ejecución del Punto 4 ----------
genoma_wuhan_p4 = leer_archivo("archivos/SARS-COV-2-MN908947.3.txt")
genoma_texas_p4 = leer_archivo("archivos/SARS-COV-2-MT106054.1.txt")

print("====================================================================")
print("PUNTO 4: Comparación Wuhan (2019) vs Texas (2020)")
print("====================================================================")
print(f"Longitudes: Wuhan={len(genoma_wuhan_p4):,} nt | Texas={len(genoma_texas_p4):,} nt")

# ¿Dónde son iguales / dónde difieren?
difs_p4, iguales_runs_p4 = comparar_nt_p4(genoma_wuhan_p4, genoma_texas_p4)
comp_len_p4 = min(len(genoma_wuhan_p4), len(genoma_texas_p4))
identidad_p4 = (comp_len_p4 - len(difs_p4)) / comp_len_p4 * 100 if comp_len_p4 else 0.0

print(f"\n== RESUMEN NUCLEÓTIDOS ==")
print(f"Comparación hasta: {comp_len_p4:,} nt")
print(f"Iguales: {comp_len_p4 - len(difs_p4):,}  |  Diferencias: {len(difs_p4):,}")
print(f"Identidad: {identidad_p4:.2f}%")

print_tramos_iguales_p4(iguales_runs_p4, max_show=10)
print_difs_p4(difs_p4, max_show=20)

# ¿Las diferencias resultan en aminoácidos diferentes?
# Por defecto, mostramos el marco 0 (el más común en estos ejercicios)
cambios_frame0 = cambios_codones_p4(difs_p4, genoma_wuhan_p4, genoma_texas_p4, frame=0)
sinonimos_0 = sum(1 for r in cambios_frame0 if r["tipo"] == "sinonimo")
no_sinonimos_0 = sum(1 for r in cambios_frame0 if r["tipo"] == "no-sinonimo")

print("\n== IMPACTO EN AMINOÁCIDOS (marco 0) ==")
print(f"Total diferencias evaluadas en codones completos: {len(cambios_frame0):,}")
print(f"  • Sinónimas (mismo AA): {sinonimos_0:,}")
print(f"  • No sinónimas (AA distinto): {no_sinonimos_0:,}")
print_cambios_codones_p4(cambios_frame0, max_show=50)

