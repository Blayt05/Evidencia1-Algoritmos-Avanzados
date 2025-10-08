
#*****************************************************************
#------- Nos sirve para leer el archivo del genoma completo ------
#*****************************************************************

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
    len = 0
    i = 1
    lps[0] = 0
    while i < M:
        if pat[i] == pat[len]:
            lps[i] = len + 1
            len += 1
            i += 1
        else:
            if len != 0:
                len = lps[len - 1]
            else:
                lps[i] = 0
                i += 1 

#********************************************************************
#--- Nos sirve para utilizar el algoritmo de string matching KMP --- 
#********************************************************************
def KMPSearch(pat, txt, nombre):
    N = len(txt)
    M = len(pat)
    lps = [0]*M
    computeLPSArray(pat, M, lps)
    i = 0
    j = 0
    while i < N:
        if txt[i] == pat[j]:
            i += 1
            j += 1
        if j == M:
            return i-j, i-1, txt[i-j:i-j + 12], nombre
            # print("El inicio del genoma ", nombre," es", i - j)
            # print(txt[i-j:i-j + 12]) #Se imprimen los primeros 12 nucleotidos del patron encontrado
            # j = lps[j-1]
        elif i < N and txt[i] != pat[j]:
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return -1, -1, -1, nombre

#**************************************************************************
#--- Nos sirve para estructurar el genoma a partir de diferentes frames --- 
#**************************************************************************   

def frames(genoma):
    proteinas_genoma = []
    #Los 3 frames disponibles
    secuencia=""
    for frame in range(3):
        for i in range(frame, len(genoma), 3):
            codon_actual = genoma[i:i+3]
            if len(codon_actual) < 3:
                break
            else:
                secuencia = secuencia + codones[codon_actual]
        proteinas_genoma.append(secuencia)
        secuencia = ""
    return proteinas_genoma

#******************************************************************************************
#--- Nos sirve para inicializar y mandar a llamar a todas nuestras funciones y archivos --- 
#******************************************************************************************
        
#********************************************************************
#---------------------- LLamadas a problema 1 ----------------------
#********************************************************************

# -------- Archivos de prueba -----------
# #Este seria el texto 
# genoma = leer_archivo("archivos/SARS-COV-2-MN908947.3.txt") # Archivo del genoma completo

# #Estos serian los patrones
# genM   = leer_archivo("archivos/gen-M.txt") # el gen M  
# genS   = leer_archivo("archivos/gen-S.txt") # el gen S
# genORF = leer_archivo("archivos/gen-ORF1AB.txt") # el gen ORF1ab
# nombre = "genM"
# nombre2 = "genS"
# nombre3 = "genORF"
# KMPSearch(genM,genoma, nombre)
# KMPSearch(genS,genoma,nombre2)
# KMPSearch(genORF,genoma,nombre3)


#********************************************************************
#---------------------- LLamadas a problema 2 ----------------------
#********************************************************************

#Tabla de aminoacidos, aqui vemos las equivalencias de tripletes de nucleotidos a aminoacidos
codones = {
    "ATG": "M",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "AAT": "N",
    "AAC": "N",
    "GAT": "D",
    "GAC": "D",
    "TGT": "C",
    "TGC": "C",
    "CAA": "Q",
    "CAG": "Q",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "CAT": "H",
    "CAC": "H",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "TTA": "L",
    "TTG": "L",
    "AAA": "K",
    "AAG": "K",
    "TTT": "F",
    "TTC": "F",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TAA": "*",
    "TGA": "*",
    "TAG": "*"
}

# #Este seria el texto 
genoma = leer_archivo("archivos/SARS-COV-2-MN908947.3.txt") # Archivo del genoma completo

resultado = frames(genoma)
i = 0
# for prot in resultado:
#     i = i + 1
#     print(i)
#     print(prot)

# for prot in resultado:
#     print("Inicio", prot["inicio"])
#     print("Final", prot["final"])
#     print("Secuencia", prot["secuencia"])
#     print("Frame", prot["frame"])
proteinas = leer_archivo_proteina("archivos/seq-proteins.txt")


# for prot_ref in proteinas:
#     i = i + 1
#     print(i)
#     print(prot_ref["nombre"])
#     print(prot_ref["secuencia"])

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
            print("")
            # print("No se encontro la proteina con el nombre: ", nombre)

#         if prot_res["secuencia"] == prot_ref["secuencia"]:
#             print("Si hay coincidencia")
#             print(prot_ref["nombre"])
#             print(prot_res["frame"])
#             print(prot_res["inicio"])
#             print(prot_res["final"])
#             print(prot_res["secuencia"])


# print(proteinas)
# print(genoma)

print(proteinas)

from algorithms import Algorithms
from pathlib import Path
from utils import Utils , FileUtils         


def run_for_file(label: str, path: Path) -> None:
    try:
        gf = FileUtils(label, str(path))        
        gf.path = str(path)       
        seq = gf.read_file()
        longest = Algorithms.mancher_algorithm(seq)
        print(f"\n[{label}] {path.name}")
        print(f"- Seq length: {len(seq)}")
        print(f"- Longest palindrome length: {len(longest)}")
        print(f"- Longest palindrome: {longest}")
    except Exception as e:
        print(f"\n[{label}] {path.name} -> ERROR: {e}")


if __name__ == "__main__":
    base = Path(__file__).parent  
    files = {
        "Gen M": base / "archivos/gen-M.txt",
        "Gen ORF1ab": base / "archivos/gen-ORF1AB.txt",
        "Gen S": base / "archivos/gen-S.txt",
    }

    for label, fpath in files.items():
        run_for_file(label, fpath)


