
#Funcion que nos ayuda a normalizar nuestro genomas y genes en una sola linea sin espacios
def leer_archivo(ruta):
    seq = []
    with open(ruta, "r", encoding="utf-8") as f: #Se abre el archivo fasta y se cierra con with
        for line in f: #Para cada linea del archivo 
            if line.startswith(">"): # Si empieza con > 
                continue  # Se ignora la cabecera
            seq.append(line.strip())  # Se quitan saltos de línea y espacios
    return "".join(seq) # Se unen todos los arreglos en uno solo sin espacios

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
            print("El inicio del genoma ", nombre," es", i - j)
            print(txt[i-j:i-j + 12]) #Se imprimen los primeros 12 nucleotidos del patron encontrado
            j = lps[j-1]
        elif i < N and txt[i] != pat[j]:
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
            

def frames(genoma):
    proteinas = []
    #Frame 0 de 3 en 3
    estado = "fuera"
    secuencia = ""
    inicio = -1
    for i in range(0, len(genoma), 3):
        codon_actual = genoma[i:i+3]
        
        if estado == "fuera":
            if codon_actual == "ATG":
                inicio = i
                estado = "dentro"
                secuencia = "M"

        if estado == "dentro":
            if codon_actual == "TAA" or codon_actual == "TGA" or codon_actual == "TAG":
                final = i
                proteinas.append({"inicio": inicio,
                                  "final": final,
                                  "frame": 0,
                                  "secuencia": secuencia
                                })
                
                estado = "fuera"
            else:
                if codon_actual in codones:
                    secuencia = secuencia + codones[codon_actual]

    #Frame 1 de 3 en 3
    estado = "fuera"
    secuencia = ""
    inicio = -1
    for i in range(1, len(genoma), 3):
        codon_actual = genoma[i:i+3]
        
        if estado == "fuera":
            if codon_actual == "ATG":
                inicio = i
                estado = "dentro"
                secuencia = "M"

        if estado == "dentro":
            if codon_actual == "TAA" or codon_actual == "TGA" or codon_actual == "TAG":
                final = i
                proteinas.append({"inicio": inicio,
                                "final": final,
                                "frame": 1,
                                "secuencia": secuencia
                                })
                
                estado = "fuera"
            else:
                if codon_actual in codones:
                    secuencia = secuencia + codones[codon_actual]

    #Frame 2 de 3 en 3
    estado = "fuera"
    secuencia = ""
    inicio = -1
    for i in range(2, len(genoma), 3):
        codon_actual = genoma[i:i+3]
        
        if estado == "fuera":
            if codon_actual == "ATG":
                inicio = i
                estado = "dentro"
                secuencia = "M"

        if estado == "dentro":
            if codon_actual == "TAA" or codon_actual == "TGA" or codon_actual == "TAG":
                final = i
                proteinas.append({"inicio": inicio,
                                "final": final,
                                "frame": 2,
                                "secuencia": secuencia
                                })
                
                estado = "fuera"
            else:
                if codon_actual in codones:
                    secuencia = secuencia + codones[codon_actual]
    return proteinas




        

# -------- Archivos de prueba -----------
#Este seria el texto 
genoma = leer_archivo("archivos/SARS-COV-2-MN908947.3.txt") # Archivo del genoma completo

#Estos serian los patrones
genM   = leer_archivo("archivos/gen-M.txt") # el gen M  
genS   = leer_archivo("archivos/gen-S.txt") # el gen S
genORF = leer_archivo("archivos/gen-ORF1AB.txt") # el gen ORF1ab
nombre = "genM"
nombre2 = "genS"
nombre3 = "genORF"
KMPSearch(genM,genoma, nombre)
KMPSearch(genS,genoma,nombre2)
KMPSearch(genORF,genoma,nombre3)

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
    "TAA": "STOP",
    "TGA": "STOP",
    "TAG": "STOP"
}

resultado = frames(genoma)
# for prot in resultado:
#     print("Inicio", prot["inicio"])
#     print("Final", prot["final"])
#     print("Secuencia", prot["secuencia"])
#     print("Frame", prot["frame"])

proteinas = leer_archivo_proteina("archivos/seq-proteins.txt")

for prot_res in resultado:
    for prot_ref in proteinas:
        if prot_res["secuencia"] == prot_ref["secuencia"]:
            print("Si hay coincidencia")
            print(prot_ref["nombre"])
            print(prot_res["frame"])
            print(prot_res["inicio"])
            print(prot_res["final"])
            print(prot_res["secuencia"])


# print(proteinas)
# print(genoma)






