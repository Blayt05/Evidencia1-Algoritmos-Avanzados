
#Funcion que nos ayuda a normalizar nuestro genomas y genes en una sola linea sin espacios
def leer_archivo(ruta):
    seq = []
    with open(ruta, "r", encoding="utf-8") as f: #Se abre el archivo fasta y se cierra con with
        for line in f: #Para cada linea del archivo 
            if line.startswith(">"): # Si empieza con > 
                continue  # Se ignora la cabecera
            seq.append(line.strip())  # Se quitan saltos de línea y espacios
    return "".join(seq) # Se unen todos los arreglos en uno solo sin espacios

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
            
                

        

# -------- Archivos de prueba -----------
#Este seria el texto 
genoma = leer_archivo("archivos/SARS-COV-2-MN908947.3.txt") # Archivo del genoma completo

#Estos serian los patrones
genM   = leer_archivo("archivos/gen-M.txt") # el gen M  
genS   = leer_archivo("archivos\gen-S.txt") # el gen S
genORF = leer_archivo("archivos\gen-ORF1AB.txt") # el gen ORF1ab
nombre = "genM"
nombre2 = "genS"
nombre3 = "genORF"
KMPSearch(genM,genoma, nombre)
KMPSearch(genS,genoma,nombre2)
KMPSearch(genORF,genoma,nombre3)

# print(genoma)




