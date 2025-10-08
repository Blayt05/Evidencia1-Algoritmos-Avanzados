
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
        print(f"- Longitud de la secuencia: {len(seq)}")
        print(f"- Longitud del palindromo mas largo: {len(longest)}")
        print(f"- Palindromo mas largo: {longest}")
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

#********************************************************************
#---------------------- LLamadas a problema 4 ----------------------
#********************************************************************
print("********************************************************************")
print("---------------------- Comienzo a problema 4 ----------------------")
print("********************************************************************")




