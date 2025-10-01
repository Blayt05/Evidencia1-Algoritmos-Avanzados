genM   = leer_archivo("archivos/gen-M.txt") # el gen M  
genS   = leer_archivo("archivos\gen-S.txt") # el gen S
genORF = leer_archivo("archivos\gen-ORF1AB.txt") # el gen ORF1ab
nombre = "genM"
nombre2 = "genS"
nombre3 = "genORF"
KMPSearch(genM,genoma, nombre)
KMPSearch(genS,genoma,nombre2)
KMPSearch(genORF,genoma,nombre3)