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