import math
import struct

def leggi_numeri_double(file_path):
    result = []
    try:
        # Apre il file binario in modalità lettura
        with open(file_path, "rb") as file:
            i = 0
            while True:
                # Legge 8 byte alla volta (un numero double)
                bytes_double = file.read(8)
                
                # Se non ci sono più dati nel file, esci dal ciclo
                if not bytes_double:
                    break
                
                # Converte i 8 byte letti in un numero double
                numero_double = struct.unpack('d', bytes_double)[0]
                i = i+1
                # Stampa il numero
                result.append(numero_double)
                print(f"Double: {i} {numero_double}")
            return result
    except Exception as e:
        print(f"Errore: {e}")
        return result

def leggi_numeri_float(file_path):
    result = []
    try:
        # Apre il file binario in modalità lettura
        with open(file_path, "rb") as file:
            i = 0
            while True:
                # Legge 4 byte alla volta (un numero float)
                bytes_float = file.read(4)
                
                # Se non ci sono più dati nel file, esci dal ciclo
                if not bytes_float:
                   break
                
                # Converte i 4 byte letti in un numero float
                numero_float = struct.unpack('f', bytes_float)[0]
                i = i+1
                # Stampa il numero
                result.append(numero_float)
                print(f"Float {i}: {numero_float}") 
            return result
    except Exception as e:
        print(f"Errore durante la lettura dei numeri float: {e}")
        return result

def rmse(vettore_a, vettore_b):
    # Verifica che i vettori abbiano la stessa lunghezza
    if len(vettore_a) != len(vettore_b):
        raise ValueError("I vettori devono avere la stessa lunghezza.")

    # Calcolo della somma dei quadrati delle differenze
    somma_quadrati = 0
    for i in range(len(vettore_a)):
        somma_quadrati += (vettore_a[i] - vettore_b[i]) ** 2

    # Calcolo del RMSE
    errore_quadratico_medio = math.sqrt(somma_quadrati / len(vettore_a))
    return errore_quadratico_medio

def sono_simili(vettore_a, vettore_b, soglia=1.3):
    # Calcola il RMSE
    errore = rmse(vettore_a, vettore_b)

    # Confronta il RMSE con la soglia
    return errore < soglia

# Esempio d'uso

file_path = "phi_256_to20_k1_alpha1_sd3.ds2" 
file_path_result = "out32_256_3_phi.ds2"
vettore_a = leggi_numeri_double(file_path)
vettore_a.pop(0)
vettore_b = leggi_numeri_float(file_path_result)
vettore_b.pop(0)
vettore_b.pop(1)

print(len(vettore_a))
print(len(vettore_b))

# Verifica se i vettori sono simili
if sono_simili(vettore_a, vettore_b):
    print("I vettori sono simili.")
else:
    print("I vettori non sono simili.")
