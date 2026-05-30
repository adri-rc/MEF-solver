import numpy as np

E = 210000.0
nu = 0.3
t = 1.0

# 1. Definir las matrices
# Usamos np.array() y le pasamos una lista de listas
Bf = np.array([[0, 0, -0.1, 0, 0, 0.1, 0, 0, 0], 
                      [0, 0.1, 0, 0, 0, 0, 0, -0.1, 0],
                      [0, 0.1,-0.1, 0, -0.1, 0, 0, 0, 0.1]])

print(Bf)
print(20*'-')
print(Bf.T)

D = (E/(1-nu**2))*np.array([[1, nu, 0], 
                            [nu, 1, 0],
                            [0, 0,(1-nu)/2]])

K_flex = (t**3/12)*50*Bf.T@D@Bf
print(K_flex)

'''
matriz_B = np.array([[5, 6], 
                     [7, 8]])

# 2. Multiplicar las matrices
# En Python 3.5+, el operador '@' se usa para la multiplicación de matrices (producto punto)
resultado = matriz_A @ matriz_B

# (Alternativa clásica: resultado = np.dot(matriz_A, matriz_B))

# 3. Imprimir el resultado en la terminal
print("Matriz Resultado (NumPy):")
print(resultado)'''