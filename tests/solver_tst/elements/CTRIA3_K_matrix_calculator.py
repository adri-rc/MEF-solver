import numpy as np
import csv
import os

class Material:
    """Clase para almacenar las propiedades del material y la sección."""
    def __init__(self, E, nu):
        self.E = E
        self.nu = nu
        # El módulo de rigidez al corte se deduce de E y nu en materiales isótropos
        self.G = E / (2.0 * (1.0 + nu))
        
def matriz_constitutiva_membrana_tension_plana(E:float, nu:float):
    """Devuelve la matriz D (3x3) para TENSIÓN PLANA."""
    
    factor = E / (1.0 - nu**2)
    return factor * np.array([
        [1.0,  nu, 0.0],
        [nu, 1.0,  0.0],
        [0.0,  0.0,  (1.0 - nu) / 2.0]
    ])


class CTRIA3:
    """Clase que define un elemento triangular de 3 nudos."""
    def __init__(self, nodos, material, t, sistema='global'):
        """
        nodos: lista o array con las coordenadas de los 3 nudos.
               Si sistema='local', puede ser 2D [[x1,y1], [x2,y2], [x3,y3]] 
               o 3D asumiendo que Z=0.
               Si sistema='global', debe ser 3D [[X1,Y1,Z1], [X2,Y2,Z2], [X3,Y3,Z3]].
        material: instancia de la clase Material.
        sistema: 'local' o 'global'.
        """
        # Aseguramos que los nodos sean un array 3x3 (rellenando con Z=0 si vienen en 2D)
        nodos_arr = np.array(nodos, dtype=float)
        if nodos_arr.shape[1] == 2:
            self.nodos = np.hstack((nodos_arr, np.zeros((3,1))))
        else:
            self.nodos = nodos_arr
            
        self.mat = material
        self.t = t
        self.sistema = sistema
        
        # Inicializamos resultados
        self.T = np.eye(18) # Matriz de transformación (Identidad por defecto)
        self.K_local = np.zeros((18, 18))
        self.K_global = np.zeros((18, 18))
        # Atributos geométricos internos que calcularemos
        self.area = 0.0
        self.x_loc = np.zeros(3)  # Coordenadas X locales de los 3 nudos
        self.y_loc = np.zeros(3)  # Coordenadas Y locales de los 3 nudos
        
        # Indices de todos los DoFs.
        # [u1, v1, w1, theta_x1, theta_y1, theta_z1,
        #  u2, v2, w2, theta_x2, theta_y2, theta_z2,
        #  u3, v3, w3, theta_x3, theta_y3, theta_z3]
        # Indices de membrana en la matriz de 18x18
        self.idx_m = [0, 1, 6, 7, 12, 13]
        
    # ==========================================================
    # 1. GEOMETRÍA Y TRANSFORMACIÓN
    # ==========================================================
    def _calcular_geometria_y_matriz_T(self):
        """Calcula coordenadas locales, área y matriz de transformación T."""
        N1, N2, N3 = self.nodos[0], self.nodos[1], self.nodos[2]
        
        if self.sistema == 'local':
            x1, y1 = N1[0], N1[1]
            x2, y2 = N2[0], N2[1]
            x3, y3 = N3[0], N3[1]
            self.T = np.eye(18)
        else:
            v12 = N2 - N1
            v13 = N3 - N1
            
            L12 = np.linalg.norm(v12)
            e_x = v12 / L12
            v_z = np.cross(e_x, v13)
            e_z = v_z / np.linalg.norm(v_z)
            e_y = np.cross(e_z, e_x)
            
            x1, y1 = 0.0, 0.0
            x2, y2 = L12, 0.0
            x3 = np.dot(v13, e_x)
            y3 = np.dot(v13, e_y)
            
            # Construcción de la matriz T (18x18)
            Lambda = np.vstack((e_x, e_y, e_z))
            self.T = np.zeros((18, 18))
            for i in range(0, 18, 3):
                self.T[i:i+3, i:i+3] = Lambda

        # Guardamos en atributos de clase para uso posterior
        self.x_loc = np.array([x1, x2, x3])
        self.y_loc = np.array([y1, y2, y3])
        self.area = 0.5 * abs(x2*y3 - x3*y2)

    # ==========================================================
    # 2. FORMULACIÓN DE MEMBRANA (CST)
    # ==========================================================
    def _calcular_matriz_B_membrana(self):
        """Calcula la matriz B (3x6) de deformación-desplazamiento para membrana."""
        x1, x2, x3 = self.x_loc
        y1, y2, y3 = self.y_loc
        
        b1, c1 = y2 - y3, x3 - x2
        b2, c2 = y3 - y1, x1 - x3
        b3, c3 = y1 - y2, x2 - x1
        
        B = (1.0 / (2.0 * self.area)) * np.array([
            [b1,  0, b2,  0, b3,  0],
            [ 0, c1,  0, c2,  0, c3],
            [c1, b1, c2, b2, c3, b3]
        ])
        return B

    def _calcular_k_membrana_local(self, B):
        """Calcula la matriz de rigidez de membrana k_m (6x6)."""
        D = matriz_constitutiva_membrana_tension_plana(self.mat.E, self.mat.nu)
        k_m = self.t * self.area * np.dot(B.T, np.dot(D, B))
        return k_m

    # ==========================================================
    # 3. ENSAMBLAJE Y TRANSFORMACIÓN GLOBAL
    # ==========================================================
    def _ensamblar_K_local(self, k_m):
        """Inserta las submatrices (membrana, flexión) en la K_local de 18x18."""
        self.K_local = np.zeros((18, 18))
        
        # Ensamblamos la parte de membrana en sus índices correspondientes
        for i in range(6):
            for j in range(6):
                self.K_local[self.idx_m[i], self.idx_m[j]] = k_m[i, j]
                
        # NOTA FUTURA: Aquí añadirás la lógica para ensamblar k_flexion (9x9) 
        # en los DoFs correspondientes a [w, th_x, th_y]

    def _transformar_a_global(self):
        """Rota la matriz K_local al sistema de coordenadas K_global."""
        self.K_global = np.dot(self.T.T, np.dot(self.K_local, self.T))

    # ==========================================================
    # 4. ORQUESTADOR PRINCIPAL
    # ==========================================================
    def calcular_matrices(self):
        """Ejecuta en orden todos los pasos para obtener K_local y K_global."""
        
        # Paso 1: Obtener ejes locales, área y matriz de rotación
        self._calcular_geometria_y_matriz_T()
        
        # Paso 2: Calcular matriz B de membrana
        B_membrana = self._calcular_matriz_B_membrana()
        
        # Paso 3: Calcular rigidez local de membrana 6x6
        k_m = self._calcular_k_membrana_local(B_membrana)
        
        # Paso 4: Ensamblar en matriz del elemento 18x18
        self._ensamblar_K_local(k_m)
        
        # Paso 5: Rotar al sistema global
        self._transformar_a_global()

# --- FUNCIONES DE UTILIDAD ---
def print_matriz(K, titulo):
    """Imprime una matriz 18x18 filtrando los ceros para mejor lectura."""
    dof_names =  ['u1', 'v1', 'w1', 'theta_x1', 'theta_y1', 'theta_z1',
                  'u2', 'v2', 'w2', 'theta_x2', 'theta_y2', 'theta_z2',
                  'u1', 'v1', 'w1', 'theta_x3', 'theta_y3', 'theta_z3']
    print(f"\n{'='*60}\n{titulo}\n{'='*60}")
    print("      " + "".join([f"{name:>9}" for name in dof_names]))
    
    for i in range(18):
        fila = f"{dof_names[i]:<4}|"
        for j in range(18):
            val = K[i, j]
            if abs(val) < 1e-6:
                fila += "      .  "
            else:
                fila += f"{val:9.1e}"
        print(fila)

# ==========================================
# CASO 1: TRABAJANDO EN COORDENADAS LOCALES
# ==========================================
'''
     --- CASO DE VALIDACIÓN 34 ---
PSHELL  1       1       0.001
MAT1    1       70.E9           0.3
$
GRID    1               0.0     0.0     0.0
GRID    2               1.0     0.0     0.0
GRID    3               0.5     0.866   0.0
$
CTRIA3  1       1       1       2       3
'''
print("\n>>> CASO 1: Elemento definido en coordenadas LOCALES (2D)")
mat = Material(E=70e9, nu=0.3)

# Nodos en 2D (x, y)
nodos_locales = [[100.0, 100.0, 0.0], [150.0, 120.0, 0.0], [80.0, 200.0, 0.0]]
elem_local = CTRIA3(nodos_locales, mat, t=0.001, sistema='local')
elem_local.calcular_matrices()

print_matriz(elem_local.K_local, "Matriz K_local (Coordenadas Locales)")


# ==========================================
# CASO 2: TRABAJANDO EN COORDENADAS GLOBALES
# ==========================================
print("\n>>> CASO DE VALIDACIÓN 34")
'''
     --- CASO DE VALIDACIÓN 34 ---
PSHELL  1       1       0.001
MAT1    1       70.E9           0.3
$
GRID    1               0.0     0.0     0.0
GRID    2               1.0     0.0     0.0
GRID    3               0.5     0.866   0.0
$
CTRIA3  1       1       1       2       3
'''
print("\n>>> CASO 1: Elemento definido en coordenadas LOCALES (2D)")
# Cogemos el mismo triángulo, pero lo giramos 90 grados alrededor del eje X
# Ahora el triángulo está en el plano XZ en lugar del XY
nodos_globales = [[100.0, 100.0, 0.0], [150.0, 120.0, 0.0], [80.0, 200.0, 0.0]]
elem_global = CTRIA3(nodos_globales, mat, t=0.001,sistema='global')
elem_global.calcular_matrices()

print_matriz(elem_global.K_local, "Matriz K_local (Interna del elemento inclinado)")
print_matriz(elem_global.T, "Matriz de Transformación (T)")
print_matriz(elem_global.K_global, "Matriz K_global (Vista desde el sistema general)")


# --- FUNCIÓN DE EXPORTACIÓN A CSV ---
def exportar_matriz_csv(K, nombre_archivo):
    """
    Exporta una matriz de 18x18 a un archivo CSV.
    Guarda el archivo en el mismo directorio donde se encuentra este script, o en la ruta indicada.
    Incluye una primera columna con el índice de Python (0-17).
    """
    dof_names =  ['u1', 'v1', 'w1', 'theta_x1', 'theta_y1', 'theta_z1',
                  'u2', 'v2', 'w2', 'theta_x2', 'theta_y2', 'theta_z2',
                  'u3', 'v3', 'w3', 'theta_x3', 'theta_y3', 'theta_z3']
    
    if os.path.isabs(nombre_archivo):
        # Si detecta que es una ruta completa (absoluta), la usa tal cual
        ruta_completa = nombre_archivo
    else:
        # Si es solo el nombre (o ruta relativa), la une al directorio del script
        directorio_script = os.path.dirname(os.path.abspath(__file__))
        ruta_completa = os.path.join(directorio_script, nombre_archivo)
    # -----------------------------------
    
    with open(ruta_completa, mode='w', newline='') as file:
        # Usamos coma ',' como separador estándar CSV. 
        writer = csv.writer(file, delimiter=',')
        
        # 1. Escribir la cabecera
        writer.writerow(['Index', 'DoF'] + dof_names)
        
        # 2. Escribir cada fila
        for i in range(18):
            fila_valores = [f"{val:.3e}" if abs(val) > 1e-10 else "0.0" for val in K[i]]
            # Escribimos: [Índice i] + [Etiqueta Fila] + [Valores de la matriz...]
            writer.writerow([i, dof_names[i]] + fila_valores)
            
    print(f"✅ Matriz exportada exitosamente a: \n   --> '{ruta_completa}'")

# 4. Exportamos las matrices a CSV
exportar_matriz_csv(elem_local.K_local, r"C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_36\draft_ej36_calculadora_K_Local.csv")
exportar_matriz_csv(elem_global.K_global, r"C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_36\draft_ej36_calculadora_K_Global.csv")
exportar_matriz_csv(elem_global.T, r"C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_36\draft_ej36_calculadora_Transformacion.csv")
