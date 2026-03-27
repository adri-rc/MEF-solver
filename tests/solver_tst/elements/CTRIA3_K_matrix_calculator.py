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
    """
    --- CHEQUEADA CON REFERENCIAS TFM (ref ec MEF-español (4.18)) ---
    Devuelve la matriz D (3x3) para TENSIÓN PLANA."""
    
    factor = E / (1.0 - nu**2)
    return factor * np.array([
        [1.0,  nu, 0.0],
        [nu, 1.0,  0.0],
        [0.0,  0.0,  (1.0 - nu) / 2.0]
    ])

def matriz_constitutiva_flexion_tension_plana(E:float, nu:float, t:float):
    """
    NO encuentro referencias a esta fórmula, según ref MEF (8.7) la expresión de D es la misma que en membrana.
    Devuelve la matriz D_f (3x3) para FLEXIÓN de placas.
    Relaciona momentos con curvaturas.
    """
    factor = (E * t**3) / (12.0 * (1.0 - nu**2))
    return factor * np.array([
        [1.0,  nu, 0.0],
        [nu, 1.0,  0.0],
        [0.0,  0.0,  (1.0 - nu) / 2.0]
    ])

def matriz_constitutiva_cortante(G:float, t:float, k_shear=5.0/6.0):
    """
    Devuelve la matriz D_c (2x2) para CORTADURA TRANSVERSAL.
    Incluye el factor de corrección de cortante de Mindlin (k_shear).
    """
    factor = k_shear * G * t
    return factor * np.array([
        [1.0, 0.0],
        [0.0, 1.0]
    ])


class CTRIA3:
    """Clase que define un elemento triangular de 3 nudos."""
    def __init__(self, nodos, material_membrana, t, material_flexion=None, material_cortante=None, sistema='global'):
        """
        nodos: lista o array con las coordenadas de los 3 nudos.
               Si sistema='local', puede ser 2D [[x1,y1], [x2,y2], [x3,y3]] 
               o 3D asumiendo que Z=0.
               Si sistema='global', debe ser 3D [[X1,Y1,Z1], [X2,Y2,Z2], [X3,Y3,Z3]].
        material: instancia de la clase Material para las propiedades de material de membrana  (MID1).
        t: espesor del elemento.
        material_flexion: Material para la flexión (MID2). Si es None, usa el de membrana.
        material_cortante: Material para el cortante transversal (MID3). Si es None, usa el de flexión.
        sistema: 'local' o 'global'.
        """
        # Aseguramos que los nodos sean un array 3x3 (rellenando con Z=0 si vienen en 2D)
        nodos_arr = np.array(nodos, dtype=float)
        if nodos_arr.shape[1] == 2:
            self.nodos = np.hstack((nodos_arr, np.zeros((3,1))))
        else:
            self.nodos = nodos_arr
            
        self.mat_m = material_membrana
        self.mat_f = material_flexion if material_flexion is not None else material_membrana
        self.mat_c = material_cortante if material_cortante is not None else self.mat_f
        
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
        # Indices de membrana en la matriz de 18x18 [u1, v1, u2, v2, u3, v3]
        self.idx_m = [0, 1, 6, 7, 12, 13]

        # Indices de flexión en la matriz de 18x18 [w1, thx1, thy1, w2, thx2, thy2, w3, thx3, thy3]
        self.idx_f = [2, 3, 4, 8, 9, 10, 14, 15, 16]
        
    # ==========================================================
    # 1. GEOMETRÍA Y TRANSFORMACIÓN
    # ==========================================================
    def _calcular_geometria_y_matriz_T(self):
        """
        Calcula coordenadas locales, área y matriz de transformación T.
        """
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
        # La fórmula del siguinete área viene de la ecuación (ref ec MEF-español (4.35)) libro MEF español.
        # A = 0.5*Det([1, x1, y1],[1, x2, y2],[1, x3, y3]) que si desarrollamos:
        # 2*A = 1*(x_2 y_3 - x_3 y_2) - 1*(x_1 y_3 - x_3 y_1) + 1*(x_1 y_2 - x_2 y_1)
        # como x1 = y1 = 0
        # 2A = (x_2 * y_3 -x_3 y_2) - 0 - 0 + 0
        # A = 0.5*(x_2 * y_3 -x_3 y_2)
        # Dado que y2 =0 -> A = 0.5 * (x_2 * y_3) ((base*altura)/2)
        # Pero se deja en la ecuación por si un usuario definen unas coordenadas locales donde y2 no sea 0.
        self.area = 0.5 * abs(x2*y3 - x3*y2)

    # ==========================================================
    # 2A. FORMULACIÓN DE MEMBRANA (CST)
    # ==========================================================
    def _calcular_matriz_B_membrana(self):
        """
        --- CHEQUEADA CON REFERENCIAS TFM (ref ec MEF-español (4.37)) ---
        Calcula la matriz B (3x6) de deformación-desplazamiento para membrana.
        """
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
        """
        --- CHEQUEADA CON REFERENCIAS TFM (ref ec MEF-español (4.41)) ---
        Calcula la matriz de rigidez de membrana k_m (6x6).
        """
        D = matriz_constitutiva_membrana_tension_plana(self.mat_m.E, self.mat_m.nu)
        k_m = self.t * self.area * np.dot(B.T, np.dot(D, B))
        return k_m

    # ==========================================================
    # 2B. FORMULACIÓN DE FLEXIÓN (MINDLIN-REISSNER)
    # ==========================================================
    def _calcular_matrices_B_flexion_cortante(self):
        """
        Calcula la matriz B_f (3x9) de flexión y B_c (2x9) de cortante.
        """
        x1, x2, x3 = self.x_loc
        y1, y2, y3 = self.y_loc
        
        b1, c1 = y2 - y3, x3 - x2
        b2, c2 = y3 - y1, x1 - x3
        b3, c3 = y1 - y2, x2 - x1
        
        # 1. Matriz de Curvaturas B_f (Evaluada para [w, th_x, th_y])
        B_f = (1.0 / (2.0 * self.area)) * np.array([
            [0,   0,  b1,  0,   0,  b2,  0,   0,  b3],
            [0, -c1,   0,  0, -c2,   0,  0, -c3,   0],
            [0, -b1,  c1,  0, -b2,  c2,  0, -b3,  c3]
        ])
        
        # 2. Matriz de Cortante Transversal B_c
        # Evaluada en el centroide (Integración reducida de 1 punto para evitar Shear Locking)
        # En el centroide, las funciones de forma valen N1 = N2 = N3 = 1/3
        Nc = 1.0 / 3.0
        A2 = 2.0 * self.area
        
        B_c = np.array([
            [b1/A2,    0,  Nc,  b2/A2,    0,  Nc,  b3/A2,    0,  Nc],
            [c1/A2,  -Nc,   0,  c2/A2,  -Nc,   0,  c3/A2,  -Nc,   0]
        ])
        
        return B_f, B_c

    def _calcular_k_flexion_local(self, B_f, B_c):
        """
        Calcula la matriz de rigidez de flexión+cortante k_f (9x9).
        """
        D_f = matriz_constitutiva_flexion_tension_plana(self.mat_f.E, self.mat_f.nu, self.t)
        D_c = matriz_constitutiva_cortante(self.mat_c.G, self.t)
        
        # Integral analítica en el área: A * B^T * D * B
        k_flex = self.area * np.dot(B_f.T, np.dot(D_f, B_f))
        k_cort = self.area * np.dot(B_c.T, np.dot(D_c, B_c))
        
        return k_flex + k_cort
    
     
    # ==========================================================
    # 3. ENSAMBLAJE Y TRANSFORMACIÓN GLOBAL
    # ==========================================================
    def _ensamblar_K_local(self, k_m, k_f):
        """Inserta las submatrices (membrana, flexión) en la K_local de 18x18."""
        self.K_local = np.zeros((18, 18))
        
        # Ensamblamos la parte de membrana en sus índices correspondientes
        for i in range(6):
            for j in range(6):
                self.K_local[self.idx_m[i], self.idx_m[j]] = k_m[i, j]
                
        # NOTA FUTURA: Aquí añadirás la lógica para ensamblar k_flexion (9x9) 
        # Ensamblamos la parte de flexión (curvatura + cortante)
        for i in range(9):
            for j in range(9):
                self.K_local[self.idx_f[i], self.idx_f[j]] = k_f[i, j]

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
        
        # Paso 2: Calcular matrices B (Membrana, Flexión y Cortante)
        B_membrana = self._calcular_matriz_B_membrana()
        B_flexion, B_cortante = self._calcular_matrices_B_flexion_cortante()
        
        # Paso 3: Calcular rigidez local de membrana 6x6 & Flexión + Cortante 9x9
        k_m = self._calcular_k_membrana_local(B_membrana)
        k_f = self._calcular_k_flexion_local(B_flexion, B_cortante)
        
        # Paso 4: Ensamblar en matriz del elemento 18x18
        self._ensamblar_K_local(k_m, k_f)
        
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
print("\n>>> CASO 1: Elemento definido en coordenadas LOCALES (2D)")
mat_m = Material(E=70e9, nu=0.3)
mat_f = Material(E=70e9, nu=0.3)
mat_c = Material(E=70e9, nu=0.3)

# Nodos en 2D (x, y)
nodos_locales = [[100.0, 100.0, 0.0], [150.0, 120.0, 0.0], [80.0, 200.0, 0.0]]
elem_local = CTRIA3(nodos_locales, material_membrana=mat_m, t=50.0, material_flexion=mat_f, material_cortante=mat_c, sistema='local')
elem_local.calcular_matrices()

print_matriz(elem_local.K_local, "Matriz K_local (Coordenadas Locales)")
'''

# ==========================================
# CASO 2: TRABAJANDO EN COORDENADAS GLOBALES
# ==========================================
print("\n>>> CASO DE VALIDACIÓN XY")

print("\n>>> CASO 1: Elemento definido en coordenadas LOCALES (2D)")
# Cogemos el mismo triángulo, pero lo giramos 90 grados alrededor del eje X
# Ahora el triángulo está en el plano XZ en lugar del XY
nodos_globales = [[100.0, 100.0, 0.0], [150.0, 120.0, 0.0], [80.0, 200.0, 0.0]]
elem_global = CTRIA3(nodos_globales, material_membrana=mat_m, t=50.0, material_flexion=mat_f, material_cortante=mat_c, sistema='global')
elem_global.calcular_matrices()

print_matriz(elem_global.K_local, "Matriz K_local (Interna del elemento inclinado)")
print_matriz(elem_global.T, "Matriz de Transformación (T)")
print_matriz(elem_global.K_global, "Matriz K_global (Vista desde el sistema general)")



# ==========================================
#  --- FUNCIÓN DE EXPORTACIÓN A CSV ---
# ==========================================
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
#exportar_matriz_csv(elem_local.K_local, r"draft_benging_calculadora_K_Local_bueno2.csv")
exportar_matriz_csv(elem_global.K_global, r"draft_benging_calculadora_K_Global_bueno3.csv")
#exportar_matriz_csv(elem_global.T, r"draft_benging_calculadora_Transformacion.csv")
