% =========================================================================
% SCRIPT MAIN: Caso 51 - CQUAD con geometría arbitraria
% =========================================================================
clear all; clc; format short e;

disp('===================================================');
disp('   EJECUTANDO CASO 51: CQUAD GEOMETRÍA ASIMÉTRICA  ');
disp('===================================================');

% --- 1. CONFIGURACIÓN DE RUTAS ---
% Ruta donde quieres guardar los CSV
output_dir = 'C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_51';
% Prefijo para los archivos generados
filename_prefix = 'ej51_pure_membrane_ferreira';

% --- 2. MATERIALES (MAT1) ---
E = 2.1e5;        
poisson = 0.30;   
C = E/(1-poisson^2) * [1 poisson 0; poisson 1 0; 0 0 (1-poisson)/2];

% --- 3. GEOMETRÍA (GRID) ---
% Extraído del archivo .dat (X, Y) ignorando Z
nodeCoordinates = [
    0.1, 0.2; 
    2.3, 0.1; 
    2.1, 2.2; 
    0.2, 1.9  
    ];

% --- 4. CONECTIVIDAD Y PROPIEDADES (CQUAD4 y PSHELL) ---
elementNodes = [1, 2, 3, 4];
thickness = 1.0; 

% --- 5. EJECUCIÓN Y EXPORTACIÓN ---
cquad_membrane_ferreira_core(nodeCoordinates, elementNodes, C, thickness, output_dir, filename_prefix);