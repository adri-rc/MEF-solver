% =========================================================================
% SCRIPT PRINCIPAL: CASO 52 - FLEXIÓN PURA (CQUAD4)
% =========================================================================
clear all; clc;

% --- 1. CONFIGURACIÓN DE RUTAS ---
% Ruta donde quieres guardar los CSV
output_dir = 'C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_52';
% Prefijo para los archivos generados
filename_prefix = 'ej52_pure_flexure_ferreira';

% -- 2. Geometría (Nodos 1 a 4) --
% Ignoramos la coordenada Z=1.0 para la matriz local 2D, ya que Z es solo 
% un offset de cuerpo rígido en el sistema de coordenadas global.
nodeCoordinates = [0.1, 0.2; 
                   2.3, 0.1; 
                   2.1, 2.2; 
                   0.2, 1.9];
elementNodes = [1, 2, 3, 4];

% -- 3. Propiedades de los Materiales (Caso 52) --
thickness = 1.0; 
I = thickness^3 / 12.0;
kappa = 5.0 / 6.0; % Factor de corrección a cortante estándar

% Material 1 (MID2 - Flexión)
E_bend = 2.1e5;        
nu_bend = 0.30;   

% Material 3 (MID3 - Cortadura Transversal)
% Valor forzado a ~0 en Nastran para anular la rigidez a cortante
E_shear = 1.0e-20; 
nu_shear = 0.30;


% =========================================================================
% =========================================================================
% -- 4. Matrices Constitutivas --
% Matriz de rigidez a flexión (MID2)
C_bending = (I * E_bend / (1 - nu_bend^2)) * ...
            [1,       nu_bend, 0; 
             nu_bend, 1,       0; 
             0,       0,       (1-nu_bend)/2];

% Matriz de rigidez a cortante (MID3)
% Al ser E_shear ~ 0, esta matriz C_shear será ~0.
C_shear = (kappa * thickness * E_shear / (2 * (1 + nu_shear))) * eye(2);

% -- 5. Ejecución del wrapper --
cquad_bending_ferreira_core(nodeCoordinates, elementNodes, C_bending, C_shear, output_dir, filename_prefix);