% =========================================================================
% CASE 51: CQUAD WITH ARBITRARY GEOMETRY
% =========================================================================
clear all; clc; format short e;

disp('===================================================');
disp('   RUNNING CASE 51: ASYMMETRIC CQUAD GEOMETRY      ');
disp('===================================================');

% --- 1. OUTPUT CONFIGURATION ---
output_dir = fileparts(mfilename('fullpath'));
% Prefix used by generated CSV files.
filename_prefix = 'ej51_pure_membrane_ferreira';

% --- 2. MATERIAL (MAT1) ---
E = 2.1e5;        
poisson = 0.30;   
C = E/(1-poisson^2) * [1 poisson 0; poisson 1 0; 0 0 (1-poisson)/2];

% --- 3. GEOMETRY (GRID) ---
% Extracted from the input deck; Z is omitted from this 2D reference.
nodeCoordinates = [
    0.1, 0.2; 
    2.3, 0.1; 
    2.1, 2.2; 
    0.2, 1.9  
    ];

% --- 4. CONNECTIVITY AND PROPERTIES (CQUAD4 AND PSHELL) ---
elementNodes = [1, 2, 3, 4];
thickness = 1.0; 

% --- 5. EXECUTION AND EXPORT ---
cquad_membrane_ferreira_core(nodeCoordinates, elementNodes, C, thickness, output_dir, filename_prefix);
