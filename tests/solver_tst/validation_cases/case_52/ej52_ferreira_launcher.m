% =========================================================================
% CASE 52: PURE CQUAD4 BENDING
% =========================================================================
clear all; clc;

% --- 1. OUTPUT CONFIGURATION ---
output_dir = fileparts(mfilename('fullpath'));
% Prefix used by generated CSV files.
filename_prefix = 'ej52_pure_flexure_ferreira';

% --- 2. GEOMETRY (NODES 1 TO 4) ---
% Z=1.0 is omitted from the local 2D matrix because it is only a rigid-body
% offset in the global coordinate system.
nodeCoordinates = [0.1, 0.2; 
                   2.3, 0.1; 
                   2.1, 2.2; 
                   0.2, 1.9];
elementNodes = [1, 2, 3, 4];

% --- 3. MATERIAL PROPERTIES ---
thickness = 1.0; 
I = thickness^3 / 12.0;
kappa = 5.0 / 6.0; % Standard transverse-shear correction factor.

% Material 1 (MID2 bending).
E_bend = 2.1e5;        
nu_bend = 0.30;   

% Material 3 (MID3 transverse shear).
% Nastran uses a near-zero value to suppress physical shear stiffness.
E_shear = 1.0e-20; 
nu_shear = 0.30;


% =========================================================================
% =========================================================================
% --- 4. CONSTITUTIVE MATRICES ---
% Bending stiffness matrix (MID2).
C_bending = (I * E_bend / (1 - nu_bend^2)) * ...
            [1,       nu_bend, 0; 
             nu_bend, 1,       0; 
             0,       0,       (1-nu_bend)/2];

% Transverse-shear stiffness matrix (MID3). With near-zero E_shear,
% C_shear is also near zero.
C_shear = (kappa * thickness * E_shear / (2 * (1 + nu_shear))) * eye(2);

% --- 5. EXECUTION ---
cquad_bending_ferreira_core(nodeCoordinates, elementNodes, C_bending, C_shear, output_dir, filename_prefix);
