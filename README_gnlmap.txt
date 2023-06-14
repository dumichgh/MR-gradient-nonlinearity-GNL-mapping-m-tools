p-code library is tested for Matlab R2011-R2015.

Functionality: 
	(a) you can use "gnlmap" script to generate 3D maps for
	    nonlinearity tensor (normalized to isocenter values) from the pre-generated
	    gradient field map structure: BzXYZ (see "gbzmap_combo" and manual);
	(b) you also have an option to generate nonlinearity maps as MHD files
	    to visualize, e.g., with the "3DSlicer", when "mhdpars" structure
	    is provided (e.g., as saved for original BzXYZ field maps)
	(c) when input gradient field map is zero,e.g., for BzXYZ.zx (GRL-coil),
	    the script will attempt generating missing field (from one 
	    that was provided) by 90-degree flip (GRL-->GAP or GAP-->GRL)
	    assuming "cylindrical" symmetry (horizontal bore gradients)	


To run "gnlmap":

(1) put p-script in Matlab directory in your path
(e.g., user\Documents\Matlab\)

(2) Start Matlab, and CD to "results" directory:
>> cd your_results_dir;

(3) load previousely saved 'BzXYZ'structure (if not yet in your workspace):
>> load('Gfield_BzXYZ_R125x125_1512041503.mat', 'BzXYZ')

(4) optional -- if desire MHD generation, load corresponding 'mhdpar' for BzXYZ:
>> load('PARS2buildMHD_4Dn_Gfield_R125x125_1512041503_SIMS.mat', 'mhdpars')

(4) Run script  by calling, e.g.:
>> Lxyz = gnlmap(BzXYZ, mhdpars); % NOTE: omitting second argument will suppress MHD generation

(5) The resulting nine GNL tensor-element maps will be saved as "GNLmap_Lxyz_...mat" 
by default: dim1=RL; dim2=AP; dim3=IS (and as MHD in LPS if 'mhdpars' provided)
NOTE: "NaN" in mat-file map (where Lxyz is not defined) is replaced by "zero" in 
the corresponding MHD

% INPUT:
% BzXYZ (required) -- three 3D Bz maps in magnet coordinates 
% for RL(x), AP(y) and SI(z) gradients
% mhdpar (opt) -- MHD param struture for BzXYZ (see, "gbzmap_combo"
% for requirements)
% 
%
% OUTPUT 
% Lxyz - 9 - 3D matrix structure for gradient correction
% corresponding map-by-map MHDs are built if "mhdpar" is provided
%
% NOTE: depends on "build_mhd_combo" library to build MHDs

--------------------------------------
With questions, contact authors:
Dariya Malyarenko : dariya@umich.edu
