Script library is tested for Matlab R2011+.
(See user-guide for functionality description)

To run "gbzmap_combo":

(1) put m/p-scripts in Matlab directory in your path
(e.g., user\Documents\Matlab\)

(2) From OS, create a separate dir to store the resulting mat/mhd output.
You may also put the "Demo...txt" (example coefficients) files in the 
same dir (for conevenience).

(3) Start Matlab, and CD to "results" directory:
>> cd your_results_dir;

(4) Run script  by calling:
>> BzXYZ = gbzmap_combo(np, sph_opt, mhd_out);

% INPUT params: 
%       np [int] -- default = 5, number of points within FOV (same 3D dimensions assumed)
%                   NOTE: odd number (preferred) will provide a point at isocenter.
%       sfopt [0/1/2] -- default = 0 (e.g., Philips convention), spherical harmonics coeff normalization 
%			use 1 for "yes" (e.g., GE convention), 
%			2 coil specific normalization (e.g., Siemens convention)
%       mhd_out ['y'/'n'] -- default = 'y', write out  MHD files for maps (use 'n' to suppress)
% 
% OUTPUT: 
% 	    BzXYZ.zx       double    3D iso volumetric Bzmap for GX gradient. 
% 	    BzXYZ.zy       double    3D iso volumetric Bzmap for GY gradient. 
% 	    BzXYZ.zz       double    3D iso volumetric Bzmap for GZ gradient. 

Output maps are also saved in the "Gfield...mat" file in your results dir.

(5) You can start by running with a small number of points (5/7/11) to compare
    with "demo" results. 

(6) For actual system map use np>=121 with characteristic system r0 ~ 300mm (FOV=600mm)

--------------------------------------
With questions, contact authors:
Dariya Malyarenko : dariya@umich.edu
Yuxi Pang: yuxipang@umich.edu