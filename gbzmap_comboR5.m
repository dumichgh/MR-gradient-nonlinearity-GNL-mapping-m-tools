
function BzXYZ = gbzmap_combo7(np, sfopt, mhd_out)
% ------------------------------------------------------------------------------
% PURPOSE:
% 	This script is used to construct Gradient coil Bz field map from input 
%   spherical harmonic (SPH) coeffients and Bzmap geom info (r0 and npts).
% NOTE1: Field maps will be calculated for full FOV within r0
% (characteristic radius) as loaded from the sys-coefficient file for AP, RL and SI gradients
%
% INPUT params: 
%       np -- default = 5, number of points within FOV (same 3D dimensions assumed)
%       sfopt -- default = 0, sph-harm coeff normalization (use 1 for "yes", 2 for SMS convention(norm(n,0)=1))
%       mhd_out -- default = 'y', write out  MHD files for maps (use 'n' for "no")
% 
% OUTPUT: 
% 	    BzXYZ.zx       double    3D iso volumetric Bzmap for GX (RL) gradient. 
% 	    BzXYZ.zy       double    3D iso volumetric Bzmap for GY (AP) gradient. 
% 	    BzXYZ.zz       double    3D iso volumetric Bzmap for GZ (SI) gradient. 
%
% NOTES:
%   This script will call subfunction "load_sph_coeffs.m" and 
%   "scale_sph_coeffs.m" (see below). The results will be saved (if 
%   selected) as MHD/RAW format for easy visulazation with 3D Slicer.
%   global field norm  hard-coded to default = 0 (code 1 for "yes")
%
% MODIFICATION HISTORY:
% 	17/09/2015: 	v1.0, initial draft
% 	23/09/2015: 	v1.1, the charateristic radius (fov/2.0) moved into SPH coeffs text file
%   24/09/2015:     v1.2, reset SPH undefined region (i.e. R > Rsph) to
%                         ZERO based on Dasha's suggestion.
%   30/09/2015:     v1.3, consolidate three outputs into one structute to
%                         align with the input for numLxyzR.m
%   02/10/2015:     v1.4, vectorize the core calculation by replacing four loops
%   05/10/2015:     v1.5, fixed a bug in v6 which wrongly scaled GzBz with sclx
%                         (should be sclz!). 
% MD 20151117: same time-stamp from MHD and mat-output
% MD 20151201: hard-code norm = 0 (do not perfrom "global" field normalization)
%               and only use 3 input args
%   06/19/2017:     v1.6, added functionality to handle different
%                         characteristics radii for gradient magnetic field map calculation. 
% MD 20180212: for R1.3(v1.7) rotate final Bz maps into "DICOM frame" X=RL, Y=AP
%		Also, change reading input demos to explicitely label gradients GRL (former GY), 
%		GAP(former GX), GSI(former GZ) 
% MD 20180215:  chnage mhd-output pars to reflect new X=RL, Y=AP  (from input X=AP, Y=RL) 
% MD 20180223:   before making MHDs, convert output from MAT(magnet)-frame to LPS (to align RL/AP)
% YP 20180821:  update max SPH order (N_MAX) from default static value (16) to dynamic,
%               see details in subfunction "load_sph_coeffs.m".
% MD 20180907: check/allow for fipped XY convention for Siemens (VIDA): dim2=horizontal for RL channel
%              do not permute dim1-dim2 (to make dim1=horiz) before saving
% MD 20180910: generalize the XY check for single channel coeff (GRL or GAP) when other is missing (e.g. "DEMO_ONE")
%
% MD 20190828: added scale-option=2 for Siemens convention norm(n,0)=1; and normSPH(n,m)=sqrt(2*pi)*normSPH(n,m)
% ------------------------------------------------------------------------------

norm = 0; % MD 20151201 : do not normalize the fields
if (nargin < 3)
  disp('(FYI) Input options not specified -- using defaults for NON-SUPPLIED arguments: ');
  disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
%  disp('Usage: BzXYZ = make_gbzmap_v7(np[int], sph_opt[int], field_opt[int]), mhd_out[char])');
  disp('Usage: BzXYZ = gbzmap_combo(np[int], sph_opt[int], mhd_out[char])');
  disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
  disp('    default np = 5;  number of points in dim1 (dim2 or dim3)');
  disp('    default sph_opt = 0;  use yes(1) or no(0) coeffs scaled by SPH NORM factor');
%  disp('    default field_opt =0; use yes(1) or no(0) BzGx field scaled by r0/C11, BzGy by r0/S11 and BzGz by r0/C10');
  disp('    default mhd_out = "y"; yse "n" to suppress MHD output              ');
  disp(' ');
  disp('   OUTPUT:  BzXYZ.zx - 3D volumetric (np*np*np) Bzmap induced by Gradient X coil');
  disp('            BzXYZ.zy - 3D volumetric (np*np*np) Bzmap induced by Gradient Y coil'); 
  disp('            BzXYZ.zz - 3D volumetric (np*np*np) Bzmap induced by Gradient Z coil'); 
%   Bx = -1;
%   By = -1;
%   Bz = -1;
%   return
if (~exist('np', 'var')) np = 5; disp(' (FYI) using default number of points'); end
if (~exist('sfopt', 'var')) sfopt = 0; disp(' (FYI) using default coeff normalization'); end
%if (~exist('norm', 'var')) norm = 0; disp(' (FYI) using default field normalization'); end % MD 20151201 
if (~exist('mhd_out', 'var')) mhd_out = 'y'; disp(' (FYI) using default MHD output'); end
%np =5; sfopt = 0; norm = 0; mhd_out = 'y';
end
    
[r0, abxyz] = load_sph_coeffs(sfopt);

% debug
disp(['   Spherical r0 (x/y/z in mm) =  ', num2str(r0.x), ' / ', num2str(r0.y), ' / ', num2str(r0.z)]);
% disp(['   SPH COEFFS AX =  ']);
% abxyz.ax
% disp(['   SPH COEFFS BX =  ']);
% abxyz.bx
% disp(['   SPH COEFFS AY =  ']);
% abxyz.ay
% disp(['   SPH COEFFS BY =  ']);
% abxyz.by
% disp(['   SPH COEFFS AZ =  ']);
% abxyz.az
% disp(['   SPH COEFFS BZ =  ']);
% abxyz.bz

% adjust X and Z map resolution so that their np are both intergers, if
% their characteristic radii are different (e.g. Siemens systems).
if (~isstruct(r0)) r00=r0; clear r0; r0.x=r00; r0.y=r00; r0.z=r00; end  % MD 20170622 for back-compatibility

if (r0.x ~= r0.z)    
    if (r0.x > r0.z) [np1 np3 r0.x r0.z] = gnl_map_findsize(r0.x, r0.z); end
    if (r0.x < r0.z) [np3 np1 r0.z r0.x] = gnl_map_findsize(r0.x, r0.z); end   
    np2 = np1; r0.y = r0.x;
    disp( ' (FYI) the number of points (np) have been adjusted in final field maps ....');
    disp([' (FYI) np1 used for computation = ',  num2str(np1)]);
    disp([' (FYI) np2 used for computation = ',  num2str(np2)]);
    disp([' (FYI) np3 used for computation = ',  num2str(np3)]);
    disp([' (FYI) the final output map with np along one dim = ',  num2str(min([np1,np2,np3]))]);
else
    np1 = np;
    np2 = np;
    np3 = np;
end

% SPH coeffs structure (abxyz) containing a_nm and b_nm for gradient X, Y and
% Z -coils, along with Bzmap's resolution info (fov and num_pts).
nar = 1:length(fieldnames(abxyz.ax));

% define a data cube (np x np x np)
a = linspace(-1.0, 1.0, np1);
b = linspace(-1.0, 1.0, np2);
c = linspace(-1.0, 1.0, np3);

[X1,Y1,Z1] = meshgrid(a, a, a);    % X/Y/Z-coordinates of cube
[X2,Y2,Z2] = meshgrid(b, b, b);    
[X3,Y3,Z3] = meshgrid(c, c, c);    

% Note, magnet coordinate system is different from math definition, i.e. 
% in math world X axis is a horizontal line from left to right; while 
% in (Philips) dictionary, gradient field X direction is from bottom to 
% up (AP equvalent to math Y axis )    
[THETA1, PHI1, Rsph1] = cart2sph(Y1(:),X1(:),Z1(:)); 
[THETA2, PHI2, Rsph2] = cart2sph(Y2(:),X2(:),Z2(:)); 
[THETA3, PHI3, Rsph3] = cart2sph(Y3(:),X3(:),Z3(:)); 

% convert to spherical coordinate: Azimuth (THETA) is the counterclockwise angle (in rad) 
% in the x-y plane measured from the positive x-axis. Elevation (PHI) is the angle from 
% the x-y plane. The ranges of THETA and PHI are redefined for Physics definations.
THETA1 = 1.0*pi + THETA1;       % [-1.0*pi, +1.0*pi] ===> [0, 2.0*pi]
PHI1   = 0.5*pi - PHI1;         % [ 0.5*pi, -0.5*pi] ===> [0, 1.0*pi]   

THETA2 = 1.0*pi + THETA2;       % [-1.0*pi, +1.0*pi] ===> [0, 2.0*pi]
PHI2   = 0.5*pi - PHI2;         % [ 0.5*pi, -0.5*pi] ===> [0, 1.0*pi]   

THETA3 = 1.0*pi + THETA3;       % [-1.0*pi, +1.0*pi] ===> [0, 2.0*pi]
PHI3   = 0.5*pi - PHI3;         % [ 0.5*pi, -0.5*pi] ===> [0, 1.0*pi]   

% get the field names
axfn = fieldnames(abxyz.ax);
bxfn = fieldnames(abxyz.by); 

% global field scaling factors
c11x = getfield(abxyz.ax, char(axfn(1))); sclx = r0.x/max([c11x(2),eps]);
s11y = getfield(abxyz.by, char(bxfn(1))); scly = r0.y/max([s11y(2),eps]);
c10z = getfield(abxyz.az, char(axfn(1))); sclz = r0.z/max([c10z(1),eps]);

% predefine 1D field maps
BzXYZ.zx = zeros(1,np1*np1*np1);  
BzXYZ.zy = zeros(1,np2*np2*np2); 
BzXYZ.zz = zeros(1,np3*np3*np3);

tic;
% summation all L degrees SPH 
for iin = 1:length(nar)      
    L = nar(iin); M = 0:L;  
    
    % retrive L+1 SPH coeffs for L order functions
    ax = getfield(abxyz.ax, char(axfn(L))); bx = getfield(abxyz.bx, char(bxfn(L)));             
    ay = getfield(abxyz.ay, char(axfn(L))); by = getfield(abxyz.by, char(bxfn(L)));             
    az = getfield(abxyz.az, char(axfn(L))); bz = getfield(abxyz.bz, char(bxfn(L)));    
   
    % SPH coeffs, matrix size = [L+1, np]
    CmnX = repmat(ax',1,np1*np1*np1); CmnY = repmat(ay',1,np2*np2*np2); CmnZ = repmat(az',1,np3*np3*np3);  
    SmnX = repmat(bx',1,np1*np1*np1); SmnY = repmat(by',1,np2*np2*np2); SmnZ = repmat(bz',1,np3*np3*np3);  
      
    % associated legendre function, matrix size = [L+1, np]
    Lmn1   = legendre(L,cos(PHI1));  
    Lmn2   = legendre(L,cos(PHI2)); 
    Lmn3   = legendre(L,cos(PHI3)); 
    
    % radial function, matrix size = [L+1, np]
    Rmn1 = repmat(transpose(Rsph1).^L, L+1, 1);
    Rmn2 = repmat(transpose(Rsph2).^L, L+1, 1);
    Rmn3 = repmat(transpose(Rsph3).^L, L+1, 1);
      
    % Azimuth angular function, matrix size = [L+1, np]
    Mmn1 = repmat(M', 1, np1*np1*np1);
    Amn1 = repmat(THETA1', L+1, 1); 
    
    Mmn2 = repmat(M', 1, np2*np2*np2);
    Amn2 = repmat(THETA2', L+1, 1); 
    
    Mmn3 = repmat(M', 1, np3*np3*np3);
    Amn3 = repmat(THETA3', L+1, 1); 
    
    % summation all M degrees SPH for one order   
    BzXYZ.zx = BzXYZ.zx + sum(Rmn1.*Lmn1.*(CmnX.*cos(Mmn1.*Amn1) + SmnX.*sin(Mmn1.*Amn1)));    
    BzXYZ.zy = BzXYZ.zy + sum(Rmn2.*Lmn2.*(CmnY.*cos(Mmn2.*Amn2) + SmnY.*sin(Mmn2.*Amn2)));
    BzXYZ.zz = BzXYZ.zz + sum(Rmn3.*Lmn3.*(CmnZ.*cos(Mmn3.*Amn3) + SmnZ.*sin(Mmn3.*Amn3)));     
    
    % update computation status
    disp(sprintf('Computation is in progress ... %5.1f %s.', 100.0*iin/length(nar), '%'));
end
toc;

% apply global normalization leading to dBxyz/dxyz dimensionless
if norm == 0
    sclx = 1.0;
    scly = 1.0;
    sclz = 1.0;
end 
BzXYZ.zx = sclx*BzXYZ.zx;
BzXYZ.zy = scly*BzXYZ.zy;
BzXYZ.zz = sclz*BzXYZ.zz;

% keep only valid SPH region, i.e. r <= r0
invalid_id1  = Rsph1' > 1.000; 
invalid_id2  = Rsph2' > 1.000; 
invalid_id3  = Rsph3' > 1.000; 
BzXYZ.zx(invalid_id1) = NaN;
BzXYZ.zy(invalid_id2) = NaN;
BzXYZ.zz(invalid_id3) = NaN;

% start to save indivual 3D volumetric map (Bx/By/Bz) in MHD/RAW format
tmp01 = reshape(BzXYZ.zx, np1, np1, np1);
tmp02 = reshape(BzXYZ.zy, np2, np2, np2); 
tmp03 = reshape(BzXYZ.zz, np3, np3, np3);

% keep the smaller volume as a reference 3D size 
npmin = min([np1,np2,np3]);
npmax = max([np1,np2,np3]);
widx  = (npmax-npmin)/2;
idst = widx + 1;
ided = widx + npmin;
fov  = 2.0*r0.x;  
dim  = np1; 
if (np1 == np3)
    allbs(:,:,:,1) = tmp01;
    allbs(:,:,:,2) = tmp02;
    allbs(:,:,:,3) = tmp03;
end
if (np1 > np3)
    allbs(:,:,:,1) = tmp01(idst:ided, idst:ided, idst:ided);
    allbs(:,:,:,2) = tmp02(idst:ided, idst:ided, idst:ided);
    allbs(:,:,:,3) = tmp03;
    fov = 2.0*r0.z;
    dim = np3; 
end
if (np1 < np3)
    allbs(:,:,:,1) = tmp01;
    allbs(:,:,:,2) = tmp02;
    allbs(:,:,:,3) = tmp03(idst:ided, idst:ided, idst:ided); 
end

% ===================================================
% First, start to save the field map in matlab format
% ===================================================
% BzXYZ_name = ['Gfield_BzXYZ_',datestr(now,'yymmddHHMM'),'.mat'];


% MD 20180212: re-assign GAP-->GY, GRL-->GX (to align output with DICOM DWI)
allbs = allbs(:,:,:,[2 1 3]); 
% size(allbs) %debug
[d1z, d2z, d3z, d4z] = size(allbs); i0 = round(d1z/2); j0 = round(d2z/2); k0 = round(d3z/2);

% MD 20180910: generalize XY permute for single gradient (RL or AP) when only one SPH provided
grld1 = abs(squeeze(allbs(end,j0,k0,1))-squeeze(allbs(1,j0,k0,1)));
grld2 = abs(squeeze(allbs(i0,end,k0,1))-squeeze(allbs(i0,1,k0,1)));
gapd1 = abs(squeeze(allbs(end,j0,k0,2))-squeeze(allbs(1,j0,k0,2)));
gapd2 = abs(squeeze(allbs(i0,end,k0,2))-squeeze(allbs(i0,1,k0,2)));
%[abs(squeeze(allbs(end,j0,k0,1))-squeeze(allbs(1,j0,k0,1))) abs(squeeze(allbs(i0,end,k0,1))-squeeze(allbs(i0,1,k0,1)))]
%[grld1 grld2 gapd1 gapd2] %debug

if ((grld1<grld2) || (gapd2<gapd1)) % MD 20180910: permute only when horiz-G along dim2 or vert-G along dim1
  %if (abs(allbs(end,j0,k0,1)-allbs(1,j0,k0,1)) < abs(allbs(i0,end,k0,1)-allbs(i0,1,k0,1)))
  % MD 20180907: permute only when (1st) horiz "gradient" is along dim2 (or leave as is for Siemnes)
allbs = permute(allbs,[2 1 3 4]);
tnp = np1; np1 = np2; np2 = tnp;
end % check permute

BzXYZ.zx = squeeze(allbs(:,:,:,1));
BzXYZ.zy = squeeze(allbs(:,:,:,2));
BzXYZ.zz = squeeze(allbs(:,:,:,3));

timestamp = datestr(now,'yymmddHHMM'); % MD 20151117 -- time-stamp for the output file names
BzXYZ_name = ['Gfield_BzXYZ_', timestamp,'.mat']; % MD 20151117: same time-stamp from MHD and mat-output
save(BzXYZ_name, 'BzXYZ');
disp(['  FYI MAT-File output saved: ', BzXYZ_name]);

% ===================================================
% Second, start to save the field map in MHD format
% ===================================================
if (isequal(mhd_out,'y'))
	allbs(isnan(allbs(:))) = 0;              % replace NaN with ZERO
	outfile_prefix = ['Gfield_' timestamp];  % MHD output labels
	xcomment = 'Simulation';
	volabels = {'BzGx','BzGy', 'BzGz'};

% start to save indivual 3D volumetric map (Bx/By/Bz) in MHD/RAW format
	%del1 = 2.0*r0.x/(np1-1.0);      % pixel resolution - distance between two pixel centers;
	%del2 = 2.0*r0.y/(np2-1.0);
%MD 20180215: change mhd-pars to reflect new output order X=RL; Y=AP
	del1 = 2.0*r0.y/(np1-1.0);      % pixel resolution - distance between two pixel centers;
	del2 = 2.0*r0.x/(np2-1.0); % np1/2 already permuted above
	fov1 = fov;
    fov2 = fov;
    fov3 = fov;
	dim1 = dim;
    dim2 = dim;
    dim3 = dim; 
    dim4 = 3;
	orax = repmat([1 0 0 0 1 0]', 1, dim3);  % AX orientation
	locs = [repmat(-fov1/2,1,dim3); repmat(-fov1/2,1,dim3); linspace(-fov3/2, fov3/2, dim3)]; % slice locations
	npmhd = mk_mininifti_params(struct('iminf', struct('PixelSpacing', [del1; del2], ...
  'Rows', dim1, 'Columns', dim2), ...
  'nipar', struct('loc', locs,'orient', orax), ...
   'geo', [fov1 fov2 dim3 dim4]));
%MD 20180215: ???? >>>> MHD vs raw should permute 1/2 dims?: permute(allbs,[2 1 3 4])
allbs = flipdim(permute(allbs, [2 1 3 4]),3); % MD20180223: convert from MAT(magnet)-frame to LPS
	build_mhd_v5(allbs, struct('nipar', npmhd, 'fipref', outfile_prefix, 'vln', 's','sidsls', ...
'n','prcn', 'MET_FLOAT','xcomnt',xcomment,'vlbs', {volabels} ));
end % if out_mhd

return  % END of "gbzmap_main"

% ------------------------------------------------------------------------------
function [n0 n1 r0 r1] = gnl_map_findsize(r0, r1)
% ------------------------------------------------------------------------------
% find integer number of points to sample large (2*r0) and small (2*r1) FOV
% with the same resolution (+/-0.01mm); and minimum number of points > 110
% (resolution < 5mm)
% [r] = mm units
% would work for reasonable chracteritsic MRI radii of r0=200-300 mm
% MD 20170622:  added soft lower point numbers fro lardest radius n0 > (110/2)*r0/r1
%		so that for smallest n1 > 110
%
r0 = round(r0);   r1 = round(r1);    
if (r1 > r0) r1t = r1; r1 = r0; r0 = r1t; end % find r0>r1;   
r0l = r0-1; r0h = r0+1; % low and high boundaries within 1mm

kr0 = floor(55*r0/r1):r0; % minimum number of points for largest radius MD 20170622: (110/2)*r0/r1
kk0 = find(abs(kr0*r1/r0 - fix(kr0*r1/r0))==0);   kk = kk0;
kk0l = find(abs(kr0*r1/r0l - fix(kr0*r1/r0l))==0); 
kk0h = find(abs(kr0*r1/r0h - fix(kr0*r1/r0h))==0);
[mxl, mxli] =  max([length(kk0l), length(kk0), length(kk0h)]); 
%[length(kk0l), length(kk0), length(kk0h)] %debug
if (mxli==1)  
    kk = kk0l; r0 = r0l;   
elseif (mxli==3) 
    kk = kk0h; r0=r0h;   
end   % adjust largest radius by 1mm if needed

if length(kk)>1 
      k0 = min(kr0(kk)); 
      k1=k0*r1/r0;  
else  % no solution -- set "res"=2mm -- > adjust radii to be "even"
     r1 = r1+1*(abs(r1/2-fix(r1/2))>0);  r0 = r0+1*(abs(r0/2-fix(r0/2))>0); 
     %[r0 r1] %debug
     k0 = r0/2; k1 = r1/2;   
end

n0 = 2*k0+1;  n1 = 2*k1+1; 

% res0 = 2*r0/(n0-1); res1 = 2*r1/(n1-1);

%[res0*(n1-1)/2 r1] % debug

%[[res0 res1];[round(200*r0/(n0-1))/100  round(200*r1/(n1-1))/100]] % debug

return



% ------------------------------------------------------------------------------
function [r0, C] = load_sph_coeffs(SCALE_OPT, FileName, N_MAX) 
% ------------------------------------------------------------------------------
% PURPOSE:
% 	load gradient coils' spherical harmonic (SPH) coefficents and global 
%   geometry information to define Bzmap resolution (FOV/NPT) from a text file
%
% INPUT:
% 	N_MAX      		int    (optinal) SPH degree, N_MAX=16 by default
%
% 	FileName   	 string    (optinal) a text file.
% 
% OUTPUT: 
%       r0        float    SPH reference radius in mm.
%        C    structure	   A MATLAB structure that contains values read from the file in  
%                          named fields, i.e. 'ax', 'bx', 'ay', 'by', 'az', 'bz'. 
%
% 
% MODIFICATION HISTORY:
% 	17/09/2015: 	v1.0, initial draft
%   23/09/2015:     v1.1, added r0 output argument
%   06/07/2017:     v1.2, change r0 to r0.x, r0.y, r0.z for various radius 
%                         per gradient channel(e.g. in Siemens case)
% MD 20170622: 	   changed default r0=0; added a check for missing RX/RY/RZ
% MD 20180212:     relabel gradient directions (in demos) X<->AP, Y<->RL, Z<->SI
% YP 20180821:     a. bug fixed with "cellstr" to convert character array into cell array.
%                  b. max SPH order (N_MAX) is determined from user's input
%                  rather than to be predefined (e.g. default N_MAX=16). 
%                
% ------------------------------------------------------------------------------
if (nargin < 2)
    [FileName,PathName] = uigetfile('*.txt','Select a text file containing SPH coefficents ...');
    FileName = [PathName, FileName];      
end
fid = fopen(FileName,'r');
if (length(fid) < 1)
    disp('Error: Input is not a valid text file!')
    return
end

% YP 20180821, determine the max SPH order: N_MAX
N_MAX = 1; 
while ~feof(fid)
    curline = strtrim(fgetl(fid));  
    if (isempty(curline) || curline(1) ~= 'G'),  continue, end
    Temp = cellstr(strtrim(strsplit1(curline,'=')));
    JUNK = cellstr(strtrim(strsplit1(Temp{1},'_')));  % retrive DEG and ORD value
    N_DEG = int8(str2num(JUNK{5}));
    if N_DEG > N_MAX 
        N_MAX =  N_DEG;
    end
end
fclose(fid);
disp([' -----------------------------------------------------'])
disp(['   INPUT SPH MAX ORDER: ', num2str(N_MAX)])
disp([' -----------------------------------------------------'])

% prepare the data structure template 
for N = 1:N_MAX
    ax.(genvarname(['a' num2str(N,'%02d')]))  = zeros(1,N+1);
    bx.(genvarname(['b' num2str(N,'%02d')]))  = zeros(1,N+1);
    ay.(genvarname(['a' num2str(N,'%02d')]))  = zeros(1,N+1);
    by.(genvarname(['b' num2str(N,'%02d')]))  = zeros(1,N+1);
    az.(genvarname(['a' num2str(N,'%02d')]))  = zeros(1,N+1);
    bz.(genvarname(['b' num2str(N,'%02d')]))  = zeros(1,N+1);
end
C = struct('ax', ax, 'bx', bx, 'ay', ay, 'by', by, 'az', az, 'bz', bz); 
% find out global and local SPH scaling factor and option

r0.x  = 0.0;  % default r0 -- MD20170622: re-assigned to "zero" to flag if "missing"
r0.y  = 0.0; 
r0.z  = 0.0;
fid = fopen(FileName,'r');
while ~feof(fid)
    curline = strtrim(fgetl(fid));  
    if (isempty(curline) || curline(1) == '#' || curline(1) == 'G'), continue, end 
    Temp = cellstr(strtrim(strsplit1(curline,'=')));     % retrive numerical data 
    rad  = Temp{1}; 
    var  = str2num(Temp{2}); 
    if strcmp(rad(end-1:end), 'MM')  % back-compatiable with uniform r0
        r0.x = var;
        r0.y = var;
        r0.z = var;
    else % MD 20180212: handle relabeled demos X<->AP, Y<->RL, Z<->SI 
        if strcmp(rad(end-2:end), '_AP'), r0.x = var;, end
        if strcmp(rad(end-2:end), '_RL'), r0.y = var;, end
        if strcmp(rad(end-2:end), '_SI'), r0.z = var;, end
    end
end
fclose(fid);

% MD 20170622: check that radii fields are not missing and warn as needed
RX = tryGetField(r0, 'x',0.0); RY = tryGetField(r0, 'y',0.0); RZ = tryGetField(r0, 'z',0.0);
r0.x = RX*(RX>=1) + RY*(RX<1); r0.y = RY*(RY>=1) + RX*(RY<1); % use RX=RY if one of them is missing

if ((RX<1) & (RY<1)) % when both RX and RY are missing
rxy = max([RZ 1.0]);
disp(['WARNING: Rxy not provided -- check your input -- "Rxy" will be set = Rz : ' num2str(rxy) 'mm ...']); 
r0.x = rxy; r0.y = rxy;
end

if (RZ<1) 
r0.z = max([r0.x r0.y 1.0]); 
disp(['WARNING: Rz not provided -- check your input -- will be set = Rxy : ' num2str(r0.z) 'mm ...']); 
end


% load in local SPH coeffs and apply global NORMfactor or local scaling if
% requested.
fid = fopen(FileName,'r');
while ~feof(fid)
    curline = strtrim(fgetl(fid));  
    if (isempty(curline) || curline(1) ~= 'G'),  continue, end 
   % Gxyz = curline(1:2);   % GX | GY | GZ  
   % CnSn = curline(10);    % C | S
% MD 20180212: relabel X<->AP, Y<->RL, Z<->SI 
    Gxyz = curline(1:3);   % GAP | GRL | GSI  
    CnSn = curline(11);    % C | S
    Temp = cellstr(strtrim(strsplit1(curline,'=')));  % retrive numerical data
    COEFFS = str2num(Temp{2});  
 
    % apply local SPH coeffs scaling factor
    JUNK = cellstr(strtrim(strsplit1(Temp{1},'_')));  % retrive DEG and ORD value
    N_DEG = int8(str2num(JUNK{5}));
    M_ORD = int8(str2num(JUNK{6}));     
    SF = scale_sph_coeffs(SCALE_OPT, N_DEG, M_ORD);    
    % if requested (i.e. SCALE_OPT == 1), order-degree specific SPH NORM 
    % factor wiil be appled; otherwise SF = 1.0
%MD20180212: get gradient label X<->AP, Y<->RL, Z<->SI 
    gl = char('x'*isequal(Gxyz(2),'A')+'y'*isequal(Gxyz(2),'R')+'z'*isequal(Gxyz(2),'S')); 
    if strcmp(CnSn,'C')
       % C.(genvarname(['a',lower(Gxyz(2))])).(genvarname(['a' num2str(N_DEG,'%02d')]))(M_ORD+1) =  COEFFS * SF; 
	C.(genvarname(['a',gl])).(genvarname(['a' num2str(N_DEG,'%02d')]))(M_ORD+1) =  COEFFS * SF;       
    else 
       % C.(genvarname(['b',lower(Gxyz(2))])).(genvarname(['b' num2str(N_DEG,'%02d')]))(M_ORD+1) =  COEFFS * SF;
	C.(genvarname(['b',gl])).(genvarname(['b' num2str(N_DEG,'%02d')]))(M_ORD+1) =  COEFFS * SF;
    end    
end

fclose(fid);
return % load coeffs

% ------------------------------------------------------------------------------
function sf = scale_sph_coeffs(scale_opt, n_deg, m_ord) 
% ------------------------------------------------------------------------------
% PURPOSE:
% 	computate SPH normalization factor, see ref. Equation (6.8.2) from Page
% 	252 in "Numerical Recipes in C: The Art of Scientific Computing, Second
%   Edition"
%
% INPUT:
% 	scale_opt      	int    toggle to turn on or off SPH norm scaling factor
%                          calculation, i.e. zero means OFF, ONE means ON,
%                          other options could be included later on if
%                          needed
%       n_deg   	int    The degree of Spherical Harmonic
%       m_ord   	int    The order of a certain degree of Spherical Harmonic
% 
% OUTPUT: 
%      sf         float    SPH coeffs scaling factor
%
% MODIFICATION HISTORY:
% 	17/09/2015: 	v1.0, initial draft
%
% MD 20190828: added scale-option=2 for Siemens convention norm(n,0)=1; norm(1,1) = 1, 
%		and normSPH(n,m)=sqrt(2*pi)*normSPH(n,m)
% ------------------------------------------------------------------------------
if (nargin == 0)
    sf = 1.0; 
    disp(['SPH coeffs scaling factor set to:  ', num2str(sf)])
    return
elseif (nargin == 3)
    if m_ord > n_deg
        sf = 1.0;
        disp(['ERROR: SPH COEFFS ORDER (M = ', num2str(m_ord),') ', ...
              'should be SMALLER than DEGREE (N = ',num2str(n_deg),').'])
        return
    end
end

n_deg = double(n_deg); %single(n_deg); % MD 20151116 -- Matlab V<14 assumes "integer" args for factorial
m_ord = double(m_ord); %single(m_ord);
%[n_deg m_ord n_deg-m_ord n_deg+m_ord]
SPH_SCALE_FACTOR  = sqrt( ((2.0*n_deg+1.0)/(4.0*pi)) * (factorial(n_deg - m_ord))/ (factorial(n_deg + m_ord)) );
% UNK_SCALE_FACTOR  = 1.0;
switch scale_opt
    case 0
        sf  = 1.0;
    case 1
        sf  = SPH_SCALE_FACTOR;
    case 2
         sf  = SPH_SCALE_FACTOR*sqrt(2*pi);
	if (m_ord == 0)  % YP 20190821
    		sf = 1.0; 
	elseif (m_ord == 1 & n_deg == 1) %MD 20190828
		sf = 1.0;
	end
%   case 3
%         sf  = SPH_SCALE_FACTOR * UNK_SCALE_FACTOR;       
    otherwise
        sf  = 1.0;       
end % end of switch

return % scale coeffs

function [ nifti_params ] = mk_mininifti_params(nipin) %imginfo,nip,fov1,fov2,dim3,dim4)
% 
% To make "nifit_params" with minimal (MHD-like) information, eg, for
% "simulation" or mat-volumes
% follows "build_nifti_params" structure (original called by readdicom7 and readdicom7_allscan to
% extract key acquistion parameters from dicom for retention and use in mhd
% and nifti files.)
% 20150629  TLChenevert.  Plucked code out of readdicom7 (identical in
% readdicom7_allscan).
% 20150630 TLC
% 20150702 TLC. Seems functional, clean-out deadwood.
%
% 20150106 Start ...
% Currently nifti_params.orientation used to create nifti will be saved in one structure, but
% these WILL BE IN ERROR FOR MULTI-ANGLE-STACKS SERIES - fortunately a rare scientific need.
% 20150706 TLC Put in more checks against missing entities:
% "MagneticFieldStrength", "ExamDescription", "SeriesDescription".
% 
% 20150916 Start
% 20150916 "nipin" input structure mimics required "imginfo" and "nip" fields:
%          'iminf' = imginfo; 'nipar' = nip; 'geo' = [fov1 fov2 dim3 dim4];
%           NOTE: 'geo' is (required) for meaningful MHD
%
%           utilised "imginfo" fields (usually from dicom-info) 
%           (optional): 
%                     'SeriesNumber', 'SeriesDescription', 'SliceThickness', 'MagneticFieldStrength'
%                     'InstitutionName', 'Manufacturer', 'SoftwareVersion','MRAcquisitionType'
%                     'ManufacturerModelName', 'DeviceSerialNumber', 'StudyDescription'
%                     'InPlanePhaseEncodingDirection', 'ImageType', 'StationName'
%           (mhd-required):
%                     'PixelSpacing', 'Rows', 'Columns'
%
%           utilised "nip" fields (see readdicom7)
%           (optional):
%                       'td', 'ppd', 'naverages', 'echotime', 'bvalue',
%                       'diffdir', 'imtyp'
%           (mhd-required):
%                       'loc', 'orient'
%
% Full info call (eg, from within "readdicom7"):
%           nifti_params = mk_mininifti_params(struct('iminf', imginfo, 'nipar', nip, 'geo', [fov1 fov2 dim3 dim4]));
%
% 20150917 MHD-sufficient call (eg, from gradient field generator):
%           nifti_params = mk_mininifti_params(struct('iminf', struct('PixelSpacing', [del1; del2], 'Rows', dim1, 'Columns', dim2), 'nipar', struct('loc', slslocs,'orient', slscos), 'geo', [fov1 fov2 dim3 dim4]))
%

axor = [1 0 0 0 1 0]'; % AX orientation

if ~exist('nipin', 'var') nipin = []; end

if isempty(nipin) nipin = struct('iminf', [], 'nipar', [], 'geo', []); end 
% init general structures
if isfield(nipin, 'iminf') imginfo = nipin.iminf; else imginfo = []; end
if isfield(nipin, 'nipar') nip = nipin.nipar; else nip = []; end
if isfield(nipin, 'geo') geo = nipin.geo; else geo = [0 0 0 0]; end

%%% fields from "imginfo" (usually DICOM-info, but can be specified
%%% independently)

% define "parallel" (opt for mgh) field names in "imginfo" versus "nifti-params"
% NOTE: keep them at the same index-position in cell-arrays to streamline
% next call
npoflds = {'seriesnumber', 'seriesdescription', 'slicethickness', 'mracqtype','inplanephencdir','imgtypstring', 'institution', 'manufacturer'};
npoflds = [npoflds, {'model','fieldstrength', 'stationname', 'serialno', 'swversion', 'studydesc'}];
dioflds = {'SeriesNumber', 'SeriesDescription', 'SliceThickness', 'MRAcquisitionType', 'InPlanePhaseEncodingDirection', 'ImageType', 'InstitutionName', 'Manufacturer'};
dioflds = [dioflds, {'ManufacturerModelName','MagneticFieldStrength','StationName', 'DeviceSerialNumber', 'SoftwareVersion', 'StudyDescription' }];

for inof = 1:length(dioflds) % assign a bunch of optional fields
    nifti_params.(npoflds{inof}) = tryGetField(imginfo,dioflds(inof),'UNK'); 
end
%     nifti_params.seriesnumber = tryGetField(imginfo,'SeriesNumber','UNK');
% %nifti_params.seriesdescription = imginfo.SeriesDescription;
%     nifti_params.seriesdescription = tryGetField(imginfo,'SeriesDescription', 'UNK');
%     nifti_params.slicethickness = tryGetField(imginfo,'SliceThickness', 'UNK'); % 20150204 returns "[]" if no field present
%     nifti_params.mracqtype = tryGetField(imginfo,'MRAcquisitionType', 'UNK'); % 20150204 returns "[]" if no field present
     nifti_params.pixelspacing = tryGetField(imginfo,'PixelSpacing', 0); % 20150204 returns "[]" if no field present
%     nifti_params.inplanephencdir = tryGetField(imginfo,'InPlanePhaseEncodingDirection', 'UNK'); % 20150204 returns "[]" if no field present
%     nifti_params.imgtypstring = tryGetField(imginfo, 'ImageType', 'UNK'); % 20150308
     nifti_params.rows = tryGetField(imginfo, 'Rows', 0); % 20150308
     nifti_params.columns = tryGetField(imginfo, 'Columns', 0); % 20150308
%     % Institutional, system, and exam demographics:
%     % nifti_params.institution = imginfo.InstitutionName;
%     nifti_params.institution = tryGetField(imginfo,'InstitutionName','UNK');
%     %nifti_params.manufacturer = imginfo.Manufacturer;
%     nifti_params.manufacturer = tryGetField(imginfo,'Manufacturer','UNK');
%     %nifti_params.model = imginfo.ManufacturerModelName;
%     nifti_params.model = tryGetField(imginfo,'ManufacturerModelName','UNK');
%     %nifti_params.fieldstrength = imginfo.MagneticFieldStrength;
%     nifti_params.fieldstrength = tryGetField(imginfo,'MagneticFieldStrength',0);
%     nifti_params.stationname = tryGetField(imginfo,'StationName','UNK');
%     %nifti_params.serialno = imginfo.DeviceSerialNumber;
%     nifti_params.serialno = tryGetField(imginfo,'DeviceSerialNumber','UNK');
%     %nifti_params.swversion = imginfo.SoftwareVersion;
%     nifti_params.swversion = tryGetField(imginfo,'SoftwareVersion','UNK');
%     %nifti_params.studydesc = imginfo.StudyDescription;
%     nifti_params.studydesc = tryGetField(imginfo,'StudyDescription','UNK');
%%% fileds from "imginfo"

if ~isfield(imginfo, 'SeriesNumber') imginfo.SeriesNumber = '000'; end % for non-dicom imginfo
%nifti_params.defname = make3TRID(imginfo,'di'); % 'di' will use de-identification xform.  Only alternative is 'id' that
nifti_params.defname = 'SIMS'; % MD 20151117 -- just name "simulations" here

%%% generally REQUIRED fields for meaningful output
if isempty(nip)
    nip = struct('loc', zeros(3,1), 'orient',zeros(6,1),'td', 0, 'ppd', 0, 'echotime', 0, 'naverages', 0, 'bvalue', 0);  
    nip.diffdir = 0; nip.imgtyp = 0;
end
if ~isfield(nip,'loc') nip.loc = zeros(3,1); end
if ~isfield(nip,'orient') nip.orient = zeros(6,1); end
optflns = {'td', 'ppd', 'echotime', 'naverages', 'bvalue', 'diffdir', 'imgtyp'};
for iof = 1:length(optflns) % assign a bunch of optional fields
    if ~isfield(nip,optflns(iof)) nip.(optflns{iof}) = 0; end
end

if isempty(geo) geo = [0 0 0 0]; end    
fov1 = geo(1); fov2 = geo(2); dim3 = geo(3); dim4 = geo(4);

% check legit MHD geo params:
[ld1 ld2] = size(nip.loc);
%[ld2 dim3]
if (sum(geo)>0) && ~isequal(ld2, dim3) && (dim3 > 1)% slice locations are not properly provided, but
    disp('  WARNING: num-slices vs loc-slices mismatch -- auto-generating "locs" from FOV3=FOV1, assuming AX orientation...');
    fov3 = fov1;
    nip.orient = repmat(axor, 1, dim3);
    slocs = linspace(-fov3/2, fov3/2, dim3);
    clocs = repmat(-fov1/2,1,dim3);
    nip.loc = [clocs; clocs; slocs];
end    
[ld1 ld2] = size(nip.loc);
if ~isequal(ld2, dim3) disp('  WARNING: num-slices vs loc-slices mismatch -- nonlegit geo-params likely...'); end 

nifti_params.locfirstslc = nip.loc(:,1);
nifti_params.loclastslc  = nip.loc(:,end);
nifti_params.slicectr2ctr = sqrt( sum( (nifti_params.locfirstslc - nifti_params.loclastslc).^2) ) / max([1 dim3-1]); % SpacingBetweenSlices
nifti_params.timefirstphase = nip.td(1,1); % Should be (close to) "td" of first phase, 1st slice in imastor structure.
nifti_params.timelastphase  = nip.td(1,end); % Should be (close to) "td" of last phase, 1st slice in imastor structure.
nifti_params.temporalsampling = (nifti_params.timelastphase - nifti_params.timefirstphase) / max([1 dim4-1]); % Meaningless for non-dynamic series, of course.
nifti_params.orientation  = squeeze(nip.orient(:,1)); % Old nifti & mhd not designed for multi-angle stacks.  Keep this info in nifti_params.allloc & .allorientation.

nifti_params.fov = [fov1 fov2];

% % Need to deal with possible replicate locs in 4D.  Select unique slice locs, but do NOT resort via 'unique':
nifti_params.allloc = nip.loc;
nifti_params.allorientation = nip.orient;

% Possible fourth-dimension arrays for retention as an "FYI" in (mhd/nifti) header comment fields.
nifti_params.td = nip.td;   % [nslc nphase] All slice & time-delay info.
nifti_params.ppd = nip.ppd; % [nslc nphase] All slice & alternative time-delay info.
nifti_params.echotime = nip.echotime; % [1 nphase] 1st slice TEs represtative of all slices
nifti_params.naverages = nip.naverages; % [1 nphase] Navs may vary with bvalue, but 1st slice represtative of all. 
nifti_params.bvalue = nip.bvalue; % [1 nphase] 1st slice represtative of all.
nifti_params.diffdir = nip.diffdir; % [4 nphase] 1st slice represtative of all.
nifti_params.imgtyp = nip.imgtyp; % [1 nphase] 1st slice represtative of all.


return % END mk_mini_nifti

function [def3TRID] = make3TRID(dicinfo,diid)
% Create 3T Research ID code for given series based on dicom info.
% diid = 'di' (default) Use de-id transformation and have "di" as 1st characters.
% diid = 'id' Do not transform and have "id" as 1st characters.
% 20150627 TLChenevert

if (nargin < 2)
    diid = 'di';
end

if (diid == 'di')
    xform = 1;
elseif (diid == 'id')
    xform = 0;
else
    disp('Aborting, Second argument should be either "di" or "id". ');
    return;
end

sdate = tryGetField(dicinfo,'SeriesDate','00000000'); % 20150204 returns "00000000" if no field present
fndum1 = str2num(sdate(3:4)) + 3*xform;
if (fndum1 < 10)
   fndum11 = ['0' num2str(fndum1)];
else
   fndum11 = num2str(fndum1);
end
fndum2 = str2num(sdate(5:6)) + 2*xform;
if (fndum2 < 10)
   fndum22 = ['0' num2str(fndum2)];
else
   fndum22 = num2str(fndum2);
end
fndum3 = str2num(sdate(7:8)) + 1*xform;
if (fndum3 < 10)
   fndum33 = ['0' num2str(fndum3)];
else
   fndum33 = num2str(fndum3);
end
% fndum4 = dicinfo.SeriesTime(1:4);
fndum4 = tryGetField(dicinfo,'SeriesTime','0000'); % 20150204 returns "0000" if no field present
fndum4 = fndum4(1:4); % 20150204
fndum5 = dicinfo.SeriesNumber;
%fndum6 = [fndum11 fndum22 fndum33 num2str(fndum4) 's' num2str(fndum5) '_0001'];
fndum6 = [diid fndum11 fndum22 fndum33 num2str(fndum4) 's' num2str(fndum5)];
def3TRID = fndum6;

%end
return % END "make3TRID"


function [ output_args ] = build_mhd_v5(imdat,mhdpars) %nifti_params,fileprefix,voln,browz,nvect,externalcomment)
% build_mhd_v5 calls write_mhd2.m and mk_ser4Dlbl.m (opt)
% script primarily based on "v4" written by TLC (as of 08/25/2015)
%      
%   Inputs:
%       0)  If no input arguments, user prompted to browse to and select a dicom series folder.
%       1)  (semi-mandatory) imdat = structure from readdicom7 or is 2D-3D-4D image array, BUT YOU MUST also provide corresponding nifti_params.
%       2)  (semi-mandatory) mhdpars = structure with build parametrs and options (eg, generated by "mk_mhd_inpars.m")
%                           can also use "nifti_params" structure (created > July3, 2015) or full "ExamSeries" (from ExamDemographics) and you will be prompted to select desired series.
%                           NOTE: if "imdat" is not specified, reading data
%                                 from DICOM will automatically assign/overwrite
%                                 nifti-params from "mhdpars"
%  MD:20150828              "mhdpars" structure contains 'fields' for:
%         (semi-mandatory) 'nipar' = nifti_params structure 
%         (optional) 'fipref' = fileprefix string as fix to default output mhd filename
%         (optional) 'vln' = 1,2,3... or 'a' or 'A' selects single
%         vector-volume of 4D image or ALL 3D volumes; 's' or 'S' to use volume-label as a suffix to file-name
%         (optional) 'brwz' = browse 0 readdicom7 search path starts from: 0=SimpleDicom, 1=Last DicSeries, 2=Within Last DicSeries, 3=Current Dir.
%         (optional) 'nvec' = nvect 'n' (def) or 'y' to write-out N separate 3D volumes (def for 3D-Slicer) or as single N-vector 4D in one file (for 4D-Slicer).
%         (optional) 'xcomnt' = externalcomment string added to end of "extracomment" cell array used as FYI-info in mhd header.
%         (optional) 'vlbs' = {'v1' 'v2' 'v3'..} 4th dim volume label cell-array, or a string 'v1' for a single volume, 
%                               also allowed to use 'i' - to ignore, 'c' - for input, 'g' - auto-generate 
%         (optional) 'crop' = 'y'
%         (optional) 'sidsls' = 'n' see-mid-slice for volume selection
%
% The header (mhd) info and binary data are extracted from the output of readdicom7.m
% Specifically, image orientation (patient) and image position (patient)
% are taken from nifit_params structure, these infomation are critical.
% 
% The binary data (raw) in floating point format is transposed for the first
% two dimensions, making it compatible with the way of how DICOM data is stored.
%
% This script serves as an example on how to use the WRITE_MHD_RAW.m.
%
% Written by Yuxi Pang, 01-26-2015, Univ of Michigan Hospital.
% 20150611 TLC Add more scaninfo-related material to comment field of MHD header
% 20150626 TLC Discovered Siemens 3D acqs may not provide
% "SpacingBetweenSlices"; therefore created work-around.
% 20150625 TLC: "v3" 3D-4D options, as well as use input args so not necessarily starting from dicom.
% 20150627 TLC: "v4" to utilize optional input arguments:
%               imdat = 2D(?), 3D, 4D data array or full imdat structure. Either way nifti_params still required.
%               If no input args, then prompt user to read in dicom.
% [imdat,dim1,dim2,dim3,dim4,fov1,fov2,fov3,fov4,tr,unique_te,ppd,flip,ii,nifti_params] = readdicom7(0,3);
% 20150723 TLC: Include "CenterOfRotation" and "AnatomicalOrientation" output by Elastix.
%               Also allow single-precision storage, default for Elastix.
%
% 20150826 MD: recast input params into "mhdpars" structure
%              NOTE: to add more parameters, add fields and INIT statements
% 20150827 MD: added "extracomment{12}" label for the 4th dim vector
% 20150828 MD: mhd-build-pars will be saved with the base-file-name suffix
%             (in addition to being provided as the output)
% 20150831 MD: enable voln='s' to add volume label to file-name suffix 
% 20150831 MD: clean special characters from "suffix"
% 20150902 MD: allow 1-line "v4"-type input params as "struct('field-name','field-value')", e.g.:
%              build_mhd_v5([], struct('fipref','1081_20','vln', 'c','xcomnt','ADCb0b2k','vlbs', 'i' ));
% 20151023 MD: using try GetFiled, assign single 3D vol label as in "vln" to avoid repeat
% 20151023 MD: reset "deafult" for "seemidslice" to 'n'
%
% (fyi): Example of calls from Matlab workspace:
%       single 3D vol series: 
%              >> build_mhd_v([], struct('fipref', subj_label, 'vln', file_label)); 
%       4D dyn series with auto-labels: 
%              >> build_mhd_v5([], struct('fipref', subj_label, 'vln', 's', 'vlbs', 'g'));
%

% Element Type Precision:
% 'MET_FLOAT' = single precision saves on disk-space; 'MET_DOUBLE' = double precion used in MatLab.
% Regardless of this setting, always convert images to double for MatLab
% calculation.  Also mhd images converted to double upon read-in in "read_mhd_tlc.m".
% precisn = 'MET_DOUBLE'; % keep it simple matlab-like, though bulky.

%precisn = 'MET_FLOAT';  % def

% define cell-array of "special" characters, not allowed in the file-name suffix (for voln = 's' option):
specchar = {':', ';', '"', ' ', '(', ')', '[', ']', '\', '/', '-', '+', '.'};

%%% if not provided, INIT mhd-build-params structure
if (~exist('imdat')) imdat = struct([]); end % no input image data

%mhdpars  %debug

if (~exist('mhdpars'))
%     mhdpars = struct('nipar', struct([]), 'fipref', 'def', 'vln', 'c', 'brwz', 3, 'nvec', 'n', 'xcomnt', [], 'prcn', 'MET_FLOAT', 'vlbs', {});
%     mhdpars.crop = 'y';
%     mhdpars.sidsls = 'y';
    mhdpars = struct([]);
else % check if this the input is straight "nifti_params" or "ExamSeries" 
    if (isstruct(mhdpars) && ~isfield(mhdpars,'nipar')) % no "nipar" field
        temp_nipar = mhdpars;
%         mhdpars.crop = 'y';
%         mhdpars.sidsls = 'y';
        if (isfield(temp_nipar, 'NiftiParams') || isfield(temp_nipar,'allloc')) %ExamSeries or "nifti_params"
            clear mhdpars;
            mhdpars = struct('nipar', [], 'fipref', 'def', 'vln', 'c', 'brwz', 3, 'nvec', 'n', 'xcomnt', [], 'prcn', 'MET_FLOAT', 'vlbs', {'c'}, 'crop', 'y', 'sidsls', 'n');
            mhdpars.nipar = temp_nipar;
            clear temp_nipar;
        else
            mhdpars.nipar = [];
            %disp('WARNING: nonlegit input for "nifti_params" -- will ask to load from file...');
        end       
    end
end

%mhdpars  %debug

if (isempty(mhdpars) || (~isstruct(mhdpars)))
    clear mhdpars;
    mhdpars = struct('nipar', [], 'fipref', 'def', 'vln', 'c', 'brwz', 3, 'nvec', 'n', 'xcomnt', [], 'prcn', 'MET_FLOAT', 'vlbs', {'c'}, 'crop', 'y', 'sidsls', 'n');
%     mhdpars.crop = 'y';
%     mhdpars.sidsls = 'y';
end

%%% INIT mhd-build-params (to provided or defaults); NOTE: add more here as needed
if (isfield(mhdpars, 'prcn')) 
    precisn = mhdpars.prcn; 
else
    precisn = 'MET_FLOAT';
    mhdpars.prcn = precisn; 
end % def data precision

if (isfield(mhdpars, 'vlbs'))
    if (iscell(mhdpars.vlbs))
        volbl = mhdpars.vlbs;
    elseif (ischar(mhdpars.vlbs))
        volbl = {mhdpars.vlbs};
        mhdpars.vlbs = volbl;
    else % non-char label
        volbl = {'i'}; 
        mhdpars.vlbs = volbl; 
    end
else  %MD 20151023: assign single 3D vol label as in "vln" to avoid repeat
    vln1 = tryGetField(mhdpars, 'vln', 'c'); 
    volbl = {vln1}; 
    mhdpars.vlbs = volbl; 
end  % 4th dim volume label(s) (cell-array)
if (~iscell(mhdpars.vlbs) && length(mhdpars.vlbs)==1) volbl = {mhdpars.vlbs}; mhdpars.vlbs = volbl; end

%mhdpars  %debug
%volbl %debug
%mhdpars.vlbs  %debug
if (isfield(mhdpars, 'nipar')) nifti_params = mhdpars.nipar; else nifti_params = []; mhdpars.nipar = nifti_params; end % nifit-parameter structure
if (isfield(mhdpars, 'fipref')) fileprefix = mhdpars.fipref; else fileprefix = 'def'; mhdpars.fipref = fileprefix; end % file prefix for the output
if (isfield(mhdpars, 'vln')) voln = mhdpars.vln; else voln = 'c'; mhdpars.vln = voln; end % default is to prompt user to 'c'hoose volume number
if (isfield(mhdpars, 'brwz')) browz = mhdpars.brwz; else browz = 3; mhdpars.brwz = browz; end % default is to look for DCM from current location
if (isfield(mhdpars, 'nvec')) nvect = mhdpars.nvec; else nvect = 'n'; mhdpars.nvec = nvect; end % 'n' (default) = write out N separate 3D files; 'y' = write out a single 4D file readable by 4D-Slicer "vv".
if (isfield(mhdpars, 'xcomnt')) externalcomment = mhdpars.xcomnt; else externalcomment = []; mhdpars.xcomnt = externalcomment; end

if (isfield(mhdpars, 'crop')) ask2crop = mhdpars.crop; else ask2crop = 'y'; mhdpars.crop = ask2crop; end %  ask user to crop or not
if (isfield(mhdpars, 'sidsls')) seemidslice = mhdpars.sidsls; else seemidslice = 'n'; mhdpars.sidsls = seemidslice; end

% ask2crop = 'y'; % Check if sagittal, then ask user to crop or not.  Could improve co-reg of sag to axials.
% seemidslice = 'y';
%mhdpars  %debug

sagv1 = [0;1;0]; sagv2 = [0;0;-1]; % Use these to ask to crop sagittal series - potentially improve co-reg to axial.
corv1 = [1;0;0]; corv2 = [0;0;-1]; % Use these to ask to crop coronal series - potentially improve co-reg to axial.
axv1 = [1;0;0]; axv2 = [0;1;0];


if (isempty(imdat)) % no image data; NOTE: will always overwrite "nifti_params" 
%    fileprefix = 'def'; % Will redefine later to default fname convention based on dicom stored in nifti_params.
%    voln = 'c'; % In the event of 4D data, prompt user to 'c'hoose which volume number.
%    browz = 3; % Lately, tend to start search from current location.
%    nvect = 'n'; % 'n' (default) = write out N separate 3D files required by 3D-Slicer; 'y' = write out a single N-vector 4D file readable by 4D-Slicer "vv".
%    externalcomment = []; % (default) = no comment.
    [imdat,dim1,dim2,dim3,dim4,fov1,fov2,fov3,fov4,tr,unique_te,ppd,flip,ii,nifti_params] = readdicom7(0,browz,'fp');    
end % no input arguments or empty "nifit_params"       

if (isempty(nifti_params)) % no nifti-params
    [fnm,path2f] = uigetfile('*.mat', 'Pick a file containing "nifti_params" (e.g., "image_order" or "ExamDemographics")'); 
    if (isequal(fnm, 'image_order.mat')) 
        load([path2f '\' fnm], 'nifti_params');
    elseif (isequal(fnm,'ExamDemographics.mat'))
        load([path2f '\' fnm], 'ExamSeries');
        nifti_params = ExamSeries;
        clear ExamSeries;
    else
        npf = load([path2f '\' fnm]); % arbitrary file that contains nifti_params
        fldn = fieldnames(npf);
        nifti_params = getfield(npf, char(fldn(1)));
        clear npf;
    end % load from file
end 

    % User can provide nifti_params (1,1) structure, or full ExamSeries (1,nser) structure.  
    % Check if input structure has dims > 1.  If yes, assume input "nifti_params" is actually "ExamSeries". 
[s1 nser] = size(nifti_params); % check if multiple series
if (nser > 1) % assume that "ExamSeries" was used
    for iser = 1:nser
        cs{iser,1} = iser;
            %c{iser,2} = nifti_params(1,iser).SeriesNumber;
        if (isfield(nifti_params,'SeriesDescription')) 
            cs{iser,2} = nifti_params(1,iser).SeriesDescription; 
        else
            cs{iser,2} = 'Not a valid ExamSeries structure' ;
        end
    end % for iser
    disp('The full ExamSeries list is: ');
        cs
    iser = input('What is index of desired series? ');
    if (isfield(nifti_params, 'NiftiParams')) 
        np = nifti_params(1,iser).NiftiParams;
        clear nifti_params; % actually
        nifti_params = np; % now re-assign nifti_params to actual nifti_params:
    else
        np = 'Not a valid "ExamSeries" structure';
        disp(['WARNING: ' np '; "nifti_params" are not assigned']);
        nifti_params = [];
    end     
end % if nser    
mhdpars.nipar = nifti_params; % update input param structure

%mhdpars  %debug

% retrive binary data
if (class(imdat) == 'struct')
    raw = getsafield(imdat);
    %[nsl nvol] = size(imdat);
    [dim1 dim2 dim3 dim4] = size(raw);
    nsl = dim3;
    nvol = dim4;
    ipp_arr = (getsafield(imdat(:,1), 'loc'))';    % Need to specify "(:,1)" in case imdat is 4D
    iop_arr = (getsafield(imdat(:,1), 'orient'))'; % Need to specify "(:,1)" in case imdat is 4D
    idz_arr = zeros(nsl,1);
else
    raw = imdat; % Then imdat is only a simple array of data. Get rest from nifti_params
    [dim1 dim2 dim3 dim4] = size(raw);
    nsl = dim3;
    nvol = dim4;
    ipp_arr = (getsafield(nifti_params,'allloc'))';
    iop_arr = (getsafield(nifti_params,'allorientation'))';
    idz_arr = zeros(nsl,1);
end % if class


if ( seemidslice == 'y' )
    if ( isa(voln,'numeric') )
        img(squeeze(raw(:,:,ceil(dim3/2),voln)));
        title([fileprefix ' ' externalcomment]);
    else
        img(squeeze(raw(:,:,ceil(dim3/2),1)));
        title([fileprefix ' ' externalcomment]);
    end % if isa
end

idz_arr = zeros(nsl,1); % Will be filled in later
loc = ipp_arr'; % in 4D case, use locs list from vector plane 1.  It should be same across all vector planes.
pixelspacing = (nifti_params.pixelspacing)';

slicectr2ctr = sqrt( sum( (loc(:,1) - loc(:,end) ).^2) ) / max([1 (dim3-1)]); % SpacingBetweenSlices

if ( (ndims(raw) < 2) || (ndims(raw) > 4) )
    disp('Failed: This is not a 2D, 3D or 4D DICOM data set!'); 
    return;
end % if ndims

% Define base output filename
%if (fileprefix == 'def')
if (strcmp(fileprefix,'def'))
    basefn = ['SID_00_' nifti_params.defname];
else
    basefn = [fileprefix '_' nifti_params.defname];
end

%%% 4th dim labeling ....
vlbls = ''; % init string ov volume labels
% volbl  % NOTE: not used??
% iscell(mhdpars.vlbs)
if (isequal(mhdpars.vlbs{1},'g')) % auto-generate 4D volume labels from "nifti_params"
    disp('(FYI) Auto-generating 4th-dim labels from "nifti_params"...');
    mhdpars.vlbs = mk_ser4Dlbl(mhdpars.nipar);
elseif (isequal(mhdpars.vlbs{1},'c')) % will work for single 3D volume only for now
    disp('');
    disp(' (FYI) Assign (single) volume label, eg, "3DSagT1w" or "phase3" (quotes are optional)');
    disp(' (FYI) OR enter (without quotes) "i" to ignore; "g" to auto-generate from "nifti_params" ');
    s1 = input('     Choose volume label : ', 's'); 
    mhdpars.vlbs{1} = s1;
    if (isequal(s1,'i'))
        mhdpars.vlbs{1} = 'IGNORED';
        vlbls = 'IGNORED';
    elseif (isequal(s1,'g'))
        disp('(FYI) Auto-generating 4th-dim labels from "nifti_params"...');
        mhdpars.vlbs = mk_ser4Dlbl(mhdpars.nipar);
    else
        mhdpars.vlbs{1} = s1;
        vlbls = ['USER INPUT: ' mhdpars.vlbs{1}];
    end   
elseif (isequal(mhdpars.vlbs{1},'i'))
   mhdpars.vlbs{1} = 'IGNORED';
   vlbls = 'IGNORED';
end

lbl1 = mhdpars.vlbs{1};
sl = 1; % start label count
if (isequal(lbl1,'g') || isequal(lbl1(1),'(')) % auto-generated labels
   sl = 2; % 2nd label is 1st phase
end 
lenlbl = length(mhdpars.vlbs);
% check raw-dim4 versus label vector length
if (~isequal(lenlbl-sl+1, dim4)) % somehow don't use for single 3D volume label
    templbl = {};
    for ii1 = 1:dim4
        templbl{ii1} = vlbls;
    end
    mhdpars.vlbs = templbl;
    vlbls = ['(FYI) dim4-label mismatch' mhdpars.vlbs{1}]; % dimention mismatch -- use first label
else
    for i1 = sl:lenlbl % label string
        vlbls = [vlbls mhdpars.vlbs{i1}]; % all 4D phases
    end
end
%vlbls %debug
%mhdpars.vlbs %debug
% Add extra comment lines into mhd for scanner demographics etc into cell array
ncommentlines = 12; %11;
extracomment = cell(1,ncommentlines+2);
extracomment{1} = '% *****************************************************';
extracomment{2} = ['% Institution = ' nifti_params.institution];
extracomment{3} = ['% Manufacturer = ' nifti_params.manufacturer];
extracomment{4} = ['% Model = ' nifti_params.model];
extracomment{5} = ['% FieldStrength = ' num2str(nifti_params.fieldstrength)];
extracomment{6} = ['% StationName = ' nifti_params.stationname];
extracomment{7} = ['% SerialNo = ' nifti_params.serialno];
extracomment{8} = ['% SoftwareVersion = ' nifti_params.swversion];
extracomment{9} = ['% StudyDescription = ' nifti_params.studydesc];
extracomment{10} = ['% SeriesDescription = ' nifti_params.seriesdescription];
extracomment{11} = ['% ImageTypeString = ' nifti_params.imgtypstring];
extracomment{12} = ['% 4th dim vector labels = ' vlbls]; % all 4th dim labels
extracomment{13} = ['% ExternalComment(opt) = ' externalcomment];
extracomment{14} = '% *****************************************************';

% Logic to select single 3D volume, or mutiple 3D volumes from 4D
% As well as filenames these mhd's will have.
if (nvol == 1)
    fnametemp = cell(1,1);
    if (ischar(voln))
        if(isequal({voln},{'c'}))
            fnametemp{1} = [basefn '_1'];
        else
            fnametemp{1} = [basefn '_' voln];
        end % if isequal
    else
        %fnametemp{1} = [basefn '_' num2str(voln)];
        % Then a specific output volume name-string was provided in "voln" for this single volume.
        % Still inside of voln being char string.
        fnametemp = cell(1,1);
        fnametemp{1} = [basefn '_' voln]; % note, "voln" was entered as a string.  
    end % if ischar
    
    % fnametemp{1} = [basefn '_1'];
    id4s = 1; % start volume number
    id4e = 1; % end volume number
    extracomment{12} = ['% 3D volume label = ' mhdpars.vlbs{1}];
else
    % The 4D image was provided
    if (ischar(voln))
        if (isequal({voln}, {'c'}) || isequal({voln}, {'C'})) % Prompt for volume choice
            disp(['      FYI: This data contains ' num2str(nvol) ' 3D volumes.']);
            vol = input('      Create mhd of which volume (1,2,...)? Or enter "A" (no quotes) to do ALL volumes. ','s');
            if ( (vol == 'a') || (vol == 'A') ) % Then do all volumes automatically using default filenames
                fnametemp = cell(1,nvol);
                for ii = 1:nvol
                    fnametemp{ii} = [basefn '_' num2str(ii)];
                end % for ii
                id4s = 1;
                id4e = nvol;
            else % then a specific volume number was requested 
               fnametemp = cell(1,1);
               fnametemp{1} = [basefn '_' vol]; % note, "vol" was entered as a string
               id4s = str2num(vol);
               extracomment{12} = ['% 3D volume label = ' mhdpars.vlbs{id4s + sl - 1}];
               %extracomment{12}
               id4e = id4s;
            end % if vol == 'a'
        elseif (isequal({voln}, {'a'}) || isequal({voln}, {'A'}))
            fnametemp = cell(1,nvol);
            for ii = 1:nvol
                fnametemp{ii} = [basefn '_' num2str(ii)];
            end % for ii
            id4s = 1;
            id4e = nvol;         
        elseif (isequal({voln}, {'s'}) || isequal({voln}, {'S'})) % MD 20150828: same as volume labels
            fnametemp = cell(1,nvol);
            for ii = 1:nvol
                fnsufx = cleanStr(mhdpars.vlbs{ii+sl-1}, specchar); % MD 20150831: clean special characteres from suffix
                if (~isempty(fnsufx))
                    fnametemp{ii} = [basefn '_' fnsufx]; % add volume labels to file name
                else
                    fnametemp{ii} = [basefn '_' num2str(ii)]; % could not retrieve legit suffix from label
                end
            end % for ii
            id4s = 1;
            id4e = nvol;         
        end
%         % Then a specific output volume name-string was provided in "voln" for presumably a single volume.
%         % Still inside of voln being char string.
%         fnametemp = cell(1,1);
%         fnametemp{1} = [basefn '_' voln]; % note, "voln" was entered as a string.
%         id4s = 1;
%         id4e = 1;
    else % Then voln was entered as a number. Use default naming based on that number:
        fnametemp = cell(1,1);
        fnametemp{1} = [basefn '_' num2str(voln)]; % note, "voln" was entered as a number
        extracomment{12} = ['% 3D volume label = ' mhdpars.vlbs{voln+sl-1}];
        id4s = voln;
        id4e = voln;        
    end % if ischar
    
end % if nvol == 1


% retrive patient image position/location and orientation
% ipp_arr = (getsafield(imdat, 'loc'))';
% iop_arr = (getsafield(imdat, 'orient'))';
% idz_arr = repmat(0, size(ipp_arr,1), 1);


% define mhd header structure template.  These are default values to be updated to dicom-derived values later.

% Need to distinguish case of 4D with nvect = 'y' from everything else:
if ( ~((ndims(raw) == 4) && (nvect == 'y')) ) % This is 3D or is 4D with nvect = 'n', so treat as multiple 3D's.
    
    % *****************************************************************
    % *****************************************************************

    mhd = struct('ObjectType', 'Image', ...
                 'NDims', 3,      ...
                 'BinaryData','True', ...
                 'BinaryDataByteOrderMSB', 'False', ...
                 'CompressedData','False', ...
                 'TransformMatrix', reshape(eye(3),1,9), ...
                 'Offset', zeros(1,3), ...
                 'CenterOfRotation', zeros(1,3), ...
                 'AnatomicalOrientation', 'RAI', ...
                 'ElementSpacing', ones(1,3), ...
                 'DimSize', [256,256,64], ...
                 'ElementType', precisn, ...
                 'ElementDataFile', 'temp.raw');
                 % mhd.TransformMatrix = [iop_arr(1,1:3), iop_arr(1,4:6), cross(iop_arr(1,1:3), iop_arr(1,4:6))];
                 
    % fill in the structure
    mhd.TransformMatrix = [iop_arr(1,1:3), iop_arr(1,4:6), cross(iop_arr(1,1:3), iop_arr(1,4:6))];
    %loc = getsafield(imdat(:,1),'loc'); % in 4D case, use locs list from vector plane 1.  It should be same across all vector planes.

    %slicectr2ctr = sqrt( sum( (loc(:,1) - loc(:,end) ).^2) ) / max([1 (dim3-1)]); % SpacingBetweenSlices
    %mhd.ElementSpacing  = [(ii.PixelSpacing)', slicectr2ctr];
    mhd.ElementSpacing  = [pixelspacing, slicectr2ctr];
    %[pixelspacing, slicectr2ctr]

    % -----------------------------------------------
    % sort image according distance along z direction
    % -----------------------------------------------
    for idx=1:dim3 
        idz_arr(idx) = sum(ipp_arr(idx,:).* cross(iop_arr(idx,1:3),iop_arr(idx,4:6)));
    end
    [temp, idx] = sort(idz_arr);
    % in case the sorting is not in accending order, e.g. original SAG scan
    if (idx(2) - idx(1)) ~= 1     
        ipp_arr = ipp_arr(idx, :);
        iop_arr = iop_arr(idx, :);
        idz_arr = idz_arr(idx);
        %raw     = raw(:,:,idx);
        raw     = raw(:,:,idx,:); % See if this works for 4D.
    end
    mhd.Offset = ipp_arr(1, :);

    % rearange data in XY plane
    mhd.DimSize = [dim2, dim1, dim3];

    % raw = permute(raw, [2,1,3]);
    raw = permute(raw, [2,1,3,4]); % Works for 3D and 4D.

    % specify raw data file base name, move up file path one layer of directories
    %[pathstr,name1,ext]  = fileparts(ii.Filename);   % dicom file list
    % pathstr
    % name1

    %[pathstr2,name2,ext] = fileparts(pathstr);     % 
    % pathstr2
    % name2

%     [pathstr3,name3,ext] = fileparts(pathstr2);
%     % pathstr3
%     % name3


    %for id4 = 1:dim4
    id4index = 1;
    for id4 = id4s:id4e
        %mhd.ElementDataFile = [name3,'_',name2,'_',num2str(id4),'.raw'];
        %fnamebase = [name2,'_',num2str(id4)];
        fnamebase = fnametemp{id4index};
        mhd.ElementDataFile = [fnamebase '.raw']; 

        % save mhd and raw files, by default data is 3D
        % Note, fnamebase is char array; raw is 3D array; mhd is a 1x1 struct; extracomment is a 1x13ish cell array.
        if ( isequal({precisn},{'MET_FLOAT'}) )
            raw = single(raw);
        end
        %mhdpars.vlbs %debug
        extracomment{12} = ['% 4th dim label = ' mhdpars.vlbs{id4+sl-1}]; %(MD)20150827: to write volume label
        
        write_mhd_raw2(fnamebase, raw(:,:,:,id4), mhd, extracomment);
        disp(['  FYI MHD/RAW saved: ' fnamebase]);
        id4index = id4index + 1;
    end % for id4

    % *****************************************************************
    % *****************************************************************
%whos
else % This is 4D with nvect = 'y'
    
    % *****************************************************************
    % *****************************************************************
    
    mhd = struct('ObjectType', 'Image', ...
                 'NDims', 4,      ...
                 'BinaryData','True', ...
                 'BinaryDataByteOrderMSB', 'False', ...
                 'CompressedData','False', ...
                 'TransformMatrix', reshape(eye(4),1,16), ...
                 'Offset', zeros(1,4), ...
                 'CenterOfRotation', zeros(1,3), ...
                 'AnatomicalOrientation', 'RAI', ...
                 'ElementSpacing', ones(1,4), ...
                 'DimSize', [256,256,64,5], ...
                 'ElementType', precisn, ...
                 'ElementDataFile', 'temp.raw');
             
     % fill in the structure
     mhd.TransformMatrix = [iop_arr(1,1:3), 0, iop_arr(1,4:6), 0, cross(iop_arr(1,1:3), iop_arr(1,4:6)), 0, 0, 0, 0, 1];
     %loc = getsafield(imdat(:,1),'loc'); % in 4D case, use locs list from vector plane 1.  It should be same across all vector planes.
     
     %slicectr2ctr = sqrt( sum( (loc(:,1) - loc(:,end) ).^2) ) / max([1 (dim3-1)]); % SpacingBetweenSlices
     %mhd.ElementSpacing  = [(ii.PixelSpacing)', slicectr2ctr, 1];
     mhd.ElementSpacing  = [pixelspacing, slicectr2ctr, 1]; % an extra 1 for 4th dimension
     
    % -----------------------------------------------
    % sort image according distance along z direction
    % -----------------------------------------------
    for idx=1:dim3 
        idz_arr(idx) = sum(ipp_arr(idx,:).* cross(iop_arr(idx,1:3),iop_arr(idx,4:6)));
    end
    [temp, idx] = sort(idz_arr);
    % in case the sorting is not in accending order, e.g. original SAG scan
    if (idx(2) - idx(1)) ~= 1     
        ipp_arr = ipp_arr(idx, :);
        iop_arr = iop_arr(idx, :);
        idz_arr = idz_arr(idx);
        %raw     = raw(:,:,idx);
        raw     = raw(:,:,idx,:); % See if this works for 4D.
    end
    mhd.Offset = [ipp_arr(1, :) 0]; % an extra 0 for 4th dimension

    % rearange data in XY plane
    mhd.DimSize = [dim2, dim1, dim3, dim4];

    % raw = permute(raw, [2,1,3]);
    raw = permute(raw, [2,1,3,4]); % Works for 3D and 4D.

    % specify raw data file base name, move up file path one layer of directories
    %[pathstr,name1,ext]  = fileparts(ii.Filename);   % dicom file list
    % pathstr
    % name1

    %[pathstr2,name2,ext] = fileparts(pathstr);     % 
    % pathstr2
    % name2

    %for id4 = 1:dim4
        %mhd.ElementDataFile = [name3,'_',name2,'_',num2str(id4),'.raw'];
        %fnamebase = [name2,'_',num2str(id4)];
        %fnamebase = [name2,'_V'];
        fnamebase = [basefn '_V']; % Special case, use "V" to indicate full 4D vector
        mhd.ElementDataFile = [fnamebase '.raw']; 

        % save mhd and raw files, by default data is 3D
        %write_mhd_raw2( fnamebase, raw(:,:,:,id4), mhd );
        if ( isequal({precisn},{'MET_FLOAT'}) )
            raw = single(raw);
        end
        write_mhd_raw2(fnamebase, raw, mhd, extracomment);
        disp(['  FYI MHD/RAW saved: ' fnamebase]);
    %end % for id4

    % *****************************************************************
    % *****************************************************************
             
end % if ndims & nvect

clear extracomment;
output_args = mhdpars;
save(['PARS2buildMHD_4D' nvect '_' basefn '.mat'], 'mhdpars');
return % END build _mhd

function [ output_args ] = write_mhd_raw2(fnamebase, raw, mhd, extracomment)
% WITE_MHD_RAW: Save data into mhd/raw format for analysis and visulization
% purposes using free-software 3D Slicer (A&V) or 4D Slicer (V).
% 
% The input data's dimension size is limited to 2,3 and 4, and by default it is 3.
% The datatype is double by default due to data type scaling when converting enhanced into simple DICOM.
% The mhd/raw file name is constructed by two upper layer directory's name.
% 
% Written by Yuxi Pang, 20150126 at Univ of Michigan Hospital
%
% 20150611 TLC Add more scaninfo-related material to comment field of MHD header
% 20150625 TLC needed alternative fname convention, esp for 4D.
% 20150723 TLC Addin CenterOfRotation and AnatomicalOrientation.  Also allow single-precision storage.

if ((ndims(raw) < 2) || (ndims(raw) > 4)) 
    disp(['ERROR: Input binary (RAW) data is not volumetric (3D or 4D) ... QUITTING!'])
    return;
end
if ~isstruct(mhd)
    disp(['ERROR: Input header (MHD) is not structue data type ... QUITTING!'])
    return;
end

% [pathstr,name,ext] = fileparts(fname);  
% if ~exist(pathstr,'dir')
%     disp(['ERROR: File directory does not exsit ... QUITTING!'])
%     return;
% end
% 
% % mhd/raw file name defined using two upper level dir names
% [pathstr2,name2,ext] = fileparts(pathstr);      % 
% name3 = [name2,'_',name]; 

% write out *.raw file
%fid=fopen([pathstr, filesep, name3,'.raw'],'w','native');

fid=fopen([fnamebase '.raw'],'w','native');
    switch mhd.ElementType
        case 'MET_DOUBLE'
            fwrite(fid,raw,'double');
        case 'MET_FLOAT' 
            fwrite(fid,raw,'float');
        case 'MET_USHORT'
            fwrite(fid,raw,'uint16');
        otherwise
            fwrite(fid,raw,'double');  
    end
fclose(fid);

[junk ncommentline] = size(extracomment);
% write *.mhd file
%fid=fopen([pathstr, filesep, name3,'.mhd'], 'w','native');
fid=fopen([fnamebase '.mhd'], 'w','native');
    for icom = 1:ncommentline
        fprintf(fid,'%s\n', extracomment{icom});
    end % for icom
    fprintf(fid,'ObjectType = %s\n', mhd.ObjectType);
    fprintf(fid,'NDims = %d\n',mhd.NDims);
    fprintf(fid,'BinaryData = %s\n', mhd.BinaryData);
    fprintf(fid,'BinaryDataByteOrderMSB = %s\n',mhd.BinaryDataByteOrderMSB);
    fprintf(fid,'CompressedData = %s\n',mhd.CompressedData);
    fprintf(fid,'TransformMatrix = %s\n',num2str(mhd.TransformMatrix, '%-12.8f'));
    fprintf(fid,'Offset = %s\n',num2str(mhd.Offset, '%-12.6f'));
    fprintf(fid,'CenterOfRotation = %s\n',num2str(mhd.CenterOfRotation, '%-12.6f'));
    fprintf(fid,'AnatomicalOrientation = %s\n',mhd.AnatomicalOrientation);
    fprintf(fid,'ElementSpacing = %s\n',num2str(mhd.ElementSpacing, '%-12.6f'));
    fprintf(fid,'DimSize = %s\n',num2str(mhd.DimSize, '%-5d'));
    fprintf(fid,'ElementType = %s\n',mhd.ElementType);
    fprintf(fid,'ElementDataFile = %s\n',mhd.ElementDataFile);
fclose(fid);

%end
return % END write_mhd

function mf = getsafield(s,f,mtcheck)
% 'getsafield' is a Matlab10 script to extract field 'f' from structure-array 's'.
% without loops (to save time).
% Currently-designed intent is to extract numerical and character arrays from structs so
% old-style mlab scripts (based on array mathematics) can be used with
% little rewrite.  Main use is to extract image data from structures
% created by readdicom4, thus 'idata' will be default for 'f'.  eg:
% x = readdicom4; % Note, x is a structure.
% xx = getsafield(x,'idata'); Or equivalent default: xx = getsafield(x);
% returns xx as a float array with dimensions like (xres,yres,nslc,ntps).
%
% xx = getsafield(x,'loc'); % Is another example.
%
% TL Chenevert UMICH
% Copyright 2005.  The Regents of the University of Michigan.
% rev 1 TLC     August 16, 2005
% rev 2 TLC     Dec 20, 2007: Add nbyte input arg
%               nbyte = 8 = same as argument not provided, Is to return double precision floating point (matlab default)
%               nbyte = 4, Is to return single precision floating point
%               nbyte = 2, Is to return signed integer
%               You need to make calling-script properly handle requested
%               data type.
%
% TLC UMich Apr30, 2012: For some reason, a few recent PHYSIOLOGS have
% returned a few "empty cells" in the waveforms that crash the script.
% FYI, when one wave (eg. ppu) has an empty cell, they all tend to have
% that same cell index points as empty.  As a workaround, lines **80-96**
% were and "check4empty" flag was added, although the default will be to
% NOT use the workaround in order to maintain consistency with all prior use of this
% script.  In the event, if you get "CAT" errors in cell2mat, try setting
% the "check4empty" to 'y'; rerun it; then return check4empty to 'n'.
% TLC 20140409: Create an optional third input argument mtcheck = 'n' (def), or 'y' filter-out empty cells. 



subsamp = 'n'; % Normally = 'n', but n some cases, may need to subsample due to memory limits

% 20140409 Start ...
if nargin < 3
    check4empty = 'n'; % the default
else
    check4empty = mtcheck;
end

% ************************************************************************
% Unblot next line for manual override:
% check4empty = 'y'; % Default is 'n'.  Switch to 'y' if sporatic empty cells in physiolog discovered.
% 20140409 ... End
% ************************************************************************

if nargin < 2
    f = 'idata'; % Default use is to extract images after readdicom(4+)
end

struc_chk = isstruct(s); dims = size(s);

if (struc_chk == 0 || length(dims) < 2)
        mf = s;
        disp('Sorry, your input is not a structure-array. ');
        return
else
        sca = struct2cell(s);
    %     disp('FYI, this structure has fields ... ');
    %     s
end

flds = fieldnames(s);
fldi = find(strcmp(f,flds) == 1); % required field index

if (isempty(fldi))
    disp('Sorry, required field cannot be found in this structure. ');
    mf = s;
        return
    else
    %     disp('FYI, this structure has fields ... ');
    %     s
end

%if (length(dims) > 3)
 %   scaf = sca(fldi,:,:,:,:); 
if (length(dims) > 2)
    scaf = sca(fldi,:,:,:);
else
    scaf = sca(fldi,:,:);
end;
%whos scaf
%scaf{1,1,1}
    % cell array for required field
    if iscell(scaf{1,1,1}) % for the case of 1D-cell (e.g., "mark" in logvalues)
        mf = squeeze(scaf);
        mf(1:length(mf)) = mf{1:length(mf)};
    else
        
        % ********************************************
        if (check4empty == 'y')
            % For some unknown reason, a few cells may turn-up empty.  Rather
            % than crash the script, notify user and assign to "zero" or the
            % last non-empty cell value.
            if (isempty(scaf{1,1,1}))
                disp('Found empty cell in point 1');
                scaf(1,1,1) = {int16(0)};
            end % 1st point is special since no prior point available
            for iisempty = 2:length(s)
                if (isempty(scaf{1,1,iisempty}))
                    disp(['Found empty cell in point ' num2str(iisempty)]);
                    scaf(1,1,iisempty) = scaf(1,1,(iisempty-1)); 
                end % end if isempty
            end % for iisempty
        end % if check4empty
        % ********************************************
  
        
        mf = cell2mat(scaf); mf=squeeze(mf); % requested "numeric" matrix-field
        dsc = size(scaf{1,1,1});
        if ((dims(1) > 1) && (dsc(2) > 1)) % reshape for 4D-data
            mf = reshape(mf, dsc(1),dsc(2),dims(1),dims(2));
        end;
        if (length(dims) > 2)
            mf = reshape(mf, dims(3)*dims(2),1); % to comply with "getafield" 1D-output
        end;
    end;

clear s sca scaf;
return % END getsfield

function val = tryGetField(s, field, dftVal)
% TLC 20150204.  Externalized from "build_nii.m"  Originally from Xiangrui Li available
% on MatLab User Community website
% 20150716 DM TC.  Make sure returned "val" is not a structure, but is a field.

if isfield(s, field)
    val = s.(field);
    
    % 20150716 Start ...
    if (isstruct(val))
            fn = fieldnames(val);
            vv = getfield(val,char(fn(1)));
            val = vv;
    end
    % ... End 20150716  
    
elseif nargin>2
    val = dftVal;
else
    val = [];
end

return % END tryGetfield

function strout = cleanStr (strin, char2clean)
% 
% clean string from characteres in the "char2clean" cell array
%
% INPUT : strin - input string
%         char2clean - cell array of characters to clean off , eg, {' ', '"', ';', ':'} 
%
% OUTPUT: string without specified characters
%
% MD 2015/08/31: 
% 

kr = [0]; % init character index array

for i1 = 1:length(char2clean)
    kr = [kr strfind(strin, char2clean{i1})];
end

if (length(kr) > 1) % found characters to remove
    kc = setdiff(1:length(strin), kr(2:end));
    if (length(kc)>0)
        strout = strin(kc);
    else
        strout=''; % empty
    end
else
    disp('  (FYI) Input string did not contain special characters: no cleaning performed');
    strout = strin;
end
    

return

% Matlab 2014 String utility functions (in case missing)
function [c, matches] = strsplit1(str, aDelim, varargin)
%STRSPLIT  Split string at delimiter
%   C = STRSPLIT(S) splits the string S at whitespace into the cell array
%   of strings C.
%
%   C = STRSPLIT(S, DELIMITER) splits S at DELIMITER into C. DELIMITER can
%   be a string or a cell array of strings. If DELIMITER is a cell array of
%   strings, STRSPLIT splits S along the elements in DELIMITER, in the
%   order in which they appear in the cell array.
%
%   C = STRSPLIT(S, DELIMITER, PARAM1, VALUE1, ... PARAMN, VALUEN) modifies
%   the way in which S is split at DELIMITER.
%   Valid parameters are:
%     'CollapseDelimiters' - If true (default), consecutive delimiters in S
%       are treated as one. If false, consecutive delimiters are treated as
%       separate delimiters, resulting in empty string '' elements between
%       matched delimiters.
%     'DelimiterType' - DelimiterType can have the following values:
%       'Simple' (default) - Except for escape sequences, STRSPLIT treats
%         DELIMITER as a literal string.
%       'RegularExpression' - STRSPLIT treats DELIMITER as a regular
%         expression.
%       In both cases, DELIMITER can include the following escape
%       sequences:
%           \\   Backslash
%           \0   Null
%           \a   Alarm
%           \b   Backspace
%           \f   Form feed
%           \n   New line
%           \r   Carriage return
%           \t   Horizontal tab
%           \v   Vertical tab
%
%   [C, MATCHES] = STRSPLIT(...) also returns the cell array of strings
%   MATCHES containing the DELIMITERs upon which S was split. Note that
%   MATCHES always contains one fewer element than C.
%
%   Examples:
%
%       str = 'The rain in Spain stays mainly in the plain.';
%
%       % Split on all whitespace.
%       strsplit(str)
%       % {'The', 'rain', 'in', 'Spain', 'stays',
%       %  'mainly', 'in', 'the', 'plain.'}
%
%       % Split on 'ain'.
%       strsplit(str, 'ain')
%       % {'The r', ' in Sp', ' stays m', 'ly in the pl', '.'}
%
%       % Split on ' ' and on 'ain' (treating multiple delimiters as one).
%       strsplit(str, {' ', 'ain'})
%       % ('The', 'r', 'in', 'Sp', 'stays',
%       %  'm', 'ly', 'in', 'the', 'pl', '.'}
%
%       % Split on all whitespace and on 'ain', and treat multiple
%       % delimiters separately.
%       strsplit(str, {'\s', 'ain'}, 'CollapseDelimiters', false, ...
%                     'DelimiterType', 'RegularExpression')
%       % {'The', 'r', '', 'in', 'Sp', '', 'stays',
%       %  'm', 'ly', 'in', 'the', 'pl', '.'}
%
%   See also REGEXP, STRFIND, STRJOIN.

%   Copyright 2012-2014 The MathWorks, Inc.

narginchk(1, Inf);

% Initialize default values.
collapseDelimiters = true;
delimiterType = 'Simple';

% Check input arguments.
if ~ischar(str)
    error(message('MATLAB:strsplit1:InvalidStringType'));
end
if nargin < 2
    delimiterType = 'RegularExpression';
    aDelim = {'\s'};
elseif ischar(aDelim)
    aDelim = {aDelim};
elseif ~iscellstr(aDelim)
    error(message('MATLAB:strsplit1:InvalidDelimiterType'));
end
if nargin > 2
    funcName = mfilename;
    p = inputParser;
    p.FunctionName = funcName;
    p.addParameter('CollapseDelimiters', collapseDelimiters);
    p.addParameter('DelimiterType', delimiterType);
    p.parse(varargin{:});
    collapseDelimiters = verifyScalarLogical(p.Results.CollapseDelimiters, ...
        funcName, 'CollapseDelimiters');
    delimiterType = validatestring(p.Results.DelimiterType, ...
        {'RegularExpression', 'Simple'}, funcName, 'DelimiterType');
end

% Handle DelimiterType.
if strcmp(delimiterType, 'Simple')
    % Handle escape sequences and translate.
    aDelim = strescape(aDelim);
    aDelim = regexptranslate('escape', aDelim);
else
    % Check delimiter for regexp warnings.
    regexp('', aDelim, 'warnings');
end

% Handle multiple delimiters.
aDelim = strjoin(aDelim, '|');

% Handle CollapseDelimiters.
if collapseDelimiters
    aDelim = ['(?:', aDelim, ')+'];
end

% Split.
[c, matches] = regexp(str, aDelim, 'split', 'match');

%end
return % END strsplit1
%--------------------------------------------------------------------------
function tf = verifyScalarLogical(tf, funcName, parameterName)

if isscalar(tf) && isnumeric(tf) && any(tf == [0, 1])
    tf = logical(tf);
else
    validateattributes(tf, {'logical'}, {'scalar'}, funcName, parameterName);
end

%end
return
%--------------------------------------------------------------------------
function escapedStr = strescape(str)
%STRESCAPE  Escape control character sequences in a string.
%   STRESCAPE(STR) converts the escape sequences in a string to the values
%   they represent.
%
%   Example:
%
%       strescape('Hello World\n')
%
%   See also SPRINTF.

%   Copyright 2012 The MathWorks, Inc.

escapeFcn = @escapeChar;                                        %#ok<NASGU>
escapedStr = regexprep(str, '\\(.|$)', '${escapeFcn($1)}');

%end
return  % END strescape
%--------------------------------------------------------------------------
function c = escapeChar(c)
    switch c
    case '0'  % Null.
        c = char(0);
    case 'a'  % Alarm.
        c = char(7);
    case 'b'  % Backspace.
        c = char(8);
    case 'f'  % Form feed.
        c = char(12);
    case 'n'  % New line.
        c = char(10);
    case 'r'  % Carriage return.
        c = char(13);
    case 't'  % Horizontal tab.
        c = char(9);
    case 'v'  % Vertical tab.
        c = char(11);
    case '\'  % Backslash.
    case ''   % Unescaped trailing backslash.
        c = '\';
    otherwise
        warning(message('MATLAB:strescape:InvalidEscapeSequence', c, c));
    end
%end
return % END escapeChar
%--------------------------------------------------------------------------
function joinedStr = strjoin(c, aDelim)
%STRJOIN  Join cell array of strings into single string
%   S = STRJOIN(C) constructs the string S by linking each string within
%   cell array of strings C together with a space.
%
%   S = STRJOIN(C, DELIMITER) constructs S by linking each element of C
%   with the elements of DELIMITER. DELIMITER can be either a string or a
%   cell array of strings having one fewer element than C.
%
%   If DELIMITER is a string, then STRJOIN forms S by inserting DELIMITER
%   between each element of C. DELIMITER can include any of these escape
%   sequences:
%       \\   Backslash
%       \0   Null
%       \a   Alarm
%       \b   Backspace
%       \f   Form feed
%       \n   New line
%       \r   Carriage return
%       \t   Horizontal tab
%       \v   Vertical tab
%
%   If DELIMITER is a cell array of strings, then STRJOIN forms S by
%   interleaving the elements of DELIMITER and C. In this case, all
%   characters in DELIMITER are inserted as literal text, and escape
%   characters are not supported.
%
%   Examples:
%
%       c = {'one', 'two', 'three'};
%
%       % Join with space.
%       strjoin(c)
%       % 'one two three'
%
%       % Join as a comma separated list.
%       strjoin(c, ', ')
%       % 'one, two, three'
%
%       % Join with a cell array of strings DELIMITER.
%       strjoin(c, {' + ', ' = '})
%       % 'one + two = three'
%
%   See also STRCAT, STRSPLIT.

%   Copyright 2012-2014 The MathWorks, Inc.

narginchk(1, 2);

% Check input arguments.
if ~iscellstr(c)
    error(message('MATLAB:strjoin:InvalidCellType'));
end

% Return early when C is empty.
numStrs = numel(c);
if numStrs < 1 && ( nargin < 2 || ischar(aDelim) )
    joinedStr = '';
    return;
end

% Allocate a cell to join into - the first row will be C and the second, D.
joinedCell = cell(2, numStrs);
joinedCell(1, :) = reshape(c, 1, numStrs);
if nargin < 2
    theDelim = {' '};
elseif ischar(aDelim)
    theDelim = {strescape(aDelim)};
elseif iscellstr(aDelim)
    if numel(aDelim) ~= numStrs - 1
        error(message('MATLAB:strjoin:WrongNumberOfDelimiterElements'));
    end
    theDelim = aDelim;
else
    error(message('MATLAB:strjoin:InvalidDelimiterType'));
end

% Join.
joinedCell(2, 1:numStrs-1) = theDelim;
joinedStr = [joinedCell{:}];

%end
return % END strjoin

