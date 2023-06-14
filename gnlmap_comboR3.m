
function [Lxyz] = gnlmap(BzXYZ, mhdpar)
%
% based on "gradLxyzR" as of 1 Dec 2015
%
% INPUT:
% BzXYZ (required) -- Bz map in 3D cartezian coordinates 
% columns are coil-fields (dBi), rows are derivatives d/drj
% mhdpar (opt) -- MHD param struture for BzXYZ (see, "gbzmap_combo" or
% "mk_mini_nifit" for requirements)
% 
%
% OUTPUT 
% Lxyz - 9 - 3D matrix structure for gradient correction
% corresponding map-by-map MHDs are built if "mhdpar" is provided
%
% NOTE: depends on "mhd_combo" library to build MHDs
%
% MD 20151207: generate missing Gx/Gy field as needed (from cylindric symmetry)
% MD 20151207: replace NaN with zeros for MHDs
% MD 20180223: before making MHDs, convert output from MAT(magnet)-frame to LPS (to align RL/AP) 
% MD 20180306: flip SI-axis to "IS" (for proper corss-term grad. "sign" in L-mat) 
%

if (nargin < 1) disp('ERROR: BzXYZ structure input is required -- existing...'); return; end
mk_mhd = 'y'; 
if (~exist('mhdpar', 'var')) mk_mhd = 'n'; mhdpar = []; end  %
    
bzz = BzXYZ.zz; [d1z, d2z, d3z] = size(bzz); % assuming the same dimensions for BzXYZ
bzx = BzXYZ.zx; i0 = round(d1z/2); j0 = round(d2z/2); k0 = round(d3z/2);
bzy = BzXYZ.zy;

% MD 20151207: generate missing Gx/Gy field as needed (from cylindrical symmetry)
bzy0 = bzy; bzy0((isnan(bzy(:))>0))= 0;
bzx0 = bzx; bzx0((isnan(bzx(:))>0))= 0;
if ((sum(abs(bzy0(:)))==0)*(sum(abs(bzx0(:)))>0))>0
    bzy = permute(bzx, [2 1 3]);
elseif ((sum(abs(bzx0(:)))==0)*(sum(abs(bzy0(:)))>0))>0
    bzx = permute(bzy, [2 1 3]);
end
clear bzx0 bzy0;

[dxy, dxx, dxz]=gradient(bzx); % "horizontal" gradient (y) is 1st direction
[dyy, dyx, dyz]=gradient(bzy);
[dzy, dzx, dzz]=gradient(bzz);
clear bzx bzy bzz;
g0x = dxx(i0,j0,k0); g0y=dyy(i0,j0,k0); g0z = dzz(i0,j0,k0);
dxx = (dxx/g0x); dyy = (dyy/g0y); dzz = (dzz/g0z);
dxy = (dxy/g0y); dxz = (dxz/g0z); dyx = (dyx/g0x); dyz = (dyz/g0z); dzx = (dzx/g0x); dzy = (dzy/g0y);
% invert X,Y-axis
% dxx=dxx(end:-1:1,end:-1:1,:);dxy=dxy(end:-1:1,end:-1:1,:);dxz=dxz(end:-1:1,end:-1:1,:);
% dyx=dyx(end:-1:1,end:-1:1,:);dyy=dyy(end:-1:1,end:-1:1,:);dyz=dyz(end:-1:1,end:-1:1,:);
% dzx=dzx(end:-1:1,end:-1:1,:);dzy=dzy(end:-1:1,end:-1:1,:);dzz=dzz(end:-1:1,end:-1:1,:);

% MD 20180306: added flip-SI to IS for correct xz,yz,zx,zy grad-sign (xx,yy,zz,xy,yx are symmertic around SI)
% NOTE: also fixes 3rd dim for LPS conversion for MHD 
dxx = dxx(:,:,end:-1:1); dxy = dxy(:,:,end:-1:1); dxz = dxz(:,:,end:-1:1); 
dyx = dyx(:,:,end:-1:1); dyy = dyy(:,:,end:-1:1); dyz = dyz(:,:,end:-1:1); 
dzx = dzx(:,:,end:-1:1); dzy = dzy(:,:,end:-1:1); dzz = dzz(:,:,end:-1:1);

Lxyz.xx = dxx;
Lxyz.xy = dyx; %X-derivative, y-field (Bammer Y-field, x-derivative?)
Lxyz.xz = dzx;
Lxyz.yx = dxy;
Lxyz.yy = dyy;
Lxyz.yz = dzy;
Lxyz.zx = dxz;
Lxyz.zy = dyz;
Lxyz.zz = dzz;
timestamp = datestr(now,'yymmddHHMM'); % time-stamp for default output file names
outfile_prefix = 'GNLmaps_Lxyz_';
if (mk_mhd == 'y') % make MHDs
    if (isfield(mhdpar,'nipar'))
        La(:,:,:,1) = dxx; La(:,:,:,2) = dyx; La(:,:,:,3) = dzx;
        La(:,:,:,4) = dxy; La(:,:,:,5) = dyy; La(:,:,:,6) = dzy;
        La(:,:,:,7) = dxz; La(:,:,:,8) = dyz; La(:,:,:,9) = dzz;
        La((isnan(La(:))>0)) = 0; % MD 20151207: replace NaN with "zeros"
        mhdpl = mhdpar;
        np = mhdpar.nipar;
        mhdpl.vlbs = {'Lxx' 'Lxy' 'Lxz' 'Lyx' 'Lyy' 'Lyz' 'Lzx' 'Lzy' 'Lzz'};
        ofipr = ''; 
        if (isfield(mhdpar, 'fipref')) ofipr = mhdpar.fipref; end
        k1 = strfind(ofipr, '_'); 
        if (((~isempty(k1))*(length(ofipr)>k1(1)))>0 ) 
            outfile_prefix = [outfile_prefix ofipr(k1(1):end)]; %use Gfield file prefix
        else
            if isfield(np, 'pixelspacing') 
               mres = ['R' num2str(round(np.pixelspacing(1))) 'x' num2str(round(np.pixelspacing(2))) '_'];
            else
               mres = ['Runk_']; % unknown resolution
            end
            outfile_prefix = [outfile_prefix mres timestamp];  % default file  prefix;
        end        
        mhdpl.fipref = outfile_prefix;
%La = flipdim(permute(La, [2 1 3 4]),3); % MD20180223: convert from MAT(magnet)-frame to LPS
La = permute(La, [2 1 3 4]); % MD20180306: convert from MAT(magnet)-frame to LPS w/o SI flip (done above)
        build_mhd_combo(La, mhdpl);
        clear La dxx dyx dzx dxy dyy dzy dxz dyz dzz;
    else
        disp('WARNING: nonlegit MHD-pars -- unable to generate MHDs ...');
        outfile_prefix = [outfile_prefix 'Runk_' timestamp];  % default file  prefix
    end
else
    outfile_prefix = [outfile_prefix 'Runk_' timestamp];  % default file  prefix
end
Lxyz_name = [outfile_prefix '.mat']; % same res/time-stamp from MHD and mat-output
save(Lxyz_name, 'Lxyz');
disp(['  FYI MAT-File output saved: ', Lxyz_name]);
return


