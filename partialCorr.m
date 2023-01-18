function c=partialCorr(varargin)

% Compute spatial correlation between two images, controlling for the third
% Input: Paths of three images, e.g. 
%c=spatCorr('/opt/rikandi/Data/DC/Results/2nd/betadifferenceISCinDCposregions.img','/opt/rikandi/Data/DC/Results/2nd/GMposDC.img')
% Images can also be given through SPM UI, in this case just use 
% c=partialCorr();
%
% Output: c - partial spatial correlation between the first two given images

if nargin<1
    input=spm_select(3,'any','Select input images');
    i1=spm_read_vols(spm_vol(input(1,:)));
    i2=spm_read_vols(spm_vol(input(2,:)));
    i3=spm_read_vols(spm_vol(input(3,:)));
end

if nargin==1 || nargin==2 || nargin>3
    disp('Input should be three images')
    return;
end

if nargin==2
    i1=spm_read_vols(spm_vol(varargin{1}));
    i2=spm_read_vols(spm_vol(varargin{2}));
    i3=spm_read_vols(spm_vol(varargin{3}));
end

c=partialcorr([nonzeros(i1(:)) nonzeros(i2(:)) nonzeros(i3(:))]);
c=c(2);