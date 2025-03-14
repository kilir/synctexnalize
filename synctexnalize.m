% --------------------------------------------------------------------------
%                         88888 8888 Yb  dP          db    8    w
%  d88b Yb  dP 8d8b. .d8b   8   8www  YbdP  8d8b.   dPYb   8    w 888P .d88b
%  `Yb.  YbdP  8P Y8 8      8   8     dPYb  8P Y8  dPwwYb  8    8  dP  8.dP'
%  Y88P   dP   8   8 `Y8P   8   8888 dP  Yb 8   8 dP    Yb 8888 8 d888 `Y88P
%        dP
% --------------------------------------------------------------------------
%
%% polefigures from 2D detector data obtained from diffraction experiments
% a script to read detector data and produce interpolated pole figures by
% single peak fitting
% - use at your own risk, if you publish non-sense, so do not blame the author
% - assume everything can be wrong unless you know what you are doing
% - feel free to enhance or suggest enhancements
% the script uses functionality provided by the mtex toolbox
% (https://mtex-toolbox.github.io/index)

% TODO:
% - undisort rings using the PONI file in case the detector normal was not
%   parallel to the beam
% - speed up peak fitting option

% written by Rüdiger Kilian (ruediger.kilian@geo.uni-halle.de)

%%
% addpath 'put mtex path here'                             % <- change here
% startup_mtex

%% 0) general set up of input, experiment parameters

% path to data
inpath = '/path/to/image/data/';                     % <- change here
% filename pattern, if just *, all the files in inpath
infiles = dir([inpath 'experiments01_-*.tif']);      % <- change here

% location of masks, backgounds etc. if needed
synchrotronpath = '/path/to/additional/data/'        % <- change here

% rotation increments between images
rotval = 10*degree;

% experimental backgrounds for subtraction if needed
% e.g. 
% bgf = '../air/air_air20_mm.tif'; 
% bgf = [synchrotronpath 'backgrounds_22/water/bkg_SeaWater_0001.tif'];

% mask for dead pixels etc...
bshm_image = [synchrotronpath 'mask_simple.tif'];    % <- change here

% center points of diffraction pattern (specific to experiment geometry)
center = [1032.42 1019.14]; % default

% wavelength
% in meter, e.g. obtained from the poni file
wavelength = 1.4234695572124025e-11;

% detector sample distance in meters
sd = 1.81849;

% pixelsize
pixSize = 0.0002;

rc = linspace(0,2*pi,360);

%% 1) read images
dback = readData(infiles);
%% 2) eventually keep a copy in case reading  is slow
data = im2double(dback);
%% 3) background correction option a) subtract experimental bg
% do you have an appropriate bg?

% function displays the result on the average of the stack
% option 'show' displays before- after comparison on average slice
% 3rd argument is the facor for eventually scaling the bg

%data = bgSubtract(data,bgf,1,'show');
%data = bgSubtract(data,bgf,0.15);


%% 3) background subtraction option b) high pass Fourier filter 
% part1: determine parameters
% works quite ok on broad backgrounds
% what do you need to do? -> find good parameters
% parameters:
% data
% lower limit
% upper limit
% optional: size of gaussian on fft mask

% test on the average of the stack
% a_slice = mean(data,3);
% lp = FTpass(a_slice,11,3600,5,'show');

%%
% part2: run the filter on all slices
%data = FTpass(data,11,3600,5);

%% 4) filter detector stripes (e.g. applies to one of the ESRF detectors)
%data = filterVerticalStripes(data,'show');
% data = filterVerticalStripes(data);

%% 5) apply masks for beamstop etc.
%data = maskData(data,bshm_image,'show');
data = maskData(data,bshm_image);

%% 6) select rings
%% option 1 part1: plot and ...

% uncomment if needed

%{
figure
a_slice = mean(data,3);
ll = quantile(a_slice(:),0.001); 
la = quantile(a_slice(:),0.999);

imagesc(a_slice,[ll la])
colormap parula
axis equal
hold on
plot(center(1),center(2),'rx')


% maybe zoom now if you want to do click
%% option 1 part 2: select rings by clicking lower / upper limit

% select rim of interest (max / min)
dvalue = ginput(2);

%compute distances rom center
r_min_max = sqrt((center(1) - dvalue(:,1)).^2 + (center(2) - dvalue(:,2)).^2);
hold on
% check
plot(center(1) +(cos(rc)*max(r_min_max)), center(2) +(sin(rc)*max(r_min_max)),'k','linewidth',1)
plot(center(1) +(cos(rc)*min(r_min_max)), center(2) +(sin(rc)*min(r_min_max)),'k','linewidth',1)
hold off

% ring in meters (2 theta)
rim = mean(r_min_max) * pixSize;
% 2theta
twothetarim = atand(rim/sd);
% wavelength
dvalue = wavelength / (2*sind(twothetarim/2));
%}

%% %% option 2: supply dvalue and range
%
% % either read CIF or define yourself
%cs = loadCIF([CIFpath 'Kaolinite.cif'])
% or
cs = crystalSymmetry('1',[5.2, 8.9, 7.4],[91.7, 104.862, 89.822]*degree,'X||a*', 'Z||c');
% Miller of interest
h = Miller(0,0,1,cs);

dvalue = h.dspacing * 1e-10;
twothetarim = 2 * asind(wavelength / dvalue / 2);
rim = tand(twothetarim)*sd;
crim = rim/pixSize;

% determine how wide a peak should be
r_min_max = [crim+5 crim-5];       % <- this value needs to be set for automation

% eventually plot
% hold on
% % check
% plot(center(1) +(cos(rc)*max(crim)), center(2) +(sin(rc)*max(crim)),'--r','linewidth',1)
% plot(center(1) +(cos(rc)*max(r_min_max)), center(2) +(sin(rc)*max(r_min_max)),'k','linewidth',1)
% plot(center(1) +(cos(rc)*min(r_min_max)), center(2) +(sin(rc)*min(r_min_max)),'k','linewidth',1)
% hold off

%% 7) read data from theta-eta maps
% here you can also specify the bin width (eta-range over which a constant theta is assumed)
% for the integration of peak areas two methods are available: 'trapz'
% (trapezoidal integration; fast) or 'peakfit' (interation of a fitted Lorentz function; slower but usually better)
% have a look at the function below for the switch between different
% detector setups
[pfvals, g] = rings2VecData(data,center,r_min_max,twothetarim,rotval,'pfmethod','peakfit','binSpec',[0:3:360]*degree);

% [pfvals, g] = rings2VecData(data,center,r_min_max,twothetarim,rotval,'pfmethod','peakfit','binSpec',60);
% [pfvals, g] = rings2VecData(data,center,r_min_max,twothetarim,rotval,'pfmethod','trapz','binSpec',[0:3:360]*degree);
% [pfvals, g] = rings2VecData(data,center,r_min_max,twothetarim,rotval,'pfmethod','trapz');

%% 8) prepare polefigure
[pfv, vu] = prepPF(pfvals,g);

% in case there are single spikes/outliers,'removeOutliers' sigma can be
% used for filtering, sigma defaults to 0.8
% [pfv, vu] = prepPF(pfvals,g,'removeOutliers',0.5)

%% 9) plot results and interpolate
% select a proper bandwidth, potentially as large as possible without
% producing artefacts

% pole figure on evaluated grid
figure
plot(vu,pfv,'MarkerSize',3,'minmax','antipodal')
mtexTitle('rawdata')
nextAxis;

% direct hamornic approximation on the grid
vu.antipodal = 1;
sfd1   = interp(vu,pfv,'harmonic','regularize',10^-8,'SobolevIndex',2.8); 
plot(sfd1,'minmax')
mtexTitle('interpolation')
nextAxis;

% density distribution
sfd2 = sfd1-(min(sfd1));
sfds = sfd2 * (4*pi/sum(sfd2));
plot(sfds,'minmax')
mtexTitle('density')



%% functions
%_-----------------------------------------------------------------------
%% reading the data
function dback = readData(infiles)
% number of files
nfiles = length(infiles);
filelist= 1:nfiles; % 1:2:nfiles % maybe one would only want to take every x-th image
for i=1:length(filelist) % do not read the last image from 0==180
    dback(:,:,i) = double(imread([infiles(i).folder '/' infiles(i).name]));
end
end

%_-----------------------------------------------------------------------
%% plotting and comparing data
function tplot(data,n)
subplot(1,2,n)
a_slice = mean(data,3);
%a_slice = a_slice - min(a_slice(:));
ll = quantile(a_slice(:),0.001); % adapt values
la = quantile(a_slice(:),0.999); % adapt values
imagesc(a_slice,[ll la])
axis equal
end

%_-----------------------------------------------------------------------
%% subtracting the background
function data = bgSubtract(data,bgf,factor,varargin)

if check_option(varargin,'show'), tplot(data,1); end

backg = double(imread(bgf));
data = data - backg*factor;

if check_option(varargin,'show'), tplot(data,2); end

end

%_-----------------------------------------------------------------------
%% Fourier Filtering
function outimage = FTpass(B,ll,ul,varargin)
% data, lower limit, upper limit, (size of gaussian on fft mask)
if check_option(varargin,'show'), tplot(B,1); end

if ndims(B) == 2

    %1) replace nan
    B(isnan(B))=0;

    %2) add padding;
    pwidth = 400;
    B = padArray(B,pwidth);

    %3) do fft
    FT_img = fft2(double(B));

    %4) make the filter
    % get domain size
    [M, N] = size(B);
    % sawtooth pos/neg
    u = 0:(M-1);
    idx = find(u>M/2);
    u(idx) = u(idx)-M;
    v = 0:(N-1);
    idy = find(v>N/2);
    v(idy) = v(idy)-N;

    % grid of u,v
    [V, U] = meshgrid(v, u);

    % distance from 0 frequency
    D = sqrt(U.^2+V.^2);

    % determien filter unsing cut-off frequencies
    H = double(D >= ll & D <= ul);
    if ~isempty(varargin) & isnumeric(varargin{1})
        H = imgaussfilt(H,varargin{1});
    end
    
    %5) convolve fft image and the filter
    G = H.*FT_img;

    %6) do the inverse fft
    outimage = real(ifft2(double(G)));

    %7) remove padding
    outimage = outimage(pwidth+1:end-pwidth-1, pwidth+1:end-pwidth-1);

else
    for i= 1:size(B,3)
        outimage(:,:,i) = FTpass(B(:,:,i),ll,ul,varargin{1});
    end
end

if check_option(varargin,'show'), tplot(outimage,2); end

    function b = padArray(a,w)
        b = [a(:,w:-1:1)  a   a(:,end:-1:end-w)];
        b = [b(w:-1:1,:); b ; b(end:-1:end-w,:)];
    end

end

%_-----------------------------------------------------------------------
%% filter vertical stripes
function imgs = filterVerticalStripes(imgs,varargin)

k = ...
    [0 0 0 0 0;...
     0 0 0 0 0;...
     1 2 4 2 1;...
     0 0 0 0 0;...
     0 0 0 0 0];
kn = k/sum(k(:));

if check_option(varargin,'show'), tplot(imgs,1); end

for i = 1:size(imgs,3)
    d_temp = imgs(:,:,i);
    d_temp = conv2(d_temp,kn,'same');
    imgs(:,:,i) = d_temp;
end

if check_option(varargin,'show'), tplot(imgs,2); end

end

%_-----------------------------------------------------------------------
%% mask with mask
function imgs = maskData(imgs,mask,varargin)

if check_option(varargin,'show'), tplot(imgs,1); end

bshmi = logical(imread(mask));
% and set to nan
for i = 1:size(imgs,3)
    d_temp = imgs(:,:,i);
    d_temp(bshmi) = nan;
    imgs(:,:,i) = d_temp;
end
if check_option(varargin,'show'), tplot(imgs,2); end
end

%_-----------------------------------------------------------------------
%% extract the data
function [pfvals, g] = rings2VecData(data,center,r_min_max,twothetarim,rotval,varargin)

% prepare
a_slice = mean(data,3);
% get coordinates of eta - theta maps
[pix_X,pix_Y] = meshgrid(1:size(a_slice,1),1:size(a_slice,2));
% shift grid values relative to center
pix_X = center(1) - pix_X;% - center(1);
pix_Y = center(2) - pix_Y;% - center(2);

pix_eta = atan2(pix_Y,pix_X)+pi; % radial distance to select ring
pix_theta = sqrt(pix_Y.^2 + pix_X.^2); % circlular map to get data

% %visualize
% subplot(1,2,1)
% imagesc(pix_eta/degree); hold on
% plot(center(1),center(2),'rx'); % hold off
% axis equal; title('eta'); xlim([0 size(a_slice,1)]); ylim([0 size(a_slice,2)])
%
% subplot(1,2,2)
% imagesc(pix_theta); hold on
% plot(center(1),center(2),'rx'); hold off
% axis equal; title('theta'); xlim([0 size(a_slice,1)]); ylim([0 size(a_slice,2)])

 
% define read-out eta-theta within rim

% those are the ids of pixels within the ring
id_theta = pix_theta < max(r_min_max) & pix_theta > min(r_min_max);

% eta_theta in the ring
eta_value = pix_eta(id_theta);
theta_value = pix_theta(id_theta);


% get bins and ids eta bins such that we can have user-defined eta-intervals
binspec = get_option(varargin,'binSpec',60);

if length(binspec) == 1 
    eta_bins = linspace(0,2*pi,binspec);
elseif length(binspec) > 1 & max(binspec) ==  2*pi & min(binspec) ==  0;
    eta_bins = binspec;
else
    disp('Please make sure your vector of bins goes from 0 to 2pi')
end
[~,~,etabinId] = histcounts(eta_value,eta_bins);


% read the data to vector3d and values
clear pfinc pvals pfazi g
% data comes in x degree steps, so we loop for each image
for i= 1:size(data,3)
    %read out data for each slice
    tempdata = data(:,:,i);
    % reduce to those in the rim
    tempdata = tempdata(id_theta);

    % centers of eta bins which will be plotted
    pfazi{i} = (eta_bins(1:end-1) + (eta_bins(2) - eta_bins(1))/2)';

    % since we need to sort values for theta for proper integration
    % we need this loop
    intG = zeros(1,length(pfazi{i}));
    for k = 1:length(pfazi{i}) % loop over every binned/rounded eta

        % what goes into the k'th bin data
        val_k= tempdata(etabinId == k);

        % theta values
        theta_k = theta_value(etabinId == k);

        % sorting for theta inside bin
        [theta_k, idtheta_s] = sort(theta_k);
        val_k = val_k(idtheta_s);

        % in case something is < 0 and not masked
        theta_k(val_k<0) = []; val_k(val_k<0) = [];

        % testing
        % plot(thetatest_s,valstest_s,'-x') %

        % decide which method to use
        pfmethod = get_option(varargin,'pfmethod','trapz');
        switch pfmethod
            case 'peakfit'
                % slow and fine in case of very few points in each interval
                intG(k) = lorentzIntegral(theta_k,val_k,r_min_max);

            case 'trapz'
                % fast and fine in case of many points in each interval
                intG(k) = trapz(theta_k,val_k);
        end

    end
    % write data into output
    pfvals{i} = intG';

    % since the entire cone is 4 theta, we need to have 2* theta 
    pfinc{i} = repmat(pi/2 - twothetarim*degree,length( pfvals{i}),1);

    % each new scan is rotated with respect to the first one; the first
    % rotation is the zero rotation
    rotinc = rotation.byAxisAngle(xvector,(rotval * (i-1))); % setup ESRF
    

    %define grid
    g{i} = (rotinc * vector3d.byPolar(pfinc{i},pfazi{i}))';

end

    function li = lorentzIntegral(theta,data,limits)
        % fit a lorentz function to theta- data values and get the
        % integral, actually this gives pretty much the
        % same as trapz as long as the points are fine, for certain
        % eta/theta ranges, it might however be beneficial

        % also lsqcurvfir does not like nan
        theta(isnan(data)) = [];
        data(isnan(data)) = [];

        % initial parameter guesses
        maxt = max(theta(:));
        mint = min(theta(:));
        p3 = ((maxt-mint)./10).^2;
        p2 = (maxt+mint)./2;
        p1 = max(data(:)).*p3;
        p0 = [p1 p2 p3];
        % compute coefficients
        p = lsqcurvefit(@lorentz,p0,theta,data,[],[],optimset('Display','off'));

        % do the integral
        li = integral(@(x) (p(1))./((x-p(2)).^2+p(3)),min(r_min_max),max(r_min_max));
    end

    %define lorentz
    function F = lorentz(p,x)
        F = p(1)./((x-p(2)).^2+p(3));
    end


end

%_-----------------------------------------------------------------------
%% prep polefigure

function [pfv, vu] = prepPF(pfvals,g,varargin)
v = [g{:}];
% just in case two measurement points exactly overlap
% we will average them
[vu,~,ir]=unique(v,'stable');
pfv = cell2mat(pfvals(:));
pfv = accumarray(ir,pfv,[],@nanmean);

% filter nan values (just in case)
idnan = isnan(pfv);
vu(idnan) = [];
pfv(idnan) = [];

if check_option(varargin,'removeOutliers')
    % remove outliers (slow) but maybe helpful
    % here numeric value is fraction of std
    olsigma = get_option(varargin,'removeOutliers',0.8);
    %%
    ido = isOutlier(vu,pfv,olsigma);
    %%
    vu(ido) = [];
    pfv(ido) = [];
end
function ind = isOutlier(v,val,facstd)
% lifted from Mtex Polefigure/isOutlier
        % find neighbours
        next = find(v,v,3*v.resolution);
        %remove diagonal
        next(speye(length(next))==1) = false;
        % compute mean
        dmean = next * val ./ sum(next,2);
        dstd = std(val);
        ind = logical(abs(dmean - val)>facstd*dstd);
 end

end
%_-----------------------------------------------------------------------
