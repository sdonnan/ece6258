function features = computeIGradientDmd(img, dmdOpts)
% The function computes BIGD features for an image img using the options
% provided in the dmdOpts structure.
% Arguments:
% img- The input image for BIGD feature computation. 
% dmdOpts.xi, dmdOpts.yi - The sampling coordinates for BIGD feature
%                           computation.
% dmdOpts.scale: The size of the microblocks considered for computing the
%                BIGD features. The micro-block size is varied from 1 to the
%                scale. The sampling points are equally divided among all
%                scales.
%
% dmdOpts.gridspace: The patches for the feature computation are extracted
%                 after every K pixel along X and Y direction. K represent
%                 the gridspace. 
% Return:
% xi, yi- Sampling coordinates for block (Dimension: 2 x nPoints). 
%
%
% Example: 
%       blkRadii = 6;
%       nPoints = 80;
%       [ xi, yi ] = get_sampling_points(blkRadii, nPoints);
%       dmdOpts.xi = xi;
%       dmdOpts.yi = yi;
%       dmdOpts.gridspace = 2;
%       dmdOpts.scale = 4;
%       img = imread('lena.jpg');
%       dmdFeature = computeDmd(img, dmdOpts);
% initialization
if size(size(img),2) == 3
    img = im2double(rgb2gray(img));
end
pts1 = dmdOpts.xi; 
pts2 = dmdOpts.yi; 
blkRadii = dmdOpts.radii; 
gridSpacing = dmdOpts.gridspace; 
maxScale = dmdOpts.scale;
numSamp = size(pts1,2);
samp_per_scale = numSamp/maxScale;

% Check boundary conidtion for the block sizes
pts1((pts1) > blkRadii - maxScale + 1) = blkRadii - maxScale +1;
pts1((pts1) < -blkRadii) = -blkRadii;
pts2((pts2) > blkRadii - maxScale+1) = blkRadii - maxScale+1;
pts2((pts2) < -blkRadii) = -blkRadii;

% Shift the sampling points 
pts1 = uint16(pts1 + blkRadii +1);
pts2 = uint16(pts2 + blkRadii +1);

% calculate block size and padding
blkSize = 2*blkRadii +1;
npts = size(pts1,2);
dimg = double((img));
[r, c] = size(img);
effr = r-blkSize;
effc = c-blkSize;


% Calculate image gradients and gradient magnitudes
% [Gx, Gy] = imgradientxy(img);
[Gx, Gy] = imgradientxy(dimg);
% Igradient_fd = 5;
Igradient_fd = 3;
featureimg_total = zeros(r,c,Igradient_fd);
% featureimg_total(:,:,1) = img;
featureimg_total(:,:,1) = dimg;
featureimg_total(:,:,2) = sqrt(Gx.^2 + Gy.^2);
featureimg_total(:,:,3) = atan2(Gx,Gy);

v_new = [];

% for each dimmension of the descriptor 
for j = 1:Igradient_fd

    % Assign feature images (I, dx, abs(dx), dy, and abs(dy))
    featureimg = featureimg_total(:,:,j);    
    dfeatureimg = double((featureimg));
    v = [];    
    
    % Compute integral images
    itimg = cumsum(featureimg,1);
    itimg = cumsum(itimg,2);
    iimg = zeros(size(itimg) +2);
    iimg(2:end-1,2:end-1) = itimg;

    % Compute the normalization constant for each block 
    % blur a squared version of the image
    normMat = sqrt( imfilter(dfeatureimg.^2, ones((blkRadii*2) +1) ) );
    % extract the part in image coordinates
    normMat = normMat(blkRadii+1:end-blkRadii,blkRadii+1:end-blkRadii);
    % any normalization constant that is zero set to a small value
    idx_zero = normMat == 0;
    normMat(idx_zero) = 1e-10;

    % Compute BIGD features using different microblock sizes for each point
    % pair as necessary to cover full range of scales
    for i = 1:npts

        % Micro-block size for the current sampling pair
        % This is assigned to distribute sampling pairs across the desired
        % scales. Smallest scale is calculated first.
        mbSize = floor((i + samp_per_scale -1)/samp_per_scale);

        % the next section takes advantace of the integrated images to
        % calculated block sums. not clear why we don't do this straight up
        % without integrating first (ASD)
        
        % Integral image coordinates for computing the sum of the pixel values
        % of size mbSize
        iiPt1 = iimg(pts1(1,i) + mbSize:pts1(1,i) + effr + mbSize, ...
                     pts1(2,i) + mbSize:pts1(2,i) + effc + mbSize);
        
        iiPt2 = iimg(pts1(1,i) + mbSize:pts1(1,i) + effr + mbSize, ...
                     pts1(2,i):pts1(2,i) + effc);
        
        iiPt3 = iimg(pts1(1,i):pts1(1,i) + effr, ...
                     pts1(2,i)+ mbSize:pts1(2,i) + effc + mbSize);
        
        iiPt4 = iimg(pts1(1,i):pts1(1,i) + effr, ...
                     pts1(2,i):pts1(2,i) + effc);

        % The block sum for the mbSize for whole image
        blockSum1 = iiPt4 + iiPt1 - iiPt2 - iiPt3;

        % Integral image coordinates for the blocks using pts2 as reference
        iiPt1 = iimg(pts2(1,i) + mbSize:pts2(1,i) + effr + mbSize, ...
                     pts2(2,i) + mbSize:pts2(2,i) + effc + mbSize);
                     
        iiPt2 = iimg(pts2(1,i) + mbSize:pts2(1,i) + effr + mbSize, ...
                     pts2(2,i):pts2(2,i) + effc);
                     
        iiPt3 = iimg(pts2(1,i):pts2(1,i) + effr, ...
                     pts2(2,i) + mbSize:pts2(2,i) + effc + mbSize);
                     
        iiPt4 = iimg(pts2(1,i):pts2(1,i) + effr, ...
                     pts2(2,i):pts2(2,i) + effc);

        % The block sum for the mbSize for whole image
        blockSum2 = iiPt4 + iiPt1 - iiPt2 - iiPt3;

        % Average intensities of micro-blocks
        blockSum1 = blockSum1/(mbSize*mbSize);
        blockSum2 = blockSum2/(mbSize*mbSize);

        % Block difference
        if i == 3
            diffImg = wrapToPi((blockSum1 - blockSum2)) ./ normMat;
        else
            diffImg = ((blockSum1 - blockSum2)) ./ normMat;
        end

        % Sample the features 
        selectedGrid = diffImg(1:gridSpacing:end,1:gridSpacing:end );
        % v is mb val per block on m for each point/mbsize on n
        v = [v selectedGrid(:)];

    end
    % v_new mb val per block on m for each point/mbsize/featureimg on n
    v_new = [v_new v];
end

% max(max(v_new))

% features= v';
% transpose to put feature source on m and feature values on n
features= v_new';
