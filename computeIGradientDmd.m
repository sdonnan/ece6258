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
%maxScale = 2 ^ dmdOpts.scale;

blkSize = 2*blkRadii +1;
npts = size(pts1,2);
dimg = double((img));
[r, c] = size(img);
effr = r-blkSize;effc = c-blkSize;


% Calculate image gradients and gradient magnitudes
% [Gx, Gy] = imgradientxy(img);
[Gx, Gy] = imgradientxy(dimg);
Igradient_fd = 5;
featureimg_total = zeros(r,c,Igradient_fd);
% featureimg_total(:,:,1) = img;
featureimg_total(:,:,1) = dimg;
featureimg_total(:,:,2) = Gx;
featureimg_total(:,:,3) = abs(Gx);
featureimg_total(:,:,4) = Gy;
featureimg_total(:,:,5) = abs(Gy);


v_new = [];

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
    normMat = sqrt(imfilter(dfeatureimg.^2,ones((blkRadii*2) +1)));
    normMat = normMat(blkRadii+1:end-blkRadii,blkRadii+1:end-blkRadii);
    idx_zero = find(normMat == 0);
    normMat(idx_zero) = 1e-10;

    % Compute BIGD features
    for i = 1:npts

        % Micro-block size for the current sampling pair
        %mbSize = 2 ^ mod(i, dmdOpts.scale + 1);
        mbSize = pts1(3, i);
        
        xy1 = pts1(1:2, i);
        xy2 = pts2(1:2, i);
        
        % Check boundary conidtion for the block sizes
%         xy1((xy1) > blkRadii - mbSize + 1) = blkRadii - mbSize + 1;
%         xy1((xy1) < -blkRadii) = -blkRadii;
%         xy2((xy2) > blkRadii - mbSize + 1) = blkRadii - mbSize + 1;
%         xy2((xy2) < -blkRadii) = -blkRadii;

        % Shift the sampling points 
        xy1 = uint16(xy1 + blkRadii +1);
        xy2 = uint16(xy2 + blkRadii +1);

        % Integral image coordinates for computing the sum of the pixel values
        % of size mbSize
        iiPt1 = iimg(xy1(1) + mbSize:xy1(1)+effr + mbSize,xy1(2) + mbSize:xy1(2)+effc + mbSize);
        iiPt2 = iimg(xy1(1)+ mbSize:xy1(1)+effr+ mbSize,xy1(2):xy1(2)+effc);
        iiPt3 = iimg(xy1(1):xy1(1)+effr,xy1(2)+ mbSize:xy1(2)+effc+ mbSize);
        iiPt4 = iimg(xy1(1):xy1(1)+effr,xy1(2):xy1(2)+effc);

        % The block sum for the mbSize for whole image
        blockSum1 = iiPt4 + iiPt1 - iiPt2 - iiPt3;


        % Integral image coordinates for the blocks using xy2 as reference
        iiPt1 = iimg(xy2(1) + mbSize:xy2(1)+effr + mbSize,xy2(2) + mbSize:xy2(2)+effc + mbSize);
        iiPt2 = iimg(xy2(1) + mbSize:xy2(1)+effr+ mbSize,xy2(2):xy2(2)+effc);
        iiPt3 = iimg(xy2(1):xy2(1)+effr,xy2(2)+ mbSize:xy2(2)+effc+ mbSize);
        iiPt4 = iimg(xy2(1):xy2(1)+effr,xy2(2):xy2(2)+effc);

        % The block sum for the mbSize for whole image
        blockSum2 = iiPt4 + iiPt1 - iiPt2 - iiPt3;


        % Average intensities of micro-blocks
        blockSum1 = blockSum1/(mbSize*mbSize);
        blockSum2 = blockSum2/(mbSize*mbSize);

        % Block difference
        diffImg = ((blockSum1 - blockSum2)) ./ normMat;

        % Sample the features 
        selectedGrid = diffImg(1:gridSpacing:end,1:gridSpacing:end );
        v = [v selectedGrid(:)];

    end
    
        v_new = [v_new v];
end

% max(max(v_new))

% features= v';
features= v_new';
