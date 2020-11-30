function [ xi, yi ] = computeCoordinates(blkRadii, nPoints, scale)
% The function generates random sampling pairs from an Isotropic Gaussian
% Distribution for the BIGD feature computation.
% Arguments:
% blkRadii- The radii of the block. The block size is computed as 2*blkRadii+1. 
% nPoints- Number of sampling points conisdered for BIGD feature
% computation.
%
% Return:
% xi, yi- Sampling coordinates for block (Dimension: 2 x nPoints). 
%
% Example: 
%       blkRadii = 6;
%       nPoints = 80;
%       [ xi, yi ] = get_sampling_points(blkRadii, nPoints);
blkSize = 2*blkRadii +1;

xi = zeros(3,nPoints);
yi = zeros(3,nPoints);

idx = 1;

while(idx <= nPoints)
    
    mbSize = 2 ^ mod(idx, scale + 1);
    
    %pts1 = round(normrnd(0,sqrt(blkSize*blkSize/25),[2 nPoints]));
    pts1 = unidrnd(blkSize - mbSize + 1, 2, 1) - blkRadii - 1;
    pts1((pts1) > blkRadii - mbSize + 1) = blkRadii - mbSize + 1;
    pts1((pts1) < -blkRadii) = -blkRadii;
    xi(1:2,idx) = pts1;
    xi(3,idx) = mbSize;
    
    pts2 = unidrnd(blkSize - mbSize + 1, 2, 1) - blkRadii - 1;
    %pts2 = round(normrnd(0,sqrt(blkSize*blkSize/25),[2 nPoints]));
    pts2((pts2) > blkRadii - mbSize + 1) = blkRadii - mbSize + 1;
    pts2((pts2) < -blkRadii) = -blkRadii;
    yi(1:2,idx) = pts2;
    yi(3,idx) = mbSize;
    
    if pts1 ~= pts2
        idx = idx + 1;
    end
    
end
