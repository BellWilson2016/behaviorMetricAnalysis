%%
%
%   Head speed in any direction. (Abs value!)
%
function out = fcnSpeed( tA, tAIX)

    samplePeriod = .05;

    dX = squeeze(tA(:,:,3));
    dY = squeeze(tA(:,:,4));
       
    out = sqrt(dX.^2 + dX.^2);
