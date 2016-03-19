function out = fcnOccupancy(tA, tAIX)

    samplePeriod = .05;

    % headPosX = bodyX + headX;
    headPosX = squeeze(tA(:,:,1));
    out = sign(headPosX);
