%%
% Calculates whether fly is stopped
%
function out = fcnStopped( tA, tAIX)

	dX = squeeze(tA(:,:,3));
	dY = squeeze(tA(:,:,4));

	speed = sqrt(dX.^2 + dY.^2);

    out = (speed < .1);
   
