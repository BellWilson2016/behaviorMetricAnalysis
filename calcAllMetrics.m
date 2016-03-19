function [scores, meanIX, perTrackScores, tAIX, metricLabels] = calcAllMetrics(expN)

	load(['/groups/wilson/derived/filtered/f',num2str(expN,'%04d'),'.mat']);
	
	tA = procTracks;
	tAIX = trackIndex;

    scores = [];
	perTrackScores = [];
    metricLabels = {};
    
    % PI
    gateWindow = [60 90];
    [gIX, stSamp, enSamp] = gateOnTime(tAIX, gateWindow);
    calcFcn = @fcnOccupancy;
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'PI';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));

    
    % decPI
    gatePosition = [-5 5];
    gateWindow = [60 90];
    firstOnly = false; timeEntries = true; timeExits = false;
    [gIX, stSamp, exSamp] = gateOnXPosition(tA, tAIX, gateWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnOccupancy;
    [perGateMeas, measIX] = calcFcnOnGate( tA, tAIX, gIX, exSamp, exSamp + 1, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'decPI';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));

    % Speed in lasered zone
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = true;
    timeExits = true;
    gatePosition = [0 30];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnSpeed;
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Speed in laser';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));
    
    % Speed in non-lasered zone
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = true;
    timeExits = true;
    gatePosition = [-30 0];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnSpeed;
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Speed in non-laser';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));

    % Pstopped in lasered zone
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = true;
    timeExits = true;
    gatePosition = [0 30];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnStopped;
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Pstop in laser';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));
    
    % Pstopped in non-lasered zone
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = true;
    timeExits = true;
    gatePosition = [-30 0];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnStopped;
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Pstop in non-laser';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));
	
    
    % Heading on Laser Entry + [.5 -> 1 sec]
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = false;
    timeExits = true;
    gatePosition = [0 30];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnCosHeadAngle;
	stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Heading after laser entry';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));


    % Heading on Laser Exit + [.5 -> 1 sec]
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = false;
    timeExits = true;
    gatePosition = [-30 0];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnCosHeadAngle;
	stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Heading after laser exit';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));
	

    % Crab speed on Laser Entry + [.5 -> 1 sec]
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = false;
    timeExits = true;
    gatePosition = [0 30];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnCrabSpeed;
	stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Crab-speed after laser entry';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));

    % Crab speed on Laser Exit + [.5 -> 1 sec]
    timeWindow = [60 90];
    firstOnly = false;
    timeEntries = false;
    timeExits = true;
    gatePosition = [-30 0];
    [gIX, stSamp, enSamp] = gateOnXPosition(tA, tAIX, timeWindow, gatePosition, firstOnly, timeEntries, timeExits);
    calcFcn = @fcnCrabSpeed;
	stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Crab-speed after laser exit';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));

    
    % Head angle on upwind laser flash on + [.5 -> 1sec]
    timeStart = 60;
    gatePosition = [5 30]; upLaser = true;
    [gIX, stSamp] = gateFlashEncounterByDirection(tA, tAIX, timeStart, gatePosition, upLaser);
    calcFcn = @fcnCosHeadAngle;
    stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Cos(heading) on upwind flash-on';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));
   
    % Head angle on downwind laser flash on + [.5 -> 1sec]
    timeStart = 60;
    gatePosition = [5 30]; upLaser = false;
    [gIX, stSamp] = gateFlashEncounterByDirection(tA, tAIX, timeStart, gatePosition, upLaser);
    calcFcn = @fcnCosHeadAngle;
    stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Cos(heading) on downwind flash-on';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));

    % Head angle on upwind laser flash off + [.5 -> 1sec]
    timeStart = 90;
    gatePosition = [5 30]; upLaser = true;
    [gIX, stSamp] = gateFlashEncounterByDirection(tA, tAIX, timeStart, gatePosition, upLaser);
    calcFcn = @fcnCosHeadAngle;
    stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Cos(heading) on upwind flash-off';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));
   
    % Head angle on downwind laser flash off + [.5 -> 1sec]
    timeStart = 90;
    gatePosition = [5 30]; upLaser = false;
    [gIX, stSamp] = gateFlashEncounterByDirection(tA, tAIX, timeStart, gatePosition, upLaser);
    calcFcn = @fcnCosHeadAngle;
    stSamp = stSamp + round(.5/.05);
    enSamp = stSamp + round(.5/.05);
    [perGateMeas, measIX] = calcFcnOnGate(tA, tAIX, gIX, stSamp, enSamp, calcFcn);
    [meanMetrics, meanIX] = stratifyByFly(perGateMeas, measIX, tAIX);
    scores = cat(2,scores,meanMetrics(:));
    metricLabels{end+1} = 'Cos(heading) on downwind flash-off';
    [perTrackMetrics, perTrackIX] = stratifyByTrack(perGateMeas, measIX, tAIX);
	perTrackScores = cat(2,perTrackScores,perTrackMetrics(:));


