function expList = applyMetricsToData()

	tSL = totalSummaryList();

	% Can apply gates to different genos here...
% Main geno group for paper
%	useGenos = [7,8,9,10,12,13,14,15,18,19,23,24,53];
%	useGenos = [9, 13, 36, 18, 75, 24];


paperGenos = [7,8,9,10,12,13,14,15,18,19,23,24,53,36,75]; % For paper!
% useGenos = [18]; % 83b only

	singleGenos = [7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 23, 24, 53, 70];
	comboGenos = [25, 28, 29, 31, 33, 34, 35, 39, 41, 43,...
				  44, 45, 48, 49, 51, 54, 55, 56, 58, 60,...
				  62, 64, 66, 67, 69, 85, 89, 93];
	useGenos = unique([singleGenos, comboGenos, paperGenos]);
	expList = [];

	% Make list of experiments
	for genoNn = 1:length(useGenos)
		row = tSL{useGenos(genoNn)};
		for n = 2:length(row)
								% GenoN, expN, repN
			expList(end+1,:) = [useGenos(genoNn),row{n},n-1];
		end
	end
	expList


	genoList   = [];
	flyNList   = [];
	powerNList = [];
	scoreList  = [];


	PTgenoList   = [];
	PTflyNList   = [];
	PTpowerNList = [];
	PTtrialNList = [];
	PTscoreList  = [];

	insList = [];
	outsList = [];
	piList = [];
	trialNList = [];
	allTrackIndex = [];

	for expNn = 1:size(expList,1)
		genoN = expList(expNn,1);
		expN  = expList(expNn,2);
		repN  = expList(expNn,3);

		disp([num2str(expNn),'/',num2str(size(expList,1))]);
		
		% MeanIX has [powerN flyN]
		[scores, meanIX, perTrackScores, tAIX, metricLabels] = calcAllMetrics(expN);
		% [scores, meanIX, perTrackScores, tAIX, metricLabels] = calcAllMetricsFigMatch(expN);

		genoList   = cat(1,genoList,  ones(size(scores,1),1)*genoN);
		flyNList   = cat(1,flyNList,  meanIX(:,2) + (repN - 1)*8);
		powerNList = cat(1,powerNList,meanIX(:,1)); 
		scoreList  = cat(1,scoreList,scores);

		PTgenoList   = cat(1,PTgenoList,  ones(size(perTrackScores,1),1)*genoN);
		PTflyNList   = cat(1,PTflyNList,  tAIX(:,3) + (repN - 1)*8);
		PTpowerNList = cat(1,PTpowerNList,tAIX(:,1)); 
		PTtrialNList = cat(1,PTtrialNList,tAIX(:,4));
		PTscoreList  = cat(1,PTscoreList,perTrackScores);
	end

	scoresByFly = scoreList;
	scoresByTrack = PTscoreList;

	save(['~/RTFW/Code/kennedyAnalysis/PTkennedyAnalysisAllGenos.mat'],...
		'scoresByFly','genoList','flyNList','powerNList',...
		'scoresByTrack','PTgenoList','PTflyNList','PTpowerNList','PTtrialNList',...
		'metricLabels');




