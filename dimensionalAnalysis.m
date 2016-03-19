function dimensionalAnalysis()


% Kennedy-style descriptors
	load('~/RTFW/Code/kennedyAnalysis/PTkennedyAnalysis.mat');
%	load('~/RTFW/Code/kennedyAnalysis/PTkennedyAnalysisAllGenos.mat');
	meanIX = [powerNList(:),flyNList(:),genoList(:)];

%% Remove Orco - #18
removeOrco = false;
if removeOrco
	ix = find(genoList == 18);
	meanIX(ix,:) = [];
	powerNList(ix) = [];
	flyNList(ix) = [];
	genoList(ix) = [];
	scoresByFly(ix,:) = [];
end

	% Take PI and decPI out
	PI = scoresByFly(:,1);
	decPI = scoresByFly(:,2);
	scoresByFly = scoresByFly(:,3:end);
	for n = 1:size(scoresByFly,2)
		mL{n} = metricLabels{n+2};
	end
	metricLabels = mL;

%%% Machine learning descriptors
%	load('/groups/wilson/gatePop/gatePopResults-totalPop150616.mat');
%	oldMeanIX = meanIX;
%	meanIX(:,1) = oldMeanIX(:,3);
%	meanIX(:,3) = oldMeanIX(:,1);
%	powerNList = meanIX(:,1);
%	genoList = meanIX(:,3);
%	for metricN = 1:size(scoresByFly,2)
%		metricLabels{metricN} = ['m',num2str(metricN,'%03d')];
%	end
%


%% Make a color list
uniqueGenos = unique(genoList);
for genoNn = 1:length(uniqueGenos)
	genoN = uniqueGenos(genoNn);
	ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == 8));
	sortVals(genoNn) = nanmean(decPI(ix));
end
[B, genoSortOrder] = sort(sortVals,'descend');
colorSpace = flipud(jet);
colorIXList = round(linspace(1,64,length(uniqueGenos)));


	% Baseline subtract each genotype
	scoresByFly = baselineSubtractScores(scoresByFly, powerNList, genoList);

	plotTitle = 'Fit Metrics';
	scoresByFly = nanZscore(scoresByFly);
	ix = find(isnan(scoresByFly)); scoresByFly(ix) = 0;
	classIX = uniqueTypes(meanIX);

	% Do the signal correlation
	[sigCorr, noiseCorr] = signalCorrelation(scoresByFly, classIX);

	% Resort by the signal auto-correlation
	[B, IX] = sort(diag(sigCorr),'descend');
	% IX = IX(1:nBest);
	sigCorr = sigCorr(IX,IX);
	noiseCorr = noiseCorr(IX,IX);
	scoresByFly = scoresByFly(:,IX);

	nBoots = 10;
	allPVals = bootSignalCorrelation( scoresByFly, classIX, nBoots);
	diagPVals = diag(allPVals);
	corrPVals = bonferroniHolm(diagPVals);
	reliability = diag(sigCorr);
	disp('Reliability | Raw P | Corrected P | Metric Name');
	disp('---------------------------------');
	for n=1:size(scoresByFly,2)
		disp([num2str(reliability(n)),'  ',...
			  num2str(diagPVals(n)),'  ',...
			  num2str(corrPVals(n)),'  ',...
			  metricLabels{IX(n)}]);
	end

	nToShow = min([100 size(scoresByFly,2)]);
	for n = 1:nToShow
		disp([num2str(n,'%02d'),'  ',num2str(sigCorr(n,n)),'   ',num2str(IX(n),'%02d'),'  ',metricLabels{IX(n)}]);
	end
	gateRawOrder = IX;

	ffsubplot(2,3,1);
	image(sigCorr,'CDataMapping','scaled');
	title('Sig. Corr.'); colorbar;
	axis square;

	ffsubplot(2,3,2);
	image(noiseCorr,'CDataMapping','scaled');
	title('Noise Corr.'); colorbar;
	axis square;

	FF = ffsubplot(2,3,3);
	nPCs = min([64 size(scoresByFly,2)-1]);
	[coeffs,scores,latent,mu,expVar] = signalPCA(scoresByFly,classIX,nPCs);
	plot(latent,'b.-'); hold on;
	xlim([0 10]);


	nPerms = 1000;
	nFlies = size(scoresByFly,1);
	allLatents = zeros(nPerms,nPCs);
	permuteScreePlot = false;
if permuteScreePlot
	for permN = 1:nPerms
		permN
		% Shuffle scores
		shuffledScores = permuteScores(scoresByFly, classIX);
		[scoeffs,sscores,slatent,smu,sexpVar] = signalPCA(shuffledScores,classIX,nPCs);
		allLatents(permN,:) = slatent(1:nPCs);
	end
	for PCn = 1:nPCs
		lIX = round(.05*nPerms);
		uIX = round(.95*nPerms);
		sortedL = sort(allLatents(:,PCn),'ascend');
		lBound(PCn) = sortedL(lIX);
		uBound(PCn) = sortedL(uIX);
	end

	% plot(mean(allLatents,1),'k');
	% plot(lBound,'k--');
	plot(uBound,'k--');
else
	disp('Skipping Scree Plot Permutations');
end
	title(['First 2 PCs: ',num2str(sum(latent(1:2)))]);
	xlim([.5 10.5]);
	xlabel('PC #');
	axis square;
	ylim([0 .65]);

	cumsum(latent(1:nPCs))

% sendKey('','r');
% FF.PDF();



%% Calculate fraction of explainable variance explained for each genotype
calcFracVar = false;
if calcFracVar
	figure;
	ffsubplot(2,1,1);
	uniqueGenos = unique(genoList(:));
	sigMatrix = signalDataMatrix(scores,classIX);
	[uniqueGenos(:) genoSortOrder(:)]
	for genoNn = 1:length(uniqueGenos)
		genoN = uniqueGenos(genoSortOrder(genoNn));
		
		ix = find(genoList == genoN);
		subScores = sigMatrix(ix,:);
		allVar = totalVariance(subScores);
		first2Var = totalVariance(subScores(:,1:2));
		
		totalExpVar = allVar;
		fracExpVar = first2Var/allVar;
		scatter(totalExpVar,fracExpVar,...
			'MarkerEdgeColor','none',...
			'MarkerFaceColor',colorSpace(colorIXList(genoNn),:)); hold on;
		text(totalExpVar+.1,fracExpVar-.01,num2str(genoNn));
	end

	xlabel('Total Explainable Variance (AU)');
	ylabel('Fraction Explained by First 2 PCs');
	return;
end
%% End Calc frac var


	for eigN = 1:4
		V = coeffs(:,eigN);
		[B,IX] = sort(abs(V),'descend');
		disp(['PC #',num2str(eigN)]);
		disp('--------------------');
		for loadN = 1:nPCs
			disp([num2str(V(IX(loadN)),'%+1.3f'),' - #',num2str(IX(loadN)),' - ',metricLabels{IX(loadN)},'  ', num2str(gateRawOrder(IX(loadN)))]);
		end
		disp(' ');
	end



	transScoresByFly = scores;
	[transSigCov, transNoiseCov] = signalCovariance(transScoresByFly, classIX);

	ffsubplot(2,3,4);
	image(transSigCov./expVar,'CDataMapping','scaled');
	title('Signal cov. after sigPCA'); colorbar;
	xlim([.5 10.5]);
	ylim([.5 10.5]);
	axis square;

	ffsubplot(2,3,5);
	image(transNoiseCov./expVar,'CDataMapping','scaled');
	title('Noise Cov. after sigPCA'); colorbar;
	axis square;

	FF = ffsubplot(2,3,6);
	image((transSigCov + transNoiseCov)./expVar,'CDataMapping','scaled');
	title('Total cov. after sigPCA'); colorbar;

	FF.setTitle(plotTitle);	
	
	[orcoAx, nullAx] = nullDimension(transScoresByFly, meanIX,scoresByFly,coeffs,sigCorr, decPI, gateRawOrder,metricLabels);

return;

	orcoAx
	nullAx

	allGenos = unique(cGenoList);
	figure;
	for genoNn = 1:length(allGenos)
		genoN = allGenos(genoNn);

		for powerN = 1:8
			ix = find((cGenoList == genoN) & (cPowerNList == powerN));
			t1(powerN) = mean(cX(ix,1));
			t2(powerN) = mean(cX(ix,2));
			t3(powerN) = mean(cX(ix,3));
		end
		
		T = [t1(:),t2(:),t3(:)];
		T = [T*orcoAx, T*nullAx];

		T
		ffsubplot(3,2,1);
		plot(T(:,1),T(:,2),'.-','Color',pretty(genoNn)); hold on;
		plot(T(1,1),T(1,2),'o','Color',pretty(genoNn)); hold on;
	end

	return;

	figure;
	scatterByPower(transScoresByFly, meanIX,scoresByFly,coeffs,sigCorr, decPI);
	
	figure;
	plotRaster(transScoresByFly, meanIX, decPI);

	function shuffledScores = permuteScores(scoresByFly, classIX)

		nClasses = length(unique(classIX));
		for dimN = 1:size(scoresByFly,2)
			classMapping = randperm(nClasses);
			for destClass = 1:nClasses

				destIXs = find(classIX == destClass);
				sourceClass = classMapping(destClass);
				sourceIXs = find(classIX == sourceClass);
				sourceIX = sourceIXs(randi(length(sourceIXs),length(destIXs),1));
				shuffledScores(destIXs(:),dimN) = scoresByFly(sourceIX(:),dimN);
			end
		end




