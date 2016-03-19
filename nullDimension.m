function [orcoAx, nullAx] = nullDimension( X, meanIX, scoresByFly, coeffs, sigCorr, decPI, gateRawOrder, metricLabels)

	laserPowers = [1,2,4,8,12,16,32,64].*(1.5/64);
	xTicks = [1.5/64, [.045:.015:.15],[.15+.15:.15:1.5]];

	genoList = unique(meanIX(:,3));
	sortVals = zeros(length(genoList),1);

	% Find an ordering for the genoList
	for genoNn = 1:length(genoList)
		genoN = genoList(genoNn);

		ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == 8));
		sortVals(genoNn) = nanmean(decPI(ix));
	end
	[B, genoSortOrder] = sort(sortVals,'descend');
	colorSpace = flipud(jet);
	colorIXList = round(linspace(1,64,length(genoList)));

	% Find the direction that controls vary in
    genoN = 24;		
	for powerN = 1:8 
		ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == powerN));
		n1(powerN) = mean(X(ix,1));
		n2(powerN) = mean(X(ix,2));
		n3(powerN) = mean(X(ix,3));
	end
%	nullDir = [n1(:),n2(:),n3(:)];
	nullDir = [n1(:),n2(:),0*n3(:)];
	[coeff, score] = princomp(nullDir);
	nullAx = coeff(:,1);
	nullAx = -nullAx./norm(nullAx);

	orcoAx = -[-nullAx(2),nullAx(1),0]';

%	% Find the direction orco varies in
%    genoN = 18;		
%	for powerN = 1:8 
%		ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == powerN));
%		o1(powerN) = mean(X(ix,1));
%		o2(powerN) = mean(X(ix,2));
%		o3(powerN) = mean(X(ix,3));
%	end
%%	orcoDir = [o1(:),o2(:),o3(:)];
%	orcoDir = [o1(:),o2(:),0*o3(:)];
%	% Remove the component of null variation
%	nonNullOrcoDir = orcoDir - ((orcoDir*nullAx)*nullAx');
%	% Find the first non-control component of orco, normalized
%	[coeff, score] = princomp(nonNullOrcoDir);
%	orcoAx = -1*coeff(:,1);
%	orcoAx = 1*orcoAx./norm(orcoAx);


%	ax2 = -cross(nullAx, orcoAx);
	ax2 = [0,0,1]';

	% nullScore = X(:,[1:3])*nullDir(:);

	allOdorAx = X(:,1:3)*orcoAx;
	allThermalAx = X(:,1:3)*nullAx;
	allGenoAx = X(:,1:3)*ax2;
	odorCorr = nancorr(decPI(:),allOdorAx);
	thermalCorr = nancorr(decPI(:),allThermalAx);
	genoCorr = nancorr(decPI(:),allGenoAx);

	figure();
	for genoNn = 1:length(genoList)
		genoN = genoList(genoSortOrder(genoNn));
		plotColor = colorSpace(colorIXList(genoNn),:);

%		disp(' ');
%		genoNn
%		genoN
%		plotColor
		
		for powerN = 1:8 
			ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == powerN));
			t1(powerN) = mean(X(ix,1));
			t2(powerN) = mean(X(ix,2));
			t3(powerN) = mean(X(ix,3));
		end

		T = [t1(:),t2(:),t3(:)];

		% Remove the null variation
		% T = T - ((T*nullAx)*nullAx');

		% Rotate s/t orcoAx is in X, and nullAx is in Y by projecting onto [orcoAx,ax2,nullAx]
		T = [T*orcoAx,T*ax2,T*nullAx];


		ffsubplot(3,2,1);
		plot(T(:,1),T(:,2),'.-','Color',plotColor); hold on;
		plot(T(1,1),T(1,2),'o','Color',plotColor);
		xlabel('Odor Axis'); ylabel('PC3');

		ffsubplot(3,2,2);
		plot(t1(:),t2(:),'.-','Color',plotColor); hold on;
		plot(t1(1),t2(1),'o','Color',plotColor); hold on;

		% Plot the axes
		plot([orcoAx(1) 0 nullAx(1)],[orcoAx(2) 0 nullAx(2)],'k');
		xlabel('PC1'); ylabel('PC2');
		axis auto;
		xlim(xlim()*1.1); ylim(ylim()*1.1);
		axis equal;

		ffsubplot(3,2,3);
		semilogx(laserPowers,T(:,3),'.-','Color',plotColor); hold on;
		semilogx(laserPowers(1),T(1,3),'o','Color',plotColor);
		xlabel('Laser power'); ylabel('Thermal Axis');
		xlim([.75 96].*1.5/64);
		set(gca,'XTick',xTicks,'XTickLabel',{});
		title(['Correlation with decPI: ',num2str(thermalCorr)]);

		ffsubplot(3,2,5);
		semilogx(laserPowers,T(:,1),'.-','Color',plotColor); hold on;
		semilogx(laserPowers(1),T(1,1),'o','Color',plotColor);
		xlabel('Laser power'); ylabel('Odor Axis');
		xlim([.75 96].*1.5/64);
		set(gca,'XTick',xTicks,'XTickLabel',{});
		title(['Correlation with decPI: ',num2str(odorCorr)]);

		ffsubplot(3,2,6);
		semilogx(laserPowers,T(:,2),'.-','Color',plotColor); hold on;
		semilogx(laserPowers(1),T(1,2),'o','Color',plotColor);
		xlabel('Laser power'); ylabel('PC3');
		xlim([.75 96].*1.5/64);
		set(gca,'XTick',xTicks,'XTickLabel',{});
		title(['Correlation with decPI: ',num2str(genoCorr)]);

%		plot3(T(:,1),T(:,2),T(:,3),'.-','Color',pretty(genoNn)); hold on;
%		plot3(T(1,1),T(1,2),T(1,3),'o','Color',pretty(genoNn));

	end

%% Do hypothesis testing on Orco Axis
	% Get Control scores for comparison first
	genoN = 24; 
	ix = find((meanIX(:,3) == genoN));
	flyNs = unique(meanIX(ix,2));
	ctrlScores = zeros(length(flyNs),8);
	for flyNn = 1:length(flyNs)
		flyN = flyNs(flyNn);
		for powerN = 1:8
			ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == powerN) &...
					  (meanIX(:,2) ==  flyN));
			ctrlScores(flyNn,powerN) = X(ix,1:3)*orcoAx;
		end
	end

	% Now get scores for each genotype
	for genoNn = 1:length(genoList)
		genoN = genoList(genoSortOrder(genoNn));
		plotColor = colorSpace(colorIXList(genoNn),:);
		
		% Find how many flies there are
		ix = find((meanIX(:,3) == genoN));
		flyNs = unique(meanIX(ix,2));
		orcoAxScores = zeros(length(flyNs),8);
		for flyNn = 1:length(flyNs)
			flyN = flyNs(flyNn);
			for powerN = 1:8
				ix = find((meanIX(:,3) == genoN) & (meanIX(:,1) == powerN) &...
						  (meanIX(:,2) ==  flyN));
				orcoAxScores(flyNn,powerN) = X(ix,1:3)*orcoAx;
			end
		end

		% Bootstrap P-values against the no-Gal4 controls
		Pvals(genoNn) = bootstrapP(ctrlScores, orcoAxScores, 10000);
%		figure;
%		plot(orcoAxScores');
%		title([num2str(genoN),' ',num2str(Pvals(genoNn))]);
	end
	% Adjust for multiple comparisons
	Padj = bonferroniHolm(Pvals);
%% End hypothesis testing on Orco Axis


	orcoAx = orcoAx;
	genoAx = ax2;
	thermalAx = nullAx;
%	save('~/gateLearn/pcAxes.mat','orcoAx','genoAx','thermalAx','coeffs','gateRawOrder');

goodMetricLabels = metricLabels;
	% Load this to get genotype names
	load('~/networkDiagrams/grandResults.mat');
	% ffsubplot(3,2,4);
	figure;
	ffsubplot(1,1,1);
	plot(0,0); 
	ylim([-(length(genoList)/2 + 2) 2 ]);
	xlim([-.25 4.25]);
	for genoNn = 1:length(genoList)
		nameString = summaryResults{genoList(genoSortOrder(genoNn))}{1};
		plotColor = colorSpace(colorIXList(genoNn),:);
		% text(2*floor((genoNn-1)/8),-mod(genoNn-1,8),nameString,'Color',plotColor);
		text(2*floor((genoNn-1)/24),-mod(genoNn-1,24),nameString,'Color',plotColor);
		hold on;

		output{genoNn} = [' ',...
			num2str(Pvals(genoNn),'%3.5f'),'   ',...
			num2str( Padj(genoNn),'%3.5f'),'   ',...
			nameString];
	end
	set(gca,'XTick',[],'YTick',[]);
	metricLabels = goodMetricLabels;

	disp({'Raw P','Adjusted P','Genotype'});
	disp('-----------------------------------------------');
	for genoNn=1:length(genoList)
		disp(output{genoNn});
	end	
	disp(' ');

	orcoAxScores = X(:,1:3)*orcoAx;
	genoAxScores = X(:,1:3)*ax2;
	thermalAxScores = X(:,1:3)*nullAx;




	ptSize = 16;

	figure;
	ffsubplot(3,2,1);
	scatter(orcoAxScores,genoAxScores,ptSize,decPI(:));
	xlabel('Odor Axis'); ylabel('PC3');
	title('Color: decPI');
	
	ffsubplot(3,2,3);
	scatter(orcoAxScores,thermalAxScores,ptSize,decPI(:));
	xlabel('Odor Axis'); ylabel('Thermal Axis');
	title('Color: decPI');

	ffsubplot(3,2,4);
	scatter(genoAxScores,thermalAxScores,ptSize,decPI(:));
	xlabel('PC3 Axis'); ylabel('Thermal Axis');
	title('Color: decPI');


	corrMatrix(:,1) = nancorr(scoresByFly,orcoAxScores);
	corrMatrix(:,2) = nancorr(scoresByFly,genoAxScores);
	corrMatrix(:,3) = nancorr(scoresByFly,thermalAxScores);
	corrMatrix(isnan(corrMatrix)) = 0;

	nLoads = 12;
	for dimN = 1:3
		[B,IX] = sort(abs(corrMatrix(:,dimN)),'descend');
		disp(['Dim #',num2str(dimN)]);
		disp('---------------------');
		for loadN = 1:nLoads
			disp(['IX: ',num2str(IX(loadN)),'   ',...
				  'G',num2str(gateRawOrder(IX(loadN))),'   ',...
				  num2str(sigCorr(IX(loadN))),'   ',...
				  num2str(corrMatrix(IX(loadN),dimN))]);
		end
		disp(' ');
		disp(' ');
	end

	figure;
	ffsubplot(1,4,1);
	nScores = min([size(scoresByFly,2) 64]);
	R = diag(sigCorr);
	h = barh(R(1:nScores),'hist');
	set(h,'EdgeColor','none');
	set(gca,'YDir','reverse');
	title('Reliability');
	ylim([.5 nScores+.5]);
	xlim([0 .55]);
	if nScores < 16
		for n = 1:nScores
			text(0, n, metricLabels{gateRawOrder(n)},...
				'HorizontalAlignment','left','VerticalAlignment','middle');
		end
	end


%	ffsubplot(5,1,2);
%	PIcorr = corr(scoresByFly,decPI(:));
%	bar(PIcorr(1:nScores),'hist');
%	title('PIcorr');
%	xlim([.5 nScores+.5]);

	ffsubplot(1,4,2);
	OrcoCorr = corr(scoresByFly,orcoAxScores);
	h = barh(OrcoCorr(1:nScores),'hist');
	set(h,'EdgeColor','none');
	set(gca,'YDir','reverse');
	title('OrcoAx Corr');
	ylim([.5 nScores+.5]);
	
	ffsubplot(1,4,3);
	GenoCorr = corr(scoresByFly,genoAxScores);
	h = barh(GenoCorr(1:nScores),'hist');
	set(h,'EdgeColor','none');
	set(gca,'YDir','reverse');
	title('GenoAx Corr');
	ylim([.5 nScores+.5]);
	
	ffsubplot(1,4,4);
	ThermalCorr = corr(scoresByFly,thermalAxScores);
	h = barh(ThermalCorr(1:nScores),'hist');
	set(h,'EdgeColor','none');
	set(gca,'YDir','reverse');
	title('ThermalAx Corr');
	ylim([.5 nScores+.5]);


	orcoSpec = abs(OrcoCorr)./(abs(GenoCorr) + abs(ThermalCorr));
	genoSpec = abs(GenoCorr)./(abs(OrcoCorr) + abs(ThermalCorr));
	thermalSpec = abs(ThermalCorr)./(abs(OrcoCorr) + abs(GenoCorr));

	ix = find(isnan(orcoSpec)); orcoSpec(ix) = 0;
	ix = find(isnan(genoSpec)); genoSpec(ix) = 0;
	ix = find(isnan(thermalSpec)); thermalSpec(ix) = 0;


	DescMatrix = [[1:length(orcoSpec)]', gateRawOrder(:), OrcoCorr(:),GenoCorr(:),ThermalCorr(:),...
						  orcoSpec(:), genoSpec(:), thermalSpec(:)];
	ix = find(isnan(DescMatrix)); DescMatrix(ix) = 0;

	[B,orcoIX] = sort(orcoSpec,'descend');
	[B,genoIX] = sort(genoSpec,'descend');
	[B,thermalIX] = sort(thermalSpec,'descend');

	disp('Orco Specificity');
	DescMatrix(orcoIX(1:12),:)

	disp('Geno Specificity');
	DescMatrix(genoIX(1:12),:)

	disp('Thermal Specificity');
	DescMatrix(thermalIX(1:12),:)




%	figure;
%	scatter3(OrcoCorr, GenoCorr, ThermalCorr,'o');
