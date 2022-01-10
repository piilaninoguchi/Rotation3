%%
%dispOM3D

%Currently analyzes 2D data.
%Given cell image (and a b/w image showing where the cell area is) and orientation map, show:
% - local orientation varianece ie degree of parallel alignment between
%   between MTs
% - alignment with neatrest orientation ot nearest point on the cell border
%   (or, optionally, center)
% - histogram of of the alignment at the different distance ranges

%code commented out contained in the actual working progrsm is removed here, so
%optional and alternative implemenatation is missing

%This actually doesn't call anny custom helpers, but the orientation map is
%prepared separately (currently using functions like BlkSVDOrient3D, all by Xiaoguang Feng)

% I'd like to imporove the terminology used in the comments to be
% consistent
%%
clear all
path(path,'/Developer/matlabdev/diffusionXiaodong/CircStat2012a');
path(path,'/Developer/matlabdev/diffusionXiaodong');
path(path,'/Developer/matlabdev/diffusionXiaodong/Data');

%LoadParamData
%dataSet = 1;%which cell's photo, mask and orientaino map files to use
noisePar=.0001;
nhoodsize=6;%radius of pixel neighborhood 
nbins=4;%in degree of alignment histogram


% 
% switch dataSet %load existing cell data from different files 
%                 %As new data files appear, add cases with different file names and.
%     case 1
%         
        
        [Idir,Ifile]=uigetfile('*','Select image');
        
        [maskdir,maskfile]=uigetfile('*','Select mask');
        
        I = imread([Idir,Ifile]);
        mask = double(imread([maskdir,maskfile]));
%         
%         
% % I=imread('/Developer/matlabdev/diffusionXiaodong/Data/2point8mMGlu_SIR_ALX.jpg');%phtoto itself, called I
% % mask = double(imread('/Developer/matlabdev/diffusionXiaodong/Data/2point8mMGlu_SIR_ALXmask.bmp'));%mask to show one specific cell in the photo
% % load('OM2point8mMGlu_SIR_ALX2.mat');%orientation map of the photo, generated beforehand, variable called OM
% end

%PrepVars
%PrepData


fprintf(1, 'Loading orientation map... ');


%I=double(imnoise(uint8(I),'gaussian',0,noisePar));
I = double(I);

OM=BlkSVDOrient(I,1,3,0);


fprintf(1, '  Done\n');





%PrepDataSamples
%involves some manual approxiamtion

%determine which pixels show an MT:
TholdToMaskAreaRatio =.3291;
npixesInMask = sum((mask(:)/255));
tholds=[0.01:.025:.5];
a1 = repmat(permute(tholds,([3 1 2])),size(mask));
a2 = repmat(I/255,[1,1,length(tholds)]) >a1;
repmask = repmat(mask,[1,1,length(tholds)])/255;
a3 = sum(sum(a2.*repmask));
a3=permute(a3,([1 3 2]));
a4= a3 ./npixesInMask - TholdToMaskAreaRatio;
[~,a5]=min(abs(a4));
thold = tholds(a5);%

npixesInData = sum(((I(:)/255)>thold).*logical(mask(:)));
 
constantNpixels = 34546;
constantNpixels = 26000;
constantNpixels = 50000;

percentpixels = constantNpixels/npixesInData;

if percentpixels>1

    percentpixels = 1;
    
    disp('warining - too few piexls in data');

    
end

I=imresize(I,size(OM),'bilinear');%all loaded data images must be the same size


mask = double(logical(mask));

mask = 1-double((imresize(mask,size(OM),'nearest')));%also, must me correct value range ()



c=repmat(mask,[1,1,3]);

%%
%RootateData
if false%Below is one way of treating orientaion analysis. Currently not used, all prepared in advance, loaded above.
I0=I;
hists=[];
bins=[-pi/2 + (pi/360)/2 : pi/360: pi/2 - (pi/360)/2];
for Actr = 0:1:359

  I= imrotate(I0,Actr,'nearest','crop');
  OM=BlkSVDOrient(double(I),1,2,0);
  angs=atan2(real(OM),imag(OM));
  ex=(angs.*((I/255)>thold).*~mask);
  hists(end+1,:)=hist(ex(ex~=0),bins);
  fprintf(1, [num2str(Actr), ' ']);
end
figure
imagesc(hists);
return
end

%%
%FindMaskBound

%determine pixels that make up the cell border:
[extmp ,wytmp] = meshgrid((1:size(I,1)),(1:size(I,2)));
maskinds1 =  wytmp(~logical(mask(:)));%coordinates of pixels that make up the cell
maskinds2 =  extmp(~logical(mask(:)));
maskindslinear = sub2ind(size(I),maskinds1,maskinds2);%1D indecies of above 2D coordinates 
boundscell=bwboundaries(mask)';%analysis to find borders in the mask image
cs=cell2mat(cellfun(@size,boundscell,'uni',false)');
[~,whichCell]=min(cs(:,1)); %bwboundaries result has two parts - actual cell border and a rectangular outline of the image
boundsinds = boundscell{2};                                         %so we want to chose the right one
boundsindsx = boundsinds(:,1);                                      %actually not guaranteed to work automatically
boundsindsy = boundsinds(:,2);                                      %check for each cell manually to be sure.
%because cell was outlained manually, the bounary is masy and crooked,
%smooth it:
smboundsx=smooth(boundsindsx,0.05,'loess');smboundsy=smooth(boundsindsy,0.05,'loess');
boundsO = atan(-diff(smboundsx)./diff(smboundsy))+pi/2;%orientation at points on border
boundsO(end+1)=boundsO(1)-boundsO(end);%add connecting point
boundsO=mod(boundsO,pi);%to wrap around connecting point;
%now create images from the above coordinate lists
boundsOmap = zeros(size(mask));
boundsmap = zeros(size(mask));
boundsindslinear=sub2ind(size(I),boundsindsx,boundsindsy);
boundsOmap(boundsindslinear)=boundsO;
boundsmap(boundsindslinear)=1;
bounDists = zeros(size(I));

%%
%deal with orientaion map, combine with mask, etc:
%OM comtains complex numbers, get the angle values from it
cm=colormap('hsv');
angs=atan2(real(OM),imag(OM));
angsNaN=angs;%remember to add to mask
angs2(isnan(angs))=0;
angs2=(angs-min(angs(:)));
angs2=(angs2/max(angs2(:)));
angs3=imresize(angs2,size(I),'nearest');
angsrgb=ind2rgb(floor(angs3.*63+1),cm);
I2 = I/255;%imresize(I,size(angs),'nearest');
I3 = repmat(I,[1 1 3])./255;
I4=((I3).*angsrgb);
ex=((I3).*angsrgb);
ex2=rgb2ind(ex,64);
 Omap = double(rgb2ind((I3).*angsrgb,64));
 
 %more set up including, neighborhood kernel
 nhoodKernel1=[ 1 2 3;1 2 3; 1 2 3];
 nhoodKernel2=[1 1 1;2 2 2;3 3 3];
 nhoodKernel12=nhoodKernel1-2;
 nhoodKernel22=nhoodKernel1-2;
 OmapPadded = padarray(2*angs,[nhoodsize nhoodsize],'circular');
 OmapPadded(isnan(OmapPadded))=0;

 [x,y]=meshgrid(-nhoodsize:nhoodsize,-nhoodsize:nhoodsize);
 circKernel = sqrt(x.^2+y.^2)<nhoodsize;
 circKernel(nhoodsize:nhoodsize+2,nhoodsize:nhoodsize+2)=0;
 
 angVarmat = zeros(size(OM));
  isnans=zeros(size(OM));
 Ithold = double(I2>thold);
 ItholdPadded = padarray(Ithold,[nhoodsize nhoodsize],'circular');
 Ipadded = padarray(I2,[nhoodsize nhoodsize],'circular');
tmpIthold = zeros(size(OM));
tmpItholdCum= zeros(size(OM));
divisorCum=zeros(size(OM));

%%
fprintf(1, 'Calculating local circular variance... ');
%angular variance is simply the sum of angles, larger means more of
%a particular direction, smaller - more evenly distributed directions
  for ctr1=-nhoodsize:nhoodsize
      for ctr2=-nhoodsize:nhoodsize
          
          %not part of clean comments, but better leave this comment here, might be important:
          % threshold values!!!!!!!!!!!!!!!!!!!!!!!!
          
        if circKernel(ctr1+nhoodsize+1,ctr2+nhoodsize+1)%shape of neighborhood may be square or circular, the latter means we skip some pixels
      %see if give neighbor pixel is part of an MT, use only those which
      %are
      tmpIthold = ItholdPadded((nhoodsize+1:end-nhoodsize)+ctr1,(nhoodsize+1:end-nhoodsize)+ctr2);
      tmpItholdCum = tmpItholdCum+tmpIthold;
      divisortmp = Ipadded((nhoodsize+1:end-nhoodsize)+ctr1,(nhoodsize+1:end-nhoodsize)+ctr2);
      divisorCum = divisorCum + divisortmp;
      tmp=OmapPadded((nhoodsize+1:end-nhoodsize)+ctr1,(nhoodsize+1:end-nhoodsize)+ctr2);
   
      angVarmattmp =  tmpIthold .* exp(1i.*tmp);
      diffTmp = angs-tmp;
      morethanpi =  diffTmp>pi/2;
      lessthanmpi = diffTmp < -pi/2;
      diffTmp(morethanpi) =diffTmp(morethanpi)- pi;
      diffTmp(lessthanmpi) =diffTmp(lessthanmpi)+ pi;
      angVarmat=angVarmat+angVarmattmp;
        end
      end
  end
  
  %now we're curious about the result, so check out the result
  %might want to look at more than one version of it
  ex=(abs(angVarmat(:))./(tmpItholdCum(:)+1));
  meanAngVarmat= mean(ex(logical(Ithold)&~mask));
  fprintf(1, '  Done\n');
  tmpIthold3=repmat(tmpItholdCum,[1,1,3]);
  angVarmat0=angVarmat;
  angVarmat=abs(angVarmat)./(tmpItholdCum+1);
  figure;imagesc(((I3).*angsrgb).*repmat(1-mask/255,[1,1,3]))
  
  %%
%FindDists
%find distances from nearest border point a given pixel 
cellcenter= [mean(maskinds2) mean(maskinds1)];
maskPlusCenter=mask;
maskPlusCenter(floor(cellcenter(2)):floor(cellcenter(2))+4,floor(cellcenter(1)):floor(cellcenter(1))+4)=255;
distsFromC= zeros(size(Omap));
distsFromC0= zeros(size(Omap));
bounDists=zeros(size(Omap))-1;
nearestOs=zeros(size(Omap));

fprintf(1, 'Calculating distance from border map... ');%... and the orientation of the neares point on boundary
 for indctr = 1:length(maskinds1)%for each pixel in cell
     
     angVSdist(1,indctr)=Omap(maskinds1(indctr), maskinds2(indctr));%not used
     repmat([maskinds1(indctr) maskinds2(indctr)],[length(boundsindsx),1]);%not used
     distsTmp = [ boundsindsx boundsindsy] - repmat([maskinds1(indctr) maskinds2(indctr)],[length(boundsindsx),1]);                 
     [minVal, minI] = min(sum(distsTmp.^2,2).^.5);%actual distance
     bounDists(maskinds1(indctr), maskinds2(indctr))=minVal;
     
     angVSdist(2,indctr)=minVal;%not used
     
     %we also look at oter distances: from cell center
     distsFromC0(maskinds1(indctr), maskinds2(indctr)) =...
         ((maskinds2(indctr)-cellcenter(1)).^2 + (maskinds1(indctr)-cellcenter(2)).^2).^.5;
     
     distsFromC(maskinds1(indctr), maskinds2(indctr)) =...
         ((boundsindsy(minI)-cellcenter(1)).^2 + (boundsindsx(minI)-cellcenter(2)).^2).^.5;
     
     %orientation at the nearest boundary pixel:
     nearestOs(maskindslinear(indctr)) = boundsOmap(boundsindslinear(minI));
 end
 
 
 fprintf(1, '  Done\n');
 
 
 
fprintf(1, 'Binning and analysis... ');
 
%finalize the result
angs= (angs+pi/2);
angs= mod(angs+pi/2,pi);
 nearestOs = mod(nearestOs+pi/2,pi);
nearestOsrgb =ind2rgb(floor(nearestOs/pi.*63+1),cm);
 ex1=angs(~mask&(I>50));
ex2=nearestOs(~mask&(I>50));
ex5=bounDists(~mask&(I>50));
exdiff=ex1-ex2+pi/2;
 diffs=angs-nearestOs;
 
 cm2=cm;
 cm2(1,:)=[0 0 0];

 
diffsrgb=ind2rgb(floor((abs(diffs))/(pi).*63+1),cm2);
%figure;

%%count # of pixels in data, those in diffsrgbBool0==1
%old syntax - needs purging
diffsrgb3=(((I3>thold).*diffsrgb).*repmat(1-mask,[1,1,3]));
diffsrgbBool=(sum(diffsrgb3,3)==0);
diffsrgbBool0=~logical(sum(diffsrgbBool,3));
diffsrgbBool=repmat(diffsrgbBool,[1,1,3]);
diffsrgb3(diffsrgbBool)=1;
% imagesc(diffsrgb3);
% colormap(cm2);
% imagesc(((I3>.3).*diffsrgb).*repmat(1-mask/255,[1,1,3]));





%Show the map of local angles relative to that of border
%parallel to border is one extreme, perpendicular the other
ex2=((diffs).*~mask.*Ithold);
ex2(ex2>pi/2)= ex2(ex2>pi/2)-pi;
ex2(ex2<-pi/2)= ex2(ex2<-pi/2)+pi;
ex2=abs(ex2);
diffsReady= ex2;
cmj=colormap('jet');
cmj(end,:)=[.5 .5 .5];
ex2(ex2==0)=max(ex2(:))+.1;
figure
colormap(cmj);
imagesc(ex2);
%%
%FindDataBins
%Find degrees of alignment for <nbins> areas of different offsets form border ranges

%Initiate variables.
binSubject=diffsReady;
bounDists(isinf(bounDists))=0;
circumbin=[];
circumbinNs=[];
maxdist=max(bounDists(:));
circboundimg=0;


binarea = sum((diffsrgbBool0(:)))/nbins;% histogram bins each contain <binarea> values

%this was in previous version... why?
%binarea = sum((diffsrgbBool0(:)))/nbins;% histogram bins each contain <binarea> values


rctr = 1;
hibound=1;lowbound=0; %values are number of pixels, initially range is between 0 and one pixel from border
binN=1;

%In order to fill the <nbins> bins neatly and equally, we iterate gradually through
%groups of pixels located wthin a range of  <hibound> to <hibound> pixels
%from border (or optionally center). Is this overkill?
while hibound<maxdist

    hibound=hibound+.1; %<hibound> is raised gradually untill the number of pixels in range
                        %matches the bin size (<binarea>).
                        %Then bin is created and <lowbound> is mover to <hibound>
    circumbintmp =bounDists<hibound&bounDists>lowbound;
    circumbintmp = circumbintmp&diffsrgbBool0;
    binareatmp =  sum(circumbintmp(:));
    if binareatmp>=binarea
        lowbound=hibound;        
        
        %Take only <percentpixels> percent of pixels from zone.
        %This is a step towards normalizing results between numerous cells,
        %which have different areas.
        bintmp = binSubject(circumbintmp(:)&diffsrgbBool0(:));
        binsampleinds = randperm(length(bintmp));
        binsampletmp=bintmp(binsampleinds(1: floor(length(bintmp)*percentpixels)));
        circumbin =[circumbin; binsampletmp];
        
        circumbinNs = [circumbinNs; ones(length(binsampletmp),1)*binN ];
        binN=binN+1;
        circboundimg = circboundimg+ circumbintmp*hibound.*binSubject.*diffsrgbBool0;
    
    end   
end
%make the final, INNERMOST (#<nbins>) bin  (special case)
        bintmp = binSubject(circumbintmp(:)&diffsrgbBool0(:));
        binsampleinds = randperm(length(bintmp));
        binsampletmp=bintmp(binsampleinds(1: floor(length(bintmp)*percentpixels)));
        circumbin =[circumbin; binsampletmp];
        circumbinNs = [circumbinNs; ones(length(binsampletmp),1)*binN ];
 circboundimg = circboundimg+ circumbintmp*hibound.*binSubject.*diffsrgbBool0;
    
 
 %%
%AnalyzeAndSave
[P,ANOVATAB,STATS] =anova1(circumbin/(pi/2),(circumbinNs));
figure
[c,m,h,nms] = multcompare(STATS);

datasetOut = struct('mask',mask,'I',I,'binSubject',binSubject,'angs',angs,'circumbin',circumbin,'circumbinNs',circumbinNs,'multOut',{c,m,h,nms});
%num2str(dataSet)
datasetName=['DatasetOutRadDiff2' num2str(dataSet) '.mat'];
disp(datasetName);
disp(size(circumbinNs));
save(datasetName,'datasetOut');

return