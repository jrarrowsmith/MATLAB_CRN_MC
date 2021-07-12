function [most_common_rate, max_rate, min_rate] = clean_and_process_pdf(h, sigmas, run_name, x)
%funcation taks input histogram of pdf of rates and normalizes it and then
%clips it by the amount in sigmas
%J R Arrowsmith, January 2021
%x is just the x location of the text for the plot

MFC=[1,1,1];
bincenters = h.BinEdges(1):h.BinWidth:(h.BinEdges(h.NumBins)+h.BinWidth);
a = size(bincenters);

wbc=h.BinCounts;
b= size(wbc);

%this is a little hack to fix why sometimes the lengths of these vectors
%don't match
if a(2)~=b(2)
    bincenters = h.BinEdges(1):h.BinWidth:(h.BinEdges(h.NumBins));
end

area_of_each_bar=wbc.*h.BinWidth;
orginal_area_under_curve=sum(area_of_each_bar);
normalized_area_of_each_bar = area_of_each_bar./orginal_area_under_curve;
plot(bincenters, normalized_area_of_each_bar,'k-')
na = sum(normalized_area_of_each_bar);

%adjust for significance
while na>sigmas
    tf=find(wbc>0);
    temp=wbc(tf)-1;
    wbc(tf)=temp;
    area_of_each_bar=wbc.*h.BinWidth;
    area_under_curve=sum(area_of_each_bar);
    normalized_area_of_each_bar = area_of_each_bar./orginal_area_under_curve;
    na=sum(normalized_area_of_each_bar);
end

%find the mean
t=0; mean_pos=0;
while t<=sigmas/2
    mean_pos=mean_pos+1;
    t=t+normalized_area_of_each_bar(mean_pos);
end
mean_rate=bincenters(mean_pos);
t/sigmas;

hold on

fill(bincenters, normalized_area_of_each_bar,[192/255 192/255 192/255])
plot(mean_rate,normalized_area_of_each_bar(mean_pos),'d',...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',MFC)
max_prob=max(normalized_area_of_each_bar);
loc=find(normalized_area_of_each_bar==max_prob);
most_common_rate=bincenters(loc);
plot(most_common_rate,max_prob, 'o',...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',MFC)
tf=normalized_area_of_each_bar>0;
min_rate = min(bincenters(tf))

%loc=find(bincenters==min_rate)-1  %this is the first non-zero but we want
%the last zero old way
loc=find(bincenters==min_rate);  %this is the first non-zero but we want the last zero
min_rate = bincenters(loc);
plot(min_rate,normalized_area_of_each_bar(loc), 's',...
    'MarkerSize',10,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',MFC)
max_rate=max(bincenters(tf));
loc=find(bincenters==max_rate)+1;%this is the las non-zero but we want the first zero on the high side
if loc>length(bincenters)
    loc=length(bincenters);
end
max_rate=bincenters(loc);
plot(max_rate,normalized_area_of_each_bar(loc), 's',...
    'MarkerSize',10,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',MFC)


atxt=[];
atxt = sprintf('Most common rate = %0.4f\nsig range = %0.4f\nmin = %0.4f to max = %0.4f\nMean rate = %0.4f',most_common_rate,sigmas,min_rate,max_rate,mean_rate);
%atxt=strcat({atxt},{run_name});
%[A]=gtext(atxt)

text(x,max_prob,atxt, 'VerticalAlignment','top')
ylabel('probability')

end