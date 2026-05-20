clear;
close all;
N=31;% To keep syncronized with the code

data = readmatrix('tbasisfun_output.dat');
idx = 0:N: size(data, 1);
idx(1)=1;

spans = cell(numel(idx) - 1, 1);
for s = 1:numel(spans)
    block = data(idx(s) : idx(s + 1), :);
    spans{s} = block;
end
figure;
for s = 1:numel(spans)
    plot(spans{s}(:,1), spans{s}(:,2), '.-');
    hold on;
end
 title("BSplines on the interval");


