function ims = imgshift(im, shift)
% shift image according to shift vector
[m, n] = size(im); % dimension of image
xs = shift(1);  % shift horizontally
ys = shift(2);  % shift vertically
ims = im;

% source, horizontal range
x1s=max(1, 1+xs);
x2s=min(n, n+xs);

% source, vertical range
y1s=max(1, 1+ys);
y2s=min(m, m+ys);

% destination, horizontal range
x1=max(1, 1-xs);
x2=min(n, n-xs);

% destination, vertical range
y1=max(1, 1-ys);
y2=min(m, m-ys);

ims(y1:y2, x1:x2)=ims(y1s:y2s, x1s:x2s);
end