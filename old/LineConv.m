function l=LineConv(pectoral)

%[y1 x1; y2 x2]=pectoral;

Point1=pectoral(1,:);
Point2=pectoral(2,:);

% [y1,x1]=Point1;
% [y2,x2]=Point2;

l=cross([Point1(2),Point1(1),1],[Point2(2),Point2(1),1]);