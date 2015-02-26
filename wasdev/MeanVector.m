function r = MeanVector(lens, angles, varargin)
%r = MeanVector(lens, angles, ...)
%returns normalized vector lenth given list of lengths, angles
%angle is in degrees by default.   
% r = MeanVector(lens, angles, 'double') doubles the angles (for axial
% rather than directional data)


dbl = 1;
rmbase = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'double',5)
        dbl = 2;
    elseif strncmpi(varargin{j},'rmbase',5) %remove spont rate/baseilne
        rmbase = 1;
    elseif strncmpi(varargin{j},'selftest',5)
        TestMeanVector(1);
        return;
    end
    j = j+1;
end


if strncmpi(lens,'selftest',5)
        TestMeanVector(1);
        return;
end
if sum(size(angles) == size(lens)) ~= length(size(angles))
    r = NaN;
    return;
end

id = find(~isinf(lens) & ~isinf(angles));
if rmbase
    zid = find(isinf(angles));
    if ~isempty(zid)
        lens(id) = lens(id) - mean(lens(zid));
    end
end

angles = dbl .* angles .* pi/180;
r = lens .* cos(angles) + i .*lens .* sin(angles);
r = mean(r(id))./mean(abs(r(id)));

if dbl == 2
    a = angle(r);
    r = abs(r) .*  (cos(a/2) + i .* sin(a/2));
end


function TestMeanVector(type)

a = [0:10:180];
sds = [10:10:100];
for j = 1:length(sds)
    y = Gauss(sds(j),a);
    r(j) = MeanVector(y,a,'double');
end
hold off;
plot(sds,abs(r));
hold on;
cv = exp(-(sds.^2)./(2 .* 45.6462^2));
cv = sd2cv(sds./1.5);
plot(sds,cv,'r');

