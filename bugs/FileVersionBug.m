function Z = FileVersionBug(varargin)
%
%Create a file, test.mat, that gnerates warnings when loaded
%by a differenv verstion of matlab

Z.x = rand(1,10);
plot(Z.x);
Z.axis = gca;
save('test.mat','Z');
fprintf('Now load test.mat in another resease of matlab\n');
