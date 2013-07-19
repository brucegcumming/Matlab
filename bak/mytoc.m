function took = mytoc(tstart)
%took mytoc(start)
%calculates elapsed time since start (start = now)
%if no return value is requested, prints result.
took = (now-tstart) * 24 * 60 * 60;
if nargout == 0
   fprintf('Took %.3f sec\n',took);
end
       
