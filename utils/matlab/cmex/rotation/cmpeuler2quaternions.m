function cmp

FileName = 'mexeuler2quaternions';

if ispc
    strcomp = [' mex -I..\lib ',FileName,'.c ', ...
 		'..\lib\gmatrix.c ..\lib\gmatrix_linalg.c ..\lib\gmatrix_matlab.c ..\lib\rotation.c',' -output ',FileName]
else
    strcomp = [' mex -I../lib ',FileName,'.c ', ...
 		'../lib/gmatrix.c ../lib/gmatrix_linalg.c ../lib/gmatrix_matlab.c ../lib/rotation.c',' -output ',FileName]
end
disp(sprintf('Gerando %s...',FileName)); eval(strcomp); disp('Concluido!');
