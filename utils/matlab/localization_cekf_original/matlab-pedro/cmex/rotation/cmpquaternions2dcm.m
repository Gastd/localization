function cmp

FileName = 'mexquaternions2dcm';

strcomp = [' mex -I..\lib ',FileName,'.c ', ...
 		'..\lib\gmatrix.c ..\lib\gmatrix_linalg.c ..\lib\gmatrix_matlab.c ..\lib\rotation.c',' -output ',FileName]
disp(sprintf('Gerando %s...',FileName)); eval(strcomp); disp('Concluido!');
