function cmp

FileName = 'mexlocalization';

strcomp = [' mex -I..\lib ',FileName,'.c ', ...
 		'..\lib\localization.c ..\lib\kalman.c ..\lib\gmatrix.c ..\lib\gmatrix_linalg.c ..\lib\gmatrix_matlab.c ..\lib\rotation.c',' -output ',FileName]
disp(sprintf('Gerando %s...',FileName)); eval(strcomp); disp('Concluido!');
