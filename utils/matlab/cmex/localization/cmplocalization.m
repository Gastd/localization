function cmp

FileName = 'mexlocalization';

if ispc
    strcomp = [' mex -I..\lib -DSYSTEM_WINDOWS ',FileName,'.c ', ...
 		'..\lib\localization.c ..\lib\kalman.c ..\lib\gmatrix.c ..\lib\gmatrix_linalg.c ..\lib\gmatrix_matlab.c ..\lib\rotation.c',' -output ',FileName]
else
    strcomp = [' mex -I../lib -DSYSTEM_LINUX ',FileName,'.c ', ...
 		'../lib/localization.c ../lib/kalman.c ../lib/gmatrix.c ../lib/gmatrix_linalg.c ../lib/gmatrix_matlab.c ../lib/rotation.c',' -output ',FileName]
end


disp(sprintf('Gerando %s...',FileName)); eval(strcomp); disp('Concluido!');
