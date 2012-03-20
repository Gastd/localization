function rn = rotation_from_quaternions(q,rb)

%%% Eh importante notar que essa funcao recebe o vetor 3x1 rb e retorno outro vetor 3x1 rn. Ou seja, temos de tirar e por os zeros. 

%%% o primeiro passo e determinar o conjugado de q
q_c = [q(1);-q(2:4)];

%%% agora temos de acrescentar o zero
rb_q = [0;rb];

%%% agora multiplicamos tudo
rn_q = multiplica_quaternion(q,rb_q);
rn_q = multiplica_quaternion(rn_q,q_c);

%%% agora soh falta tirar o zero de rn_q
rn = rn_q(2:4);


