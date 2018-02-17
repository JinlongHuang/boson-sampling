% Classical simulation for 8 photons and 8 modes using Glynn's algorithm
% and Gurvits's algorithm. 
% 
% Both the initial state and final state are
% standard state, i.e. [1,1,1,1,1,1,1,1]. Transition amplitudes are
% calculated. 
% 
% Total Glynn's algorithm takes 45.8 minutes, so only 1/8 fraction of it is 
% executed. If number fo samples is chosen to be 5*10^5, total Gurvits's 
% algorithm takes 92s.
% 
% ©2018 Jin-Long Huang All Right Reserved
% -------------------------------------------------------------------------

% random unitary matrix in [0,1]
C = rand(8) + 1i * rand(8);
randH = (C + C')/2;
randU = expm(1i*randH)

file1 = fopen('test8.txt', 'w');

% Glyn's Algorithm (deterministic algorithm)
i0 = 1;
i1 = 1;
i2 = 1;
i3 = 1;
i4 = 1;
i5 = 1;
i6 = 1;
i7 = 1;
sum = 0+0*1i;
mGenGly = 0+0*1i;
Z = [1 + 0*1i, exp(1i*2*pi/8), exp(1i*4*pi/8), exp(1i*6*pi/8)...
    , exp(1i*8*pi/8), exp(1i*10*pi/8), exp(1i*12*pi/8),exp(1i*14*pi/8)];
while i0 <= 8
    stic0 = tic;
    while i1 <= 8
        stic1 = tic;
        while i2 <= 8
            while i3 <= 8
                while i4 <= 8
                    while i5 <= 8
                        while i6 <= 8    
                            while i7 <= 8  
                            ZVec = [Z(i0);Z(i1);Z(i2);Z(i3);Z(i4);Z(i5);Z(i6);Z(i7)];
                            mGenGly = Z(i0)'*Z(i1)'*Z(i2)'*Z(i3)'*Z(i4)'...
                            *Z(i5)'*Z(i6)'*Z(i7)'*dot(randU(1,:),ZVec)...
                            *dot(randU(2,:),ZVec)...
                            *dot(randU(3,:),ZVec)*dot(randU(4,:),ZVec)...
                            *dot(randU(5,:),ZVec)*dot(randU(6,:),ZVec)...
                            *dot(randU(7,:),ZVec)*dot(randU(8,:),ZVec);                            
                            sum = sum + mGenGly;  
                            i7 = i7 + 1;
                            end
                        i7 = 1;    
                        i6 = i6 + 1;
                        end
                        i6 = 1;
                        i5 = i5 + 1;
                    end
                    i5 = 1;    
                    i4 = i4 + 1;
                end                    
                i4 = 1;
                i3 = i3 + 1;
            end
            i3 = 1;
            i2 = i2 + 1;
        end
        i2 = 1;            
        stoc1 = toc(stic1);
        fprintf(file1, 'i1   %3d   %14.8e\n', i1, stoc1);
        transAmp = abs(sum)/(i1*8^6);
        fprintf(file1, 'abs(sum) is %14.8e\n', abs(sum));
        fprintf(file1, 'transAmp average over (%2d * 8^6) is %14.8e\n', i1, transAmp);
        i1 = i1 + 1;
    end
    i1 = 1;    
    break;
    stoc0 = toc(stic0);
%     fprintf(file1, 'i0   %3d   %14.8e\n', i0, stoc0);
%     transAmp = abs(sum)/(i0*8^7);
%     fprintf(file1, 'abs(sum) is %14.8e\n', abs(sum));
%     fprintf(file1, 'transAmp average over (%2d * 8^7) is %14.8e\n \n ', i0, transAmp);
    i0 = i0 + 1;
end
transAmp2 = abs(sum)/(8^7);
fprintf(file1, '\nTransition amplitude for Glynn from 8^7 samples is %14.8e\n', transAmp2);
fprintf(file1, '\n1/8 time for Glynn is   %14.8e\n', stoc0);


% Gurvits' algorithm (Sampling algorithm)
tic2 = tic;
T = 5*10^5;
count = 1;
sum = 0 + 0*1i;
while count <= T
    r = randi([1,8],8,1);
    Z = [1 + 0*1i, exp(1i*2*pi/8), exp(1i*4*pi/8), exp(1i*6*pi/8) ...
        , exp(1i*8*pi/8), exp(1i*10*pi/8), exp(1i*12*pi/8), ...
        exp(1i*14*pi/8),exp(1i*16*pi/8)];
    ZVec = [Z(r(1));Z(r(2));Z(r(3));Z(r(4));Z(r(5));Z(r(6));Z(r(7));Z(r(8))];
    mGenGly = Z(r(1))'*Z(r(2))'*Z(r(3))'*Z(r(4))'*Z(r(5))'*Z(r(6))'...
            *Z(r(7))'*Z(r(8))' ...
            *dot(randU(1,:),ZVec) ...
            *dot(randU(2,:),ZVec) ...
            *dot(randU(3,:),ZVec) ...
            *dot(randU(4,:),ZVec) ...
            *dot(randU(5,:),ZVec) ...
            *dot(randU(6,:),ZVec) ...
            *dot(randU(7,:),ZVec) ...
            *dot(randU(8,:),ZVec);   
    sum = sum + mGenGly;
    count = count + 1;
end
transAmp1 = sum/T;
fprintf(file1, 'Transition amplitude for Gurvits is %14.8e\n', abs(transAmp1));
toc2 = toc(tic2);
fprintf(file1, 'Total time for Gurvits is %14.8e\n\n', toc2);


fclose(file1);
