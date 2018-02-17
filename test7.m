% Classical simulation for 7 photons and 7 modes using Glynn's algorithm
% and Gurvits's algorithm. 
% 
% Both the initial state and final state are
% standard state, i.e. [1,1,1,1,1,1,1]. Transition amplitudes are
% calculated. 
% 
% Total Glynn's algorithm takes 17s. 
% If number fo samples is chosen to be 5*10^5, total Gurvits's algorithm takes 73s.
% 
% ©2018 Jin-Long Huang All Right Reserved
% -------------------------------------------------------------------------

% random unitary matrix in [0,1]
C = rand(7) + 1i * rand(7);
randH = (C + C')/2;
randU = expm(1i*randH)

file1 = fopen('test7.txt', 'w');

% Glyn's Algorithm (deterministic algorithm)
i0 = 1;
i1 = 1;
i2 = 1;
i3 = 1;
i4 = 1;
i5 = 1;
i6 = 1;
sum = 0+0*1i;
mGenGly = 0+0*1i;
Z = [1 + 0*1i, exp(1i*2*pi/7), exp(1i*4*pi/7), exp(1i*6*pi/7)...
    , exp(1i*8*pi/7), exp(1i*10*pi/7), exp(1i*12*pi/7),exp(1i*14*pi/7)];
while i0 <= 7
    stic0 = tic;
    while i1 <= 7
        stic1 = tic;
        while i2 <= 7
            while i3 <= 7
                while i4 <= 7
                    while i5 <= 7
                        while i6 <= 7
                        ZVec = [Z(i0);Z(i1);Z(i2);Z(i3);Z(i4);Z(i5);Z(i6)];
                        mGenGly = Z(i0)'*Z(i1)'*Z(i2)'*Z(i3)'*Z(i4)'...
                        *Z(i5)'*Z(i6)'*dot(randU(1,:),ZVec)...
                        *dot(randU(2,:),ZVec)...
                        *dot(randU(3,:),ZVec)*dot(randU(4,:),ZVec)...
                        *dot(randU(5,:),ZVec)*dot(randU(6,:),ZVec)...
                        *dot(randU(7,:),ZVec);                            
                        sum = sum + mGenGly;                              
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
        i1 = i1 + 1;
    end
    i1 = 1;    
    stoc0 = toc(stic0);
    fprintf(file1, 'i0   %3d   %14.8e\n', i0, stoc0);
    transAmp = abs(sum)/(i0*7^6);
    fprintf(file1, 'abs(sum) is %14.8e\n', abs(sum));
    fprintf(file1, 'TransAmp from (%2d * 7^6) samples is %14.8e\n \n ', i0, transAmp);
    i0 = i0 + 1;
end
transAmp2 = abs(sum)/(7^7);
fprintf(file1, '\nTransition amplitude for Glynn from 7^7 samples is %14.8e\n', transAmp2);
fprintf(file1, '\nTotal time for Glynn is   %14.8e\n', stoc0);


% Gurvits' algorithm (Sampling algorithm)
tic2 = tic;
T = 5*10^5;
count = 1;
sum = 0 + 0*1i;
while count <= T
    r = randi([1,7],7,1);
    Z = [1 + 0*1i, exp(1i*2*pi/7), exp(1i*4*pi/7), exp(1i*6*pi/7) ...
        , exp(1i*8*pi/7), exp(1i*10*pi/7), exp(1i*12*pi/7), ...
        exp(1i*14*pi/7)];
    ZVec = [Z(r(1));Z(r(2));Z(r(3));Z(r(4));Z(r(5));Z(r(6));Z(r(7))];
    mGenGly = Z(r(1))'*Z(r(2))'*Z(r(3))'*Z(r(4))'*Z(r(5))'*Z(r(6))'...
            *Z(r(7))' ...
            *dot(randU(1,:),ZVec) ...
            *dot(randU(2,:),ZVec) ...
            *dot(randU(3,:),ZVec) ...
            *dot(randU(4,:),ZVec) ...
            *dot(randU(5,:),ZVec) ...
            *dot(randU(6,:),ZVec) ...
            *dot(randU(7,:),ZVec);   
    sum = sum + mGenGly;
    count = count + 1;
end
transAmp1 = sum/T;
fprintf(file1, 'Transition amplitude for Gurvits is %14.8e\n', abs(transAmp1));
toc2 = toc(tic2);
fprintf(file1, 'Total time for Gurvits is %14.8e\n\n', toc2);


fclose(file1);
