% Classical simulation for 10 photons and 10 modes using Glynn's algorithm
% and Gurvits's algorithm. 
% 
% Both the initial state and final state are
% standard state, i.e. [1,1,1,1,1,1,1,1,1,1]. Transition amplitudes are
% calculated. 
% 
% Total Glynn's algorithm takes 27.7 days, so only 1/10^4 fraction of it is 
% executed. If number fo samples is chosen to be 5*10^6, total Gurvits's 
% algorithm takes 17 minutes.
% 
% ©2018 Jin-Long Huang All Right Reserved
% -------------------------------------------------------------------------

file1 = fopen('test10.txt', 'w');

randU = RandomUnitary(10)

% Glynn's Algorithm (deterministic algorithm)
i0 = 1;
i1 = 1;
i2 = 1;
i3 = 1;
i4 = 1;
i5 = 1;
i6 = 1;
i7 = 1;
i8 = 1;
i9 = 1;
sum = 0+0*1i;
mGenGly = 0+0*1i;
Z = [1 + 0*1i, exp(1i*2*pi/10), exp(1i*4*pi/10), exp(1i*6*pi/10)...
    , exp(1i*8*pi/10), exp(1i*10*pi/10), exp(1i*12*pi/10), ...
    exp(1i*14*pi/10),exp(1i*16*pi/10),exp(1i*18*pi/10),];
tic0 = tic;
while i0 <= 10
    while i1 <= 10
        while i2 <= 10
            while i3 <= 10
                while i4 <= 10
                    tic1 = tic;
                    while i5 <= 10                          
                        while i6 <= 10                               
                            while i7 <= 10                                
                                while i8 <= 10
                                    while i9 <= 10
                                        ZVec = [Z(i0);Z(i1);Z(i2);...
                                                Z(i3);Z(i4);Z(i5);...
                                                Z(i6);Z(i7);Z(i8);Z(i9)];
                                        mGenGly = Z(i0)'*Z(i1)'*Z(i2)' ...
                                                *Z(i3)'*Z(i4)'*Z(i5)'...
                                                *Z(i6)'*Z(i7)'*Z(i8)'...
                                                *Z(i9)' ...
                                                *dot(randU(1,:),ZVec)...
                                                *dot(randU(2,:),ZVec) ...
                                                *dot(randU(3,:),ZVec) ...
                                                *dot(randU(4,:),ZVec) ...
                                                *dot(randU(5,:),ZVec) ...
                                                *dot(randU(6,:),ZVec) ...
                                                *dot(randU(7,:),ZVec) ...
                                                *dot(randU(8,:),ZVec) ...
                                                *dot(randU(9,:),ZVec) ...
                                                *dot(randU(10,:),ZVec);                            
                                        sum = sum + mGenGly; 
                                        i9 = i9 + 1;
                                    end
                                    i9 = 1;
                                    i8 = i8 + 1;
                                end                              
                                i8 = 1;
                                i7 = i7 + 1;
                            end                            
                            i7 = 1;
                            i6 = i6 + 1;
                        end
                        i6 = 1;                      
                        i5 = i5 + 1;
                    end
                    toc1 = toc(tic1);
                    sum1 = abs(sum)/(10^5*i4);
                    fprintf(file1, 'Transition amplitude from %d*10^5 samples is %14.8e\n', i4, sum1);
                    if(i4==1)
                        fprintf(file1, '1/10^5 time for Glynn is %14.8e \n\n', toc1);  
                    end
                    i5 = 1;
                    i4 = i4 + 1;
                end
                break;
                i4 = 1;
                i3 = i3 + 1;
            end
            break;
            i3 = 1;
            i2 = i2 + 1;
        end
        break;
        i2 = 1;            
        i1 = i1 + 1;
    end
    break;
    i1 = 1;        
    i0 = i0 + 1;          
end
transAmp = abs(sum)/(10^6);
fprintf(file1, 'Transition amplitude from 10^6 samples is %14.8e\n', transAmp);


% Gurvits' algorithm (Sampling algorithm)
tic2 = tic;
T = 5*10^6;
count = 1;
sum = 0 + 0*1i;
while count <= T
    r = randi([1,10],10,1);
    Z = [1 + 0*1i, exp(1i*2*pi/10), exp(1i*4*pi/10), exp(1i*6*pi/10) ...
        , exp(1i*8*pi/10), exp(1i*10*pi/10), exp(1i*12*pi/10), ...
        exp(1i*14*pi/10),exp(1i*16*pi/10),exp(1i*18*pi/10),];
    ZVec = [Z(r(1));Z(r(2));Z(r(3));Z(r(4));Z(r(5));Z(r(6));Z(r(7));Z(r(8));Z(r(9));Z(r(10))];
    mGenGly = Z(r(1))'*Z(r(2))'*Z(r(3))'*Z(r(4))'*Z(r(5))'*Z(r(6))'...
            *Z(r(7))'*Z(r(8))'*Z(r(9))'*Z(r(10))' ...
            *dot(randU(1,:),ZVec)...
            *dot(randU(2,:),ZVec) ...
            *dot(randU(3,:),ZVec) ...
            *dot(randU(4,:),ZVec) ...
            *dot(randU(5,:),ZVec) ...
            *dot(randU(6,:),ZVec) ...
            *dot(randU(7,:),ZVec) ...
            *dot(randU(8,:),ZVec) ...
            *dot(randU(9,:),ZVec) ...
            *dot(randU(10,:),ZVec);   
    sum = sum + mGenGly;
    count = count + 1;
end
transAmp1 = sum/T;
fprintf(file1, 'Transition amplitude for Gurvits is %14.8e\n', abs(transAmp1));
toc2 = toc(tic2);
fprintf(file1, 'Total time for Gurvits is %14.8e\n\n', toc2);


fclose(file1);