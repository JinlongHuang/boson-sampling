% Classical simulation for 6 photons and 6 modes using Glynn's algorithm
% and Gurvits's algorithm. 
% 
% Both the initial state and final state are
% standard state, i.e. [1,1,1,1,1,1]. Transition amplitudes are
% calculated. 
% 
% Total Glynn's algorithm takes 7s. If number of samples is 
% chosen to be 5*10^2, total Gurvits's algorithm takes 0.07s.
% If number of samples is chosen to be 5*10^3, total Gurvits's algorithm takes 0.7s.
% 
% ©2018 Jin-Long Huang All Right Reserved
% -------------------------------------------------------------------------


file1 = fopen('test6.txt', 'w');

% Glyn's Algorithm (deterministic algorithm)
randU = RandomUnitary(6)

i0 = 1;
i1 = 1;
i2 = 1;
i3 = 1;
i4 = 1;
i5 = 1;
sum = 0+0*1i;
mGenGly = 0+0*1i;
Z = [1 + 0*1i, exp(1i*2*pi/6), exp(1i*4*pi/6), exp(1i*6*pi/6)...
    , exp(1i*8*pi/6), exp(1i*10*pi/6)];
tic0 = tic;
while i0 <= 6
    tic1 = tic;
    while i1 <= 6
%             stic1 = tic;
        while i2 <= 6
            while i3 <= 6
                while i4 <= 6
                    while i5 <= 6
                        ZVec = [Z(i0);Z(i1);Z(i2);Z(i3);Z(i4);Z(i5)];
                        mGenGly = Z(i0)'*Z(i1)'*Z(i2)' ...
                            *(Z(i3)')*Z(i4)'*Z(i5)'...
                            *dot(randU(1,:),ZVec)*dot(randU(2,:),ZVec) ...
                            *dot(randU(3,:),ZVec)*dot(randU(4,:),ZVec) ...
                            *dot(randU(5,:),ZVec)*dot(randU(6,:),ZVec);                            
                        sum = sum + mGenGly; 
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
    toc1 = toc(tic1);
    transAmp = abs(sum)/(i0*6^5);
    fprintf(file1, '\n 1/6 Time for Glynn is   %14.8e\n', toc1);
    fprintf(file1, 'Transition amplitude from (%2d * 6^5) samples is %14.8e\n', i0, transAmp);        
    i0 = i0 + 1;
%     break;
end
transAmp2 = abs(sum)/(6^6);
toc0 = toc(tic0);
fprintf(file1, 'Total time for Glynn is   %14.8e\n', toc0);
fprintf(file1, 'Transition amplitude for Glynn is %14.8e\n', transAmp2);


% Gurvits' algorithm (Sampling algorithm)
% Parameters
TOL0 = 0.050;
TOL2 = 0.005;
STEP = 5;
maxExp = 1000;

% T = 5*10^2
TOL1 = TOL0;
FID0 = 1 - TOL0;
FID2 = 1 - TOL2;
INTE = (TOL1-TOL2)/(STEP-1);
numExp = 1;
suc = 0;
datay = zeros(STEP,1);
roww = 1;
T = 5*10^2;

while TOL1 >= TOL2-0.0001
    tic2 = tic;   
    while numExp <= maxExp   
        count = 1;
        sum = 0 + 0*1i;
        while count <= T
            r = randi([1,6],6,1);
            Z = [1 + 0*1i, exp(1i*2*pi/6), exp(1i*4*pi/6), exp(1i*6*pi/6) ...
                , exp(1i*8*pi/6), exp(1i*10*pi/6)];
            ZVec = [Z(r(1));Z(r(2));Z(r(3));Z(r(4));Z(r(5));Z(r(6))];
            mGenGly = Z(r(1))'*Z(r(2))'*Z(r(3))' ...
                *Z(r(4))'*Z(r(5))'*Z(r(6))'...
                *dot(randU(1,:),ZVec)*dot(randU(2,:),ZVec) ...
                *dot(randU(3,:),ZVec)*dot(randU(4,:),ZVec) ...
                *dot(randU(5,:),ZVec)*dot(randU(6,:),ZVec);   
            sum = sum + mGenGly;
            count = count + 1;
        end
        transAmp1 = sum/T;
%         fprintf(file1, '\ntransAmp1 is %14.8e\n', abs(transAmp1));

%         dif = abs(transAmp1-transAmp2)/abs(transAmp2);
        if abs(transAmp1-transAmp2) <= TOL1
            suc = suc + 1;
        end 
        pro = suc/maxExp;            
%         fprintf(file1, '\nTime for Gurvits is   %14.8e\n', toc2);
        numExp = numExp + 1;
    end
    toc2 = toc(tic2)
    fprintf(file1, '\n Time is   %14.8e\n', toc2);
%     break;
    FID1 = 1-TOL1;
    fprintf(file1,'%4.3f \t %4.3f \n', FID1, pro);
    datay(roww) = pro;
    roww = roww + 1;
    numExp = 1;
    suc = 0;
    TOL1 = TOL1 - INTE;
end

% Gurvits' algorithm (Sampling algorithm)
T = 5*10^3
TOL1 = TOL0;
INTE = (TOL1-TOL2)/(STEP-1);
numExp = 1;
suc = 0;
datay2 = zeros(STEP,1);
roww = 1;

while TOL1 >= TOL2-0.0001
    tic2 = tic;   
    while numExp <= maxExp   
        T = 5*10^3;
        count = 1;
        sum = 0 + 0*1i;
        while count <= T
            r = randi([1,6],6,1);
            Z = [1 + 0*1i, exp(1i*2*pi/6), exp(1i*4*pi/6), exp(1i*6*pi/6) ...
                , exp(1i*8*pi/6), exp(1i*10*pi/6)];
            ZVec = [Z(r(1));Z(r(2));Z(r(3));Z(r(4));Z(r(5));Z(r(6))];
            mGenGly = Z(r(1))'*Z(r(2))'*Z(r(3))' ...
                *Z(r(4))'*Z(r(5))'*Z(r(6))'...
                *dot(randU(1,:),ZVec)*dot(randU(2,:),ZVec) ...
                *dot(randU(3,:),ZVec)*dot(randU(4,:),ZVec) ...
                *dot(randU(5,:),ZVec)*dot(randU(6,:),ZVec);   
            sum = sum + mGenGly;
            count = count + 1;
        end
        transAmp1 = sum/T;
%         fprintf(file1, '\ntransAmp1 is %14.8e\n', abs(transAmp1));

%         dif = abs(transAmp1-transAmp2)/abs(transAmp2);
        if abs(transAmp1-transAmp2) <= TOL1
            suc = suc + 1;
        end 
        pro = suc/maxExp;            
%         fprintf(file1, '\nTime for Gurvits is   %14.8e\n', toc2);
        numExp = numExp + 1;
    end
    toc2 = toc(tic2)
    fprintf(file1, '\n Time is   %14.8e\n', toc2);
%     break;
    FID1 = 1-TOL1;
    fprintf(file1,'%4.3f \t %4.3f \n', FID1, pro);
    datay2(roww) = pro;
    roww = roww + 1;
    numExp = 1;
    suc = 0;
    TOL1 = TOL1 - INTE;
end


% Theoretical bound
xthe = linspace(0.95,0.995,20);
ythe = zeros(20,1);
count = 1;
for i = 0.05:-0.0025:0.005-0.0025
    ythe(count) = 1 - 4*exp(-(5*10^3*i^2)/4);
    count = count + 1;
end


x = linspace(FID0, FID2, STEP);
plot(x,datay, '-*' , x, datay2, '-^', xthe, ythe,'--','LineWidth', 1.5); 
axis([FID0 FID2 0.0 1.0]);
ax = gca;
ax.FontSize = 18;
xlabel('Fidelity','FontSize',18);
ylabel('Success Probability','FontSize',18);
legend({'T1 = 5x10^2','T2 = 5x10^3','Chernoff Bound for T2'},'FontSize',18);

fclose(file1);