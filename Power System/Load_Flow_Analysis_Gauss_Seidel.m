%% Load Flow Analysis using Gauss-Seidel Method

a = [0 0.02 0.01;
    0.02 0 0.0125;
    0.01 0.0125 0];

b = [0 0.04 0.03;
    0.04 0 0.025;
    0.03 0.025 0;];

LineData = complex(a,b); 

Ybus = calculateYbus(LineData , 3);

disp('Ybus=');
disp(Ybus);

BusTypes = [1 2 3;]; 
Base_mva = 100;

Genload = [0 0 0 0;
    50 30 305.6 140.2;
    0.0 0.0 138.6 45.2;];

P = zeros(size(Ybus,1),1);
Q = zeros(size(Ybus,1),1);

for i = 2:3
    P(i) = (Genload(i,1)- Genload(i,3))/Base_mva;
    Q(i) = (Genload(i,2)- Genload(i,4))/Base_mva;
end

PQ = complex(P,Q);

V = [1.05+0i;
    1+0i;
    1+0i;];

V_old = V;


%% (a) phasor values of the voltage at buses 2 and 3
for i = 1:2
    for k = 2:3
        temp = [0;
            0;
            0;];
 
        for n = 1:3
            if k ~= n
                temp(k,1) = temp(k,1) +(Ybus(k,n) * V(n,1)); 
            end
        
        end
        V(k,1) = (1/Ybus(k,k))*((conj(PQ(k,1))/V(k,1)) - temp(k,1));
    end
end
disp("V=");
disp(V);

r2 = abs(V(2,1));
theta2 = rad2deg(angle(V(2,1)));
r3 = abs(V(3,1));
theta3 = rad2deg(angle(V(3,1)));
disp("(a) After 2 iterations");
fprintf("V2 = %f<%f\n",r2,theta2);
fprintf("V3 = %f<%f\n\n\n",r3,theta3);


%% (b)Slack bus real and reactive power after second iteration
SlackBusIndex = 1;
temp = 0;
for i = 1:3
    temp = temp + (Ybus(SlackBusIndex,i)* V(i,1)* conj(V(SlackBusIndex,1)));
end

disp('(b) Slack bus real and reactive power after second iteration:');
fprintf('\nP_slack = %f MW\n', (real(temp))*100);
fprintf('Q_slack = %f MVAR\n\n\n', -(imag(temp))*100);

%% (c)line flows and line losses after second iteration.
I = zeros(3,3);
I(1,2) = Ybus(1,2)*(V(1,1)-V(2,1));
I(2,1) = -I(1,2);
I(1,3) = Ybus(1,3)*(V(1,1)-V(3,1));
I(3,1) = -I(1,3);
I(3,2) = Ybus(3,2)*(V(3,1)-V(2,1));
I(2,3) = -I(3,2);

S = zeros(3,3);
for i= 1:3
    for j = 1:3
        if i~=j
            S(i,j) = V(i,1)* conj(I(i,j));
        end
    end
end

disp("(c) Flow and Line Losses after 2nd iteration S = ");
disp(S);

PLoss12 = real(S(1,2)) + real(S(2,1));
PLoss13 = real(S(1,3)) + real(S(3,1));
PLoss23 = real(S(2,3)) + real(S(3,2));

fprintf('\nPLoss12 = %f MW\n', (PLoss12*100));
fprintf('PLoss13 = %f MW\n', (PLoss13*100));
fprintf('PLoss23 = %f MW\n\n', (PLoss23*100));


QLoss12 = imag(S(1,2)) + imag(S(2,1));
QLoss13 = imag(S(1,3)) + imag(S(3,1));
QLoss23 = imag(S(2,3)) + imag(S(3,2));

fprintf('\nQLoss12 = %f MVAR\n', (QLoss12*100));
fprintf('QLoss13 = %f MVAR\n', (QLoss13*100));
fprintf('QLoss23 = %f MVAR\n\n\n', (QLoss23*100));

%% function to calculate ybus
function Ybus = calculateYbus(LineData, numBuses)

    Ybus = zeros(numBuses, numBuses);
    s = size(Ybus, 1);
    for i = 1:s
        for j = 1: size(Ybus,2)
            data = 0;
            if i ~=j
                data = -(1/LineData(i,j));
            else
                for k = 1:3
                    if i ~= k
                        data = data + 1/(LineData(i,k));
                    end
                end

            end
                
            Ybus(i, j) = data;
        end
    end
end
