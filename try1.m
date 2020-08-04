bps=4; %bits per signal
M=2^bps;
nBits=1e3*bps;
ebno=10; %SNR

%constellation
symMap       = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
bitTable     = de2bi(symMap,bps,'left-msb');

%symbolset modulation and reshaping
sym          = qammod(symMap(1:M),M,symMap,'UnitAveragePower',false,'PlotConstellation',true); 
sym          = reshape(sym,[M,1]);
averagePower = sum(abs(sym).^2)/16;


%awgnChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',10,'SignalPower',averagePower);

%built-in error determination function
berRate  = comm.ErrorRate;

%random input datastream
data            = randi([0 1],nBits,1);
%modulation
modData         = qammod(data,M,symMap,'InputType','bit','UnitAveragePower',false,'PlotConstellation',true);


h11_r = normrnd(0,(1/sqrt(2)));
h11_i = normrnd(0,(1/sqrt(2)));
h12_r = normrnd(0,(1/sqrt(2)));
h12_i = normrnd(0,(1/sqrt(2)));
h21_r = normrnd(0,(1/sqrt(2)));
h21_i = normrnd(0,(1/sqrt(2)));
h22_r = normrnd(0,(1/sqrt(2)));
h22_i = normrnd(0,(1/sqrt(2)));
h11 = h11_r + 1j*h11_i;
h12 = h12_r + 1j*h12_i;
h21 = h21_r + 1j*h21_i;
h22 = h22_r + 1j*h22_i;

H = [[h11,h12];[h21,h22]]; %matrix H
len=nBits/bps; %Number of signals
x = reshape(modData,[2,len/2]);
y1 = H*x;
y1 = y1(:);

rxsig           = step(awgnChan,y1);
rxsig_real      = real(rxsig);
rxsig_img       = imag(rxsig);


for i=1:len
    rxsig_2(2*i-1)=rxsig_real(i);
    rxsig_2(2*i)=rxsig_img(i);
end

rxsig_2=reshape(rxsig_2,[4,len/2]);

%mapping and manipulation

B = [[h11_r,-h11_i,h12_r,-h12_i];[h11_i,h11_r,h12_i,h12_r];[h21_r,-h21_i,h22_r,-h22_i];[h21_i,h21_r,h22_i,h22_r]];
shift = zeros(4,1)+3;

y_c = rxsig_2 + B*shift;
[Q,R] = qr(2*B);
D = diag(sign(diag(R)));
Qunique = Q*D; 
Runique = D*R;
y_dash = Qunique'*y_c;

decodedData = zeros(4,len/2);

for z = 1:(len/2)
    x_mat = sphere_dec(1/2,Runique,y_dash(:,z));
    tempmat = x_mat-y_dash(:,z);  
    dist = vecnorm(tempmat);
    [M1,I] = min(dist);
    decodedData(:,z) = 2*x_mat(:,I) - shift;
   
end

decodedData = reshape(decodedData,[2,len])';
decodedData = decodedData(:,1)+1j*decodedData(:,2);

%demodulation
rxData = qamdemod(decodedData,M,symMap,'OutputType','bit','UnitAveragePower',false);

errorStats = step(berRate,data,rxData);
errorStats;