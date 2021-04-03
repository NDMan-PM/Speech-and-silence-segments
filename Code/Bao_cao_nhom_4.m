clear
%====================PHAN CHUONG TRINH CHINH===============================

% Doc tin hieu tu file audio
% x: mang chua cac gia tri bien do
% Fs: Tan so lay mau
 [x,Fs] = audioread('lab_female.wav');
% [x,Fs] = audioread('lab_male.wav');
% [x,Fs] = audioread('studio_male.wav');
% [x,Fs] = audioread('studio_female.wav');
fr_time = 0.02;         %fr_time(s): Do dai moi khung hinh
fr_len = fr_time*Fs;    %fr_len: So mau tren moi khung hinh
N = length(x);          %N: tong so mau cua tin hieu
% Ham tinh tong so khung hinh(frame)
fr_num = NumberOfFrames(N,fr_len);
% Ham tinh nang luong cua moi khung hinh
Power = SumPower(x, fr_num, fr_len);
min_P = min(Power);
max_P = max(Power);

%==================== Chuong trinh thuat toan ham STE ====================

% Ham chuan hoa ve nguong [0;1]
Standard = Standard_function(Power, fr_num, min_P, max_P);
% Ham xac dinh nguong thuat toan STE
Thres_STE = Ethres_function(Power,10);
% Ham phat hien tieng noi va khoang lang
Voiced_STE = DetecVoiced_function(Standard, fr_num, Thres_STE);
% Ham xac dinh vi tri bien dau-cuoi cua tieng noi
Position_STE = Position_function(Voiced_STE, fr_num, fr_len);
% Ham truoc khi xu ly nhung doan khoang lang < 200ms
Pre_Voiced_STE = Pre_DetecVoiced_function(Power, fr_num, Thres_STE);
% Ham xac dinh vi tri bien dau-cuoi cua tieng noi
Pre_Position_STE = Position_function(Pre_Voiced_STE, fr_num, fr_len);

%=========== Chuong trinh thuat toan ham ZCR + STE ========================

% Chuan hoa lai tin hieu de tim ham zcr
X_ZCR = X_ZCR_function(x);
% Dai luong can thiet de tinh ham zcr (-1,1)
sgn = sgn_function(X_ZCR, fr_num, fr_len);
% Ham zcr
ZCR = ZCR_function(sgn, fr_num, fr_len);
% Ham thuat toan ket hop giua ZCR va STE
W = W_function(Power, ZCR, fr_len, fr_num);
% Ham xac dinh nguong thuat toan ZCR ket hop STE
Thres_ZS = Ethres_function(W,30);
% Ham phat hien tieng noi va khoang lang
Voiced_ZS = DetecVoiced_function(W, fr_num, Thres_ZS);
% Ham xac dinh vi tri bien dau-cuoi cua tieng noi
Position_ZS = Position_function(Voiced_ZS, fr_num, fr_len);
% Ham truoc khi xu ly nhung doan khoang lang < 200ms
Pre_Voiced_ZS = Pre_DetecVoiced_function(W, fr_num, Thres_ZS);
% Ham xac dinh vi tri bien dau-cuoi cua tieng noi
Pre_Position_ZS = Position_function(Pre_Voiced_ZS, fr_num, fr_len);

%=============== Chuong trinh thuat toan MA ===============================

MA = MA_function(x, fr_num, fr_len);
min_MA = min(MA);
max_MA = max(MA);
%Ham chuan hoa ve nguong [0;1]
Standard_MA = Standard_function(MA, fr_num, min_MA, max_MA);
%Ham xac dinh nguong thuat toan MA
Thres_MA = Ethres_function(Standard_MA,10);
%Ham phat hien tieng noi va khoang lang
Voiced_MA = DetecVoiced_function(Standard_MA, fr_num, Thres_MA);
%Ham xac dinh vi tri bien dau - cuoi cua tieng noi
Position_MA = Position_function(Voiced_MA, fr_num, fr_len);
% Ham truoc khi xu ly nhung doan khoang lang < 200ms
Pre_Voiced_MA = Pre_DetecVoiced_function(Standard_MA, fr_num, Thres_MA);
% Ham xac dinh vi tri bien dau-cuoi cua tieng noi
Pre_Position_MA = Position_function(Pre_Voiced_MA, fr_num, fr_len);

% Bien tieng noi- khoang lang khi xac dinh thu cong cua 4 file
Lab_female = [11066, 19634, 31334, 39429, 55010, 61614, 74171, 82260,  97346, 103315, 113887, 125102];
Lab_male = [14519, 21359, 37435, 43339, 65582, 70570, 107175, 112666, 128879, 134030, 152803, 159493];
studio_male = [11201, 47415, 52431, 71211, 76935, 91042, 96296, 116301];
studio_female = [8996, 34395, 39570, 68135];

subplot(411), plot(x);
title('Do thi phan doan thu cong file Labfemale.wav');
hold on;
for i =1:length(Lab_female)/2
    xline(Lab_female(2*i-1), 'r','LineWidth',1);
    xline(Lab_female(2*i), 'g','LineWidth',1);
end
legend('Tin hieu mau','Vi tri bat dau co tieng noi','Vi tri ket thuc tieng noi','location','northeast');
hold off;

subplot(412), plot(Standard);
title('Do thi ham nang luong ngan han (STE)');
hold on;
for i =1:length(Pre_Position_STE)/2
    xline(Pre_Position_STE(2*i-1)/fr_len, 'r','LineWidth',1);
    xline(Pre_Position_STE(2*i)/fr_len, 'g','LineWidth',1);
end
yline(Thres_STE,'-r','Thres');
hold off;


subplot(413), plot(x);
title('Phan doan tieng noi va khoang lang (STE)');
hold on;
for i =1:length(Position_STE)/2
    xline(Position_STE(2*i-1), 'r','LineWidth',1);
    xline(Position_STE(2*i), 'g','LineWidth',1);
end
hold off;
legend('Tin hieu mau','Vi tri bat dau co tieng noi','Vi tri ket thuc tieng noi','location','northeast');

figure;
subplot(411);
plot(W);
title('Do thi ham W = P*(1-ZCR)*1000 (ZCR + STE)');
hold on;
for i =1:length(Pre_Position_ZS)/2
    xline(Pre_Position_ZS(2*i-1)/fr_len, 'r','LineWidth',1);
    xline(Pre_Position_ZS(2*i)/fr_len, 'g','LineWidth',1);
end
yline(Thres_ZS,'-r','Thres');
hold off;
subplot(412), plot(x);
title('Phan doan tieng noi va khoang lang (ZCR + STE)');
hold on;
for i =1:length(Position_ZS)/2
    xline(Position_ZS(2*i-1), 'r','LineWidth',1);
    xline(Position_ZS(2*i), 'g','LineWidth',1);
end
hold off;
legend('Tin hieu mau','Vi tri bat dau co tieng noi','Vi tri ket thuc tiengnoi','location','northeast');
subplot(413);
plot(W);
title('Do thi ham bien do trung binh (MA)');
hold on;
for i =1:length(Pre_Position_MA)/2
    xline(Pre_Position_MA(2*i-1)/fr_len, 'r','LineWidth',1);
    xline(Pre_Position_MA(2*i)/fr_len, 'g','LineWidth',1);
end
yline(Thres_ZS,'-r','Thres');
hold off;
subplot(414), plot(x);
title('Phan doan tieng noi va khoang lang (MA)');
hold on;
for i =1:length(Position_MA)/2
    xline(Position_MA(2*i-1), 'r','LineWidth',1);
    xline(Position_MA(2*i), 'g','LineWidth',1);
end
hold off;
legend('Tin hieu mau','Vi tri bat dau co tieng noi','Vi tri ket thuc tiengnoi','location','northeast');

% % Xac dinh sai so giua xac dinh thu cong va thuat toan tren file lab_female
% RMSE_STE_lab_female = sqrt(mean((Position_STE - Lab_female).^2))
% RMSE_ZS_lab_female = sqrt(mean((Position_ZS - Lab_female).^2))
% RMSE_MA_lab_female = sqrt(mean((Position_MA - Lab_female).^2))
 
% % Xac dinh sai so giua xac dinh thu cong va thuat toan tren file lab_male
% RMSE_STE_lab_male = sqrt(mean((Position_STE - Lab_male).^2))
% RMSE_ZS_lab_male = sqrt(mean((Position_ZS - Lab_male).^2))
% RMSE_MA_lab_male = sqrt(mean((Position_MA - Lab_male).^2))
 
% % Xac dinh sai so giua xac dinh thu cong va thuat toan tren file studio_male
% RMSE_STE_studio_male = sqrt(mean((Position_STE - studio_male).^2))
% RMSE_ZS_studio_male = sqrt(mean((Position_ZS - studio_male).^2))
% RMSE_MA_studio_male = sqrt(mean((Position_MA - studio_male).^2))

% % Xac dinh sai so giua xac dinh thu cong va thuat toan tren file studio_male
% RMSE_STE_studio_female = sqrt(mean((Position_STE - studio_female).^2))
% RMSE_ZS_studio_female = sqrt(mean((Position_ZS - studio_female).^2))
% RMSE_MA_studio_female = sqrt(mean((Position_MA - studio_female).^2))

%============================PHAN DINH NGHIA HAM===========================
% Ham tinh so khung hinh(frame) cua tin hieu
% Ham tra ve tong so frame cua tin hieu
% N: Tong so mau cua tin hieu
% fr_len: so mau cua tin hieu tren mot frame
% su dung ham floor de lay nguyen phep chia
function fr_num = NumberOfFrames(N,fr_len)
    fr_num = floor(N/fr_len);          
end 

% Ham tinh nang luong tren moi khung
% Ham tra ve: mang nang luong cua moi khung hinh
% x: mang gia tri bien do cua tin hieu am thanh
% fr_num: tong so frame cua tin hieu
% fr_len: so mau tren mot frame
function Power = SumPower(x, fr_num, fr_len)
    Power = zeros(1,fr_num); % Tao mang chua Power
    for k = 1:fr_num        % Duyet tat ca cac khung
        tempPower = 0; 
        for j =(k-1)*fr_len +1 : (fr_len*k -1)  % duyet cac mau co trong khung
            tempPower = tempPower + abs(x(j)^2); % tinh nang luong cua tung khung
        end
        Power(k) = tempPower;  % Luu vao mang Power chua cac gia tri nang luong 
    end
end


% Ham chuan hoa vector Power ve [0;1]
% fr_num: tong so frame cua tin hieu
% min, max la gia tri min va max cua vector Power
function Standard = Standard_function(Power, fr_num, min, max)
    Standard = zeros(1,fr_num);     % Tao mang chua cac gia tri sau khi chuan hoa ve nguong 0,1
    for k = 1:fr_num                % Duyet tat ca cac khung
        %cong thuc chuan hoa 
        temp = (Power(k)-min)/(max-min);    % Chuan hoa bien do ve 0-1
        Standard(k) = temp;         % luu du lieu vao mang Standard sau khi chuan hoa
    end
end

% Ham dua tin hieu ve [0;1]
% Ham tra ve mot vector chua 0,1
% Voi 0 duoc hieu la khoang lang
%     1 duoc hieu la tieng noi
% Ham tu xac dinh khoang lang, neu phat hien khoang lang 
% co thoi gian < 200ms thi bien doi thanh tieng noi
function [Voiced] = DetecVoiced_function(Standard, fr_num, Thres)
    Voiced = zeros(1,fr_num); 
    for k = 1:fr_num                % Duyet het tat cac khung
        if (Standard(k) > Thres)    % So sanh voi nguong phan tach tieng noi khoang lang
            Voiced(k) = 1;          % Gan Voiced(i) = 1 neu khung lon hon gia tri nguong
        end
    end
    for k = 1:fr_num                % Duyet het tat cac khung
        if (Standard(k)< Thres)     % So sanh voi nguong phan tach tieng noi khoang lang
        Voiced(k) = 0;              % Gan Voiced(i) = 0 neu khung be hon gia tri nguong
        end 
    end
    for k = 1: fr_num                % Duyet het tat cac khung
        if (Voiced(k) == 1)          % Kiem tra neu khung thu k co gia tri bang 1 ?
            for i = k:k+9            % Duyet tiep them 9 khung ke tu khung thu k
                if Voiced(i) == 1    % Kiem tra neu cac khung tiep theo co gia tri bang 1
                    for j=k:i        % Duyet cac khung tu k den k +9
                        Voiced(j) = 1; % Gan gia tri cac khung tu k den k +9 = 1;
                    end
                end
            end        
        end
    end
end

% Ham xac dinh vi tri bien cua tieng noi-khoang lang
% Ham tra ve vi tri bien dau, bien cuoi cua tieng noi
% duoc luu vao vector Position
function [Position] = Position_function(Voiced, fr_num, fr_len)
    Position = [];
    j = 1;  
    for i = 2:fr_num      % Duyet tat ca cac khung
        if(Voiced(i) == 1 && Voiced(i-1)==0)
            Position(2*j - 1)=  fr_len/2 +(i-1)*fr_len; % Luu vi tri bien dau cua tieng noi 
            j = j + 1;
        end 
    end 
    j = 1;
    for i = 2:fr_num   
        if (Voiced(i) == 0 && Voiced(i-1) == 1)
            Position(2*j) = fr_len/2 + (i-1)*fr_len;  % Luu vi tri bien cuoi cua tieng noi
            j = j + 1;
        end
    end 
end 

% Ham xac dinh nguong phan tach tieng noi - khoang lang
% Ham tra ve gia tri cua nguong phan tach tieng noi - khoang lang
function [Ethres] = Ethres_function(Power,S)
 avg = 0;
 for i= 1:10                %
    avg = avg + Power(i)/10; %Tim gia tri trung binh cua khung trong 10 khung dau tien
 end                        % 
 vari = 0;
 for i= 1:10
    vari = vari + ((Power(i) - (vari))^2)/10;  % Tim phuong sai cua 10 khung dau tien
 end
    Ethres = S*(avg + vari*0.2*vari^0.8);  % Ap dung cong thuc de tim nguong
end

% Ham chuyen doi dai luong can thiet de tim ham zcr
% Ham tra ve ma tran cac tin hieu (-1,1)
function [sgn] = sgn_function(X_ZCR, fr_num, fr_len)
    sgn = zeros(fr_num,fr_len);
    k = 1;
    for i = 1:fr_num        % Duyet tat ca cac khung
        for j = 1:fr_len    % Duyet tat ca cac mau co trong khung
            x(i,j) = X_ZCR(k); 
            if(X_ZCR(k)>0)     % Kiem tra xem khung thu k cua ZCR lon hon 0
                sgn(i,j) = 1;   % Gan cho mang sgn gia tri 1
            else
                sgn(i,j) = -1; % Gan cho mang sgn gia tri -1
            end
            k = k + 1 ;
        end
    end
end
 
% Ham tim ZCR
% Ham tra ve mang cac gia tri ZCR cua moi khung
function [ZCR] = ZCR_function(sgn,fr_num,fr_len)
ZCR = zeros(1,fr_num);
    for i = 1 : fr_num
        ZCR(1) = abs(sgn(1,1))/2;
        for j = 2:fr_len
            zcr = abs(sgn(i,j)- sgn(i,j-1))/2;
            ZCR(i) = ZCR(i) + zcr;
        end
    end
end

% Ham chuan hoa lai tin hieu ban dau de tim ham zcr
% Ham tra ve tin hieu sau khi da chuan hoa
function [X_ZCR] = X_ZCR_function(x)
avg = 0;
for i = 1:length(x)
    avg = avg + x(i);       % Tim tong tat cac cac gia tri bien do cua tin hieu ban dau
end
avg = avg/length(x);        % Lay trung binh cong gan vao bien avg
x1 = zeros(1,length(x));
for i = 1:length(x)
    x1(i) = x(i) - avg;      % Dich chuyen vi tri cua x avg don vi
end
X_ZCR = x1;                  
end

% Ham xac dinh W ket hop giua ham STE va ZCR
% Ham tra ve mang gia tri cua W khi ket hop hai ham STE va ZCR
function [W] = W_function(Power, ZCR, fr_len, fr_num)
    P = Power/(fr_len);
    ZCR = ZCR/(fr_len);
    W = zeros(1,fr_num);
    for i = 1:fr_num
        W(i) = P(i)*(1-ZCR(i))*1000;
    end
end


% Ham truoc khi xu ly nhung doan khoang lang < 200ms
function [Pre_Voiced] = Pre_DetecVoiced_function(Standard, fr_num, Thres)
    Pre_Voiced = zeros(1,fr_num);
    for k = 1:fr_num
        if (Standard(k) > Thres)
            Pre_Voiced(k) = 1;
        end
    end
    for k = 1:fr_num
        if (Standard(k)< Thres)
        Pre_Voiced(k) = 0;
        end 
    end
end
% Ham tinh trung binh do lon tren moi khung
% Tra ve: Gia tri trung binh tren moi khung
% x: mang gia tri bien do cua tin hieu am thanh
% fr_num: tong so frame cua tin hieu
% fr_len: so mau tren mot frame
function [MA] = MA_function(x, fr_num, fr_len)
    MA = zeros(1, fr_num); % Tao mang de chua cac gia tri trung binh
    for k = 1:fr_num    % Duyet tat ca cac khung 
        temp = 0;       % Bien tam dung luu gia tri trung binh cua moi khung khi chay vong lap
        % Duyet tat ca cac mau cua tin hieu trong 1 khung
        for j = (k-1)*fr_len + 1 : k*fr_len - 1
            temp = temp + abs(x(j));
        end
        MA(k) = temp; % Gia tri trung binh tren moi khung duoc dua vao mang MA
    end
end





