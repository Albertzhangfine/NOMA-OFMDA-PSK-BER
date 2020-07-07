clc
close 
clear 

%% 参数
order_psk = [2,4,8,16,64];      % PSK符号阶数，2表示BPSK，4表示4PSK，以此类推
for order = 1:length(order_psk)
    sym_total = 2.048e4*52; % 总符号数
    len_turbo = 1280;   % Turbo码长
    bit_total = sym_total * log2(order_psk(order));
    N_ofdma_u1 = 26;    % OFDMA用户1的子载波数
    N_ofdma_u2 = 26;    % OFDMA用户2的子载波数
    p1 = 0.1;           % NOMA用户1的功率
    p2 = 0.9;           % NOMA用户2的功率

    N_ofdm = 64;        % OFDM子载波数
    N_data = 52;        % 数据载波数
    N_GB = [4; 3];      % 载波间隔参数
    N_P = [12; 26; 40; 54]; % 导频参数
    CP = 1/4;           % CP占比

    Ts = 1/10000;
    FD = 500;           % 瑞利信道参数
    SNR = 25:1:60;           % 信噪比 

    %% 参数合理性判断
    if (mod(bit_total, len_turbo)~=0)
        error('总比特数必须是Turbo编码长度的整数倍');
    end
    if (mod((3*len_turbo+12)*bit_total/len_turbo, N_data)~=0)
        error('总比特数必须保证OFDM符号完整');
    end
    if (N_ofdma_u1+N_ofdma_u2~=N_data)
        error('OFDMA中两用户的子载波数目之和必须等于数据载波总数');
    end
    if (p1+p2~=1)
        error('功率系数p1与p2之和必须等于1');
    end

    %% 生成两个用户的发射符号序列，加入Turbo编码
    [sym_seq_u1, bit_seq_u1] = data_gen(bit_total, len_turbo, order_psk(order));
    [sym_seq_u2, bit_seq_u2] = data_gen(bit_total, len_turbo, order_psk(order));
    % mean(abs(sym_seq_u1).^2)
    % mean(abs(sym_seq_u2).^2)

    %% NOMA和OFDMA
    sym_seq_noma = noma_enc(sym_seq_u1, sym_seq_u2, p1, p2);
    % mean(abs(sym_seq_noma).^2)
    % 对两用户数据进行NOMA编码
    sym_seq_ofdma = ofdma_enc(sym_seq_u1, sym_seq_u2, N_ofdma_u1, N_ofdma_u2);
    % mean(abs(sym_seq_ofdma).^2)
    % 对两用户数据进行OFDMA编码

    %% OFDM调制
    num_ofdmsym_noma = length(sym_seq_noma)/N_data;
    mod_ofdm_noma = comm.OFDMModulator(...
    'FFTLength',N_ofdm,...
    'NumGuardBandCarriers',N_GB,...
    'PilotInputPort',true,...
    'PilotCarrierIndices',N_P,...
    'NumSymbols',num_ofdmsym_noma,...
    'CyclicPrefixLength',N_ofdm*CP,...
    'InsertDCNull',true);
    % 构造NOMA的OFDM调制器
    num_ofdmsym_ofdma = length(sym_seq_ofdma)/N_data;
    mod_ofdm_ofdma = comm.OFDMModulator(...
    'FFTLength',N_ofdm,...
    'NumGuardBandCarriers',N_GB,...
    'PilotInputPort',true,...
    'PilotCarrierIndices',N_P,...
    'NumSymbols',num_ofdmsym_ofdma,...
    'CyclicPrefixLength',N_ofdm*CP,...
    'InsertDCNull',true);
    % 构造OFDMA的OFDM调制器

    tx_noma = ofdm_tx(sym_seq_noma, mod_ofdm_noma);
    % tx_noma = sym_seq_noma;
    tx_ofdma = ofdm_tx(sym_seq_ofdma, mod_ofdm_ofdma);
    % OFDM调制
    % mean(abs(tx_noma).^2)
    % mean(abs(tx_ofdma).^2)

    for snr = 1:length(SNR)
        %% 瑞利信道
    %     crl = rayleighchan(Ts, FD);
    %     tx_noma = filter(crl, tx_noma);
    %     tx_ofdma = filter(crl, tx_ofdma);
        cawgn = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
        cawgn.SNR = SNR(snr);
        rx_noma = step(cawgn, tx_noma);
        rx_ofdma = step(cawgn, tx_ofdma);
        % mean(abs(rx_noma).^2)
        % mean(abs(rx_ofdma).^2)
        % 接收信号经过瑞利信道与加入高斯噪声

        %% OFDM解调
        sym_seq_noma_mat = ofdm_rx(rx_noma, mod_ofdm_noma);
        sym_seq_noma = reshape(sym_seq_noma_mat, numel(sym_seq_noma_mat), 1);
    %     sym_seq_noma = rx_noma;
        sym_seq_ofdma_mat = ofdm_rx(rx_ofdma, mod_ofdm_ofdma);
        sym_seq_ofdma = reshape(sym_seq_ofdma_mat, numel(sym_seq_ofdma_mat), 1);
        % mean(abs(sym_seq_noma).^2)
        % mean(abs(sym_seq_ofdma).^2)

        %% NOMA和OFDMA解码
        if (p1>p2)
            [bit_u1, bit_u2] = noma_dec(sym_seq_noma, len_turbo, order_psk(order), p1, p2);
        else
            [bit_u2, bit_u1] = noma_dec(sym_seq_noma, len_turbo, order_psk(order), p2, p1);
        end

        [n1,r1(order, snr)] = biterr(bit_u1,bit_seq_u1);
        [n2,r2(order, snr)] = biterr(bit_u2,bit_seq_u2);

        % 按照功率顺序解码用户1和用户2的符号序列
        [bit_u11, bit_u22] = ofdma_dec(sym_seq_ofdma, len_turbo, order_psk(order), N_ofdma_u1, N_ofdma_u2);

        [n11,r11(order,snr)] = biterr(bit_u11,bit_seq_u1(1:length(bit_u11)));
        [n22,r22(order,snr)] = biterr(bit_u22,bit_seq_u2(1:length(bit_u22)));
        [r1(order,snr), r2(order,snr), r11(order,snr), r22(order,snr)]
    end
    
end
%% 作图
save('ber.mat','r1','r2','r11','r22');
plot_ber(order_psk, SNR, N_data, p1, p2, N_ofdma_u1, N_ofdma_u2);

%%
function [ symbol_seq, bit_seq ] = data_gen( bit_num, len_turbo, M )
% [ symbol_seq, bit_seq ] = data_gen( bit_num, len_turbo, M )
% 输入总比特数、Turbo码长和QAM阶数，生成符号序列，同时返回原始比特序列。

bit_seq = round(rand(bit_num, 1));
% 生成随机比特序列

[turbo_seq] = turbo_enc(bit_seq, len_turbo);
% size(turbo_seq)
% Turbo编码，编码后长度为3*len+12

hModulator = comm.PSKModulator(M,'BitInput',true);
hModulator.PhaseOffset = pi/M;
symbol_seq = step(hModulator, turbo_seq);
% size(symbol_seq)
% 将比特序列转化为QAM符号并进行功率归一化
end
%%
function [ bit_u1, bit_u2 ] = noma_dec( sym_seq_noma, len, M, p1, p2 )
% [ seq_u1, seq_u2 ] = noma_dec( sym_seq_noma, M, p1, p2 )
% 按信噪比依次解码用户1序列和用户2序列
% mean(abs(sym_seq_noma).^2)

sym_seq_u1 = sym_seq_noma/sqrt(p1);
% size(sym_seq_u1)
hDemod = comm.PSKDemodulator(M, 'BitOutput', true, 'PhaseOffset',pi/M);
pskdemod_u1 = step(hDemod, sym_seq_u1);
% size(pskdemod_u1)
% 用户1PSK解调

bit_u1 = turbo_dec(pskdemod_u1, len);
% 用户1Turbo解码

turbo_seq_u1r = turbo_enc(bit_u1, len);
hMod_u1 = comm.PSKModulator(M,'BitInput',true);
hMod_u1.PhaseOffset = pi/M;
sym_seq_u1r = step(hMod_u1, turbo_seq_u1r);
% 用户1信号重构

% size(sym_seq_noma);
% size(sym_seq_u1r);

sym_seq_u2 = (sym_seq_noma - sym_seq_u1r * sqrt(p1))/sqrt(p2);
% mean(abs(sym_seq_u2).^2)
% 用户2信号
pskdemod_u2 = step(hDemod, sym_seq_u2);
% 用户2QAM解调
bit_u2 = turbo_dec(pskdemod_u2, len);
% 用户2Turbo解码

end

%%
function [ seq_noma ] = noma_enc( seq_u1, seq_u2, p1, p2 )
% [ seq_noma ] = noma_co( seq_u1, seq_u2, p1, p2 )
% 将用户1和用户2的数据序列按照p1和p2的功率分配合成NOMA序列

seq_noma = sqrt(p1) * seq_u1 + sqrt(p2) * seq_u2;
end

%%
function [ sym_seq ] = ofdm_rx( rx_data, mod)
% [ sym_seq ] = ofdm_rx( rx_data, mod )
% OFDM解调，mod为调制器

demod = comm.OFDMDemodulator(mod);
[sym_seq, ~] = step(demod,rx_data);

end

%%
function [ tx_data ] = ofdm_tx( sym_seq, mod_ofdm )
% [ tx_data ] = ofdm_tx( sym_seq, mod )
% OFDM调制，mod为调制器

N_data = mod_ofdm.FFTLength - sum(mod_ofdm.NumGuardBandCarriers) - 1 - length(mod_ofdm.PilotCarrierIndices);
N_P = length(mod_ofdm.PilotCarrierIndices);
num_ofdmsym = mod_ofdm.NumSymbols;

pilot_seq = complex(rand(N_P, num_ofdmsym),rand(N_P, num_ofdmsym));
sym_mat = reshape(sym_seq, N_data, num_ofdmsym);
% 构造OFDM符号和导频分组

tx_data = step(mod_ofdm,sym_mat,pilot_seq);

end

%%

function [ bit_u1, bit_u2 ] = ofdma_dec( sym_seq, len, M, N1, N2 )
% [ bit_u1, bit_u2 ] = ofdma_dec( sym_seq, len, M, N1, N2 )
% 根据用户1和用户2的子载波数进行OFDMA解码

sym_mat = reshape(sym_seq, N1+N2, numel(sym_seq)/(N1+N2));
% size(sym_mat)
sym_mat1 = sym_mat(1:N1,:);
sym_mat2 = sym_mat(N1+1:end,:);

sym_seq1 = reshape(sym_mat1, numel(sym_mat1), 1);
sym_seq2 = reshape(sym_mat2, numel(sym_mat2), 1);

hDemod = comm.PSKDemodulator(M, 'BitOutput', true, 'PhaseOffset',pi/M);
pskdemod_u1 = step(hDemod, sym_seq1);
bit_u1 = turbo_dec(pskdemod_u1, len);

pskdemod_u2 = step(hDemod, sym_seq2);
bit_u2 = turbo_dec(pskdemod_u2, len);
end

%%
function [ seq_ofdma ] = ofdma_enc( seq_u1, seq_u2, N_u1, N_u2 )
% [ seq_ofdma ] = ofdma_co( seq_u1, seq_u2, N_u1, N_u2, p1, p2 )
% 将用户1和用户2的数据分别按照N_u1和N_u2的子载波数合成OFDMA序列，只保留部分序列。

sym_total = length(seq_u1);
sym_u1 = sym_total * (N_u1/(N_u1+N_u2));
sym_u2 = sym_total * (N_u2/(N_u1+N_u2));
seq_u1_r = seq_u1(1:sym_u1);
seq_u2_r = seq_u2(1:sym_u2);
% 截取用户1和用户2的部分序列，保证序列总长度和NOMA的一致

mat_u1_r = reshape(seq_u1_r, N_u1, numel(seq_u1_r)/N_u1);
mat_u2_r = reshape(seq_u2_r, N_u2, numel(seq_u1_r)/N_u2);

mat_ofdma = [mat_u1_r; mat_u2_r];
seq_ofdma = reshape(mat_ofdma, numel(mat_ofdma), 1);
end

%%
function [ flag ] = plot_ber( M, SNR, N_data, p1,p2,N1,N2 )
% function [ flag ] = plot_ber( M, SNR, N_data, p1,p2,N_ofdma_u1,N_ofdma_u2 )
% 画BER的图

load('ber.mat');
SNR_ad = SNR - 10*log10(sqrt(N_data));

%% NOMA bpsk
figure;

SNR_u1_noma = SNR_ad + 10*log10(p1/log2(M(1)));
SNR_u2_noma = SNR_ad + 10*log10(p2/log2(M(1)));
SNR_u1_ofdma = SNR_ad + 10*log10((N1+N2)/N1/log2(M(2)));
SNR_u2_ofdma = SNR_ad + 10*log10((N1+N2)/N2/log2(M(2)));

semilogy(SNR_u1_noma, r1(1,:),'-ob');
hold on;
semilogy(SNR_u2_noma, r2(1,:),'-vb');
semilogy(SNR_u1_ofdma, r11(2,:), '-or');
semilogy(SNR_u2_ofdma, r22(2,:), '-vr');
box on;
grid on;

legend ('NOMA用户1','NOMA用户2','OFDMA用户1','OFDMA用户2','Location','southwest');
title('NOMA-BPSK');
xlim([10 25]);
ylim([1e-4 1]);

xlabel('Eb/N0');
ylabel('BER');

%% NOMA-4psk
figure;

SNR_u1_noma = SNR_ad + 10*log10(p1/log2(M(2)));
SNR_u2_noma = SNR_ad + 10*log10(p2/log2(M(2)));
SNR_u1_ofdma = SNR_ad + 10*log10((N1+N2)/N1/log2(M(4)));
SNR_u2_ofdma = SNR_ad + 10*log10((N1+N2)/N2/log2(M(4)));

semilogy(SNR_u1_noma, r1(2,:),'-ob');
hold on;
semilogy(SNR_u2_noma, r2(2,:),'-vb');
semilogy(SNR_u1_ofdma, r11(4,:), '-or');
semilogy(SNR_u2_ofdma, r22(4,:), '-vr');
box on;
grid on;

legend ('NOMA用户1','NOMA用户2','OFDMA用户1','OFDMA用户2','Location','southwest');
title('NOMA-4PSK');
xlim([10 32]);
ylim([1e-4 1]);

xlabel('Eb/N0');
ylabel('BER');

%% NOMA 8psk
figure;

SNR_u1_noma = SNR_ad + 10*log10(p1/log2(M(3)));
SNR_u2_noma = SNR_ad + 10*log10(p2/log2(M(3)));
SNR_u1_ofdma = SNR_ad + 10*log10((N1+N2)/N1/log2(M(5)));
SNR_u2_ofdma = SNR_ad + 10*log10((N1+N2)/N2/log2(M(5)));

semilogy(SNR_u1_noma, r1(3,:),'-ob');
hold on;
semilogy(SNR_u2_noma, r2(3,:),'-vb');
semilogy(SNR_u1_ofdma, r11(5,:), '-or');
semilogy(SNR_u2_ofdma, r22(5,:), '-vr');
box on;
grid on;

legend ('NOMA用户1','NOMA用户2','OFDMA用户1','OFDMA用户2','Location','southwest');
title('NOMA-8PSK');
xlim([15 45]);
ylim([1e-4 1]);

xlabel('Eb/N0');
ylabel('BER');

end

%%

function [ data_seq ] = turbo_dec( turbo_seq, len )
% [ data_seq ] = turbo_dec( turbo_seq, len )
% Turbo解码

len_p = 3*len+12;
turbo_mat = reshape(turbo_seq, len_p, numel(turbo_seq)/len_p);
data_seq = [];
for i = 1:numel(turbo_seq)/len_p
    data_seq0 = lteTurboDecode(turbo_mat(:,i));
    data_seq = [data_seq; data_seq0];
end

end


%%
function [ turbo_seq ] = turbo_enc( bit_seq, len )
% [ turbo_seq ] = turbo_enc( bit_seq, len )
% 对输入比特序列进行Turbo编码

bit_mat = reshape(bit_seq, len, numel(bit_seq)/len);
turbo_seq = [];
for i = 1:numel(bit_seq)/len
    turbo_seq0 = lteTurboEncode(bit_mat(:,i));
    turbo_seq = [turbo_seq; turbo_seq0];
end
end

