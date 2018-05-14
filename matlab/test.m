function [theta_search, Pmusic] = test(SNRdB, snapshot, M)
    %----------------------------------参数设置----------------------------------------                                                              
    f0=200;                                        %工作频率
    fs=600;                                        %采样率
    c=3e8;                                         %传播速度
    lamda=c/f0;                                 %工作波长
    d=lamda/2;                                  %阵元间距，半波长
    % M=8;                                            %阵元数
    theta=[-20 0 20]*pi/180;         %入射角度-90~90度
    K=length(theta);                         %信号个数
    % SNRdB=20;                               %信噪比
    % snapshot=200;                          %快拍数
    %方向矩阵
    A=zeros(M,K);
    for ii=1:K
        A(:,ii)=exp(-1j*2*pi*d*(0:M-1)*sin(theta(ii))/lamda).';
    end

    %入射信号
    S=zeros(K,snapshot);
    % for ii=1:K
    %     for jj=1:snapshot
    %         S(ii,jj)=(10^(SNRdB/20))*exp(j*2*pi*(f0/fs*jj+rand));%产生非相干信号
    %     end
    % end
    for ii=1:K
            S(ii,:)=(10^(SNRdB/20))*exp(j*2*pi*(f0/fs*(1:snapshot)+rand));%产生相干信号
    end

    Noise = (randn(M, snapshot) + 1i * randn(M, snapshot)) / sqrt(2);
    %均匀面阵接收信号
    X=A*S+Noise;
    %数据协方差矩阵
    R=X*X'/snapshot;%数据协方差矩阵
    %% 均匀线阵空间平滑
    m=6;                    %每个子阵阵元数
    p=M+1-m;            %划分子阵个数
    Rfk=zeros(m,m);
    Rbk=zeros(m,m);

    for k=1:p
        Zk=[zeros(m,k-1),eye(m),zeros(m,p-k)];%前向空间平滑选择矩阵
        Rfk=Rfk+Zk*R*Zk';

        Qk=[zeros(m,k-1),fliplr(eye(m)),zeros(m,p-k)];%后向空间平滑选择矩阵
        Rbk=Rbk+Qk*conj(R)*Qk';
    end
    Rf=Rfk/p;                    %前向空间平滑
    Rb=Rbk/p;                  %后向空间平滑
    Rfb=(Rf+Rb)/2;          %前后向空间平滑。使用前后向空间平滑能够增加子阵数。

    [EV,D]=eig(Rfb);
    [EVA,IBX]=sort(diag(D));%升序排列
    EV=EV(:,IBX);
    EN=EV(:,1:m-K);            %噪声子空间的提取

    theta_search=(-90:90)*pi/180;
    for ii=1:length(theta_search)
            a=exp(-1j*2*pi*d*(0:m-1)*sin(theta_search(ii))/lamda).';
            Pmusic(ii)=1/(a'*EN*EN'*a);
    end
end
