clc; clear all; close all;
n=input('Enter the number Users (Signals): ');
nseconds=input('Enter the number of seconds for each User: ');
fprintf('Start recording Users masseges \n');
recording_FS=16000;
resampling_fs=8000;
recObj = audiorecorder(recording_FS,16,1); 
MUsersMsg=zeros(n,nseconds*recording_FS);   
for i= 1:n
    disp(['Press "Enter" then start speaking for ', nseconds,' seconds']);
    pause;
    recordblocking(recObj,nseconds);
    disp('The recording has ended.');
    myRecording = getaudiodata(recObj);
    sound(myRecording,recording_FS);
    MUsersMsg(i,:)=myRecording;   
end

%% Arguments of qauntization function qauntize
partition=-1:2/255:1-2/255;
par=[partition 1];
%%

%% Resampling,Quantization,MSQE audio signal for each User 
ind_mat=zeros(n,nseconds*resampling_fs);
quan_mat=zeros(n,nseconds*resampling_fs);
MSQE=zeros(1,n);
for i= 1:n
    k1(i,:)=resample(MUsersMsg(i,:),resampling_fs,recording_FS); 
    [index,quantized_sig,distor] = quantiz(k1(i,:),partition,par);
    ind_mat(i,:)=index;
    quan_mat(i,:)=quantized_sig;
    MSQE(1,i)=distor;
end
%%

%% PCM 
pcm_mat=zeros(n,nseconds*resampling_fs*8);
for i=1:n
    x1=de2bi(ind_mat(i,:),'left-msb');
    x2=transpose(x1);
    x3=reshape(x2,1,nseconds*resampling_fs*8);
    pcm_mat(i,:)=x3(1,:);
end
%%

%% TDM signal generation 
%framing of TDMsignal
fram_bit_vactor=[1 1 1 1 0 1 0 0 0 0];
nb=nseconds*resampling_fs;
fram_mat=zeros(nb,8*n+1);
temp=1;
temp1=1;
for i=1:nb
    fram_mat(i,1)=fram_bit_vactor(temp1);
    for j=1:n
        start=(temp-1)*8+1;
        ed=start+7;
        f=pcm_mat(j,start:ed);
        start2=(j-1)*8+1;
        ed2=start2+7;
        fram_mat(i,start2+1:ed2+1)=f;
    end
    temp =temp+1;
    temp1 =temp1+1;
    if (temp1==11)
        temp1=1;
    end
end
%%
%% Adding Noise
y = awgn(fram_mat,15,'measured') ;
%%

%% Detection
for i=1:nb
    for j=1:8*n+1;
        if (y(i,j)>=.5)
            y(i,j)=1;
        else
            y(i,j)=0;
        end
    end
end

%%
%% Receiving
rec_pcm=zeros(n,nseconds*resampling_fs*8);

tem=1;
for i=1:n
    temp=1;
    for j=1:nb
        start=(tem-1)*8+1;
        ed=start+7 ;
        f=y(j,start+1:ed+1);
        start2=(temp-1)*8+1;
        ed2=start2+7;
        rec_pcm(i,start2:ed2)=f;
        temp=temp+1;
    end
    tem=tem+1;
end

map=zeros(1,n);
s='User number ';
s2=' will call User number =';
for i=1:n
    MUsersMsg=[s num2str(i)];
    MUsersMsg=[MUsersMsg s2];  
    map(i)=input(MUsersMsg);
end
%%
%% Maping each couple of users
final_pcm=zeros(n,nseconds*resampling_fs*8);
for i=1:n
    final_pcm(i,:)=rec_pcm(map(i),:);
end


rec_ind_mat=zeros(n,nseconds*resampling_fs);
for i=1:n
    x=final_pcm(i,:);
    x=reshape(x,8,[]);
    x=transpose(x);
    x=bi2de(x,'left-msb');
    x=transpose(x);
    rec_ind_mat(i,:)=x;
end
    
rec_qants=rec_ind_mat*(2/255)-1;
filename1=input('Write the extension between single quotations to save the recorded files \n');
s='\user';
s3='.wav';
for i=1:n
    filename=filename1;
    filename=[filename s];
    filename=[filename num2str(i)];
    filename=[filename s3];
    audiowrite(filename,rec_qants(i,:),resampling_fs);
end

%%
%% Plot the spectrum of the TDM signal. 
% By taking sample of each user (quan_mat) we can construct the TDM signal
TDMsignal1=reshape(quan_mat,n*nseconds*resampling_fs,1);
N=length(TDMsignal1);
X = fftshift(fft(TDMsignal1));
dF = resampling_fs/N;                           
f = -resampling_fs/2:dF:resampling_fs/2-dF;    
figure;
plot(f,abs(X)/N);
xlabel('Frequency [Hz]');
title('Magnitude Response');
%%
%% The BER for each user.
    
ber=[]
for i=1:n
    b1=rec_pcm(i,:);
    b2=pcm_mat(i,:);
    [number,ratio] = biterr(b2,b1);
    ber=[ber ratio];
end
%%

%% The mean square quantization error
s='The mean square quantization error and the BER of user ';
s2=' = ';
for i=1:n
    s1=[s num2str(i)];
    s1=[s1 s2];
    disp(s1);
    disp(MSQE(1,i));
    disp(ber(i));
end
%%
