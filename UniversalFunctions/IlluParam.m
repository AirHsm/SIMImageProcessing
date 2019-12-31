% IlluParam
function [Phase, K, MfList] = IlluParam(Data,OTF,WF)
  [M,N,H] = size(Data);

  K = zeros(H,2);
  Phase = zeros(H,1);

  %hwait = waitbar(0, 'Estimating frequency and phase...');
  parfor i = 1:H
    img = jRL_deconvolve(Data(:,:,i), ones(M,N), OTF, 20);
    K(i,:) = get_frequency(img, OTF);
    Phase(i) = get_phase(img, K(i,:));
    %K(i,:) = get_frequency(Data(:,:,i), OTF);
    %Phase(i) = get_phase(Data(:,:,i), K(i,:));
    %str=['estimating frequencies and phases...',num2str(i), ' images'];
    %waitbar(i/H,hwait,str);
    disp(['Estimation of Image:', int2str(i),' complete!'])
  end
  %close(hwait);
  MfList = MfEstimate2(Data, WF, OTF, K, Phase);
  %Mf = MfEstimate(Data, OTF, K, Phase, ones(H,1));
  %MfList = zeros(H,2);
  %MfList(:,1) = 1;
  %MfList(:,2) = Mf;

  disp('Frequency:')
  disp(K);
  Kv = (K(:,1)+1i*K(:,2));

  disp('Km:');
  disp(abs(Kv));

  disp('Angle:');
  disp(mod(angle(Kv)/pi*180,180));
  p = Phase / pi * 180;
  p = mod(p,360);

  disp('Phase:');
  disp(p);

  disp('Modulation factor:');
  disp(MfList);
