function [Pattern, K, Phase, MfList] = Estimate(Data, OTF, WF)
  Pattern = zeros(size(Data));
  M = size(Data,1);
  % [phase, k, MfList, I0List] = get_illumination_parament(Data, OTF);
  [Phase, K, MfList] = IlluParam(Data,OTF,WF);
  % Phase = Phase + pi;
  % Mf2 = MfEstimate(Data, OTF, K, Phase, ones(size(Phase)));
  for i = 1:size(Data,3)
    % Pattern(:,:,i) = get_pattern(k(i,:), phase(i), m, MfList(i), I0List(i));
    k = K(i,:);
    phase = Phase(i);
    Mf = MfList(2,i);
    % Pattern(:,:,i) = GenPattern(k, phase, M, Mf);
    % Pattern(:,:,i) = get_pattern(k,phase,M,1,1);
    Pattern(:,:,i) = get_pattern(k,phase,M,Mf,1);
  end
