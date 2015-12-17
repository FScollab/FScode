%NOTES.  Created by Jessica Nadalin & Mark Kramer, December 2015.

%  From [Traub J Neurophys 2003].
%  x. Use fast Na activation & inactivation.

%  From [Cunningham PNAS 2004 SI].
%  x.  No persistent Na
%  x.  No h-current
%  x.  Updated fast Na activation (done)        <--- REPLACED [12/17/15]
%  x.  Updated fast Na inactivation (done)      <--- REPLACED [12/17/15]
%  x.  Updated KDR (done)

% INPUTS:
%  T  = #steps in simulation
%  I0 = Initial current
%  gL = Leak conductance
%  gNaF = Fast sodium max conductance
%  gKDR = Fast potassium max conductance
%  gCaH = High-threshold calcium max conductance
%  gKM = M-current max conductance
%  gKv3 = Kv3.1 potassium max conductance
%  EK = Reversal potential for potassium currents
%  C = Capacitance
%  sigma = Noise level
%  ic = Initial conditions

% OUTPUTS:
%  V = Voltage
%  t = Simulation time
%  mNaF = Sodium channel activation
%  hNaF = Sodium channel inactivation
%  mKDR = Potassium delayed rectifier channel activation
%  mCaH = Calcium channel activation
%  mkV = Potassium Kv3 channel activation
%  mKM = Potassium muscarinic receptor-supressed channel activation
%  ic = Initial Conditions

function [V,t,mNaF,hNaF,mKDR,mCaH,mkV,mKM,ic] = traub_edit(T, I0, gL, gNaF, gKDR, gCaH, gKM, gKv3, EK, C, sigma, ic)

  dt = 0.005;               %Time step, milliseconds
 
  V = zeros(T,1);
  mNaF = zeros(T,1);
  hNaF = zeros(T,1);
  mKDR = zeros(T,1);
  mCaH = zeros(T,1);
  mKM = zeros(T,1);
  mkV = zeros(T,1);
  
  %Set the initial conditions (if passed as input).
  if isstruct(ic)
      V(1) = ic.V;
      mNaF(1) = ic.mNaF;
      hNaF(1) = ic.hNaF;
      mKDR(1) = ic.KDR;
      mCaH(1) = ic.mCaH;
      mKM(1)  = ic.mKM;
      mkV(1)   = ic.kV;
  end
      
  for i=1:T-1                                                   %For each time step,
      
      noise   = sqrt(dt).*randn().*sigma;                       %Update the noise,
      
      V(i+1) = V(i) + dt.*C*( -I0(i)  ...                       %Update the voltage,
          - gL*(70.0 + V(i)) + ...                              %... leak current,
          - gNaF.*mNaF(i).^3.0.*hNaF(i).*(-50.0 + V(i)) ...     %... NaF current,
          - (  gKDR.*mKDR(i).^4.0 ...                           %... KDR current,
             + gKM.*  mKM(i) ...                                %... M-current,
             + gKv3*   mkV(i)).*(-EK(i) + V(i)) ...              %... Kv3.1 current
          - gCaH.*mCaH(i).^2.0.*(-125 + V(i))) ...              %... high-threshold calcium current.
          + noise;
      
      mNaF(i+1) = mNaF(i) + dt*(alpha_mNaF(V(i))*(1-mNaF(i))  - beta_mNaF(V(i))*mNaF(i));       %Update fast sodium activation variable.
      hNaF(i+1) = hNaF(i) + dt*(alpha_hNaF(V(i))*(1-hNaF(i))  - beta_hNaF(V(i))*hNaF(i));       %Update fast sodium inactivation variable.
      mKDR(i+1) = mKDR(i) + dt*(alpha_mKDR(V(i))*(1-mKDR(i))  - beta_mKDR(V(i))*mKDR(i));       %Update fast potassium activation variable.
      mCaH(i+1) = mCaH(i) + dt*((1.6.*(1 - mCaH(i)))./(1 + exp(-0.072.*(-5 + V(i)))) - ...      %Update high threshold calcium activation variable.
                                (0.02.*mCaH(i).*(8.9 + V(i)))./(-1 + exp((8.9 + V(i))/5.)));
      mkV(i+1)   = mkV(i)   + dt*(alpha_KV(  V(i))*(1-  mkV(i))  - beta_KV(  V(i))*mkV(i));         %Update Kv3.1 activation variable.
      mKM(i+1)  = mKM(i)  + dt*(alpha_mKM( V(i))*(1- mKM(i))  - beta_mKM( V(i))*mKM(i));        %Update M-current activation variable.

  end
  
  %Set theinitial conditions to the last value of the output.
  ic = {};
  t = dt*(1:T);
  ic.V = V(end);
  ic.mNaF = mNaF(end);
  ic.hNaF = hNaF(end);
  ic.KDR = mKDR(end);
  ic.mCaH = mCaH(end);
  ic.mKM = mKM(end);
  ic.kV = mkV(end);

end
    
function res = alpha_mNaF(V)    % NaF activation [Traub J Neurophysiol 2003]. <--- REPLACED [12/17/15]
  minf = 1.0 ./ (1.0 + exp( (-V-34.5)/10 ));
  if V < -26.5
      taum = 0.0225 + 0.1425 * exp((V+26.5)/10);
  else
      taum = 0.0225 + 0.1425 * exp((-V-26.5)/10);
  end
  res = minf./taum;
end

function res = beta_mNaF(V)    % NaF activation  [Traub J Neurophysiol 2003]. <--- REPLACED [12/17/15]
  minf = 1.0 ./ (1.0 + exp( (-V-34.5)/10 ));
  if V < -26.5
      taum = 0.0225 + 0.1425 * exp((V+26.5)/10);
  else
      taum = 0.0225 + 0.1425 * exp((-V-26.5)/10);
  end
  res = (1.0 - minf)./taum;
end

function res = alpha_hNaF(V)    % NaF inactivation [Cunningham SI 2004].
  hinf = 1.0 ./ (1.0 + exp( (V+58.3)/6.7 ));
  tauh = 0.225 + 1.125 / (1.0+exp((V+37.0)/15.0));
  res = hinf./tauh;
end

function res = beta_hNaF(V)    % NaF inctivation  [Cunningham SI 2004].
  hinf = 1.0 ./ (1.0 + exp( (V+58.3)/6.7 ));
  tauh = 0.225 + 1.125 / (1.0+exp((V+37.0)/15.0));
  res = (1.0 - hinf)./tauh;
end

function res = alpha_mKDR(V)  % KDR activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-27.0)/11.5 ));
  if V < -10
      taum = 0.25 + 4.35 * exp((V+10)/10);
  else
      taum = 0.25 + 4.35 * exp((-V-10)/10);
  end
  res = minf./taum;
end

function res = beta_mKDR(V)  % KDR activation [Cunningham SI 2004].
  minf = 1.0 ./ (1.0 + exp( (-V-27.0)/11.5 ));
  if V < -10
      taum = 0.25 + 4.35 * exp((V+10)/10);
  else
      taum = 0.25 + 4.35 * exp((-V-10)/10);
  end
  res = (1.0 - minf)./taum;
end

function aK = alpha_KV(V)  % K-current forward rate function [Lien, 2003].
a = 0.0189324;
b = -4.18371;
c = 6.42606;
aK = a*(-((V)+b)) / (exp(-((V)+b)/c)-1);   
end

function bK = beta_KV(V)  % K-current backward rate function [Lien, 2003].
d = 0.015857;
e = 25.4834;
bK = d*exp(-(V)/e);
end

function alpha = alpha_mKM(V)  % 	  ;M-current forward rate function [Traub, 2003].
  alpha = 0.02./(1.0 + exp((-20 - V)/5.));
end

function beta = beta_mKM(V)  % 	  ;M-current backward rate function [Traub, 2003].
  beta = 0.01.*exp((-43 - V)/18.);
end