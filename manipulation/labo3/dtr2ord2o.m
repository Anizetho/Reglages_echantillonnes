function [po,m,p,wo,wr,num,den] = dtr2ord2o(trx,x,d);
%   Function [po,m,p,wo,wr,num,den]=dtr2ord2(trx,x,d);
%
%  [num,den]=dtr2ord2(trx,x,d);
%
%   Pour un système du second ordre défini
%
%      trx : le temps de réponse à x % de la valeur asymptotique.
%      d   : le dépassement en %
%      x   : fourchette en % ou tr est défini.
%   Calcule
%
%      po  :  vecteur colonne des pôles
%      m   :  le coefficient d'amortissement réduit
%      wo  :  la pulsation naturelle non-amortie
%      wr  :  la pulsation de résonnance
%       p  :  le facteur de surtension
%

for k=1:length(trx),
   if ((d(k)<=0)|(d(k)>=100)),
      error('Le dépassement doit être supérieur à 0 % et inférieur à 100 %');
   end
   mk=log(d(k)/100);mk=sqrt(mk^2/(mk^2+pi^2));
   m=[];
   m=[m;mk];
   wok=-log(x(k)/100)/(mk*trx(k));
   wo=[];
   wo=[wo;wok];
   pk=1/(2*mk*sqrt(1-mk^2));
   p=[];
   p=[p;pk];
   wrk=wok*sqrt(1-mk^2);
   wr=[];
   wr=[wr;wrk];
   pok=[-wok*mk+j*wrk -wok*mk-j*wrk];
   po=[];
   po=[po;pok];
   denk=poly(pok');
   den=[];
   den=[den;denk];
   num=[];
   num=[num;denk(length(denk))];
end


if nargout==2,
   po=num;
   m=den;
elseif nargout==0,

   figh=findobj('name','Second Ordre fondamental');

   if isempty(figh)
     figh=figure;
     set(figh,'name','Second Ordre fondamental');
   else
     delete(figh);
     figh=figure;
     set(figh,'name','Second Ordre fondamental');
   end


   tmax=max(-ones(size(trx'))./real(po(:,1)))*5;
   tmin=min(-ones(size(trx'))./real(po(:,1)))/20;   
   t=[0:tmin:tmax]';
   wmax=max(wr);
   wmin=min(wr);
   w=logspace(floor(log10(wmin)),ceil(log10(wmax))+1,500);
   for k=1:length(trx),
     [am,ph]=bode(num(k,:),den(k,:),w);
     Am=[Am am];
     Ph=[Ph ph];
     y=[y step(num(k,:),den(k,:),t)];
   end
     
   ht=axes('position',[0.1 0.62 0.5 0.35],...
                    'box','on','FontSize',8,...
                    'TickLength',[0.0 0.0],'xscale','log','yscale','linear');
   line(w,20*log(Am));
   grid
   ylabel('Magnitude (dB)');

   ht=axes('position',[0.1 0.35 0.5 0.20],...
                    'box','on','FontSize',8,...
                    'TickLength',[0.0 0.0],'xscale','log','yscale','linear');
   line(w,Ph);
   grid
   xlabel('Pulsation (rad/s)');
   ylabel('Phase (deg)');

   ht=axes('position',[0.67 0.62 0.3 0.35],...
                    'box','on','FontSize',8,...
                    'TickLength',[0.0 0.0],'xscale','linear','yscale','linear');
   line(t,y);
   grid
   xlabel('time (s)');
   ylabel('Output');

   ht=axes('position',[0.67 0.35 0.3 0.2],...
                    'box','on','FontSize',8,...
                    'TickLength',[0.0 0.0],'xscale','linear','yscale','linear');
   line(real(po'),imag(po'),'linestyle','+');
   grid

   for k=1:length(trx),
        S = str2mat(S,[sprintf('%5.2e',trx(k)) ...
                       sprintf('        %5.3f',x(k)) ...
                       sprintf('       %5.2e',d(k)) ...
                       sprintf('    %4.2f',m(k)) ...
                       sprintf('     %5.2e',wo(k)) ...
                       sprintf('    %5.2e',real(po(k,1))) ...
                       sprintf('    %5.2e',imag(po(k,1))) ...
                       sprintf('    %5.2e',wr(k)) ...
                       sprintf('    %5.2f',p(k))]);
   end
   S = str2mat('  à x ( % )          ( x %)            en %        amort     naturelle           R(p)             +-I(p)                Wr          Réson',S);
   S = str2mat('Temps rep     fourchette     dépasse       coeff        Wo                          Pôles                    Puls Réson      Pic ',S);
   htext=axes('position',[0 0 1 0.25],'visible','off');
   St='Second Ordre Fondamental';
   clk=fix(clock);
   St=[ St ' Le ' date '  ' sprintf('%2.0i h %2.0i min',clk(1,4),clk(1,5)) ' arn@ecam.be'];
   h=text(0.1,1,St,'FontSize',10);
   [l,c]=size(S);
   ypos=(([-1:-1:-l]')*0.1)+0.9;
   xpos=zeros(size(ypos))+0.05;
   h=text(xpos,ypos,S,'FontSize',8);
 end