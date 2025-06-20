%----
R=1000;%ohm
C=100*10^(-9);%nf
V0=1;
%---------
fc=1/(R*C*2*pi);
%w=@(a)a;
w=2*pi*fc;
%----


fig=figure;

I0=V0./(R-(1i/(w*C)));
VR=I0*R;
VC=-(1i./(C.*w))*I0;
c=compass([V0,VR,VC]);
c1=c(1);
c1.LineWidth=2;

c2=c(2);
c2.LineWidth=2;
c2.Color='r';

c3=c(3);
c3.LineWidth=2;
c3.Color='g';

legend('V0','VR','VC');
set(gca,'FontSize',18)

%-----Slider---------
b = uicontrol('Parent',fig,'Style','slider','Position',[81,54,419,23],...
              'value',w, 'min',0.01*2*pi*fc, 'max',15*2*pi*fc);
b.Callback = @myFcn;

bgcolor = c.Color;
bl1 = uicontrol('Parent',fig,'Style','text','Position',[50,54,50,50],...
                'String','0.01fc','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',fig,'Style','text','Position',[500,54,50,50],...
                'String','15fc','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',fig,'Style','text','Position',[240,25,100,23],...
                'String','Frecuencia','BackgroundColor',bgcolor);


function myFcn(src,event)
w = src.Value;

R=1000;%ohm
C=100*10^(-9);%nf
V0=1;

I0=V0./(R-(1i/(w*C)));
VR=I0*R;
VC=-(1i./(C.*w))*I0;

c=compass(gca,[V0,VR,VC]);
c(1).LineWidth=2;
c(2).LineWidth=2;
c(2).Color='r';
c(3).LineWidth=2;
c(3).Color='g';
end