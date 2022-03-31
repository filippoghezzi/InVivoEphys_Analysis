figure
hold on
plot(LFP(ops.L4best,ops.fs*60*12:ops.fs*60*17),'k')
plot(LFP(ops.L4best-6,ops.fs*60*12:ops.fs*60*17)+6000,'k')
plot(LFP(ops.L4best+8,ops.fs*60*12:ops.fs*60*17)-4000,'k')
plot(LFP(ops.L4best+16,ops.fs*60*12:ops.fs*60*17)-8000,'k')

U43=s.st(s.sclu==43);
U43=U43((U43>ops.fs*60*12) & (U43<ops.fs*60*17))-ops.fs*60*12;
plot([U43,U43],[1500,2000],'r-')

U164=s.st(s.sclu==164);
U164=U164((U164>ops.fs*60*12) & (U164<ops.fs*60*17))-ops.fs*60*12;
plot([U164,U164],[-2700,-2200],'r-')

U137=s.st(s.sclu==137);
U137=U137((U137>ops.fs*60*12) & (U137<ops.fs*60*17))-ops.fs*60*12;
plot([U137,U137],[-6300,-5800],'r-')

U162=s.st(s.sclu==163);
U162=U162((U162>ops.fs*60*12) & (U162<ops.fs*60*17))-ops.fs*60*12;
plot([U162,U162],[-10000,-9500],'r-')

plot([ops.fs*60*2,ops.fs*60*3],[-11,-11]*10^3,'k','LineWidth',1)% 1 min
plot([ops.fs*60*2,ops.fs*60*2],[-11,-10]*10^3,'k','LineWidth',1)% 1 mV

ax=gca;
ax.Box='off';
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 4\4.6\SB and units')

ax.XLim=[8.5,11.5]*10^5;
plot([0.9,0.95]*10^6,[-11,-11]*10^3,'k','LineWidth',1) % 1.5 s
print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 4\4.6\SB and units_Zoom')
