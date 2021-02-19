function filteredData=BandPassButerworthFilter(data,sr)

Wn=[1.5/(sr/2),100/(sr/2)];
[b,a] = butter(1,Wn,'bandpass');
filteredData=filtfilt(b,a,data); 

end
