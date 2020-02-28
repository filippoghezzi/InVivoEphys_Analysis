function findEventTime(eventIdx)


LED=diff(LFP.eventArray(3,:))==1;

LED_ON=(LFP.eventArray(1,LED)==1)'; %1 is with laser, 0 is without

laser=diff(LFP.eventArray(1,:))==1;

laser_ON=(LFP.eventArray(3,laser)==1)';