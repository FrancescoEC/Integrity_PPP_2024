
Crit=sqrt(satdataV1A2.Sat_Vel_X.^2+satdataV1A2.Sat_Vel_Y.^2+satdataV1A2.Sat_Vel_Z.^2);

Index=satdataV1A2.Sat_type=='GPS';

Crit_1=Crit<4000;
Pos=Crit_1 & Index;
SetGPS=unique(satdataV1A2.Sat_ID(Pos))

Index=satdataV1A2.Sat_type=='GALILEO';

Crit_2=Crit<4000;
Pos=Crit_2 & Index;
SetGAL=unique(satdataV1A2.Sat_ID(Pos))

Index=satdataV1A2.Sat_type=='GPS';


Crit_3=Crit>4000;
Pos=Crit_3 & Index;
SetGPS_LEO=unique(satdataV1A2.Sat_ID(Pos))

Index=satdataV1A2.Sat_type=='GALILEO';

Crit_4=Crit>4000;
Pos=Crit_4 & Index;
SetGAL_LEO=unique(satdataV1A2.Sat_ID(Pos));


