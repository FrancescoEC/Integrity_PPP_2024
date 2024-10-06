
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
SetGAL_LEO=unique(satdataV1A2.Sat_ID(Pos))


SetGPS=[4   7   9  16  19  20  27  30]
SetGAL=[1   2   3   8  11  12  13  14  18  19  20  22  25  26  29  31  32  33  34  35  36]


% SetGPS_LEO=[1   2   3   5   6   8  10  11  12  13  14  15  17  18  21  22  23  24  25  26  28  29  31  32 ]
% SetGAL_LEO=[4   5   6   7   9  10  15  16  17  21  23  24  27  28  30 ]
