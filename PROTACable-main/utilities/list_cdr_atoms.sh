#!/bin/sh

awk 'BEGIN{
       for(i=24;i<=34;i++)L[i]++; for(i=50;i<=56;i++)L[i]++; for(i=89;i<=97;i++)L[i]++;
       for(i=26;i<=32;i++)H[i]++; for(i=52;i<=56;i++)H[i]++; for(i=95;i<=102;i++)H[i]++;
     };
     /^@/{ina=ins=0};
     {if(ina) {id[++natm]=$1; res[natm]=$7; isH[natm]=(($6=="H")?1:0); }
      else if(ins){ split(substr($2,4),a,"[a-zA-Z]");
         if( ($6=="L" && (a[1] in L)) || ($6=="H" && (a[1] in H)) )hl[$1]++ }
     };
     /^@<TRIPOS>ATOM/{ina=1}; /^@<TRIPOS>SUBST/{ins=1};
     END{ for(i=1;i<=natm;i++)if(!isH[i] && (res[i] in hl))print id[i]}' $*
