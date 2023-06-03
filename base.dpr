program U2343ptnov; {ðàñ÷åò îáîáùàåò ïðîãðàììó U232GI3 èñïîëüçóÿ íîâûé ìåòîä ñøèâêè ìèíîðíûõ èçîòîïîâ, íóìåðàöèÿ 1-232, 2-234, 3-235 è ò.ä., èñïîëüçóåòñÿ ìåòîä Õóêà-Äæèâñà äëÿ ìèíèìèçàöèè ðàçíîñòåé çàäàâàåìûõ è ðàññ÷èòûâàåìûõ çíà÷åíèé êîíöåíòðàöèé èçîòîïîâ 2-4 íà îòâàëå êàñêàäà.
Ìèíèìèçàöèÿ ñîäåðæèò ñèñòåìàòè÷åñêóþ îøèáêó}

{$APPTYPE CONSOLE}

uses
  SysUtils;

type mas2 = array [1..4,1..200] of extended;
       mas1 = array [1..200] of extended;
       mas =  array [1..5] of extended;
    var i,i1,iis,ip,j,j1,j2,k,p,m,jc,jj,sdv,pp,ppp,jn,jnn:integer;
        f,ff,f1,f2,f3,f4,f5,f11,f22,f33,f44,f55,ko1,ko2,ko3,ko4,hcm,hcp,snn,sn1,km,kp,ff3,ff4,
        aakv,bbkv,cckv,kkor,cp2,cp3,cp4,kpp1,koo1,koo2,koo3,koo4,kooo1,kooo2,kooo3,kooo4,TL,TLL,TL_1,TLL_1,
        kp11,kp1,kp2,kp3,kp4,km11,km1,km2,km3,km4,cmn3,cmv3,cmn2,cmv2,d,sg,sg1,cmn4,cmv4,smmm,s1,s2,s3,s4,
        T0,T,T1,gn,gv,gn1,gv1,teta,m1,m2,m3,m4,m5,eu,bet1,bet2,bet3,r1,r2,r3,sn,sig,smm,smm1,cz,TT,T00,T000:extended;
        ci,cn,cv,ggg,gg,g,gggp,ggp,gp,gggm,ggm,gm,tau1,tttet,ttet,tet,e,ee,r,rp,rm,n,n1,nnn,nn,zm,zp:mas1;
        xp,x1,x2,yp,y1,y2,ssx,sx,hhx,hx:mas1;
        ccc,cc,c,cccp,ccp,cp,cccm,ccm,cm,sm,ssm,sp,xi,xxi,taui1,ft:mas2;
        bl,bl1:boolean;
        a,l,fi,x1r,x2r,xpr,sxr,hxr,fa:mas;
  function potCm(c1,c2,c3,c4,c5:extended):extended;
     var f,f12,f13,f14,f15,f23,f24,f25,f34,f35,f45,
         b12,b13,b14,b15,b21,b23,b24,b25,b31,b32,b34,b35,
         b41,b42,b43,b45,b51,b52,b53,b54:extended;
     begin
        b15:=1; b21:=1; b13:=0; b14:=0; b12:=0; b23:=0; b24:=0; b25:=0;
        b34:=0; b43:=0; b35:=0; b53:=0; b45:=0; b54:=0;
        f:=(m2-m1)*(m2+m1-2*m3); b31:=1{-sqr(m3-m1)/f}; b32:=0{sqr(m3-m2)/f};
        f:=(m2-m1)*(m2+m1-2*m4); b41:=1{-sqr(m4-m1)/f}; b42:=0{sqr(m4-m2)/f};
        f:=(m2-m1)*(m2+m1-2*m5); b51:=1{-sqr(m5-m1)/f}; b52:=0{sqr(m5-m2)/f};
        f12:=1/sqr(m2-m1)*(b12*c1-b21*c2)*ln(c1/c2); f13:=1/sqr(m3-m1)*(b13*c1-b31*c3)*ln(c1/c3);
        f14:=1/sqr(m4-m1)*(b14*c1-b41*c4)*ln(c1/c4); f15:=1/sqr(m5-m1)*(b15*c1-b51*c5)*ln(c1/c5);
        f23:=1/sqr(m3-m2)*(b23*c2-b32*c3)*ln(c2/c3); f24:=1/sqr(m4-m2)*(b24*c2-b42*c4)*ln(c2/c4);
        f25:=1/sqr(m5-m2)*(b25*c2-b52*c5)*ln(c2/c5); f34:=1/sqr(m4-m3)*(b34*c3-b43*c4)*ln(c3/c4);
        f35:=1/sqr(m5-m3)*(b35*c3-b53*c5)*ln(c3/c5); f45:=1/sqr(m5-m4)*(b45*c4-b54*c5)*ln(c4/c5);
        potCM:=sqr(m5-m1)*(f12+f13+f14+f15+f23+f24+f25+f34+f35+f45)
     end;
 procedure soxr;
    begin
       for i:=1 to k do
          begin
             nn[i]:=n[i];
             gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
             cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
             ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
             ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
             cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
          end;
    end;
 procedure xxis;
    begin
       for i:=1 to 4 do
         for j:=1 to k do xxi[i,j]:=xi[i,j]
    end;
 procedure xis;
    begin
       for i:=1 to k do
         begin
            xi[3,i]:=exp((m5-m3)/(m5-m1)*ln(xi[1,i])); xi[2,i]:=exp((m5-m2)/(m5-m1)*ln(xi[1,i]));
            xi[4,i]:=exp((m5-m4)/(m5-m1)*ln(xi[1,i]))
         end;
    end;
 procedure rasn;
    var akv,bkv,ckv,kor,a0,a1,a2,a3:extended;
     begin
        f4:=0; a0:=1.0; a1:=1.0; a2:=1.0; a3:=0.2;
        for i:=1 to k do
           begin
              kor:=ln(xi[1,i]);
              n1[i]:=g[i]*exp(-1/a3*(a0+a1*tet[i]-a2*tet[i]*tet[i]-kor));
              f4:=f4+{(n1[i]-n[i])*(n1[i]-n[i])}g[i];
           end;
     end;
 procedure xar;
    begin
       for i:=1 to k do
        begin
           xi[1,i]:=exp(1+1*tet[i]-1*sqr(tet[i])-0.2*ln(g[i]/n[i]));
           xi[3,i]:=exp((m5-m3)/(m5-m1)*ln(xi[1,i])); xi[2,i]:=exp((m5-m2)/(m5-m1)*ln(xi[1,i]));
           xi[4,i]:=exp((m5-m4)/(m5-m1)*ln(xi[1,i]))
        end;
    end;
  procedure tran;
    begin
      for i:=1 to k do if i<=p then
                 tau1[i]:=T1 else tau1[i]:=T1-T0;
      for i:=1 to k do if i<=p then
                 taui1[3,i]:=T1*cm[3,1] else taui1[3,i]:=T1*cm[3,1]-T0*ko3;
      for i:=1 to k do if i<=p then
                 taui1[2,i]:=T1*cm[2,1] else taui1[2,i]:=T1*cm[2,1]-T0*ko2;
      for i:=1 to k do if i<=p then
                 taui1[1,i]:=T1*cm[1,1] else taui1[1,i]:=T1*cm[1,1]-T0*ko1;
      for i:=1 to k do if i<=p then
                 taui1[4,i]:=T1*cm[4,1] else taui1[4,i]:=T1*cm[4,1]-T0*ko4;
    end;

   procedure raspre;
    var i,j:integer;
    begin
           j2:=0; jn:=-999; jj:=0;
           {writeln(' raspre '); readln;}
           {writeln(' sm[2,1]= ',sm[2,1]:12:8,' sm[3,1]= ',sm[3,1]:12:8,' sm[4,1]= ',sm[4,1]:12:8);}
          { for i:=1 to k do if i<=p then
                 taui1[2,i]:=T1*sm[2,1] else taui1[2,i]:=T1*sm[2,1]-T0*ko2;
           for i:=1 to k do if i<=p then
                 taui1[3,i]:=T1*sm[3,1] else taui1[3,i]:=T1*sm[3,1]-T0*ko3;
           for i:=1 to k do if i<=p then
                 taui1[4,i]:=T1*sm[4,1] else taui1[4,i]:=T1*sm[4,1]-T0*ko4; }

           {cp2:=(T0*ko2-T1*sm[2,1])/T; cp3:=(T0*ko3-T1*sm[3,1])/T; cp4:=(T0*ko4-T1*sm[4,1])/T;
           f3:=1-(xi[1,k]-1)/xi[1,k]*kp1-(xi[2,k]-1)/xi[2,k]*cp2-(xi[3,k]-1)/xi[3,k]*cp3-(xi[4,k]-1)/xi[4,k]*cp4;
           sm[1,k]:=1/xi[1,k]*kp1/f3;}
            {writeln(' raspre 1'); readln;
            writeln(' sm[1,k]= ',sm[1,k]:12:8);}
           for i:=2 to k do
              begin
                 f3:=1+(xi[1,i-1]-1)*sm[1,i-1]+(xi[2,i-1]-1)*sm[2,i-1]+(xi[3,i-1]-1)*sm[3,i-1]+(xi[4,i-1]-1)*sm[4,i-1];
                 sp[1,i-1]:=xi[1,i-1]*sm[1,i-1]/f3; sp[2,i-1]:=xi[2,i-1]*sm[2,i-1]/f3;
                 sp[3,i-1]:=xi[3,i-1]*sm[3,i-1]/f3; sp[4,i-1]:=xi[4,i-1]*sm[4,i-1]/f3;

                 if sp[1,i-1]<kpp1 then
                          begin
{                              writeln(' i= ',i:4);
                             writeln(' sp1= ',sp[1,i-1]:12:9);  }
                             tau1[i]:=T1; taui1[1,i]:=T1*sm[1,1];
                             taui1[2,i]:=T1*sm[2,1]; taui1[3,i]:=T1*sm[3,1]; taui1[4,i]:=T1*sm[4,1];
                             if i>ppp then
                                begin
                                   tau1[i]:=tau1[i]-T000; taui1[1,i]:=taui1[1,i]-T000*kooo1;
                                   taui1[2,i]:=taui1[2,i]-T000*kooo2;  taui1[3,i]:=taui1[3,i]-T000*kooo3;
                                   taui1[4,i]:=taui1[4,i]-T000*kooo4;
                                end;
                              if i>pp then
                                 begin
                                   tau1[i]:=tau1[i]-T00; taui1[1,i]:=taui1[1,i]-T00*koo1;
                                   taui1[2,i]:=taui1[2,i]-T00*koo2; taui1[3,i]:=taui1[3,i]-T00*koo3;
                                   taui1[4,i]:=taui1[4,i]-T00*koo4;
                                 end;
                              if i>p then
                                 begin
                                    tau1[i]:=tau1[i]-T0; taui1[1,i]:=taui1[1,i]-T0*ko1;
                                    taui1[2,i]:=taui1[2,i]-T0*ko2; taui1[3,i]:=taui1[3,i]-T0*ko3;
                                    taui1[4,i]:=taui1[4,i]-T0*ko4;
                                 end;
                              cp2:=(T0*ko2+T00*koo2+T000*kooo2-T1*sm[2,1])/T;
                              cp3:=(T0*ko3+T00*koo3+T000*kooo3-T1*sm[3,1])/T;
                              cp4:=(T0*ko4+T00*koo4+T000*kooo4-T1*sm[4,1])/T;
                              f3:=1-(xi[1,k]-1)/xi[1,k]*kp1-(xi[2,k]-1)/xi[2,k]*cp2-(xi[3,k]-1)/xi[3,k]*cp3-(xi[4,k]-1)/xi[4,k]*cp4;
                              sm[1,k]:=1/xi[1,k]*kp1/f3;

                          end;
                       if sp[1,i-1]>=kpp1 then
                          begin
                             jj:=jj+1;
 {                            writeln(' i= ',i:4,' jj= ',jj:4);
                             writeln(' sp1= ',sp[1,i-1]:12:9);}
                             if jj=1 then
                                begin
                                   TL:=TT*(kpp1-sp[1,i-2])/(sp[1,i-1]-sp[1,i-2]);
                                   TL_1:=TT-TL; jn:=i-1;
                                   tau1[i]:=T1+TL_1; taui1[1,i]:=T1*sm[1,1]+TL_1*sp[1,i-2];
                                   taui1[2,i]:=T1*sm[2,1]+TL_1*sp[2,i-2]; taui1[3,i]:=T1*sm[3,1]+TL_1*sp[3,i-2];
                                   taui1[4,i]:=T1*sm[4,1]+TL_1*sp[4,i-2];

                                   f:=(T000*kooo2+T00*koo2+T0*ko2-T1*sm[2,1]-TL_1*sp[2,i-2]-TL*sp[2,i-1])/T;
                                   f3:=(T000*kooo3+T00*koo3+T0*ko3-T1*sm[3,1]-TL_1*sp[3,i-2]-TL*sp[3,i-1])/T;
                                   f4:=(T000*kooo4+T00*koo4+T0*ko4-T1*sm[4,1]-TL_1*sp[4,i-2]-TL*sp[4,i-1])/T;
                                   f5:=1-(xi[1,k]-1)/xi[1,k]*kp1-(xi[2,k]-1)/xi[2,k]*f-(xi[3,k]-1)/xi[3,k]*f3-(xi[4,k]-1)/xi[4,k]*f4;
                                   sm[1,k]:=1/xi[1,k]*kp1/f5; f:=1/xi[2,k]*f/f5; f3:=1/xi[3,k]*f3/f5; f4:=1/xi[4,k]*f4/f5;
                                   ff:=T1*sm[2,1]+TL_1*sp[2,i-2]+TL*sp[2,i-1];  ff3:=T1*sm[3,1]+TL_1*sp[3,i-2]+TL*sp[3,i-1];
                                   ff4:=T1*sm[4,1]+TL_1*sp[4,i-2]+TL*sp[4,i-1];
                                   if i>ppp then
                                       begin
                                          tau1[i]:=tau1[i]-T000; taui1[1,i]:=taui1[1,i]-T000*kooo1;
                                          taui1[2,i]:=taui1[2,i]-T000*kooo2; taui1[3,i]:=taui1[3,i]-T000*kooo3;
                                          taui1[4,i]:=taui1[4,i]-T000*kooo4;
                                       end;
                                   if i>pp then
                                       begin
                                          tau1[i]:=tau1[i]-T00; taui1[1,i]:=taui1[1,i]-T00*koo1;
                                          taui1[2,i]:=taui1[2,i]-T00*koo2; taui1[3,i]:=taui1[3,i]-T00*koo3;
                                          taui1[4,i]:=taui1[4,i]-T00*koo4;
                                       end;
                                    if i>p then
                                       begin
                                          tau1[i]:=tau1[i]-T0; taui1[1,i]:=taui1[1,i]-T0*ko1;
                                          taui1[2,i]:=taui1[2,i]-T0*ko2; taui1[3,i]:=taui1[3,i]-T0*ko3;
                                          taui1[4,i]:=taui1[4,i]-T0*ko4;
                                       end;
                                end;
                             if (jj>=2) and (i<=ppp) then
                                begin
                                   tau1[i]:=-T+T0+T00+T000; taui1[1,i]:=-T*kp1+T0*ko1+T00*koo1+T000*kooo1;
                                   taui1[2,i]:=ff;  taui1[3,i]:=ff3; taui1[4,i]:=ff4;
                                end;
                             if (jj>=2) and (i<=pp) and (i>ppp) then
                                begin
                                   tau1[i]:=-T+T0+T00; taui1[1,i]:=-T*kp1+T0*ko1+T00*koo1;
                                   taui1[2,i]:=ff-T000*kooo2;  taui1[3,i]:=ff3-T000*kooo3;
                                   taui1[4,i]:=ff4-T000*kooo4;
                                end;
                             if (jj>=2) and (i>pp) and (i<=p) then
                                begin
                                   tau1[i]:=-T+T0; taui1[1,i]:=-T*kp1+T0*ko1;
                                   taui1[2,i]:=ff-T000*kooo2-T00*koo2;  taui1[3,i]:=ff3-T000*kooo3-T00*koo3;
                                   taui1[4,i]:=ff4-T000*kooo4-T00*koo4;
                                end;
                             if (jj>=2) and (i>p) then
                                begin
                                   tau1[i]:=-T; taui1[1,i]:=-T*kp1;
                                   taui1[2,i]:=ff-T000*kooo2-T00*koo2-T0*ko2; taui1[3,i]:=ff3-T000*kooo3-T00*koo3-T0*ko3;
                                   taui1[4,i]:=ff4-T000*kooo4-T00*koo4-T0*ko4;
                                end;
                          end;

                 {writeln(' raspre 3',' i= ',i:3,' f3= ',f3:12:8); readln;
                 writeln(' sm[1,i-1]= ',sm[1,i-1]:14:10,' sm[2,i-1]= ',sm[2,i-1]:12:8,' sm[3,i-1]= ',sm[3,i-1]:12:8,' sm[4,i-1]= ',sm[4,i-1]:12:8);
                 writeln(' sp[1,i-1]= ',sp[1,i-1]:14:10,' sp[2,i-1]= ',sp[2,i-1]:12:8,' sp[3,i-1]= ',sp[3,i-1]:12:8,' sp[4,i-1]= ',sp[4,i-1]:12:8);
                 writeln(' sm[1,i]= ',sm[1,i]:14:10, 'sp[1,i-1]-sm[1,i]= ',(sp[1,i-1]-sm[1,i]):12:8); readln;}
                 gm[i]:=(tau1[i]*sp[1,i-1]-taui1[1,i])/(sp[1,i-1]-sm[1,i]);
                 {writeln(' gm[i]= ',gm[i]:12:8); readln;}
                 if i>1 then gp[i-1]:=gm[i]-tau1[i];
                 sm[2,i]:=(gp[i-1]*sp[2,i-1]+taui1[2,i])/gm[i];
                 sm[3,i]:=(gp[i-1]*sp[3,i-1]+taui1[3,i])/gm[i];
                 sm[4,i]:=(gp[i-1]*sp[4,i-1]+taui1[4,i])/gm[i];
                 {writeln(' i= ',i:3,' sm[2,i]= ',sm[2,i]:12:8,' sm[3,i]= ',sm[3,i]:12:8,' sm[4,i]= ',sm[4,i]:12:8);}
              end;
           {writeln(' raspre 2'); readln;}
           {cp2:=(T0*ko2-T1*sm[2,1])/T; cp3:=(T0*ko3-T1*sm[3,1])/T; cp4:=(T0*ko4-T1*sm[4,1])/T;}
           {f:=1-(xi[1,k]-1)/xi[1,k]*kp1-(xi[2,k]-1)/xi[2,k]*cp2-(xi[3,k]-1)/xi[3,k]*cp3-(xi[4,k]-1)/xi[4,k]*cp4;
           f2:=1/xi[2,k]*cp2/f; f3:=1/xi[3,k]*cp3/f; f4:=1/xi[2,k]*cp4/f;}

                       {begin
                          writeln(' sm[2,1]= ',sm[2,1]:12:10,' sm[2,k]à= ',sm[2,k]:12:10,' sm[2,k]¡= ',f1:12:10);
                          writeln('  sm[3,1]= ',sm[3,1]:12:10,' sm[3,k]à= ',sm[3,k]:12:10,' sm[3,k]¡= ',f2:12:10);
                          writeln('  gm[2]= ',gm[2]:15:7,' gm[k-1]= ',gm[k-1]:15:7,' gm[k]= ',gm[k]:15:7);
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;}



            f5:=1+(xi[1,k]-1)*sm[1,k]+(xi[2,k]-1)*sm[2,k]+(xi[3,k]-1)*sm[3,k]+(xi[4,k]-1)*sm[4,k];
            sp[1,k]:=xi[1,k]*sm[1,k]/f5; sp[2,k]:=xi[2,k]*sm[2,k]/f5;
            sp[3,k]:=xi[3,k]*sm[3,k]/f5; sp[4,k]:=xi[4,k]*sm[4,k]/f5;

            cp2:=(T000*kooo2+T00*koo2+T0*ko2-T1*sm[2,1]-TL_1*sp[2,jn-1]-TL*sp[2,jn])/T;
            cp3:=(T000*kooo3+T00*koo3+T0*ko3-T1*sm[3,1]-TL_1*sp[3,jn-1]-TL*sp[3,jn])/T;
            cp4:=(T000*kooo4+T00*koo4+T0*ko4-T1*sm[4,1]-TL_1*sp[4,jn-1]-TL*sp[4,jn])/T;

            f:=(sp[1,k]-kp1)*(sp[1,k]-kp1)+(sp[2,k]-cp2)*(sp[2,k]-cp2)+(sp[3,k]-cp3)*(sp[3,k]-cp3)+(sp[4,k]-cp4)*(sp[4,k]-cp4);

     end;
  procedure hukr;
     var bulr,bl1r:boolean;
         mr,j1r,mhr,iir,ir,i1r:integer; kpr:double;
     begin
        for ir:=1 to 4 do sm[ir,1]:=ssm[ir,1];
        x1r[2]:=ssm[2,1]; x1r[3]:=ssm[3,1]; x1r[4]:=ssm[4,1]; smm:=1e+18;
        bulr:=true; mhr:=50; kpr:=2; smm1:=1e+18;
        {writeln(' hukr ');}
        raspre;
        smm:=f;
        for mr:=1 to 20 do if bulr then
           begin
              for ir:=2 to 4 do
                 begin
                    sxr[ir]:=1;
                    hxr[ir]:=sm[ir,1]/1000;
                 end;
              for j1r:=2 to 4 do
                 begin
                    i1r:=0; iir:=0;
                    repeat
                       for ir:=2 to 4 do xpr[ir]:=ssm[ir,1];
                       {writeln(' j1r= ',j1r:3,' isx sm[j1r,1]= ',sm[j1r,1]:16:14);
                       writeln(' hxr[j1r]= ',hxr[j1r]:12:10);}
                       sm[j1r,1]:=sm[j1r,1]+sxr[j1r]*hxr[j1r];
                       {writeln(' j1r= ',j1r:3,' new sm[j1r,1]= ',sm[j1r,1]:16:14);}

                       bl1r:=true;
                       raspre; if iir>100 then bl1r:=false;
                       {for i:=1 to 4 do for j:=1 to k do
                          if (cm[i,j]<0) or (cp[i,j]<0) or (c[i,j]<0) then bl1r:=false;}
                       {writeln(' f= ',f:12:10,' smm= ',smm:12:10); readln;}
                       if bl1r then
                          begin

                             if f<smm then
                                begin
                                   iir:=iir+1;
                                   smm:=f;
                                   ssm[2,1]:=sm[2,1]; ssm[3,1]:=sm[3,1]; ssm[4,1]:=sm[4,1];
                                end
                             else
                                begin
                                  i1r:=i1r+1; sm[j1r,1]:=xpr[j1r];
                                  sxr[j1r]:=-sxr[j1r]; hxr[j1r]:=hxr[j1r]/2;
                                end
                          end
                       else
                          begin
                             i1r:=i1r+1; sm[j1r,1]:=xpr[j1r];
                             sxr[j1r]:=-sxr[j1r]; hxr[j1r]:=hxr[j1r]/2;
                          end;

                    until i1r > mhr;
                 end;
              for ir:=2 to 4 do x2r[ir]:=ssm[ir,1];
              for ir:=2 to 4 do hxr[ir]:=kpr*(x2r[ir]-x1r[ir]);
              i1r:=0; iir:=0;
              repeat
                 for ir:=2 to 4 do xpr[i]:=ssm[ir,1];
                 for ir:=2 to 4 do
                    begin
                       sm[ir,1]:=sm[ir,1]+hxr[ir];
                       {writeln(' ir= ',ir:3,' sm[ir]= ',sm[ir,1]:12:8);}
                    end;
                 bl1r:=true;
                 raspre;  if iir>100 then bl1r:=false;
                 {for i:=1 to 4 do for j:=1 to k do
                          if (cm[i,j]<0) or (cp[i,j]<0) or (c[i,j]<0) then bl1r:=false; }
                 if bl1r then
                    begin

                      {writeln('strelka f= ',f:12:8,' smm= ',smm:12:8);}
                       if f<smm then
                          begin
                             iir:=iir+1;
                             smm:=f;
                             ssm[2,1]:=sm[2,1]; ssm[3,1]:=sm[3,1]; ssm[4,1]:=sm[4,1];
                          end
                       else
                          begin
                             i1r:=i1r+1;
                             for ir:=2 to 4 do
                                begin
                                   sm[ir,1]:=xpr[ir]; hxr[ir]:=hxr[ir]/2;
                                end
                          end
                    end
                 else
                    begin
                       i1r:=i1r+1;
                       for ir:=2 to 4 do
                          begin
                             sm[ir,1]:=xpr[ir]; hxr[ir]:=hxr[ir]/2;
                          end
                    end;

              until i1r > mhr;
             { writeln(' smm = ',smm:15:10,' smm1 = ',smm1:15:10);}
              if abs(smm-smm1)< 1e-7 then
                begin
                  bulr:=false; bl:=true; iis:=1;
                  {writeln(' hukr iis= ',iis:5);}
                end;
              smm1:=smm;
              for ir:=2 to 4 do x1r[ir]:=x2r[ir];
              for ir:=1 to k do if (gm[ir]<0) then
                begin
                   bl:=false; iis:=0;
                end;
              {writeln(' mr = ',mr:3);}
           end;
            {writeln(' hukr iis= ',iis:5);}

     end;

  procedure ras;
     var f:extended;
     begin
        bl:=false;
        for i:=1 to 4 do for j:=1 to k do sm[i,j]:=cm[i,j];
        {writeln(' sm[1,1]= ',sm[1,1]:14:10,' sm[2,1]= ',sm[2,1]:12:8,' sm[3,1]= ',sm[3,1]:12:8,' sm[4,1]= ',sm[4,1]:12:8);
        writeln(' sm[1,2]= ',sm[1,2]:14:10,' sm[2,2]= ',sm[2,2]:12:8,' sm[3,2]= ',sm[3,2]:12:8,' sm[4,2]= ',sm[4,2]:12:8);
        writeln(' sm[1,3]= ',sm[1,3]:14:10,' sm[2,3]= ',sm[2,3]:12:8,' sm[3,3]= ',sm[3,3]:12:8,' sm[4,3]= ',sm[4,3]:12:8);}
        for i:=1 to 4 do for j:=1 to k do ssm[i,j]:=sm[i,j];
        hukr;
        cm[1,k]:=sm[1,k];
        for i:=1 to 4 do sm[i,1]:=ssm[i,1];
        raspre;
        for i:=1 to k do
           begin
              cp[1,i]:=sp[1,i]; cm[2,i]:=sm[2,i]; cp[2,i]:=sp[2,i];
              cm[3,i]:=sm[3,i]; cp[3,i]:=sp[3,i];
              cm[4,i]:=sm[4,i]; cp[4,i]:=sp[4,i];
           end;

        for i:=1 to k-1 do gp[i]:=gm[i+1]-tau1[i+1];
        for i:=1 to k do g[i]:=gp[i]+gm[i];

     for i:=1 to k do
         begin
            if (gm[i]<0) or (gp[i]<0) then iis:=0;
            if 1-cm[1,i]-cm[2,i]-cm[3,i]-cm[4,i]<0 then iis:=0;
         end;
     {writeln(' ras iis= ',iis:5);}
     if abs(cp[1,k]-kp1)>1e-3 then iis:=0;
     {writeln(' ras iis= ',iis:5);}
     end;
    procedure err;
     begin
        if (i<>1) and (i<>p) and (i<>k) then
           begin
              ee[i]:=ggp[i]*potCM(ccp[1,i],ccp[2,i],ccp[3,i],ccp[4,i],f1);
              ee[i]:=ee[i]+ggm[i]*potCM(ccm[1,i],ccm[2,i],ccm[3,i],ccm[4,i],f2);
              f11:=1-ccp[1,i-1]-ccp[2,i-1]-ccp[3,i-1]-ccp[4,i-1]; if f11<1e-1000 then f11:=1e-1000;
              f22:=1-ccm[1,i+1]-ccm[2,i+1]-ccm[3,i+1]-ccm[4,i+1]; if f22<1e-1000 then f22:=1e-1000;
              ee[i]:=ee[i]-ggp[i-1]*potCM(ccp[1,i-1],ccp[2,i-1],ccp[3,i-1],ccp[4,i-1],f11);
              ee[i]:=ee[i]-ggm[i+1]*potCM(ccm[1,i+1],ccm[2,i+1],ccm[3,i+1],ccm[4,i+1],f22);
           end;
           if i=1 then
              begin
                 ee[i]:=ggp[i]*potCM(ccp[1,i],ccp[2,i],ccp[3,i],ccp[4,i],f1);
                 ee[i]:=ee[i]+ggm[i]*potCM(ccm[1,i],ccm[2,i],ccm[3,i],ccm[4,i],f2);
                 f11:=1-ccm[1,i+1]-ccm[2,i+1]-ccm[3,i+1]-ccm[4,i+1]; if f11<1e-1000 then f11:=1e-1000;
                 ee[i]:=ee[i]-ggm[i+1]*potCM(ccm[1,i+1],ccm[2,i+1],ccm[3,i+1],ccm[4,i+1],f11);
                 if i=p then ee[i]:=ee[i]-T0*potCM(ko1,ko2,ko3,ko4,1-ko1-ko2-ko3-ko4);
               end;
           if i=k then
              begin
                 ee[i]:=ggp[i]*potCM(ccp[1,i],ccp[2,i],ccp[3,i],ccp[4,i],f1);
                 ee[i]:=ee[i]+ggm[i]*potCM(ccm[1,i],ccm[2,i],ccm[3,i],ccm[4,i],f2);
                 f11:=1-ccp[1,i-1]-ccp[2,i-1]-ccp[3,i-1]-ccp[4,i-1]; if f11<1e-1000 then f11:=1e-1000;
                 ee[i]:=ee[i]-ggp[i-1]*potCM(ccp[1,i-1],ccp[2,i-1],ccp[3,i-1],ccp[4,i-1],f11);
                 if i=p then ee[i]:=ee[i]-T0*potCM(ko1,ko2,ko3,ko4,1-ko1-ko2-ko3-ko4);
               end;
            if (i=p) and (i<>1) and (i<>k) then
               begin
                 ee[i]:=ggp[i]*potCM(ccp[1,i],ccp[2,i],ccp[3,i],ccp[4,i],f1);
                 ee[i]:=ee[i]+ggm[i]*potCM(ccm[1,i],ccm[2,i],ccm[3,i],ccm[4,i],f2);
                 f11:=1-ccp[1,i-1]-ccp[2,i-1]-ccp[3,i-1]-ccp[4,i-1]; if f11<1e-1000 then f11:=1e-1000;
                 f22:=1-ccm[1,i+1]-ccm[2,i+1]-ccm[3,i+1]-ccm[4,i+1]; if f22<1e-1000 then f22:=1e-1000;
                 ee[i]:=ee[i]-ggp[i-1]*potCM(ccp[1,i-1],ccp[2,i-1],ccp[3,i-1],ccp[4,i-1],f11);
                 ee[i]:=ee[i]-ggm[i+1]*potCM(ccm[1,i+1],ccm[2,i+1],ccm[3,i+1],ccm[4,i+1],f22);
                 ee[i]:=ee[i]-T0*potCM(ko1,ko2,ko3,ko4,1-ko1-ko2-ko3-ko4);
               end;
     end;
   procedure fun;
     begin
           bl1:=true;
           for i:=1 to k do
              begin
                 f:=1+(xi[1,i]-1)*cm[1,i]+(xi[2,i]-1)*ccm[2,i]+(xi[3,i]-1)*ccm[3,i]+(xi[4,i]-1)*ccm[4,i];
                 cp[1,i]:=xi[1,i]*cm[1,i]/f;
                 if i>1 then
                   gm[i]:=(tau1[i]*cp[1,i-1]-taui1[1,i])/(cp[1,i-1]-cm[1,i]);
                 if gm[i]<=0 then
                    begin
                       bl1:=false;
                     {  writeln(' fun - ¯®â®ª¨ ®âà¨æ â¥«ì­ë¥')}
                    end;
              end;
           if bl1 then
              begin
                 for i:=1 to k-1 do gp[i]:=gm[i+1]-tau1[i+1];
                 for i:=1 to k do
                    begin
                       g[i]:=gp[i]+gm[i]; tet[i]:=gp[i]/g[i];
                    end;
                 ip:=1; iis:=0;
                 ras;
                 if iis=1 then
                    begin
                       for i:=1 to k do
                          begin
                             c[1,i]:=(cm[1,i]*gm[i]+cp[1,i]*gp[i])/g[i];
                             c[2,i]:=(cm[2,i]*gm[i]+cp[2,i]*gp[i])/g[i];
                             c[3,i]:=(cm[3,i]*gm[i]+cp[3,i]*gp[i])/g[i];
                             c[4,i]:=(cm[4,i]*gm[i]+cp[4,i]*gp[i])/g[i];
                             tet[i]:=gp[i]/g[i];
                          end;
                    end
                   else bl1:=false;
              end
        end;
   procedure hukk;
     var bul:boolean;
         mh,ii:integer; kp:double;
     begin
        for i:=1 to 4 do for j:=1 to k do
         begin
            cm[i,j]:=ccm[i,j]; cp[i,j]:=ccp[i,j]; c[i,j]:=cc[i,j];
            ccm[i,j]:=ccm[i,j]; ccp[i,j]:=ccp[i,j]; cc[i,j]:=cc[i,j];
         end;
        for i:=1 to k do x1[i]:=ccm[1,i];
        bul:=true; mh:=20; kp:=2; sg1:=1e+18;

        fun; {sn:=1e+18;} cz:=cp[3,k];
        if bl1 then
           begin
              rasn; sn:=f4; smmm:=smm; cz:=cp[3,k];
              for i:=1 to k do
                 begin
                    nn[i]:=n1[i]; TLL:=TL; TLL_1:=TL_1;  jnn:=jn;
                    gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                    cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                    ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                    ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                    cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                 end;
           end;
        for m:=1 to 10 do if bul then
           begin
              for i:=2 to k-1 do
                 begin
                    sx[i]:=1;
                    hx[i]:=1e-10;
                 end;
              for j1:=2 to k-1 do
                 begin
                    i1:=0; ii:=0;
                    repeat
                       for i:=2 to k-1 do xp[i]:=ccm[1,i];
                       cm[1,j1]:=cm[1,j1]+sx[j1]*hx[j1];
                       if cm[1,j1] >= cm[1,j1+1] then cm[1,j1]:=cm[1,j1+1]-1e-12;
                       if cm[1,j1] <= cm[1,j1-1] then cm[1,j1]:=cm[1,j1-1]+1e-12;
                       fun; if ii>100 then bl1:=false;
                       {writeln('huk f4= ',f4:12:8,' sn= ',sn:12:8,' iis= ',iis:3);}
                       f1:=T0*ko1+T00*koo1+T000*kooo1-T*cp[1,k]-T1*cm[1,1]-cp[1,jnn-1]*TL_1-cp[1,jnn]*TL;
                       f2:=T0*ko2+T00*koo2+T000*kooo2-T*cp[2,k]-T1*cm[2,1]-cp[2,jnn-1]*TL_1-cp[2,jnn]*TL;
                       f3:=T0*ko3+T00*koo3+T000*kooo3-T*cp[3,k]-T1*cm[3,1]-cp[3,jnn-1]*TL_1-cp[3,jnn]*TL;
                       f4:=T0*ko4+T00*koo4+T000*kooo4-T*cp[4,k]-T1*cm[4,1]-cp[4,jnn-1]*TL_1-cp[4,jnn]*TL;
                       {writeln('huk f4= ',f4:12:8,' sn= ',sn:12:8,' iis= ',iis:3);}
                       bl1:=(f1<1e-12)and(f2<1e-12)and(f3<1e-13)and(f4<1e-12);
                       if bl1 then
                          begin
                             rasn;
                            { writeln(' fun - ãá¯¥å:áã¬¬  ¯®â®ª®¢ ',f4:15:10,' cm[2,1]= ',cm[2,1]:12:9);}
                             jj:=0;
                             if (f4<sn)and(cp[3,k]<cz)then
                                begin
                                   ii:=ii+1;
                                   sn:=f4; smmm:=smm; cz:=cp[3,k];
                                   for i:=1 to k do
                                      begin

                                         nn[i]:=n1[i]; TLL:=TL; TLL_1:=TL_1;  jnn:=jn;
                                         gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                         cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                         ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                         ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                         cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                      end;
                                end
                             else
                                begin
                                  i1:=i1+1; cm[1,j1]:=xp[j1];
                                  sx[j1]:=-sx[j1]; hx[j1]:=hx[j1]/2;
                                end
                          end
                       else
                          begin
                             i1:=i1+1; cm[1,j1]:=xp[j1];
                             sx[j1]:=-sx[j1]; hx[j1]:=hx[j1]/2;
                          end;
{                       writeln(' i1= ',i1:4,' cm[',j1:4,']*100=',cm[1,j1]*100:15:12,' sn= ',sn:12:7,' m= ',m:2);}
                    until i1 > mh;
                 end;
              for i:=1 to k do x2[i]:=ccm[1,i];
              for i:=1 to k do hx[i]:=kp*(x2[i]-x1[i]);
              i1:=0; ii:=0;
              repeat
                 for i:=1 to k do xp[i]:=ccm[1,i];
                 for i:=1 to k do
                    begin
                       cm[1,i]:=cm[1,i]+hx[i];
                       {if cm[1,i] >= cm[1,i+1] then cm[1,i]:=cm[1,i+1]-1e-12;
                       if cm[1,i] <= cm[1,i-1] then cm[1,i]:=cm[1,i-1]+1e-12;}
                    end;
                 fun;  if ii>100 then bl1:=false;
                 f1:=T0*ko1+T00*koo1+T000*kooo1-T*cp[1,k]-T1*cm[1,1]-cp[1,jnn-1]*TL_1-cp[1,jnn]*TL;
                 f2:=T0*ko2+T00*koo2+T000*kooo2-T*cp[2,k]-T1*cm[2,1]-cp[2,jnn-1]*TL_1-cp[2,jnn]*TL;
                 f3:=T0*ko3+T00*koo3+T000*kooo3-T*cp[3,k]-T1*cm[3,1]-cp[3,jnn-1]*TL_1-cp[3,jnn]*TL;
                 f4:=T0*ko4+T00*koo4+T000*kooo4-T*cp[4,k]-T1*cm[4,1]-cp[4,jnn-1]*TL_1-cp[4,jnn]*TL;
                 {writeln('huk f4= ',f4:12:8,' sn= ',sn:12:8,' iis= ',iis:3);}
                 bl1:=(f1<1e-12)and(f2<1e-12)and(f3<1e-13)and(f4<1e-12);
                 if bl1 then
                    begin
                       rasn;
                       jj:=0;
                       {writeln('huk strela ', f4:12:8,' sn= ',sn:12:8);}
                       if (f4<sn)and(cp[3,k]<cz)then
                          begin
                             ii:=ii+1;
                             sn:=f4; smmm:=smm; cz:=cp[3,k];
                             for i:=1 to k do
                                begin
                                   nn[i]:=n1[i]; TLL:=TL; TLL_1:=TL_1;  jnn:=jn;
                                   gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                   cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                   ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                   ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                   cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                end;
                          end
                       else
                          begin
                             i1:=i1+1;
                             for i:=1 to k do
                                begin
                                   cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                                end
                          end
                    end
                 else
                    begin
                       i1:=i1+1;
                       for i:=1 to k do
                          begin
                             cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                          end
                    end;
{                 writeln(' i1= ',i1:4,' áâà¥«ì¡  cm[2]*100=',cm[1,2]*100:15:12,' sn= ',sn:12:7);}
              until i1 > mh/2;
              if abs(sn-sg1) < 1e-10 then bul:=false;
              sg1:=sn;
              for i:=1 to k do x1[i]:=x2[i];
           end;
     end;
  procedure huk;
     var bul:boolean;
         mh,ii:integer; kp:double;
     begin
        for i:=1 to 4 do for j:=1 to k do
         begin
            cm[i,j]:=ccm[i,j]; cp[i,j]:=ccp[i,j]; c[i,j]:=cc[i,j];
            ccm[i,j]:=ccm[i,j]; ccp[i,j]:=ccp[i,j]; cc[i,j]:=cc[i,j];
         end;
        for i:=1 to k do x1[i]:=ccm[1,i];
        bul:=true; mh:=20; kp:=2; sg1:=1e+18;

        fun; {sn:=1e+18;}
        if bl1 then
           begin
              rasn; sn:=f4; smmm:=smm; cz:=cp[3,k];
              for i:=1 to k do
                 begin
                    nn[i]:=n1[i];
                    gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                    cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                    ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                    ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                    cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                 end;
           end;
        for m:=1 to 10 do if bul then
           begin
              for i:=2 to k-1 do
                 begin
                    sx[i]:=1;
                    hx[i]:=1e-10;
                 end;
              for j1:=2 to k-1 do
                 begin
                    i1:=0; ii:=0;
                    repeat
                       for i:=2 to k-1 do xp[i]:=ccm[1,i];
                       cm[1,j1]:=cm[1,j1]+sx[j1]*hx[j1];
                       if cm[1,j1] >= cm[1,j1+1] then cm[1,j1]:=cm[1,j1+1]-1e-12;
                       if cm[1,j1] <= cm[1,j1-1] then cm[1,j1]:=cm[1,j1-1]+1e-12;
                       fun; if ii>100 then bl1:=false;
                       {writeln('huk f4= ',f4:12:8,' sn= ',sn:12:8,' iis= ',iis:3);}
                       if bl1 then
                          begin
                             rasn;
                            { writeln(' fun - ãá¯¥å:áã¬¬  ¯®â®ª®¢ ',f4:15:10,' cm[2,1]= ',cm[2,1]:12:9);}
                             jj:=0;
                             if (f4<sn)then
                                begin
                                   ii:=ii+1;
                                   sn:=f4; smmm:=smm;
                                   for i:=1 to k do
                                      begin
                                         nn[i]:=n1[i];
                                         gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                         cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                         ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                         ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                         cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                      end;
                                end
                             else
                                begin
                                  i1:=i1+1; cm[1,j1]:=xp[j1];
                                  sx[j1]:=-sx[j1]; hx[j1]:=hx[j1]/2;
                                end
                          end
                       else
                          begin
                             i1:=i1+1; cm[1,j1]:=xp[j1];
                             sx[j1]:=-sx[j1]; hx[j1]:=hx[j1]/2;
                          end;
{                       writeln(' i1= ',i1:4,' cm[',j1:4,']*100=',cm[1,j1]*100:15:12,' sn= ',sn:12:7,' m= ',m:2);}
                    until i1 > mh;
                 end;
              for i:=1 to k do x2[i]:=ccm[1,i];
              for i:=1 to k do hx[i]:=kp*(x2[i]-x1[i]);
              i1:=0; ii:=0;
              repeat
                 for i:=1 to k do xp[i]:=ccm[1,i];
                 for i:=1 to k do
                    begin
                       cm[1,i]:=cm[1,i]+hx[i];
                       {if cm[1,i] >= cm[1,i+1] then cm[1,i]:=cm[1,i+1]-1e-12;
                       if cm[1,i] <= cm[1,i-1] then cm[1,i]:=cm[1,i-1]+1e-12;}
                    end;
                 fun;  if ii>100 then bl1:=false;
                 if bl1 then
                    begin
                       rasn;
                       jj:=0;
                       {writeln('huk strela ', f4:12:8,' sn= ',sn:12:8);}
                       if (f4<sn)then
                          begin
                             ii:=ii+1;
                             sn:=f4; smmm:=smm;
                             for i:=1 to k do
                                begin
                                   nn[i]:=n1[i];
                                   gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                   cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                   ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                   ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                   cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                end;
                          end
                       else
                          begin
                             i1:=i1+1;
                             for i:=1 to k do
                                begin
                                   cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                                end
                          end
                    end
                 else
                    begin
                       i1:=i1+1;
                       for i:=1 to k do
                          begin
                             cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                          end
                    end;
{                 writeln(' i1= ',i1:4,' áâà¥«ì¡  cm[2]*100=',cm[1,2]*100:15:12,' sn= ',sn:12:7);}
              until i1 > mh/2;
              if abs(sn-sg1) < 1e-10 then bul:=false;
              sg1:=sn;
              for i:=1 to k do x1[i]:=x2[i];
           end;
     end;
  procedure huk0;
     var bul:boolean;
         mh,k1,k2,ii:integer; kp,h:double;
     begin
        k1:=trunc((k-1)/2);
        for i:=1 to k do x1[i]:=ccm[1,i];
        bul:=true; mh:=20; kp:=2; sg1:=1e+18;
        for m:=1 to 50 do if bul then
           begin
              for i:=2 to k-1 do
                 begin
                    sx[i]:=1;
                    hx[i]:=1e-10;
                 end;
              for j1:=1 to k1 do
                 begin
                    i1:=0; ii:=0;
                    repeat
                       for i:=2 to k-1 do xp[i]:=ccm[1,i];
                       if j1>1 then
                          begin
                             cm[1,2*j1-1]:=cm[1,2*j1-1]+sx[2*j1-1]*hx[2*j1-1];
                             if cm[1,2*j1-1] >= cm[1,2*j1] then cm[1,2*j1-1]:=cm[1,2*j1]-1e-12;
                             if cm[1,2*j1-1] <= cm[1,2*j1-2] then cm[1,2*j1-1]:=cm[1,2*j1-2]+1e-12;
                          end;
                       cm[1,2*j1]:=cm[1,2*j1]+sx[2*j1]*hx[2*j1];
                       if cm[1,2*j1] >= cm[1,2*j1+1] then cm[1,2*j1]:=cm[1,2*j1+1]-1e-12;
                       if cm[1,2*j1] <= cm[1,2*j1-1] then cm[1,2*j1]:=cm[1,2*j1-1]+1e-12;
                       fun; if ii>100 then bl1:=false;
                       if bl1 then
                          begin
                             rasn;
                            { writeln(' fun - ãá¯¥å:áã¬¬  ¯®â®ª®¢ ',f4:15:10,' cm[2,1]= ',cm[2,1]:12:9);}
                               jj:=0;
                             if (f4<sn) then
                                begin
                                   ii:=ii+1;
                                   sn:=f4;
                                   for i:=1 to k do
                                      begin
                                         nn[i]:=n1[i];
                                         gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                         cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                         ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                         ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                         cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                      end;
                                end
                             else
                                begin
                                  i1:=i1+1;
                                  if j1>1 then
                                     begin
                                        cm[1,2*j1-1]:=xp[2*j1-1];
                                        sx[2*j1-1]:=-sx[2*j1-1]; hx[2*j1-1]:=hx[2*j1-1]/2;
                                     end;
                                  cm[1,2*j1]:=xp[2*j1];
                                  sx[2*j1]:=-sx[2*j1]; hx[2*j1]:=hx[2*j1]/2;
                                end
                          end
                       else
                          begin
                             i1:=i1+1;
                             if j1>1 then
                                begin
                                   cm[1,2*j1-1]:=xp[2*j1-1];
                                   sx[2*j1-1]:=-sx[2*j1-1]; hx[2*j1-1]:=hx[2*j1-1]/2;
                                end;
                              cm[1,2*j1]:=xp[2*j1];
                              sx[2*j1]:=-sx[2*j1]; hx[2*j1]:=hx[2*j1]/2;
                          end;
{                       if j1=1 then  writeln(' i1= ',i1:4,' cm[',2*j1:4,']*100=',cm[1,2*j1]*100:15:12,' sn= ',sn:12:7)
                       else writeln(' i1= ',i1:4,' cm[',2*j1-1:4,']*100=',cm[1,2*j1-1]*100:15:12,' sn= ',sn:12:7,' m=  ',m:3);}
                    until i1 > mh;
                 end;
              for i:=1 to k do x2[i]:=ccm[1,i];
              for i:=1 to k do hx[i]:=kp*(x2[i]-x1[i]);
              i1:=0; ii:=0;
              repeat
                 for i:=1 to k do xp[i]:=ccm[1,i];
                 for i:=1 to k do
                    begin
                       cm[1,i]:=cm[1,i]+hx[i];
                       if cm[1,i] >= cm[1,i+1] then cm[1,i]:=cm[1,i+1]-1e-12;
                       if cm[1,i] <= cm[1,i-1] then cm[1,i]:=cm[1,i-1]+1e-12;
                    end;
                 fun; if ii>100 then bl1:=false;
                 if bl1 then
                    begin
                       rasn;
                          jj:=0;
                          if (f4<sn) then
                          begin
                             ii:=ii+1;
                             sn:=f4;
                             for i:=1 to k do
                                begin
                                   nn[i]:=n1[i];
                                   gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                    cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                    ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                    ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                    cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                end;
                          end
                       else
                          begin
                             i1:=i1+1;
                             for i:=1 to k do
                                begin
                                   cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                                end
                          end
                    end
                 else
                    begin
                       i1:=i1+1;
                       for i:=1 to k do
                          begin
                             cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                          end
                    end;
{                 writeln(' i1= ',i1:4,' áâà¥«ì¡  cm[2]*100=',cm[1,2]*100:15:12,' sn= ',sn:12:7);}
              until i1 > mh/2;
              if abs(sn-sg1) < 1e-10 then bul:=false;
              sg1:=sn;
              for i:=1 to k do x1[i]:=x2[i];
           end;
     end;
  procedure huk1;
     var bul:boolean;
         mh,k1,k2,ii:integer; kp:double;
     begin
        k1:=trunc((k-1)/3);
        for i:=1 to k do x1[i]:=ccm[1,i];
        bul:=true; mh:=20; kp:=2; sg1:=1e+18;
        for m:=1 to 50 do if bul then
           begin
              for i:=2 to k-1 do
                 begin
                    sx[i]:=1;
                    hx[i]:=1e-10;
                 end;
              for j1:=1 to k1 do
                 begin
                    i1:=0; ii:=0;
                    repeat
                       for i:=2 to k-1 do xp[i]:=ccm[1,i];
                       if j1>1 then
                          begin
                             cm[1,3*j1-2]:=cm[1,3*j1-2]+sx[3*j1-2]*hx[3*j1-2];
                             if cm[1,3*j1-2] >= cm[1,3*j1-1] then cm[1,3*j1-2]:=cm[1,3*j1-1]-1e-12;
                             if cm[1,3*j1-2] <= cm[1,3*j1-3] then cm[1,3*j1-2]:=cm[1,3*j1-3]+1e-12;
                          end;
                       if j1<34 then cm[1,3*j1-1]:=cm[1,3*j1-1]+sx[3*j1-1]*hx[3*j1-1];
                       if cm[1,3*j1-1] >= cm[1,3*j1] then cm[1,3*j1-1]:=cm[1,3*j1]-1e-12;
                       if cm[1,3*j1-1] <= cm[1,3*j1-2] then cm[1,3*j1-1]:=cm[1,3*j1-2]+1e-12;
                       if j1<34 then cm[1,3*j1]:=cm[1,3*j1]+sx[3*j1]*hx[3*j1];
                       if cm[1,3*j1] >= cm[1,3*j1+1] then cm[1,3*j1]:=cm[1,3*j1+1]-1e-12;
                       if cm[1,3*j1] <= cm[1,3*j1-1] then cm[1,3*j1]:=cm[1,3*j1-1]+1e-12;
                       fun;  if ii>100 then bl1:=false;
                       if bl1 then
                          begin
                             rasn;
                            { writeln(' fun - ãá¯¥å:áã¬¬  ¯®â®ª®¢ ',f4:15:10,' cm[2,1]= ',cm[2,1]:12:9);}
                                jj:=0;
                             if (f4<sn) then
                                begin
                                   ii:=ii+1;
                                   sn:=f4;
                                   for i:=1 to k do
                                      begin
                                         nn[i]:=n1[i];
                                         gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                         cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                         ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                         ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                         cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                      end;
                                end
                             else
                                begin
                                  i1:=i1+1;
                                  if j1>1 then
                                     begin
                                        cm[1,3*j1-2]:=xp[3*j1-2];
                                        sx[3*j1-2]:=-sx[3*j1-2]; hx[3*j1-2]:=hx[3*j1-2]/2;
                                     end;
                                  cm[1,3*j1-1]:=xp[3*j1-1];
                                  sx[3*j1-1]:=-sx[3*j1-1]; hx[3*j1-1]:=hx[3*j1-1]/2;
                                  cm[1,3*j1]:=xp[3*j1];
                                  sx[3*j1]:=-sx[3*j1]; hx[3*j1]:=hx[3*j1]/2;
                                end
                          end
                       else
                          begin
                             i1:=i1+1;
                             if j1>1 then
                                begin
                                   cm[1,3*j1-2]:=xp[3*j1-2];
                                   sx[3*j1-2]:=-sx[3*j1-2]; hx[3*j1-2]:=hx[3*j1-2]/2;
                                end;
                             cm[1,3*j1-1]:=xp[3*j1-1];
                             sx[3*j1-1]:=-sx[3*j1-1]; hx[3*j1-1]:=hx[3*j1-1]/2;
                             cm[1,3*j1]:=xp[3*j1];
                             sx[3*j1]:=-sx[3*j1]; hx[3*j1]:=hx[3*j1]/2;
                          end;
{                       if j1=1 then  writeln(' i1= ',i1:4,' cm[',3*j1-1:4,']*100=',cm[1,3*j1-1]*100:15:12,' sn= ',sn:12:7)
                       else writeln(' i1= ',i1:4,' cm[',3*j1-2:4,']*100=',cm[1,3*j1-2]*100:15:12,' sn= ',sn:12:7,' m=  ',m:3);}
                    until i1 > mh;
                 end;
              for i:=1 to k do x2[i]:=ccm[1,i];
              for i:=1 to k do hx[i]:=kp*(x2[i]-x1[i]);
              i1:=0; ii:=0;
              repeat
                 for i:=1 to k do xp[i]:=ccm[1,i];
                 for i:=1 to k do
                    begin
                       cm[1,i]:=cm[1,i]+hx[i];
                       if cm[1,i] >= cm[1,i+1] then cm[1,i]:=cm[1,i+1]-1e-12;
                       if cm[1,i] <= cm[1,i-1] then cm[1,i]:=cm[1,i-1]+1e-12;
                    end;
                 fun;  if ii>100 then bl1:=false;


                 if bl1 then
                    begin
                       rasn;
                          jj:=0;
                          if (f4<sn) then
                          begin
                             ii:=ii+1;
                             sn:=f4;
                             for i:=1 to k do
                                begin
                                   nn[i]:=n1[i];
                                   gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                   cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                   ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                   ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                   cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                                end;
                          end
                       else
                          begin
                             i1:=i1+1;
                             for i:=1 to k do
                                begin
                                   cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                                end
                          end
                    end
                 else
                    begin
                       i1:=i1+1;
                       for i:=1 to k do
                          begin
                             cm[1,i]:=xp[i]; hx[i]:=hx[i]/2;
                          end
                    end;
{                 writeln(' i1= ',i1:4,' áâà¥«ì¡  cm[2]*100=',cm[1,2]*100:15:12,' sn= ',sn:12:7);}
              until i1 > mh/2;
              if abs(sn-sg1) < 1e-10 then bul:=false;
              sg1:=sn;
              for i:=1 to k do x1[i]:=x2[i];
           end;
     end;
  procedure  con;
     begin
    writeln(' ª®­æ¥­âà æ¨¨ 1-© ª®¬¯®­¥­âë:  U-235 ');
    for i:=1 to k do
           begin
              f:=ccm[1,i]; f1:=ccp[1,i]; f2:=cc[1,i];
              writeln(' cm[1,',i:2,']= ',f:13:11,' cp[1,',i:2,']= ',f1:13:11,' c[1,',i:2,']= ',f2:13:11);
              if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
           end;
    writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
    writeln(' ª®­æ¥­âà æ¨¨ 2-© ª®¬¯®­¥­âë:  U-234 ');
    for i:=1 to k do
       begin
          f:=ccm[2,i]; f1:=ccp[2,i]; f2:=cc[2,i];
          writeln(' cm[2,',i:2,']= ',f:12:10,' cp[2,',i:2,']= ',f1:12:10,' c[2,',i:2,']= ',f2:12:10);
          if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
       end;
    write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
    writeln(' ª®­æ¥­âà æ¨¨ 3-© ª®¬¯®­¥­âë:  U-232 %');
           for i:=1 to k do
              begin
                 f:=ccm[3,i]*100; f1:=ccp[3,i]*100; f2:=cc[3,i]*100;
                 writeln(' cm[3,',i:2,']= ',f:17:14,' cp[3,',i:2,']= ',f1:15:12,' c[3,',i:2,']= ',f2:14:11);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
    write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
   writeln(' ª®­æ¥­âà æ¨¨ 4-© ª®¬¯®­¥­âë:  U-236 ');
           for i:=1 to k do
              begin
                 f:=ccm[4,i]; f1:=ccp[4,i]; f2:=cc[4,i];
                 writeln(' cm[4,',i:2,']= ',f:12:10,' cp[4,',i:2,']= ',f1:12:10,' c[4,',i:2,']= ',f2:12:10);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
    write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
    writeln(' ª®­æ¥­âà æ¨¨ 5-© ª®¬¯®­¥­âë U-238');
           for i:=1 to k do
              begin
                 f1:=1-ccm[1,i]-ccm[2,i]-ccm[3,i]-ccm[4,i]; f2:=1-ccp[1,i]-ccp[2,i]-ccp[3,i]-ccp[4,i];
                 f:=1-cc[1,i]-cc[2,i]-cc[3,i]-cc[4,i];
                 writeln(' cm[5,',i:2,']= ',f1:12:10,' cp[5,',i:2,']= ',f2:12:10,' c[5,',i:2,']= ',f:12:10);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
    write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
     writeln(' sigma ');
    for i:=1 to k do
       begin
          f:=ggp[i]*ccp[1,i]/gg[i]/cc[1,i];
          f1:=f/xi[1,i]/(1-f);
          f:=ggp[i]*ccp[2,i]/gg[i]/cc[2,i];
          f2:=f/xi[2,i]/(1-f);

          writeln(' sig[1,',i:2,']= ',f1:12:9,' sig[2,',i:2,']= ',f2:12:9);
          if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
       end;
    write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
     end;
  procedure potok;
     begin
       sg:=0;
    writeln(' ¯®â®ª¨ áâã¯¥­¥© ');
    for i:=1 to k do
       begin
          f:=ggm[i]; f1:=ggp[i]; f2:=gg[i]; sg:=sg+gg[i];
          writeln(' gm[',i:2,']= ',f:12:10,' gp[',i:2,']= ',f1:12:10,' g[',i:2,']= ',f2:12:10);
          if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
       end;
     writeln('  «ãçè ï áã¬¬  ¯®â®ª®¢ ¯¨â ­¨ï - ', sg:15:10);
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
    {writeln(' ¤ ¢«¥­¨ï ¢ áâã¯¥­ïå ');
    for i:=1 to k do
       begin
          f:=(0.84*gg[i]-ggm[i])/nn[i]/0.00299; if f>0 then f:=sqrt(f) else f:=-999;
          f1:=gg[i]/nn[i]/3.5;
          writeln(' pt[',i:2,']= ',f:12:7,' po[',i:2,']= ',f1:12:7);
          if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
       end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
      writeln(' ª®«¨ç¥áâ¢  æ¥­âà¨äã£ ¢ áâã¯¥­ïå ');
    for i:=1 to k do
       begin
          n1[i]:=nn[i];
          f:=n1[i]; f1:=gg[i]/n1[i]; f2:=ggp[i]/gg[i];
          writeln(' n[',i:2,']= ',f:12:7,' gæ[',i:2,']= ',f1:12:7,' tet[',i:2,']= ',f2:6:4,' xi1[',i:2,']= ',xi[1,i]:6:4);
          if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
       end;
        writeln('  «ãçè ï áã¬¬  ¯®â®ª®¢ ¯¨â ­¨ï - ', sg:15:10);
     if TT<>0 then
        begin
           s1:=(ccp[1,jnn-1]*TLL_1+ccp[1,jnn]*TLL)/TT;  s2:=(ccp[2,jnn-1]*TLL_1+ccp[2,jnn]*TLL)/TT;
           s3:=(ccp[3,jnn-1]*TLL_1+ccp[3,jnn]*TLL)/TT;  s4:=(ccp[4,jnn-1]*TLL_1+ccp[4,jnn]*TLL)/TT;
           writeln(' jn= ',jnn:3,' TL1= ',TLL_1:10:7,' TL= ',TLL:10:7,' s1= ',s1:12:10,' s2 (%)= ',s2*100:14:12);
           writeln(' s3 (%)= ',s3*100:17:14,' s4 (%)= ',s4*100:14:12);
        end;

     {writeln('  «ãçè ï áã¬¬  ª¢ ¤à â®¢ à §­®áâ¥© ª®«¨ç¥áâ¢ æ¥­âà¨äã£ - ', sn:16);}
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
     end;
  procedure hukx;
     var bula,bb:boolean;
         mhh,i3,i2,j2,m1,iii:integer; kpp:double;
     begin
        for i3:=1 to k do y1[i3]:=xxi[1,i3];
        bula:=true; mhh:=6; kpp:=2; sn1:=1e+18;
        for m1:=1 to 8 do if bula then
           begin
              for i3:=1 to k do
                 begin
                    ssx[i3]:=1;
                    hhx[i3]:=y1[i3]/10;
                 end;
              for j2:=1 to k do
                 begin
                    i2:=0; iii:=0; bb:=false;
                    repeat
                       for i3:=1 to k do yp[i3]:=xxi[1,i3];
                       xi[1,j2]:=xi[1,j2]+ssx[j2]*hhx[j2];
                       if xi[1,j2]<1 then xi[1,j2]:=1; if xi[1,j2]>7 then xi[1,j2]:=7;
               {        writeln(' sx[',j2:3,']= ',ssx[j2]:12:7,' hx[',j2:3,']= ',hhx[j2]:12:7);}
                       xis; huk;
                          writeln(' j2=',j2:3,' xi[1,',j2:3,']= ',xi[1,j2]:7:5,' sn= ',sn:12:8,' snn= ',snn:12:8,' smmm=',smmm:12:10);
                          {writeln(' cm[5,k-1]=',1-cm[1,k-1]-cm[2,k-1]-cm[3,k-1]-cm[4,k-1]:15:13,' cm[5,k]=',1-cm[1,k]-cm[2,k]-cm[3,k]-cm[4,k]:15:13);}
                          begin
                             if sn<snn then
                                begin
                                   iii:=iii+1;  if iii>10 then bb:=true;
                                   snn:=sn;
                                   soxr; xxis;
                                   {f:=0; for i:=1 to k do f:=f+(nnn[i]-n[i])*(nnn[i]-n[i]);
                                   writeln('                  snn=',snn:12:8,' summa kvadratov N ',f:12:8);}
                                end
                             else
                                begin
                                  i2:=i2+1; xi[1,j2]:=yp[j2];
                                  ssx[j2]:=-ssx[j2]; hhx[j2]:=hhx[j2]/2;
                                end
                          end;
{                    writeln(' i2= ',i2:4,' xi1[',j2:2,']= ',xi[1,j2]:12:7,' sn= ',snn:16:12,' m1= ',m1:3);}
                    until (i2 > mhh) or bb;
                 end;

              for i3:=1 to k do y2[i3]:=xxi[1,i3];
              for i3:=1 to k do hhx[i3]:=kpp*(y2[i3]-y1[i3]);
              i2:=0; iii:=0; bb:=false;
              repeat
                 for i3:=1 to k do yp[i3]:=xxi[1,i3];
                 for i3:=1 to k do
                   begin
                      xi[1,i3]:=xi[1,i3]+hhx[i3];
                      if xi[1,i3]<1 then xi[1,i3]:=1; if xi[1,i3]>7 then xi[1,j2]:=7;
                   end;
                 xis; huk;
                    begin
                        if sn<snn then
                           begin
                              iii:=iii+1; if iii>10 then bb:=true;
                              snn:=sn;
                              soxr; xxis;
                          {f:=0; for i:=1 to k do f:=f+(nnn[i]-n[i])*(nnn[i]-n[i]);
                                   writeln(' strelka   snn=',snn:12:8,' summa kvadratov N ',f:12:8);}
                           end
                        else
                           begin
                              i2:=i2+1;
                              for i3:=1 to k do
                                 begin
                                    xi[1,i3]:=yp[i3];
                                    hhx[i3]:=hhx[i3]/2;
                                 end;
                           end
                    end;
              {writeln(' áâà¥«ª   i2= ',i2:4,' sn= ',snn:16:12);}
              until (i2 > mhh/2) or bb;
              if abs(snn-sn1) < 1e-10 then bula:=false;
              sn1:=snn;
              writeln('  sn= ',snn:16:12,' m1= ',m1:3);
              for i3:=1 to k do y1[i3]:=y2[i3];
           end;
     end;

begin
  { TODO -oUser -cConsole Main : Insert code here }
   k:=28; p:=14;

      ko1:=0.0085; ko2:=0.00016; ko3:=1.5e-9;  ko4:=0.0035;

     ppp:=1; T000:=0/8.76/3.6; kooo1:=0.003; kooo2:=0.000015; kooo3:=1e-15; kooo4:=1e-15;
     pp:=13; T00:=0/8.76/3.6; koo1:=0.003; koo2:=0.000015;  koo3:=1e-15; koo4:=1e-15;

     m1:=235; m2:=234; m3:=232; m4:=236; m5:=238;
     {m1:=70; m2:=72; m3:=73; m4:=74; m5:=76; iz:=3;} eu:=0.0141;

     {km3:=6e-11; km2:=3e-5; km1:=0.003; km4:=0.002072;
     kp3:=1.08e-8; kp2:=0.001; kp1:=0.0495; kp4:=0.01272;}
     {T0:=10; T:=T0*(ko1-km1)/(kp1-km1); T1:=T0-T;
     writeln(' T1= ',T1:12:8,' T0= ',T0:12:8,' T= ',T:12:8);
     km:=km1; kp:=kp1;}
     Randomize;

     f:=1.3;
     for i:=1 to k do
        begin
           xi[1,i]:=f; xi[2,i]:=exp((m5-m2)/(m5-m1)*ln(xi[1,i]));
           xi[3,i]:=exp((m5-m3)/(m5-m1)*ln(xi[1,i])); xi[4,i]:=exp((m5-m4)/(m5-m1)*ln(xi[1,i]))
        end;
     writeln(' xi1= ',xi[1,1]:12:8,' eps15= ',ln(xi[1,1]):12:8);
     writeln(' xi2= ',xi[2,1]:12:8,' eps25= ',ln(xi[2,1]):12:8,' xi3= ',xi[3,1]:12:8,' eps35= ',ln(xi[3,1]):12:8);
     writeln(' xi4= ',xi[4,1]:12:8,' eps35= ',ln(xi[4,1]):12:8);
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;

     f:=xi[1,1];{á¨¬¬¥âà¨ç­ë© ¯® 1,5} sig:=1/sqrt(f);
     writeln(' sig= ',sig:15:8);

     fi[1]:=sig*xi[1,1]/(1+sig*xi[1,1]);
     fi[5]:=sig/(1+sig);
     fi[3]:=sig*xi[3,1]/(1+sig*xi[3,1]); fi[4]:=sig*xi[4,1]/(1+sig*xi[4,1]);
     fi[2]:=sig*xi[2,1]/(1+sig*xi[2,1]);
     writeln('fi[1]= ',fi[1]:12:8,' fi[2]= ',fi[2]:12:8,' fi[3]= ',fi[3]:12:8,' fi[4]= ',fi[4]:12:8,' fi[5]= ',fi[5]:12:8);
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;


     for i:=1 to 5 do if abs(fi[i]-0.5)>1e-15 then a[i]:=1/(2*fi[i]-1) else a[i]:=1e+15;
     writeln(' a[1]= ',a[1]:15:12,' a[2]= ',a[2]:15:12,' a[3]= ',a[3]:15:12,' a[4]= ',a[4]:15:12,' a[5]= ',a[5]:15:12);
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;

     for i:=1 to 5 do l[i]:=fi[i]/(1-fi[i]);
     writeln(' l[1]= ',l[1]:15:12,' l[2]= ',l[2]:15:12,' l[3]= ',l[3]:15:12,' l[4]= ',l[4]:15:12,' l[5]= ',l[5]:15:12);
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;

     if abs(l[1]-1)>1e-10 then f1:=ko1/(exp(p*ln(l[1]))-exp((p-k-1)*ln(l[1])))*(exp(p*ln(l[1]))-1)
     else f1:=ko1*p/(k+1);
     if abs(l[2]-1)>1e-10 then f2:=ko2/(exp(p*ln(l[2]))-exp((p-k-1)*ln(l[2])))*(exp(p*ln(l[2]))-1)
     else f2:=ko2*p/(k+1);
     if abs(l[3]-1)>1e-10 then f3:=ko3/(exp(p*ln(l[3]))-exp((p-k-1)*ln(l[3])))*(exp(p*ln(l[3]))-1)
     else f3:=ko3*p/(k+1);
     if abs(l[4]-1)>1e-10 then f4:=ko4/(exp(p*ln(l[4]))-exp((p-k-1)*ln(l[4])))*(exp(p*ln(l[4]))-1)
      else f4:=ko4*p/(k+1);
     f5:=(1-ko1-ko2-ko3-ko4)/(exp(p*ln(l[5]))-exp((p-k-1)*ln(l[5])))*(exp(p*ln(l[5]))-1);
     f:=f1+f2+f3+f4+f5; cp[1,k]:={kp1}f1/f; cp[2,k]:=f2/f; cp[3,k]:=f3/f; cp[4,k]:=f4/f;






     if abs(l[1]-1)>1e-10 then f1:=ko1/(exp(p*ln(l[1]))-exp((p-k-1)*ln(l[1])))*(1-exp((p-k-1)*ln(l[1])))
     else f1:=ko1*(k-p+1)/(k+1);
     if abs(l[2]-1)>1e-10 then f2:=ko2/(exp(p*ln(l[2]))-exp((p-k-1)*ln(l[2])))*(1-exp((p-k-1)*ln(l[2])))
     else f2:=ko2*(k-p+1)/(k+1);
     if abs(l[3]-1)>1e-10 then f3:=ko3/(exp(p*ln(l[3]))-exp((p-k-1)*ln(l[3])))*(1-exp((p-k-1)*ln(l[3])))
     else f3:=ko3*(k-p+1)/(k+1);
     if abs(l[4]-1)>1e-10 then f4:=ko4/(exp(p*ln(l[4]))-exp((p-k-1)*ln(l[4])))*(1-exp((p-k-1)*ln(l[4])))
      else f4:=ko4*(k-p+1)/(k+1);
     f5:=(1-ko1-ko2-ko3-ko4)/(exp(p*ln(l[5]))-exp((p-k-1)*ln(l[5])))*(1-exp((p-k-1)*ln(l[5])));
     f:=f1+f2+f3+f4+f5; cm[1,1]:={km1}f1/f; cm[2,1]:=f2/f; cm[3,1]:=f3/f; cm[4,1]:=f4/f;
     kp1:=cp[1,k]; km1:=cm[1,1];
     teta:=(ko1-km1)/(kp1-km1); T0:=1/teta; T:=T0*teta; T1:=T0-T;
     T0:=100/8.76/3.6; T:=T0*(ko1-km1)/(kp1-km1); T1:=T0-T;

     writeln('  teta= ',teta:8:6,'  T= ',T:15:12,'  T1= ',T1:15:12,' T0= ',T0:15:12);
     writeln(' cm[1,1]= ',cm[1,1]:22:19,' cm[2,1]= ',cm[2,1]:22:19,' cm[3,1]= ',cm[3,1]:22:19);
     writeln(' cp[1,k]= ',cp[1,k]:22:19,' cp[2,k]= ',cp[2,k]:22:19,' cp[3,k]= ',cp[3,k]:22:19);
     writeln(' cm[4,1]= ',cm[4,1]:22:19,' cp[4,k]= ',cp[4,k]:22:19);
     f1:=T0*ko1-T*cp[1,k]-T1*cm[1,1]; f2:=T0*ko2-T*cp[2,k]-T1*cm[2,1]; f3:=T0*ko3-T*cp[3,k]-T1*cm[3,1];
     f4:=T0*ko4-T*cp[4,k]-T1*cm[4,1];

     f5:=T0*(1-ko1-ko2-ko3-ko4)-T*(1-cp[1,k]-cp[2,k]-cp[3,k]-cp[4,k])-T1*(1-cm[1,1]-cm[2,1]-cm[3,1]-cm[4,1]);

     writeln('¡ « ­áë: 1 ',f1:15:12,' 2 ',f2:15:12,' 3 ',f3:15:12,' 4 ',f4:15:12,' 5 ',f5:15:12);
     f1:=1-ko1-ko2-ko3-ko4; f2:=1-cp[1,k]-cp[2,k]-cp[3,k]-cp[4,k]; f3:=1-cm[1,1]-cm[2,1]-cm[3,1]-cm[4,1];
     writeln(' ª®­æ¥­âà æ¨ï 5: ¯¨â ­. ',f1:15:12,' ®â¡. ',f2:15:12,' ®â¢. ',f3:15:12);
     {f1:=1-cp[1,k]-cp[2,k]-cp[3,k]-cp[4,k]; if f1<=1e-18 then f1:=1e-12;
     f2:=1-cm[1,1]-cm[2,1]-cm[3,1]-cm[4,1]; if f2<= 1e-18 then f2:= 1e-12;
     f3:=1-ko1-ko2-ko3-ko4;
     f1:=potCM(cp[1,k],cp[2,k],cp[3,k],cp[4,k],f1); f2:=potCM(cm[1,1],cm[2,1],cm[3,1],cm[4,1],f2);
     f3:=potCM(ko1,ko2,ko3,ko4,f3);
     f4:=T*f1+T1*f2-T0*f3;
     writeln('  ¨¤¥ «ì­ ï áã¬¬  ¯®â®ª®¢ ',f4/eu:15:10); }
     {for i:=1 to k do
        if i<=p then
           begin
              tau1[i]:=T1; taui1[1,i]:=T1*km1;
           end
         else
           begin
              tau1[i]:=-T; taui1[1,i]:=-T*kp1;
           end; }
     for i:=1 to k do
        begin
         {  if i>=p then f1:=T*cp[1,k]*a[1]*(1-exp((i-k-1)*ln(l[1])));
           if i<p then f1:=T1*cm[1,1]*a[1]*(exp(i*ln(l[1]))-1);
           if i>=p then f2:=T*cp[2,k]*a[2]*(1-exp((i-k-1)*ln(l[2])));
           if i<p then f2:=T1*cm[2,1]*a[2]*(exp(i*ln(l[2]))-1);
           if i>=p then f3:=T*cp[3,k]*a[3]*(1-exp((i-k-1)*ln(l[3])));
           if i<p then f3:=T1*cm[3,1]*a[3]*(exp(i*ln(l[3]))-1);}
           for j:=1 to 4 do
              begin
                 if abs(l[j]-1)>1e-10 then
                   begin
                      if i>=p then fa[j]:=T*cp[j,k]*a[j]*(1-exp((i-k-1)*ln(l[j])));
                      if i<p then fa[j]:=T1*cm[j,1]*a[j]*(exp(i*ln(l[j]))-1)
                   end
                  else
                   begin
                      if i>=p then fa[j]:=2*T*cp[j,k]*(k-i+1);
                      if i<p then fa[j]:=2*T1*cm[j,1]*i
                   end;
              end;
           f1:=fa[1]; f2:=fa[2]; f3:=fa[3]; f4:=fa[4];
           if i>=p then
                 f5:=T*(1-cp[1,k]-cp[2,k]-cp[3,k]-cp[4,k])*a[5]*(1-exp((i-k-1)*ln(l[5])));
           if i<p then f5:=T1*(1-cm[1,1]-cm[2,1]-cm[3,1]-cm[4,1])*a[5]*(exp(i*ln(l[5]))-1);
           f:=f1+f2+f3+f4+f5;
           c[1,i]:=f1/f; c[2,i]:=f2/f; c[3,i]:=f3/f; c[4,i]:=f4/f;
           f11:=fi[1]*f1; f22:=fi[2]*f2; f33:=fi[3]*f3;
           f44:=fi[4]*f4; f55:=fi[5]*f5;
           f:=f11+f22+f33+f44+f55;
           gp[i]:=f;
           cp[1,i]:=f11/f; cp[2,i]:=f22/f; cp[3,i]:=f33/f; cp[4,i]:=f44/f;

           f11:=f1-f11; f22:=f2-f22; f33:=f3-f33; f44:=f4-f44; f55:=f5-f55;
           f:=f11+f22+f33+f44+f55;
           gm[i]:=f;
           cm[1,i]:=f11/f; cm[2,i]:=f22/f; cm[3,i]:=f33/f; cm[4,i]:=f44/f;
           g[i]:=gp[i]+gm[i]; tet[i]:=gp[i]/g[i];
        end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
{     writeln(' 1-© à áç¥â ¯®â®ª®¢ ');
     for i:=1 to k do
        begin
           g[i]:=gp[i]+gm[i]; tet[i]:=gp[i]/g[i];
           writeln(' gm[',i:2,']= ',gm[i]:12:8,' g[',i:2,']= ',g[i]:12:8,' tet[',i:2,']= ',tet[i]:12:8);
           if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
        end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
{     writeln(' ª®­æ¥­âà æ¨¨ 5-© ª®¬¯®­¥­âë ');
           for i:=1 to k do
              begin
                 f1:=1-cm[1,i]-cm[2,i]-cm[3,i]-cm[4,i]; f2:=1-cp[1,i]-cp[2,i]-cp[3,i]-cp[4,i];
                 f:=(f1*gm[i]+f2*gp[i])/g[i];
                 writeln(' cm[5,',i:2,']= ',f1:12:10,' cp[5,',i:2,']= ',f2:12:10,' c[5,',i:2,']= ',f:12:10);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end; }
     for i:=1 to 4 do for j:=1 to k do
         begin
            ccm[i,j]:=cm[i,j]; ccp[i,j]:=cp[i,j]; cc[i,j]:=c[i,j];
         end;
     for i:=1 to k do
         begin
            ggm[i]:=gm[i]; ggp[i]:=gp[i]; gg[i]:=g[i];
         end;

     con;

     f4:=0; for i:=1 to k do f4:=f4+g[i];
     writeln(' áã¬¬  ¯®â®ª®¢ ¯¨â ­¨ï ¯® ¨áå®¤­ë¬ ª®­æ¥­âà æ¨ï¬ - ', f4:15:10); sg:=f4;
     {f1:=1-cp[1,k]-cp[2,k]-cp[3,k]-cp[4,k]; if f1<=1e-18 then f1:=1e-12;
     f2:=1-cm[1,1]-cm[2,1]-cm[3,1]-cm[4,1]; if f2<= 1e-18 then f2:= 1e-12;
     f3:=1-ko1-ko2-ko3-ko4;
     f1:=potCM(cp[1,k],cp[2,k],cp[3,k],cp[4,k],f1); f2:=potCM(cm[1,1],cm[2,1],cm[3,1],cm[4,1],f2);
     f3:=potCM(ko1,ko2,ko3,ko4,f3);
     f4:=T*f1+T1*f2-T0*f3;
     writeln('  ¨¤¥ «ì­ ï áã¬¬  ¯®â®ª®¢ ',f4/eu:15:10,' ---- ',f4/eu/sg*100:6:3,'%');}
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;

     writeln('  ¨§¬¥­¥­¨¥ ­  § ¤ ­­ë¥ ª®­æ¥­âà æ¨¨ ¯® 1-¬ã ª®¬¯®­¥­âã ');

     km11:=0.001; kp11:=0.055; jc:=1; kpp1:=0.0085;
     hcm:=(km11-km1)/10; hcp:=(kp11-kp1)/10;

    { cgrn:=0.012; cgrv:=0.01266;}
     repeat
       km1:=km1+hcm; kp1:=kp1+hcp;
       {T:=8/8.76/3.6; TT:=2.5/8.76/3.6;
       T1:=(T00*(ko1-koo1)+TT*(kpp1-ko1)+T*(kp1-ko1))/(ko1-km1); T0:=T+TT+T1-T00;}
       T0:=200/8.76/3.6; TT:=150/8.76/3.6; T:=(T0-TT)*(ko1-km1)/(kp1-km1); T1:=T0-TT-T;
       {f:=200/167.5; f1:=-f/(1+f)*kpp1-1/(1+f)*km1+kp1; f2:=-f/(1+f)*kpp1-1/(1+f)*km1+ko1;
       f3:=-f/(1+f)*kpp1-1/(1+f)*km1+koo1; f4:=-f/(1+f)*kpp1-1/(1+f)*km1+kooo1;
       T:=(T0*f2+T00*f3+T000*f4)/f1; T1:=(T0+T00+T000-T)/(1+f); TT:=f*T1;}
       writeln('TT= ',TT:8:3,' T= ',T:8:3,' T1= ',T1:8:3,' T000= ',T000:8:3,' T00= ',T00:8:3,' T0= ',T0:8:3);
      { writeln(' T000= ',T000:12:8,' T00= ',T00:12:8,' T0= ',T0:12:8);}
       gm[1]:=T1; gp[k]:=T;
      { for i:=1 to k do
        if i<=pp then
           begin
              tau1[i]:=T1; taui1[1,i]:=T1*km1
           end; }
    {   for i:=1 to k do
        if (i>pp) and (i<=p) then
           begin
              tau1[i]:=T1-T00; taui1[1,i]:=T1*km1-T00*koo1
           end;  }
 {        else
           begin
              tau1[i]:=-T; taui1[1,i]:=-T*kp1;
           end};
       cp[1,k]:=kp1;  cm[1,1]:=km1;
       for i:=2 to k-1 do
          if i<p then
             cm[1,i]:=cm[1,i]+hcm*(ko1-cm[1,i])/(ko1-cm[1,1])
            else
             cm[1,i]:=cm[1,i]+hcp*(cp[1,i]-ko1)/(cp[1,k]-ko1);
       ip:=1; iis:=0;
       ras; jc:=jc+1
     until jc>10;

   for jc:=1 to 30 do
    begin
     if jc>=2 then
        begin
           p:=p+1; writeln(' p= ',p:3);
           {for i:=1 to k do
              if i<=p then
                begin
                   tau1[i]:=T1; taui1[1,i]:=T1*km1;
                end
              else
                begin
                   tau1[i]:=-T; taui1[1,i]:=-T*kp1;
                end;}
        end;
     ip:=1; iis:=0; ras;
     writeln(' iis= ',iis:3,' smm = ',smm:17:14); readln;
     if iis=1 then
        begin
           writeln(' ­ ©¤¥­® 1-¥ à¥è¥­¨¥: ');
           writeln('  teta= ',teta:8:6,'  T= ',T:15:12,'  T1= ',T1:15:12,' T0= ',T0:15:12);
           writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
{           writeln(' ª®­æ¥­âà æ¨¨ 1-© ª®¬¯®­¥­âë ');
           for i:=1 to k do
              begin c[1,i]:=(cm[1,i]*gm[i]+cp[1,i]*gp[i])/g[i];
                    writeln(' cm[1,',i:2,']= ',cm[1,i]:12:10,' cp[1,',i:2,']= ',cp[1,i]:12:10,' c[1,',i:2,']= ',c[1,i]:12:10);
                     if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
           writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
{           writeln(' ª®­æ¥­âà æ¨¨ 2-© ª®¬¯®­¥­âë ');
           for i:=1 to k do
              begin
                 c[2,i]:=(cm[2,i]*gm[i]+cp[2,i]*gp[i])/g[i];
                 writeln(' cm[2,',i:2,']= ',cm[2,i]:12:10,' cp[2,',i:2,']= ',cp[2,i]:12:10,' c[2,',i:2,']= ',c[2,i]:12:10);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
           write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
{    writeln(' ®â­®á¨â¥«ì­ë¥ ª®­æ¥­âà æ¨¨ 2-© ª®¬¯®­¥­âë ');
    for i:=1 to k do
       begin
          f:=cm[2,i]/cm[1,i]; f1:=cp[2,i]/cp[1,i];
          f2:=c[2,i]/c[1,i];
          writeln(' rm[2,',i:2,']= ',f:12:10,' rp[2,',i:2,']= ',f1:12:10,' r[2,',i:2,']= ',f2:12:10);
          if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
       end;
    write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
{          writeln(' ª®­æ¥­âà æ¨¨ 3-© ª®¬¯®­¥­âë ');
           for i:=1 to k do
              begin
                 c[3,i]:=(cm[3,i]*gm[i]+cp[3,i]*gp[i])/g[i];
                 writeln(' cm[3,',i:2,']= ',cm[3,i]:12:10,' cp[3,',i:2,']= ',cp[3,i]:12:10,' c[3,',i:2,']= ',c[3,i]:12:10);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
           write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
           writeln(' ª®­æ¥­âà æ¨¨ 4-© ª®¬¯®­¥­âë ');
           for i:=1 to k do
              begin
                 f1:=1-cm[1,i]-cm[2,i]-cm[3,i]; f2:=1-cp[1,i]-cp[2,i]-cp[3,i];
                 f:=(f1*gm[i]+f2*gp[i])/g[i];
                 writeln(' cm[4,',i:2,']= ',f1:12:10,' cp[4,',i:2,']= ',f2:12:10,' c[4,',i:2,']= ',f:12:10);
                      if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
           write(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
           writeln(' ¯®â®ª¨ áâã¯¥­¥© ');
           for i:=1 to k do
              begin
                 tet[i]:=gp[i]/g[i];
                 writeln(' g[',i:2,']= ',g[i]:12:8,' tet[',i:2,']= ',tet[i]:12:8);
                    if (i=22) or (i=44) or (i=66) or (i=88) or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
           f4:=0; for i:=1 to k do f4:=f4+g[i];
           writeln(' áã¬¬  ¯®â®ª®¢ ¯¨â ­¨ï ¯® ¨áå®¤­ë¬ ª®­æ¥­âà æ¨ï¬ - ', f4:15:10); sg:=f4;
           rasn; sn:=f4;
           writeln(' ª®«¨ç¥áâ¢  ¨ ¯®â®ª¨ æ¥­âà¨äã£ ');
           for i:=1 to k do
              begin
                 writeln(' n[',i:2,']= ',n1[i]:12:8,' gæ[',i:2,']= ',g[i]/n1[i]:12:8,' tet[',i:2,']= ',tet[i]:12:8);
                    if (i=22) or (i=44) or (i=66) or (i=88) or (i=110) or (i=132) then
                       begin
                          writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                       end;
              end;
           {writeln(' áã¬¬  ª¢ ¤à â®¢ à §­®áâ¥© ª®«¨ç¥áâ¢ æ¥­âà¨äã£ ', sn:12:7);}
        end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
     f1:=(cp[1,k]-ko1)/(cp[1,k]-cm[1,1]);
     writeln(' teta 1: ',f1:15:12,' 2: ',(cp[2,k]-ko2)/(cp[2,k]-cm[2,1]):15:12,' 3: ',(cp[3,k]-ko3)/(cp[3,k]-cm[3,1]):15:12);
     writeln(' ª®«¨ç¥áâ¢® á«ãç ©­ëå £¥­¥à æ¨© ¨ ¤®«ï ®â ¨­â¥à¢ «  ¢®§¬®¦­ëå ®âª«®­¥­¨©');
     readln(m,d);
     for i:=1 to k do
        begin
           nn[i]:=n1[i]; TLL:=TL; TLL_1:=TL_1;  jnn:=jn;
           gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
           ci[i]:=cm[1,i];
           cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
           ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
           ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
           cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
          if (i<>1) and (i<>k) then
              begin
                 cn[i]:=ci[i]-d*(cm[1,i]-cm[1,i-1])/2;
                 cv[i]:=ci[i]+d*(cm[1,i+1]-cm[1,i])/2;
              end
            else
              begin
                 cn[i]:=cm[1,i]; cv[i]:=cm[1,i];
              end;
        end;
     sm[1,1]:=cm[1,1]; sm[1,k]:=cm[1,k];
     gn:=1; gv:=0; gn1:=1; gv1:=0;
     for i1:=1 to m do
        begin
           bl1:=true;
           for i:=1 to k do
              begin
                 sm[2,i]:=ccm[2,i]; sm[3,i]:=ccm[3,i];
                 sm[4,i]:=ccm[4,i];
                 if (i<>1) and (i<>k) then
                    begin
                       sm[1,i]:=cn[i]+random*(cv[i]-cn[i]);
                    end;
                 f:=1+(xi[1,i]-1)*sm[1,i]+(xi[2,i]-1)*sm[2,i]+(xi[3,i]-1)*sm[3,i]+(xi[4,i]-1)*sm[4,i];
                 sp[1,i]:=xi[1,i]*sm[1,i]/f;
                if i>1 then
                   gm[i]:=(tau1[i]*sp[1,i-1]-taui1[1,i])/(sp[1,i-1]-sm[1,i]);
                   if gm[i]<=0 then bl1:=false
              end;
           if bl1 then
              begin
                 for i:=1 to k-1 do gp[i]:=gm[i+1]-tau1[i+1];
                 for i:=1 to k do
                    begin
                       g[i]:=gp[i]+gm[i]; tet[i]:=gp[i]/g[i];
                    end;
                 for i:=1 to k do
                    begin
                       cm[1,i]:=sm[1,i]; cp[1,i]:=sp[1,i];
                       cm[2,i]:=sm[2,i]; cm[3,i]:=sm[3,i];
                       cm[4,i]:=sm[4,i];
                    end;
                 ip:=1; iis:=0;
                 ras;
                 if iis=1 then
                    begin
                       for i:=1 to k do
                          begin
                             c[1,i]:=(cm[1,i]*gm[i]+cp[1,i]*gp[i])/g[i];
                             c[2,i]:=(cm[2,i]*gm[i]+cp[2,i]*gp[i])/g[i];
                             c[3,i]:=(cm[3,i]*gm[i]+cp[3,i]*gp[i])/g[i];
                             c[4,i]:=(cm[4,i]*gm[i]+cp[4,i]*gp[i])/g[i];
                             tet[i]:=gp[i]/g[i];
                          end;
                       rasn;
                       jj:=0;
                             for i:=1 to k do if jj=0 then
                                begin
                                   f:=(0.84*g[i]-gm[i])/0.00299/n1[i]; if f>0 then f:=sqrt(f);
                                   if f>36.5 then jj:=-1;
                                end;
                             if (f4<sn) and (jj=0) then
                         begin
                            sn:=f4;
                            for i:=1 to k do
                               begin
                                  nn[i]:=n1[i]; TLL:=TL; TLL_1:=TL_1;
                                  gg[i]:=g[i]; ggm[i]:=gm[i]; ggp[i]:=gp[i]; ttet[i]:=tet[i];
                                  cc[1,i]:=c[1,i]; cc[2,i]:=c[2,i]; cc[3,i]:=c[3,i];
                                  ccm[1,i]:=cm[1,i]; ccm[2,i]:=cm[2,i]; ccm[3,i]:=cm[3,i];
                                  ccp[1,i]:=cp[1,i]; ccp[2,i]:=cp[2,i]; ccp[3,i]:=cp[3,i];
                                  cc[4,i]:=c[4,i]; ccm[4,i]:=cm[4,i]; ccp[4,i]:=cp[4,i];
                               end;
                         end
                    end;
              end
        end;
   for i:=1 to k do cm[1,i]:=ccm[1,i];
   soxr;
    {huk0; huk; huk1; huk;
    huk0; huk; huk1;} hukk;
    writeln(' ­ ©¤¥­® à¥è¥­¨¥ ¯®á«¥ á«ãç ©­ëå ¨â¥à æ¨© ¨ ¬¥â®¤  •ãª  ¨ „¦¨¢á : ');
    con; rasn; snn:=sn; potok; {soxr; xxis;} {hukx};
   { for i:=1 to k do
       begin
          nn[i]:=nnn[i]; xi[1,i]:=xxi[1,i];
          gg[i]:=ggg[i]; ggm[i]:=gggm[i]; ggp[i]:=gggp[i]; ttet[i]:=tttet[i];
          cc[1,i]:=ccc[1,i]; cc[2,i]:=ccc[2,i]; cc[3,i]:=ccc[3,i];
          ccm[1,i]:=cccm[1,i]; ccm[2,i]:=cccm[2,i]; ccm[3,i]:=cccm[3,i];
          ccp[1,i]:=cccp[1,i]; ccp[2,i]:=cccp[2,i]; ccp[3,i]:=cccp[3,i];
          cc[4,i]:=ccc[4,i]; ccm[4,i]:=cccm[4,i]; ccp[4,i]:=cccp[4,i];
       end;
    con; sn:=snn; f:=0;
    for i:=1 to k do f:=f+(nn[i]-n[i])*(nn[i]-n[i]);
    writeln(' snn=',snn:12:8,' summa kvadratov N ',f:12:8); readln;

    potok;}
{     for i:=1 to 3 do
         for j:=1 to k do ft[i,j]:=ggp[j]*ccp[i,j]/gg[j]/cc[i,j];
     writeln(' áà¥§ë ¯® ª®¬¯®­¥­â ¬ ');
     for i:=1 to k do
        begin
           writeln(' ft[1,',i:2,']= ',ft[1,i]:12:10,' ft[2,',i:2,']= ',ft[2,i]:12:10,' ft[3,',i:2,']= ',ft[3,i]:12:10);
           if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
        end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
    { writeln(' ª®íää¨æ¨¥­âë à §¤¥«¥­¨ï alfa12, beta12, alfa12, beta12  ');
     for i:=1 to k do
        begin
           f1:=ccp[1,i]/ccp[2,i]/cc[1,i]*cc[2,i];  f2:=cc[1,i]/cc[2,i]/ccm[1,i]*ccm[2,i];
           f3:=ccp[1,i]/(1-ccp[1,i]-ccp[2,i])/cc[1,i]*(1-cc[1,i]-cc[2,i]);
           f4:=cc[1,i]/(1-cc[1,i]-cc[2,i])/ccm[1,i]*(1-ccm[1,i]-ccm[2,i]);
           writeln(' a12[',i:2,']= ',f1:8:5,' b12[',i:2,']= ',f2:8:5,'    a13[',i:2,']= ',f3:8:5,' b13[',i:2,']= ',f4:8:5);
           if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
        end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;}
 writeln(' ¡ « ­áë ¯® áâã¯¥­ï¬ ');
     for i:=1 to k do
              begin
                 if i=1 then
                    begin
                       f1:=ggm[2]*ccm[1,2]-ggp[1]*ccp[1,1]-ggm[1]*ccm[1,1];
                       f2:=ggm[2]*ccm[2,2]-ggp[1]*ccp[2,1]-ggm[1]*ccm[2,1];
                       f3:=ggm[2]*ccm[3,2]-ggp[1]*ccp[3,1]-ggm[1]*ccm[3,1];
                       f4:=ggm[2]*ccm[4,2]-ggp[1]*ccp[4,1]-ggm[1]*ccm[4,1];
                       f5:=ggm[2]*(1-ccm[1,2]-ccm[2,2]-ccm[3,2]-ccm[4,2])-ggp[1]*(1-ccp[1,1]-ccp[2,1]-ccp[3,1]-ccp[4,1])
                       -ggm[1]*(1-ccm[1,1]-ccm[2,1]-ccm[3,1]-ccm[4,1]);
                    end;
                 if i=k then
                    begin
                       f1:=ggp[k-1]*ccp[1,k-1]-ggp[k]*ccp[1,k]-ggm[k]*ccm[1,k];
                       f2:=ggp[k-1]*ccp[2,k-1]-ggp[k]*ccp[2,k]-ggm[k]*ccm[2,k];
                       f3:=ggp[k-1]*ccp[3,k-1]-ggp[k]*ccp[3,k]-ggm[k]*ccm[3,k];
                       f4:=ggp[k-1]*ccp[4,k-1]-ggp[k]*ccp[4,k]-ggm[k]*ccm[4,k];
                       f5:=ggp[k-1]*(1-ccp[1,k-1]-ccp[2,k-1]-ccp[3,k-1]-ccp[4,k-1])
                       -ggp[k]*(1-ccp[1,k]-ccp[2,k]-ccp[3,k]-ccp[4,k])
                       -ggm[k]*(1-ccm[1,k]-ccm[2,k]-ccm[3,k]-ccm[4,k]);
                    end;
                 if (i<>1) and (i<>k) then
                    begin
                       f1:=ggp[i-1]*ccp[1,i-1]+ggm[i+1]*ccm[1,i+1]-ggp[i]*ccp[1,i]-ggm[i]*ccm[1,i];
                       f2:=ggp[i-1]*ccp[2,i-1]+ggm[i+1]*ccm[2,i+1]-ggp[i]*ccp[2,i]-ggm[i]*ccm[2,i];
                       f3:=ggp[i-1]*ccp[3,i-1]+ggm[i+1]*ccm[3,i+1]-ggp[i]*ccp[3,i]-ggm[i]*ccm[3,i];
                       f4:=ggp[i-1]*ccp[4,i-1]+ggm[i+1]*ccm[4,i+1]-ggp[i]*ccp[4,i]-ggm[i]*ccm[4,i];
                       f5:=ggp[i-1]*(1-ccp[1,i-1]-ccp[2,i-1]-ccp[3,i-1]-ccp[4,i-1])
                       -ggp[i]*(1-ccp[1,i]-ccp[2,i]-ccp[3,i]-ccp[4,i])
                       -ggm[i]*(1-ccm[1,i]-ccm[2,i]-ccm[3,i]-ccm[4,i]);
                       f5:=f5+ggm[i+1]*(1-ccm[1,i+1]-ccm[2,i+1]-ccm[3,i+1]-ccm[4,i+1])
                    end;
                 if i=ppp then
                    begin
                       f1:=f1+T000*kooo1; f2:=f2+T000*kooo2; f3:=f3+T000*kooo3; f4:=f4+T000*kooo4;
                       f5:=f5+T000*(1-kooo1-kooo2-kooo3-kooo4)
                    end;
                 if i=pp then
                    begin
                       f1:=f1+T00*koo1; f2:=f2+T00*koo2; f3:=f3+T00*koo3; f4:=f4+T00*koo4;
                       f5:=f5+T00*(1-koo1-koo2-koo3-koo4)
                    end;
                 if i=p then
                    begin
                       f1:=f1+T0*ko1; f2:=f2+T0*ko2; f3:=f3+T0*ko3; f4:=f4+T0*ko4; f5:=f5+T0*(1-ko1-ko2-ko3-ko4)
                    end;
                  if i=jnn then
                    begin
                       f1:=f1-TLL_1*ccp[1,i-1]; f2:=f2-TLL_1*ccp[2,i-1]; f3:=f3-TLL_1*ccp[3,i-1]; f4:=f4-TLL_1*ccp[4,i-1];
                       f5:=f5-TLL_1*(1-ccp[1,i-1]-ccp[2,i-1]-ccp[3,i-1]-ccp[4,i-1])
                    end;
                if i=jnn+1 then
                    begin
                       f1:=f1-TLL*ccp[1,i-1]; f2:=f2-TLL*ccp[2,i-1]; f3:=f3-TLL*ccp[3,i-1]; f4:=f4-TLL*ccp[4,i-1];
                       f5:=f5-TLL*(1-ccp[1,i-1]-ccp[2,i-1]-ccp[3,i-1]-ccp[4,i-1])
                    end;
                 writeln(' ¡1[',i:2,']= ',f1:17:14,' ¡2[',i:2,']= ',f2:17:14,' ¡3[',i:2,']= ',f3:17:14,
                 ' ¡4[',i:2,']= ',f4:17:14,' ¡5[',i:2,']= ',f5:17:14);
                 if (i=11) or (i=22) or (i=33) or (i=44) or (i=55) or (i=66) or (i=77) or (i=88)  or (i=99) then
                    begin
                       writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
                    end;
              end;
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
     f1:=1-ccp[1,k]-ccp[2,k]-ccp[3,k]-ccp[4,k]; f2:=1-ccm[1,1]-ccm[2,1]-ccm[3,1]-ccm[4,1]; f3:=1-ko1-ko2-ko3-ko4;
     writeln('ª®­æ¥­âà æ¨¨ 5-© ( ¡á.) ',f1:15:12,f2:15:12,f3:15:12); if f1<0 then f1:=1e-12;
     f1:=potCM(ccp[1,k],ccp[2,k],ccp[3,k],ccp[4,k],f1); f2:=potCM(ccm[1,1],ccm[2,1],ccm[3,1],ccm[4,1],f2);
     f3:=potCM(ko1,ko2,ko3,ko4,f3);
     f4:=T*f1+T1*f2-T0*f3;
     writeln(' T0 g/c= ',T0:12:6,' T= ',T:12:6,' TT= ',TT:12:6,' T1= ',T1:12:6);
     writeln('  ¨¤¥ «ì­ ï áã¬¬  ¯®â®ª®¢ ',f4/eu:15:10,' ---- ',f4/eu/sg*100:6:3,'%');
     writeln('  ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦ âì ENTER  '); readln;
     writeln('  à §¤¥«¨â¥«ì­ ï á¯®á®¡­®áâì áâã¯¥­¥© (­®¬.¢¥«., ã¤¥«.æ¥­âà.,íää.,Š„)');
     f5:=0;
     for i:=1 to k do
        begin
           f1:=1-ccp[1,i]-ccp[2,i]-ccp[3,i]-ccp[4,i]; if f1<0 then f1:=1e-12;
           f2:=1-ccm[1,i]-ccm[2,i]-ccm[3,i]-ccm[4,i]; if f2<0 then f2:=1e-12;
           f3:=1-cc[1,i]-cc[2,i]-cc[3,i]-cc[4,i]; if f3<0 then f3:=1e-12;
           e[i]:=ggp[i]*potCM(ccp[1,i],ccp[2,i],ccp[3,i],ccp[4,i],f1);
           e[i]:=e[i]+ggm[i]*potCM(ccm[1,i],ccm[2,i],ccm[3,i],ccm[4,i],f2);
           e[i]:=e[i]-gg[i]*potCM(cc[1,i],cc[2,i],cc[3,i],cc[4,i],f3); f5:=f5+e[i];
           err;
           writeln('i= ',i:3,e[i]:15:10,e[i]/nn[i]:15:10,ee[i]:15:10,ee[i]/e[i]*100:15:10);
           if (i=22) or (i=44) or (i=66) or (i=88)  or (i=110) or (i=132) then
              begin
                 writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
              end;
        end;
     writeln('  áã¬¬  à §¤¥«. á¯®á®¡. áâã¯¥­¥© ',f5:15:10);
     writeln('  à §¤¥«¨â¥«ì­ ï á¯®á®¡­®áâì ª áª ¤  ', f4:15:10,' ---- ',f4/f5*100:6:3,'%');
     writeln(' ¤«ï ¯à®¤®«¦¥­¨ï ­ ¦¬¨â¥ ENTER '); readln;
     f1:=T0*ko1+T00*koo1+T000*kooo1-T*ccp[1,k]-T1*ccm[1,1]-ccp[1,jnn-1]*TLL_1-ccp[1,jnn]*TLL;
     f2:=T0*ko2+T00*koo2+T000*kooo2-T*ccp[2,k]-T1*ccm[2,1]-ccp[2,jnn-1]*TLL_1-ccp[2,jnn]*TLL;
     f3:=T0*ko3+T00*koo3+T000*kooo3-T*ccp[3,k]-T1*ccm[3,1]-ccp[3,jnn-1]*TLL_1-ccp[3,jnn]*TLL;
     f4:=T0*ko4+T00*koo4+T000*kooo4-T*ccp[4,k]-T1*ccm[4,1]-ccp[4,jnn-1]*TLL_1-ccp[4,jnn]*TLL;
     f:=T0*(1-ko1-ko2-ko3-ko4)-T*(1-ccp[1,k]-ccp[2,k]-ccp[3,k]-ccp[4,k])-T1*(1-ccm[1,1]-ccm[2,1]-ccm[3,1]-ccm[4,1]);
     f:=f+T00*(1-koo1-koo2-koo3-koo4)+T000*(1-kooo1-kooo2-kooo3-kooo4);
     f:=f-(1-ccp[1,jnn-1]-ccp[2,jnn-1]-ccp[3,jnn-1]-ccp[4,jnn-1])*TLL_1-(1-ccp[1,jnn]-ccp[2,jnn]-ccp[3,jnn]-ccp[4,jnn])*TLL;
     writeln('  ¡ « ­áë ¯® 1,2,3,4,5,-© ª®¬¯®­¥­â ¬: ',f1:17:14,f2:17:14,f3:19:16,f4:17:14,f:17:14);

     writeln('  ¤«ï ®ª®­ç ­¨ï ­ ¦ âì ENTER  '); readln;


    end;

end.
