#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define L2 L1*57	    
int main() {
        FILE *f4;
        long int i=0,i1,i2,j,k,posit,ks[2],l=147,L=L1,dw0,sze,sze1,soft[1000][3],nuc,ksm,ks1,t,jfn,jb,jbxn,h,w,b1,b2;
        double p[L1],p1,p0,pot[L1],Es[L1],Va,gamaN=1,CN,mu,E0,Zf[L1],Zb[L1],Zfor,Zbor,Pn[L1],O[L1],Zfsum,Zbsum,pdO[57];
        double Osum,Osum1,Osum2,Osum3,nucavg,nucavg1,nucavg2,nucavg3,On[L1],On1[L1],On2[L1],On3[L1],Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b;
        double Omd[57],sO,Pnsum,sumOmd,pO,d1,d2,sumd,sump;
        static double Eij[L1][57];
                     
        
        f4=fopen("./d1_mu_E0d.txt", "w");
        if(f4==NULL) {
          printf("quit \n");
          exit(0);
        }
             Va=0;ks[0]=147;ks[1]=17; nuc=56+1; ksm=91;                   
             for(i=0;i<nuc;i++) {
                soft[i][0]=91+i; soft[i][1]=147-soft[i][0]; soft[i][2]=0;
             }
              
       for(h=-250;h<250;h++) {  
          E0=(double)h/(10*(147-91)); 
          for(w=0;w<500;w++) {
             mu=-(double)(300-w)/5;CN=exp(mu); 
             Osum=0; Osum1=0; Osum2=0; Osum3=0;       //b1=(4000-floor(w/2)); b2=(4000+floor(w/2));                          
             for(i=0;i<L1;i++) {
                for(t=0;t<nuc;t++) {
                   Eij[i][t]=-E0*soft[t][0]; //Eij[i][t]=gamaN*pot[i]+Va;
                }
                Es[i]=Eij[i][56]; Zf[i]=0; Zb[i]=0; Pn[i]=0; O[i]=0; On[i]=0; On1[i]=0; On2[i]=0; On3[i]=0;
             }
            
             for(i=0;i<(ksm-1);i++) {
                Zf[i]=0;
                Zb[L1-ksm+1+i]=0;
             }
             for(i=0;i<(L1-ksm+1);i++) { // DNA length
                Zbsum=0; Zfsum=0; 
                for(t=0;t<(nuc-1);t++) { // sub nucleosome
                   jfn=i+ksm-soft[t][0]; // position of TF, nuc{
                   if(jfn>=0) {                      // forward nuc
                     if(jfn==0) {
                       Zfor=exp(-Zf[i+ksm-2])*CN*exp(-Eij[jfn][t]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
                     }
                     else {
                       Zfor=exp(Zf[jfn-1]-Zf[i+ksm-2])*CN*exp(-Eij[jfn][t]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
                     }
                   }  
                   else {
                     Zfor=0;
                   }
                   
                   jb=L1-ksm-i;                      // backward nuc
                   jbxn=jb+soft[t][0];
                   if((ksm+i)>=soft[t][0]) {
                     if(jbxn>(L1-1)) {
                       Zbor=exp(-Zb[L1-ksm-i+1])*CN*exp(-Eij[jb][t]);
                     }
                     else {
                       Zbor=exp(Zb[jbxn]-Zb[L1-ksm-i+1])*CN*exp(-Eij[jb][t]);
                     }
                   }
                   else {
                     Zbor=0;
                   }
                   Zfsum=Zfsum+Zfor;
                   Zbsum=Zbsum+Zbor;
                }
                                
                if((i+ksm)>=l) {                  // forwards nucleosome
                  if((i+ksm)==l) {
                    Zfor=exp(-Zf[i+ksm-2])*CN*exp(-Es[i+ksm-l]); 
                  }
                  else {
                    Zfor=exp(Zf[i+ksm-l-1]-Zf[i+ksm-2])*CN*exp(-Es[i+ksm-l]); 
                  }
                }
                else {
                  Zfor=0;
                }          
           
               if((ksm+i)>=l) {               // backwards nucleosome
                 if((L1-ksm-i+l)>(L1-1)) {
                   Zbor=exp(-Zb[L1-ksm-i+1])*CN*exp(-Es[L1-ksm-i]); 
                 }
                 else {
                   Zbor=exp(Zb[L1-ksm-i+l]-Zb[L1-ksm-i+1])*CN*exp(-Es[L1-ksm-i]);
                 }
               }
               else {
                 Zbor=0;
               }
               Zfsum=Zfsum+Zfor;
               Zbsum=Zbsum+Zbor;
               Zf[i+ksm-1]=Zf[i+ksm-2] + log(1+Zfsum); // forward sums
               Zb[L1-ksm-i]=Zb[L1-ksm-i+1] + log(1+Zbsum); // backward sums
             }
             sO=0; pO=0;
             for(t=0;t<(nuc-1);t++) {  // subnucleosome
                for(i=0;i<L1;i++) {
                   Pn[i]=0;O[i]=0;
                }
               Pnsum=0; sumOmd=0;          
               Pn[0]=exp(Zb[0+soft[t][0]]-Zf[L1-1])*CN*exp(-Eij[0][t]); Pnsum=Pnsum+Pn[0];                 
               for(i=1;i<=(L1-soft[t][0]-1);i++) {
                  Pn[i]=exp(Zf[i-1]+Zb[i+soft[t][0]]-Zf[L1-1])*CN*exp(-Eij[i][t]); Pnsum=Pnsum+Pn[i];            
               }
               Pn[L1-soft[t][0]]=exp(Zf[L1-soft[t][0]-1]-Zf[L1-1])*CN*exp(-Eij[L1-soft[t][0]][t]); Pnsum=Pnsum+Pn[L1-soft[t][0]]; pdO[t]=Pnsum; pO=pO+pdO[t];
               ks1=soft[t][0];
               if(ks1>0 && ks1<110) {
                 for(i=0;i<L1;i++) {
                    if(i<ks1-1) {
                      if(i==0) {
                        O[i]=Pn[i];
                      }
                      else {
                        for(j=0;j<=i;j++) {
                           O[i]=O[i]+Pn[j];
                        }
                      }
                    }
                    else if(i>(L1-ks1)) {
                      if(i==L1-1) {
                        O[i]=Pn[L1-ks1];
                      }
                      else {
                        for(j=i-ks1+1;j<=L1-ks1;j++) {
                           O[i]=O[i]+Pn[j];
                        }
                      }
                    }
                    else {
                      for(j=i-ks1+1;j<=i;j++) {
                         O[i]=O[i]+Pn[j];
                      }
                    }
                 }
                 for(i=0;i<L1;i++) {       
                    On[i]=On[i]+O[i]; On1[i]=On1[i]+O[i]; sumOmd=sumOmd+O[i]; 
                 }                 
               }
               
               if(ks1>109 && ks1<129) {
                 for(i=0;i<L1;i++) {
                    if(i<ks1-1) {
                      if(i==0) {
                        O[i]=Pn[i];
                      }
                      else {
                        for(j=0;j<=i;j++) {
                           O[i]=O[i]+Pn[j];
                        }
                      }
                    }
                    else if(i>(L1-ks1)) {
                      if(i==L1-1) {
                        O[i]=Pn[L1-ks1];
                      }
                      else {
                        for(j=i-ks1+1;j<=L1-ks1;j++) {
                           O[i]=O[i]+Pn[j];
                        }
                      }
                    }
                    else {
                      for(j=i-ks1+1;j<=i;j++) {
                         O[i]=O[i]+Pn[j];
                      }
                    }
                 }
                 for(i=0;i<L1;i++) {       
                    On[i]=On[i]+O[i]; On2[i]=On2[i]+O[i]; sumOmd=sumOmd+O[i];
                 }                 
               } 
               
               if(ks1>128 && ks1<147) {
                 for(i=0;i<L1;i++) {
                    if(i<ks1-1) {
                      if(i==0) {
                        O[i]=Pn[i];
                      }
                      else {
                        for(j=0;j<=i;j++) {
                           O[i]=O[i]+Pn[j];
                        }
                      }
                    }
                    else if(i>(L1-ks1)) {
                      if(i==L1-1) {
                        O[i]=Pn[L1-ks1];
                      }
                      else {
                        for(j=i-ks1+1;j<=L1-ks1;j++) {
                           O[i]=O[i]+Pn[j];
                        }
                      }
                    }
                    else {
                      for(j=i-ks1+1;j<=i;j++) {
                         O[i]=O[i]+Pn[j];
                      }
                    }
                 }
                 for(i=0;i<L1;i++) {       
                    On[i]=On[i]+O[i]; On3[i]=On3[i]+O[i]; sumOmd=sumOmd+O[i];
                 }                 
               }  
               Omd[t]=sumOmd/L1; sO=sO+Omd[t];                
             }
        
             for(i=0;i<L1;i++) { // full nuc
                Pn[i]=0;O[i]=0;
             }
             Pnsum=0; sumOmd=0;        
             Pn[0]=exp(Zb[0+l]-Zf[L1-1])*CN*exp(-Es[0]); Pnsum=Pnsum+Pn[0];
             Pn[L1-l]=exp(Zf[L1-l-1]-Zf[L1-1])*CN*exp(-Es[L1-l]); Pnsum=Pnsum+Pn[L1-l];
             for(i=1;i<=(L1-l-1);i++) {
                Pn[i]=exp(Zf[i-1]+Zb[i+l]-Zf[L1-1])*CN*exp(-Es[i]); Pnsum=Pnsum+Pn[i];
             }
             for(i=0;i<L1;i++) {
                if(i<ks[0]-1) {
                  if(i==0) {
                    O[i]=Pn[i];
                  }
                  else {
                    for(j=0;j<=i;j++) {
                       O[i]=O[i]+Pn[j];
                    }
                  }
                }
                else if(i>(L1-ks[0])) {
                  if(i==L1-1) {
                    O[i]=Pn[L1-ks[0]];
                  }
                  else {
                    for(j=i-ks[0]+1;j<=L1-ks[0];j++) {
                       O[i]=O[i]+Pn[j];
                    }
                  }
                }
                else {
                  for(j=i-ks[0]+1;j<=i;j++) {
                     O[i]=O[i]+Pn[j];
                  }
                }
             }            
             for(i=0;i<L1;i++) {       
                On[i]=On[i]+O[i]; On3[i]=On3[i]+O[i]; Osum=Osum+On[i]; Osum1=Osum1+On1[i]; Osum2=Osum2+On2[i]; Osum3=Osum3+On3[i]; sumOmd=sumOmd+O[i];
             }
             Omd[56]=sumOmd/L1; sO=sO+Omd[56]; pdO[56]=Pnsum; pO=pO+pdO[56]; 
             sumd=0; sump=0;
             for(t=0;t<nuc;t++) {
                sumd=sumd+pow(Omd[t]-(sO/nuc),2.0); sump=sump+pow(pdO[t]-(pO/nuc),2.0);
             }
             d1=sqrt(sumd/nuc); d2=sqrt(sump/nuc);
                                      
             fprintf(f4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,Osum3/L1,Osum2/L1,Osum1/L1); 

         } 
         printf("ZfL %lf nuc %ld ksm %ld h %ld w %ld \n",Zf[L-1],nuc,ksm,h,w);  
      }

      fclose(f4); return(0); 
}


