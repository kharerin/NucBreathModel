#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define L2 L1*57	    
int main() {
        FILE *f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9;
        char c[10];
        long int i=0,i1,i2,j,k,posit,ks[2],l=147,L=L1,dw0,sze,sze1,amin,amax,lavg,dle,nuc,ksm,ks1,t,jfn,jb,jbxn,h,w,b1,b2;
        double p[L1],p1,p0,pot[L1],Es[L1],Va,gamaN=1,CN,mu,E0,Zf[L1],Zb[L1],Zfor,Zbor,Pn[L1],O[L1],Zfsum,Zbsum,pdO,lavgx,le;
        double Osum,Osum1,Osum2,Osum3,nucavg,nucavg1,nucavg2,nucavg3,On[L1],On1[L1],On2[L1],On3[L1],Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b;
        double Omd,sO,Pnsum,sumOmd,pO,d1,d2,sumd,sump,lx,lm;
        static double Eij[L1];
                     
        
        f9=fopen("./occ_0p525.txt", "w");
        if(f9==NULL) {
          printf("quit \n");
          exit(0);
        } 
        f0=fopen("./occ_0p575.txt", "w");
        if(f0==NULL) {
          printf("quit \n");
          exit(0);
        }        
        f1=fopen("./occ_0p625.txt", "w");
        if(f1==NULL) {
          printf("quit \n");
          exit(0);
        }
        f2=fopen("./occ_0p675.txt", "w");
        if(f2==NULL) {
          printf("quit \n");
          exit(0);
        }
        f3=fopen("./occ_0p725.txt", "w");
        if(f3==NULL) {
          printf("quit \n");
          exit(0);
        }
        
        f4=fopen("./d1_mu_E0d.txt", "w");
        if(f4==NULL) {
          printf("quit \n");
          exit(0);
        }
        f5=fopen("./occ_0p775.txt", "w");
        if(f5==NULL) {
          printf("quit \n");
          exit(0);
        }
        f6=fopen("./occ_0p825.txt", "w");
        if(f6==NULL) {
          printf("quit \n");
          exit(0);
        }
        f7=fopen("./occ_0p875.txt", "w");
        if(f7==NULL) {
          printf("quit \n");
          exit(0);
        }
        f8=fopen("./rms_0.txt", "w");
        if(f8==NULL) {
          printf("quit \n");
          exit(0);
        }
       
       amax=147; amin=91; nuc=1; dle=amax-amin; le=0; i=0;  
      /* for(h=-250;h<250;h=h+2) {  
          E0=(double)h/(10*(147-91)); 
          printf("le %lf i %ld E0 %lf \n",le,i,E0);
          if(E0>0 || E0<0) {
            le=le+amin+dle*(exp(-E0*dle)/(exp(-E0*dle)-1))+1/E0;
            printf("le %lf i %ld \n",le,i); 
            i=i+1;
          }
       }
       lavgx=le/(double)i; 
       if((lavgx-floor(lavgx))>0.5) {
         lavg=(int)(floor(lavgx)+1);
       }
       else {
         lavg=(int)(floor(lavgx));
       } */  
       Va=0; ksm=lavg;                     
       for(h=-250;h<250;h=h+1) {  
          E0=(double)h/(10*(147-91)); 
          if(E0>0 || E0<0) {
            le=amin+dle*(exp(-(-1)*E0*dle)/(exp(-(-1)*E0*dle)-1))+1/((-1)*E0);
            printf("le %lf h %ld E0 %lf \n",le,h,(-1)*E0);
          }
          lavgx=le; 
          if((lavgx-floor(lavgx))>0.5) {
            lavg=(int)(floor(lavgx)+1);
          }
          else {
            lavg=(int)(floor(lavgx));
          } 
          ksm=lavg;
          for(w=0;w<500;w=w+1) {
             mu=-(double)(300-w)/5;CN=exp(mu); 
             Osum=0; Osum1=0; Osum2=0; Osum3=0;       //b1=(4000-floor(w/2)); b2=(4000+floor(w/2));                          
             for(i=0;i<L1;i++) {
                Eij[i]=-E0*lavg; //Eij[i][t]=gamaN*pot[i]+Va;                
                Zf[i]=0; Zb[i]=0; Pn[i]=0; O[i]=0; On[i]=0; On1[i]=0; On2[i]=0; On3[i]=0;
             }
            
             for(i=0;i<(ksm-1);i++) {
                Zf[i]=0;
                Zb[L1-ksm+1+i]=0;
             }
             for(i=0;i<(L1-ksm+1);i++) { // DNA length
                Zbsum=0; Zfsum=0; 
                for(t=0;t<1;t++) { // sub nucleosome
                   jfn=i+ksm-lavg; // position of TF, nuc{
                   if(jfn>=0) {                      // forward nuc
                     if(jfn==0) {
                       Zfor=exp(-Zf[i+ksm-2])*CN*exp(-Eij[jfn]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
                     }
                     else {
                       Zfor=exp(Zf[jfn-1]-Zf[i+ksm-2])*CN*exp(-Eij[jfn]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
                     }
                   }  
                   else {
                     Zfor=0;
                   }
                   
                   jb=L1-ksm-i;                      // backward nuc
                   jbxn=jb+lavg;
                   if((ksm+i)>=lavg) {
                     if(jbxn>(L1-1)) {
                       Zbor=exp(-Zb[L1-ksm-i+1])*CN*exp(-Eij[jb]);
                     }
                     else {
                       Zbor=exp(Zb[jbxn]-Zb[L1-ksm-i+1])*CN*exp(-Eij[jb]);
                     }
                   }
                   else {
                     Zbor=0;
                   }
                   Zfsum=Zfsum+Zfor;
                   Zbsum=Zbsum+Zbor;
                }
               Zf[i+ksm-1]=Zf[i+ksm-2] + log(1+Zfsum); // forward sums
               Zb[L1-ksm-i]=Zb[L1-ksm-i+1] + log(1+Zbsum); // backward sums
             }
             sO=0; pO=0; lx=0;
             for(t=0;t<1;t++) {  // subnucleosome
                for(i=0;i<L1;i++) {
                   Pn[i]=0;O[i]=0;
                }
               Pnsum=0; sumOmd=0;          
               Pn[0]=exp(Zb[0+lavg]-Zf[L1-1])*CN*exp(-Eij[0]); Pnsum=Pnsum+Pn[0];                 
               for(i=1;i<=(L1-lavg-1);i++) {
                  Pn[i]=exp(Zf[i-1]+Zb[i+lavg]-Zf[L1-1])*CN*exp(-Eij[i]); Pnsum=Pnsum+Pn[i];            
               }
               Pn[L1-lavg]=exp(Zf[L1-lavg-1]-Zf[L1-1])*CN*exp(-Eij[L1-lavg]); Pnsum=Pnsum+Pn[L1-lavg]; pdO=Pnsum; pO=pO+pdO; 
               /*if((pdO[t]-floor(pdO[t]))>0.5) {
                 lx=lx+(floor(pdO[t])+1)*soft[t][0];
               }
               else {
                 lx=lx+floor(pdO[t])*soft[t][0];
               }*/
               lx=lx+pdO*lavg;
               ks1=lavg;
               if(ks1>0) {
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
                  
               Omd=sumOmd/L1; sO=sO+Omd;                
             }        
             
             lm=0;             
             if(pO>0) {
               lm=lx/pO;
             }
             for(i=0;i<L1;i++) {
                Osum=Osum+On[i]; 
             }
             
             sumd=0; sump=0;
             for(t=0;t<1;t++) {
                sumd=sumd+pow(Omd-(sO/nuc),2.0); sump=sump+pow(pdO-(pO/nuc),2.0);
             }
             d1=sqrt(sumd/nuc); d2=sqrt(sump/nuc);
             
             if((Osum/L1)>0.524 && (Osum/L1)<0.526) {
                fprintf(f9,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             if((Osum/L1)>0.574 && (Osum/L1)<0.576) {
                fprintf(f0,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             if((Osum/L1)>0.624 && (Osum/L1)<0.626) {
                fprintf(f1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             if((Osum/L1)>0.674 && (Osum/L1)<0.676) {
                fprintf(f2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.724 && (Osum/L1)<0.726) {
                fprintf(f3,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola2 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.774 && (Osum/L1)<0.776) {
                fprintf(f5,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola3 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.824 && (Osum/L1)<0.826) {
                fprintf(f6,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.874 && (Osum/L1)<0.876) {
                fprintf(f7,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola2 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((d1>0 && d1<0.0001) && ((Osum/L1)>0.1)){
                fprintf(f8,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola4 E0d %lf mu %lf \n",E0*(147-91),mu);
             }

             fprintf(f4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); 
             //fprintf(f5,"%ld %ld %lf \n",h,w,den1b); fprintf(f6,"%ld %ld %lf \n",h,w,den2b); fprintf(f7,"%ld %ld %lf \n",h,w,den3b);  
         } 
         printf("ZfL %lf nuc %ld ksm %ld lx %lf pO %lf lm %lf h %ld w %ld lavg %ld \n",Zf[L-1],nuc,ksm,lx,pO,lm,h,w,lavg);  
      }
       /* f8=fopen("./E_map.txt", "w");
        if(f8==NULL) {
          printf("quit \n");
          exit(0);
        }
        for(i=0;i<L1;i++) { 
          for(j=0;j<nuc;j++) {
             fprintf(f8,"%ld %ld %lf \n",i,j,Eij[i][j]);
          }
        }*/
        fclose(f0); fclose(f1); fclose(f2); fclose(f3); fclose(f4); fclose(f5); fclose(f6); fclose(f7); fclose(f8); fclose(f9); return(0); 
}


