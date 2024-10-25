#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define Ns 57 // 147-91+1, min size 10
#define L2 L1*Ns
#define L3 545*Ns
#define n1 3
   
double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1]);
double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1]); 
int main() {
        FILE *f4;
        char c[10];
        long int i=0,i1,i2,j,k,posit,ks[2],l=147,L=L1,dw0,sze,sze1,amin,amax,nuc,ksm,lavg,dle,ks1,t,jfn,jb,jbxn,h=1,w,b1,b2,js,a,cnt;
        double p[L1],p1,p0,pot[L1],Va,gamaN=1,CN,mu,E0,Zf[L1],Zb[L1],Pn[L1],O[L1],Zfsum,Zbsum,pOx,pdO,Eh[6],lavgx,le;
        double Osum,Osum1,Osum2,Osum3,nucavg,nucavg1,nucavg2,nucavg3,On[L1],On1[L1],On2[L1],On3[L1],Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b;
        double Omd,sO,Pnsum,sumOmd,pO,d1,d2,sumd,sump,lx,lm;
        static double Eij[L1];
        
        
        f4=fopen("./den_mu_E0d.txt", "w");
        if(f4==NULL) {
          printf("quit \n");
          exit(0);
        }

       amax=147; amin=91; nuc=1; dle=amax-amin; le=0; i=0; a=0;
       /*for(h=-250;h<250;h=h+2) {  
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
       }*/   
       Va=0; ksm=lavg;                       
       Eh[0]=-10; Eh[1]=0.001; Eh[2]=5; Eh[3]=20; Eh[4]=40; Eh[5]=60;
       cnt=0;                 
       for(h=0;h<6;h++) {  
          E0=Eh[h]/(147-1); 
          if(E0>0 || E0<0) {
            le=amin+dle*(exp(-(-1)*E0*dle)/(exp(-(-1)*E0*dle)-1))+1/((-1)*E0);
            printf("le %lf cnt %ld E0 %lf \n",le,cnt,(-1)*E0); cnt=cnt+1;
          }
          lavgx=le; 
          if((lavgx-floor(lavgx))>0.5) {
            lavg=(int)(floor(lavgx)+1);
          }
          else {
            lavg=(int)(floor(lavgx));
          }
          ksm=lavg;   
          for(w=0;w<575;w++) {
             mu=-(double)(375-w)/5;CN=exp(mu);
             Osum=0; Osum1=0; Osum2=0; Osum3=0;
             for(i=0;i<L1;i++) {
                for(t=0;t<nuc;t++) {
                   Eij[i]=-E0*lavg;//Eij[i][t]=gamaN*pot[i]+Va;
                }
                Zf[i]=0; Zb[i]=0; Pn[i]=0; O[i]=0; On[i]=0; On1[i]=0; On2[i]=0; On3[i]=0; //Es[i]=Eij[i][nuc-1];
             }
             
             for(i=0;i<(L1-ksm+1);i++) { // DNA length
                Zbsum=0; Zfsum=0; 
                for(t=0;t<nuc;t++) { // sub nucleosome
                   sze=lavg+a;
                   Zfsum=Zfsum+Zfor(i,ksm,sze,t,CN,Zf,Eij); // forward nuc                    
                   Zbsum=Zbsum+Zbor(i,ksm,sze,t,CN,Zb,Eij); // backward nuc                                                  
                }                                       
                Zf[i+ksm-1]=Zf[i+ksm-2] + log(1+Zfsum); // forward sums
                Zb[L1-ksm-i]=Zb[L1-ksm-i+1] + log(1+Zbsum); // backward sums
             }
                          
             sO=0; pO=0; lx=0;
             for(t=0;t<nuc;t++) {  // subnucleosome
                for(i=0;i<L1;i++) {
                   Pn[i]=0;O[i]=0;
                }
               sze=lavg+a; Pnsum=0; sumOmd=0;          
               Pn[0]=exp(Zb[0+sze]-Zf[L1-1])*CN*exp(-Eij[0]); Pnsum=Pnsum+Pn[0];          
               for(i=1;i<=(L1-sze-1);i++) {
                  Pn[i]=exp(Zf[i-1]+Zb[i+sze]-Zf[L1-1])*CN*exp(-Eij[i]); Pnsum=Pnsum+Pn[i]; 
               }
               Pn[L1-sze]=exp(Zf[L1-sze-1]-Zf[L1-1])*CN*exp(-Eij[L1-lavg]); Pnsum=Pnsum+Pn[L1-sze]; pdO=Pnsum; pO=pO+pdO; 
              /* if((pdO[t]-floor(pdO[t]))>0.5) {
                 lx=lx+(floor(pdO[t])+1)*soft[t][0];
               }
               else {
                 lx=lx+floor(pdO[t])*soft[t][0];
               }*/
               lx=lx+pdO*lavg;
               ks1=lavg;
               //if(ks1>0 && ks1<110) {
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
              // }
              
               Omd=sumOmd/L1; sO=sO+Omd;                
             }
             
             lm=0;
             /*if((pO-floor(pO))>0.5) {
               lm=lx/(floor(pO)+1);
             }
             else {
               if(floor(pO)>0) {
                 lm=lx/floor(pO);
               }
             }*/
             if(pO>0) {
               lm=lx/pO;
             }
                
             for(i=0;i<L1;i++) {       
                Osum=Osum+On[i]; Osum1=Osum1+On1[i]; Osum2=Osum2+On2[i]; Osum3=Osum3+On3[i];  
             }        
             sumd=0; sump=0;
             for(t=0;t<nuc;t++) {
                sumd=sumd+pow(Omd-(sO/nuc),2.0); sump=sump+pow(pdO-(pO/nuc),2.0);
             }
             d1=sqrt(sumd/nuc); d2=sqrt(sump/nuc);           

             fprintf(f4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd-Omd),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); 
             //fprintf(f5,"%ld %ld %lf \n",h,w,den1b); fprintf(f6,"%ld %ld %lf \n",h,w,den2b); fprintf(f7,"%ld %ld %lf \n",h,w,den3b);  
             
             printf("ZfL %lf nuc %ld ksm %ld lx %lf pO %lf lm %lf h %ld w %ld \n",Zf[L-1],nuc,ksm,lx,pO,lm,h,w);
         } 
         printf("h %ld done! \n",h);  
      }
      fclose(f4); return(0); 
}


double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1]) { //i,ksm,tfsize,t,CN,(Zf or Zb),E1
       double Zx; long int jfn;
       jfn=num0+num1-num2;               // jfn=i+ksm-soft[t][0]; // position of TF, nuc 
       if(jfn>=0) {                      // forward nuc
         if(jfn==0) {
           Zx=exp(-num5[num0+num1-2])*num4*exp(-num6[jfn]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
         else {
           Zx=exp(num5[jfn-1]-num5[num0+num1-2])*num4*exp(-num6[jfn]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
       }  
       else {
         Zx=0;
       }
       return Zx;
}

double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1]) { //i,ksm,tfsize,t,CN,(Zf or Zb),E1
       double Zx; long int jfn,jb,jbxn;
       jfn=num0+num1-num2;               // jfn=i+ksm-soft[t][0]; // position of TF, nuc        
       jb=L1-num1-num0;                      // backward nuc
       jbxn=jb+num2;
       if(jfn>=0) {                      
         if(jfn==0) {
           Zx=exp(-num5[L1-num1-num0+1])*num4*exp(-num6[jb]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
         else {
           Zx=exp(num5[jbxn]-num5[L1-num1-num0+1])*num4*exp(-num6[jb]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
       }  
       else {
         Zx=0;
       }
       return Zx;
}


