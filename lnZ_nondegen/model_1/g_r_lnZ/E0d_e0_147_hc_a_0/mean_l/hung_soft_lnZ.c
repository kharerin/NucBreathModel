#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define Ns 57 // 147-91+1, min size 10
#define L2 L1*Ns
#define L3 545*Ns
#define n1 3
   
double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]);
double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]); 
int main() {
        FILE *f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9,*f15,*f16,*f17;
        char c[10];
        long int i=0,i1,i2,dx,j,k,posit,ks[2],l=147,L=L1,dw0,sze,sze1,soft[1000][3],nuc,ksm,ks1,t,jfn,jb,jbxn,h=1,w,b1,b2,js,a,amin,pk,tn,hx;
        double p[L1],p1,p0,pot[L1],Es[L1],Va,gamaN=1,CN,mu,E0,Zf[L1],Zb[L1],Pn[L1],O[L1],Zfsum,Zbsum,pOx,pdO[Ns],GN0,sumG,dg;
        double Osum,Osum1,Osum2,Osum3,nucavg,nucavg1,nucavg2,nucavg3,On[L1],On1[L1],On2[L1],On3[L1],Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b;
        double Omd[Ns],sO,Pnsum,Pnsum1,sumOmd,pO,pO1,d1,d2,Zm,gsum,gsum1,sumd,sump,lx,lm;
        static double Eij[L1][Ns],gn[2000][Ns],Gn[2000],GN[2000],Gp[5][2];
        
        
        
        f4=fopen("./d1_mu_E0d.txt", "w");
        if(f4==NULL) {
          printf("quit \n");
          exit(0);
        }
        f6=fopen("./occ_0p7.txt", "w");
        if(f6==NULL) {
          printf("quit \n");
          exit(0);
        }
        f7=fopen("./occ_0p75.txt", "w");
        if(f7==NULL) {
          printf("quit \n");
          exit(0);
        }
        f8=fopen("./rms_0.txt", "w");
        if(f8==NULL) {
          printf("quit \n");
          exit(0);
        }
        
        f15=fopen("./occ_0p8.txt", "w");
        if(f15==NULL) {
          printf("quit \n");
          exit(0);
        }
        f16=fopen("./occ_0p85.txt", "w");
        if(f16==NULL) {
          printf("quit \n");
          exit(0);
        }
        f17=fopen("./occ_0p9.txt", "w");
        if(f17==NULL) {
          printf("quit \n");
          exit(0);
        }
           amin=91;     
           Va=0;ks[0]=147;ks[1]=17; nuc=147-amin+1; a=0; ksm=amin+a; nuc=1;                
           for(i=0;i<nuc;i++) {
              soft[i][0]=amin+i; soft[i][1]=147-soft[i][0]; soft[i][2]=0;
           }
           soft[0][0]=147; ksm=147;          
       for(h=-250;h<250;h=h+2) {  
          E0=(double)h/(10*(147-91));
          for(w=0;w<500;w=w+2) {
             mu=-(double)(300-w)/5;CN=exp(mu);
             Osum=0; Osum1=0; Osum2=0; Osum3=0;
             for(i=0;i<L1;i++) {
                for(t=0;t<nuc;t++) {
                   Eij[i][t]=-E0*soft[t][0];//Eij[i][t]=gamaN*pot[i]+Va;
                }
                Zf[i]=0; Zb[i]=0; Pn[i]=0; O[i]=0; On[i]=0; On1[i]=0; On2[i]=0; On3[i]=0; //Es[i]=Eij[i][nuc-1];
             }
           
             for(i=0;i<(L1-ksm+1);i++) { // DNA length
                Zbsum=0; Zfsum=0; 
                for(t=0;t<nuc;t++) { // sub nucleosome
                   sze=soft[t][0]+a;
                   Zfsum=Zfsum+Zfor(i,ksm,sze,t,CN,Zf,Eij); // forward nuc                    
                   Zbsum=Zbsum+Zbor(i,ksm,sze,t,CN,Zb,Eij); // backward nuc                                                  
                }                                       
                Zf[i+ksm-1]=Zf[i+ksm-2] + log(1+Zfsum); // forward sums
                Zb[L1-ksm-i]=Zb[L1-ksm-i+1] + log(1+Zbsum); // backward sums
             }
                          
             sO=0; pO=0; lx=0; //printf("t %ld dx %ld gsum %lf w %ld \n",t,dx,gsum1,w);
             for(t=0;t<nuc;t++) {  // subnucleosome
                for(i=0;i<L1;i++) {
                   Pn[i]=0;O[i]=0;
                }
               sze=soft[t][0]+a; Pnsum=0; sumOmd=0;          
               Pn[0]=exp(Zb[0+sze]-Zf[L1-1])*CN*exp(-Eij[0][t]); Pnsum=Pnsum+Pn[0];          
               for(i=1;i<=(L1-sze-1);i++) {
                  Pn[i]=exp(Zf[i-1]+Zb[i+sze]-Zf[L1-1])*CN*exp(-Eij[i][t]); Pnsum=Pnsum+Pn[i]; 
               }
               Pn[L1-sze]=exp(Zf[L1-sze-1]-Zf[L1-1])*CN*exp(-Eij[L1-soft[t][0]][t]); Pnsum=Pnsum+Pn[L1-sze]; pdO[t]=Pnsum; pO=pO+pdO[t]; 
               /*if((pdO[t]-floor(pdO[t]))>0.5) {
                 lx=lx+(floor(pdO[t])+1)*soft[t][0];
               }
               else {
                 lx=lx+floor(pdO[t])*soft[t][0];
               }*/
               lx=lx+pdO[t]*soft[t][0];
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
               
               if(ks1>128 && ks1<148) {
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
                sumd=sumd+pow(Omd[t]-(sO/nuc),2.0); sump=sump+pow(pdO[t]-(pO/nuc),2.0);
             }
             d1=sqrt(sumd/nuc); d2=sqrt(sump/nuc);             

             
             if((Osum/L1)>0.699 && (Osum/L1)<0.701) {
                fprintf(f6,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.749 && (Osum/L1)<0.751) {
                fprintf(f7,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola2 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((d1>0 && d1<0.0001) && ((Osum/L1)>0.1)){
                fprintf(f8,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola4 E0d %lf mu %lf \n",E0*(147-91),mu);
             }

             fprintf(f4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); 
             //fprintf(f5,"%ld %ld %lf \n",h,w,den1b); fprintf(f6,"%ld %ld %lf \n",h,w,den2b); fprintf(f7,"%ld %ld %lf \n",h,w,den3b);  
             
             if((Osum/L1)>0.799 && (Osum/L1)<0.801) {
                fprintf(f15,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola3 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.849 && (Osum/L1)<0.851) {
                fprintf(f16,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola1 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
             
             if((Osum/L1)>0.899 && (Osum/L1)<0.901) {
                fprintf(f17,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1); //printf("ola2 E0d %lf mu %lf \n",E0*(147-91),mu);
             }
         } 
         printf("ZfL %lf nuc %ld ksm %ld lx %lf pO %lf lm %lf h %ld w %ld \n",Zf[L-1],nuc,ksm,lx,pO,lm,h,w);  
      }
        fclose(f4); fclose(f6); fclose(f7); fclose(f8); fclose(f15); fclose(f16); fclose(f17); return(0); 
}


double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]) { //i,ksm,tfsize,t,CN,(Zf or Zb),E1
       double Zx; long int jfn;
       jfn=num0+num1-num2;               // jfn=i+ksm-soft[t][0]; // position of TF, nuc 
       if(jfn>=0) {                      // forward nuc
         if(jfn==0) {
           Zx=exp(-num5[num0+num1-2])*num4*exp(-num6[jfn][num3]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
         else {
           Zx=exp(num5[jfn-1]-num5[num0+num1-2])*num4*exp(-num6[jfn][num3]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
       }  
       else {
         Zx=0;
       }
       return Zx;
}

double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]) { //i,ksm,tfsize,t,CN,(Zf or Zb),E1
       double Zx; long int jfn,jb,jbxn;
       jfn=num0+num1-num2;               // jfn=i+ksm-soft[t][0]; // position of TF, nuc        
       jb=L1-num1-num0;                      // backward nuc
       jbxn=jb+num2;
       if(jfn>=0) {                      
         if(jfn==0) {
           Zx=exp(-num5[L1-num1-num0+1])*num4*exp(-num6[jb][num3]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
         else {
           Zx=exp(num5[jbxn]-num5[L1-num1-num0+1])*num4*exp(-num6[jb][num3]); //Eij[jfn][t]=Es[jfn-soft[t][2]]+soft[t][1]*E0+log(CN);
         }
       }  
       else {
         Zx=0;
       }
       return Zx;
}

