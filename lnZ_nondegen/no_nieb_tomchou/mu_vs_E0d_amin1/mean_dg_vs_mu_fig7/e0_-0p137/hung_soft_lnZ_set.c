#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define Ns 147 // 147-91+1, min size 10
#define L2 L1*Ns
#define L3 545*Ns
#define n1 3
   
double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]);
double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]); 
int main() {
        FILE *f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9,*f15,*f16,*f17;
        char c[10];
        long int i=0,i1,i2,dx,j,k,posit,ks[2],l=147,L=L1,dw0,sze,sze1,soft[1000][3],nuc,ksm,ks1,t,jfn,jb,jbxn,h=1,w,b1,b2,js,a,amin,pk,tn,hx;
        double p[L1],p1,p0,pot[L1],Es[L1],Va,gamaN=1,CN,mu,E0,Zf[L1],Zb[L1],Pn[L1],O[L1],Zfsum,Zbsum,pOx,pdO[Ns],GN0,sumG,dg,Eh[6],mh[9];
        double Osum,Osum1,Osum2,Osum3,nucavg,nucavg1,nucavg2,nucavg3,On[L1],On1[L1],On2[L1],On3[L1],Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b;
        double Omd[Ns],sO,Pnsum,Pnsum1,sumOmd,pO,pO1,d1,d2,Zm,gsum,gsum1,sumd,sump,lx,lm;
        static double Eij[L1][Ns],gn[2000][Ns],Gn[2000],GN[2000],Gp[5][2];
        
        tn=1; hx=tn-1;
        
        f4=fopen("./d1_mu_E0d.txt", "w");
        if(f4==NULL) {
          printf("quit \n");
          exit(0);
        }

        amin=1;     
        Va=0;ks[0]=147;ks[1]=17; nuc=147-amin+1; a=20; ksm=amin+a;                   
        for(i=0;i<nuc;i++) {
            soft[i][0]=amin+i; soft[i][1]=147-soft[i][0]; soft[i][2]=0;
        }
        //mh[0]=-8.6; mh[1]=-8.25; mh[2]=-7.9; mh[3]=-7.55; mh[4]=-7.2; mh[5]=-6.85; mh[6]=-6.5; mh[7]=-6.15; mh[8]=-5.8;
        //Eh[0]=-10; Eh[1]=0; Eh[2]=5; Eh[3]=20; Eh[4]=40; Eh[5]=60;
                            
        for(h=hx;h<(hx+1);h=h+1) {  
           E0=0.137;
           for(w=0;w<20;w=w+1) {
              mu=-22+(double)w*1.1;CN=exp(mu);
              for(i=0;i<L1;i++) {
                 for(t=0;t<nuc;t++) {
                    Eij[i][t]=-E0*soft[t][0];//Eij[i][t]=gamaN*pot[i]+Va;
                 }
                 Zf[i]=0; Zb[i]=0; Pn[i]=0; O[i]=0; On[i]=0; On1[i]=0; On2[i]=0; On3[i]=0; //Es[i]=Eij[i][nuc-1];
              }
             
              for(i=0;i<2000;i++) {
                 Gn[i]=0;
                 for(j=0;j<nuc;j++) {
                    gn[i][j]=0;
                 }          
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
                       
              pO=0; pO1=0; lx=0; //printf("t %ld dx %ld gsum %lf w %ld \n",t,dx,gsum1,w);
              for(t=0;t<nuc;t++) {  // g_r subnucleosome 
                 sze=soft[t][0]+a;        
                 for(dx=0;dx<(1169-sze);dx++) {
                    if(dx==0) {
                      Zm=0;
                    }
                    else {
                      Zm=Zf[dx-1]; 
                    }
                    gsum1=0;
                    for(j=0;j<nuc;j++) {
                       sze1=soft[j][0]+a; gsum=0;
                       //gsum=gsum+exp(Zm+Zb[0+dx+sze+sze1]-Zf[L1-1])*CN*exp(-Eij[0][t])*CN*exp(-Eij[0+dx+sze][j]);            
                       for(i=4500;i<5000;i++) {  //for(i=0;i<=(L1-(dx+sze+sze1)-1);i++) {              
                          gsum=gsum+exp(Zf[i-1]+Zm+Zb[i+dx+sze+sze1]-Zf[L1-1])*CN*exp(-Eij[i][t])*CN*exp(-Eij[i+dx+sze][j]);
                       }
                       // gsum=gsum+exp(Zf[L1-(dx+sze+sze1)-1]+Zm-Zf[L1-1])*CN*exp(-Eij[L1-(dx+sze+sze1)][t])*CN*exp(-Eij[L1-(dx+sze+sze1)+dx+sze][j]);
                       gsum1=gsum1+gsum; //gsum1=gsum1+gsum/(2+L1-(dx+sze+sze1));
                    }                 
                    gn[sze+dx][t]=gsum1; //gn[soft[t][0]+dx][t]=gsum/(2+L1-(dx+2*soft[t][0]));
                    //printf("t %ld dx %ld gsum %lf %lf \n",t,dx,gsum1,Zf[L-1]);
                 }                                         
                 for(i=0;i<L1;i++) { // subnucleosome
                    Pn[i]=0;O[i]=0;
                 }
                 Pnsum=0;          
                 Pn[0]=exp(Zb[0+sze]-Zf[L1-1])*CN*exp(-Eij[0][t]); Pnsum=Pnsum+Pn[0];          
                 for(i=1;i<=(L1-sze-1);i++) {
                    Pn[i]=exp(Zf[i-1]+Zb[i+sze]-Zf[L1-1])*CN*exp(-Eij[i][t]); Pnsum=Pnsum+Pn[i]; 
                 }
                 Pn[L1-sze]=exp(Zf[L1-sze-1]-Zf[L1-1])*CN*exp(-Eij[L1-soft[t][0]][t]); Pnsum=Pnsum+Pn[L1-sze]; pdO[t]=Pnsum; pO=pO+pdO[t]; 
                 Pnsum1=0;
                 for(i=4500;i<5000;i++) {
                    Pnsum1=Pnsum1+exp(Zf[i-1]+Zb[i+sze]-Zf[L1-1])*CN*exp(-Eij[i][t]);
                 }
                 pO1=pO1+Pnsum1;  
                 lx=lx+pdO[t]*soft[t][0]; printf("i...%ld...t...%ld \n",i,t);                                       
              }
             
              for(i=0;i<5;i++) {
                 Gp[i][0]=0; Gp[i][1]=0;
              }
              pk=0; i1=0; i2=-1; GN0=0; //sumG=0;
              for(i=0;i<1169;i++) {
                 for(j=0;j<nuc;j++) {
                    Gn[i]=Gn[i]+gn[i][j];
                 }
                 GN[i]=(Gn[i]/pO1)/(pO/L1); //sumG=sumG+GN[i]*((double)i);
                 if(i>0) {
                   GN0=GN[i+i2];
                 }
                 if(pk<5 && i<1169) {
                   if(GN[i+i1]>=GN0) {
                     j=j+1; 
                   }
                   else {
                     j=0;  
                     if(i1>-1 && GN[i-1]>1) {
                       Gp[pk][0]=i-1; Gp[pk][1]=GN[i-1]; i1=-1; i2=0; pk=pk+1;
                     }
                     else {
                       i1=0; i2=-1;
                     }
                   }
                 }
              }
              if(pk>1) {
                dg=(Gp[pk-1][0]-Gp[0][0])/(double)(pk-1);
              }
              else {
                dg=1169;
              }
              lm=0;
              if(pO>0) {
                lm=lx/pO;
              }
              /*for(i=0;i<5;i++) {
                 printf("Gp %lf %lf %ld %lf h %ld w %ld \n",Gp[i][0],Gp[i][1],pk-1,dg,h,w);
              }*/

              fprintf(f4,"%lf %lf %lf %lf %lf \n",mu,E0,pO,dg,lm);  
              //fprintf(f5,"%ld %ld %lf \n",h,w,den1b); fprintf(f6,"%ld %ld %lf \n",h,w,den2b); fprintf(f7,"%ld %ld %lf \n",h,w,den3b);  
                                      
              printf("mu loop....pO %lf lm %lf dg %lf h %ld w %ld \n",pO,lm,dg,h,w); 
           } 
           printf("ZfL %lf nuc %ld ksm %ld pO %lf dg %lf lm %lf pk %ld h %ld w %ld \n",Zf[L-1],nuc,ksm,pO,dg,lm,pk,h,w);  
        }
        fclose(f4); return(0); 
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


