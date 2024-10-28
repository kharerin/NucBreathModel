#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define L2 L1*57
#define L3 545*57
#define n1 3	
float ran2(long *idum);    
double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][57]);
double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][57]); 
int main() {
        FILE *f2,*f10,*f11;
        long int i=0,i1,i2,i3,j,k,nTF=2,posit,ks[2],l=147,L=L1,dw0,sze,soft[1000][3],nuc,ksm,ks1,t,w[n1],b1[n1],b2[n1],b3,hc1,hc2,cen[n1],h[n1],Nsim,b1in,b1out,b2in,b2out;
        double p[L1],p1,p0,pot[L1],Es[L1],Va,gamaN=1,CN,mu,E0=0.05,Zf[L1],Zb[L1],Pn[L1],O[L1],Zfsum,Zbsum,Esum,cy,cy1,hmx1,hyx,hmn1,hmx,hmn,hy,sd,px[57],Px[57],hy1[2];
        double Dy1,Dy2,Dy3,Dy4,sumdy,score,subf,subf1,Dyt,sum1,sum2,sum3,d,dkl,dyt[L3],pp3[L3],py2[L3];;
        double Osum=0,Osum1=0,Osum2=0,Osum3=0,nucavg,nucavg1,nucavg2,nucavg3,On1,On2,On3,On4,Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b,frag[57],frag1[57],frag4[57];
        static double Eij[L1][57],subheat[L1][57],den5[L1][57],den3[L1][57],dya[L1][57],Subheat[2000][57],Den5[2000][57],Den3[2000][57],Dya[2000][57],Eseq[L1][57],niebE[L1][2],NiebE[2000][2];
                      
                     
        Va=0;ks[0]=147;ks[1]=17; nuc=56+1; ksm=91;                   
        for(i=0;i<nuc;i++) {
           soft[i][0]=91+i; soft[i][1]=147-soft[i][0]; soft[i][2]=0; frag[i]=0; frag1[i]=0; frag4[i]=0;
        }             
        for(i=0;i<2000;i++) {
           NiebE[i][0]=0; NiebE[i][1]=0;
           for(j=0;j<nuc;j++) {
              Subheat[i][j]=0; Den5[i][j]=0; Den3[i][j]=0; Dya[i][j]=0;
           }          
        }
       
        E0=5.0/(147-91); mu=-25.0; CN=exp(mu); Nsim=1;     
        for(i1=0;i1<Nsim;i1++) {                                                      
           for(i=0;i<L1;i++) {
              for(t=0;t<nuc;t++) { 
                 Eij[i][t]=0; 
              }
           }
           Osum=0; Osum1=0; Osum2=0; Osum3=0; h[0]=10; j=0; 
           w[0]=50; w[1]=200; w[2]=50; cy=5; hmn=5; 

           for(i=0;i<L1;i++) { 
              for(t=0;t<nuc;t++) {
                 Eseq[i][t]=-E0;                                                                                                                   
              }
           }                                   
           for(i=0;i<L1;i++) {
              for(t=0;t<nuc;t++) {
                 if(i<=(L1-soft[t][0])) {
                   Esum=0;
                   for(j=0;j<soft[t][0];j++) {
                      Esum=Esum+Eseq[i+j][t];                      
                   }
                   Eij[i][t]=Esum;
                 }
              }
              Es[i]=Eij[i][nuc-1]; Zf[i]=0; Zb[i]=0; Pn[i]=0; O[i]=0;
           }

           for(i=0;i<(L1-ksm+1);i++) { // DNA length
              Zbsum=0; Zfsum=0; 
              for(t=0;t<nuc;t++) { // sub nucleosome
                 sze=soft[t][0];
                 Zfsum=Zfsum+Zfor(i,ksm,sze,t,CN,Zf,Eij); // forward nuc                    
                 Zbsum=Zbsum+Zbor(i,ksm,sze,t,CN,Zb,Eij); // backward nuc                                                  
              }                                       
              Zf[i+ksm-1]=Zf[i+ksm-2] + log(1+Zfsum); // forward sums
              Zb[L1-ksm-i]=Zb[L1-ksm-i+1] + log(1+Zbsum); // backward sums
           }        
           for(t=0;t<nuc;t++) {  // subnucleosome
              for(i=0;i<L1;i++) {
                 Pn[i]=0;O[i]=0; den5[i][t]=0; den3[i][t]=0; dya[i][t]=0;
              }           
              Pn[0]=exp(Zb[0+soft[t][0]]-Zf[L1-1])*CN*exp(-Eij[0][t]); 
              den5[0][t]=Pn[0]; den3[0+soft[t][0]-1][t]=Pn[0]; dya[0+(long int)floor(soft[t][0]/2)-1][t]=Pn[0];
              Pn[L1-soft[t][0]]=exp(Zf[L1-soft[t][0]-1]-Zf[L1-1])*CN*exp(-Eij[L1-soft[t][0]][t]); 
              den5[L1-soft[t][0]][t]=Pn[L1-soft[t][0]]; den3[L1-1][t]=Pn[L1-soft[t][0]]; dya[L1-soft[t][0]+(long int)floor(soft[t][0]/2)-1][t]=Pn[L1-soft[t][0]];             
              for(i=1;i<=(L1-soft[t][0]-1);i++) {
                 Pn[i]=exp(Zf[i-1]+Zb[i+soft[t][0]]-Zf[L1-1])*CN*exp(-Eij[i][t]); 
                 den5[i][t]=Pn[i]; den3[i+soft[t][0]-1][t]=Pn[i]; dya[i+(long int)floor(soft[t][0]/2)-1][t]=Pn[i];
              }               
              ks1=soft[t][0];
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
                   subheat[i][t]=O[i]; 
                }                                                    
            }

            b3=5000+100;
            for(i=0;i<2000;i++) {
               NiebE[i][0]=NiebE[i][0]+niebE[b3-1000+i][0]; NiebE[i][1]=NiebE[i][1]+niebE[b3-1000+i][1];
               for(j=0;j<nuc;j++) {
                  Subheat[i][j]=Subheat[i][j]+subheat[b3-1000+i][j]; Den5[i][j]=Den5[i][j]+den5[b3-1000+i][j]; Den3[i][j]=Den3[i][j]+den3[b3-1000+i][j]; Dya[i][j]=Dya[i][j]+dya[b3-1000+i][j];
               }
            }
            
            Osum1=0;
            for(j=0;j<nuc;j++) {
               Osum=0;
               for(i=0;i<1000;i++) {
                  Osum=Osum+dya[b3+i][j];
               }
               frag[j]=Osum;
               Osum1=Osum1+Osum;
            }
            for(j=0;j<nuc;j++) {
               frag[j]=frag[j]/Osum1;
            }
                       
            Osum1=0;
            for(j=0;j<nuc;j++) {
               Osum=0;
               for(i=0;i<100;i++) {
                  Osum=Osum+dya[b3+i][j];
               }
               frag1[j]=Osum;
               Osum1=Osum1+Osum;
            }
            for(j=0;j<nuc;j++) {
               frag1[j]=frag1[j]/Osum1;
            }
            
            Osum1=0;
            for(j=0;j<nuc;j++) {
               Osum=0;
               for(i=0;i<100;i++) {
                  Osum=Osum+dya[b3+400+i][j];
               }
               frag4[j]=Osum;
               Osum1=Osum1+Osum;
            }
            for(j=0;j<nuc;j++) {
               frag4[j]=frag4[j]/Osum1;
            }
            
            printf("i1 %ld ZfL %lf nuc %ld ksm %ld b3 %ld \n",i1,Zf[L-1],nuc,ksm,b3);   
        }
     
        f2=fopen("./frag_nie2_occ_0pxx.txt", "w"); //dendya_h10_flat2_nie3.txt
        if(f2==NULL) {
          printf("quit \n");
          exit(0);
        }
        
        for(i=0;i<57;i++) {
           fprintf(f2,"%lf %lf %lf \n",frag1[i],frag4[i],frag[i]);
        }                
        
        f10=fopen("./dyamapfull_h10_sigm8_nie2.txt", "w");
        if(f10==NULL) {
           printf("quit \n");
           exit(0);
        }
        
        f11=fopen("./Emapfull_h10_sigm8_nie2.txt", "w");
        if(f11==NULL) {
           printf("quit \n");
           exit(0);
        }        
        
        for(i=0;i<L1;i++) {
           for(j=0;j<nuc;j++) {
              fprintf(f10,"%ld %ld %lf \n",i,91+j,dya[i][j]);
              fprintf(f11,"%ld %ld %lf \n",i,91+j,Eij[i][j]);
           }
        }
              
        printf("E0d %lf mu %lf frag[0] %lf frag[56] %lf\n",E0*56,mu,frag[0],frag[56]); 
        fclose(f2); fclose(f10); fclose(f11); return(0); 
}




double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][57]) { //i,ksm,tfsize,t,CN,(Zf or Zb),E1
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

double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][57]) { //i,ksm,tfsize,t,CN,(Zf or Zb),E1
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


