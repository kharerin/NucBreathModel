#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L1 10000        //length of DNA 
#define Ns 57 // 147-10+1, min size 10
#define L2 L1*Ns
#define L3 545*Ns
#define n1 3	
 
double Zfor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]);
double Zbor(long int num0,long int num1,long int num2,long int num3,double num4,double num5[L1],double num6[L1][Ns]); 
int main() {
        FILE *f2,*f10,*f11,*f12;
        long int i,i1,i2,i3,j,k,dx,ks[2],l=147,L=L1,dw0,sze,sze1,soft[1000][3],nuc,ksm,ks1,t,w[n1],b1[n1],b2[n1],b3,hc1,hc2,cen[n1],h[n1],Nsim,b1in,b1out,b2in,b2out,a,amin;
        double p[L1],p1,p0,Es[L1],Va,CN,mu,E0,Zf[L1],Zb[L1],Pn[L1],O[L1],Zfsum,Zbsum,Esum,cy,cy1,hmx1,hyx,hmn1,hmx,hmn,hy,sd,px[Ns],Px[Ns],hy1[2];
        double Dy1,Dy2,Dy3,Dy4,sumdy,subf,subf1,Dyt,d,dkl,dyt[L3],pp3[L3],py2[L3],Zm,gsum,Pnsum,pO,gsum1,pO1,Pnsum1;
        double Osum=0,Osum1=0,Osum2=0,Osum3=0,nucavg,nucavg1,nucavg2,nucavg3,On1,On2,On3,On4,Osumb,Osum1b,Osum2b,Osum3b,denb,den1b,den2b,den3b,frag[Ns],frag1[Ns],frag4[Ns];
        static double Eij[L1][Ns],subheat[L1][Ns],den5[L1][Ns],den3[L1][Ns],dya[L1][Ns],Subheat[2000][Ns],Den5[2000][Ns],Den3[2000][Ns],Dya[2000][Ns],Eseq[L1][Ns],niebE[L1][2],NiebE[2000][2],gn[2000][Ns],Gn[2000];
                                
        i=0; amin=91;               
        Va=0;ks[0]=147;ks[1]=17; nuc=147-amin+1; a=0; ksm=amin+a;                  
        for(i=0;i<nuc;i++) {
           soft[i][0]=amin+i; soft[i][1]=147-soft[i][0]; soft[i][2]=0; frag[i]=0; frag1[i]=0; frag4[i]=0;
        }            
        for(i=0;i<2000;i++) {
           NiebE[i][0]=0; NiebE[i][1]=0; Gn[i]=0;
           for(j=0;j<nuc;j++) {
              Subheat[i][j]=0; Den5[i][j]=0; Den3[i][j]=0; Dya[i][j]=0; gn[i][j]=0;
           }          
        }
             
        E0=0.0/(147-amin); mu=5.0; CN=exp(mu); Nsim=1;  //E0=10.0/(147-91);  
        for(i1=0;i1<Nsim;i1++) {
           hmx=10; hmx1=10;                                                        
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
                 sze=soft[t][0]+a;
                 Zfsum=Zfsum+Zfor(i,ksm,sze,t,CN,Zf,Eij); // forward nuc                    
                 Zbsum=Zbsum+Zbor(i,ksm,sze,t,CN,Zb,Eij); // backward nuc                                                  
              }                                       
              Zf[i+ksm-1]=Zf[i+ksm-2] + log(1+Zfsum); // forward sums
              Zb[L1-ksm-i]=Zb[L1-ksm-i+1] + log(1+Zbsum); // backward sums
           }   
           
           for(t=0;t<nuc;t++) {  // g_r subnucleosome }  
              sze=soft[t][0]+a;        
              for(dx=0;dx<(2000-sze);dx++) {
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
                    for(i=2500;i<7500;i++) {  //for(i=0;i<=(L1-(dx+sze+sze1)-1);i++) {              
                       gsum=gsum+exp(Zf[i-1]+Zm+Zb[i+dx+sze+sze1]-Zf[L1-1])*CN*exp(-Eij[i][t])*CN*exp(-Eij[i+dx+sze][j]);
                    }
                   // gsum=gsum+exp(Zf[L1-(dx+sze+sze1)-1]+Zm-Zf[L1-1])*CN*exp(-Eij[L1-(dx+sze+sze1)][t])*CN*exp(-Eij[L1-(dx+sze+sze1)+dx+sze][j]);
                    gsum1=gsum1+gsum; //gsum1=gsum1+gsum/(2+L1-(dx+sze+sze1));
                 }                 
                 gn[sze+dx][t]=gsum1; //gn[soft[t][0]+dx][t]=gsum/(2+L1-(dx+2*soft[t][0]));
                 printf("t %ld dx %ld gsum %lf \n",t,dx,gsum1);
              }
           } 
           pO=0; pO1=0;  
           for(t=0;t<nuc;t++) {  // subnucleosome
              for(i=0;i<L1;i++) {
                 Pn[i]=0;O[i]=0; den5[i][t]=0; den3[i][t]=0; dya[i][t]=0;
              }   
              sze=soft[t][0]+a; Pnsum=0;        
              Pn[0]=exp(Zb[0+sze]-Zf[L1-1])*CN*exp(-Eij[0][t]); Pnsum=Pnsum+Pn[0];
              den5[0][t]=Pn[0]; den3[0+soft[t][0]-1][t]=Pn[0]; dya[0+(long int)floor(soft[t][0]/2)-1][t]=Pn[0]; 
              Pn[L1-sze]=exp(Zf[L1-sze-1]-Zf[L1-1])*CN*exp(-Eij[L1-soft[t][0]][t]); Pnsum=Pnsum+Pn[L1-sze];
              den5[L1-sze][t]=Pn[L1-sze]; den3[L1-1-a][t]=Pn[L1-sze]; dya[L1-sze+(long int)floor(soft[t][0]/2)-1][t]=Pn[L1-sze];             
              for(i=1;i<=(L1-sze-1);i++) {
                 Pn[i]=exp(Zf[i-1]+Zb[i+sze]-Zf[L1-1])*CN*exp(-Eij[i][t]); Pnsum=Pnsum+Pn[i]; 
                 den5[i][t]=Pn[i]; den3[i+soft[t][0]-1][t]=Pn[i]; dya[i+(long int)floor(soft[t][0]/2)-1][t]=Pn[i];
              }
              Pnsum1=0;
              for(i=2500;i<7500;i++) {
                 Pnsum1=Pnsum1+exp(Zf[i-1]+Zb[i+sze]-Zf[L1-1])*CN*exp(-Eij[i][t]);
              }                
              ks1=soft[t][0]; pO=pO+Pnsum; pO1=pO1+Pnsum1;               
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
                  Subheat[i][j]=Subheat[i][j]+subheat[b3-1000+i][j]; Den5[i][j]=Den5[i][j]+den5[b3-1000+i][j]; Den3[i][j]=Den3[i][j]+den3[b3-1000+i][j]; Dya[i][j]=Dya[i][j]+dya[b3-1000+i][j]; Gn[i]=Gn[i]+gn[i][j];
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
        
        f12=fopen("./gr.txt", "w");
        if(f12==NULL) {
           printf("quit \n");
           exit(0);
        }
        
        for(i=0;i<2000;i++) {
           fprintf(f12,"%ld %lf \n",i,(Gn[i]/pO1)/(pO/L1));
        }
        
        for(i=0;i<L1;i++) {
           for(j=0;j<nuc;j++) {
              fprintf(f10,"%ld %ld %lf \n",i,amin+j,dya[i][j]);
              fprintf(f11,"%ld %ld %lf \n",i,amin+j,Eij[i][j]);
           }
        }
               
        printf("E0d %lf mu %lf Gn[1999] %lf density %lf done! \n",E0,mu,Gn[1999],pO/L1);  
        fclose(f2); fclose(f10); fclose(f11); fclose(f12); return(0); 
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




