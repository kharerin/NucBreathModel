#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
        FILE *f4;
        char c[10];
        long int i=0,amin,amax,h,h1;
        double px,mu,E0,E01;
        
        f4=fopen("./d1_mu_E0d.txt", "w");
        if(f4==NULL) {
          printf("quit \n");
          exit(0);
        }
        
        amin=91; amax=147;     
                            
        for(h=-250;h<0;h++) {  
           E0=(double)h/(10*(147-91));
           px=-E0/(1-((amax+amin)*E0)/2);
           mu=log(-E0/(amax-amin+1));
           fprintf(f4,"%lf %lf %lf \n",mu,E0,px);                       
        } 
        printf("h %ld done! \n",h);  
        fclose(f4); return(0); 
}





