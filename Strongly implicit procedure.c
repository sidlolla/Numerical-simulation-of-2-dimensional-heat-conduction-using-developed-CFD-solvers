#include<stdio.h>
#include<math.h>
#define pi 3.141592654
#define nx 51
#define ny 51
#define len_x 2
#define len_y 2
#define Tm 2
main()
{
	int N=nx*ny,i,j,k=0;
	double dx,dy,beta,e=1;
	double x,y,alpha=0.5;
	double AP[N+1],AN[N+1],AE[N+1],AS[N+1],AW[N+1],b[N+1],b_new[N+1];
	double LP[N+1],LS[N+1],LW[N+1];
	double UP[N+1],UN[N+1],UE[N+1];
	double MSE[N+1],MNW[N+1];
	double T[N+1],T_new[N+1],z[N+1];
	double xpos,ypos,s1,s2,s3,s4,T_ana[nx+1][ny+1];
	
	dx = (double)len_x/(nx-1);
	dy = (double)len_y/(ny-1);
	beta = (double)dx/dy;
	
	//initial
	for(i=0;i<=N;i++){
		AP[i] = 0; AN[i] = 0; AE[i] = 0; AS[i] = 0; AW[i] = 0;b[i] = 0;
		LP[i] = 0; LS[i] = 0; LW[i] = 0;
		UP[i] = 1; UN[i] = 0; UE[i] = 0;
		MSE[i] = 0; MNW[i] = 0;
		T[i] = 0; T_new[i] = 0; b[i] = 0;
	}
	
	//A coefficients at internal nodes
	for(i=1;i<=N;i++){
		AP[i] = -2*(1+beta*beta);
	}
	for(i=2;i<=N;i++){
		AS[i] = (beta*beta);
	}	
	for(i=1;i<=N-1;i++){
		AN[i] = (beta*beta);
	}
	for(i=1+ny;i<=N;i++){
		AW[i] = 1;
	}
	for(i=1;i<=N-ny;i++){
		AE[i] = 1;
	}
	//at boundaries
	for(i=1;i<=ny;i++){
		AP[i]=1 ; AS[i]=0; AN[i]=0; AE[i]=0; AW[i]=0; b[i]=0; 	//left
	}
	for(i=N;i>=N-ny;i--){
		AP[i]=1 ; AS[i]=0; AN[i]=0; AE[i]=0; AW[i]=0; b[i]=0;	//right
	}
	for(i=1;i<=N;i=i+ny){
		AP[i]=1 ; AS[i]=0; AN[i]=0; AE[i]=0; AW[i]=0; b[i]=0;	//bottom
	}
	for(i=ny;i<=N;i=i+ny){
		x=((i/ny)-1)*dx;
		AP[i]=1 ; AS[i]=0; AN[i]=0; AE[i]=0; AW[i]=0; b[i]=Tm*sin((double)pi*x/len_x);	//top  Tm*sin(pi*x/L)
	}
	
	
	//creating L and U matrices
	for(i=1;i<=N;i++){
		if(i<ny){
			LW[i] = AW[i];
		}
		else{
			LW[i] = AW[i]/(1+alpha*UN[i-ny]);
		}
		
		LS[i] = AS[i]/(1+alpha*UE[i-1]);
		
		if(i<ny){
			LP[i] = AP[i] - LS[i]*UN[i-1] + alpha*(LS[i]*UE[i-1]);
		}
		else{
			LP[i] = AP[i] - LW[i]*UE[i-ny] - LS[i]*UN[i-1] + alpha*(LW[i]*UN[i-ny] + LS[i]*UE[i-1]);
		}
		
		if(i<ny){
			UN[i] = (double)(AN[i])/LP[i];
		}
		else{
			UN[i] = (double)(AN[i] - alpha*LW[i]*UN[i-ny])/LP[i];
		}
	
		UE[i] = (double)(AE[i] - alpha*LS[i]*UE[i-1])/LP[i];
	//	printf("%lf	%lf	%lf	%lf	%lf\n",LW[i],LS[i],LP[i],UN[i],UE[i]);		
	}
	for(i=ny;i<=N;i++){
		MNW[i] = LW[i]*UN[i-ny];
	}
	for(i=1;i<=N;i++){
		MSE[i] = LS[i]*UE[i-1];
	//	printf("%lf	%lf	\n",MNW[i],MSE[i]);	
	}
	
	while(e>0.000001){
	e=0;
		
		for(i=1;i<=N;i++){
			if(i<=ny-1){
				b_new[i] = b[i] - MNW[i]*(alpha*(T[i+1]-T[i])) + MSE[i]*(T[i+ny-1] -alpha*(T[i-1]+T[i+ny]-T[i]));
			}
			else if(i>=N-ny+1 && i<N){
				b_new[i] = b[i] + MNW[i]*(T[i-ny+1] -alpha*(T[i+1]+T[i-ny]-T[i])) - MSE[i]*(alpha*(T[i-1]-T[i]));
			}			
			else{
				b_new[i] = b[i] + MNW[i]*(T[i-ny+1] -alpha*(T[i+1]+T[i-ny]-T[i])) + MSE[i]*(T[i+ny-1] -alpha*(T[i-1]+T[i+ny]-T[i]));
			}			
		}
		
		for(i=1;i<=N;i++){
			if(i<ny){
				z[i] = (b_new[i]-LS[i]*z[i-1])/LP[i];
			}
			else{
				z[i] = (b_new[i]-LS[i]*z[i-1]-LW[i]*z[i-ny])/LP[i];
			}		
		}
		
		for(i=N;i>=1;i--){
			if(i==N){
				T_new[i] = z[i];
			}
			else if(i>N-ny&&i<N){
				T_new[i] = z[i]-UN[i]*T_new[i+1];
			}
			else{
				T_new[i] = z[i]-UN[i]*T_new[i+1]-UE[i]*T_new[i+ny];
			}
			e+=fabs(T_new[i]-T[i]);
			T[i] = T_new[i];
		}
		
	}
		
		
		//ANALYTIC SOLUTION
	for(i=1; i<=nx; i++){
    	for(j=1; j<=ny; j++){
    		T_ana[i][j]=0;
    		xpos=(i-1)*dx;
    		ypos=(j-1)*dy;
				s2=(double)pi*ypos/len_x;
				s3=(double)pi*len_y/len_x;
				s1 = (double)sinh(s2)/(double)sinh(s3);  
				s4 = (double)sin(pi*xpos/len_x); 					
    			T_ana[i][j] = Tm*s1*s4;    		
			    		
		}
	}

	
			
	FILE *fa,*fb,*fc,*fd,*fe,*ff;
	
	fa = fopen("part(2a).txt", "w+"); 
		fprintf(fa,"Temperature distribution\n\n");
		for(i=1;i<=ny;i++){
			for(j=i;j<=N;j=j+ny){
				fprintf(fa,"%lf	",T[j]);
			}
			fprintf(fa,"\n");
		}
	fd = fopen("part(2a) analytic.txt", "w+"); 
		fprintf(fd,"Analytic Temperature distribution\n\n");
		for(j=1; j<=ny; j++){
			for(i=1; i<=nx; i++){    		
				fprintf(fd,"%lf	",T_ana[i][j]);
			}
			fprintf(fd,"\n");
		}	
		
	fb = fopen("part(2b).txt", "w+");  
		fprintf(fb,"mid vertical line , x=0.5\n\n");
		i=(int)(nx)/2;
		i=i*ny+1;
			for(j=i;j<i+ny;j++){
				fprintf(fb,"%lf\n",T[j]);
		}
	
	fe = fopen("part(2b) analytic.txt", "w+");  
		fprintf(fe,"analytic mid vertical line , x=0.5\n\n");
		for(j=1; j<=ny; j++){ 
			i= (int)(nx-1)/2;  		
				fprintf(fe,"%lf\n",T_ana[i][j]);
			}

	
	fc = fopen("part(2c).txt", "w+");  
		fprintf(fc,"mid horizontal line , y=0.5\n\n");
		i=(int)(ny+1)/2;
			for(j=i;j<N;j=j+ny){
				fprintf(fc,"%lf\n",T[j]);
		}
	
	ff = fopen("part(2c) analytic.txt", "w+");  
		fprintf(ff,"analytic mid horizontal line , y=0.5\n\n");
		for(j=1; j<=ny; j++){ 
			i= (int)(ny+1)/2;  		
				fprintf(ff,"%lf\n",T_ana[j][i]);
			}
		
	
	printf("\n\nRESULTS HAVE BEEN PRINTED IN OUTPUT TEXT FILES\n\n");
}
