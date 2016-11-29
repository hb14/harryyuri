# harryyuri

//
//  main.c
//  RNG copy
//
//  Use ran2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Here we define ran2

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
    //Long period (> 2×10^18) random number generator of L’Ecuyer with Bays-Durham shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.
    
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0)                             //Initialize.
    {
        if (-(*idum) < 1) *idum=1;              //Be sure to prevent idum = 0.
        else *idum = -(*idum);
        
        idum2 = (*idum);
        
        for (j = NTAB+7; j>=0; j--)             //Load the shuffle table (after 8 warm-ups).
        {
            k = (*idum)/IQ1;
            *idum = IA1*(*idum-k*IQ1)-k*IR1;
            
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        
        iy = iv[0];
    }
    
    k = (*idum)/IQ1;                              //Start here when not initializing.
    *idum = IA1*(*idum-k*IQ1)-k*IR1;              //Compute idum=(IA1*idum) % IM1 without overflows by Schrage’s method.
    if (*idum < 0) *idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2)-k*IR2;              //Compute idum2=(IA2*idum) % IM2 likewise.
    if (idum2 < 0) idum2 += IM2;
    j = iy/NDIV;                                  //Will be in the range 0..NTAB-1.
    iy = iv[j]-idum2;                             //Here idum is shuffled, idum and idum2 are combined to generate output.
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM*iy) > RNMX)
        return RNMX;                              //Because users don’t expect endpoint values.
    
    else
        return temp;
}

//Here we generate a random number

float gasdev(long *idum)                          //Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) as the source of uniform deviates.
{
    float ran2(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    if(*idum < 0 ) iset = 0;                      //Reinitialize
    if (iset == 0)                                //We don’t have an extra deviate handy, so
    {
        do
        {
            v1=2.0*ran2(idum)-1.0;                //pick two uniform numbers in the square
            v2=2.0*ran2(idum)-1.0;                //extending from -1 to +1 in each direction,
            rsq=v1*v1+v2*v2;                      //see if they are in the unit circle,
        }
        while (rsq >= 1.0 || rsq == 0.0);         //and if they are not, try again.
        fac=sqrt(-2.0*log(rsq)/rsq);
        //Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
        gset=v1*fac;
        iset=1;                                   //Set flag.
        
        return v2*fac;
    }
    else                                          //We have an extra deviate handy,
    {
        iset=0;                                   //so unset the flag,
        
        return gset;                              //and return it.
    }
}

int main()
{
    int i;
    long seed = (-1);
    
    double n; //Random number & potentialy a vector
    double KB=1.38064852*pow(10,-23), pi=3.141592653589793; //Boltzman, pi
    double T=300, t, eta=8.9*pow(10,-4), a=pow(10,-8), D, Z; //Temperature, time, time-step, viscosity, radius, friction
    double r2;
    double r[10000]; //Array of r^2 values
    
    FILE *fh;
    fh = fopen ("tryagain.csv","w"); //Creates a file called RNGdistribution.csv
    if (fh == NULL)
    {
        puts("Can't open that file!!\n"); //File error
        exit(1);
    }
    
    Z = (6*pi*eta*a); //Calculating the friction coefficient
    D = KB*T/Z;
    t = i;
    
    printf("D = %.15lf\n", D);
    
    n = gasdev(&seed); //Defining n as the random number
    
    r[0] = 2*sqrt((6*KB*T*i)/(Z))*n;
    
    for(i=1; i<10000; i++)
    {
        n = gasdev(&seed);
        
        r[i+1] = 2*sqrt((6*D*i))*n + r[i]; //Finding each r
    }
    
    for(i=0; i<10000; i++)
    {
        r2 = r[i]*r[i];
        fprintf(fh, "%d, %.15lf\n", i, r2);
    }
    
    fclose (fh);
    
    return 0;
}
