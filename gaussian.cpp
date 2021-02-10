#include <math.h>
//#include "random.cpp"
extern double genrand_res53();
extern long unsigned int genrand_int32();
extern void init_genrand(unsigned long);
// interface file

class CDataDistrib
{
public: 
CDataDistrib(){m_dblVar = 0;};

double GetNormVar(double flMean, double flStDev, bool bAllowSubZero);

private:
double GetNormVar();
double m_dblVar;
};

// implentation file...
double CDataDistrib::GetNormVar(double flMean, double flStDev, bool bAllowSubZero)
{
// call private func.
m_dblVar = GetNormVar(); // returns a value - btwn -1 and +1 ...if it's -ve - then the value is under mean
//return 
m_dblVar = (m_dblVar*flStDev)+flMean; // make it within suitable range

// if below zero recursively call this func
if(!bAllowSubZero && m_dblVar < 0)
m_dblVar = GetNormVar(flMean,flStDev,bAllowSubZero);

return m_dblVar; // default return
}

double CDataDistrib::GetNormVar()
{
static int iset=0;
static double gset;
double fac,rsq,v1,v2;

iset =0;//defo called first time round

if (iset ==0)
{ 
do 
{
v1= 2.0*genrand_res53()-1.0; 
v2= 2.0*genrand_res53()-1.0;
rsq=v1*v1+v2*v2;
} 
while (rsq >= 1.0 || rsq ==0.0); // repeats the call of the above code ubtil the value is &lt; 1

fac=sqrt(-2.0*log(rsq)/rsq); 
gset=v1*fac; 
iset=1; // iset = 1

return v2*fac; // 0.583208 = 0.784079 * 0.743813
} 
else 
{
iset=0;
return gset;
}
}
