#include"burgers.hpp"
#include<math.h>
double u0(double x)
{
    return 0.25+0.5*sin(M_PI*x);
}
double du0(double x)
{
    return 0.5*M_PI*cos(M_PI*x);
}
double f(double theta,double x,double t)
{
    return theta+u0(theta)*t-x;
}
double df(double theta,double x,double t)
{
    return 1.0+du0(theta)*t;
}
double exactBurgers(double x,double t,double theta0)
{
    
    double theta=theta0;
    double delta=100;
    double eps=1e-10;
    int tt=0;
    do
    {
        delta=-f(theta,x,t)/df(theta,x,t);
        double thetaPre=theta;
        theta=theta+0.5*delta;
        tt++;
        //std::cout<<times<<' '<<theta<<'\n';
    } while (abs(delta)>eps);
    return u0(theta);
}