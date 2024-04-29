#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <random>

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
int main(int argc, char** argv){
    double tend,dt;
    int n;
    std::cout<<argc<<std::endl;
    if (argc==4) 
    {
        n=std::stoi(argv[1]);
        tend=std::stod(argv[2]);
        dt=std::stod(argv[3]);
        std::cout<<"success running"<<std::endl;
    }
    else 
    {
        n=200,tend=10,dt=0.1;
    }
        std::cout<<"n=  "<<n<<"   t=  "<<tend<<"   dt=  "<<dt<<std::endl;
    std::vector<double> u;
    u.resize(n,-2.0);
    double h=2.0/n,eps=1e-8;
    for(double t=0;t<tend+dt/2;t=t+dt)
    {
        std::string fname="BurgersT"+std::to_string(t)+".dat";
        std::fstream file(fname,std::ios::out);
        file<<"Title=\"exact solution of burgers equation at t = "<<t<<"\""<<'\n';
        file<< "Variables=\"x\",\"theta\",\"u\" \n";

        for(int i =0;i<n;i++)
        {
            double theta=u[i];
            if (i>0) theta=u[i-1];
            double x=h/2.0+i*h-1.0;
            double delta=100;
            int tt=0;
            do
            {
                delta=-f(theta,x,t)/df(theta,x,t);
                double thetaPre=theta;
                theta=theta+0.5*delta;
                tt++;
                if (tt==50) 
                { 
                    theta=u[i];
                    continue;
                }
                //std::cout<<times<<' '<<theta<<'\n';
            } while (abs(delta)>eps);
            u[i]=theta;
            file<<x<<' '<<theta<<' '<<u0(theta)<<'\n';
        }
        file.close();
        std::cout<<"t= "<<t<<"  finished"<<'\n';
    }

    


    
    
}
