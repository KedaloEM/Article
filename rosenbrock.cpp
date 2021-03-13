#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <stdio.h>
#include <typeinfo>
#include <stdlib.h>
#include <vector>
#include "units.hpp"
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;
namespace phoenix = boost::phoenix;

typedef vector<double> state_type;
typedef matrix<double> matrix_type;

double kvv (int i, int i_f, int j, int j_f)
{
    int max1 = std::max(i,i_f);
    int max2 = std::max(j,j_f);
    return (kVV0*max1*max2*Q01_VV*exp(-deltaVT*abs(max1-max2))*(2.0/3.0-0.5*exp(-deltaVT*abs(max1-max2))));
}


double kvt (int j)
{
    return (kVT0*(j+1)*P01_VT*exp(deltaVT*j));
}


double kev (int i, int j)
{
    if (i==0&&j==1)
    {
        return kev0*exp(-0.7*fabs(j-i))/(1+0.05*i);
    }
    else
    {
        return 0;
    }

}

double Macheret(int j)
{
    return exp(MF*h*c*w*j*(1-xe*j)/(kb*T));
}


struct stiff_system
{
void operator() ( const state_type &x_inp, state_type &dxdt, const double /* t */ )
{
    state_type x = x_inp;
    x *= Ng;
    state_type subrange1;
    state_type subrange2;
    subrange1 = subrange(x,nvibr,2*nvibr);
    subrange2 = subrange(x,2*nvibr,3*nvibr);
    double N1_tot = sum(subrange1) + x(N-10)+x(N-6)+x(N-4)+x(N-2)+x(N-7)*surf_rate_1;//fast H diffusion is taken for account
    double N2_tot = sum(subrange2) + x(N-9)+x(N-5)+x(N-3)+x(N-1)+x(N-7)*surf_rate_2;

    for (int v = 0; v < nvibr; v++) //gas equations
    {
        dxdt(v) = surf_rate_1*(-x(v)*st*wD*(1.0 - N1_tot/a_max) + width_r*kdes1*x(v+nvibr)) +
                  surf_rate_2*(-x(v)*st*wD*(1.0 - N2_tot/a_max) + width_r*kdes2*x(v+2*nvibr));
        if (v != nvibr-1)
        {
            dxdt(v) += kvt(v+1)*x(v+1)*Ng;
        }
        if (v != 0)
        {
            dxdt(v) += -kvt(v)*x(v)*Ng;
        }
        for (int u = 0; u < nvibr-1; u++)
        {
            if (v!=0 && v!=nvibr-1)
            {
                dxdt(v) += -kvv(v,v+1,u+1,u)*x(v)*x(u+1) - kvv(v,v-1,u,u+1)*x(v)*x(u) + kvv(v-1,v,u+1,u)*x(v-1)*x(u+1) + kvv(v+1,v,u,u+1)*x(v+1)*x(u);
                if (u > v)
                {
                    dxdt(v) += -kev(v,u)*ne*x(v) + kev(v,u)*ne*x(u);
                }
                else if (u < v)
                {
                    dxdt(v) += kev(u,v)*ne*x(u) - kev(u,v)*ne*x(v); //???
                }
            }
            else if (v == 0)
            {
                dxdt(v) += -kvv(v,v+1,u+1,u)*x(v)*x(u+1)  + kvv(v+1,v,u,u+1)*x(v+1)*x(u) - kev(v,u)*ne*x(v) + kev(v,u)*ne*x(u);
            }
            else if (v == nvibr-1)
            {
                dxdt(v) +=  -kvv(v,v-1,u,u+1)*x(v)*x(u) + kvv(v-1,v,u+1,u)*x(v-1)*x(u+1) + kev(u,v)*ne*x(u) - kev(u,v)*ne*x(v);
            }
        }

    }

   // dxdt(N-2) = 0; //surface atomic N equation
   // dxdt(N-1) = 0;

    for (int v = nvibr; v < 2*nvibr; v++) //surface1 equations
    {
        dxdt(v) = (wD*x(v-nvibr)*st/width_r)*(1-N1_tot/a_max) - kdes1*x(v) - ka_1*Macheret(v-nvibr)*x(v)*(1-N1_tot/a_max);
        if (v!= 2*nvibr-1)
        {
            dxdt(v) += kvt_s*x(v+1);
        }
        if (v!= nvibr)
        {
            dxdt(v) += -kvt_s*x(v);
        }
        dxdt(N-10) += ka_1*Macheret(v-nvibr)*x(v)*(1-N1_tot/a_max);//atomic surf 1
    }

    for (int v = 2*nvibr; v < 3*nvibr; v++) //surface2 equations
    {
        dxdt(v) = (wD*x(v-2*nvibr)*st/width_r)*(1-N2_tot/a_max) - kdes2*x(v) - ka_2*Macheret(v-2*nvibr)*x(v)*(1-N2_tot/a_max) ;
        if (v!= 3*nvibr-1)
        {
            dxdt(v) += kvt_s*x(v+1);
        }
        if (v!= 2*nvibr)
        {
            dxdt(v) += -kvt_s*x(v);
        }
        dxdt(N-9) += ka_2*Macheret(v-2*nvibr)*x(v)*(1-N2_tot/a_max);//atomic_surf 2
    }

    dxdt(N-10) += -kNH_1*x(N-10)*surf_rate_1*x(N-7);//surf 1 chemistry
    dxdt(N-6) += kNH_1*x(N-10)*x(N-7)*surf_rate_1 - kNH2_1*x(N-6)*x(N-7)*surf_rate_1;
    dxdt(N-4) += kNH2_1*x(N-6)*x(N-7)*surf_rate_1 - kNH3_1*x(N-4)*x(N-7)*surf_rate_1;
    dxdt(N-2) += kNH3_1*x(N-4)*x(N-7)*surf_rate_1 - kNH3_des1*x(N-2)*surf_rate_1;

    dxdt(N-9) += -kNH_2*x(N-9)*surf_rate_2*x(N-7);//surf 2 chemistry
    dxdt(N-5) += kNH_2*x(N-9)*x(N-7)*surf_rate_2 - kNH2_2*x(N-5)*x(N-7)*surf_rate_2;
    dxdt(N-3) += kNH2_2*x(N-5)*x(N-7)*surf_rate_2 - kNH3_2*x(N-3)*x(N-7)*surf_rate_2;
    dxdt(N-1) += kNH3_2*x(N-3)*x(N-7)*surf_rate_2 - kNH3_des2*x(N-1)*surf_rate_2;

    dxdt(N-8) += surf_rate_1*(-x(N-8)*st*wDH2*(1.0 - N1_tot/a_max) + width_r*kH2_des*x(N-7)*x(N-7)) + //molecular hydrogen
        surf_rate_2*(-x(N-8)*st*wDH2*(1.0 - N2_tot/a_max) + width_r*kH2_des*x(N-7)*x(N-7));
    dxdt(N-7) +=  2*((wDH2*x(N-8)*st/width_r)*(surf_rate_1*(1-N1_tot/a_max)+surf_rate_2*(1-N2_tot/a_max))
    - kH2_des*x(N-7)*x(N-7)) - surf_rate_1*x(N-7)*(kNH_1*x(N-10)+kNH2_1*x(N-6)+kNH3_1*x(N-4)) -
    surf_rate_2*x(N-7)*(kNH_2*x(N-9)+kNH2_2*x(N-5)+kNH3_2*x(N-3)); //atomic hydrogen
    dxdt(N) = kNH3_des2*x(N-1)*surf_rate_2+ kNH3_des1*x(N-2)*surf_rate_1;
    dxdt /= Ng;
    std::cout<<x(N)/Ng<<"\n";
}
};



struct stiff_system_jacobi
{
    void operator() (const state_type &x_inp, matrix_type &J, const double /* t*/, state_type &dfdt )
    {
        state_type x = x_inp;
        x *= Ng;
        state_type subrange1;
        state_type subrange2;
        subrange1 = subrange(x,nvibr,2*nvibr);
        subrange2 = subrange(x,2*nvibr,3*nvibr);
        double N1_tot = sum(subrange1) + x(N-10)+x(N-6)+x(N-4)+x(N-2)+x(N-7)*surf_rate_1;//fast H diffusion is taken for account
        double N2_tot = sum(subrange2) + x(N-9)+x(N-5)+x(N-3)+x(N-1)+x(N-7)*surf_rate_2;

        for (int i = 0; i < N+1; i++)
        {
            dfdt(i) = 0; //explicit dependence of right hand side
            for (int j = 0; j < N+1; j++)
            {
                J(i,j) = 0;
            }
        }
        //gas
        for (int v = 0; v < nvibr; v++)
        {
            J(v,v) += -wD*surf_rate_1*st*(1-N1_tot/a_max) - wD*surf_rate_2*st*(1-N2_tot/a_max);
            if (v != nvibr-1)
            {
                J(v,v+1) += kvt(v+1)*Ng;
            }
            if (v != 0)
            {
                J(v,v) += -kvt(v)*Ng;
            }
            J(v,v+nvibr) += width_r*surf_rate_1*kdes1;
            J(v,v+2*nvibr) += width_r*surf_rate_2*kdes2;
            J(v,N-10) += wD*surf_rate_1*st*x(v)/(a_max);
            J(v,N-9) += wD*surf_rate_2*st*x(v)/(a_max);
            J(v,N-6) += wD*surf_rate_1*st*x(v)/(a_max);
            J(v,N-4) += J(v,N-6);
            J(v,N-2) += J(v,N-6);
            J(v,N-5) += wD*surf_rate_2*st*x(v)/(a_max);
            J(v,N-3) += J(v,N-5);
            J(v,N-1) += J(v,N-5);
            J(v,N-7) += J(v,N-5)+J(v,N-6);
            for (int u = 0; u<nvibr-1; u++)
            {
                if (v != nvibr-1)
                {
                    J(v,v) += -kvv(v,v+1,u+1,u)*x(u+1);
                    J(v,u+1) += -kvv(v,v+1,u+1,u)*x(v);
                    J(v,v+1) += kvv(v+1,v,u,u+1)*x(u);
                    J(v,u) += kvv(v+1,v,u,u+1)*x(v+1);
                }
                if (v != 0)
                {
                    J(v,v-1) += kvv(v-1,v,u+1,u)*x(u+1);
                    J(v,u+1) += kvv(v-1,v,u+1,u)*x(v-1);
                    J(v,v) += -kvv(v,v-1,u,u+1)*x(u);
                    J(v,u) += -kvv(v,v-1,u,u+1)*x(v);
                }

                if (u < v)
                {
                    J(v,u) += kev(u,v)*ne;
                    J(v,v) += -kev(u,v)*ne;
                }
                else if (u > v)
                {
                    J(v,v) += -kev(v,u)*ne;
                    J(v,u) += kev(v,u)*ne;
                }

                J(v,u+nvibr) += wD*surf_rate_1*st*x(v)/(a_max);
                J(v,u+2*nvibr) += wD*surf_rate_2*st*x(v)/(a_max);

            }
            J(v,2*nvibr-1) += wD*surf_rate_1*st*x(v)/(a_max);
            J(v,3*nvibr-1) += wD*surf_rate_2*st*x(v)/(a_max);

        }

        //molecular surface1
        for (int v = nvibr; v < 2*nvibr; v++)
        {
            J(v,v) += -kdes1 - ka_1*Macheret(v-nvibr)*(1-N1_tot/a_max);
            if (v != 2*nvibr - 1)
            {
                J(v,v+1) += kvt_s;
            }
            if (v != nvibr)
            {
                J(v,v) += -kvt_s;
            }
            J(v,v-nvibr) += (wD*st/width_r)*(1-N1_tot/a_max);
            J(v,N-10) += -wD*st*x(v-nvibr)/(width_r*a_max) + ka_1*Macheret(v-nvibr)*x(v)/(a_max);
            J(v,N-6) += J(v,N-10);
            J(v,N-4) += J(v,N-10);
            J(v,N-2) += J(v,N-10);
            J(v,N-7) += J(v,N-10)*surf_rate_1;
            for (int u = nvibr; u<2*nvibr; u++)
            {
                J(v,u) += -wD*st*x(v-nvibr)/(width_r*a_max)+ka_1*Macheret(v-nvibr)*x(v)/a_max;
            }
        }

        //molecular surface2
        for (int v = 2*nvibr; v < 3*nvibr; v++)
        {
            J(v,v) += -kdes2 - ka_2*Macheret(v-2*nvibr)*(1-N2_tot/a_max);
            if (v != 3*nvibr - 1)
            {
                J(v,v+1) += kvt_s;
            }
            if (v != 2*nvibr)
            {
                J(v,v) += -kvt_s;
            }
            J(v,v-2*nvibr) += wD*st/width_r*(1-N2_tot/a_max);
            J(v,N-9) += (-wD*st*x(v-2*nvibr)/(width_r*a_max) + ka_2*Macheret(v-2*nvibr)*x(v)/(a_max));
            J(v,N-5) += J(v,N-9);
            J(v,N-3) += J(v,N-9);
            J(v,N-1) += J(v,N-9);
            J(v,N-7) += J(v,N-9)*surf_rate_2;
            for (int u = 2*nvibr; u<3*nvibr; u++)
            {
                J(v,u) += -wD*st*x(v-2*nvibr)/(width_r*a_max)+ka_2*Macheret(v-2*nvibr)*x(v)/(a_max);
            }
        }

        //atomic surf1
        for (int v = nvibr; v < 2*nvibr; v++)
        {
            J(N-10,v) += ka_1*Macheret(v-nvibr);
            double add = -ka_1*Macheret(v-nvibr)*x(v)/a_max;
            J(N-10,N-10) += add;
            J(N-10,N-7) += add*surf_rate_1;
            J(N-10,N-6) += add;
            J(N-10,N-4) += add;
            J(N-10,N-2) += add;
            J(N-8,v) += surf_rate_1*x(N-8)*st*wDH2/a_max;
            J(N-7,v) += -2*surf_rate_1*st*wDH2*x(N-8)/(width_r*a_max);
            for (int u = nvibr; u< 2*nvibr; u++)
            {
                if (u == v)
                {
                    J(N-10,u) += 2*add;
                }
                else
                {
                    J(N-10,u) += add;
                }
            }
        }

        //atomic surf2
        for (int v = 2*nvibr; v < 3*nvibr; v++)
        {
            J(N-9,v) += ka_2*Macheret(v-2*nvibr);
            double add2 = -ka_2*Macheret(v-2*nvibr)*x(v)/a_max;
            J(N-9,N-9) += add2;
            J(N-9,N-7) += add2*surf_rate_2;
            J(N-9,N-5) += add2;
            J(N-9,N-3) += add2;
            J(N-9,N-1) += add2;
            J(N-8,v) +=  surf_rate_2*x(N-8)*st*wDH2/a_max;
            J(N-7,v) +=  -2*surf_rate_2*st*wDH2*x(N-8)/(width_r*a_max);
            for (int u = 2*nvibr; u< 3*nvibr; u++)
            {
                if (u == v)
                {
                    J(N-9,u) += 2*add2;
                }
                else
                {
                    J(N-9,u) += add2;
                }
            }
        }

        J(N-10,N-10) += -kNH_1*surf_rate_1*x(N-7);
        J(N-10,N-7) += -kNH_1*surf_rate_1*x(N-10);
        J(N-6,N-10) += kNH_1*x(N-7)*surf_rate_1;
        J(N-6,N-6) += -kNH2_1*x(N-7)*surf_rate_1;
        J(N-6,N-7) += kNH_1*x(N-10)*surf_rate_1 - kNH2_1*x(N-6)*surf_rate_1;
        J(N-4,N-7) += kNH2_1*x(N-6)*surf_rate_1 -kNH3_1*x(N-4)*surf_rate_1;
        J(N-4,N-6) += kNH2_1*x(N-7)*surf_rate_1;
        J(N-4,N-4) += -kNH3_1*x(N-7)*surf_rate_1;
        J(N-2,N-4) += kNH3_1*x(N-7)*surf_rate_1;
        J(N-2,N-7) += kNH3_1*x(N-4)*surf_rate_1;
        J(N-2,N-2) += -kNH3_des1*surf_rate_1;

        J(N-9,N-9) += -kNH_2*surf_rate_2*x(N-7);
        J(N-9,N-7) += -kNH_2*surf_rate_2*x(N-9);
        J(N-5,N-9) += kNH_2*x(N-7)*surf_rate_2;
        J(N-5,N-5) += -kNH2_2*x(N-7)*surf_rate_2;
        J(N-5,N-7) += kNH_2*x(N-9)*surf_rate_2 - kNH2_2*x(N-5)*surf_rate_2;
        J(N-3,N-7) += kNH2_2*x(N-5)*surf_rate_2 -kNH3_2*x(N-3)*surf_rate_2;
        J(N-3,N-5) += kNH2_2*x(N-7)*surf_rate_2;
        J(N-3,N-3) += -kNH3_2*x(N-7)*surf_rate_2;
        J(N-1,N-3) += kNH3_2*x(N-7)*surf_rate_2;
        J(N-1,N-7) += kNH3_2*x(N-3)*surf_rate_2;
        J(N-1,N-1) += -kNH3_des2*surf_rate_2;

        J(N-8,N-8) += -surf_rate_1*st*wDH2*(1-N1_tot/a_max) - surf_rate_2*st*wDH2*(1-N2_tot/a_max);
        J(N-8,N-6) += surf_rate_1*st*wDH2*x(N-8)/a_max;
        J(N-8,N-4) += J(N-8,N-6);
        J(N-8,N-2) += J(N-8,N-6);
        J(N-8,N-10) += J(N-8,N-6);
        J(N-8,N-5) += surf_rate_2*st*wDH2*x(N-8)/a_max;
        J(N-8,N-3) += J(N-8,N-5);
        J(N-8,N-1) += J(N-8,N-5);
        J(N-8,N-9) += J(N-8,N-5);
        J(N-8,N-7) += surf_rate_1*J(N-8,N-6)+surf_rate_2*J(N-8,N-5)+2*width_r*kH2_des*x(N-7);

        J(N-7,N-7) += -2*x(N-8)*st*wDH2/(a_max*width_r) - 4*x(N-7)*kH2_des-surf_rate_1*(kNH_1*x(N-10)+kNH2_1*x(N-6)+kNH3_1*x(N-4)) -
    surf_rate_2*(kNH_2*x(N-9)+kNH2_2*x(N-5)+kNH3_2*x(N-3));
        J(N-7,N-8) += 2*wDH2*st/width_r*(surf_rate_1*(1-N1_tot/a_max)+surf_rate_2*(1-N2_tot/a_max));
        J(N-7,N-6) += -2*surf_rate_1*x(N-8)*wDH2*st/(width_r*a_max)-surf_rate_1*x(N-7)*kNH2_1;
        J(N-7,N-4) += -2*surf_rate_1*x(N-8)*wDH2*st/(width_r*a_max)-surf_rate_1*x(N-7)*kNH3_1;
        J(N-7,N-10) += -2*surf_rate_1*x(N-8)*wDH2*st/(width_r*a_max)-surf_rate_1*x(N-7)*kNH_1;
        J(N-7,N-5) += -2*surf_rate_2*x(N-8)*wDH2*st/(width_r*a_max)-surf_rate_2*x(N-7)*kNH2_2;
        J(N-7,N-3) += -2*surf_rate_2*x(N-8)*wDH2*st/(width_r*a_max)-surf_rate_2*x(N-7)*kNH3_2;
        J(N-7,N-9) += -2*surf_rate_2*x(N-8)*wDH2*st/(width_r*a_max)-surf_rate_2*x(N-7)*kNH_2;
        J(N,N-1) += kNH3_des2*surf_rate_2;
        J(N,N-2) += kNH3_des1*surf_rate_1;

        //J/=Ng;
    }
};

void calc_jacobi(const state_type &x) //numerical Jacobian
{
    double del = 1e-9;
    double t = 0;
    int col = 0;
    state_type dxdt (N+1);
    state_type dxdt1 (N+1);
    state_type dfdt (N+1);
    matrix_type J (N+1,N+1);
    matrix_type J_num (N+1,N+1);
    state_type x1 = x;
    for(int j=0 ; j<N+1 ; j++){
        double dx = x1(j)*del;
        if(dx == 0.0)
            dx = 1e-9;
        x1(j) -= dx;
        stiff_system()(x1,dxdt,t);
        x1(j) += 2.0*dx;
        stiff_system()(x1,dxdt1,t);
        state_type der = (-dxdt+dxdt1)/(2.0*dx);
        for(int i=0 ; i<N+1 ; i++){
            J_num(i,j) = der(i);
        }
    }
    stiff_system_jacobi()(x,J,t,dfdt);
    double err_max = 0.0;
    int i_max = 0;
    int j_max = 0;
    for (int i=0; i<N+1; i++)
    {
        for (int j=0; j<N+1; j++){
            double err = fabs(J(i,j) - J_num(i,j));
            double numer = fabs(J(i,j) + J_num(i,j));
            if(numer > 0)
                err /= numer;
            if(err > err_max ){
                err_max = err;
                i_max = i;
                j_max = j;
            }
        }
    }
    //std::cout << J_num <<std::endl;
    std::cout<< "Maximum J diff = " << err_max << " at i = " << i_max << "  j = " << j_max << " ( J = " << J(i_max, j_max) <<
            ", J_num  = " << J_num(i_max, j_max) << ")" <<std::endl;
    for (int i = 0;i<N+1;i++){
    for (int j = 0;j<N+1;j++){
    std::cout<<i<<"  "<<j<<" "<<J(i,j)<<" "<<J_num(i,j)<<std::endl;
    }
    }
};


struct streaming_observer
{
    std::ostream &m_out;
    streaming_observer( std::ostream &out ) : m_out( out ) { }
    void operator()( const state_type &x , double t ) const
    {

        m_out << t;
        for( size_t i=0 ; i<x.size() ; ++i )
        { m_out <<" "<< x[i];
        }
        m_out << "\n";
    }
};


int main(){
    state_type x (N+1);
    for (int i=0;i<N+1;i++)
    {
    x(i) = 1e-12*((double) rand() / (RAND_MAX));
    }

    x(0) = 0.75;
    x(N-8) = 0.25;
    //calc_jacobi(x);
    //exit(0);
    std::ofstream out2("treanor.dat");

    size_t num_of_steps =  integrate_adaptive( make_dense_output< rosenbrock4< double > >(1e-9,1e-3) ,
                                           std::make_pair( stiff_system() , stiff_system_jacobi() ) ,
                                           x ,0.0 , 1e-2 , 1e-9 ,
                                           streaming_observer(std::cout));

    //Treanor distribution
    for (int i = 0; i<nvibr; i++)
    {
        out2 << i << " "<<log10(x(i)*Ng)<<"\n";
    }
    out2.close();
    return 0;
}
