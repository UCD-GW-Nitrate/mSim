//--------------------------------------------------------
// -----------define structures and classes
//--------------------------------------------------------
// ind takes i ,j and number of rows and returns the linear index
#define ind(i,j,ni) (i-1+ni*(j-1))
const double PI = 4*atan(1);

struct partoptions{
	int Nal;// = 1000; 		1 
	int search;// = 50; 	2
	bool bed_corr;// = 0; 	3
	int Knodes;// = 0; 		4
	double step;// = 0.5;	5
 	double minstep;// = 0.1;6
	double errmin;// = 0.01;7
	double errmax;// = 0.1; 8
    int method;// = 1->Euler, 2->RK2, 3->RK4, 4->RK45; 9
    int pornodes;// =0;		10
    int porinterp;	//		11
    int Ngen;// = 15;		12
    double maxstep;// = 0.5;//times the maximum diagonal of the current element 13
    int stall_times;// = 100;     				14
    int el_type; // = 1-> prism, 2-> hex		15
    int el_order; // = 1;						16
    bool valid;		//							17
    int Nel; //Number of elements				18
    int Np; //Number of points					19
    int Nsides; // Number of element faces		20	
    int Nsh; // Number of shape functions		21
    int Nlay; // Number of layers				22
};

class particle_point
{
	public:
	double coord[3]; //coordinates
	int idel; //element id
	int lay;// layer id
	double t[3]; // parametric coordinates in xy
	double u; // parametric coordinates in z
	double vel[3]; // velocity 
	double oldvel[3];// previous velocity
	double step;
	double maxstep; // maximum allowable step for current element
	
	void initialize(double stp){
        coord[0]=0; coord[1]=0; coord[2]=0;
        idel=-88;
        lay=-88;
        t[0]=-88;t[1]=-88;t[2]=-88;
        u=-88;
        vel[0]=0; vel[1]=0; vel[2]=0;
        oldvel[0]=0; oldvel[1]=0; oldvel[2]=0;
        step=stp;
    }
    
    void setcoord(double CRD[]){
		coord[0]=CRD[0];
		coord[1]=CRD[1];
		coord[2]=CRD[2];
	}
	
	void setoldvel(double vel[]){
		oldvel[0]=vel[0];
		oldvel[1]=vel[1];
		oldvel[2]=vel[2];
	}
	
	void show(){
        std::cout <<"---------------------------------------------------------------\n";
        std::cout << "coords: " << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
        std::cout << "idel:" << idel << "\n";
        std::cout << "lay:" << lay << "\n";
        std::cout << "t: " << t[0] << " " << t[1] << " " << t[2] << "\n";
        std::cout << "u:" << u << "\n";
        std::cout << "vel (x1000): " << vel[0]*1000 << " " << vel[1]*1000 << " " << vel[2]*1000 << "\n";
        std::cout << "prevel: " << oldvel[0] << " " << oldvel[1] << " " << oldvel[2] << "\n";
        std::cout << "step:" << step << "\n";
        std::cout << "maxstep:" << maxstep << "\n";
    }
    
    bool isPntinElem(int elementtype){
		if (elementtype==1){
			if (t[0] <0 || t[0] >1 || t[1] <0 || t[1] >1 || t[2] <0 || t[2] >1)
				return false;
			else
				return true;
		}
		return false; // Just to remove the warning
    }
    int checkpoint(){
        int status = 0;
        if (idel == -88)
			status = -8; // initial position is outside of the domain
        if (idel ==0)
            status=2; //Particle exit from the side
        if (lay==0)
            status=1; //Particle exits from the top
        if (lay==-9)
            status=-9; //Particle exits from the bottom
        return status;
    }
};
/*
    The class my_comp and function comp_func are used for elemet sorting
    They are used to get a similar output to the matlab command [c d]=sort(x)
*/
class my_comp
{
    public:
    double num;
    int ind;
};
bool comp_func (const my_comp& a, const my_comp& b)
{
	return a.num < b.num;
}

//--------------------------------------------------------
// ---------End of Definitions
//--------------------------------------------------------


//----------------------------------------------
// --------Utility functions
//---------------------------------------
void check_sVolume(int& cnt_stall_times, double sV[], particle_point pnt){
	bool expand=false;
	if (pnt.coord[0] > sV[1]){
		sV[1] = pnt.coord[0];
		expand=true;
	}
	if (pnt.coord[0] < sV[0]){
		sV[0] = pnt.coord[0];
		expand=true;
	}
	if (pnt.coord[1] > sV[3]){
		sV[3] = pnt.coord[1];
		expand=true;
	}
	if (pnt.coord[1] < sV[2]){
		sV[2] = pnt.coord[1];
		expand=true;
	}
	if (pnt.coord[2] > sV[5]){
		sV[5] = pnt.coord[2];
		expand=true;
	}
	if (pnt.coord[2] < sV[4]){
		sV[4] = pnt.coord[2];
		expand=true;
	}
	if (expand==true)
		cnt_stall_times=0;
	else
		cnt_stall_times++;
}


double check_velAngle(particle_point pnt){
	double angle;
	angle = (180.0/PI)*acos((pnt.vel[0]*pnt.oldvel[0]+
					        pnt.vel[1]*pnt.oldvel[1]+
					        pnt.vel[2]*pnt.oldvel[2])/
			(sqrt(pow(pnt.vel[0],2)+pow(pnt.vel[1],2)+pow(pnt.vel[2],2))*
			 sqrt(pow(pnt.oldvel[0],2)+pow(pnt.oldvel[1],2)+pow(pnt.oldvel[2],2))));
	return angle;
}

