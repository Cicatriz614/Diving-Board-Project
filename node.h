/****************************NODE CLASS****************************************/
#define TRUE 1
#define FALSE 0

struct Connection
{
      public:
        int tonode;
        double area;
        double length;
};

class Node {
      private:
              bool *Fixed; //only boundary condition used in part 4
              double *initialstate; //Initial position in global space
              double *curstate, *paststate; //displacement or delta-charge, current and previous
              double *selffactors;
              double momentofinertia;
              double areamoment;
//              double nodelength;
              double nodeheight;
              int DegreesFreedom;
              int conncount;
              Connection *conn; //will hold which nodes are connected to each other
      public:
              void init(int,int); //this function will have the nodes either hard-coded or read from a CSV
              void initProb4(int);
              void initProb(Node*, int, double, double, double, double, double, double, double, double,
                            int, double, double, double, double,int*);
              void initNodeCell(int, bool, bool, double, int, double, double, double, double, double,
                                double, double, double, int, double, int, double);
              void UpdateDelta(double*);
              double GetCurState(int);
              double GetInitState(int);
              double GetPosition(int);
              double GetMoment();
              double GetAreaMoment();
              bool isFixed();
              bool isFixed(int);
              int getconn();
              int connto(int);
              double get_2D_Distance(Node,int,int);
              double get_init_length(Node);
              double getFrac(Node,int,int);
              double getselffactor(int);
              //double getconnfactor(int,int);
              double GetConArea(int);
              void initProb3(int);
              void gettransformation(Node,double**);
              double tanInv(Node);
              void formDampMatrix(double**,double*,Node,int,double,double);
              double getInitialstate(int);
              double getCurstate(int);
              double getconnlength(int);
};

/***********************CLASS FUNCTIONS***********************************/
double Node::GetMoment()
{
    return momentofinertia;
}

double Node::getconnlength(int index)
{
    return conn[index].length;
}
double Node::GetAreaMoment()
{
    return areamoment;
}
void Node::init(int index, int DoF)
{
    DegreesFreedom = DoF;
    curstate = new double[DoF];
    paststate = new double[DoF];
    initialstate = new double[DoF];
    Fixed = new bool[DoF];
    selffactors = new double[3]; //proportional, derivative and second derivative term
    std::cout<< "Enter initial conditions for Node " << index << std::endl;
    for(int i = 0;i<DegreesFreedom;i++)
    {
        std::cout << "Initial Value (ex. charge in circuit, position in mds system) for Degree " << (i+1) << ": ";
        std::cin >> initialstate[i] ;
        paststate[i] = 0;
        curstate[i] = 0;
        std::cout<<std::endl<<"Is this degree fixed (1 if true, 0 if false)?";
        std::cin >> Fixed[i];
    }
    std::cout<<"Factor for self-referential proportional   ";
    std::cin>>selffactors[0];
    std::cout<<"Factor for self-referential derivative     ";
    std::cin>>selffactors[1];
    std::cout<<"Factor for self-referential 2nd derivative ";
    std::cin>>selffactors[2];
    std::cout<<"How many higher indexed nodes is this node connected to?"<<std::endl;
    std::cin>>conncount;
    conn = new Connection[conncount];
    /*for(int i = 0;i<conncount;i++)
    {
        std::cout<<"Input node index for connection"<<(i+1)<<std::endl;
        std::cin>>conn[i].tonode;
        std::cout<<"Factor for coupled proportional   ";
        std::cin>>conn[i].factors[0];
        std::cout<<"Factor for coupled derivative     ";
        std::cin>>conn[i].factors[1];
        std::cout<<"Factor for coupled 2nd derivative ";
        std::cin>>conn[i].factors[2];
    }*/
    std::cout<<std::endl;
}

void Node:: UpdateDelta(double *newPos)
{
     for(int i = 0;i<DegreesFreedom;i++)
     {
         paststate[i] =  curstate[i]; //pos[][0] is the old displacement
         curstate[i] = newPos[i];
     }
     return;
}

double Node::GetInitState(int index)
{
    return initialstate[index];
}

double Node::GetCurState(int index)
{
    return curstate[index];
}

double Node::GetPosition(int index)
{
    return curstate[index]+initialstate[index];
}
bool Node::isFixed(int index)
{
     return Fixed[index];
}
bool Node::isFixed()
{
    bool temp = TRUE;
    for(int i = 0;i<DegreesFreedom;i++)
    {
        temp = temp&Fixed[i];
    }
    return temp;
}
//to get distance if degrees describe physical orientation.
double Node::get_2D_Distance(Node other, int index1, int index2)
{
    double sum;
    sum = (((*this).GetPosition(index1))-(other.GetPosition(index1)))*(((*this).GetPosition(index1))-(other.GetPosition(index1)));
    sum += (((*this).GetPosition(index2))-(other.GetPosition(index2)))*(((*this).GetPosition(index2))-(other.GetPosition(index2)));
    return sqrt(sum);
}

double Node::get_init_length(Node other)
{
    double sum;
    sum = (initialstate[0]-other.initialstate[0])*(initialstate[0]-other.initialstate[0]);
    sum += (initialstate[1]-other.initialstate[1])*(initialstate[1]-other.initialstate[1]);
    return sqrt(sum);
}

//this gets the fraction that one degree affects another
//If the system being solved was 2-D. X and Y would be index
//0 and 1 respectively. Cos is delta-X over length.
//This would be specified by calling Node1.getFrac(Node2,0,1)
//In a similar manner, sin would be  Node1.getFrac(Node2,1,0)
//main code needs to check degree of freedom to decide whether or not this
//function is necessary
double Node::getFrac(Node other, int index1, int index2)
{
    double delta = other.GetPosition(index1)-(*this).GetPosition(index1);
    double dis = (*this).get_2D_Distance(other,index1,index2);
    return delta/dis;
}

int Node::getconn()
{
    return conncount;
}

int Node::connto(int index)
{
    return conn[index].tonode;
}

double Node::GetConArea(int index)
{
    return conn[index].area;
}
/*
double Node:: getconnfactor(int index,int select)
{
//    return conn[index].factors[select];
}
*/
double Node:: getselffactor(int index)
{
    return  selffactors[index];
}

void Node::initProb(Node Nodes[], int ribs, double twall, double trib, double l1,
                    double l2, double h1, double h2, double h3, double width, int noofnodes,
                    double fulcrum, double stiff, double damp, double density, int*nodesfixed)
{
    double sl = (l1 + l2)/(noofnodes-1);
         //  sl = floor(sl*1000000)/1000000;
    double slope1 = (h2-h1)/l1;
         //  slope1 = floor(slope1*1000000)/1000000;
    double slope2 = (h3-h2)/l2;
         //  slope2 = floor(slope2*1000000)/1000000;
    double vm = 0;
    double hm = 0;
    double furthest_edge = 0;
    double moment = 0;
    double height = 0;
    double x = 0;
    double y = 0;
    int current_node = 1;
    bool fulcrum_done = 0;

    //First node (fixed node)
    DegreesFreedom = 3;
    curstate = new double[3];
    paststate = new double[3];
    initialstate = new double[3];
    Fixed = new bool[3];
    selffactors = new double[3]; //proportional, derivative and second derivative term
    for(int i = 0;i<DegreesFreedom;i++)
    {
        initialstate[i] = 0;
        paststate[i] = 0;
        curstate[i] = 0;
    }
    Fixed[0] = 1;
    Fixed[1] = 1;
    Fixed[2] = 1;
    conncount = 1;
    conn = new Connection[1];
    nodeheight = (h1+slope1*sl/2);
    conn[0].area = (2*twall*width) + (2*twall+ribs*trib)*(nodeheight-2*twall);
    conn[0].tonode = 1;
    conn[0].length = sl;
    //nodeheight = floor(nodeheight*1000000)/1000000;
    vm = ((nodeheight-2*twall)*(ribs*trib+2*twall)*(sl/2) + (nodeheight-2*twall)*(width-ribs*trib-2*twall)*(twall)) *(density);
    hm = 2*(width)*(twall)*(sl/2)*(density);
    selffactors[2] = vm+hm;
    areamoment = ((width*nodeheight*nodeheight*nodeheight)-(width-ribs*trib-2*twall)*(nodeheight-2*twall)*(nodeheight-2*twall)*(nodeheight-2*twall))/12;
    momentofinertia = (hm/12)*width*width + (hm*(nodeheight-twall)*(nodeheight-twall)/4) + (vm/12)*(2*twall+ribs*trib)*(2*twall+ribs*trib) + (width-2*twall-ribs*trib)*(width-2*twall-ribs*trib)*(width-2*twall-ribs*trib)*(height-2*twall)*twall*density/12;
    furthest_edge += sl/2;

    while(furthest_edge < l1)
    {
        if(furthest_edge >= fulcrum && fulcrum_done == 0)
        {
             height = Nodes[current_node-1].nodeheight + slope1*sl;
             x = Nodes[current_node-1].initialstate[0] + sl;
             y = (h1-height)/2;
             vm = (height-2*twall)*(ribs*trib+2*twall)*(sl)*(density);
             hm = 2*(width)*(twall)*(sl)*(density);

             moment = (hm/12)*width*width + (hm*(height-twall)*(height-twall)/4) + (vm/12)*(2*twall+ribs*trib)*(2*twall+ribs*trib);
             Nodes[current_node].initNodeCell(3,0,1,vm+hm,1,stiff,damp,moment,x,y,height,twall,trib,ribs,width,current_node,sl);
             nodesfixed[3] = 3*(current_node-1) + 1;
             furthest_edge += sl;
             current_node++;
             fulcrum_done = 1;
             continue;
        }
        height = Nodes[current_node-1].nodeheight + slope1*sl;
        x = Nodes[current_node-1].initialstate[0] + sl;
        y = (h1-height)/2;
        vm = (height-2*twall)*(ribs*trib+2*twall)*(sl)*(density);
        hm = 2*(width)*(twall)*(sl)*(density);

        moment = (hm/12)*width*width + (hm*(height-twall)*(height-twall)/4) + (vm/12)*(2*twall+ribs*trib)*(2*twall+ribs*trib);
        Nodes[current_node].initNodeCell(3,0,0,vm+hm,1,stiff,damp,moment,x,y,height,twall,trib,ribs,width,current_node,sl);

        furthest_edge += sl;
        current_node++;
    }

    while(furthest_edge < l2 + l1)
    {
        if(furthest_edge >= fulcrum && fulcrum_done == 0)
        {
             height = Nodes[current_node-1].nodeheight + slope2*sl;
             x = Nodes[current_node-1].initialstate[0] + sl;
             y = (h1-height)/2;
             vm = (height-2*twall)*(ribs*trib+2*twall)*(sl)*(density);
             hm = 2*(width)*(twall)*(sl)*(density);

             moment = (hm/12)*width*width + (hm*(height-twall)*(height-twall)/4) + (vm/12)*(2*twall+ribs*trib)*(2*twall+ribs*trib);
             Nodes[current_node].initNodeCell(3,0,1,vm+hm,1,stiff,damp,moment,x,y,height,twall,trib,ribs,width,current_node,sl);
             nodesfixed[3] = 3*(current_node-1) + 1;
             furthest_edge += sl;
             current_node++;
             fulcrum_done = 1;
             continue;
        }
        if((furthest_edge + sl) > (l2 + l1))
        {
             height = ((Nodes[current_node-1].nodeheight + slope2*sl/2) + h3)/2;
             x = Nodes[current_node-1].initialstate[0] + sl;
             y = (h1-height)/2;
             vm = ((height-2*twall)*(ribs*trib+2*twall)*(sl/2) + (height-2*twall)*(width-ribs*trib-2*twall)*(twall)) *(density);
             hm = 2*(width)*(twall)*(sl)*(density);

             moment = (hm/12)*width*width + (hm*(height-twall)*(height-twall)/4) + (vm/12)*(2*twall+ribs*trib)*(2*twall+ribs*trib) + (1/12)*(width-2*twall-ribs*trib)*(width-2*twall-ribs*trib)*(width-2*twall-ribs*trib)*(height-2*twall)*twall*density;
             Nodes[current_node].initNodeCell(3,0,0,vm+hm,1,stiff,damp,moment,x,y,height,twall,trib,ribs,width,current_node,sl);
             Nodes[current_node].conncount = 0;
             break;
        }
        height = Nodes[current_node-1].nodeheight + slope2*sl;
        x = Nodes[current_node-1].initialstate[0] + sl;
        y = (h1-height)/2;
        vm = (height-2*twall)*(ribs*trib+2*twall)*(sl)*(density);
        hm = 2*(width)*(twall)*(sl)*(density);

        moment = (hm/12)*width*width + (hm*(height-twall)*(height-twall)/4) + (vm/12)*(2*twall+ribs*trib)*(2*twall+ribs*trib);
        Nodes[current_node].initNodeCell(3,0,0,vm+hm,1,stiff,damp,moment,x,y,height,twall,trib,ribs,width,current_node,sl);
        furthest_edge += sl;
        current_node++;
    }
}



void Node::initNodeCell(int DoF, bool fixed_x, bool fixed_y, double mass, int NoC, double stiff,
                        double damp, double moment, double position_x, double position_y, double height,
                        double wall_thickness, double rib_thickness, int ribs, double width, int currentnode, double sl)
{
    DegreesFreedom = DoF;
    curstate = new double[DoF];
    paststate = new double[DoF];
    initialstate = new double[DoF];
    Fixed = new bool[DoF];
    selffactors = new double[3]; //proportional, derivative and second derivative term
    initialstate[0] = position_x;
    paststate[0] = 0;
    curstate[0] = 0;
    initialstate[1] = position_y;
    paststate[1] = 0;
    curstate[1] = 0;
    Fixed[0] = fixed_x;
    Fixed[1] = fixed_y;
    Fixed[2] = 0;
    selffactors[2] = mass;
    conncount = NoC;
    conn = new Connection[conncount];
    momentofinertia = moment;
    areamoment = (width*nodeheight*nodeheight*nodeheight-(width-ribs*rib_thickness-2*wall_thickness)*(nodeheight-2*wall_thickness)*(nodeheight-2*wall_thickness)*(nodeheight-2*wall_thickness))/12;

    nodeheight = height;
    conncount = 1;
    conn = new Connection[1];
    conn[0].area = (2*wall_thickness*width) + (2*wall_thickness+ribs*rib_thickness)*(nodeheight-2*wall_thickness);
    conn[0].tonode = currentnode+1;
    conn[0].length = sl;
}

double Node::getInitialstate(int index)
{
    return initialstate[index];
}

double Node::getCurstate(int index)
{
    return curstate[index];
}

double Node::tanInv(Node next)
{
    double x = (next.getInitialstate(0)+next.getCurstate(0)) - (initialstate[0]+curstate[0]);
    double y = (next.getInitialstate(1)+next.getCurstate(1)) - (initialstate[1]+curstate[1]);
    if(x < 0 && y > 0)
    {
        return pi + atan(y/x);
    }
    if(x < 0 && y < 0)
    {
        return atan(y/x) - pi;
    }
    if(x == 0 && y > 0)
    {
        return pi/2;
    }
    if(x == 0 && y < 0)
    {
        return -pi/2;
    }
    if(x == 0 && y == 0)
    {
        return NAN;
    }
    return atan(y/x);
}

void Node::formDampMatrix(double**damp, double*CurDisplacement, Node next, int node, double width, double linfag)
{
    double c = cos(((*this).tanInv(next)) + CurDisplacement[3*node + 2]);
    double c2 = c*c;
    double s = sin((*this).tanInv(next) + CurDisplacement[3*node + 2]);
    double s2 = s*s;
    double b = 1.33917*(*this).getconnlength(0)*width;
    damp[0][0] = linfag*b*c2;
    damp[0][1] = linfag*b*c*s;
    damp[0][2] = 0;
    damp[1][0] = -linfag*b*c*s;
    damp[1][1] = linfag*b*s2;
    damp[1][2] = 0;
    damp[2][0] = 0;
    damp[2][1] = 0;
    if(node == 0)
        damp[2][2] = 100;
    else
        damp[2][2] = 0;

    return;
}
