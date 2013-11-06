/****************************NODE CLASS****************************************/
#define TRUE 1
#define FALSE 0

struct Connection
{
      public:
        int tonode;
        double area;
};

class Node {
      private:
              bool *Fixed; //only boundary condition used in part 4
              double *initialstate; //Initial position in global space
              double *curstate, *paststate; //displacement or delta-charge, current and previous
              double *selffactors;
              double momentofinertia;
              double nodelength;
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
                                double, double, double, int, double, int);
              void UpdateDelta(double*);
              double GetCurState(int);
              double GetInitState(int);
              double GetPosition(int);
              double GetMoment();
              bool isFixed();
              bool isFixed(int);
              int getconn();
              int connto(int);
              double get_2D_Distance(Node,int,int);
              double get_init_length(Node,int,int);
              double getFrac(Node,int,int);
              double getselffactor(int);
              double getconnfactor(int,int);
              double GetConArea(int);
              void initProb3(int);
              void gettransformation(Node,double**);
};

/***********************CLASS FUNCTIONS***********************************/
double Node::GetMoment()
{
    return momentofinertia;
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

double Node::get_init_length(Node other, int index1, int index2)
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

double Node:: getconnfactor(int index,int select)
{
//    return conn[index].factors[select];
}

double Node:: getselffactor(int index)
{
    return  selffactors[index];
}

void Node::initProb(Node Nodes[], int ribs, double wall_thickness, double rib_thickness, double length1,
                    double length2, double height1, double height2, double height3, double width, int noofnodes,
                    double fulcrum, double stiff, double damp, double density, int*nodesfixed)
{
    double standard_length = (length1 + length2)/(noofnodes-1);
         //  standard_length = floor(standard_length*1000000)/1000000;
    double slope1 = (height2-height1)/length1;
         //  slope1 = floor(slope1*1000000)/1000000;
    double slope2 = (height3-height2)/length2;
         //  slope2 = floor(slope2*1000000)/1000000;
    double vertical_mass = 0;
    double horizontal_mass = 0;
    double furthest_edge = 0;
    double moment = 0;
    double height = 0;
    double position_x = 0;
    double position_y = 0;
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
    Fixed[2] = 0;
    conncount = 1;
    conn = new Connection[1];
    nodeheight = (height1+slope1*standard_length/2);
    conn[0].area = (2*wall_thickness*width) + (2*wall_thickness+ribs*rib_thickness)*(nodeheight-2*wall_thickness);
    conn[0].tonode = 1;
    //nodeheight = floor(nodeheight*1000000)/1000000;
    vertical_mass = ((nodeheight-2*wall_thickness)*(ribs*rib_thickness+2*wall_thickness)*(standard_length/2) + (nodeheight-2*wall_thickness)*(width-ribs*rib_thickness-2*wall_thickness)*(wall_thickness)) /(density);
    horizontal_mass = 2*(width)*(wall_thickness)*(standard_length/2)/(density);
    momentofinertia = (1/12)*(vertical_mass)*((2*wall_thickness+ribs*rib_thickness)*(2*wall_thickness+ribs*rib_thickness)+(nodeheight-2*wall_thickness)*(nodeheight-2*wall_thickness)) + (1/12)*(wall_thickness*(width-2*wall_thickness-ribs*rib_thickness)*(nodeheight-2*wall_thickness)/density)*((width-2*wall_thickness-ribs*rib_thickness)*(width-2*wall_thickness-ribs*rib_thickness)+(nodeheight-2*wall_thickness)*(nodeheight-2*wall_thickness)) + (1/6)*(horizontal_mass)*((width)*(width)+(wall_thickness)*(wall_thickness))+(horizontal_mass)*((nodeheight+wall_thickness)/2);
    selffactors[2] = vertical_mass + horizontal_mass + wall_thickness*(width-2*wall_thickness-ribs*rib_thickness)*(nodeheight-2*wall_thickness)/density;
    furthest_edge += standard_length/2;

    while(furthest_edge < length1)
    {
        if(furthest_edge >= fulcrum && fulcrum_done == 0)
        {
             height = Nodes[current_node-1].nodeheight + slope1*standard_length;
             position_x = Nodes[current_node-1].initialstate[0] + standard_length;
             position_y = (height1-height)/2;
             vertical_mass = (height-2*wall_thickness)*(ribs*rib_thickness+2*wall_thickness)*(standard_length)/(density);
             horizontal_mass = 2*(width)*(wall_thickness)*(standard_length)/(density);
             moment = (1/12)*(vertical_mass)*((2*wall_thickness+ribs*rib_thickness)*(2*wall_thickness+ribs*rib_thickness)+(height-2*wall_thickness)*(height-2*wall_thickness)) + (1/6)*(horizontal_mass)*((width)*(width)+(wall_thickness)*(wall_thickness))+(horizontal_mass)*((height+wall_thickness)/2);
             Nodes[current_node].initNodeCell(3,0,1,vertical_mass+horizontal_mass,1,stiff,damp,moment,position_x,position_y,height,wall_thickness,rib_thickness,ribs,width,current_node);
             nodesfixed[2] = current_node;
             furthest_edge += standard_length;
             current_node++;
             fulcrum_done = 1;
             continue;
        }
        height = Nodes[current_node-1].nodeheight + slope1*standard_length;
        position_x = Nodes[current_node-1].initialstate[0] + standard_length;
        position_y = (height1-height)/2;
        vertical_mass = (height-2*wall_thickness)*(ribs*rib_thickness+2*wall_thickness)*(standard_length)/(density);
        horizontal_mass = 2*(width)*(wall_thickness)*(standard_length)/(density);
        moment = (1/12)*(vertical_mass)*((2*wall_thickness+ribs*rib_thickness)*(2*wall_thickness+ribs*rib_thickness)+(height-2*wall_thickness)*(height-2*wall_thickness)) + (1/6)*(horizontal_mass)*((width)*(width)+(wall_thickness)*(wall_thickness))+(horizontal_mass)*((height+wall_thickness)/2);
        Nodes[current_node].initNodeCell(3,0,0,vertical_mass+horizontal_mass,1,stiff,damp,moment,position_x,position_y,height,wall_thickness,rib_thickness,ribs,width,current_node);

        furthest_edge += standard_length;
        current_node++;
    }

    while(furthest_edge < length2 + length1)
    {
        if(furthest_edge >= fulcrum && fulcrum_done == 0)
        {
             height = Nodes[current_node-1].nodeheight + slope2*standard_length;
             position_x = Nodes[current_node-1].initialstate[0] + standard_length;
             position_y = (height1-height)/2;
             vertical_mass = (height-2*wall_thickness)*(ribs*rib_thickness+2*wall_thickness)*(standard_length)/(density);
             horizontal_mass = 2*(width)*(wall_thickness)*(standard_length)/(density);
             moment = (1/12)*(vertical_mass)*((2*wall_thickness+ribs*rib_thickness)*(2*wall_thickness+ribs*rib_thickness)+(height-2*wall_thickness)*(height-2*wall_thickness)) + (1/6)*(horizontal_mass)*((width)*(width)+(wall_thickness)*(wall_thickness))+(horizontal_mass)*((height+wall_thickness)/2);
             Nodes[current_node].initNodeCell(3,0,1,vertical_mass+horizontal_mass,1,stiff,damp,moment,position_x,position_y,height,wall_thickness,rib_thickness,ribs,width,current_node);
             nodesfixed[2] = current_node;
             furthest_edge += standard_length;
             current_node++;
             fulcrum_done = 1;
             continue;
        }
        height = Nodes[current_node-1].nodeheight + slope2*standard_length;
        position_x = Nodes[current_node-1].initialstate[0] + standard_length;
        position_y = (height1-height)/2;
        vertical_mass = (height-2*wall_thickness)*(ribs*rib_thickness+2*wall_thickness)*(standard_length)/(density);
        horizontal_mass = 2*(width)*(wall_thickness)*(standard_length)/(density);
        moment = (1/12)*(vertical_mass)*((2*wall_thickness+ribs*rib_thickness)*(2*wall_thickness+ribs*rib_thickness)+(height-2*wall_thickness)*(height-2*wall_thickness)) + (1/6)*(horizontal_mass)*((width)*(width)+(wall_thickness)*(wall_thickness))+(horizontal_mass)*((height+wall_thickness)/2);
        Nodes[current_node].initNodeCell(3,0,0,vertical_mass+horizontal_mass,1,stiff,damp,moment,position_x,position_y,height,wall_thickness,rib_thickness,ribs,width,current_node);
        furthest_edge += standard_length;
        current_node++;
        if(furthest_edge + standard_length > length2 + length1)
        {
             height = ((Nodes[current_node-1].nodeheight + slope2*standard_length/2) + height3)/2;
             position_x = Nodes[current_node-1].initialstate[0] + ((length1+length2)-furthest_edge)/2;
             position_y = (height1-height)/2;
             vertical_mass = ((height-2*wall_thickness)*(ribs*rib_thickness+2*wall_thickness)*(standard_length/2) + (height-2*wall_thickness)*(width-ribs*rib_thickness-2*wall_thickness)*(wall_thickness)) /(density);
             horizontal_mass = 2*(width)*(wall_thickness)*(standard_length)/(density);
             moment = (1/12)*(vertical_mass)*((2*wall_thickness+ribs*rib_thickness)*(2*wall_thickness+ribs*rib_thickness)+(height-2*wall_thickness)*(height-2*wall_thickness)) + (1/6)*(horizontal_mass)*((width)*(width)+(wall_thickness)*(wall_thickness))+(horizontal_mass)*((height+wall_thickness)/2);
             Nodes[current_node].initNodeCell(3,0,0,vertical_mass+horizontal_mass,1,stiff,damp,moment,position_x,position_y,height,wall_thickness,rib_thickness,ribs,width,current_node);
             Nodes[current_node].conncount = 0;
             break;
        }
    }
}



void Node::initNodeCell(int DoF, bool fixed_x, bool fixed_y, double mass, int NoC, double stiff,
                        double damp, double moment, double position_x, double position_y, double height,
                        double wall_thickness, double rib_thickness, int ribs, double width, int currentnode)
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
    nodeheight = height;
    conncount = 1;
    conn = new Connection[1];
    conn[0].area = (2*wall_thickness*width) + (2*wall_thickness+ribs*rib_thickness)*(nodeheight-2*wall_thickness);
    conn[0].tonode = currentnode+1;
}
