int ncount;
int mcount;
Node[] nodes;
int Emin_ID;
int Emax_ID;

double Emax;
double Emin;

int garea_x0=0;
int garea_y0=0;
int garea_w=width;
int garea_h=height;
int[][] m_order;
int main_row=0, copy_row=1;

int Eaxis_wd = (width / 100) * 5;
int Eaxis_hd = (height / 100) * 5;

boolean enable_mulwidth=true;  // draw basin nodes 
boolean enable_mulcolor=false; // color the graph according to multiplicity
boolean enable_minprob=false;  // show probabilities of reaching each minumum

/* DO NOT TOUCH THIS! */
boolean nodraw=false;
boolean debug_mode=false;

int screen_width=1024;
int screen_height=700;

String input_path="/tmp/AAT.emap";
String imout_path="/tmp/output.png";

PFont font;
PFont fontE;

/**********************/

void setParams() {
    input_path="data/AAT.3.emap";
    imout_path=input_path + ".png";
    screen_width=1024;
    screen_height=700;
}

void setup() {
    setParams();

    font  = createFont("SansSerif", 10);
    fontE = createFont("SansSerif", 8);
    textFont(font);

    loadData(input_path);
    prepData();
   
    size(screen_width, screen_height);
    background(color(0,0,0));
  
    Eaxis_wd = (width / 100)  * 8;
    Eaxis_hd = (height / 100) * 5;
  
    garea_x0 = Eaxis_wd + 15;
    garea_y0 = Eaxis_hd;
    garea_w  = width  - (2 * Eaxis_wd);
    garea_h  = height - (2 * Eaxis_hd);
  
    drawAxis();

    for (int id = Emin_ID; nodes[id].PID() != id; id = nodes[id].PID()) {
	if (nodes[id].isBasin())
	    drawScale(nodes[id].getEnergy()); 
    }
   
    stroke(0);
   
    firstboot();
   
    smooth();
}

void loadData(String path) {
    String[] lines;
    String s;
    char[] sa;
    Node n;
    int offset;
 
    lines  = loadStrings(path);
    nodes  = new Node[lines.length];
    ncount = 0;
    offset = 0;
 
    while((s = lines[offset].trim()).length() > 0) {
	sa = s.toCharArray();
	if (sa[0] == '#') {
	    if(sa[1] != '!') {
		String[] hdr = splitTokens(s, "=");
         
		if (hdr.length == 2) {
		    if (hdr[0].equals("#EmaxID")) {
			Emax_ID=Integer.parseInt(hdr[1]);
			if (debug_mode) println("EmaxID=" + Emax_ID); 
		    } else if (hdr[0].equals("#EminID")) {
			Emin_ID=Integer.parseInt(hdr[1]);
			if (debug_mode) println("EminID=" + Emin_ID);
		    } else if (hdr[0].equals("#mcount")) {
			mcount=Integer.parseInt(hdr[1]);
			if (debug_mode) println("mcount=" + mcount);
		    } else if (hdr[0].equals("#Emin")) {
			Emin=Float.parseFloat(hdr[1]) - 1.0;
			if (debug_mode) println("Emin=" + Emin);
		    } else if (hdr[0].equals("#Emax")) {
			Emax=Float.parseFloat(hdr[1]) + 1.0;
			if (debug_mode) println("Emax=" + Emax);
		    }
		}
	    }
	}
	offset++;
    }
 
    if (lines[offset].length() == 0)
	offset++;
    else {
	println(path + ": Syntax error: '\\n' expected on line " + (offset + 1));
	return;
    }
 
    for (int i = offset; i < lines.length; i++) {
	if (lines[i].trim().length() > 1) {
	    n = new Node(lines[i].trim());
	    nodes[n.ID()] = n;
	    ncount++;
	}
    }
 
    compWeight(Emax_ID);
}

void compWeight_recurse(int id, double w, int l) {
    nodes[id].msumWeight   = nodes[id].msum / (double)nodes[nodes[id].PID()].msum;
    nodes[id].globalWeight = ((nodes[id].msumWeight + 1/l) * w); ///(double)l;
  
    if (nodes[id].basin) {
	for (int i = 0; i < nodes[id].multiplicity(); i++) {
	    compWeight_recurse(nodes[id].CID[i], nodes[id].globalWeight, l + 1); 
	}
    }
    nodes[id].level   = l;
    nodes[id].precomp = true;
}

void compWeight(int id) {
    compWeight_recurse(id, 1.0, 1); 
}

void prepData() {
    int last_basin=0, new_lb;
    int last_col=0, new_lc;
    boolean lb_found=false;
   
    m_order = new int[2][nodes.length];
    m_order[main_row][last_basin] = Emax_ID;
  
    while(last_basin != mcount) {
	// copy all cols up to last_basin - 1
	for (int i = 0; i < last_basin; i++) {
	    m_order[copy_row][i] = m_order[main_row][i];
	}

	// substitute basin with children
	new_lb = last_basin;
      
	if (!nodes[m_order[main_row][last_basin]].isBasin()) {
	    println("Internal error; m_order[" + main_row + "][" + last_basin + "] is not a basin!");
	    if (debug_mode) println("last_col=" + last_col);
	    if (debug_mode) println("last_basin=" + last_basin);
	    return; 
	}
      
	new_lc = nodes[m_order[main_row][last_basin]].multiplicity() + last_basin - 1;
      
	for (int i = 0; i < nodes[m_order[main_row][last_basin]].multiplicity(); i++) {
	    m_order[copy_row][last_basin + i] = nodes[m_order[main_row][last_basin]].cid(i);

	    if (!lb_found) {
		if (!nodes[m_order[copy_row][new_lb]].isBasin())
		    new_lb++;
		else
		    lb_found=true;
	    }
	}
      
	if (last_basin > last_col) {
	    if (debug_mode) println("wuut?");
	    return;
	}
      
	if (last_basin < last_col) {
	    for (int i = 0; i < (last_col - last_basin); i++) {
		m_order[copy_row][new_lc + 1 + i] = m_order[main_row][last_basin + 1 + i];
	    }
         
	    new_lc = new_lc + (last_col - last_basin);
	}  
       
	// update last_basin
	last_col   = new_lc;
	last_basin = new_lb;
     
	if (!lb_found)  {
	    while (!nodes[m_order[copy_row][last_basin]].isBasin() && last_basin < last_col) {
		++last_basin;
	    }
	}
      
	// switch main row
	main_row = (main_row + 1) % 2;
	copy_row = (copy_row + 1) % 2;
	lb_found = false;
      
	//println(last_basin);
	if (last_basin > mcount) {
	    if (debug_mode) println("fuu!\n");
	    return; 
	}
    }
   
    for (int i = 0; i < mcount; i++) {
	println(m_order[main_row][i]); 
    }
}

void drawAxis() {
    stroke(color(190,190,190));

    textAlign(CENTER, CENTER);
    text("log(E)", Eaxis_wd/2, Eaxis_hd - Eaxis_hd/3);
   
    // top scale line
    line(Eaxis_wd, Eaxis_hd,
	 Eaxis_wd + 10, Eaxis_hd);
       
    // axis line
    line(Eaxis_wd, Eaxis_hd,
	 Eaxis_wd, height - 1 - Eaxis_hd);
       
    // bottom scale line     
    line(Eaxis_wd, height - 1 - Eaxis_hd,
	 Eaxis_wd + 10, height - 1 - Eaxis_hd);
       
    stroke(0);
}

void drawScale(double E) {
    stroke(color(20, 20, 20));
    strokeWeight(0.1);

    textFont(fontE);
    textAlign(CENTER, CENTER);
    text("" + E, Eaxis_wd/2, Eaxis_hd + EtoY(E));
  
    line(Eaxis_wd, Eaxis_hd + EtoY(E),
	 Eaxis_wd + garea_w, Eaxis_hd + EtoY(E));
         
    strokeWeight(1.0);

    stroke(color(190,190,190));  
    line(Eaxis_wd, Eaxis_hd + EtoY(E),
	 Eaxis_wd + 10, Eaxis_hd + EtoY(E));
       
    stroke(0);
}


int EtoY(double E) {
    double y = garea_h * (E - Emin)/(Emax - Emin);
    
    if (E < Emin || E > Emax) {
        if (debug_mode) println("FUCK! E=" + E + "; Emin=" + Emin + "; Emax=" + Emax);
        return -100;  
    }
    
    return (int)(garea_h - y);
}

void draw() {
    if (!nodraw) {
	for (int i = 0; i < ncount; i++) {
	    nodes[i].drawAt(nodes[i].getX(), garea_y0 + EtoY(nodes[i].getEnergy()));
	} 
	smooth();
    }
    
    //if (!nodraw)
    //save("/tmp/output.png");
    
    nodraw = true;
}

void keyPressed() {
    if (key == 'S') {
	save(imout_path); 
	println("Screenshot saved at /tmp/output.png");
    }
}

void firstboot() {
    double x;
    int px;
    boolean no_unset;
    int[][] nlayer;
    int lcount;
    int lmain, lnext;
    
    if (nodraw)
	return;
   
    //background(102);
   
    for (int n = 0; n < mcount; n++) {
	x = garea_x0 + n * garea_w/mcount;
	//nodes[m_order[main_row][n]].drawAt((int)x, garea_y0 + EtoY(nodes[m_order[main_row][n]].getEnergy()));
	nodes[m_order[main_row][n]].setX((int)x);
    }
    
    nlayer = new int[2][nodes.length];
    lcount = mcount;
    lmain  = 0;
    lnext  = 1;
   
    for (int i = 0; i < mcount; i++) {
	nlayer[lmain][i] = m_order[main_row][i]; 
    }
    
    while (lcount > 0) {
	int newcnt = 0;
	int pid;
      
	for (int i = 0; i < lcount; ++i) {
	    pid = nodes[nlayer[lmain][i]].PID();
          
	    // sanity check
	    if (nodes[nlayer[lmain][i]].getX() == -1) {
		println("Internal error; getX() == -1");
		return;
	    }
          
	    if (nodes[pid].getX() == -1) {
		nodes[pid].setX(nodes[nlayer[lmain][i]].getX());
             
		//, garea_y0 + EtoY(nodes[pid].getEnergy()));
             
		nlayer[lnext][newcnt] = pid;
		newcnt++;
	    } else {
		boolean found;
             
		nodes[pid].setX((nodes[nlayer[lmain][i]].getX() + nodes[pid].getX()) / 2);
             
		//, garea_y0 + EtoY(nodes[pid].getEnergy()));
             
		found = false;
             
		for (int a = 0; a < newcnt; a++) {
		    if (nlayer[lnext][a] == pid) {
			found = true;
			break;
		    }
		}
             
		if (!found && nodes[pid].ID() != nodes[pid].PID()) {
		    if (debug_mode) println("warning: pid="+pid+" not found in the node queue!");
		    nlayer[lnext][newcnt] = pid;
		    newcnt++;
		}
	    }
	}

	lmain  = (lmain + 1) % 2;
	lnext  = (lnext + 1) % 2;
	lcount = newcnt;
    } 
}

class Node {
    boolean basin;
    int ID;
    int PID;
    int[] CID;
    double energy;
    int msum;
    int x, y;
  
    double   globalWeight;
    double   msumWeight;
    int     level;
    boolean precomp;
  
    Node(String definition) {
	String[] args = splitTokens(definition, " \t");
     
	x = y = -1;
	precomp = false;
     
	if (args[0].equals("B")) {
	    this.basin = true;
	    this.ID     = Integer.parseInt(args[1]);
	    this.energy = Double.parseDouble(args[2]);
	    this.PID    = Integer.parseInt(args[3]);
	    this.msum   = Integer.parseInt(args[4]);
	    this.CID    = new int[Integer.parseInt(args[5])];

	    for(int i = 0; i < this.CID.length; i++) {
		this.CID[i] = Integer.parseInt(args[6+i]);
	    }
	} else if (args[0].equals("M")) {
	    this.basin  = false;
	    this.ID     = Integer.parseInt(args[1]);
	    this.energy = Double.parseDouble(args[2]);
	    this.PID    = Integer.parseInt(args[3]);
	    this.msum   = 1;
	} else {
	    println("fuuuu!");
	    return; 
	}
    }
 
    int getX() {
	return this.x;
    }
 
    int getY() {
	return this.y;
    }
 
    double getEnergy() {
	return this.energy; 
    }
 
    void setX(int _x)
    {
	this.x = _x; 
    }
 
    void setY(int _y)
    {
	this.y = _y; 
    }
 
    void drawAt(int _x, int _y) {
	this.x = _x;
	this.y = _y;
	this.drawNode();
    }
  
    void drawNode() {
	//strokeWeight(0.1);
	//noStroke();
	//fill(0);
          
	//stroke(1);
    
	if (this.basin) {
	    for(int i = 0; i < this.CID.length; i++) {

		if (nodes[this.CID[i]].precomp) { 
		    //println("w=" + nodes[this.CID[i]].globalWeight);
		    strokeWeight((float)nodes[this.CID[i]].globalWeight);
		    stroke(color((int)(nodes[this.CID[i]].globalWeight * 255), 100, 50));
		} else {
		    if (debug_mode) println("precomp==false!");
		}
		//strokeWeight(nodes[this.CID[i]].msum/(float)this.msum);
              
		//println("sw=" + nodes[this.CID[i]].globalWeight);
		line(this.x,
		     garea_y0 + EtoY(this.energy),
		     nodes[this.CID[i]].getX(),
		     garea_y0 + EtoY(nodes[this.CID[i]].getEnergy()));
	    }
	    strokeWeight(1.0);
	} else {
      
	    noStroke();
	    float r = pathSum(this.ID())/(float)(nodes[Emax_ID].msum * this.level); //(float)(this.globalWeight);
	    if (debug_mode) println("r=" + r);
	    //int c = 100; //(int)(255 * this.globalWeight);
	    fill(color(255.0 * r, 0, 255.0 * (1.0 - r)));
	    //println("r=" + r);
	    ellipse(this.x, garea_y0 + EtoY(this.energy), 2, 2);
	    stroke(0);
	    fill(color(255,255,255));
	}
    }
 
    int pathSum(int id) {
	if (nodes[id].PID() == id)
	    return nodes[id].msum;
	else
	    return nodes[id].msum + pathSum(nodes[id].PID()); 
    }
 
    int ID() {
	return this.ID; 
    }
 
    int PID() {
	return this.PID; 
    }
 
    public int multiplicity() {
	if (this.basin) {
	    return this.CID.length; 
	} else
	    return 0;
    }
 
    public int cid(int i) {
	return this.CID[i]; 
    }
 
    public boolean isBasin() {
	return this.basin; 
    }
}
