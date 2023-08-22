package svm;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.net.URL;
import gui.*;
import alg.*;
import tools.*;
import io.*;


public class SVM extends Frame {
	public Toolkit tool;
	public MenuBar mb;
	public Dimension res;	      
	public Image ico, bkg, color, calculates;
	
	public Design design;
	public Settings settings;
	public SimulationControl control;
	public About about;
	public Options options;
		
	public OutputData outd;
	public InputData ind;
	
	public Algorithm algorithm;
         
public static void main (String args[]){new SVM();}


public SVM(){
	tool=getToolkit(); 
    res=tool.getScreenSize();
	loadImages();
	setIconImage(ico);
    setTitle("SVM Simulator"); 
	adaugaMenuBar(); 	
	
	design = new Design(this);
	add("Center", design);
	
	settings = new Settings(this);
	settings.resize(376, 600);
	settings.move((res.width-376)/2,(res.height-600)/2);	
	
	about = new About(this);
	about.resize(712, 410);
	about.move((res.width-712)/2,(res.height-410)/2);	
	
	control = new SimulationControl(this, 400, res.height-80);
	control.resize(400, res.height-80);
	control.move(res.width-405,35);	
	
	options = new Options(this);

	outd = new OutputData(this);	
	ind = new InputData(this);
	
    setResizable(false);	
	setBackground(settings.background_color);
	resize(res.width,res.height-40);	
    move(0,0);	
    show();   	
}

void adaugaMenuBar(){
    mb=new MenuBar();    
    Menu file = new Menu("File");
    file.add("Load Input Data");
    file.add("-");
    file.add("Exit");    
    mb.add(file);
    Menu algorithms = new Menu("Algorithms");
	algorithms.add("Median");
	Menu perceptron = new Menu("Perceptron Algorithms");
    perceptron.add("Perceptron"); 
	perceptron.add("Perceptron Dual");
	perceptron.add("MMMPerceptron");
	perceptron.add("Median > Perceptron");
    algorithms.add(perceptron);
	algorithms.add("Adaline");
	algorithms.add("MaxMarginRandom");
	algorithms.add("MaxMarginWolfe");
	algorithms.add("SMO");
	algorithms.add("Adatron");
	algorithms.add("OneVersusRest");
	algorithms.add("OneVersusOne");
	algorithms.add("Project");
	mb.add(algorithms);
    Menu view = new Menu("View");
	view.add("Show Simulation Control");
	view.add("Show Input Data");
	view.add("Show Output Data");
	view.add("-");	
	view.add("Show Cursor Coordinates"); 
	mb.add(view);
    Menu tools = new Menu("Tools");
    tools.add("Input Data Generator"); 
	tools.add("-");	
	tools.add("Settings");	
    mb.add(tools);
	Menu help = new Menu("Help");
	help.add("Help");	
	help.add("About");
	mb.add(help);
    setMenuBar(mb);    
}

public URL getResources(String s) {return this.getClass().getResource(s);}

public void loadImages(){
    try {                                 
		bkg = tool.getImage(getResources("res/bkg.jpg"));        
		ico = tool.getImage(getResources("res/ico.png"));  
		color = tool.getImage(getResources("res/color.png")); 
		calculates = tool.getImage(getResources("res/calculates.gif")); 	
    }
	catch(Throwable e) {System.out.println("Eroare la incarcarea imaginilor!");}
}

public boolean handleEvent(Event e){
    if(e.id==Event.WINDOW_DESTROY){
    	System.exit(0);
    }else if(e.id==Event.ACTION_EVENT && e.target instanceof MenuItem){
        if("Exit".equals(e.arg)){
            System.exit(0);
        }else if("Load Input Data".equals(e.arg)){
            ind.loadInputData();			
            return true; 
        }else if("Median".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new Median(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
        }else if("Perceptron".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new Perceptron(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
		}else if("MMMPerceptron".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new MMMPerceptron(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
        }else if("Perceptron Dual".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new PerceptronDual(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;			
        }else if("Median > Perceptron".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new MPerceptron(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;			
        }else if("Adaline".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}				
				//algorithm = new Adaline(this);				
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");	
			}
            return true;
        }else if("Median > Adaline".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;	
        }else if("MaxMarginRandom".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new MaxMarginRandom(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;			
		}else if("MaxMarginWolfe".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new MaxMarginWolfe(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
		}else if("SMO".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new SMO(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
		}else if("Adatron".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new Adatron(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
		}else if("OneVersusRest".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new OneVersusRest(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
		}else if("OneVersusOne".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new OneVersusOne(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
		}else if("Project".equals(e.arg)){
			if(ind.V != null){
				if(algorithm !=null){algorithm.stop(); algorithm = null; init2();}	
				algorithm = new Project(this);
				control.show();
				mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
			}
            return true;
	}else if("Show Simulation Control".equals(e.arg)){
			control.show();
			mb.getMenu(2).getItem(0).setLabel("Hide Simulation Control");
            return true;
        }else if("Hide Simulation Control".equals(e.arg)){
			control.hide();
			mb.getMenu(2).getItem(0).setLabel("Show Simulation Control");
            return true;	
		}else if("Show Input Data".equals(e.arg)){
			ind.show();
			mb.getMenu(2).getItem(1).setLabel("Hide Input Data");
            return true;
        }else if("Hide Input Data".equals(e.arg)){
			ind.hide();
			mb.getMenu(2).getItem(1).setLabel("Show Input Data");
            return true;
		}else if("Show Output Data".equals(e.arg)){
			outd.show();
			mb.getMenu(2).getItem(2).setLabel("Hide Output Data");
            return true;
        }else if("Hide Output Data".equals(e.arg)){
			outd.hide();
			mb.getMenu(2).getItem(2).setLabel("Show Output Data");
            return true;						
        }else if("Show Cursor Coordinates".equals(e.arg)){
			design.show_coords = true;
			design.repaint();
			mb.getMenu(2).getItem(4).setLabel("Hide Cursor Coordinates");
            return true;
        }else if("Hide Cursor Coordinates".equals(e.arg)){
			design.show_coords = false;
			design.repaint();
			mb.getMenu(2).getItem(4).setLabel("Show Cursor Coordinates");
            return true;			
        }else if("Input Data Generator".equals(e.arg)){
            InputDataGenerator inputDataGenerator = new InputDataGenerator(this);
            return true; 
        }else if("Settings".equals(e.arg)){
			settings.loadSettings();
			settings.show();
            return true; 
        }else if("Help".equals(e.arg)){
			File help = new File("svm/SVM.pdf");
			try{
				if (help.toString().endsWith(".pdf")) 
					Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + help);
				else {
					Desktop desktop = Desktop.getDesktop();
					desktop.open(help);
				}	
			}
			catch(IOException ex){
				System.out.println("No application registered for PDFs !");
			}
            return true; 
		}else if("About".equals(e.arg)){
			about.show();
            return true; 
        }            	
    }else return false;	
    return super.handleEvent(e);
}

public void init(){
	if(algorithm !=null){
		algorithm.stop();
		algorithm = null;
	}		
	ind.input_file = null;
	ind.V = null;
	design.show_line = false;
	design.calculates = false;
	control.init = true;
	control.start.setLabel("Start Simulation");		
	design.repaint();
}

public void init2(){
	design.show_line = false;
	design.calculates = false;
	control.init = true;
	control.start.setLabel("Start Simulation");		
	design.repaint();	
	ind.init();
}	

}
