package gui;

import java.awt.*;
import svm.SVM;

public class Slider extends Panel implements Runnable {
    About about;
    Thread t;
    Image im;
    Graphics img;
    int w1, h1;
    int delay = 50;
    String[] students;
    Font fnt = new Font("Arial", 0, 12);
    int y, w = 10;
    boolean init = true, control;
    int ww, hh;

    public Slider(About about, int ww, int hh) {
    	this.about = about;
		this.ww = ww;
		this.hh = hh;    		
    	setStudents();
		setFont(fnt); 
		start();
    }

	public void setStudents(){
		String[] student = new String[39];
		student[0] = "ABABEI PETRU";
		student[1] = "ANDONE BOGDAN";
		student[2] = "ARMENIA ANDREEA";
		student[3] = "B\u0102L\u0102CEANU VASILE";
		student[4] = "BENESCU MATEI";
		student[5] = "BERLAN ELENA";
		student[6] = "BOGDAN GEORGIANA";
		student[7] = "BOUR MARIANA";
		student[8] = "CIUBEREA ANDRA";
		student[9] = "C\u00CERJ\u0102 ELENA";
		student[10] = "COMORA\u0218U DIANA";
		student[11] = "CO\u0218ERU BIANCA";
		student[12] = "CONU\u021A DIANA";
		student[13] = "DAVIDESCU MARIA";
		student[14] = "DROBOT\u0102 PARASCHIV";
		student[15] = "DUGHIL\u0102 ANDREI";
		student[16] = "FRUNZ\u0102 \u0218TEFAN";
		student[17] = "GAVRIL\u0102 ANDREEA";
		student[18] = "GHIMICI IONU\u021A";
		student[19] = "GRASU VALENTIN";
		student[20] = "GRIGORIU DENISA";		
		student[21] = "HANDRAGEL DIANA";
		student[22] = "HOAMEA TEODORA";
		student[23] = "HOLCA ROBERT";
		student[24] = "HUC ALEXANDRU";
		student[25] = "IFTINCHI IOANA";
		student[26] = "IONI\u021A\u0102 ALEXANDRU";
		student[27] = "JESCU CRISTINA";
		student[28] = "LUNGU MARINA";
		student[29] = "MUTESCU ANDREI";		
		student[30] = "NECULAU SILVIU";
		student[31] = "PAVEL ANDREI";
		student[32] = "POPA ANA";
		student[33] = "POPA IOANA";
		student[34] = "PURICE ELENA";
		student[35] = "ROTAR IULIA";
		student[36] = "TUMURUC ALEXANDRU";
		student[37] = "V\u00C2RG\u0102 SABINA";
		student[38] = "VOINEA DANIEL";
		
		
		students = new String[student.length*100];
		for(int i = 0; i < students.length; i++)
			students[i] = student[i % student.length];
	}
	
    public void start() {
    	if(t == null){
    		t = new Thread(this); 
    		t.start();
	        try{Thread.sleep(1000);}
			catch(InterruptedException e) { }    			
    	}
    }
    	
    public void stop() {if(t != null){ t.stop(); t = null;}}	
	
    public void run() {
	    do {
	        repaint();
			try {Thread.sleep(delay);}
			catch(InterruptedException e) {return;}
	    } while(true);
    }  	
	
	public void reset(){
		y = hh + 10;	
		repaint();
		stop();
	}

    public final void paint(Graphics g) {
    	if(init){
			im = createImage(ww, hh);
			img = im.getGraphics();	
			for(int i = 0; i < students.length; i++) {
				FontMetrics fm = img.getFontMetrics(fnt);
				h1 += fm.getHeight();
				if(fm.stringWidth(students[i]) > w1) w1 = fm.stringWidth(students[i]);
			}
			y = hh + 10; 	
			init = false;
		}
		Color color = null;
		for(int l = 0; l < students.length; l++) {
			float f = (float)hh / 4.0F;
			float f1 = 1.0F, f2 = 1.0F, f3 = 1.0F;
			int i1 = y + (int)(1.5 * getFont().getSize() * l);
			if(i1 >= 0 && i1 <= hh) {
				float ff = 0;
				if((float)i1 <= f)
					ff = (float)i1 / f;
				else if((float)i1 >= (float)hh - f)
					ff = ((float)hh - (float)i1) / f;
				else
					ff = 1.0F;
				color = new Color((int)((float)255 * ff), (int)((float)255 * ff), (int)((float)255 * ff));
			}else color = new Color(0, 0, 0);
			img.setColor(color);
			img.setFont(fnt);
			img.drawString(students[l], w, i1);
	    }
	    g.drawImage(im, 0, 0, this);
    }

    public final void update(Graphics g) {   
    	if(img!=null){	
			img.setColor(Color.black);
			img.fillRect(0,0,ww,hh);		
		}
		if(y < -(int)((float)h1*0.75f))
			y = hh + 10;
		else
			y --;
		paint(g);
	}		

}
