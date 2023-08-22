package alg;

import svm.SVM;
import io.Vector;

public class MMMPerceptron extends Algorithm{
	float[] mw = new float[dim+1];
	float mm, b0, b1;
	int S = 1000;

public MMMPerceptron(SVM svm){
	super(svm);
	if(svm.ind.V != null){
		name = "MMMPerceptron";
		svm.outd.algorithm = name;
		svm.outd.showInputData();
	}
}

public void run(){
	t = System.currentTimeMillis();
	mm = Float.MIN_VALUE;
	for(long st = 1; st <= S; st++){
		float[] w = new float[dim+1];
		for(int j = 0; j < w.length; j++) w[j] = -0.5f + (float)Math.random();
		for(long p = 1; p <= P; p++){
			boolean erori = false;
			for(int i = 0; i < N; i++) {
				float s = 0;
				for(int j = 0; j < dim; j++)
					s += w[j]*svm.ind.V[i].X[j];
				s += w[dim];
				int y = s < 0 ? 0 : 1;
				int e = svm.ind.V[i].cl.Y - y;
				if(e != 0){
					erori = true;
					for(int j = 0; j < dim; j++)
					w[j] += eta*svm.ind.V[i].X[j]*e;
					w[dim] += eta*e;
				}
			}
			if(!erori)break;
		}

		float min0 = Float.MAX_VALUE, min1 = Float.MAX_VALUE;
		int i0=-1, i1=-1;
		float B0, B1;
		for(int i = 0; i < N; i++) {
			float d = dist(w, svm.ind.V[i]);
			if(svm.ind.V[i].cl.Y == 0){
				if(d < min0){
					min0 = d;
					i0 = i;
				}
			}else{
				if(d < min1){
					min1 = d;
					i1 = i;
				}
			}
		}
		float ps = 0;
		for(int j = 0; j < dim; j++)
			ps += w[j]*svm.ind.V[i0].X[j];
		B0 = ps;
		ps = 0;
		for(int j = 0; j < dim; j++)
			ps += w[j]*svm.ind.V[i1].X[j];
		B1 = ps;
		
		float[] w0 = new float[dim+1];
		for(int j=0; j<dim; j++) w0[j] = w[j];
		w0[dim] = -B0;
		float[] w1 = new float[dim+1];
		for(int j=0; j<dim; j++) w1[j] = w[j];
		w1[dim] = -B1;	
		
		w[dim] = -B1;
		float mg = dist(w, svm.ind.V[i0]);
		if(mg > mm){
			mm = mg;
			b0 = B0;
			b1 = B1;
			System.arraycopy(w,0,mw,0,w.length);
			mw[dim] = -(b0+b1)/2;
			
			svm.control.ta.append("Stage " + st + "\n");
			svm.control.ta.append("Margin " + mg + "\n");
			String s = "";
			for(int j = 0; j < mw.length; j++) s += "w["+j+"] = " + mw[j] + "; ";
			svm.control.ta.append(s + "\n");			
			if(dim==2) svm.design.setPointsOfMaxLine(mw,w0,w1);
		}
		
	}
	
	svm.outd.stages_count = S;
	svm.outd.computing_time = System.currentTimeMillis() - t;
	svm.outd.w = mw;
	svm.outd.accuracy = getAccuracy(mw);
	svm.outd.margin = Tools.translate(mw, svm.ind.V)[3];
	svm.outd.showInputData();
	svm.outd.showOutputData();
	svm.design.calculates = false;
	svm.design.repaint();
	svm.control.start.enable(true);	

}

public float dist(float[] w, Vector V){
	float dist = 0;
	float norm = 0;
	for(int j = 0; j < w.length-1; j++){
		dist += w[j]*V.X[j];
		norm += w[j]*w[j];
	}
	dist += w[w.length-1];
	dist = Math.abs(dist);
	norm = (float)Math.sqrt(norm);
	dist /= norm;
	return dist;
}

}