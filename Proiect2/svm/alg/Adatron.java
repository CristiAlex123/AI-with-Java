package alg;

import svm.SVM;
import io.*;
import java.util.*;

public class Adatron extends Algorithm{
	SVM svm;
	
	public Adatron(SVM svm){
		super(svm);
		this.svm=svm;
		if(svm.ind.V!=null){
			name="Adatron";
			svm.outd.algorithm = name;
			svm.outd.max_stages_count =1;
			svm.outd.showInputData();
		}
	}
	
	public float dotProduct(float[] x,float[] y){
		float sum=0;
		for(int i=0;i<x.length;i++){
			sum+=x[i]*y[i];
		}
		return sum;
	}
	
	public void run(){
		t = System.currentTimeMillis();
		float Z=0;
		int i,j,k;
		float y_changed_i;
		float y_changed_j;
		float niu=0.1f;
		float C=100;
		float[] w3 = new float[dim];
		float[] w2 = new float[dim+1];
		float eta=0;
		float G=0;
		float tol=0.01f;
		float beta=0;
		boolean ok=true;
		float[] alpha=new float[N];
		for(i=0;i<N;i++){
			alpha[i]=0;
		}
	
		do{
			for(i=0;i<N;i++){
				beta=alpha[i];
				y_changed_i=svm.ind.V[i].cl.Y==0?-1:1;
				Z=0;
				for(j=0;j<N;j++){
					y_changed_j=svm.ind.V[j].cl.Y==0?-1:1;
					Z+=alpha[j]*y_changed_j*dotProduct(svm.ind.V[i].X,svm.ind.V[j].X);
				}
				G=1-y_changed_i*Z;
				if(Math.abs(G)<tol)break;
				G=G/Math.abs(G);
				eta=niu;
				while(alpha[i]+eta*G<alpha[i])eta=eta/10;
				alpha[i]=alpha[i]+eta*G;
				if(alpha[i]<0)
					alpha[i]=0;
				if(alpha[i]>C)
					alpha[i]=C;
				if(Math.abs(alpha[i]-beta)<tol)
					ok=false;
			}

		}while(ok==true);
		// Representing the support vectors
		
		ArrayList<Integer> vS = new ArrayList<Integer>();
		for(i = 0; i < N; i++){
			if(alpha[i] == Float.NaN)
				continue;
			if(alpha[i] != 0)
				vS.add(i);
		}
		for(i = 0; i < vS.size(); i++){
			int index = vS.get(i);
			int y_i = svm.ind.V[index].cl.Y == 0 ? -1 : 1;
			for(k = 0; k < dim; k++)
				w3[k] += alpha[index] * svm.ind.V[index].X[k] * y_i;
		}
		for(i = 0; i < dim; i++)
			w2[i] = w3[i];
		w2[dim]=beta;
		svm.outd.computing_time = System.currentTimeMillis() - t;
		
		float[] b = translate(w2,svm.ind.V);
		w2[dim] = -b[2];
		float[] w0 = new float[dim+1];
		for(j=0; j<dim; j++) w0[j] = w2[j];
		w0[dim] = -b[0];
		float[] w1 = new float[dim+1];
		for(j=0; j<dim; j++) w1[j] = w2[j];
		w1[dim] = -b[1];			
		if(dim==2) svm.design.setPointsOfMaxLine(w2,w0,w1);	
		
		svm.outd.w = w2;
		svm.outd.accuracy = getAccuracy(w2);
		svm.outd.margin = b[3];
		svm.outd.stages_count = 1;
		svm.outd.max_stages_count = 1;
		svm.outd.showInputData();
		svm.outd.showOutputData();
		svm.design.calculates = false;
		svm.design.repaint();
		svm.control.start.enable(false);
		
	}


 public static float[] translate(float[] w, io.Vector[] V){
		int dim = w.length - 1;
		int N = V.length;
		float max = Float.MIN_VALUE, min;
		int imax = -1, imin;
		for(int i = 0; i < N; i++){
			float d = distFromHiperplanToVector(w, V[i]);
			if(d > max){max = d; imax = i;}
		}
		int y = V[imax].cl.Y;
		float s = 0;
		for(int j = 0; j < dim; j++) s += w[j]*V[imax].X[j];
		w[dim] = -s;
		min = Float.MAX_VALUE; max = Float.MIN_VALUE;
		imax = -1; imin = -1;
		for(int i = 0; i < N; i++){
			float d = distFromHiperplanToVector(w, V[i]);
			if(V[i].cl.Y == y){
				if(d > max){max = d; imax = i;}
			}else{
				if(d < min){min = d; imin = i;}
			}
		}
		float bmax=0, bmin=0;
		s = 0;
		for(int j = 0; j < dim; j++) s += w[j]*V[imax].X[j];
		bmax = s;
		s = 0;
		for(int j = 0; j < dim; j++) s += w[j]*V[imin].X[j];
		bmin = s;
		float[] b = new float[4];
		b[0] = bmin; 
		b[1] = bmax;
		b[2] = (bmin + bmax)/2;
		w[dim] = -bmax;
		b[3] = distFromHiperplanToVector(w, V[imin]); 
		return b;
	}

	public static float distFromHiperplanToVector(float[] w, io.Vector V){
		float dist = 0, norm = 0;
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