package alg;
import svm.SVM;
import io.*;

public class Adaline extends Algorithm
{
	public Adaline(SVM svm)
	{
		super(svm);
		if(svm.ind.V != null) 
		{
			name = "Adaline";
			svm.outd.algorithm = name;
			svm.outd.showInputData();
		}
	}
	
	public void run()
	{
		t = System.currentTimeMillis();
		boolean flag = false;
		float b = 0, emax = 0.001;
		float[] w = new float[dim+1];
		for(int k = 0; k <w.length; k++) 
			w1[j]= w[k] = -0.5f+(float)Math.random();
			
		for(int p=1; p<=P; p++){
			double epm = 0;
			for(int i=1; i<=N;i++){
				float s = 0;
				for(int j = 0; j < dim; j++)
					s += w[j]*svm.ind.V[i].X[j];
				s+= w[dim];
				double y = s -b;
				double e = svm.ind.V[i].cl.Y - y;
				epm = epm + e*e;
				for(int j = 0; j < dim; j++)
					w[i] = w[i]+ eta*svm.ind.V[i].cl.Y * svm.ind.V[i].X[j];
				b = b - eta*svm.ind.V[i].cl.Y;
			}
			epm = epm/(2*N);
			if(epm < emax) //return w,b
		}
	}
}