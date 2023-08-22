package alg;

import svm.SVM;

public class Project extends Algorithm{

    public Project(SVM svm){
        super(svm);

        
        if(svm.ind.V != null){
            name = "Project";
            svm.outd.algorithm = name;
            svm.outd.max_stages_count = 1;
            svm.outd.showInputData();
        }
    }

    public void run(){
        t = System.currentTimeMillis();

        
        float[] w = new float[dim+1];
        int[] n = new int[2];
        float[][] m = new float[2][dim];

        for(int i = 0; i < N; i++){
            int c = svm.ind.V[i].cl.Y;
            n[c]++;
            for(int j = 0; j < dim; j++){
                m[c][j] += svm.ind.V[i].X[j];
            }
        }
        for(int c = 0; c < 2; c++){
            for(int j = 0; j < dim; j++){
                m[c][j] /= n[c];
            }
        }

     
        for(int j = 0; j < dim; j++){
            w[j] = m[1][j] - m[0][j];
        }
        w[dim] = -(w[0] * (m[0][0] + m[1][0])/2 + w[1] * (m[0][1] + m[1][1])/2);

        
        if(dim==2){
            svm.design.setPointsOfLine(w);
        }

        svm.outd.stages_count = 1;
        svm.outd.max_stages_count = 1;
        svm.outd.computing_time = System.currentTimeMillis() - t;
        svm.outd.w = w;
        svm.outd.accuracy = getAccuracy(w);
        svm.outd.showInputData();
        svm.outd.showOutputData();

        svm.design.calculates = false;
        svm.design.repaint();
        svm.control.start.enable(false);        
    }   
}
