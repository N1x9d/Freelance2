
import java.util.concurrent.CountDownLatch;

//import java.util.Arrays;

public class NewThreadE implements Runnable {

    float[][] Ez;
    float[][] Hx;
    float[][] Hy;
    float[][] e;
    int iBegin;
    int iEnd;
    float tt;
    int nt = 1;
    
    CountDownLatch cdl;

    public NewThreadE(float[][] E, float[][] H, float[][] H2, float[][] eps,
            int iBegi, int iEn, float ttt, CountDownLatch cdl) {
        //Ez = Arrays.copyOf(E, E.length);
        //Hx = Arrays.copyOf(H, H.length);
        //Hy = Arrays.copyOf(H2, H2.length);
        this.cdl = cdl;
        
        Ez = E;
        Hx = H;
        Hy = H2;
        e = eps;
        iBegin = iBegi;
        iEnd = iEn;
        tt = ttt;

        }

    @Override
    public void run() {

        
        for (int i = iBegin; i <= iEnd; i++) {        // main Ez  
            for (int j = 2; j <= tt; j++) {
                Ez[i][j] += e[i][j] * ((Hx[i][j - 1] - Hx[i][j] + Hy[i][j] - Hy[i - 1][j]));
            }
        }
        cdl.countDown();

    }
}
